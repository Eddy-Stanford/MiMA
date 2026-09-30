      module rrtm_radiation
!
!    Modeling an idealized Moist Atmosphere (MiMA)
!    Copyright (C) 2015  Martin Jucker
!    See https://mjucker.github.io/MiMA for documentation
!    Main reference is Jucker and Gerber, JClim 2017, doi: http://dx.doi.org/10.1175/JCLI-D-17-0127.1
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!
!   RRTM_VARS:
!   Contains all variables needed to
!   run the RRTM code, version for GCMs (hence the 'G'),
!   other than astronomy, i.e. all variables needed
!   for radiation that are not within astro.f90
!
!   external modules
        use parkind, only         : im => kind_im, rb => kind_rb
        use mima_interpolator_mod, only: interpolate_type
        use mpp_domains_mod, only: domain2d
        use restart_file_mod, only: restart_file_type, open_restart_read, open_restart_write, &
                                    close_restart, read_restart_field, write_restart_field, &
                                    restart_field_exists
!
!  rrtm_radiation variables
!
        implicit none
        private

        public :: rrtm_radiation_init, interp_temp, run_rrtmg, &
                  rrtm_precip_accum, rrtm_radiation_end

        interface read_saved
           module procedure read_saved_2d, read_saved_3d
        end interface

        logical                                    :: rrtm_init=.false.    ! has radiation been initialized?
        type(interpolate_type),save                :: o3_interp            ! use external file for ozone
        integer(kind=im)                           :: ncols_rrt,nlay_rrt   ! RRTM field sizes
                                                                           ! ncols_rrt = (size(lon)/lonstep*
                                                                           !             size(lat)
                                                                           ! nlay_rrt  = size(pfull)
        ! gas volume mixing ratios, dimensions (ncols_rrt,nlay=nlevels)
        ! vmr = mass mixing ratio [mmr,kg/kg] scaled with molecular weight [g/mol]
        ! any modification of the units must be accounted for in rrtm_xw_rad.nomica.f90, where x=(s,l)
        real(kind=rb),allocatable,dimension(:,:)   :: h2o                  ! specific humidity [kg/kg]
                                                                           ! dimension (ncols_rrt x nlay_rrt)
        real(kind=rb),allocatable,dimension(:,:)   :: o3                   ! ozone [mmr]
                                                                           ! dimension (ncols_rrt x nlay_rrt)
        real(kind=rb),allocatable,dimension(:,:)   :: co2                  ! CO2 [vmr]
                                                                           ! dimension (ncols_rrt x nlay_rrt)
        real(kind=rb),allocatable,dimension(:,:)   :: zeros                ! place holder for any species set
                                                                           !  to zero
        real(kind=rb),allocatable,dimension(:,:)   :: ones                 ! place holder for secondary species
        ! the following species are only set if use_secondary_gases=.true.
        real(kind=rb),allocatable,dimension(:,:)   :: ch4                  ! CH4 [vmr]
                                                                           ! dimension (ncols_rrt x nlay_rrt)
        real(kind=rb),allocatable,dimension(:,:)   :: n2o                  ! N2O [vmr]
                                                                           ! dimension (ncols_rrt x nlay_rrt)
        real(kind=rb),allocatable,dimension(:,:)   :: o2                   ! O2  [vmr]
                                                                           ! dimension (ncols_rrt x nlay_rrt)
        real(kind=rb),allocatable,dimension(:,:)   :: cfc11                ! CFC11 [vmr]
                                                                           ! dimension (ncols_rrt x nlay_rrt)
        real(kind=rb),allocatable,dimension(:,:)   :: cfc12                ! CFC12 [vmr]
                                                                           ! dimension (ncols_rrt x nlay_rrt)
        real(kind=rb),allocatable,dimension(:,:)   :: cfc22                ! CFC22 [vmr]
                                                                           ! dimension (ncols_rrt x nlay_rrt)
        real(kind=rb),allocatable,dimension(:,:)   :: ccl4                 ! CCL4 [vmr]
                                                                           ! dimension (ncols_rrt x nlay_rrt)
        real(kind=rb),allocatable,dimension(:,:)   :: emis                 ! surface LW emissivity per band
                                                                           ! dimension (ncols_rrt x nbndlw)
                                                                           ! =1 for black body
        ! clouds stuff
        !  cloud & aerosol optical depths, cloud and aerosol specific parameters. Set to zero
        real(kind=rb),allocatable,dimension(:,:,:) :: taucld,tauaer, sw_zro, zro_sw
        ! heating rates and fluxes, zenith angle when in-between radiation time steps
        real(kind=rb),allocatable,dimension(:,:)   :: sw_flux,lw_flux,zencos! surface fluxes, cos(zenith angle)
                                                                            ! dimension (lon x lat)
        real(kind=rb),allocatable,dimension(:,:,:) :: tdt_rad               ! heating rate [K/s]
                                                                            ! dimension (lon x lat x pfull)
        real(kind=rb),allocatable,dimension(:,:,:) :: tdt_sw_rad,tdt_lw_rad ! SW, LW radiation heating rates,
                                                                            ! diagnostics only [K/s]
                                                                            ! dimension (lon x lat x pfull)
        real(kind=rb),allocatable,dimension(:,:,:) :: t_half                ! temperature at half levels [K]
                                                                            ! dimension (lon x lat x phalf)
        real(kind=rb),allocatable,dimension(:,:)   :: rrtm_precip           ! total time of precipitation
                                                                            ! between radiation steps to
                                                                            ! determine precip_albedo
                                                                            ! dimension (lon x lat)
        real(kind=rb),allocatable,dimension(:,:)   :: olr,isr               ! Outgoing LW and net SW radiation at TOA
                                                                            ! diagnostics only [W/m2]
                                                                            ! dimension (lon x lat)
        real(kind=rb),allocatable,dimension(:,:)   :: swdn_toa,lwup_sfc     ! SW flux down at TOA and LW flux up at
                                                                            ! the surface, diagnostics only [W/m2]
                                                                            ! dimension (lon x lat)
        integer(kind=im)                           :: num_precip            ! number of times precipitation
                                                                            ! has been summed in rrtm_precip
        integer(kind=im)                           :: dt_last               ! time of last radiation calculation
                                                                            ! used for alarm
        type(domain2d)                             :: domain                ! grid domain, for the restart file
                                                                            ! RESTART/rrtm_radiation.res.nc, which
                                                                            ! holds dt_last and the fields kept
                                                                            ! between radiation steps
!---------------------------------------------------------------------------------------------------------------
! some constants
        real(kind=rb)      :: daypersec=1./86400,deg2rad
! no clouds in the radiative scheme
        integer(kind=im) :: icld=0,idrv=0, &
             inflglw=0,iceflglw=0,liqflglw=0, &
             iaer=0
!---------------------------------------------------------------------------------------------------------------
!                                namelist values
!---------------------------------------------------------------------------------------------------------------
! input files: file names are always given without '.nc', which is always assumed
!  the field to be read within the file needs to have the same name as the file
        logical            :: do_read_ozone=.true.           ! read ozone from an external file?
                                                              !  this is the only way to get ozone into the model
        character(len=256) :: ozone_file='ozone_1990'              !  file name of ozone file to read
        real(kind=rb)      :: scale_ozone = 1.0               ! scale the ozone values in the file by this factor
        real(kind=rb)      :: o3_val = 0.0                    ! if do_read_ozone = .false., give ozone this constant value
! secondary gases (CH4,N2O,O2,CFC11,CFC12,CFC22,CCL4)
        logical            :: include_secondary_gases=.false. ! non-zero values for above listed secondary gases?
        real(kind=rb)      :: ch4_val  = 0.                   !  if .true., value for CH4
        real(kind=rb)      :: n2o_val  = 0.                   !                       N2O
        real(kind=rb)      :: o2_val   = 0.                   !                       O2
        real(kind=rb)      :: cfc11_val= 0.                   !                       CFC11
        real(kind=rb)      :: cfc12_val= 0.                   !                       CFC12
        real(kind=rb)      :: cfc22_val= 0.                   !                       CFC22
        real(kind=rb)      :: ccl4_val = 0.                   !                       CCL4
! some safety boundaries
        real(kind=rb)      :: h2o_lower_limit = 2.e-7         ! never use smaller than this in radiative scheme
        real(kind=rb)      :: temp_lower_limit = 100.         ! never go below this in radiative scheme
        real(kind=rb)      :: temp_upper_limit = 370.         ! never go above this in radiative scheme
! primary gases: CO2 and H2O
        real(kind=rb)      :: co2ppmv=390.                    ! CO2 ppmv concentration
        logical            :: do_zm_tracers=.false.           ! Feed only the zonal mean of tracers to radiation

! radiation time stepping and spatial sampling
        integer(kind=im)   :: dt_rad=4500                        ! Radiation time step - every step if dt_rad<dt_atmos
        logical            :: store_intermediate_rad =.true.  ! Keep rad constant over entire dt_rad?
                                                              ! Else only heat radiatively at every dt_rad
        logical            :: do_rad_time_avg =.true.         ! Average coszen for SW radiation over dt_rad?
        integer(kind=im)   :: dt_rad_avg = 4500             ! If averaging, over what time?
                                                              !  no averaging if dt_rad_avg = 0. (equivalent to do_rad_time_avg=.false.)
                                                              !  dt_rad_avg=dt_rad if dt_rad_avg < 0
                                                              !  Default is to average over the whole day, i.e. remove diurnal  cycle.
                                                              !  This seems safest as the diurnal cycle has been observed
                                                              !  to create strong atmospheric tides with topography.
        integer(kind=im)   :: lonstep=4                       ! Subsample fields along longitude
                                                              !  for faster radiation calculation
! some fancy radiation tweaks
        real(kind=rb)      :: slowdown_rad = 1.0              ! factor do simulate slower seasonal cycle: >1 means faster, <1 slower
        logical            :: do_zm_rad=.false.               ! Only compute zonal mean radiation
        logical            :: do_precip_albedo=.false.        ! Modify albedo depending on large scale
                                                              !  precipitation (crude cloud parameterization)
        real(kind=rb)      :: precip_albedo=0.35              ! If so, what's the cloud albedo?
        real(kind=rb)      :: precip_lat = 0.0                ! If so, poleward of which latitude should it be applied?
        character(len=14)  :: precip_albedo_mode = 'full'     ! If so, use
                                                              !  full precipitation ('full')
                                                              !  only large scale condensation ('lscale')
                                                              !  only convection ('conv')
!---------------------------------------------------------------------------------------------------------------
!
!-------------------- diagnostics fields -------------------------------

        integer :: id_tdt_rad,id_tdt_sw,id_tdt_lw,id_coszen,id_flux_sw,id_flux_lw,id_albedo,id_ozone,id_thalf
        integer :: id_olr,id_isr,id_swdn_toa,id_lwup_sfc
        character(len=9), parameter :: mod_name = 'radiation'
        real :: missing_value = -999.

!---------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------

        namelist/rrtm_radiation_nml/ include_secondary_gases, do_read_ozone, ozone_file, scale_ozone, o3_val, &
             &ch4_val, n2o_val, o2_val, cfc11_val, cfc12_val, cfc22_val, ccl4_val, &
             &h2o_lower_limit,temp_lower_limit,temp_upper_limit,co2ppmv, &
             &slowdown_rad, &
             &store_intermediate_rad, do_rad_time_avg, dt_rad, dt_rad_avg, &
             &lonstep, do_zm_tracers, do_zm_rad, &
             &do_precip_albedo, precip_albedo_mode, precip_albedo, precip_lat


      contains

!*****************************************************************************************
        subroutine rrtm_radiation_init(axes,Time,ncols,nlay,lonb,latb,domain_in)
!
! Initialize diagnostics, allocate variables, set constants
!
! Modules
          use rrtm_astro, only:       astro_init,solday
          use parrrtm, only:          nbndlw
          use parrrsw, only:          nbndsw
          use diag_manager_mod, only: register_diag_field, send_data
          use mima_interpolator_mod, only: interpolate_type, interpolator_init, &
                                      &CONSTANT, ZERO,INTERP_WEIGHTED_P
          use fms_mod, only:          input_nml_file, check_nml_error,  &
                                      &mpp_pe, mpp_root_pe, &
                                      &write_version_number, stdlog, &
                                      &error_mesg, NOTE, WARNING
          use time_manager_mod, only: time_type
! Local variables
          implicit none

          integer, intent(in), dimension(4) :: axes
          type(time_type), intent(in)       :: Time
          integer(kind=im),intent(in)       :: ncols,nlay
          real(kind=rb),dimension(:),intent(in) :: lonb,latb
          type(domain2d),intent(in)         :: domain_in   ! domain decomposition of the model grid

          integer :: i,k,seconds

          integer :: ierr, io, unit


! read namelist and copy to logfile
          read (input_nml_file, nml=rrtm_radiation_nml, iostat=io)
          ierr = check_nml_error (io, 'rrtm_radiation_nml')

          !call write_version_number ( version, tagname )
          if ( mpp_pe() == mpp_root_pe() ) then
             write (stdlog(), nml=rrtm_radiation_nml)
          endif
!----
!------------ initialize diagnostic fields ---------------

          id_tdt_rad = &
               register_diag_field ( mod_name, 'tdt_rad', axes(1:3), Time, &
                 'Temperature tendency due to radiation', &
                 'K/s', missing_value=missing_value               )
          id_tdt_sw = &
               register_diag_field ( mod_name, 'tdt_sw', axes(1:3), Time, &
                 'Temperature tendency due to SW radiation', &
                 'K/s', missing_value=missing_value               )
          id_tdt_lw = &
               register_diag_field ( mod_name, 'tdt_lw', axes(1:3), Time, &
                 'Temperature tendency due to LW radiation', &
                 'K/s', missing_value=missing_value               )
          id_coszen  = &
               register_diag_field ( mod_name, 'coszen', axes(1:2), Time, &
                 'cosine of zenith angle', &
                 'none', missing_value=missing_value               )
          id_flux_sw = &
               register_diag_field ( mod_name, 'swnet_sfc', axes(1:2), Time, &
                 'Net SW flux at surface (positive down)', &
                 'W/m2', missing_value=missing_value               )
          id_flux_lw = &
               register_diag_field ( mod_name, 'lwdn_sfc', axes(1:2), Time, &
                 'LW flux down at surface', &
                 'W/m2', missing_value=missing_value               )
          id_olr     = &
               register_diag_field ( mod_name, 'olr', axes(1:2), Time, &
                 'Outgoing longwave radiation at TOA', &
                 'W/m2', missing_value=missing_value               )
          id_isr     = &
               register_diag_field ( mod_name, 'swnet_toa', axes(1:2), Time, &
                 'Net SW flux at TOA (positive down)', &
                 'W/m2', missing_value=missing_value               )
          id_swdn_toa = &
               register_diag_field ( mod_name, 'swdn_toa', axes(1:2), Time, &
                 'SW flux down at TOA', &
                 'W/m2', missing_value=missing_value               )
          id_lwup_sfc = &
               register_diag_field ( mod_name, 'lwup_sfc', axes(1:2), Time, &
                 'LW flux up at surface', &
                 'W/m2', missing_value=missing_value               )
          id_albedo  = &
               register_diag_field ( mod_name, 'albedo_rad', axes(1:2), Time, &
                 'Surface albedo seen by the radiation', &
                 'none', missing_value=missing_value               )
          id_ozone   = &
               register_diag_field ( mod_name, 'ozone', axes(1:3), Time, &
                 'Ozone', &
                 'mmr', missing_value=missing_value               )
          id_thalf   = &
               register_diag_field ( mod_name, 'thalf', (/axes(1),axes(2),axes(4)/), Time, &
                 'Temperature on half levels', &
                 'K', missing_value=missing_value               )
!
!------------ make sure namelist choices are consistent -------
! this does not work at the moment, as dt_atmos from coupler_mod induces a circular dependency at compilation
!          if(dt_rad .le. dt_atmos .and. store_intermediate_rad)then
!             call error_mesg ( 'rrtm_gases_init', &
!                  ' dt_rad <= dt_atmos, for conserving memory, I am setting store_intermediate_rad=.false.', &
!                  WARNING)
!             store_intermediate_rad = .false.
!          endif
!          if(dt_rad .gt. dt_atmos .and. .not.store_intermediate_rad)then
!             call error_mesg( 'rrtm_gases_init', &
!                  ' dt_rad > dt_atmos, but store_intermediate_rad=.false. might cause time steps with zero radiative forcing!', &
!                  WARNING)
!          endif

!------------ set some constants and parameters -------

          deg2rad = acos(0.)/90.

          dt_last = -dt_rad !make sure we are computing radiation at the first time step

          ncols_rrt = ncols/lonstep
          nlay_rrt  = nlay

          if(dt_rad_avg .eq. 0.0) then
             do_rad_time_avg = .false.
          else if(dt_rad_avg .lt. 0) then
             dt_rad_avg = dt_rad
          endif

!------------ allocate arrays to be used later  -------
          allocate(t_half(size(lonb,1)-1,size(latb)-1,nlay+1))

             allocate(h2o(ncols_rrt,nlay_rrt),o3(ncols_rrt,nlay_rrt), &
                  co2(ncols_rrt,nlay_rrt))
             allocate(ones(ncols_rrt,nlay_rrt), &
                  zeros(ncols_rrt,nlay_rrt))
             allocate(emis(ncols_rrt,nbndlw))
             allocate(taucld(nbndlw,ncols_rrt,nlay_rrt), &
                  tauaer(ncols_rrt,nlay_rrt,nbndlw))
             allocate(sw_zro(nbndsw,ncols_rrt,nlay_rrt), &
                  zro_sw(ncols_rrt,nlay_rrt,nbndsw))
             if(id_coszen > 0)allocate(zencos (size(lonb,1)-1,size(latb,1)-1))

             ! gases
             h2o   = 0.     !this will be set by the water vapor tracer
             o3    = o3_val !this will be set by an input file if do_read_ozone=.true.
             co2   = co2ppmv*1.e-6 ! convert ppmv
             zeros = 0. ! gases and clouds
             ones  = 1. ! gases and clouds

             emis  = 1. !black body: 1.0

             ! absorption
             taucld = 0.
             tauaer = 0.
             ! clouds
             sw_zro = 0.
             zro_sw = 0.

          if(do_read_ozone)then
             call interpolator_init (o3_interp, trim(ozone_file)//'.nc', lonb, latb, data_out_of_bounds=(/ZERO/))
          endif

          if(store_intermediate_rad .or. id_flux_sw > 0) &
               allocate(sw_flux(size(lonb,1)-1,size(latb,1)-1))
          if(store_intermediate_rad .or. id_flux_lw > 0) &
               allocate(lw_flux(size(lonb,1)-1,size(latb,1)-1))
          if(do_precip_albedo)allocate(rrtm_precip(size(lonb,1)-1,size(latb,1)-1))
          if(store_intermediate_rad .or. id_tdt_rad > 0)&
               allocate(tdt_rad(size(lonb,1)-1,size(latb,1)-1,nlay))
          if(id_tdt_sw .gt. 0) allocate(tdt_sw_rad(size(lonb,1)-1,size(latb,1)-1,nlay))
          if(id_tdt_lw .gt. 0) allocate(tdt_lw_rad(size(lonb,1)-1,size(latb,1)-1,nlay))
          if(id_isr .gt. 0) allocate(isr(size(lonb,1)-1,size(latb,1)-1))
          if(id_olr .gt. 0) allocate(olr(size(lonb,1)-1,size(latb,1)-1))
          if(id_swdn_toa .gt. 0) allocate(swdn_toa(size(lonb,1)-1,size(latb,1)-1))
          if(id_lwup_sfc .gt. 0) allocate(lwup_sfc(size(lonb,1)-1,size(latb,1)-1))

          if(do_precip_albedo)then
             rrtm_precip = 0.
             num_precip  = 0
          endif

          domain = domain_in
          call read_restart_rrtm

          call astro_init

          if(solday .gt. 0)then
             call error_mesg( mod_name, &
                  ' running perpetual simulation', NOTE)
          endif

          rrtm_init=.true.

        end subroutine rrtm_radiation_init
!*****************************************************************************************
        subroutine interp_temp(z_full,z_half,t_surf_rad,t)
          implicit none

          real(kind=rb),dimension(:,:,:),intent(in)  :: z_full,z_half,t
          real(kind=rb),dimension(:,:)  ,intent(in)  :: t_surf_rad

          integer i,j,k,kend
          real dzk,dzk1,dzk2

! also, for some reason, z_half(k=1)=0. so we need to deal with k=1 separately
          kend=size(z_full,3)
          do k=2,kend
             do j=1,size(t,2)
                do i=1,size(t,1)
                   dzk2 = 1./( z_full(i,j,k-1)   - z_full(i,j,k) )
                   dzk  = ( z_half(i,j,k  )   - z_full(i,j,k) )*dzk2
                   dzk1 = ( z_full(i,j,k-1)   - z_half(i,j,k) )*dzk2
                   t_half(i,j,k) = t(i,j,k)*dzk1 + t(i,j,k-1)*dzk
                enddo
             enddo
          enddo
! top of the atmosphere: need to extrapolate. z_half(1)=0, so need to use values on full grid
          do j=1,size(t,2)
             do i=1,size(t,1)
                !standard linear extrapolation
                !top: use full points, and distance is 1.5 from k=2
                t_half(i,j,1) = 0.5*(3*t(i,j,1)-t(i,j,2))
                !bottom: z=0 => distance is -z_full(kend-1)/(z_full(kend)-z_full(kend-1))
!!$                t_half(i,j,kend+1) = t(i,j,kend-1) &
!!$                     + (z_half(i,j,kend+1) - z_full(i,j,kend-1))&
!!$                     * (t     (i,j,kend  ) - t     (i,j,kend-1))&
!!$                     / (z_full(i,j,kend  ) - z_full(i,j,kend-1))
                !bottom: t_half = t_surf
                t_half(i,j,kend+1) = t_surf_rad(i,j)
             enddo
          enddo


        end subroutine interp_temp
!*****************************************************************************************
!*****************************************************************************************
        subroutine run_rrtmg(is,js,Time,Time_diag,lat,lon,p_full,p_half,albedo,q,t,t_surf_rad,tdt,coszen,flux_sw,flux_lw)
!
! Driver for RRTMG radiation scheme.
! Prepares all inputs, calls SW and LW radiation schemes,
!  transforms outputs back into FMS form
!
! Modules
          use fms_mod, only:         error_mesg, FATAL
          use mpp_mod, only:         mpp_pe,mpp_root_pe
          use rrtmg_lw_rad, only:    rrtmg_lw
          use rrtmg_sw_rad, only:    rrtmg_sw
          use rrtm_astro, only:      compute_zenith,use_dyofyr,solr_cnst,&
                                     solrad,solday,equinox_day
          use time_manager_mod,only: time_type,get_time,set_time
          use mima_interpolator_mod,only: interpolator
!---------------------------------------------------------------------------------------------------------------
! In/Out variables
          implicit none

          integer, intent(in)                               :: is, js          ! index range for each CPU
          type(time_type),intent(in)                        :: Time            ! global time in calendar
          type(time_type),intent(in)                        :: Time_diag       ! time the diagnostics are sent at
                                                                               ! (Time_next, as for the other physics)
          real(kind=rb),dimension(:,:,:),intent(in)         :: p_full,p_half   ! pressure, full and half levels
                                                                               ! dimension (lat x lon x p*)
          real(kind=rb),dimension(:,:,:),intent(in)         :: q               ! water vapor mixing ratio [g/g]
                                                                               ! dimension (lat x lon x pfull)
          real(kind=rb),dimension(:,:,:),intent(in)         :: t               ! temperature [K]
                                                                               ! dimension (lat x lon x pfull)
          real(kind=rb),dimension(:,:),intent(in)           :: lat,lon         ! latitude, longitude
                                                                               ! dimension (lat x lon)
          real(kind=rb),dimension(:,:),intent(in)           :: albedo          ! surface albedo
                                                                               ! dimension (lat x lon)
          real(kind=rb),dimension(:,:),intent(in)           :: t_surf_rad      ! surface temperature [K]
                                                                               ! dimension (lat x lon)
          real(kind=rb),dimension(:,:,:),intent(inout)      :: tdt             ! heating rate [K/s]
                                                                               ! dimension (lat x lon x pfull)
          real(kind=rb),dimension(:,:),intent(out)          :: coszen          ! cosine of zenith angle
                                                                               ! dimension (lat x lon)
          real(kind=rb),dimension(:,:),intent(out),optional :: flux_sw,flux_lw ! surface fluxes [W/m2]
                                                                               ! dimension (lat x lon)
                                                                               ! need to have both or none!
!---------------------------------------------------------------------------------------------------------------
! Local variables
          integer k,j,i,ij,j1,i1,ij1,kend,dyofyr,seconds,days
          integer si,sj,sk,locmin(3)
          real(kind=rb),dimension(size(q,1),size(q,2),size(q,3)) :: o3f
          real(kind=rb),dimension(ncols_rrt,nlay_rrt) :: pfull,tfull&
               , hr,hrc, swhr, swhrc
          real(kind=rb),dimension(size(tdt,1),size(tdt,2),size(tdt,3)) :: tdt_rrtm
          real(kind=rb),dimension(ncols_rrt,nlay_rrt+1) :: uflx, dflx, uflxc, dflxc&
               ,swuflx, swdflx, swuflxc, swdflxc
          real(kind=rb),dimension(size(q,1)/lonstep,size(q,2),size(q,3)  ) :: swijk,lwijk
          real(kind=rb),dimension(size(q,1)/lonstep,size(q,2)) :: swflxijk,lwflxijk,olrijk,isrijk
          real(kind=rb),dimension(size(q,1)/lonstep,size(q,2)) :: swdntoaijk,lwupsfcijk
          real(kind=rb),dimension(ncols_rrt,nlay_rrt+1):: phalf,thalf
          real(kind=rb),dimension(ncols_rrt)   :: tsrf,cosz_rr,albedo_rr
          real(kind=rb) :: dlon,dlat,dj,di
          type(time_type) :: Time_loc
          real(kind=rb),dimension(size(q,1),size(q,2)) :: albedo_loc
          real(kind=rb),dimension(size(q,1),size(q,2),size(q,3)) :: q_tmp
! debug
          integer :: indx2(2),indx(3),ii,ji,ki
!---------------------------------------------------------------------------------------------------------------

          if(.not. rrtm_init)&
               call error_mesg('run_rrtm','module not initialized', FATAL)

!check if we really want to recompute radiation (alarm,input file(s))
! alarm
          call get_time(Time,seconds,days)
          if(seconds < dt_last) dt_last=dt_last-86400 !it's a new day
          if(seconds - dt_last .ge. dt_rad) then
             dt_last = seconds
          else
             if(store_intermediate_rad)then
                tdt_rrtm = tdt_rad
                flux_sw = sw_flux
                flux_lw = lw_flux
             else
                tdt_rrtm = 0.
                flux_sw  = 0.
                flux_lw  = 0.
             endif
             tdt = tdt + tdt_rrtm
             call write_diag_rrtm(Time_diag,is,js)
             return !not time yet
          endif
!make sure we run perpetual when solday > 0)
          if(solday > 0)then
             Time_loc = set_time(seconds,solday)
          elseif(slowdown_rad .ne. 1.0)then
             seconds = days*86400 + seconds
             Time_loc = set_time(int(seconds*slowdown_rad))
          else
             Time_loc = Time
          endif
!
! compute zenith angle
!  this is also an output, so need to compute even if we read radiation from file
          if(do_rad_time_avg) then
             call compute_zenith(Time_loc,equinox_day,dt_rad_avg,lat,lon,coszen,dyofyr)
          else
             call compute_zenith(Time_loc,equinox_day,0     ,lat,lon,coszen,dyofyr)
          end if
!---------------------------------------------------------------------------------------------
! we know now that we want to run radiation

          si=size(tdt,1)
          sj=size(tdt,2)
          sk=size(tdt,3)

          if(.not. use_dyofyr) dyofyr=0 !use solrad instead of day of year

          !get ozone
          if(do_read_ozone)then
             call interpolator( o3_interp, Time_loc, p_half, o3f, trim(ozone_file))
             o3f = o3f*scale_ozone
             !due to interpolation, some values might be negative
             o3f = max(0.0,o3f)
          else
             o3f = o3_val   ! only used for the ozone diagnostic; RRTM uses o3 = o3_val
          endif

          !interactive albedo: zonal mean of precipitation
          if(do_precip_albedo .and. num_precip>0)then
             where ( abs(lat) < precip_lat*3.14159265/180. ) rrtm_precip = 0.
             do i=1,size(albedo,1)
                albedo_loc(i,:) = albedo(i,:) + (precip_albedo - albedo(i,:))&
                     &*sum(rrtm_precip,1)/size(rrtm_precip,1)/num_precip
             enddo
             rrtm_precip = 0.
             num_precip  = 0
          else
             albedo_loc = albedo
          endif
!---------------------------------------------------------------------------------------------------------------
          !Compute zonal means if that's what we want to feed to RRTM
          if(do_zm_tracers)then
             do i=1,size(q,1)
                q_tmp(i,:,:) = sum(q,1)/size(q,1)
             enddo
          else
             q_tmp = q
          endif
!---------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------
          !RRTM's first pressure level is at the surface - need to inverse order
          !also, RRTM's pressures are in hPa
          !reshape arrays
          pfull = reshape(p_full(1:si:lonstep,:,sk  :1:-1),(/ si*sj/lonstep,sk   /))*0.01
          phalf = reshape(p_half(1:si:lonstep,:,sk+1:1:-1),(/ si*sj/lonstep,sk+1 /))*0.01
          !for RRTM, we need the top level to be greater than 0
          if(minval(phalf(:,sk+1)) .le. 0.) &
               &phalf(:,sk+1) = pfull(:,sk)*0.5
          tfull = reshape(t     (1:si:lonstep,:,sk  :1:-1),(/ si*sj/lonstep,sk   /))
          thalf = reshape(t_half(1:si:lonstep,:,sk+1:1:-1),(/ si*sj/lonstep,sk+1 /))
          h2o   = reshape(q_tmp (1:si:lonstep,:,sk  :1:-1),(/ si*sj/lonstep,sk   /))
          if(do_read_ozone)o3 = reshape(o3f(1:si:lonstep,:,sk :1:-1),(/ si*sj/lonstep,sk  /))

          cosz_rr   = reshape(coszen    (1:si:lonstep,:),(/ si*sj/lonstep /))
          albedo_rr = reshape(albedo_loc(1:si:lonstep,:),(/ si*sj/lonstep /))
          tsrf      = reshape(t_surf_rad(1:si:lonstep,:),(/ si*sj/lonstep /))

!---------------------------------------------------------------------------------------------------------------
! now actually run RRTM
!
          swhr = 0.
          swdflx = 0.
          swuflx = 0.
          ! make sure we don't go beyond 'dangerous' values for radiation. Same as Merlis spectral_am2rad
          h2o   = max(h2o  , h2o_lower_limit)
          tfull = max(tfull, temp_lower_limit)
          tfull = min(tfull, temp_upper_limit)
          thalf = max(thalf, temp_lower_limit)
          thalf = min(thalf, temp_upper_limit)


          !h2o=h2o_lower_limit
          ! SW seems to have a problem with too small coszen values.
          ! anything lower than 0.01 (about 15min) is set to zero
          !where(cosz_rr < 1.e-2)cosz_rr=0.

          if(include_secondary_gases)then
             call rrtmg_sw &
                  (ncols_rrt, nlay_rrt , icld     , iaer         , &
                  pfull     , phalf    , tfull    , thalf        , tsrf         , &
                  h2o       , o3       , co2      , ch4_val*ones , n2o_val*ones , o2_val*ones , &
                  albedo_rr , albedo_rr, albedo_rr, albedo_rr, &
                  cosz_rr   , solrad   , dyofyr   , solr_cnst, &
                  inflglw   , iceflglw , liqflglw , &
                  ! cloud parameters
                  zeros     , taucld   , sw_zro   , sw_zro   , sw_zro , &
                  zeros     , zeros    , 10*ones  , 10*ones  , &
                  tauaer    , zro_sw   , zro_sw   , zro_sw    , &
                  ! output
                  swuflx    , swdflx   , swhr     , swuflxc  , swdflxc, swhrc)
          else
             call rrtmg_sw &
                  (ncols_rrt, nlay_rrt , icld     , iaer     , &
                  pfull     , phalf    , tfull    , thalf    , tsrf , &
                  h2o       , o3       , co2      , zeros    , zeros, zeros, &
                  albedo_rr , albedo_rr, albedo_rr, albedo_rr, &
                  cosz_rr   , solrad   , dyofyr   , solr_cnst, &
                  inflglw   , iceflglw , liqflglw , &
                  ! cloud parameters
                  zeros     , taucld   , sw_zro   , sw_zro   , sw_zro , &
                  zeros     , zeros    , 10*ones  , 10*ones  , &
                  tauaer    , zro_sw   , zro_sw   , zro_sw   , &
                  ! output
                  swuflx    , swdflx   , swhr     , swuflxc  , swdflxc, swhrc)
          endif

          swijk   = reshape(swhr(:,sk:1:-1),(/ si/lonstep,sj,sk /))*daypersec
          isrijk  = reshape(swdflx(:,sk+1)-swuflx(:,sk+1),(/ si/lonstep,sj /))
          swdntoaijk = reshape(swdflx(:,sk+1),(/ si/lonstep,sj /))

          hr = 0.
          dflx = 0.
          uflx = 0.
          if(include_secondary_gases)then
             call rrtmg_lw &
                  (ncols_rrt     , nlay_rrt       , icld           , idrv , &
                  pfull          , phalf          , tfull          , thalf, tsrf  , &
                  h2o            , o3             , co2            , &
                  ! secondary gases
                  ch4_val*ones   , n2o_val*ones   , o2_val*ones    , &
                  cfc11_val*ones , cfc12_val*ones , cfc22_val*ones , ccl4_val*ones , &
                  ! emissivity and cloud composition
                  emis           , inflglw        , iceflglw       , liqflglw      ,  &
                  ! cloud parameters
                  zeros          , taucld         , zeros          , zeros         , 10*ones, 10*ones, &
                  tauaer         , &
                  ! output
                  uflx           , dflx           , hr             , uflxc         , dflxc  , hrc)
          else
             call rrtmg_lw &
                  (ncols_rrt, nlay_rrt, icld    , idrv , &
                  pfull     , phalf   , tfull   , thalf, tsrf , &
                  h2o       , o3      , co2     , zeros, zeros, zeros, &
                  zeros     , zeros   , zeros   , zeros, &
                  ! emissivity and cloud composition
                  emis      , inflglw , iceflglw, liqflglw, &
                  ! cloud parameters
                  zeros     , taucld  , zeros   , zeros, 10*ones, 10*ones, &
                  tauaer    , &
                  ! output
                  uflx      , dflx    , hr      , uflxc, dflxc  , hrc)
          endif

          lwijk   = reshape(hr(:,sk:1:-1),(/ si/lonstep,sj,sk /))*daypersec
          olrijk  = reshape(uflx(:,sk+1),(/ si/lonstep,sj /))
          lwupsfcijk = reshape(uflx(:,1),(/ si/lonstep,sj /))

!---------------------------------------------------------------------------------------------------------------
          ! get radiation
! interpolate back onto GCM grid (latitude is kept the same due to parallelisation)
             dlon=1./lonstep
             do i=1,size(swijk,1)
                i1 = i+1
                ! close toroidally
                if(i1 > size(swijk,1)) i1=1
                do ij=1,lonstep
                   di = (ij-1)*dlon
                   ij1 = (i-1)*lonstep + ij
                   if(do_zm_rad) then
                      tdt_rrtm(ij1,:,:) = sum(swijk+lwijk,1)/max(1,size(swijk,1))
                   else
                      tdt_rrtm(ij1,:,:) =  &
                           +       di *(swijk(i1,:,:) + lwijk(i1,:,:)) &
                           +   (1.-di)*(swijk(i ,:,:) + lwijk(i ,:,:))
                   endif
                   if(id_tdt_sw .gt. 0)tdt_sw_rad(ij1,:,:)=di*swijk(i1,:,:)+(1.-di)*swijk(i,:,:)
                   if(id_tdt_lw .gt. 0)tdt_lw_rad(ij1,:,:)=di*lwijk(i1,:,:)+(1.-di)*lwijk(i,:,:)
                   if(id_olr    .gt. 0)olr(ij1,:)         =di*olrijk(i1,:) +(1.-di)*olrijk(i,:)
                   if(id_isr    .gt. 0)isr(ij1,:)         =di*isrijk(i1,:) +(1.-di)*isrijk(i,:)
                   if(id_swdn_toa .gt. 0)swdn_toa(ij1,:)  =di*swdntoaijk(i1,:)+(1.-di)*swdntoaijk(i,:)
                   if(id_lwup_sfc .gt. 0)lwup_sfc(ij1,:)  =di*lwupsfcijk(i1,:)+(1.-di)*lwupsfcijk(i,:)
                enddo
             enddo
          tdt = tdt + tdt_rrtm
          ! store radiation between radiation time steps
          if(store_intermediate_rad .or. id_tdt_rad > 0) tdt_rad = tdt_rrtm


          ! get the surface fluxes
          if(present(flux_sw).and.present(flux_lw))then
             !only surface fluxes are needed
             swflxijk = reshape(swdflx(:,1)-swuflx(:,1),(/ si/lonstep,sj /)) ! net down SW flux
             lwflxijk = reshape(  dflx(:,1)            ,(/ si/lonstep,sj /)) ! down LW flux
             dlon=1./lonstep
             do i=1,size(swijk,1)
                i1 = i+1
                ! close toroidally
                if(i1 > size(swijk,1)) i1=1
                do ij=1,lonstep
                   di = (ij-1)*dlon
                   ij1 = (i-1)*lonstep + ij
                   if(do_zm_rad) then
                      flux_sw(ij1,:) = sum(swflxijk,1)/max(1,size(swflxijk,1))
                      flux_lw(ij1,:) = sum(lwflxijk,1)/max(1,size(lwflxijk,1))
                   else
                      flux_sw(ij1,:) = di*swflxijk(i1,:) + (1.-di)*swflxijk(i ,:)
                      flux_lw(ij1,:) = di*lwflxijk(i1,:) + (1.-di)*lwflxijk(i ,:)
                   endif
                enddo
             enddo
             ! store between radiation steps
             if(store_intermediate_rad)then
                sw_flux = flux_sw
                lw_flux = flux_lw
             else
                if(id_flux_sw > 0)sw_flux = flux_sw
                if(id_flux_lw > 0)lw_flux = flux_lw
             endif
             if(id_coszen  > 0)zencos  = coszen
          endif

          ! check if we want surface albedo as a function of precipitation
          !  call diagnostics accordingly
          if(do_precip_albedo)then
             call write_diag_rrtm(Time_diag,is,js,o3f,t_half,albedo_loc)
          else
             call write_diag_rrtm(Time_diag,is,js,o3f,t_half)
          endif
        end subroutine run_rrtmg

!*****************************************************************************************
!*****************************************************************************************
        subroutine write_diag_rrtm(Time,is,js,ozone,thalf,albedo_loc)
!
! write out diagnostics fields
!
! Modules
          use diag_manager_mod, only: register_diag_field, send_data
          use time_manager_mod,only:  time_type
! Input variables
          implicit none
          type(time_type)               ,intent(in)          :: Time
          integer                       ,intent(in)          :: is, js
          real(kind=rb),dimension(:,:,:),intent(in),optional :: ozone,thalf
          real(kind=rb),dimension(:,:  ),intent(in),optional :: albedo_loc
! Local variables
          logical :: used

!------- temperature tendency due to radiation ------------
          if ( id_tdt_rad > 0 ) then
             used = send_data ( id_tdt_rad, tdt_rad, Time, is, js, 1 )
          endif
!------- temperature tendency due to SW radiation ---------
          if ( id_tdt_sw > 0 ) then
             used = send_data ( id_tdt_sw, tdt_sw_rad, Time, is, js, 1 )
          endif
!------- temperature tendency due to LW radiation ---------
          if ( id_tdt_lw > 0 ) then
             used = send_data ( id_tdt_lw, tdt_lw_rad, Time, is, js, 1 )
          endif
!------- cosine of zenith angle                ------------
          if ( id_coszen > 0 ) then
             used = send_data ( id_coszen, zencos, Time, is, js )
          endif
!------- Net SW surface flux                   ------------
          if ( id_flux_sw > 0 ) then
             used = send_data ( id_flux_sw, sw_flux, Time, is, js )
          endif
!------- Net LW surface flux                   ------------
          if ( id_flux_lw > 0 ) then
             used = send_data ( id_flux_lw, lw_flux, Time, is, js )
          endif
!------- Net SW radiation at TOA               ------------
          if ( id_isr > 0 ) then
             used = send_data ( id_isr, isr, Time, is, js )
          endif
!------- Incoming SW radiation at TOA          ------------
          if ( id_swdn_toa > 0 ) then
             used = send_data ( id_swdn_toa, swdn_toa, Time, is, js )
          endif
!------- Upward LW flux at the surface         ------------
          if ( id_lwup_sfc > 0 ) then
             used = send_data ( id_lwup_sfc, lwup_sfc, Time, is, js )
          endif
!------- Outgoing LW radiation                   ------------
          if ( id_olr > 0 ) then
             used = send_data ( id_olr, olr, Time, is, js )
          endif
!------- Interactive albedo                    ------------
          if ( present(albedo_loc)) then
             used = send_data ( id_albedo, albedo_loc, Time, is, js )
          endif
!------- Ozone                                 ------------
          if ( present(ozone) .and. id_ozone > 0 ) then
             used = send_data ( id_ozone, ozone, Time, is, js, 1 )
          endif
!------- Half grid point temperature                                ------------
          if ( present(thalf) .and. id_thalf > 0 ) then
             used = send_data ( id_thalf, thalf, Time, is, js, 1 )
          endif
        end subroutine write_diag_rrtm
!*****************************************************************************************

        subroutine rrtm_precip_accum(precip, rain, snow)
!
! Count where it precipitates, for the precipitation-dependent albedo.
! precip is the total precipitation; rain and snow the large-scale parts.
!
          implicit none
          real(kind=rb),dimension(:,:),intent(in) :: precip, rain, snow

          if(do_precip_albedo)then
             if( trim(precip_albedo_mode) .eq. 'full' )then
                where(precip > 0.) rrtm_precip = rrtm_precip + 1.
             elseif( trim(precip_albedo_mode) .eq. 'lscale')then
                where(rain+snow > 0.) rrtm_precip = rrtm_precip + 1. !precip -> total precip, rain+snow -> lscale
             elseif( trim(precip_albedo_mode) .eq. 'conv') then
                where(precip-rain-snow > 0.) rrtm_precip = rrtm_precip + 1.
             endif
             num_precip = num_precip + 1
          endif
        end subroutine rrtm_precip_accum
!*****************************************************************************************

        subroutine rrtm_radiation_end
          use mima_interpolator_mod, only: interpolator_end
          implicit none

          if(do_read_ozone)call interpolator_end(o3_interp)

          call write_restart_rrtm

          if(allocated(t_half))      deallocate(t_half)
          if(allocated(h2o))         deallocate(h2o, o3, co2, ones, zeros, emis, &
                                                taucld, tauaer, sw_zro, zro_sw)
          if(allocated(zencos))      deallocate(zencos)
          if(allocated(sw_flux))     deallocate(sw_flux)
          if(allocated(lw_flux))     deallocate(lw_flux)
          if(allocated(rrtm_precip)) deallocate(rrtm_precip)
          if(allocated(tdt_rad))     deallocate(tdt_rad)
          if(allocated(tdt_sw_rad))  deallocate(tdt_sw_rad)
          if(allocated(tdt_lw_rad))  deallocate(tdt_lw_rad)
          if(allocated(isr))         deallocate(isr)
          if(allocated(olr))         deallocate(olr)
          if(allocated(swdn_toa))    deallocate(swdn_toa)
          if(allocated(lwup_sfc))    deallocate(lwup_sfc)
          rrtm_init = .false.

        end subroutine rrtm_radiation_end

!*****************************************************************************************
        subroutine read_restart_rrtm
!
! Read the state kept between radiation steps, if INPUT/rrtm_radiation.res.nc exists.
! If a field in use is missing (e.g. a diagnostic was added), recompute radiation at
! the first step, as without a restart file.
!
          use fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, NOTE
          implicit none
          type(restart_file_type) :: rst
          real(kind=rb) :: x
          logical :: found

          if(.not. open_restart_read(rst, 'INPUT/rrtm_radiation.res.nc', domain)) return
          if(mpp_pe() == mpp_root_pe()) call error_mesg(mod_name, &
               'Reading NetCDF formatted restart file: INPUT/rrtm_radiation.res.nc', NOTE)
          found = .true.
          if(allocated(tdt_rad))    call read_saved(rst, 'tdt_rad',    tdt_rad,    found)
          if(allocated(sw_flux))    call read_saved(rst, 'sw_flux',    sw_flux,    found)
          if(allocated(lw_flux))    call read_saved(rst, 'lw_flux',    lw_flux,    found)
          if(allocated(zencos))     call read_saved(rst, 'zencos',     zencos,     found)
          if(allocated(tdt_sw_rad)) call read_saved(rst, 'tdt_sw_rad', tdt_sw_rad, found)
          if(allocated(tdt_lw_rad)) call read_saved(rst, 'tdt_lw_rad', tdt_lw_rad, found)
          if(allocated(olr))        call read_saved(rst, 'olr',        olr,        found)
          if(allocated(isr))        call read_saved(rst, 'isr',        isr,        found)
          if(allocated(swdn_toa))   call read_saved(rst, 'swdn_toa',   swdn_toa,   found)
          if(allocated(lwup_sfc))   call read_saved(rst, 'lwup_sfc',   lwup_sfc,   found)
          if(found)then
             call read_restart_field(rst, 'dt_last', x)
             dt_last = nint(x)
          elseif(mpp_pe() == mpp_root_pe())then
             call error_mesg(mod_name, 'INPUT/rrtm_radiation.res.nc lacks fields now in use;'// &
                  ' radiation is recomputed at the first time step', NOTE)
          endif
          if(do_precip_albedo)then
             if(restart_field_exists(rst, 'rrtm_precip'))then
                call read_restart_field(rst, 'rrtm_precip', rrtm_precip)
                call read_restart_field(rst, 'num_precip', x)
                num_precip = nint(x)
             endif
          endif
          call close_restart(rst)

        end subroutine read_restart_rrtm

        subroutine read_saved_2d(rst, name, data, found)
          implicit none
          type(restart_file_type),intent(inout)    :: rst
          character(len=*),intent(in)              :: name
          real(kind=rb),dimension(:,:),intent(out) :: data
          logical,intent(inout)                    :: found

          if(restart_field_exists(rst, name))then
             call read_restart_field(rst, name, data)
          else
             found = .false.
          endif
        end subroutine read_saved_2d

        subroutine read_saved_3d(rst, name, data, found)
          implicit none
          type(restart_file_type),intent(inout)      :: rst
          character(len=*),intent(in)                :: name
          real(kind=rb),dimension(:,:,:),intent(out) :: data
          logical,intent(inout)                      :: found

          if(restart_field_exists(rst, name))then
             call read_restart_field(rst, name, data)
          else
             found = .false.
          endif
        end subroutine read_saved_3d
!*****************************************************************************************

        subroutine write_restart_rrtm
!
! Write the state kept between radiation steps (only the fields in use).
!
          use fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, NOTE
          implicit none
          type(restart_file_type) :: rst

          if(mpp_pe() == mpp_root_pe()) call error_mesg(mod_name, &
               'Writing NetCDF formatted restart file: RESTART/rrtm_radiation.res.nc', NOTE)
          call open_restart_write(rst, 'RESTART/rrtm_radiation.res.nc', domain)
          call write_restart_field(rst, 'dt_last', real(dt_last))
          if(allocated(tdt_rad))    call write_restart_field(rst, 'tdt_rad',    tdt_rad)
          if(allocated(sw_flux))    call write_restart_field(rst, 'sw_flux',    sw_flux)
          if(allocated(lw_flux))    call write_restart_field(rst, 'lw_flux',    lw_flux)
          if(allocated(zencos))     call write_restart_field(rst, 'zencos',     zencos)
          if(allocated(tdt_sw_rad)) call write_restart_field(rst, 'tdt_sw_rad', tdt_sw_rad)
          if(allocated(tdt_lw_rad)) call write_restart_field(rst, 'tdt_lw_rad', tdt_lw_rad)
          if(allocated(olr))        call write_restart_field(rst, 'olr',        olr)
          if(allocated(isr))        call write_restart_field(rst, 'isr',        isr)
          if(allocated(swdn_toa))   call write_restart_field(rst, 'swdn_toa',   swdn_toa)
          if(allocated(lwup_sfc))   call write_restart_field(rst, 'lwup_sfc',   lwup_sfc)
          if(do_precip_albedo)then
             call write_restart_field(rst, 'rrtm_precip', rrtm_precip)
             call write_restart_field(rst, 'num_precip',  real(num_precip))
          endif
          call close_restart(rst)

        end subroutine write_restart_rrtm
!*****************************************************************************************
      end module rrtm_radiation
