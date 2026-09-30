
                    module moist_processes_mod

!-----------------------------------------------------------------------
!
!         interface module for moisture processes
!         ---------------------------------------
!             Betts-Miller convective adjustment
!             moist convective adjustment
!             large-scale condensation
!
!-----------------------------------------------------------------------

use    betts_miller_mod, only: betts_miller, betts_miller_init

use      moist_conv_mod, only: moist_conv, moist_conv_init
use     lscale_cond_mod, only: lscale_cond, lscale_cond_init
use  sat_vapor_pres_mod, only: lookup_es

use    time_manager_mod, only: time_type

use    diag_manager_mod, only: register_diag_field, send_data

use             fms_mod, only: check_nml_error,    &
                               input_nml_file, &
                               write_version_number,           &
                               mpp_pe, mpp_root_pe, stdlog,    &
                               error_mesg, FATAL, NOTE

use mima_diag_integral_mod, only:     diag_integral_field_init, &
                             sum_diag_integral_field

use       constants_mod, only: CP_AIR, GRAV, RDGAS, RVGAS, HLV, KAPPA, ES0

use  field_manager_mod, only: MODEL_ATMOS
use tracer_manager_mod, only: get_tracer_index,&
                              get_number_tracers, &
                              get_tracer_names, &
                              query_method, &
                              NO_TRACER
use atmos_tracer_utilities_mod, only : wet_deposition
use            fms_mod, only : mpp_clock_id, mpp_clock_begin, &
                               mpp_clock_end, CLOCK_MODULE, &
                               MPP_CLOCK_SYNC

! start chemistry modules
!use    chem_interface, only : id_wet
!use           mo_hook, only : moz_hook
! end chemistry modules
!
!------------------mj pass precipitation to rrtm for albedo-------------
use radiation_mod, only: radiation_precip_accum
implicit none
private

!-----------------------------------------------------------------------
!-------------------- public data/interfaces ---------------------------

   public   moist_processes, moist_processes_init, moist_processes_end

!-----------------------------------------------------------------------
!-------------------- private data -------------------------------------


      real,parameter :: epst=200.

   real, parameter :: d622 = RDGAS/RVGAS
   real, parameter :: d378 = 1.-d622


!--------------------- version number ----------------------------------
   character(len=128) :: version = '$Id: moist_processes.f90,v 12.0.4.3 2005/05/21 02:02:03 pjp Exp $'
   character(len=128) :: tagname = '$Name:  $'
   logical            :: module_is_initialized = .false.
!-----------------------------------------------------------------------
!-------------------- namelist data (private) --------------------------


   logical :: do_mca=.false., do_lsc=.true.,  &
              use_tau=.false., do_gust_cv = .false., &
              do_bm=.true.

   real :: pdepth = 150.e2
   real :: tfreeze = 273.16
   real :: gustmax = 3.             ! maximum gustiness wind (m/s)
   real :: gustconst = 10./86400.   ! constant in kg/m2/sec, default =
                                    ! 1 cm/day = 10 mm/day

!---------------- namelist variable definitions ------------------------
!
!   do_mca   = switch to turn on/off moist convective adjustment;
!                [logical, default: do_mca=false ]
!   do_lsc   = switch to turn on/off large scale condensation
!                [logical, default: do_lsc=true ]
!   use_tau  = switch to determine whether current time level (tau)
!                will be used or else future time level (tau+1).
!                if use_tau = true then the input values for t,q, and r
!                are used; if use_tau = false then input values
!                tm+tdt*dt, etc. are used.
!                [logical, default: use_tau=false ]
!
!   pdepth   = boundary layer depth in pascals for determining mean
!                temperature tfreeze (used for snowfall determination)
!   tfreeze  = mean temperature used for snowfall determination (deg k)
!                [real, default: tfreeze=273.16]
!
!  do_gust_cv = switch to use convective gustiness (default = false)
!  gustmax    = maximum convective gustiness (m/s)
!  gustconst  = precip rate which defines precip rate which begins to
!               matter for convective gustiness (kg/m2/sec)
!
!   do_bm    = switch to turn on/off betts-miller scheme
!                [logical, default: do_bm=true ]
!
!   notes: pdepth and tfreeze are used to determine liquid vs. solid
!          precipitation for the mca and lsc schemes.
!
!-----------------------------------------------------------------------

namelist /moist_processes_nml/ do_mca, do_lsc,  &
                               pdepth, tfreeze,        &
                               use_tau, &
                               do_gust_cv, &
                               gustmax, gustconst, &
                               do_bm

!-----------------------------------------------------------------------
!-------------------- diagnostics fields -------------------------------

integer :: id_tdt_conv, id_qdt_conv, id_prec_conv, id_snow_conv, &
           id_tdt_ls  , id_qdt_ls  , id_prec_ls  , id_snow_ls  , &
           id_precip  , id_WVP, id_gust_conv, id_rh, &
           id_q_conv_col, id_q_ls_col, id_t_conv_col, id_t_ls_col, &
           id_cape, id_cin, id_tref, id_qref, id_rhsurf, &
           id_bmflag, id_klzbs, id_invtaubmt, id_invtaubmq, &
           id_entrop_ls

integer, dimension(:), allocatable :: id_tracerdt_conv,  &
                                      id_tracerdt_conv_col, &
                                      id_conv_tracer,  &
                                      id_conv_tracer_col
character(len=5) :: mod_name = 'moist'

real :: missing_value = -999.
integer :: convection_clock, largescale_clock

logical :: do_tracers_in_mca = .false.
logical, dimension(:), allocatable :: tracers_in_mca
integer :: num_mca_tracers=0
integer :: num_tracers=0

!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine moist_processes (is, ie, js, je, Time, dt, land,            &
                            phalf, pfull, zhalf, zfull, omega, diff_t, &
                            t, q, r, u, v, tm, qm, rm, um, vm,         &
                            tdt, qdt, rdt, udt, vdt,                   &
                            lprec, fprec, gust_cv, area,               &
                            lat, mask, kbot)

!-----------------------------------------------------------------------
!
!    in:  is,ie      starting and ending i indices for window
!
!         js,je      starting and ending j indices for window
!
!         Time       time used for diagnostics [time_type]
!
!         dt         time step (from t(n-1) to t(n+1) if leapfrog)
!                    in seconds   [real]
!
!         land        fraction of surface covered by land
!                     [real, dimension(nlon,nlat)]
!
!         phalf      pressure at half levels in pascals
!                      [real, dimension(nlon,nlat,nlev+1)]
!
!         pfull      pressure at full levels in pascals
!                      [real, dimension(nlon,nlat,nlev)]
!
!         omega      omega (vertical velocity) at full levels
!                    in pascals per second
!                      [real, dimension(nlon,nlat,nlev)]
!
!         diff_t     vertical diffusion coefficient for temperature
!                    and tracer (m*m/sec) on half levels
!                      [real, dimension(nlon,nlat,nlev)]
!
!         t, q       temperature (t) [deg k] and specific humidity
!                    of water vapor (q) [kg/kg] at full model levels,
!                    at the current time step if leapfrog scheme
!                      [real, dimension(nlon,nlat,nlev)]
!
!         r          tracer fields at full model levels,
!                    at the current time step if leapfrog
!                      [real, dimension(nlon,nlat,nlev,ntrace)]
!
!         u, v,      zonal and meridional wind [m/s] at full model levels,
!                    at the current time step if leapfrog scheme
!                      [real, dimension(nlon,nlat,nlev)]
!
!         tm, qm     temperature (t) [deg k] and specific humidity
!                    of water vapor (q) [kg/kg] at full model levels,
!                    at the previous time step if leapfrog scheme
!                      [real, dimension(nlon,nlat,nlev)]
!
!         rm         tracer fields at full model levels,
!                    at the previous time step if leapfrog
!                      [real, dimension(nlon,nlat,nlev,ntrace)]
!
!         um, vm     zonal and meridional wind [m/s] at full model levels,
!                    at the previous time step if leapfrog
!                      [real, dimension(nlon,nlat,nlev)]
!
!         area       grid box area (in m2)
!                      [real, dimension(nlon,nlat)]
!
!         lat        latitude in radians
!                      [real, dimension(nlon,nlat)]
!
! inout:  tdt, qdt   temperature (tdt) [deg k/sec] and specific
!                    humidity of water vapor (qdt) tendency [1/sec]
!                      [real, dimension(nlon,nlat,nlev)]
!
!         rdt        tracer tendencies
!                      [real, dimension(nlon,nlat,nlev,ntrace)]
!
!         udt, vdt   zonal and meridional wind tendencies [m/s/s]
!
!   out:  lprec      liquid precipitiaton rate (rain) in kg/m2/s
!                      [real, dimension(nlon,nlat)]
!
!         fprec      frozen precipitation rate (snow) in kg/m2/s
!                      [real, dimension(nlon,nlat)]
!
!         gust_cv    gustiness from convection  in m/s
!                      [real, dimension(nlon,nlat)]
!
!       optional
!  -----------------
!
!    in:  mask       mask (1. or 0.) for grid boxes above or below
!                    the ground   [real, dimension(nlon,nlat,nlev)]
!
!         kbot       index of the lowest model level
!                      [integer, dimension(nlon,nlat)]
!
!
!-----------------------------------------------------------------------
integer,         intent(in)              :: is,ie,js,je
type(time_type), intent(in)              :: Time
   real, intent(in)                      :: dt
   real, intent(in) , dimension(:,:)     :: land
   real, intent(in) , dimension(:,:,:)   :: phalf, pfull, zhalf, zfull,&
                                            omega, diff_t,             &
                                            t, q, u, v, tm, qm, um, vm
   real, intent(in) , dimension(:,:,:,:) :: r, rm
   real, intent(inout),dimension(:,:,:)  :: tdt, qdt, udt, vdt
   real, intent(inout),dimension(:,:,:,:):: rdt
   real, intent(out), dimension(:,:)     :: lprec, fprec, gust_cv
   real, intent(in) , dimension(:,:)     :: area
   real, intent(in) , dimension(:,:)     :: lat


   real, intent(in) , dimension(:,:,:), optional :: mask
integer, intent(in) , dimension(:,:),   optional :: kbot

!-----------------------------------------------------------------------
real, dimension(size(t,1),size(t,2),size(t,3)) :: tin,qin,ttnd,qtnd, &
                 rlin, riin, rnew, rlnew, rinew, qnew, qinew, qlnew, &
                 rltnd, ritnd, rin, rtnd, ahuco, qrat, entrop_ls
real, dimension(size(t,1),size(t,2),size(t,3)) :: ttnd_save,qtnd_save, &
                                                  qitnd_save,   &
                                                  qltnd_save, qatnd_save
real, dimension(size(t,1),size(t,2))           :: tsnow,snow
real, dimension(size(t,1),size(t,2))           :: snow_save, rain_save
logical,dimension(size(t,1),size(t,2))         :: coldT
real, dimension(size(t,1),size(t,2),size(t,3)) :: utnd,vtnd,uin,vin

real, dimension(size(t,1),size(t,2),size(t,3)) :: qltnd,qitnd,qatnd, &
                                                  qlin, qiin, qain
real, dimension(size(r,1),size(r,2),size(r,3),size(r,4)) :: tracer, tracertnd
real, dimension(size(t,1),size(t,2),size(t,3)) :: RH, pmass, wetdeptnd, q_ref, t_ref
real, dimension(size(t,1),size(t,2))           :: rain, precip, cape, cin
real, dimension(size(t,1),size(t,2))           :: wvp
real, dimension(size(t,1),size(t,2))           :: tempdiag, bmflag, &
                                                  klzbs, invtaubmt, invtaubmq
real, dimension(size(t,1),size(t,2),size(t,3)) :: tempdiag1
integer n

integer :: i, j, k, ix, jx, kx, nt, ip, tr
real    :: dtinv
logical :: use_mask, used

real, dimension(size(t,1),size(t,2),size(t,3)+1) :: press
real, dimension(size(t,1),size(t,2),size(t,3)) :: tin_in

real :: tinsave

! tracer code:
integer :: nn
real, dimension (size(t,1), size(t,2), &
                 size(t,3), num_mca_tracers) :: qtrmca , &
                                                mca_tracers
logical :: alpha, alphc

! end tracer code


real, dimension(size(rdt,1),size(rdt,2),size(rdt,3),size(rdt,4)) :: wet_data
!-----------------------------------------------------------------------

      if (.not. module_is_initialized) call error_mesg ('moist_processes',  &
                     'moist_processes_init has not been called.', FATAL)


!-------- input array size and position in global storage --------------

      ix=size(t,1); jx=size(t,2); kx=size(t,3); nt=size(rdt,4)

!-----------------------------------------------------------------------

      use_mask=.false.
      if (present(mask) .and. present(kbot))  use_mask=.true.

      dtinv=1./dt
      lprec=0.0; fprec=0.0; precip=0.0; rain=0.0; snow=0.0

!------------------ setup input data -----------------------------------

   if (use_tau) then
      tin (:,:,:)=t(:,:,:)
      qin (:,:,:)=q(:,:,:)
      uin (:,:,:)=u(:,:,:)
      vin (:,:,:)=v(:,:,:)
   else
      tin (:,:,:)=tm(:,:,:)+tdt(:,:,:)*dt
      qin (:,:,:)=qm(:,:,:)+qdt(:,:,:)*dt
      uin (:,:,:)=um(:,:,:)+udt(:,:,:)*dt
      vin (:,:,:)=vm(:,:,:)+vdt(:,:,:)*dt
   endif

   if (use_mask) then
      tin (:,:,:)=mask(:,:,:)*tin(:,:,:)+(1.0-mask(:,:,:))*epst
      qin (:,:,:)=mask(:,:,:)*qin(:,:,:)
      uin (:,:,:)=mask(:,:,:)*uin(:,:,:)
      vin (:,:,:)=mask(:,:,:)*vin(:,:,:)
   endif

! Make local copies of all the tracers
      if (use_tau) then
        do tr = 1, size(r,4)
          tracer (:,:,:,tr)=r(:,:,:,tr)
        enddo
      else
        do tr = 1, size(r,4)
          tracer (:,:,:,tr)=rm(:,:,:,tr)+rdt(:,:,:,tr)*dt
        enddo
      endif

      if (use_mask) then
        do tr = 1, size(r,4)
          tracer (:,:,:,tr)=mask(:,:,:)*tracer(:,:,:,tr)
        enddo
      endif


!   initialize qrat and ahuco

   qrat = 1.0
   ahuco = 0.0

!---compute mass in each layer if needed by any of the diagnostics -----

      alpha = any (id_tracerdt_conv_col > 0)
      alphc = any (id_conv_tracer_col > 0)

      if ( id_q_conv_col  > 0 .or. id_t_conv_col  > 0 .or. &
           id_q_ls_col    > 0 .or. id_t_ls_col    > 0 .or. &
           id_WVP         > 0 .or. alpha .or. alphc) then
        do k=1,kx
          pmass(:,:,k) = (phalf(:,:,k+1)-phalf(:,:,k))/GRAV
        end do
      end if

!----------------- mean temp in lower atmosphere -----------------------
!----------------- used for determination of rain vs. snow -------------
!----------------- use input temp ? ------------------------------------

   call tempavg (pdepth, phalf, t, tsnow, mask)

   where (tsnow <= tfreeze)
          coldT=.TRUE.
   elsewhere
          coldT=.FALSE.
   endwhere

call mpp_clock_begin( convection_clock )

!---------------------------------------------------------------------
!    initialize an array to hold tracer tendencies due to moist
!    processes.
!---------------------------------------------------------------------
      tracertnd = 0.0


      !------- diagnostics for tracers from convection -------
!  allow any tracer to be activated here (allows control cases)
      do n=1,num_tracers
        if ( id_conv_tracer(n) > 0 ) then
          used = send_data ( id_conv_tracer(n), tracer(:,:,:,n), Time, is, js, 1, &
                            rmask=mask )
         endif
!------- diagnostics for tracers column integral tendency ------
         if ( id_conv_tracer_col(n) > 0 ) then
           tempdiag(:,:)=0.
           do k=1,kx
             tempdiag(:,:) = tempdiag(:,:) + tracer   (:,:,k,n)*pmass(:,:,k)
           end do
           used = send_data ( id_conv_tracer_col(n), tempdiag, Time, is, js )
        end if
      enddo

!-----------------------------------------------------------------------
!***********************************************************************
!----------------- moist convective adjustment -------------------------

                        if (do_mca) then
!-----------------------------------------------------------------------

!---------------------------------------------------------------------
!    if any tracers are to be transported by moist_conv_mod,
!    check each active tracer to find those to be transported and fill
!    the mca_tracers array with these fields.
!---------------------------------------------------------------------
      if (num_mca_tracers > 0) then
        nn = 1
        do n=1, num_tracers
          if (tracers_in_mca(n)) then
            mca_tracers(:,:,:,nn) = tracer(:,:,:,n)
            nn = nn + 1
          endif
        end do
         call moist_conv (tin,qin,pfull,phalf,coldT,&
                          ttnd,qtnd,rain,snow,kbot,&
                          dtinv, Time, mask, is, js, &
                          tracers=mca_tracers, qtrmca=qtrmca )
      else
         call moist_conv (tin,qin,pfull,phalf,coldT,&
                          ttnd,qtnd,rain,snow,kbot,&
                          dtinv, Time, mask, is, js )
      endif

!------- add on tendency ----------
     tdt=tdt+ttnd; qdt=qdt+qtnd

!------- save total precip and snow ---------
      lprec=lprec+rain
      fprec=fprec+snow
      precip=precip+rain+snow

!-----------------------------------------------------------------------
                           endif ! if (do_mca)

  if (do_bm) then ! betts-miller cumulus param scheme

    call betts_miller (dt,tin,qin,pfull,phalf,coldT,rain,snow,ttnd,qtnd,&
                      q_ref,bmflag,klzbs,cape,cin,t_ref,invtaubmt,&
                      invtaubmq, mask=mask)

!------- (update input values and) compute tendency -----
    tin=tin+ttnd;    qin=qin+qtnd
    ttnd=ttnd*dtinv; qtnd=qtnd*dtinv
    rain=rain*dtinv; snow=snow*dtinv

!-------- add on tendency ----------
    tdt=tdt+ttnd; qdt=qdt+qtnd
!------- save total precip and snow ---------
    lprec=lprec+rain
    fprec=fprec+snow
    precip=precip+rain+snow

!-----------------------------------------------------------------------
  endif ! if (do_bm)

! Do wet deposition for the convective routine here
! qdt should give the vertical distribution here
      wet_data = 0.0
      do n = 1, size(tracertnd,4)
        wetdeptnd = 0.0
        call wet_deposition(n, t, pfull, phalf, rain, snow, qtnd, tracer(:,:,:,n), &
                            wetdeptnd, Time, 'convect', is, js, dt)
        tracertnd(:,:,:,n) = tracertnd(:,:,:,n) - wetdeptnd
!RSH9-11-03
! added here, previously effect not added to rdt
        rdt (:,:,:,n) = rdt(:,:,:,n) - wetdeptnd
        wet_data(:,:,:,n) = wetdeptnd
      enddo

!-----------------------------------------------------------------------
! do convective gustiness


  gust_cv = 0.0
  if (do_gust_cv) then
     where((rain+snow).gt.0.)
        gust_cv = gustmax * sqrt( (rain+snow)/(gustconst+(rain+snow)) )
     endwhere
  end if

!-----------------------------------------------------------------------
!***********************************************************************
!--------------- DIAGNOSTICS FOR CONVECTIVE SCHEME ---------------------
!-----------------------------------------------------------------------

 if (do_bm) then
     if ( id_tref > 0 ) then
       used = send_data ( id_tref, t_ref, Time, is, js, 1, &
                          rmask=mask )
     end if
     if ( id_qref > 0 ) then
       used = send_data ( id_qref, q_ref, Time, is, js, 1, &
                          rmask=mask )
     end if
     if (id_bmflag > 0 ) then
       used = send_data ( id_bmflag, bmflag, Time, is, js)
     end if
     if (id_klzbs > 0 ) then
       used = send_data (id_klzbs, klzbs, Time, is, js)
     end if
 endif

 if ( do_bm ) then
     if (id_invtaubmt > 0) then
       used = send_data (id_invtaubmt, invtaubmt, Time, is, js)
     end if
     if (id_invtaubmq > 0) then
       used = send_data (id_invtaubmq, invtaubmq, Time, is, js)
     end if
 end if

!------- diagnostics for dt/dt_conv -------
      if ( id_tdt_conv > 0 ) then
        used = send_data ( id_tdt_conv, ttnd, Time, is, js, 1, &
                           rmask=mask )
      endif
!------- diagnostics for dq/dt_conv -------
      if ( id_qdt_conv > 0 ) then
        used = send_data ( id_qdt_conv, qtnd, Time, is, js, 1, &
                           rmask=mask )
      endif
!------- diagnostics for precip_conv -------
      if ( id_prec_conv > 0 ) then
        used = send_data ( id_prec_conv, rain+snow, Time, is, js )
      endif
!------- diagnostics for snow_conv -------
      if ( id_snow_conv > 0 ) then
        used = send_data ( id_snow_conv, snow, Time, is, js )
      endif

!------- diagnostics for gust_conv -------
      if ( id_gust_conv > 0 ) then
        used = send_data ( id_gust_conv, gust_cv, Time, is, js )
      endif

!------- diagnostics for water vapor path tendency ----------
      if ( id_q_conv_col > 0 ) then
        tempdiag(:,:)=0.
        do k=1,kx
          tempdiag(:,:) = tempdiag(:,:) + qtnd(:,:,k)*pmass(:,:,k)
        end do
        used = send_data ( id_q_conv_col, tempdiag, Time, is, js )
      end if

!------- diagnostics for dry static energy tendency ---------
      if ( id_t_conv_col > 0 ) then
        tempdiag(:,:)=0.
        do k=1,kx
          tempdiag(:,:) = tempdiag(:,:) + ttnd(:,:,k)*CP_AIR*pmass(:,:,k)
        end do
        used = send_data ( id_t_conv_col, tempdiag, Time, is, js )
      end if

!------- diagnostics for tracers from convection -------
      do n = 1, size(tracertnd,4)
        if (tracers_in_mca(n)) then
          if ( id_tracerdt_conv(n) > 0 ) then
            used = send_data ( id_tracerdt_conv(n), tracertnd(:,:,:,n), Time, is, js, 1, &
                               rmask=mask )
          endif

!------- diagnostics for tracers column integral tendency ------
          if ( id_tracerdt_conv_col(n) > 0 ) then
            tempdiag(:,:)=0.
            do k=1,kx
              tempdiag(:,:) = tempdiag(:,:) + tracertnd(:,:,k,n)*pmass(:,:,k)
            end do
            used = send_data ( id_tracerdt_conv_col(n), tempdiag, Time, is, js )
          end if
        end if
      enddo


! convection diagnostics
  call mpp_clock_end ( convection_clock )
!-----------------------------------------------------------------------
  call mpp_clock_begin ( largescale_clock )

!-----------------------------------------------------------------------
!***********************************************************************
!----------------- large-scale condensation ----------------------------

                      if (do_lsc) then
!-----------------------------------------------------------------------

   call lscale_cond (tin,qin,pfull,phalf,coldT,rain,snow,ttnd,qtnd,&
                     mask=mask)

!------- (update input values and) compute tendency -----
      tin=tin+ttnd;    qin=qin+qtnd

      ttnd=ttnd*dtinv; qtnd=qtnd*dtinv
      rain=rain*dtinv; snow=snow*dtinv

!------- add on tendency ----------
     tdt=tdt+ttnd; qdt=qdt+qtnd

!------- save total precip and snow ---------
      lprec=lprec+rain
      fprec=fprec+snow
      precip=precip+rain+snow

!-----------------------------------------------------------------------
                           endif
! Do wet deposition for the large scale routine here
! qdt should give the vertical distribution here
do n = 1, size(tracertnd,4)
wetdeptnd = 0.0
call wet_deposition(n, t, pfull, phalf, rain, snow, qtnd, tracer(:,:,:,n), &
                    wetdeptnd, Time, 'lscale', is, js, dt)
tracertnd(:,:,:,n) = tracertnd(:,:,:,n) - wetdeptnd
!RSH9-11-03
! added here, previously effect not added to rdt
          rdt (:,:,:,n) = rdt(:,:,:,n) - wetdeptnd
wet_data(:,:,:,n) = wet_data(:,:,:,n) + wetdeptnd

! start chemistry
!if (id_wet(n) /= 0 ) then
!! --------send wet deposition data to diag ----------
!used = send_data(id_wet(n), wet_data(:,:,:,n), Time,is_in=is,js_in=js)
!endif
! end chemistry
enddo

!-----------------------------------------------------------------------
!***********************************************************************
!--------------- DIAGNOSTICS FOR LARGE-SCALE SCHEME --------------------
!-----------------------------------------------------------------------

 if ( do_lsc ) then
   if ( id_entrop_ls > 0 ) then
     do k=1,kx
        entrop_ls(:,:,k) = ttnd(:,:,k) / t(:,:,k) * phalf(:,:,kx) /1.e5
     end do
     used = send_data ( id_entrop_ls, entrop_ls, Time, is, js, 1, &
                        rmask=mask )
   endif
 endif

 if ( do_lsc ) then
!------- diagnostics for dt/dt_strat -------
      if ( id_tdt_ls > 0 ) then
        used = send_data ( id_tdt_ls, ttnd, Time, is, js, 1, &
                           rmask=mask )
      endif
!------- diagnostics for dq/dt_strat -------
      if ( id_qdt_ls > 0 ) then
        used = send_data ( id_qdt_ls, qtnd, Time, is, js, 1, &
                           rmask=mask )
      endif
!------- diagnostics for precip_strat -------
      if ( id_prec_ls > 0 ) then
        used = send_data ( id_prec_ls, rain+snow, Time, is, js )
      endif
!------- diagnostics for snow_strat -------
      if ( id_snow_ls > 0 ) then
        used = send_data ( id_snow_ls, snow, Time, is, js )
      endif

!------- diagnostics for water vapor path tendency ----------
      if ( id_q_ls_col > 0 ) then
        tempdiag(:,:)=0.
        do k=1,kx
          tempdiag(:,:) = tempdiag(:,:) + qtnd(:,:,k)*pmass(:,:,k)
        end do
        used = send_data ( id_q_ls_col, tempdiag, Time, is, js )
      end if

!------- diagnostics for dry static energy tendency ---------
      if ( id_t_ls_col > 0 ) then
        tempdiag(:,:)=0.
        do k=1,kx
          tempdiag(:,:) = tempdiag(:,:) + ttnd(:,:,k)*CP_AIR*pmass(:,:,k)
        end do
        used = send_data ( id_t_ls_col, tempdiag, Time, is, js )
      end if

 endif  !end large-scale or strat diagnostics
  call mpp_clock_end ( largescale_clock )

!-----------------------------------------------------------------------
!------------------mj pass precipitation to rrtm for albedo-------------
  call radiation_precip_accum(precip, rain, snow)
!-----------------------------------------------------------------------
!***********************************************************************
!--------------------- GENERAL DIAGNOSTICS -----------------------------
!-----------------------------------------------------------------------

!------- diagnostics for total precip -------
   if ( id_precip > 0 ) then
        used = send_data ( id_precip, precip, Time, is, js )
   endif


!-----------------------------------------------------------------------
!------- diagnostics for column water vapor, liquid water path and
!------- ice water path


!-- compute and write out water vapor path --
   if ( id_WVP > 0 ) then
        wvp(:,:)=0.
        do k=1,kx
          wvp(:,:) = wvp(:,:) + qin(:,:,k)*pmass(:,:,k)
        end do

        used = send_data ( id_WVP, wvp, Time, is, js )
   end if




!-----------------------------------------------------------------------
!------- diagnostics for relative humidity -------

   if ( id_rh > 0 .or. id_rhsurf > 0 ) then
      call rh_calc (pfull,tin,qin,RH,mask)
      if ( id_rh     > 0 ) used = send_data ( id_rh, RH*100., Time, is, js, 1, rmask=mask )
      if ( id_rhsurf > 0 ) used = send_data ( id_rhsurf, RH(:,:,kx)*100., Time, is, js )
   endif

!-----------------------------------------------------------------------
!------- diagnostics for CAPE and CIN

!!-- compute and write out CAPE and CIN--
   if ( id_cape > 0 .or. id_cin > 0) then
!! calculate r
         rin = qin/(1.0 - qin)
         do j=js,je
            do i=is,ie
               call capecalcnew( kx, pfull(i,j,:), phalf(i,j,:), CP_AIR, RDGAS, RVGAS, &
                         HLV, KAPPA, tin(i,j,:), rin(i,j,:), cape(i,j), cin(i,j))
            end do
         end do
        if (id_cape > 0) then
             used = send_data ( id_cape, cape, Time, is, js )
        end if

        if ( id_cin > 0 ) then
             used = send_data ( id_cin, cin, Time, is, js )
        end if
   end if

!-----------------------------------------------------------------------
!---- accumulate global integral of precipiation (mm/day) -----

call sum_diag_integral_field ('prec', precip*86400., is, js)

!-----------------------------------------------------------------------
!    print *, ' end moist_processes             ', mpp_pe()

end subroutine moist_processes

!#######################################################################

subroutine moist_processes_init ( id, jd, kd, lonb, latb, pref, &
!                                 axes, Time, doing_strat)
                                  axes, Time)

!-----------------------------------------------------------------------
integer,         intent(in) :: id, jd, kd, axes(4)
real, dimension(:), intent(in) :: lonb, latb, pref
type(time_type), intent(in) :: Time
!logical,         intent(out) :: doing_strat
!-----------------------------------------------------------------------
!
!      input
!     --------
!
!      id, jd        number of horizontal grid points in the global
!                    fields along the x and y axis, repectively.
!                      [integer]
!
!      kd            number of vertical points in a column of atmosphere
!-----------------------------------------------------------------------

integer :: unit,io,ierr, n, nt, ntprog
character(len=32) :: tracer_units, tracer_name
character(len=80)  :: scheme
!-----------------------------------------------------------------------

       if ( module_is_initialized ) return

       read (input_nml_file, nml=moist_processes_nml, iostat=io)
       ierr = check_nml_error(io,'moist_processes_nml')

!--------- write version and namelist to standard log ------------

      call write_version_number ( version, tagname )
      if ( mpp_pe() == mpp_root_pe() ) &
      write ( stdlog(), nml=moist_processes_nml )

!------------------- dummy checks --------------------------------------

         if ( do_mca .and. do_bm ) call error_mesg   &
                   ('moist_processes_init',  &
                    'both do_mca and do_bm cannot be specified', FATAL)


!------------ initialize various schemes ----------
      if (do_bm) call betts_miller_init ()
      if (do_lsc)    call lscale_cond_init ()


!----- initialize quantities for global integral package -----

   call diag_integral_field_init ('prec', 'f6.3')


!----- initialize clocks -----

   convection_clock =     &
       mpp_clock_id( '   Physics_up: Moist Proc: Conv', &
           grain=CLOCK_MODULE, flags = MPP_CLOCK_SYNC )
   largescale_clock =     &
       mpp_clock_id( '   Physics_up: Moist Proc: LS', &
           grain=CLOCK_MODULE, flags = MPP_CLOCK_SYNC )

!  doing_strat = do_strat

!---------------------------------------------------------------------
!    retrieve the number of registered tracers in order to determine
!    which tracers are to be convectively transported.
!---------------------------------------------------------------------
      call get_number_tracers (MODEL_ATMOS, num_tracers= num_tracers)

!---------------------------------------------------------------------
!    allocate logical arrays to indicate the tracers which are to be
!    transported by the various available convective schemes.
!    initialize these arrays to .false..
!---------------------------------------------------------------------
      allocate (tracers_in_mca(num_tracers))
      tracers_in_mca = .false.

!----------------------------------------------------------------------
!    for each tracer, determine if it is to be transported by convect-
!    ion, and the convection schemes that are to transport it. set a
!    logical flag to .true. for each tracer that is to be transported by
!    each scheme and increment the count of tracers to be transported
!    by that scheme.
!----------------------------------------------------------------------
      do n=1, num_tracers
        if (query_method ('convection', MODEL_ATMOS, n, scheme)) then
          select case (scheme)
            case ("none")
            case ("mca")
               num_mca_tracers = num_mca_tracers + 1
               tracers_in_mca(n) = .true.
            case ("donner_and_mca", "mca_and_ras", "all")
               num_mca_tracers = num_mca_tracers + 1
               tracers_in_mca(n) = .true.
            case default  ! corresponds to "none"
          end select
        endif
      end do

!--------------------------------------------------------------------
!    set a logical indicating if any tracers are to be transported by
!    each of the available convection parameterizations.
!--------------------------------------------------------------------
      if (num_mca_tracers > 0) then
        do_tracers_in_mca = .true.
      else
        do_tracers_in_mca = .false.
      endif

!--------------------------------------------------------------------
!    initialize the convection scheme modules.
!--------------------------------------------------------------------
      if (do_mca)  then
        call  moist_conv_init (axes,Time, tracers_in_mca)
      endif


!----- initialize quantities for diagnostics output -----
      call diag_field_init ( axes, Time )


       module_is_initialized = .true.

!-----------------------------------------------------------------------

end subroutine moist_processes_init

!#######################################################################

subroutine moist_processes_end

      if( .not.module_is_initialized ) return

!----------------close various schemes-----------------

      module_is_initialized = .false.

!-----------------------------------------------------------------------

end subroutine moist_processes_end

!#######################################################################
!#######################################################################

      subroutine tempavg (pdepth,phalf,temp,tsnow,mask)

!-----------------------------------------------------------------------
!
!    computes a mean atmospheric temperature for the bottom
!    "pdepth" pascals of the atmosphere.
!
!   input:  pdepth     atmospheric layer in pa.
!           phalf      pressure at model layer interfaces
!           temp       temperature at model layers
!           mask       data mask at model layers (0.0 or 1.0)
!
!   output:  tsnow     mean model temperature in the lowest
!                      "pdepth" pascals of the atmosphere
!
!-----------------------------------------------------------------------
      real, intent(in)  :: pdepth
      real, intent(in) , dimension(:,:,:) :: phalf,temp
      real, intent(out), dimension(:,:)   :: tsnow
      real, intent(in) , dimension(:,:,:), optional :: mask
!-----------------------------------------------------------------------
 real, dimension(size(temp,1),size(temp,2)) :: prsum, done, pdel, pdep
 real  sumdone
 integer  k
!-----------------------------------------------------------------------

      tsnow=0.0; prsum=0.0; done=1.0; pdep=pdepth

      do k=size(temp,3),1,-1

         if (present(mask)) then
           pdel(:,:)=(phalf(:,:,k+1)-phalf(:,:,k))*mask(:,:,k)*done(:,:)
         else
           pdel(:,:)=(phalf(:,:,k+1)-phalf(:,:,k))*done(:,:)
         endif

         where ((prsum(:,:)+pdel(:,:))  >  pdep(:,:))
            pdel(:,:)=pdepth-prsum(:,:)
            done(:,:)=0.0
            pdep(:,:)=0.0
         endwhere

         tsnow(:,:)=tsnow(:,:)+pdel(:,:)*temp(:,:,k)
         prsum(:,:)=prsum(:,:)+pdel(:,:)

         sumdone=sum(done(:,:))
         if (sumdone < 1.e-4) exit

      enddo

         tsnow(:,:)=tsnow(:,:)/prsum(:,:)

!-----------------------------------------------------------------------

      end subroutine tempavg

!#######################################################################

      subroutine rh_calc(pfull,T,qv,RH,MASK)

        IMPLICIT NONE


        REAL, INTENT (IN),    DIMENSION(:,:,:) :: pfull,T,qv
        REAL, INTENT (OUT),   DIMENSION(:,:,:) :: RH
        REAL, INTENT (IN), OPTIONAL, DIMENSION(:,:,:) :: MASK

        REAL, DIMENSION(SIZE(T,1),SIZE(T,2),SIZE(T,3)) :: esat

!-----------------------------------------------------------------------
!       Calculate RELATIVE humidity.
!       This is calculated according to the formula:
!
!       RH   = qv / (epsilon*esat/ [pfull  -  (1.-epsilon)*esat])
!
!       Where epsilon = RDGAS/RVGAS = d622
!
!       and where 1- epsilon = d378
!
!       Note that RH does not have its proper value
!       until all of the following code has been executed.  That
!       is, RH is used to store intermediary results
!       in forming the full solution.

        !calculate water saturated vapor pressure from table
        !and store temporarily in the variable esat
        CALL LOOKUP_ES(T,esat)

        !calculate denominator in qsat formula
        RH(:,:,:) = pfull(:,:,:)

        !limit denominator to esat, and thus qs to epsilon
        !this is done to avoid blow up in the upper stratosphere
        !where pfull ~ esat
        RH(:,:,:) = MAX(RH(:,:,:),esat(:,:,:))

        !calculate RH
        RH(:,:,:)=qv(:,:,:)/(d622*esat(:,:,:)/RH(:,:,:))

        !IF MASK is present set RH to zero
        IF (present(MASK)) RH(:,:,:)=MASK(:,:,:)*RH(:,:,:)


END SUBROUTINE rh_calc

!#######################################################################

!all new cape calculation.

subroutine capecalcnew(kx,p,phalf,cp,rdgas,rvgas,hlv,kappa,tin,rin,&
                                cape,cin)

!
!    Input:
!
!    kx          number of levels
!    p           pressure (index 1 refers to TOA, index kx refers to surface)
!    phalf       pressure at half levels
!    cp          specific heat of dry air
!    rdgas       gas constant for dry air
!    rvgas       gas constant for water vapor (used in Clausius-Clapeyron,
!                not for virtual temperature effects, which are not considered)
!    hlv         latent heat of vaporization
!    kappa       the constant kappa
!    tin         temperature of the environment
!    rin         specific humidity of the environment
!
!    Output:
!    cape        Convective available potential energy
!    cin         Convective inhibition (if there's no LFC, then this is set
!                to zero)
!
!    Algorithm:
!    Start with surface parcel.
!    Calculate the lifting condensation level (uses an analytic formula and a
!       lookup table).
!    Average under the LCL if desired, if this is done, then a new LCL must
!       be calculated.
!    Calculate parcel ascent up to LZB.
!    Calculate CAPE and CIN.
      implicit none
      integer, intent(in)                    :: kx
      real, intent(in), dimension(:)         :: p, phalf, tin, rin
      real, intent(in)                       :: rdgas, rvgas, hlv, kappa, cp
      real, intent(out)                      :: cape, cin

      integer            :: k, klcl, klfc, klzb, klcl2
      logical            :: nocape
      real, dimension(kx)   :: theta, tp, rp
      real                  :: t0, r0, es, rs, theta0, pstar, value, tlcl, &
                               a, b, dtdlnp, d2tdlnp2, thetam, rm, tlcl2, &
                               plcl2, plcl, plzb

      pstar = 1.e5

      nocape = .true.
      cape = 0.
      cin = 0.
      plcl = 0.
      plzb = 0.
      klfc = 0
      klcl = 0
      klzb = 0
      tp(1:kx) = tin(1:kx)
      rp(1:kx) = rin(1:kx)

! start with surface parcel
      t0 = tin(kx)
      r0 = rin(kx)
! calculate the lifting condensation level by the following:
! are you saturated to begin with?
      call lookup_es(t0,es)
      rs = rdgas/rvgas*es/p(kx)
      if (r0.ge.rs) then
! if you're already saturated, set lcl to be the surface value.
         plcl = p(kx)
! the first level where you're completely saturated.
         klcl = kx
! saturate out to get the parcel temp and humidity at this level
! first order (in delta T) accurate expression for change in temp
         tp(kx) = t0 + (r0 - rs)/(cp/hlv + hlv*rs/rvgas/t0**2.)
         call lookup_es(tp(kx),es)
         rp(kx) = rdgas/rvgas*es/p(kx)
      else
! if not saturated to begin with, use the analytic expression to calculate the
! exact pressure and temperature where you?re saturated.
         theta0 = tin(kx)*(pstar/p(kx))**kappa
! the expression that we utilize is
! log(r/theta**(1/kappa)*pstar*rvgas/rdgas/es00) = log(es/T**(1/kappa))
! The right hand side of this is only a function of temperature, therefore
! this is put into a lookup table to solve for temperature.
         if (r0.gt.0.) then
            value = log(theta0**(-1/kappa)*r0*pstar*rvgas/rdgas/es0)
            call lcltabl(value,tlcl)
            plcl = pstar*(tlcl/theta0)**(1/kappa)
! just in case plcl is very high up
            if (plcl.lt.p(1)) then
               plcl = p(1)
               tlcl = theta0*(plcl/pstar)**kappa
               write (*,*) 'hi lcl'
            end if
            k = kx
         else
! if the parcel sp hum is zero or negative, set lcl to 2nd to top level
            plcl = p(2)
            tlcl = theta0*(plcl/pstar)**kappa
!            write (*,*) 'zero r0', r0
            do k=2,kx
               tp(k) = theta0*(p(k)/pstar)**kappa
               rp(k) = 0.
! this definition of CIN contains everything below the LCL
               cin = cin + rdgas*(tin(k) - tp(k))*log(phalf(k+1)/phalf(k))
            end do
            go to 11
         end if
! calculate the parcel temperature (adiabatic ascent) below the LCL.
! the mixing ratio stays the same
         do while (p(k).gt.plcl)
            tp(k) = theta0*(p(k)/pstar)**kappa
            call lookup_es(tp(k),es)
            rp(k) = rdgas/rvgas*es/p(k)
! this definition of CIN contains everything below the LCL
            cin = cin + rdgas*(tin(k) - tp(k))*log(phalf(k+1)/phalf(k))
            k = k-1
         end do
! first level where you're saturated at the level
         klcl = k
	 if (klcl.eq.1) klcl = 2
! do a saturated ascent to get the parcel temp at the LCL.
! use your 2nd order equation up to the pressure above.
! moist adaibat derivatives: (use the lcl values for temp, humid, and
! pressure)
         a = kappa*tlcl + hlv/cp*r0
         b = hlv**2.*r0/cp/rvgas/tlcl**2.
         dtdlnp = a/(1. + b)
! first order in p
!         tp(klcl) = tlcl + dtdlnp*log(p(klcl)/plcl)
! second order in p (RK2)
! first get temp halfway up
         tp(klcl) = tlcl + dtdlnp*log(p(klcl)/plcl)/2.
         if ((tp(klcl).lt.173.16).and.nocape) go to 11
         call lookup_es(tp(klcl),es)
         rp(klcl) = rdgas/rvgas*es/(p(klcl) + plcl)*2.
         a = kappa*tp(klcl) + hlv/cp*rp(klcl)
         b = hlv**2./cp/rvgas*rp(klcl)/tp(klcl)**2.
         dtdlnp = a/(1. + b)
! second half of RK2
         tp(klcl) = tlcl + dtdlnp*log(p(klcl)/plcl)
!         d2tdlnp2 = (kappa + b - 1. - b/tlcl*(hlv/rvgas/tlcl - &
!                   2.)*dtdlnp)/ (1. + b)*dtdlnp - hlv*r0/cp/ &
!                   (1. + b)
! second order in p
!         tp(klcl) = tlcl + dtdlnp*log(p(klcl)/plcl) + .5*d2tdlnp2*(log(&
!             p(klcl)/plcl))**2.
         if ((tp(klcl).lt.173.16).and.nocape) go to 11
         call lookup_es(tp(klcl),es)
         rp(klcl) = rdgas/rvgas*es/p(klcl)
!         write (*,*) 'tp, rp klcl:kx, new', tp(klcl:kx), rp(klcl:kx)
! CAPE/CIN stuff
         if ((tp(klcl).lt.tin(klcl)).and.nocape) then
! if you're not yet buoyant, then add to the CIN and continue
            cin = cin + rdgas*(tin(klcl) - &
                 tp(klcl))*log(phalf(klcl+1)/phalf(klcl))
         else
! if you're buoyant, then add to cape
            cape = cape + rdgas*(tp(klcl) - &
                  tin(klcl))*log(phalf(klcl+1)/phalf(klcl))
! if it's the first time buoyant, then set the level of free convection to k
            if (nocape) then
               nocape = .false.
               klfc = klcl
            endif
         end if
      end if
! then, start at the LCL, and do moist adiabatic ascent by the first order
! scheme -- 2nd order as well
      do k=klcl-1,1,-1
         a = kappa*tp(k+1) + hlv/cp*rp(k+1)
         b = hlv**2./cp/rvgas*rp(k+1)/tp(k+1)**2.
         dtdlnp = a/(1. + b)
! first order in p
!         tp(k) = tp(k+1) + dtdlnp*log(p(k)/p(k+1))
! second order in p (RK2)
! first get temp halfway up
         tp(k) = tp(k+1) + dtdlnp*log(p(k)/p(k+1))/2.
         if ((tp(k).lt.173.16).and.nocape) go to 11
         call lookup_es(tp(k),es)
         rp(k) = rdgas/rvgas*es/(p(k) + p(k+1))*2.
         a = kappa*tp(k) + hlv/cp*rp(k)
         b = hlv**2./cp/rvgas*rp(k)/tp(k)**2.
         dtdlnp = a/(1. + b)
! second half of RK2
         tp(k) = tp(k+1) + dtdlnp*log(p(k)/p(k+1))
!         d2tdlnp2 = (kappa + b - 1. - b/tp(k+1)*(hlv/rvgas/tp(k+1) - &
!               2.)*dtdlnp)/(1. + b)*dtdlnp - hlv/cp*rp(k+1)/(1. + b)
! second order in p

!         tp(k) = tp(k+1) + dtdlnp*log(p(k)/p(k+1)) + .5*d2tdlnp2*(log( &
!             p(k)/p(k+1)))**2.
! if you're below the lookup table value, just presume that there's no way
! you could have cape and call it quits
         if ((tp(k).lt.173.16).and.nocape) go to 11
         call lookup_es(tp(k),es)
         rp(k) = rdgas/rvgas*es/p(k)
         if ((tp(k).lt.tin(k)).and.nocape) then
! if you're not yet buoyant, then add to the CIN and continue
            cin = cin + rdgas*(tin(k) - tp(k))*log(phalf(k+1)/phalf(k))
         elseif((tp(k).lt.tin(k)).and.(.not.nocape)) then
! if you have CAPE, and it's your first time being negatively buoyant,
! then set the level of zero buoyancy to k+1, and stop the moist ascent
            klzb = k+1
            go to 11
         else
! if you're buoyant, then add to cape
            cape = cape + rdgas*(tp(k) - tin(k))*log(phalf(k+1)/phalf(k))
! if it's the first time buoyant, then set the level of free convection to k
            if (nocape) then
               nocape = .false.
               klfc = k
            endif
         end if
      end do
 11   if(nocape) then
! this is if you made it through without having a LZB
! set LZB to be the top level.
         plzb = p(1)
         klzb = 0
         klfc = 0
         cin = 0.
         tp(1:kx) = tin(1:kx)
         rp(1:kx) = rin(1:kx)
      end if
! if the parcel is still buoyant at the top level, take that as the LZB
      if (.not. nocape .and. klzb == 0) klzb = 1
!      write (*,*) 'plcl, klcl, tlcl, r0 new', plcl, klcl, tlcl, r0
!      write (*,*) 'tp, rp new', tp, rp
!       write (*,*) 'tp, new', tp
!       write (*,*) 'tin new', tin
!       write (*,*) 'klcl, klfc, klzb new', klcl, klfc, klzb
      end subroutine capecalcnew

!#######################################################################

! lookup table for the analytic evaluation of LCL
      subroutine lcltabl(value,tlcl)
!
! Table of values used to compute the temperature of the lifting condensation
! level.
!
! the expression that we utilize is
! log(r/theta**(1/kappa)*pstar*rvgas/rdgas/es00) = log(es/T**(1/kappa))
! the RHS is tabulated for the control amount of moisture, hence the
! division by es00 on the LHS

! Gives the values of the temperature for the following range:
!   starts with -23, is uniformly distributed up to -10.4.  There are a
! total of 127 values, and the increment is .1.
!
      implicit none
      real, intent(in)     :: value
      real, intent(out)    :: tlcl

      integer              :: ival
      real, dimension(127) :: lcltable
      real                 :: v1, v2

      data lcltable/   1.7364512e+02,   1.7427449e+02,   1.7490874e+02, &
      1.7554791e+02,   1.7619208e+02,   1.7684130e+02,   1.7749563e+02, &
      1.7815514e+02,   1.7881989e+02,   1.7948995e+02,   1.8016539e+02, &
      1.8084626e+02,   1.8153265e+02,   1.8222461e+02,   1.8292223e+02, &
      1.8362557e+02,   1.8433471e+02,   1.8504972e+02,   1.8577068e+02, &
      1.8649767e+02,   1.8723077e+02,   1.8797006e+02,   1.8871561e+02, &
      1.8946752e+02,   1.9022587e+02,   1.9099074e+02,   1.9176222e+02, &
      1.9254042e+02,   1.9332540e+02,   1.9411728e+02,   1.9491614e+02, &
      1.9572209e+02,   1.9653521e+02,   1.9735562e+02,   1.9818341e+02, &
      1.9901870e+02,   1.9986158e+02,   2.0071216e+02,   2.0157057e+02, &
      2.0243690e+02,   2.0331128e+02,   2.0419383e+02,   2.0508466e+02, &
      2.0598391e+02,   2.0689168e+02,   2.0780812e+02,   2.0873335e+02, &
      2.0966751e+02,   2.1061074e+02,   2.1156316e+02,   2.1252493e+02, &
      2.1349619e+02,   2.1447709e+02,   2.1546778e+02,   2.1646842e+02, &
      2.1747916e+02,   2.1850016e+02,   2.1953160e+02,   2.2057364e+02, &
      2.2162645e+02,   2.2269022e+02,   2.2376511e+02,   2.2485133e+02, &
      2.2594905e+02,   2.2705847e+02,   2.2817979e+02,   2.2931322e+02, &
      2.3045895e+02,   2.3161721e+02,   2.3278821e+02,   2.3397218e+02, &
      2.3516935e+02,   2.3637994e+02,   2.3760420e+02,   2.3884238e+02, &
      2.4009473e+02,   2.4136150e+02,   2.4264297e+02,   2.4393941e+02, &
      2.4525110e+02,   2.4657831e+02,   2.4792136e+02,   2.4928053e+02, &
      2.5065615e+02,   2.5204853e+02,   2.5345799e+02,   2.5488487e+02, &
      2.5632953e+02,   2.5779231e+02,   2.5927358e+02,   2.6077372e+02, &
      2.6229310e+02,   2.6383214e+02,   2.6539124e+02,   2.6697081e+02, &
      2.6857130e+02,   2.7019315e+02,   2.7183682e+02,   2.7350278e+02, &
      2.7519152e+02,   2.7690354e+02,   2.7863937e+02,   2.8039954e+02, &
      2.8218459e+02,   2.8399511e+02,   2.8583167e+02,   2.8769489e+02, &
      2.8958539e+02,   2.9150383e+02,   2.9345086e+02,   2.9542719e+02, &
      2.9743353e+02,   2.9947061e+02,   3.0153922e+02,   3.0364014e+02, &
      3.0577420e+02,   3.0794224e+02,   3.1014515e+02,   3.1238386e+02, &
      3.1465930e+02,   3.1697246e+02,   3.1932437e+02,   3.2171609e+02, &
      3.2414873e+02,   3.2662343e+02,   3.2914139e+02,   3.3170385e+02 /

      v1 = value
      if (value.lt.-23.0) v1 = -23.0
      if (value.gt.-10.4) v1 = -10.4
      ival = floor(10.*(v1 + 23.0))
      v2 = -230. + ival
      v1 = 10.*v1
! the table has 127 entries; at the upper clamp ival+2 would be 128
      tlcl = (v2 + 1.0 - v1)*lcltable(ival+1) + (v1 - v2)*lcltable(min(ival+2,127))

      end subroutine lcltabl

!#######################################################################

subroutine diag_field_init ( axes, Time )

  integer,         intent(in) :: axes(4)
  type(time_type), intent(in) :: Time

  character(len=32) :: tracer_units, tracer_name
  character(len=128) :: diaglname
  integer, dimension(3) :: half = (/1,2,4/)
  integer   :: n, nn

!------------ initializes diagnostic fields in this module -------------

if (do_bm) then
   id_qref = register_diag_field ( mod_name, &
     'qref', axes(1:3), Time, &
     'Adjustment reference specific humidity profile', &
     'kg/kg',  missing_value=missing_value               )
   id_tref = register_diag_field ( mod_name, &
     'tref', axes(1:3), Time, &
     'Adjustment reference temperature profile', &
     'K',  missing_value=missing_value                   )
   id_bmflag = register_diag_field (mod_name, &
      'bmflag', axes(1:2), Time, &
      'Betts-Miller flag', &
      'no units', missing_value=missing_value            )
   id_klzbs  = register_diag_field  (mod_name, &
      'klzbs', axes(1:2), Time, &
      'klzb', &
      'no units', missing_value=missing_value            )
endif

if ( do_bm ) then
   id_invtaubmt  = register_diag_field  (mod_name, &
      'invtaubmt', axes(1:2), Time, &
      'Inverse temperature relaxation time', &
      '1/s', missing_value=missing_value            )
   id_invtaubmq = register_diag_field  (mod_name, &
      'invtaubmq', axes(1:2), Time, &
      'Inverse humidity relaxation time', &
      '1/s', missing_value=missing_value            )
end if  ! if ( do_bm )

   id_tdt_conv = register_diag_field ( mod_name, &
     'tdt_conv', axes(1:3), Time, &
     'Temperature tendency',    'deg_K/s',  &
                        missing_value=missing_value               )

   id_qdt_conv = register_diag_field ( mod_name, &
     'qdt_conv', axes(1:3), Time, &
     'Spec humidity tendency',  'kg/kg/s',  &
                        missing_value=missing_value               )

   id_q_conv_col = register_diag_field ( mod_name, &
     'q_conv_col', axes(1:2), Time, &
    'Water vapor path tendency',   'kg/m2/s' )

   id_t_conv_col = register_diag_field ( mod_name, &
     't_conv_col', axes(1:2), Time, &
    'Column static energy tendency','W/m2' )

   id_prec_conv = register_diag_field ( mod_name, &
     'prec_conv', axes(1:2), Time, &
    'Precipitation rate',       'kg/m2/s' )

   id_snow_conv = register_diag_field ( mod_name, &
     'snow_conv', axes(1:2), Time, &
    'Frozen precip rate',       'kg/m2/s' )

   id_gust_conv = register_diag_field ( mod_name, &
     'gust_conv', axes(1:2), Time, &
    'Gustiness from deep convection ',       'm/s' )


if ( do_lsc ) then

   id_entrop_ls = register_diag_field ( mod_name, &
      'entrop_ls', axes(1:3), Time, &
      'Entropy tendency from large-scale cond',    '1/s', &
                        missing_value=missing_value               )

   id_tdt_ls = register_diag_field ( mod_name, &
     'tdt_ls', axes(1:3), Time, &
       'Temperature tendency from large-scale cond',   'deg_K/s',  &
                        missing_value=missing_value               )

   id_qdt_ls = register_diag_field ( mod_name, &
     'qdt_ls', axes(1:3), Time, &
     'Spec humidity tendency from large-scale cond', 'kg/kg/s',  &
                        missing_value=missing_value               )

   id_prec_ls = register_diag_field ( mod_name, &
     'prec_ls', axes(1:2), Time, &
    'Precipitation rate from large-scale cond',     'kg/m2/s' )

   id_snow_ls = register_diag_field ( mod_name, &
     'snow_ls', axes(1:2), Time, &
    'Frozen precip rate from large-scale cond',     'kg/m2/s' )

   id_q_ls_col = register_diag_field ( mod_name, &
     'q_ls_col', axes(1:2), Time, &
    'Water vapor path tendency from large-scale cond','kg/m2/s' )

   id_t_ls_col = register_diag_field ( mod_name, &
     't_ls_col', axes(1:2), Time, &
    'Column static energy tendency from large-scale cond','W/m2' )

endif

   id_cape = register_diag_field ( mod_name, &
     'cape', axes(1:2), Time, &
     'Convectively available potential energy',      'J/Kg')

   id_cin = register_diag_field ( mod_name, &
     'cin', axes(1:2), Time, &
     'Convective inhibition',                        'J/Kg')

   id_precip = register_diag_field ( mod_name, &
     'precip', axes(1:2), Time, &
     'Total precipitation rate',                     'kg/m2/s' )

   id_WVP = register_diag_field ( mod_name, &
     'WVP', axes(1:2), Time, &
        'Column integrated water vapor',                'kg/m2'  )

   id_rh = register_diag_field ( mod_name, &
     'rh', axes(1:3), Time, &
         'relative humidity',                            'percent',  &
                        missing_value=missing_value               )

   id_rhsurf = register_diag_field ( mod_name, &
     'rhsurf', axes(1:2), Time, &
         'Surface relative humidity',                     'percent',  &
                        missing_value=missing_value               )

!-----------------------------------------------------------------------
!---------------------------------------------------------------------
!    register the diagnostics associated with convective tracer
!    transport.
!---------------------------------------------------------------------
      allocate (id_tracerdt_conv    (num_tracers))
      allocate (id_tracerdt_conv_col(num_tracers))
      allocate (id_conv_tracer           (num_tracers))
      allocate (id_conv_tracer_col(num_tracers))

      do n = 1,num_tracers
        call get_tracer_names (MODEL_ATMOS, n, name = tracer_name,  &
                               units = tracer_units)
        if (tracers_in_mca(n)) then
          diaglname = trim(tracer_name)//  &
                        ' total tendency from moist convection'
          id_tracerdt_conv(n) =    &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'dt_conv',  &
                         axes(1:3), Time, trim(diaglname), &
                         TRIM(tracer_units)//'/s',  &
                         missing_value=missing_value)

          diaglname = trim(tracer_name)//  &
                       ' total path tendency from moist convection'
          id_tracerdt_conv_col(n) =  &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'dt_conv_col', &
                         axes(1:2), Time, trim(diaglname), &
                         TRIM(tracer_units)//'/s',   &
                         missing_value=missing_value)
         endif

         diaglname = trim(tracer_name)
         id_conv_tracer(n) =    &
                        register_diag_field ( mod_name, &
                        TRIM(tracer_name),  &
                        axes(1:3), Time, trim(diaglname), &
                        TRIM(tracer_units)      ,  &
                        missing_value=missing_value)
         diaglname =  ' column integrated' // trim(tracer_name)
         id_conv_tracer_col(n) =  &
                        register_diag_field ( mod_name, &
                        TRIM(tracer_name)//'_col', &
                        axes(1:2), Time, trim(diaglname), &
                        TRIM(tracer_units)      ,   &
                        missing_value=missing_value)
      end do


end subroutine diag_field_init

!#######################################################################


                 end module moist_processes_mod
