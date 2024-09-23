module topo_drag_mod
!=======================================================================
! TOPOGRAPHIC DRAG CLOSURE -- Garner (2005)
!=======================================================================

!-----------------------------------------------------------------------
!  Calculates horizontal velocity tendency due to topographic drag
!-----------------------------------------------------------------------


   use         fms_mod, only: mpp_npes, field_size, file_exist,  field_exist, write_version_number, stdlog, &
      mpp_pe, mpp_root_pe, error_mesg, FATAL, NOTE, read_data, write_data,  &
      open_namelist_file, close_file, check_nml_error, open_restart_file, mpp_error
   use      fms_io_mod, only: get_restart_io_mode
   use   constants_mod, only: Grav, RDgas, cp_air, PI, radian
   use   horiz_interp_mod, only: horiz_interp
   ! Added by Aman Gupta 19 Jan 2023
   use diag_manager_mod,       only:  diag_manager_init, register_diag_field, send_data
   use time_manager_mod,       only:  time_manager_init, time_type
   implicit none

   private

   character(len=128) :: version = '$Id$'
   character(len=128) :: tagname = '$Name$'

   ! Added by Aman Gupta, 19 Jan 2023
   ! Also, damping_driver namelist from Chaim does not have dd_drag in it
   integer          :: id_gwfx_ogwd, id_gwfy_ogwd
   real             :: missing_value = -999
   character(len=9) :: mod_name = 'topo_drag'


   logical :: module_is_initialized = .false.

   ! horizontal array size

   integer :: nlon, nlat
   integer :: kd=0

   ! arrays defined by topo_drag_init:

   real, allocatable, dimension(:,:) :: t11, t21, t12, t22    ! drag tensor
   real, allocatable, dimension(:,:) :: hmin, hmax
   real, allocatable, dimension(:,:) :: lon, lat

   ! parameters:

   real, parameter :: u0=1.0       ! arbitrary velocity scale for diagnostics
   real, parameter :: xl=80.0e3    ! arbitrary horiz length scale for diagnostics
   real, parameter :: ro=1.2       ! arbitrary density scale for diagnostics
   real, parameter :: lapse=Grav/cp_air ! adiabatic temperature lapse rate
   real, parameter :: tiny=1.0e-20

   real, parameter :: resolution=30.0 ! # of points per degree in topo datasets
   real, parameter :: frint=0.5

   integer, parameter :: ipts=360*resolution
   integer, parameter :: jpts=180*resolution


! parameters in namelist (topo_drag_nml):

   real :: &
      frcrit=0.7   &      ! critical value of Froude number for nonlinear flow
      ,alin=1.0     &      ! amplitude of propagating drag
      ,anonlin=5.0  &      ! amplitude of nonpropagating drag
      ,gamma=0.4    &      ! exponent in aspect ratio power law
      ,epsi=0.0     &      ! exponent in distribution power law
      ,beta=0.5     &      ! bluntness of topographic features
      ,h_frac=0.0   &      ! ratio of min to max subgrid mountain height
      ,zref_fac=1.0 &      ! adjusts level separating breaking/laminar flow
      ,tboost=1.0   &      ! surface T boost to improve PBL height estimate
      ,pcut=0.0     &      ! high-level cutoff pressure for momentum forcing
      ,samp=1.0     &      ! correction for coarse sampling of d2v/dz2
      ,max_udt=3.e-3     & ! upper bound on acceleration [m/s2]
      ,no_drag_frac=0.05 & ! fraction of lower atmosphere with no breaking
      ,max_pbl_frac=0.50  ! max fraction of lower atmosphere in PBL
   logical :: &
      do_conserve_energy=.true. & ! conserve total energy?
      ,keep_residual_flux=.true. & ! redistribute residual pseudomomentum?
      ,do_pbl_average=.false.    & ! average u,rho,N over PBL for baseflux?
      ,use_mg_scaling=.false.    & ! base flux saturates with value 'usat'?
      ,use_mask_for_pbl=.false.  & ! use bottom no_drag_layer as pbl?
      ,use_pbl_from_lock=.false.    ! use pbl height from Lock boundary scheme

   ! Aman Gupta 20 Jan 2023, changed from 64 to 128
   character(len=128)  :: topography_file='INPUT/poztopog.nc'
   character(len=128)  :: dragtensor_file='INPUT/dragelements.nc'

   logical :: do_netcdf_restart = .true. ! write to restart file

   NAMELIST /topo_drag_nml/                                               &
      frcrit, alin, anonlin, beta, gamma, epsi,                            &
      h_frac, zref_fac, tboost, pcut, samp, max_udt,                       &
      no_drag_frac, max_pbl_frac,                                          &
      do_conserve_energy, keep_residual_flux, do_pbl_average,              &
      use_mg_scaling, use_mask_for_pbl, use_pbl_from_lock,                 &    !stg
      topography_file, dragtensor_file, do_netcdf_restart


   public topo_drag, topo_drag_init, topo_drag_end


contains

!#######################################################################
! Aman Gupta 19 Jan 2023 - Added Time as an input argument
   subroutine topo_drag (                                                 &
      is, js, delt, uwnd, vwnd, atmp, &
      pfull, phalf, zfull, zhalf, &
      lat,  z_pbl, & !bqx+ z_pbl
      dtaux, dtauy,  dtemp, taux, tauy, taus, kbot, Time )

      integer, intent(in) :: is, js
      real,    intent(in) :: delt
      integer, intent(in), optional, dimension(:,:) :: kbot
      type(time_type),        intent(in)      :: Time ! Aman Gupta

! INPUT
! -----

! UWND     Zonal wind (dimensioned IDIM x JDIM x KDIM)
! VWND     Meridional wind (dimensioned IDIM x JDIM x KDIM)
! ATMP     Temperature at full levels (IDIM x JDIM x KDIM)
! PFULL    Pressure at full levels (IDIM x JDIM x KDIM)
! PHALF    Pressure at half levels (IDIM x JDIM x KDIM+1)
! ZFULL    Height at full levels (IDIM x JDIM x KDIM)
! ZHALF    Height at half levels (IDIM x JDIM x KDIM+1)

      real, intent(in), dimension(:,:,:) :: uwnd, vwnd, atmp
      real, intent(in), dimension(:,:)   :: lat,  z_pbl  !bqx+
      real, intent(in), dimension(:,:,:) :: pfull, phalf, zfull, zhalf

! OUTPUT
! ------

! DTAUX,DTAUY  Tendency of the vector wind in m/s^2 (IDIM x JDIM x KDIM)
! DTEMP        Tendency of the temperature in K/s (IDIM x JDIM x KDIM)
! TAUX,TAUY    Base momentum flux in kg/m/s^2 (IDIM x JDIM) for diagnostics
! TAUS         clipped saturation momentum flux (IDIM x JDIM x KDIM) for diagnostics

      real, intent(out), dimension(:,:)   :: taux, tauy
      real, intent(out), dimension(:,:,:) :: dtaux, dtauy, dtemp, taus

      integer, dimension(size(zfull,1),size(zfull,2)) :: kpbl, knod, kcut
      real,    dimension(size(zhalf,1),size(zhalf,2),size(zhalf,3)) :: tausat

! work arrays

      real, dimension(size(zfull,1),size(zfull,2)) :: taub, taul, taup, taun
      real, dimension(size(zfull,1),size(zfull,2)) :: frulo, fruhi, frunl, rnorm

      integer :: idim
      integer :: jdim
      integer :: i,j, k, kdim, km
      real    :: dz

! Aman Gupta 19 Jan 2023
      logical :: used


      idim = size(uwnd,1)
      jdim = size(uwnd,2)
      kdim = size(uwnd,3)

! estimate height of pbl

      call get_pbl ( atmp, zfull, pfull, phalf, kpbl, knod, kcut )
!
      if (use_pbl_from_lock) then
         kpbl = kdim
         do k = kdim, 2, -1
            where ( zfull(:,:,k) < zhalf(:,:,kdim+1) + z_pbl(:,:) )
               kpbl(:,:) = k - 1           ! the first full model level above PBL
            endwhere
         enddo
      endif

! calculate base flux

      call base_flux (                                                     &
         is, js, uwnd, vwnd, atmp, &
         lat,  z_pbl, &  !bqx+
         taux, tauy, dtaux, dtauy, &
         taub, taul, taup, taun, &
         frulo, fruhi, frunl, rnorm, &
         zfull, zhalf, pfull, phalf, kpbl )


! calculate saturation flux profile

      call satur_flux (                                                    &
         uwnd, vwnd, atmp, &
         taup, taub, tausat, &
         frulo, fruhi, frunl, &
         dtaux, dtauy, zfull, pfull, phalf, kpbl, kcut )

! calculate momentum tendency

      call topo_drag_tend (                                                &
         delt, uwnd, vwnd, atmp, &
         taux, tauy, taul, taun, tausat, &
         dtaux, dtauy,  dtemp, zfull, zhalf, pfull, phalf, kpbl ) !bqx+ dtaux_np, dtauy_np

! put saturation flux profile into 'taus' for diagnostics

      do k=1,kdim
!     taus(:,:,k) = 0.5*rnorm(:,:)*(tausat(:,:,k) + tausat(:,:,k+1))
         taus(:,:,k) = tausat(:,:,k+1)*taub(:,:)/taul(:,:)
      enddo

! put total drag into 'taux,tauy' for diagnostics

      taup = taup - tausat(:,:,1)
      taub = (taup + taun)/taul
      taux = taux*taub
      tauy = tauy*taub


! Added by Aman Gupta - writing to netCDF output files
      if (id_gwfx_ogwd > 0) then
         used = send_data (id_gwfx_ogwd, dtaux, Time, is, js, 1)
      endif
      if (id_gwfy_ogwd > 0) then
         used = send_data (id_gwfy_ogwd, dtauy, Time, is, js, 1)
      endif

   end subroutine topo_drag

!=======================================================================

   subroutine base_flux (                                                 &
      is, js, uwnd, vwnd, atmp, &
      lat,  z_pbl, & !bqx
      taux, tauy, dtaux, dtauy, &
      taub, taul, taup, taun, &
      frulo, fruhi, frunl, rnorm, &
      zfull, zhalf, pfull, phalf, kpbl )

      integer, intent(in) :: is, js
      real, intent(in),  dimension(:,:,:) :: uwnd, vwnd, atmp
      real, intent(in),  dimension(:,:)   :: lat,  z_pbl
      real, intent(in),  dimension(:,:,:) :: zfull, zhalf, pfull, phalf
      real, intent(out), dimension(:,:)   :: taux, tauy
      real, intent(out), dimension(:,:,:) :: dtaux, dtauy
      real, intent(out), dimension(:,:)   :: taub, taul, taup, taun
      real, intent(out), dimension(:,:)   :: frulo, fruhi, frunl, rnorm
      integer, intent(in), dimension(:,:) :: kpbl

      real, dimension(size(uwnd,1),size(uwnd,2)) :: ubar, vbar

      integer :: i, idim, id
      integer :: j, jdim, jd
      integer :: k, kdim, kb, kbp, kt, km

      real :: usat, bfreq2, bfreq, dphdz, vtau, d2udz2, d2vdz2
      real :: dzfull, dzhalf, dzhalf1, dzhalf2, density
      real :: frmin, frmax, frmed, frumin, frumax, frumed, fruclp, fruclm
      real :: rnormal, gterm, hterm, fru0, frusat
      real :: usum, vsum, n2sum, delp

      idim = size(uwnd,1)
      jdim = size(uwnd,2)
      kdim = size(uwnd,3)

! compute base flux

      do j=1,jdim
         do i=1,idim
            usum = 0.
            vsum = 0.
            kt = kpbl(i,j)
            kb = max(kd,kt)
            do k=kt,kb
               delp = phalf(i,j,k+1) - phalf(i,j,k)
               usum = usum + uwnd(i,j,k)*delp
               vsum = vsum + vwnd(i,j,k)*delp
            enddo
            ubar(i,j) = usum/(phalf(i,j,kb+1) - phalf(i,j,kt))
            vbar(i,j) = vsum/(phalf(i,j,kb+1) - phalf(i,j,kt))
         enddo
      enddo



      do j=1,jdim
         jd = js+j-1
         do i=1,idim
            id = is+i-1
            kt = kpbl(i,j)
            kb = max(kd,kt)
            kbp = min(kdim, kb+1) !bqx+
            dzfull = zhalf(i,j,kt) - zhalf(i,j,kb+1)
            density = (phalf(i,j,kb+1) - phalf(i,j,kt))/(Grav*dzfull)
            dzfull = zfull(i,j,kt-1) - zfull(i,j,kbp)
            bfreq2 = Grav*((atmp(i,j,kt-1) - atmp(i,j,kbp))/dzfull+lapse)/&
               (0.5*(atmp(i,j,kt-1) + atmp(i,j,kbp)))
!
            bfreq = sqrt(max(tiny, bfreq2))

!       included 'alin' 4/2015

            taux(i,j) = (ubar(i,j)*t11(id,jd) + vbar(i,j)*t21(id,jd))      &
               *bfreq*density
            tauy(i,j) = (ubar(i,j)*t12(id,jd) + vbar(i,j)*t22(id,jd))      &
               *bfreq*density

            taub(i,j) = max(tiny, sqrt(taux(i,j)**2 + tauy(i,j)**2))

!       min/max Froude numbers based on low-level flow

            vtau = max(tiny, -(ubar(i,j)*taux(i,j)                         &
               + vbar(i,j)*tauy(i,j))/taub(i,j))
            frmax = hmax(id,jd)*bfreq / vtau
            frmin = hmin(id,jd)*bfreq / vtau
            frmed = frcrit + frint

!       linear momentum flux associated with min/max Froude numbers

            dphdz = bfreq / vtau
            usat = sqrt(density/ro) * vtau / sqrt(dphdz*xl)
            frusat = frcrit*usat

            frumin = frmin*usat
            frumax = frmax*usat
            frumed = frmed*usat

            frumax = max(frumax,frumin + tiny)
            fruclp = min(frumax,max(frumin,frusat))
            fruclm = min(frumax,max(frumin,frumed))
            fru0 = (u0/vtau)*usat

!       total drag in linear limit

            rnormal =                                                      &
               (frumax**(2.0*gamma - epsi)                         &
               - frumin**(2.0*gamma - epsi))/(2.0*gamma - epsi)
            rnormal = fru0**gamma * ro/rnormal

            taul(i,j) =                                                    &
               (frumax**(2.0 + gamma - epsi)                       &
               - frumin**(2.0 + gamma - epsi))/(2.0 + gamma - epsi)

!       separate propagating and nonpropagating parts of total drag

            gterm = frusat**(beta + 1.0)*                                  &
               (frumax**(gamma - epsi - beta)                         &
               - fruclp**(gamma - epsi - beta))/(gamma - epsi - beta)

            taup(i,j) =  alin *                                             &
               ( (fruclp**(2.0 + gamma - epsi)                       &
               - frumin**(2.0 + gamma - epsi))/(2.0 + gamma - epsi) &
               + frusat*gterm )

            taun(i,j) = anonlin*usat/(1.0 + beta) *                        &
               ( (frumax**(1.0 + gamma - epsi)                       &
               - fruclp**(1.0 + gamma - epsi))/(1.0 + gamma - epsi) &
               - gterm )

!       5/2015 mg option: depth of blocking ~ U/N, not h

            if (use_mg_scaling) taun(i,j) = taun(i,j)/max(frmax,frcrit)

            fruhi(i,j) = frumax
            frulo(i,j) = frumin
            frunl(i,j) = frusat
            rnorm(i,j) = rnormal

         enddo
      enddo

! wind component opposite the drag at full levels (stored as 'dtaux')

      do k=1,kdim
         do j=1,jdim
            do i=1,idim
               dtaux(i,j,k) =                                              &
                  -(uwnd(i,j,k)*taux(i,j) + vwnd(i,j,k)*tauy(i,j))/taub(i,j)
            enddo
         enddo
      enddo

! curvature of wind at full levels (stored as 'dtauy')

      dtauy = 0.

      do j=1,jdim
         do i=1,idim
!        kt = kpbl(i,j)
            kt = min(kdim-1, kpbl(i,j)) !bqx+
            do k=2,kt
               dzfull = zhalf(i,j,k) - zhalf(i,j,k+1)
               dzhalf1 = zfull(i,j,k-1) - zfull(i,j,k)
               dzhalf2 = zfull(i,j,k) - zfull(i,j,k+1)
               d2udz2 = ((uwnd(i,j,k-1) - uwnd(i,j,k  ))/dzhalf1           &
                  - (uwnd(i,j,k  ) - uwnd(i,j,k+1))/dzhalf2)/dzfull
               d2vdz2 = ((vwnd(i,j,k-1) - vwnd(i,j,k  ))/dzhalf1           &
                  - (vwnd(i,j,k  ) - vwnd(i,j,k+1))/dzhalf2)/dzfull
               dtauy(i,j,k) = -(d2udz2*taux(i,j) + d2vdz2*tauy(i,j))/      &
                  taub(i,j)
            enddo

         enddo
      enddo

   end subroutine base_flux

!=======================================================================

   subroutine satur_flux (                                                &
      uwnd, vwnd, atmp, &
      taup, taub, tausat, &
      frulo, fruhi, frunl, &
      dtaux, dtauy, zfull, pfull, phalf, kpbl, kcut )

      real, intent(in),  dimension (:,:,:) :: uwnd, vwnd, atmp
      real, intent(in),  dimension (:,:,:) :: dtaux, dtauy
      real, intent(in),  dimension (:,:,:) :: zfull, pfull, phalf
      real, intent(in),  dimension (:,:)   :: taup
      real, intent(out), dimension (:,:)   :: taub
      real, intent(out), dimension (:,:,:) :: tausat
      real, intent(in),  dimension (:,:)   :: frulo, fruhi, frunl
      integer, intent(in), dimension (:,:) :: kpbl, kcut

      real, dimension(size(zfull,1),size(zfull,2)) :: usat

      real :: dzhalf, gterm, gterm0, density
      real :: bfreq2, bfreq, vtau, d2vtau, dphdz, xl1
      real :: frumin, frumax, fruclp, frusat, frusat0, fruclp0

      integer :: i, idim
      integer :: j, jdim
      integer :: k, kdim, k1

      idim = size(uwnd,1)
      jdim = size(uwnd,2)
      kdim = size(uwnd,3)

! get vertical profile of propagating part of momentum flux

      usat = frunl/frcrit

      do k=kdim,2,-1
         do j=1,jdim
            do i=1,idim

!          buoyancy frequency, velocity and density at half levels

               dzhalf = zfull(i,j,k-1) - zfull(i,j,k)
               density = (pfull(i,j,k) - pfull(i,j,k-1))/(Grav*dzhalf)
               bfreq2 = Grav*                                              &
                  ((atmp(i,j,k-1) - atmp(i,j,k))/dzhalf + lapse)/ &
                  (0.5*(atmp(i,j,k-1) + atmp(i,j,k)))
               bfreq = sqrt(max(tiny, bfreq2))

               vtau = max(tiny, 0.5*(dtaux(i,j,k-1) + dtaux(i,j,k)))

!          WKB correction of vertical wavelength

               d2vtau = 0.5*(dtauy(i,j,k-1) + dtauy(i,j,k))
               xl1 = xl*max(0.5, min(2.0, 1.0 - samp*vtau*d2vtau/(bfreq*bfreq)))

!          min/max and critical momentum flux values at half levels

               dphdz = bfreq / vtau
               usat(i,j) = min(usat(i,j),sqrt(density/ro) * vtau/sqrt(dphdz*xl1))
               frusat = frcrit*usat(i,j)

               frumin = frulo(i,j)
               frumax = fruhi(i,j)
               fruclp = min(frumax,max(frumin,frusat))
               frusat0 = frunl(i,j)
               fruclp0 = min(frumax,max(frumin,frusat0))

!          propagating part of momentum flux (from WKB or EP)

               gterm0 = (frumax**(gamma - epsi - beta)                     &
                  - fruclp0**(gamma - epsi - beta))/(gamma - epsi - beta)
               gterm = (fruclp0**(gamma - epsi)                            &
                  - fruclp**(gamma - epsi))/(gamma - epsi)

               tausat(i,j,k) = alin *                                      &
                  ( (fruclp**(2.0 + gamma - epsi)                       &
                  - frumin**(2.0 + gamma - epsi))/(2.0 + gamma - epsi) &
                  + frusat**2.0*(gterm0*frusat0**beta + gterm) )
            enddo
         enddo
      enddo

! make propagating flux constant with height in zero-drag top layer
! changed 5/2014

      k1 = maxval(kcut)
      do k=k1,1,-1
         where (k <= kcut)
            tausat(:,:,k) = tausat(:,:,k+1)
         endwhere
      enddo

! make propagating flux constant with height in zero-drag surface layer

      k1 = minval(kpbl)
      do k=kdim+1,k1+1,-1
         where (k > kpbl)
            tausat(:,:,k) = taup
         endwhere
      enddo

! redistribute residual forcing

      if ( keep_residual_flux ) then
         taub(:,:) = tausat(:,:,1)/(phalf(:,:,kdim+1) - phalf(:,:,1))
         do k=1,kdim
            tausat(:,:,k) = tausat(:,:,k)                                  &
               - taub(:,:)*(phalf(:,:,kdim+1) - phalf(:,:,k))
         enddo
      endif

   endsubroutine satur_flux

!=======================================================================

   subroutine topo_drag_tend (                                            &
      delt, uwnd, vwnd, atmp, &
      taux, tauy, taul, taun, tausat, &
      dtaux, dtauy, dtemp, zfull, zhalf, &
      pfull, phalf, kpbl )

      real, intent(in) :: delt
      real, intent(in), dimension(:,:,:)    :: uwnd, vwnd, atmp
      real, intent(in), dimension(:,:,:)    :: zfull, zhalf, pfull, phalf
      real, intent(in), dimension(:,:)      :: taux, tauy, taul, taun
      real, intent(inout), dimension(:,:,:) :: tausat
      real, intent(inout), dimension(:,:,:) :: dtaux, dtauy, dtemp
      integer, intent(in), dimension (:,:)  :: kpbl

      real, parameter :: bfmin=0.7e-2, bfmax=1.7e-2  ! min/max buoyancy freq [1/s]
      real, parameter :: vvmin=1.0                   ! minimum surface wind [m/s]

      integer,dimension(size(zfull,1),size(zfull,2)) :: kref
      real :: dzhalf, zlast, rscale, phase, bfreq, bfreq2, vtau
      real :: gfac, gfac1, dp, weight, wtsum, taunon, taunon1

      integer :: i, idim
      integer :: j, jdim
      integer :: k, kdim, kr, kt

      real,dimension(size(zfull,1),size(zfull,2)) :: dx, dy

      idim = size(uwnd,1)
      jdim = size(uwnd,2)
      kdim = size(uwnd,3)

! find reference level for non-propagating drag (z ~ pi U/N)

!the following re-orients the drag to align with low-level wind
!do j=1,jdim
!do i=1,idim
!k = knod(i,j)
!gfac = sqrt( (taux(i,j)**2 + tauy(i,j)**2) / max (tiny, uwnd(i,j,k)**2 + vwnd(i,j,k)**2) )
!dx(i,j) = gfac * uwnd(i,j,k)
!dy(i,j) = gfac * vwnd(i,j,k)
!enddo
!enddo

      do j=1,jdim
         do i=1,idim
            k = kpbl(i,j)
!stg        k = kdim
            phase = 0.0
            zlast = zhalf(i,j,k)
            do while (phase <= PI*zref_fac .and. k > 1)
               k = k-1
               vtau = 0.5*(dtaux(i,j,k-1) + dtaux(i,j,k))
               dzhalf = zfull(i,j,k-1) - zfull(i,j,k)
               bfreq2 = Grav*                                              &
                  ((atmp(i,j,k-1) - atmp(i,j,k))/dzhalf + lapse)/ &
                  (0.5*(atmp(i,j,k-1) + atmp(i,j,k)))
               bfreq = sqrt(max(tiny, bfreq2))
               rscale = max(bfmin, min(bfmax, bfreq))/max(vvmin, vtau)
               dzhalf = zfull(i,j,k-1) - zlast
               phase = phase + dzhalf*rscale
               zlast = zfull(i,j,k-1)
            enddo
            kref(i,j) = k
         enddo
      enddo

! CALCULATE DECELERATION DUE TO PROPAGATING DRAG (~-rho^-1 dtau/dz)

      do k=1,kdim
         do j=1,jdim
            do i=1,idim
               dp = phalf(i,j,k+1) - phalf(i,j,k)
               gfac = tausat(i,j,k+1) - tausat(i,j,k)
               gfac1 = gfac*Grav/(dp*taul(i,j))
               dtaux(i,j,k) = gfac1*taux(i,j)
               dtauy(i,j,k) = gfac1*tauy(i,j)
!dtaux(i,j,k) = -gfac1*dx(i,j)
!dtauy(i,j,k) = -gfac1*dy(i,j)
            enddo
         enddo
      enddo

! CALCULATE DECELERATION DUE TO NON-PROPAGATING DRAG

      do j=1,jdim
         do i=1,idim
            kr = kref(i,j)
            kt = kpbl(i,j)
!stg        kt = kdim
            gfac = taun(i,j)/taul(i,j) * Grav
            wtsum = 0.0
            do k=kr,kt
               dp = phalf(i,j,k+1) - phalf(i,j,k)
               weight = pfull(i,j,k) - phalf(i,j,kr)
               wtsum = wtsum + dp*weight
            enddo
            taunon = 0.0
            do k=kr,kt
               weight = pfull(i,j,k) - phalf(i,j,kr)
               gfac1 = gfac*weight/wtsum
               dtaux(i,j,k) = dtaux(i,j,k) + gfac1*taux(i,j)
               dtauy(i,j,k) = dtauy(i,j,k) + gfac1*tauy(i,j)
!bqx
!dtaux(i,j,k) = dtaux(i,j,k) - gfac1*dx(i,j)
!dtauy(i,j,k) = dtauy(i,j,k) - gfac1*dy(i,j)
               dp = phalf(i,j,k+1) - phalf(i,j,k)
               taunon = taunon + gfac1*dp
               taunon1 = taunon*taul(i,j)/Grav
               tausat(i,j,k) = tausat(i,j,k) + taunon1
            enddo
            do k=kt+1,kdim
               tausat(i,j,k) = tausat(i,j,k) + taunon1
            enddo
         enddo
      enddo

      dtaux = max(-max_udt, min(max_udt, dtaux))   !stg
      dtauy = max(-max_udt, min(max_udt, dtauy))

! CALCULATE HEATING TO CONSERVE TOTAL ENERGY

      if (do_conserve_energy) then
         dtemp = -((uwnd + 0.5*delt*dtaux)*dtaux                           &
            + (vwnd + 0.5*delt*dtauy)*dtauy)/cp_air
      else
         dtemp = 0.0
      endif

   end subroutine topo_drag_tend

!=======================================================================

   subroutine get_pbl ( atmp, zfull, pfull, phalf, kpbl, knod, kcut )

      integer, intent(out), dimension(:,:) :: kpbl, knod, kcut
      real, intent(in), dimension(:,:,:)   :: atmp
      real, intent(in), dimension(:,:,:)   :: zfull, pfull, phalf

      real, dimension(size(pfull,1),size(pfull,2)) :: ppbl, pbot
      real, dimension(size(pfull,1),size(pfull,2)) :: tbot, zbot

      integer :: i, idim
      integer :: j, jdim
      integer :: k, kdim

      idim = size(atmp,1)
      jdim = size(atmp,2)
      kdim = size(atmp,3)

      do j=1,jdim
         do i=1,idim
            ppbl(i,j) = (1.0 - max_pbl_frac)*phalf(i,j,kdim+1)
            pbot(i,j) = (1.0 - no_drag_frac)*phalf(i,j,kdim+1)
            tbot(i,j) = atmp(i,j,kdim) + tboost
            zbot(i,j) = zfull(i,j,kdim)
         enddo
      enddo

! find highest model level in no-drag surface layer

      knod = kdim-1

      do k=kdim-2,2,-1
         where ( pfull(:,:,k) >= pbot(:,:) )
            knod = k
         endwhere
      enddo

! find lowest model level in no-drag top layer

      kcut = 1

      do k=2,kdim
         where ( pfull(:,:,k) <= pcut )
            kcut = k
         endwhere
      enddo

      if (kd == 0 .and. do_pbl_average) kd = kdim-1

      if ( use_mask_for_pbl ) then
         kpbl = knod
         return
      endif

! find the first layer above PBL

      kpbl = kdim-1

      do k=kdim-2,2,-1
         where ( pfull(:,:,k) >= ppbl(:,:) .and.                          &
            tbot(:,:) - atmp(:,:,k) > lapse*(zfull(:,:,k) - zbot(:,:)) )
            kpbl = k -1
         endwhere
      enddo

   end subroutine get_pbl

!=======================================================================
! function definition changed by Aman Gupta, 19 Jan 2023
!subroutine topo_drag_init (lonb, latb)
   subroutine topo_drag_init (lonb, latb, Time, axes)

      real, intent(in), dimension(:) :: lonb, latb
! Added by Aman Gupta
      integer, dimension(4),   intent(in)      :: axes
      type(time_type),         intent(in)      :: Time

      character(len=128) :: msg
      character(len=64)  :: restart_fname='INPUT/topo_drag.res.nc'
      character(len=3)   :: tensornames(4) = (/ 't11', 't21', 't12', 't22' /)

      logical :: found_field(4)

      real, parameter :: bfscale=1.0e-2      ! buoyancy frequency scale [1/s]

      real :: xdat(ipts+1), ydat(jpts+1)
      real :: zdat(ipts,jpts)
      real :: zout(size(lonb,1)-1,size(latb,1)-1)

      real :: exponent, hmod

      integer :: n
      integer :: io, ierr,  unit
      integer :: i, j
      integer :: siz(4)

      if (module_is_initialized) return
! Added by Aman Gupta 19 Jan 2023
      call diag_manager_init
      call time_manager_init

      nlon = size(lonb,1)-1
      nlat = size(latb,1)-1
!  print *,'nlon=', nlon
!  print *,'nlat=', nlat
!---------------------------------------------------------------------
! --- Read namelist
!---------------------------------------------------------------------
      if( file_exist( 'input.nml' ) ) then
! -------------------------------------
         unit = open_namelist_file()
         ierr = 1
         do while( ierr .ne. 0 )
            read ( unit,  nml = topo_drag_nml, iostat = io, end = 10 )
            ierr = check_nml_error(io,'topo_drag_nml')
         end do
10       continue
         call close_file ( unit )
         call get_restart_io_mode(do_netcdf_restart)

! -------------------------------------
      end if

! write version number and namelist to logfile

      call write_version_number(version, tagname)
      if(mpp_pe() == mpp_root_pe()) write (stdlog(), nml=topo_drag_nml)


      allocate (t11(nlon,nlat))
      allocate (t21(nlon,nlat))
      allocate (t12(nlon,nlat))
      allocate (t22(nlon,nlat))
      allocate (hmin(nlon,nlat))
      allocate (hmax(nlon,nlat))

      if (gamma == beta + epsi) gamma = gamma + tiny

! read restart file

      if ( file_exist( 'INPUT/topo_drag.res.nc' ) ) then
         if (mpp_pe() == mpp_root_pe()) call mpp_error ('topo_drag_mod', &
            'Reading NetCDF formatted restart file: INPUT/topo_drag.res.nc', NOTE)
         call read_data ('INPUT/topo_drag.res.nc', 't11', t11)
         call read_data ('INPUT/topo_drag.res.nc', 't12', t12)
         call read_data ('INPUT/topo_drag.res.nc', 't21', t21)
         call read_data ('INPUT/topo_drag.res.nc', 't22', t22)
         call read_data ('INPUT/topo_drag.res.nc', 'hmin', hmin)
         call read_data ('INPUT/topo_drag.res.nc', 'hmax', hmax)
      else if ( file_exist( 'INPUT/topo_drag.res' ) ) then
         if (mpp_pe() == mpp_root_pe()) call mpp_error ('topo_drag_mod', &
            'Reading native formatted restart file.', NOTE)
         unit = open_restart_file('INPUT/topo_drag.res','read')
         call read_data(unit, t11)
         call read_data(unit, t12)
         call read_data(unit, t21)
         call read_data(unit, t22)
         call read_data(unit, hmin)
         call read_data(unit, hmax)
         call close_file(unit)
      else if (file_exist(topography_file) .and.                           &
         file_exist(dragtensor_file)) then

!    read and interpolate topography datasets

         if (mpp_pe() == mpp_root_pe()) then
            write ( msg, '("Reading topography file: ",a)')              &
               trim(topography_file)
            call error_mesg('topo_drag_mod', msg, NOTE)
         endif

! check for correct field size in topography
         call field_size (topography_file, 'hpoz', siz)


         if (siz(1) /= ipts .or. siz(2) /= jpts) then
            call error_mesg('topo_drag_mod', 'Field \"hpoz\" in file '// &
               trim(topography_file)//' has the wrong size', FATAL)
         endif


         do i=1,ipts+1
            xdat(i) = (i-1)/resolution / Radian
         enddo
         do j=1,jpts+1
            ydat(j) = (-90.0 + (j-1)/resolution) / Radian
         enddo

! initialize horizontal interpolation
! Note: interp_method will be conservative for lat/lon grid
!       and bilinear for all other grids

         call read_data (topography_file, 'hpoz', zdat, no_domain=.true.)


         exponent = 2.0 - gamma
         zdat = max(0.0, zdat)**exponent



         call horiz_interp ( zdat, xdat, ydat, lonb, latb, zout )


         hmax = abs(zout)**(1.0/exponent) * sqrt(2.0/exponent)
         hmin = hmax*h_frac

         if (mpp_pe() == mpp_root_pe()) then
            write ( msg, '("Reading drag tensor file: ",a)')             &
               trim(dragtensor_file)
            call error_mesg('topo_drag_mod', msg, NOTE)
         endif

! check for correct field size in tensor file
         call field_size (dragtensor_file, tensornames(1), siz)
         if (siz(1) /= ipts .or. siz(2) /= jpts) then
            call error_mesg('topo_drag_mod', 'Field \"'//tensornames(1)// &
               '\" in file '//trim(dragtensor_file)//' has the wrong size', FATAL)
         endif

         do n=1,4
            found_field(n) = field_exist(dragtensor_file, tensornames(n))
            if (.not. found_field(n)) cycle
            call read_data (dragtensor_file, tensornames(n), zdat, no_domain=.true.)

            call horiz_interp ( zdat, xdat, ydat, lonb, latb, zout )
!call horiz_interp ( zdat, xdat, ydat, lonb, latb, zout,        &
!                                 interp_method='conservative' )
            if ( tensornames(n) == 't11' ) then
               t11 = zout/bfscale
            else if ( tensornames(n) == 't21' ) then
               t21 = zout/bfscale
            else if ( tensornames(n) == 't12' ) then
               t12 = zout/bfscale
            else if ( tensornames(n) == 't22' ) then
               t22 = zout/bfscale
            endif
         enddo
         if (.not. found_field(3)) t12 = t21


      else
         print *,topography_file
         print *,dragtensor_file

         call error_mesg  ('topo_drag_init',                                &
            'Aman: No sub-grid orography available for topo_drag', FATAL)

      endif

! Writing to the file, added by Aman Gupta, 19 Jan 2023
!--------------------------------------------------------------------
!    initialize netcdf diagnostic fields.
!-------------------------------------------------------------------
      id_gwfx_ogwd =  &
         register_diag_field (mod_name, 'gwfu_ogwd', axes(1:3), Time, &
         'gravity wave forcing on mean zonal flow', &
         'm/s^2',  missing_value=missing_value)
      id_gwfy_ogwd =  &
         register_diag_field (mod_name, 'gwfv_ogwd', axes(1:3), Time, &
         'gravity wave forcing on mean meridional flow', &
         'm/s^2',  missing_value=missing_value)


      module_is_initialized = .true.

   end subroutine topo_drag_init

!=======================================================================

   subroutine topo_drag_end

      integer :: unit

      if(.not.module_is_initialized) return
      if(do_netcdf_restart) then
         if (mpp_pe() == mpp_root_pe()) call mpp_error ('topo_drag_mod', &
            'Writing NetCDF formatted restart file: RESTART/topo_drag.res.nc', NOTE)
         call write_data('RESTART/topo_drag.res.nc', 't11', t11)
         call write_data('RESTART/topo_drag.res.nc', 't12', t12)
         call write_data('RESTART/topo_drag.res.nc', 't21', t21)
         call write_data('RESTART/topo_drag.res.nc', 't22', t22)
         call write_data('RESTART/topo_drag.res.nc', 'hmin', hmin)
         call write_data('RESTART/topo_drag.res.nc', 'hmax', hmax)
      else
         if (mpp_pe() == mpp_root_pe()) call mpp_error ('topo_drag_mod', &
            'Writing native formatted restart file.', NOTE)
         unit = open_restart_file('RESTART/topo_drag.res','write')
         call write_data(unit, t11)
         call write_data(unit, t12)
         call write_data(unit, t21)
         call write_data(unit, t22)
         call write_data(unit, hmin)
         call write_data(unit, hmax)
         call close_file(unit)
      endif
      deallocate(t11)
      deallocate(t12)
      deallocate(t21)
      deallocate(t22)
      deallocate(hmin)
      deallocate(hmax)
      module_is_initialized = .false.


   end subroutine topo_drag_end

endmodule topo_drag_mod
