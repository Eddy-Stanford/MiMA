
!> Upper-boundary damping and gravity-wave drag.
!>
!> Adds to the tendencies of the winds (and temperature) up to four optional terms:
!>
!> * `do_rayleigh`: Rayleigh friction that damps the winds towards zero at pressures below
!>   `sponge_pbottom`, with a rate that increases quadratically from zero at `sponge_pbottom`
!>   to the rate of time scale `trayfric` at zero pressure;
!> * `do_mg_drag`: orographic gravity-wave drag (`mg_drag_mod`);
!> * `do_cg_drag`: non-orographic gravity-wave drag after Alexander and Dunkerton (1999)
!>   (`cg_drag_mod`);
!> * `do_const_drag`: an idealized, time-independent "gravity-wave" drag on u at pressures
!>   below e (about 2.7) hPa, modelled on the Alexander-Dunkerton winter average, linear in ln(p),
!>   with a cubic latitudinal profile and a seasonal cosine in time.
!>
!> With `do_conserve_energy`, the kinetic energy removed by the Rayleigh friction heats
!> the air.
!>
!> Namelist: `damping_driver_nml`
!> ([namelist reference](https://eddy-stanford.github.io/MiMA/Parameters/#damping_driver_nml)).
!>
!> References:
!>
!> * Alexander, M. J., and T. J. Dunkerton, 1999: A spectral parameterization
!>   of mean-flow forcing due to breaking gravity waves. J. Atmos. Sci., 56,
!>   4167-4180.
module damping_driver_mod

  use mg_drag_mod, only: mg_drag, mg_drag_init, mg_drag_end
  use cg_drag_mod, only: cg_drag_init, cg_drag_calc, cg_drag_end
  use fms_mod, only: mpp_pe, mpp_root_pe, stdlog, &
                     write_version_number, &
                     input_nml_file, error_mesg, &
                     check_nml_error, &
                     FATAL
  use diag_manager_mod, only: register_diag_field, &
                              register_static_field, send_data
  use time_manager_mod, only: time_type, get_time, length_of_year !mj
  use constants_mod, only: cp_air, grav, PI
  use mpp_domains_mod, only: domain2d

  implicit none
  private

  public damping_driver, damping_driver_init, damping_driver_end

!-----------------------------------------------------------------------
!---------------------- namelist ---------------------------------------

  real     :: trayfric = -0.5  !! Rayleigh friction time scale: [s] if > 0, [days] if < 0
! mj pk02-like sponge   integer  :: nlev_rayfric = 1
  integer  :: nlev_rayfric
  ! number of levels at the top of the model where Rayleigh friction is applied; set in
  ! damping_driver_init from sponge_pbottom
  real :: sponge_pbottom = 50. !! [Pa] bottom of the Rayleigh sponge
  logical  :: do_mg_drag = .false.  !! orographic gravity-wave drag (`mg_drag_nml`)
!epg: Use cg_drag.f90, GFDL's version of the Alexander and Dunkerton 1999
!     Non-orographic gravity wave parameterization, updated as for Cohen et al. 2013
! mj actively choose rayleigh friction
  logical  :: do_rayleigh = .false.  !! Rayleigh friction (sponge) at the top of the model
  logical  :: do_cg_drag = .true.
  !! non-orographic (convective) gravity-wave drag (`cg_drag_nml`)
  logical  :: do_const_drag = .false.
  !! idealized seasonal "gravity-wave" drag in the stratosphere
  real     :: const_drag_amp = 3.e-04  !! [m/s2] amplitude of the constant drag
  real     :: const_drag_off = 0.  !! offset of its latitudinal profile
  logical  :: do_conserve_energy = .true.  !! heat the air by the momentum lost to the damping

  namelist /damping_driver_nml/ trayfric, &
    do_rayleigh, sponge_pbottom, & ! mj
    do_cg_drag, &
    do_mg_drag, do_conserve_energy, &
    do_const_drag, const_drag_amp, const_drag_off    !mj

!-----------------------------------------------------------------------
!----- id numbers for diagnostic fields -----

  integer :: id_udt_rdamp, id_vdt_rdamp, &
             id_udt_gwd, id_vdt_gwd, &
             id_sgsmtn, &
             id_udt_cgwd, id_taus, &
             id_udt_cnstd                 !mj

  integer :: id_tdt_diss_rdamp, id_diss_heat_rdamp, &
             id_tdt_diss_gwd, id_diss_heat_gwd

  integer :: id_taubx, id_tauby

!----- missing value for all fields ------

  real :: missing_value = -999.

  character(len=7) :: mod_name = 'damping'

!-----------------------------------------------------------------------
!mj actively choose rayleigh - is now in namelist
! logical :: do_rayleigh

  real, parameter ::  daypsec = 1./86400.
  logical :: module_is_initialized = .false.

  real :: rfactr

!   note:
!     rfactr = coeff. for damping momentum at the top level

  character(len=128) :: version = '$Id: damping_driver.f90,v 10.0 2003/10/24 22:00:25 fms Exp $'
  character(len=128) :: tagname = '$Name: lima $'

!mj cg_drag alarm
  integer :: Time_lastcall, dt_integer, days, seconds
!-----------------------------------------------------------------------

contains

!#######################################################################

  !> Adds the tendencies from the Rayleigh sponge, the orographic and non-orographic
  !> gravity-wave drag and the constant drag (whichever are switched on) to `udt`, `vdt`
  !> and `tdt`, and sends their diagnostics.
  subroutine damping_driver(is, js, lat, Time, delt, pfull, phalf, zfull, zhalf, &
                            u, v, t, q, r, udt, vdt, tdt, qdt, rdt, &
                            mask, kbot)

!-----------------------------------------------------------------------
    integer, intent(in)                :: is, js  !! starting i,j indices of the physics window
    real, dimension(:, :), intent(in)           :: lat  !! latitudes [rad]
    type(time_type), intent(in)                :: Time  !! current time
    real, intent(in)                :: delt  !! physics time step [s]
    real, intent(in), dimension(:, :, :)   :: pfull, phalf, &
                                              zfull, zhalf, &
                                              u, v, t, q
    !! `pfull`, `phalf`: pressure at full and half levels [Pa]; `zfull`, `zhalf`: height at
    !! full and half levels [m]; `u`, `v`: zonal and meridional wind [m/s]; `t`: temperature
    !! [K]; `q`: specific humidity [kg/kg] (not used)
    real, intent(in), dimension(:, :, :, :) :: r  !! tracers (not used)
    real, intent(inout), dimension(:, :, :)   :: udt, vdt, tdt, qdt
    !! tendencies of u [m/s2], v [m/s2], temperature [K/s] and specific humidity [kg/kg/s]
    !! (`qdt` is not changed)
    real, intent(inout), dimension(:, :, :, :) :: rdt  !! tracer tendencies (not changed)
    real, intent(in), dimension(:, :, :), optional :: mask  !! mask for the diagnostics
    integer, intent(in), dimension(:, :), optional :: kbot  !! index of the lowest model level

!-----------------------------------------------------------------------
    real, dimension(size(udt, 1), size(udt, 2))             :: diag2
    real, dimension(size(udt, 1), size(udt, 2))             :: taubx, tauby
    real, dimension(size(udt, 1), size(udt, 2), size(udt, 3)) :: taus
    real, dimension(size(udt, 1), size(udt, 2), size(udt, 3)) :: utnd, vtnd, &
                                                                 ttnd, pmass
    logical :: used

    real, dimension(size(udt, 1), size(udt, 2), size(udt, 3) + 1) :: p_pass, &
                                                                     t_pass
    integer :: k, j, i
!-----------------------------------------------------------------------
!mj constant drag TOA
    real :: minp, cosday
    integer :: seconds, days, daysperyear
!-----------------------------------------------------------------------

    if (.not. module_is_initialized) call error_mesg('damping_driver', &
                                                     'damping_driver_init must be called first', FATAL)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!----------------- r a y l e i g h   d a m p i n g ---------------------
!-----------------------------------------------------------------------
    if (do_rayleigh) then

! mj pk02-like sponge
      call rayleigh(delt, pfull, u, v, utnd, vtnd, ttnd)
      udt = udt + utnd
      vdt = vdt + vtnd
      tdt = tdt + ttnd

!----- diagnostics -----

      if (id_udt_rdamp > 0) then
        used = send_data(id_udt_rdamp, utnd, Time, is, js, 1, &
                         rmask=mask)
      end if

      if (id_vdt_rdamp > 0) then
        used = send_data(id_vdt_rdamp, vtnd, Time, is, js, 1, &
                         rmask=mask)
      end if

      if (id_tdt_diss_rdamp > 0) then
        used = send_data(id_tdt_diss_rdamp, ttnd, Time, is, js, 1, &
                         rmask=mask)
      end if

      if (id_diss_heat_rdamp > 0) then
        do k = 1, size(u, 3)
          pmass(:, :, k) = phalf(:, :, k + 1) - phalf(:, :, k)
        end do
        diag2 = cp_air/grav*sum(ttnd*pmass, 3)
        used = send_data(id_diss_heat_rdamp, diag2, Time, is, js)
      end if

    end if
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!--------- m t n   g r a v i t y   w a v e   d r a g -------------------
!-----------------------------------------------------------------------
    if (do_mg_drag) then

      call mg_drag(is, js, delt, u, v, t, pfull, phalf, zfull, zhalf, &
                   utnd, vtnd, ttnd, taubx, tauby, taus, kbot)
      udt = udt + utnd
      vdt = vdt + vtnd
      tdt = tdt + ttnd

!----- diagnostics -----

      if (id_udt_gwd > 0) then
        used = send_data(id_udt_gwd, utnd, Time, is, js, 1, &
                         rmask=mask)
      end if

      if (id_vdt_gwd > 0) then
        used = send_data(id_vdt_gwd, vtnd, Time, is, js, 1, &
                         rmask=mask)
      end if

      if (id_taubx > 0) then
        used = send_data(id_taubx, taubx, Time, is, js)
      end if

      if (id_tauby > 0) then
        used = send_data(id_tauby, tauby, Time, is, js)
      end if

      if (id_taus > 0) then
        used = send_data(id_taus, taus, Time, is, js, 1, &
                         rmask=mask)
      end if

      if (id_tdt_diss_gwd > 0) then
        used = send_data(id_tdt_diss_gwd, ttnd, Time, is, js, 1, &
                         rmask=mask)
      end if

      if (id_diss_heat_gwd > 0) then
        do k = 1, size(u, 3)
          pmass(:, :, k) = phalf(:, :, k + 1) - phalf(:, :, k)
        end do
        diag2 = cp_air/grav*sum(ttnd*pmass, 3)
        used = send_data(id_diss_heat_gwd, diag2, Time, is, js)
      end if

    end if

!   Alexander-Dunkerton gravity wave drag

    if (do_cg_drag) then
!mj updating call to riga version of cg_drag
      !call cg_drag_calc (is, js, lat, pfull, zfull, t, u, Time,    &
      !                  delt, utnd)
      call cg_drag_calc(is, js, lat, pfull, zfull, t, u, v, Time, delt, utnd, vtnd)
      udt = udt + utnd
      vdt = vdt + vtnd !mj

!----- diagnostics -----

      if (id_udt_cgwd > 0) then
        used = send_data(id_udt_cgwd, utnd, Time, is, js, 1, &
                         rmask=mask)
      end if

    end if

! constant drag, modeled on Alexander-Dunkerton winter average
    if (do_const_drag) then
      ! get time of the year for seasonal cycle
      call get_time(length_of_year(), seconds, daysperyear)
      call get_time(Time, seconds, days)
      cosday = cos(2*PI*days/daysperyear)
      utnd = 0.
      minp = log(minval(pfull*0.01)) - 1.
      where (pfull*0.01 < exp(1.))
        ! vertical: linear in ln(p)
        utnd = -const_drag_amp*((log(pfull*0.01) - 1.)/minp)**1.
      end where
      ! latitudinal: 3rd order polynomial, and cosine in time
      do k = 1, size(utnd, 3)
        where (pfull(:, :, k)*0.01 < exp(1.))
          utnd(:, :, k) = utnd(:, :, k)*sign(1., lat)*cosday &
                          *(-1.65*abs(lat)**3 + 2.5*lat**2 + 0.17*abs(lat) + const_drag_off)
        end where
      end do
      udt = udt + utnd

!----- diagnostics -----

      if (id_udt_cnstd > 0) then
        used = send_data(id_udt_cnstd, utnd, Time, is, js, 1, &
                         rmask=mask)
      end if
    end if

!-----------------------------------------------------------------------

  end subroutine damping_driver

!#######################################################################

  !> Initializes the module: reads `damping_driver_nml`, sets up the Rayleigh sponge,
  !> initializes `mg_drag_mod` and `cg_drag_mod` if they are used and registers the
  !> diagnostics.
  subroutine damping_driver_init(lonb, latb, domain, pref, axes, Time, sgsmtn)

    real, intent(in) :: lonb(:), latb(:), pref(:)
    !! `lonb`, `latb`: longitudes and latitudes of the grid box edges [rad]; `pref`: reference
    !! pressures at full levels (plus the surface value at nlev+1) [Pa]
    type(domain2d), intent(in) :: domain  !! domain decomposition of the model grid
    integer, intent(in) :: axes(4)  !! diagnostic axes (lon, lat, pfull, phalf)
    type(time_type), intent(in) :: Time  !! current time
    real, dimension(:, :), intent(out) :: sgsmtn
    !! sub-grid scale topography variance (from `mg_drag_init`; only set if `do_mg_drag`) [m]
    integer :: unit, ierr, io
    logical :: used
!mj
    integer :: raylev(1)
!-----------------------------------------------------------------------
!----------------- namelist (read & write) -----------------------------

    read (input_nml_file, nml=damping_driver_nml, iostat=io)
    ierr = check_nml_error(io, 'damping_driver_nml')

    call write_version_number(version, tagname)
    if (mpp_pe() == mpp_root_pe()) then
      write (stdlog(), nml=damping_driver_nml)
    end if

!-----------------------------------------------------------------------
!--------- rayleigh friction ----------

!mj actively choose rayleigh friction
!   do_rayleigh=.false.

!   if (abs(trayfric) > 0.0001 .and. nlev_rayfric > 0) then
    if (do_rayleigh) then
! mj automatically determine nlev_rayfric
      raylev = minloc(abs(pref(:) - 2*sponge_pbottom))
      nlev_rayfric = raylev(1)
      if (trayfric > 0.0) then
        rfactr = (1./trayfric)
      else
        rfactr = (1./abs(trayfric))*daypsec
      end if
    end if
!      do_rayleigh=.true.
!   else
!      rfactr=0.0
!   endif

!-----------------------------------------------------------------------
!----- mountain gravity wave drag -----

    if (do_mg_drag) call mg_drag_init(lonb, latb, domain, sgsmtn)

!--------------------------------------------------------------------
!----- Alexander-Dunkerton gravity wave drag -----

    if (do_cg_drag) then
      call cg_drag_init(lonb, latb, domain, pref, Time=Time, axes=axes)
    end if

!-----------------------------------------------------------------------
!----- initialize diagnostic fields -----

    if (do_rayleigh) then

      id_udt_rdamp = &
        register_diag_field(mod_name, 'udt_rdamp', axes(1:3), Time, &
                            'u wind tendency for Rayleigh damping', 'm/s2', &
                            missing_value=missing_value)

      id_vdt_rdamp = &
        register_diag_field(mod_name, 'vdt_rdamp', axes(1:3), Time, &
                            'v wind tendency for Rayleigh damping', 'm/s2', &
                            missing_value=missing_value)

      id_tdt_diss_rdamp = &
        register_diag_field(mod_name, 'tdt_diss_rdamp', axes(1:3), Time, &
                            'Dissipative heating from Rayleigh damping', &
                            'K/s', missing_value=missing_value)

      id_diss_heat_rdamp = &
        register_diag_field(mod_name, 'diss_heat_rdamp', axes(1:2), Time, &
                            'Integrated dissipative heating from Rayleigh damping', &
                            'W/m2')
    end if

    if (do_mg_drag) then

      ! register and send static field
      id_sgsmtn = &
        register_static_field(mod_name, 'sgsmtn', axes(1:2), &
                              'sub-grid scale topography for gravity wave drag', 'm')
      if (id_sgsmtn > 0) used = send_data(id_sgsmtn, sgsmtn)

      ! register non-static field
      id_udt_gwd = &
        register_diag_field(mod_name, 'udt_gwd', axes(1:3), Time, &
                            'u wind tendency for gravity wave drag', 'm/s2', &
                            missing_value=missing_value)

      id_vdt_gwd = &
        register_diag_field(mod_name, 'vdt_gwd', axes(1:3), Time, &
                            'v wind tendency for gravity wave drag', 'm/s2', &
                            missing_value=missing_value)

      id_taubx = &
        register_diag_field(mod_name, 'taubx', axes(1:2), Time, &
                            'x base flux for grav wave drag', 'N/m2', &
                            missing_value=missing_value)

      id_tauby = &
        register_diag_field(mod_name, 'tauby', axes(1:2), Time, &
                            'y base flux for grav wave drag', 'N/m2', &
                            missing_value=missing_value)

      id_taus = &
        register_diag_field(mod_name, 'taus', axes(1:3), Time, &
                            'saturation flux for gravity wave drag', 'N/m2', &
                            missing_value=missing_value)

      id_tdt_diss_gwd = &
        register_diag_field(mod_name, 'tdt_diss_gwd', axes(1:3), Time, &
                            'Dissipative heating from gravity wave drag', &
                            'K/s', missing_value=missing_value)

      id_diss_heat_gwd = &
        register_diag_field(mod_name, 'diss_heat_gwd', axes(1:2), Time, &
                            'Integrated dissipative heating from gravity wave drag', &
                            'W/m2')
    end if

    if (do_cg_drag) then

      id_udt_cgwd = &
        register_diag_field(mod_name, 'udt_cgwd', axes(1:3), Time, &
                            'u wind tendency for cg gravity wave drag', 'm/s2', &
                            missing_value=missing_value)
    end if

    if (do_const_drag) then

      id_udt_cnstd = &
        register_diag_field(mod_name, 'udt_cnstd', axes(1:3), Time, &
                            'u wind tendency for constant drag', 'm/s2', &
                            missing_value=missing_value)
    end if

!-----------------------------------------------------------------------

    module_is_initialized = .true.

!******************** end of initialization ****************************
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

  end subroutine damping_driver_init

!#######################################################################

  !> Finalizes `mg_drag_mod` and `cg_drag_mod` (which write their restart files).
  subroutine damping_driver_end

    if (do_mg_drag) call mg_drag_end
    if (do_cg_drag) call cg_drag_end

    module_is_initialized = .false.

  end subroutine damping_driver_end

!#######################################################################

  !> Computes the Rayleigh sponge tendencies of u and v and, if `do_conserve_energy`, the
  !> heating by the dissipated kinetic energy.
  subroutine rayleigh(dt, pres, u, v, udt, vdt, tdt)

    real, intent(in)                      :: dt
    real, intent(in), dimension(:, :, :)   :: pres, u, v
    real, intent(out), dimension(:, :, :)   :: udt, vdt, tdt

    real, dimension(size(u, 1), size(u, 2)) :: fact
    integer :: k
!-----------------------------------------------------------------------
!--------------rayleigh damping of momentum (to zero)-------------------

    udt = 0.
    vdt = 0.
    do k = 1, nlev_rayfric
      where (pres(:, :, k) < sponge_pbottom)
        fact(:, :) = rfactr*(sponge_pbottom - pres(:, :, k))**2/(sponge_pbottom)**2
        udt(:, :, k) = -u(:, :, k)*fact(:, :)
        vdt(:, :, k) = -v(:, :, k)*fact(:, :)
      end where
    end do
!   do k = nlev_rayfric+1, size(u,3)
!     udt(:,:,k) = 0.0
!     vdt(:,:,k) = 0.0
!   enddo

!  total energy conservation
!  compute temperature change loss due to ke dissipation

    tdt = 0. !mj
    if (do_conserve_energy) then
      do k = 1, nlev_rayfric
        tdt(:, :, k) = -((u(:, :, k) + .5*dt*udt(:, :, k))*udt(:, :, k) + &
                         (v(:, :, k) + .5*dt*vdt(:, :, k))*vdt(:, :, k))/cp_air
      end do
!       do k = nlev_rayfric+1, size(u,3)
!          tdt(:,:,k) = 0.0
!       enddo
!   else
!       tdt = 0.0
    end if
!-----------------------------------------------------------------------

  end subroutine rayleigh

!#######################################################################

end module damping_driver_mod
