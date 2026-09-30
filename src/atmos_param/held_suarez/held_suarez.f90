module held_suarez_mod

!-----------------------------------------------------------------------
!
!   Held and Suarez (1994) idealized forcing:
!     - Newtonian relaxation of temperature towards the zonally symmetric
!       equilibrium profile Teq(lat, p), and
!     - Rayleigh friction of the horizontal wind in the boundary layer.
!
!   Held, I. M. and M. J. Suarez, 1994: A proposal for the intercomparison
!   of the dynamical cores of atmospheric general circulation models.
!   Bull. Amer. Meteor. Soc., 75, 1825-1830.
!
!   The default parameters are those of HS94. The forcing depends only on
!   latitude and sigma = p/p_surf, so it applies at any resolution and
!   with any vertical levels.
!
!-----------------------------------------------------------------------

  use fms_mod, only: input_nml_file, check_nml_error, &
                     mpp_pe, mpp_root_pe, stdlog, &
                     error_mesg, FATAL, write_version_number
  use time_manager_mod, only: time_type
  use constants_mod, only: kappa, cp_air, seconds_per_day
  use diag_manager_mod, only: register_diag_field, send_data

  implicit none
  private

  public :: held_suarez_init, held_suarez_forcing, held_suarez_end

  character(len=128) :: version = '$Id: held_suarez.f90 $'
  character(len=128) :: tagname = '$Name: $'

!-------------------- namelist -----------------------------------------

  real    :: t_zero = 315.     ! surface equilibrium temperature at the equator [K]
  real    :: t_strat = 200.     ! minimum (stratospheric) equilibrium temperature [K]
  real    :: delh = 60.      ! equator-to-pole temperature difference [K]
  real    :: delv = 10.      ! vertical potential temperature difference [K]
  real    :: p_ref = 1.e5     ! reference pressure [Pa]
  real    :: sigma_b = 0.7      ! top of the frictional boundary layer [sigma]
  real    :: ka = 40.      ! free-atmosphere relaxation time [days]
  real    :: ks = 4.       ! surface relaxation time at the equator [days]
  real    :: kf = 1.       ! boundary-layer Rayleigh friction time [days]
  logical :: do_rayleigh_friction = .true.   ! apply the boundary-layer friction
  logical :: do_conserve_energy = .false.  ! heat the air by the frictional dissipation

  namelist /held_suarez_nml/ t_zero, t_strat, delh, delv, p_ref, sigma_b, ka, ks, kf, &
    do_rayleigh_friction, do_conserve_energy

!-------------------- diagnostics --------------------------------------

  integer :: id_teq, id_tdt_hs, id_udt_hs, id_vdt_hs, id_diss_heat_hs
  character(len=11), parameter :: mod_name = 'held_suarez'
  real :: missing_value = -999.

  logical :: module_is_initialized = .false.

contains

!#######################################################################

  subroutine held_suarez_init(axes, Time)

    integer, intent(in), dimension(4) :: axes
    type(time_type), intent(in)               :: Time

    integer :: unit, ierr, io

    if (module_is_initialized) return

    read (input_nml_file, nml=held_suarez_nml, iostat=io)
    ierr = check_nml_error(io, 'held_suarez_nml')

    call write_version_number(version, tagname)
    if (mpp_pe() == mpp_root_pe()) write (stdlog(), nml=held_suarez_nml)

    if (ka <= 0. .or. ks <= 0. .or. (do_rayleigh_friction .and. kf <= 0.)) &
      call error_mesg('held_suarez_init', 'ka, ks and kf must be positive (days)', FATAL)
    if (sigma_b <= 0. .or. sigma_b >= 1.) &
      call error_mesg('held_suarez_init', 'sigma_b must lie between 0 and 1', FATAL)

    id_teq = register_diag_field(mod_name, 'teq', axes(1:3), Time, &
                                 'Held-Suarez equilibrium temperature', 'K', missing_value=missing_value)
    id_tdt_hs = register_diag_field(mod_name, 'tdt_hs', axes(1:3), Time, &
                                    'Temperature tendency from Held-Suarez relaxation', 'K/s', missing_value=missing_value)
    id_udt_hs = register_diag_field(mod_name, 'udt_hs', axes(1:3), Time, &
                                    'Zonal wind tendency from Held-Suarez friction', 'm/s2', missing_value=missing_value)
    id_vdt_hs = register_diag_field(mod_name, 'vdt_hs', axes(1:3), Time, &
                                    'Meridional wind tendency from Held-Suarez friction', 'm/s2', missing_value=missing_value)
    id_diss_heat_hs = register_diag_field(mod_name, 'diss_heat_hs', axes(1:3), Time, &
                                          'Heating from Held-Suarez frictional dissipation', 'K/s', missing_value=missing_value)

    module_is_initialized = .true.

  end subroutine held_suarez_init

!#######################################################################

  subroutine held_suarez_forcing(is, js, Time, lat, p_full, p_half, u, v, t, udt, vdt, tdt)

!-----------------------------------------------------------------------
!   Adds the Held-Suarez tendencies to udt, vdt and tdt. u, v and t are
!   the fields the tendencies are computed from (the previous time level
!   in the leapfrog scheme), p_full and p_half the pressures [Pa].
!-----------------------------------------------------------------------

    integer, intent(in)                      :: is, js
    type(time_type), intent(in)                      :: Time
    real, intent(in), dimension(:, :)   :: lat
    real, intent(in), dimension(:, :, :) :: p_full, p_half, u, v, t
    real, intent(inout), dimension(:, :, :) :: udt, vdt, tdt

    real, dimension(size(t, 1), size(t, 2), size(t, 3)) :: teq, tdt_hs, udt_hs, vdt_hs, diss
    real, dimension(size(t, 1), size(t, 2))           :: sin2, cos2, cos4, p_surf
    real    :: sigma, sigma_fac, rka, rks, rkf
    integer :: k, nlev
    logical :: used

    if (.not. module_is_initialized) &
      call error_mesg('held_suarez_forcing', 'held_suarez_init has not been called', FATAL)

    nlev = size(t, 3)
    rka = 1./(ka*seconds_per_day)
    rks = 1./(ks*seconds_per_day)
    rkf = 1./(kf*seconds_per_day)
    sin2 = sin(lat)**2
    cos2 = 1.-sin2
    cos4 = cos2**2
    p_surf = p_half(:, :, nlev + 1)

    udt_hs = 0.
    vdt_hs = 0.
    do k = 1, nlev
      ! equilibrium temperature
      teq(:, :, k) = (t_zero - delh*sin2 - delv*log(p_full(:, :, k)/p_ref)*cos2) &
                     *(p_full(:, :, k)/p_ref)**kappa
      teq(:, :, k) = max(t_strat, teq(:, :, k))

      ! Newtonian relaxation, faster near the surface in the tropics
      tdt_hs(:, :, k) = -(rka + (rks - rka)*max(0., (p_full(:, :, k)/p_surf - sigma_b)/(1.-sigma_b))*cos4) &
                        *(t(:, :, k) - teq(:, :, k))

      ! Rayleigh friction in the boundary layer
      if (do_rayleigh_friction) then
        udt_hs(:, :, k) = -rkf*max(0., (p_full(:, :, k)/p_surf - sigma_b)/(1.-sigma_b))*u(:, :, k)
        vdt_hs(:, :, k) = -rkf*max(0., (p_full(:, :, k)/p_surf - sigma_b)/(1.-sigma_b))*v(:, :, k)
      end if
    end do

    if (do_conserve_energy) then
      diss = -(u*udt_hs + v*vdt_hs)/cp_air
    else
      diss = 0.
    end if

    tdt = tdt + tdt_hs + diss
    udt = udt + udt_hs
    vdt = vdt + vdt_hs

    if (id_teq > 0) used = send_data(id_teq, teq, Time, is, js, 1)
    if (id_tdt_hs > 0) used = send_data(id_tdt_hs, tdt_hs, Time, is, js, 1)
    if (id_udt_hs > 0) used = send_data(id_udt_hs, udt_hs, Time, is, js, 1)
    if (id_vdt_hs > 0) used = send_data(id_vdt_hs, vdt_hs, Time, is, js, 1)
    if (id_diss_heat_hs > 0) used = send_data(id_diss_heat_hs, diss, Time, is, js, 1)

  end subroutine held_suarez_forcing

!#######################################################################

  subroutine held_suarez_end

    module_is_initialized = .false.

  end subroutine held_suarez_end

!#######################################################################

end module held_suarez_mod
