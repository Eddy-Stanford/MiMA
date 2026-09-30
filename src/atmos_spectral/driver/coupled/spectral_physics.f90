!> Interface between the spectral atmosphere (`atmosphere_mod`) and the physics driver
!> (`physics_driver_mod`).
!>
!> Passes the grid-point fields at the current and previous time levels of the leapfrog
!> scheme, and the latitude, longitude and area of each grid point, to the down and up
!> parts of the physics. Diagnostic tracers, and tracers initialized by the physics, are
!> not supported.
module spectral_physics_mod

  use fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL, write_version_number, fms_init

  use constants_mod, only: grav, pi

  use time_manager_mod, only: time_type, set_time, get_time, operator(-), operator(/=), time_manager_init

  use press_and_geopot_mod, only: pressure_variables, compute_pressures_and_heights

  use transforms_mod, only: get_grid_boundaries, get_deg_lon, get_deg_lat, get_wts_lat, &
                            get_grid_domain, get_lon_max, get_lat_max

  use spec_mpp_mod, only: grid_domain

  use spectral_dynamics_mod, only: get_reference_sea_level_press, get_num_levels

  use physics_driver_mod, only: physics_driver_init, physics_driver_down, physics_driver_up, physics_driver_end, &
                                surf_diff_type

  use tracer_type_mod, only: tracer_type

  use field_manager_mod, only: MODEL_ATMOS

  use tracer_manager_mod, only: get_number_tracers

  implicit none
  private

  public :: spectral_physics_init, spectral_physics_down, spectral_physics_up, &
            spectral_physics_end, surf_diff_type

  character(len=128), parameter :: version = &
                                   '$Id: spectral_physics.f90,v 12.0 2005/04/14 15:53:05 fms Exp $'

  character(len=128), parameter :: tagname = &
                                   '$Name: lima $'

  integer, parameter :: num_time_levels = 2

  real, allocatable, dimension(:, :) :: rad_lon_2d, rad_lat_2d, area_2d
  real, allocatable, dimension(:, :, :, :) :: diag_tracers
  integer :: num_levels, num_tracers, nhum
  integer :: is, ie, js, je
  logical :: module_is_initialized = .false.

contains

!------------------------------------------------------------------------------------------------

  !> Sets up the grid-point latitudes, longitudes and areas and the reference pressure
  !> profiles (for surface pressures of 1013.25 and 810.6 hPa), and initializes the physics
  !> driver.
  subroutine spectral_physics_init(Time, axes, Surf_diff, nhum_in, p_half)

    type(time_type), intent(in) :: Time  !! current time
    integer, intent(in), dimension(:) :: axes  !! diagnostic axes (lon, lat, pfull, phalf)
    type(surf_diff_type), intent(inout) :: Surf_diff  !! surface data of the implicit vertical diffusion
    integer, intent(in) :: nhum_in  !! tracer index of the humidity
    real, intent(in), dimension(:, :, :) :: p_half  !! pressure at half levels [Pa]
    real, allocatable, dimension(:, :, :, :) :: grid_tracers

    real, allocatable, dimension(:) :: rad_lon, rad_lat, wts_lat, lon_boundaries, lat_boundaries
    real, dimension(2) :: radiation_ref_press_surf = (/101325., 81060./)

    real, allocatable, dimension(:, :) :: radiation_ref_press
    real, allocatable, dimension(:)   :: p_half_1d, ln_p_half_1d, shalf
    real, allocatable, dimension(:)   :: p_full_1d, ln_p_full_1d, sfull

    real :: reference_sea_level_press
    integer :: i, j, tracer_number, unit, num_diag, ntr, nmix_rat, nsphum, lon_max, lat_max

    if (module_is_initialized) return

    call write_version_number(version, tagname)

    call fms_init
    call time_manager_init

    nhum = nhum_in

    call get_grid_domain(is, ie, js, je)
    allocate (rad_lon(is:ie), rad_lat(js:je), wts_lat(js:je))
    allocate (lon_boundaries(ie - is + 2), lat_boundaries(je - js + 2))

    call get_num_levels(num_levels)
    allocate (radiation_ref_press(num_levels + 1, 2))
    allocate (p_half_1d(num_levels + 1), ln_p_half_1d(num_levels + 1), shalf(num_levels + 1))
    allocate (p_full_1d(num_levels), ln_p_full_1d(num_levels), sfull(num_levels))

    allocate (rad_lon_2d(is:ie, js:je))
    allocate (rad_lat_2d(is:ie, js:je))
    allocate (area_2d(is:ie, js:je))

    call pressure_variables(p_half_1d, ln_p_half_1d, radiation_ref_press(1:num_levels, 1), ln_p_full_1d, &
                            radiation_ref_press_surf(1))
    call pressure_variables(p_half_1d, ln_p_half_1d, radiation_ref_press(1:num_levels, 2), ln_p_full_1d, &
                            radiation_ref_press_surf(2))
    radiation_ref_press(num_levels + 1, :) = radiation_ref_press_surf

    call get_reference_sea_level_press(reference_sea_level_press)
    call pressure_variables(p_half_1d, ln_p_half_1d, p_full_1d, ln_p_full_1d, reference_sea_level_press)
    shalf = p_half_1d/reference_sea_level_press
    sfull = p_full_1d/reference_sea_level_press

    call get_lon_max(lon_max)
    call get_lat_max(lat_max)
    call get_deg_lon(rad_lon)
    call get_deg_lat(rad_lat)
    call get_wts_lat(wts_lat)
    rad_lon = pi*rad_lon/180.
    rad_lat = pi*rad_lat/180.
    do j = js, je
      rad_lat_2d(:, j) = rad_lat(j)
      area_2d(:, j) = wts_lat(j)/(2.*lon_max)
    end do
    do i = is, ie
      rad_lon_2d(i, :) = rad_lon(i)
    end do

    call get_grid_boundaries(lon_boundaries, lat_boundaries)

    call get_number_tracers(MODEL_ATMOS, num_diag=num_diag)
    if (num_diag > 0) then
      call error_mesg('spectral_physics_init', &
                      'This version of the spectral atmospheric model not coded to handle diagnostic tracers', FATAL)
    end if
    allocate (diag_tracers(is:ie, js:je, num_levels, num_diag))
    diag_tracers = 0.

    call get_number_tracers(MODEL_ATMOS, num_prog=num_tracers)

    allocate (grid_tracers(is:ie, js:je, num_levels, num_tracers))
    grid_tracers = 0.

    call physics_driver_init(Time, lon_boundaries, lat_boundaries, grid_domain, axes, radiation_ref_press, grid_tracers, &
                             Surf_diff, p_half)

    if (sum(grid_tracers) /= 0.) then
      call error_mesg('spectral_physics_init', 'This version of the spectral atmospheric model not coded to handle'// &
                      ' initialization of tracer fields by physics_driver_init', FATAL)
    end if

    deallocate (rad_lon, rad_lat, wts_lat, lon_boundaries, lat_boundaries, grid_tracers)
    deallocate (radiation_ref_press, p_half_1d, ln_p_half_1d, shalf, p_full_1d, ln_p_full_1d, sfull)

    module_is_initialized = .true.

    return
  end subroutine spectral_physics_init
!------------------------------------------------------------------------------------------------

  !> Calls `physics_driver_down` with the fields at the current and previous time levels.
  subroutine spectral_physics_down(Time_prev, Time, Time_next, previous, current, &
                                   p_half, p_full, z_half, z_full, psg, ug, vg, tg, grid_tracers, &
                                   frac_land, rough_mom, albedo, t_surf, u_star, b_star, q_star, dtau_du, dtau_dv, tau_x, tau_y, &
                                   dt_ug, dt_vg, dt_tg, dt_tracers, flux_sw, flux_lw, gust, Surf_diff)

    type(time_type), intent(in) :: Time_prev, Time, Time_next
    !! times of the previous, current and next time levels
    integer, intent(in)         :: previous, current  !! indices of the previous and current time levels
    real, intent(in), dimension(:, :, :) :: p_full, z_full
    !! pressure [Pa] and height [m] at full levels
    real, intent(in), dimension(:, :, :) :: p_half, z_half
    !! pressure [Pa] and height [m] at half levels
    real, intent(in), dimension(:, :, :) :: psg  !! surface pressure at the two time levels [Pa] (not used)
    real, intent(in), dimension(:, :, :, :) :: ug, vg, tg
    !! zonal and meridional wind [m/s] and temperature [K] at the two time levels
    real, intent(inout), dimension(:, :, :, :, :) :: grid_tracers  !! tracers at the two time levels
    real, intent(in), dimension(:, :) :: frac_land, rough_mom, albedo, t_surf, u_star, b_star, q_star, dtau_du, dtau_dv
    !! surface fields: land fraction; roughness length for momentum [m]; albedo; surface
    !! temperature [K]; friction velocity [m/s]; buoyancy scale [m/s2]; moisture scale
    !! [kg/kg]; derivatives of the zonal and meridional surface stress with respect to the
    !! lowest-level wind [kg/m2/s]
    real, intent(inout), dimension(:, :) :: tau_x, tau_y  !! zonal and meridional surface stress [Pa]
    real, intent(inout), dimension(:, :, :) :: dt_ug, dt_vg, dt_tg
    !! tendencies of the zonal and meridional wind [m/s2] and temperature [K/s]
    real, intent(inout), dimension(:, :, :, :) :: dt_tracers  !! tracer tendencies
    real, intent(out), dimension(:, :) :: flux_sw, flux_lw, gust
    !! `flux_sw`: net downward shortwave flux at the surface [W/m2]; `flux_lw`: downward
    !! longwave flux at the surface [W/m2]; `gust`: gustiness [m/s]
    type(surf_diff_type), intent(inout) :: Surf_diff  !! surface data of the implicit vertical diffusion

!**************************************************************************************

    if (.not. module_is_initialized) then
      call error_mesg('spectral_physics_down', 'spectral_physics module is not initialized', FATAL)
    end if

    call physics_driver_down(1, ie - is + 1, 1, je - js + 1, Time_prev, Time, Time_next, &
                             rad_lat_2d, rad_lon_2d, area_2d, &
                             p_half, p_full, z_half, z_full, &
                             ug(:, :, :, current), vg(:, :, :, current), &
                             tg(:, :, :, current), grid_tracers(:, :, :, current, nhum), &
                             grid_tracers(:, :, :, current, :), ug(:, :, :, previous), vg(:, :, :, previous), &
                             tg(:, :, :, previous), grid_tracers(:, :, :, previous, nhum), &
                             grid_tracers(:, :, :, previous, :), &
                             frac_land, rough_mom, albedo, &
                             t_surf, u_star, b_star, q_star, &
                             dtau_du, dtau_dv, tau_x, tau_y, &
                             dt_ug, dt_vg, dt_tg, &
                             dt_tracers(:, :, :, nhum), dt_tracers, flux_sw(:, :), &
                             flux_lw, gust, Surf_diff)

    return
  end subroutine spectral_physics_down
!------------------------------------------------------------------------------------------------
  !> Calls `physics_driver_up` with the fields at the current and previous time levels.
  subroutine spectral_physics_up(Time_prev, Time, Time_next, previous, current, p_half, p_full, &
                                 z_half, z_full, wg_full, ug, vg, tg, grid_tracers, &
                                 frac_land, dt_ug, dt_vg, dt_tg, dt_tracers, Surf_diff, lprec, fprec, gust)

    type(time_type), intent(in) :: Time_prev, Time, Time_next
    !! times of the previous, current and next time levels
    integer, intent(in) :: previous, current  !! indices of the previous and current time levels
    real, intent(in), dimension(:, :, :) :: p_full, z_full, wg_full
    !! pressure [Pa], height [m] and vertical pressure velocity [Pa/s] at full levels
    real, intent(in), dimension(:, :, :) :: p_half, z_half
    !! pressure [Pa] and height [m] at half levels
    real, intent(in), dimension(:, :, :, :) :: ug, vg, tg
    !! zonal and meridional wind [m/s] and temperature [K] at the two time levels
    real, intent(in), dimension(:, :, :, :, :) :: grid_tracers  !! tracers at the two time levels
    real, intent(in), dimension(:, :) :: frac_land  !! land fraction
    real, intent(inout), dimension(:, :, :) :: dt_ug, dt_vg, dt_tg
    !! tendencies of the zonal and meridional wind [m/s2] and temperature [K/s]
    real, intent(inout), dimension(:, :, :, :) :: dt_tracers  !! tracer tendencies
    type(surf_diff_type), intent(inout) :: Surf_diff  !! surface data of the implicit vertical diffusion
    real, intent(out), dimension(:, :) :: lprec, fprec
    !! liquid and frozen precipitation rate [kg/m2/s]
    real, intent(inout), dimension(:, :) :: gust  !! gustiness [m/s]

    if (.not. module_is_initialized) then
      call error_mesg('spectral_physics_up', 'spectral_physics module is not initialized', FATAL)
    end if

    call physics_driver_up(1, ie - is + 1, 1, je - js + 1, Time_prev, Time, Time_next, &
                           rad_lat_2d, rad_lon_2d, area_2d, &
                           p_half, p_full, z_half, z_full, &
                           wg_full, ug(:, :, :, current), vg(:, :, :, current), &
                           tg(:, :, :, current), grid_tracers(:, :, :, current, nhum), &
                           grid_tracers(:, :, :, current, :), &
                           ug(:, :, :, previous), vg(:, :, :, previous), &
                           tg(:, :, :, previous), grid_tracers(:, :, :, previous, nhum), &
                           grid_tracers(:, :, :, previous, :), &
                           frac_land, dt_ug, dt_vg, dt_tg, &
                           dt_tracers(:, :, :, nhum), dt_tracers, Surf_diff, &
                           lprec, fprec, gust)

    return
  end subroutine spectral_physics_up
!------------------------------------------------------------------------------------------------

  !> Terminates the physics driver.
  subroutine spectral_physics_end(Time)
    type(time_type), intent(in) :: Time  !! current time

    if (.not. module_is_initialized) return

    deallocate (rad_lon_2d, rad_lat_2d, area_2d)
    call physics_driver_end(Time)
    module_is_initialized = .false.

    return
  end subroutine spectral_physics_end
!------------------------------------------------------------------------------------------------

end module spectral_physics_mod
