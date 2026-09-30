!> Spectral atmosphere: the interface between the atmosphere driver (`atmos_model_mod`)
!> and the spectral dynamical core and physics.
!>
!> Holds the grid-point state at two time levels for the leapfrog scheme. Each time step
!> is split in two calls, which correspond to the down and up sweeps of the implicit
!> vertical diffusion: `atmosphere_down` calls the first part of the physics (radiation,
!> vertical diffusion down to the surface), and `atmosphere_up` finishes the physics
!> (vertical diffusion up from the surface, moist processes), then steps the dynamics.
!> The state is saved in the restart file `RESTART/atmosphere.res.nc`.
module atmosphere_mod

  use mpp_mod, only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, MPP_CLOCK_SYNC

  use fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL, WARNING, write_version_number

  use restart_file_mod, only: restart_file_type, open_restart_read, open_restart_write, close_restart, &
                              read_restart_field, write_restart_field, get_restart_field_size

  use field_manager_mod, only: MODEL_ATMOS

  use tracer_manager_mod, only: get_number_tracers, get_tracer_index

  use spectral_physics_mod, only: spectral_physics_init, spectral_physics_down, spectral_physics_up, &
                                  spectral_physics_end, surf_diff_type

  use constants_mod, only: grav

  use transforms_mod, only: trans_grid_to_spherical, trans_spherical_to_grid, get_deg_lon, get_deg_lat, &
                            get_wts_lat, get_grid_boundaries, compute_ucos_vcos, divide_by_cos, get_lon_max, &
                            get_lat_max, get_grid_domain

  use spec_mpp_mod, only: grid_domain

  use press_and_geopot_mod, only: compute_pressures_and_heights, compute_z_bot

  use time_manager_mod, only: time_type, set_time, get_time, operator(+), operator(<), operator(-)

  use spectral_dynamics_mod, only: spectral_dynamics_init, spectral_dynamics, spectral_dynamics_end, get_num_levels, &
                                   complete_robert_filter, &
                                   get_axis_id, spectral_diagnostics, get_initial_fields

  use mpp_domains_mod, only: domain2d

  use tracer_type_mod, only: tracer_type

  implicit none
  private

  character(len=128), parameter :: version = &
                                   '$Id: atmosphere.f90,v 12.0 2005/04/14 15:52:50 fms Exp $'

  character(len=128), parameter :: tagname = &
                                   '$Name: lima $'

  public :: atmosphere_init, atmosphere_down, atmosphere_up, atmosphere_end, atmosphere_domain
  public :: atmosphere_resolution, atmosphere_boundary, get_bottom_mass, get_bottom_wind, get_atmosphere_axes
  public :: surf_diff_type
  integer :: seconds, days, num_tracers, num_levels, nhum

  integer, parameter :: num_time_levels = 2
  integer :: phyclock, dynclock

  real, allocatable, dimension(:, :, :) :: p_half, p_full, z_half, z_full, wg_full
  type(tracer_type), allocatable, dimension(:) :: tracer_attributes
  real, allocatable, dimension(:, :, :, :, :) :: grid_tracers
  real, allocatable, dimension(:, :, :) :: psg
  real, allocatable, dimension(:, :, :, :) :: ug, vg, tg

  real, allocatable, dimension(:, :) :: dt_psg
  real, allocatable, dimension(:, :, :) :: dt_ug, dt_vg, dt_tg
  real, allocatable, dimension(:, :, :, :) :: dt_tracers

! dt_real is the atmospheric time step, converted to a real number.
! delta_t is passed to physics. It is twice dt_real, except for
! the first time step of a cold start, when it equals dt_real.
  real :: delta_t, dt_real

  integer :: is, ie, js, je
  integer :: previous, current, future
  logical :: module_is_initialized = .false., atmos_domain_is_computed = .false.

  type(time_type) :: Time_step, Time_prev, Time_next

!------------------------------------------------------------------------------------------------

contains

!####################################################################################################################

  !> Initializes the dynamics and the physics, and reads the state from
  !> `INPUT/atmosphere.res.nc` if it exists (otherwise the cold-start fields of the dynamics
  !> are used).
  subroutine atmosphere_init(Time_init, Time, Time_step_in, Surf_diff)

    type(time_type), intent(in)    :: Time_init, Time, Time_step_in
    !! `Time_init`: initial time of the experiment (not used); `Time`: current time;
    !! `Time_step_in`: atmospheric time step
    type(surf_diff_type), intent(inout) :: Surf_diff  !! surface data of the implicit vertical diffusion

    integer :: j, k, time_level, lon_max, lat_max, ntr, nt
    integer, dimension(4) :: siz
    real :: level
    character(len=64) :: tr_name
    character(len=4) :: ch1, ch2, ch3, ch4, ch5, ch6
    type(restart_file_type) :: rst

    if (module_is_initialized) return

    dynclock = mpp_clock_id('Dynamics', flags=MPP_CLOCK_SYNC)
    phyclock = mpp_clock_id('Physics', flags=MPP_CLOCK_SYNC)

    call write_version_number(version, tagname)
!-----------------------------------------------------------------------------------------

!  because the time step is used in different ways,
!  it must exist as a real and time_type variable

    call get_time(Time_step_in, seconds, days)
    dt_real = float(86400*days + seconds)
    Time_step = Time_step_in

    call get_number_tracers(MODEL_ATMOS, num_prog=num_tracers)
    allocate (tracer_attributes(num_tracers))

    call spectral_dynamics_init(Time, Time_step, tracer_attributes, nhum)
    atmos_domain_is_computed = .true.
    call get_grid_domain(is, ie, js, je)
    call get_num_levels(num_levels)

    allocate (p_half(is:ie, js:je, num_levels + 1))
    allocate (z_half(is:ie, js:je, num_levels + 1))
    allocate (p_full(is:ie, js:je, num_levels))
    allocate (z_full(is:ie, js:je, num_levels))
    allocate (wg_full(is:ie, js:je, num_levels))
    allocate (psg(is:ie, js:je, num_time_levels))
    allocate (ug(is:ie, js:je, num_levels, num_time_levels))
    allocate (vg(is:ie, js:je, num_levels, num_time_levels))
    allocate (tg(is:ie, js:je, num_levels, num_time_levels))
    allocate (grid_tracers(is:ie, js:je, num_levels, num_time_levels, num_tracers))
    allocate (dt_psg(is:ie, js:je))
    allocate (dt_ug(is:ie, js:je, num_levels))
    allocate (dt_vg(is:ie, js:je, num_levels))
    allocate (dt_tg(is:ie, js:je, num_levels))
    allocate (dt_tracers(is:ie, js:je, num_levels, num_tracers))

    p_half = 0.; z_half = 0.; p_full = 0.; z_full = 0.; wg_full = 0.
    psg = 0.; ug = 0.; vg = 0.; tg = 0.; grid_tracers = 0.
    dt_psg = 0.; dt_ug = 0.; dt_vg = 0.; dt_tg = 0.; dt_tracers = 0.

    if (open_restart_read(rst, 'INPUT/atmosphere.res.nc', grid_domain)) then
      call get_lon_max(lon_max)
      call get_lat_max(lat_max)
      call get_restart_field_size(rst, 'ug', siz)
      if (lon_max /= siz(1) .or. lat_max /= siz(2)) then
        write (ch1, '(i4)') siz(1)
        write (ch2, '(i4)') siz(2)
        write (ch3, '(i4)') lon_max
        write (ch4, '(i4)') lat_max
        call error_mesg('atmosphere_init', 'Resolution of restart data does not match resolution specified on namelist.'// &
                        ' Restart data: lon_max='//ch1//', lat_max='//ch2//'  Namelist: lon_max='//ch3//', lat_max='//ch4, FATAL)
      end if
      call read_restart_field(rst, 'previous', level)
      previous = nint(level)
      call read_restart_field(rst, 'current', level)
      current = nint(level)
      do nt = 1, num_time_levels
        call read_restart_field(rst, 'ug', ug(:, :, :, nt), nt)
        call read_restart_field(rst, 'vg', vg(:, :, :, nt), nt)
        call read_restart_field(rst, 'tg', tg(:, :, :, nt), nt)
        call read_restart_field(rst, 'psg', psg(:, :, nt), nt)
        do ntr = 1, num_tracers
          tr_name = trim(tracer_attributes(ntr)%name)
          call read_restart_field(rst, trim(tr_name), grid_tracers(:, :, :, nt, ntr), nt)
        end do ! end loop over tracers
      end do ! end loop over time levels
      call read_restart_field(rst, 'wg_full', wg_full)
      call close_restart(rst)
    else
      previous = 1; current = 1
      call get_initial_fields(ug(:, :, :, 1), vg(:, :, :, 1), tg(:, :, :, 1), psg(:, :, 1), grid_tracers(:, :, :, 1, :))
    end if

    call spectral_physics_init(Time, get_axis_id(), Surf_diff, nhum, p_half)

    call compute_pressures_and_heights( &
      tg(:, :, :, current), psg(:, :, current), z_full, z_half, p_full, p_half, grid_tracers(:, :, :, current, nhum))

    module_is_initialized = .true.

    return
  end subroutine atmosphere_init
!#################################################################################################################################

  !> Computes the first part of the physics tendencies (`physics_driver_down`): radiation,
  !> and vertical diffusion down to the surface.
  !>
  !> The physics time step is twice the model time step (leapfrog), except on the first
  !> step of a cold start.
  subroutine atmosphere_down(Time, frac_land, t_surf, albedo, &
                             rough_mom, u_star, b_star, q_star, dtau_du, dtau_dv, tau_x, tau_y, &
                             gust, flux_sw, flux_lw, Surf_diff)

    type(time_type), intent(in) :: Time  !! current time
    real, intent(in), dimension(:, :) :: frac_land, t_surf, albedo
    !! `frac_land`: land fraction; `t_surf`: surface temperature [K]; `albedo`: surface albedo
    real, intent(in), dimension(:, :) :: rough_mom, u_star, b_star, q_star, dtau_du, dtau_dv
    !! `rough_mom`: roughness length for momentum [m]; `u_star`: friction velocity [m/s];
    !! `b_star`: buoyancy scale [m/s2]; `q_star`: moisture scale [kg/kg]; `dtau_du`, `dtau_dv`:
    !! derivatives of the zonal and meridional surface stress with respect to the
    !! lowest-level wind [kg/m2/s]
    real, intent(inout), dimension(:, :) :: tau_x, tau_y  !! zonal and meridional surface stress [Pa]
    real, intent(out), dimension(:, :) :: flux_sw, flux_lw, gust
    !! `flux_sw`: net downward shortwave flux at the surface [W/m2]; `flux_lw`: downward
    !! longwave flux at the surface [W/m2]; `gust`: gustiness [m/s]
    type(surf_diff_type), intent(inout)                 :: Surf_diff  !! surface data of the implicit vertical diffusion

    integer :: days, seconds

    if (.not. module_is_initialized) then
      call error_mesg('atmosphere_down', 'atmosphere module has not been initialized.', FATAL)
    end if

    dt_psg = 0.
    dt_ug = 0.
    dt_vg = 0.
    dt_tg = 0.
    dt_tracers = 0.

    if (current == previous) then
      delta_t = dt_real
      Time_prev = Time
    else
      delta_t = 2*dt_real
      Time_prev = Time - Time_step
    end if
    Time_next = Time + Time_step

    call mpp_clock_begin(phyclock)
    call spectral_physics_down(Time_prev, Time, Time_next, previous, current, p_half, p_full, z_half, z_full, psg, &
                               ug, vg, tg, grid_tracers, frac_land, rough_mom, albedo, t_surf, u_star, b_star, q_star, &
                               dtau_du, dtau_dv, tau_x, tau_y, &
                               dt_ug, dt_vg, dt_tg, dt_tracers, flux_sw, flux_lw, gust, Surf_diff)
    call mpp_clock_end(phyclock)

    return
  end subroutine atmosphere_down
!#################################################################################################################################

  !> Finishes the physics (`physics_driver_up`: vertical diffusion up from the surface and
  !> moist processes), steps the dynamics, completes the Robert filter and sends the
  !> dynamics diagnostics.
  subroutine atmosphere_up(Time, frac_land, Surf_diff, lprec, fprec, gust)

    type(time_type), intent(in)                          :: Time  !! current time
    real, intent(in), dimension(is:ie, js:je) :: frac_land  !! land fraction
    type(surf_diff_type), intent(inout)                       :: Surf_diff
    !! surface data of the implicit vertical diffusion, including the changes of the
    !! lowest-level temperature and humidity from the surface
    real, intent(out), dimension(is:ie, js:je) :: lprec, fprec, gust
    !! `lprec`, `fprec`: liquid and frozen precipitation rate [kg/m2/s]; `gust`: gustiness [m/s]

    if (.not. module_is_initialized) then
      call error_mesg('atmosphere_up', 'atmosphere module has not been initialized.', FATAL)
    end if

    call mpp_clock_begin(phyclock)
    call spectral_physics_up(Time_prev, Time, Time_next, previous, current, p_half, p_full, z_half, z_full, wg_full, ug, vg, tg, &
                             grid_tracers, frac_land, dt_ug, dt_vg, dt_tg, dt_tracers, Surf_diff, lprec, fprec, gust)
    call mpp_clock_end(phyclock)

    if (previous == current) then
      future = num_time_levels + 1 - current
    else
      future = previous
    end if

    call mpp_clock_begin(dynclock)
    call spectral_dynamics(Time, psg(:, :, future), ug(:, :, :, future), vg(:, :, :, future), &
                           tg(:, :, :, future), tracer_attributes, grid_tracers(:, :, :, future, :), &
                           dt_psg, dt_ug, dt_vg, dt_tg, dt_tracers, wg_full, p_full, p_half, z_full)
    call mpp_clock_end(dynclock)

    call complete_robert_filter(tracer_attributes)

    call spectral_diagnostics(Time_next, psg(:, :, future), ug(:, :, :, future), vg(:, :, :, future), &
                              tg(:, :, :, future), wg_full, grid_tracers(:, :, :, future, :))

    previous = current
    current = future

    call compute_pressures_and_heights( &
      tg(:, :, :, current), psg(:, :, current), z_full, z_half, p_full, p_half, grid_tracers(:, :, :, current, nhum))

    return
  end subroutine atmosphere_up
!####################################################################################################################

  !> Returns the temperature, humidity, pressure and height of the lowest model level, and
  !> the surface pressure, at the previous time level (the level the implicit vertical
  !> diffusion steps from).
  subroutine get_bottom_mass(t_bot, q_bot, p_bot, z_bot_out, p_surf)

    real, intent(out), dimension(:, :) :: t_bot, q_bot, p_bot, z_bot_out, p_surf
    !! `t_bot`: temperature [K]; `q_bot`: specific humidity [kg/kg]; `p_bot`: pressure [Pa];
    !! `z_bot_out`: height above the surface [m], all at the lowest level; `p_surf`: surface
    !! pressure [Pa]

    real, dimension(size(t_bot, 1), size(t_bot, 2), num_levels) :: p_full_prev, z_full_prev
    real, dimension(size(t_bot, 1), size(t_bot, 2), num_levels + 1) :: p_half_prev, z_half_prev

    if (.not. module_is_initialized) then
      call error_mesg('get_bottom_mass', 'atmosphere module has not been initialized.', FATAL)
    end if

! All bottom-level fields are taken at the previous time level, which is
! the level the implicit vertical diffusion steps from. The pressure and
! height of the lowest level are recomputed for that level here (p_full
! holds the values for the current level).
    t_bot = tg(:, :, num_levels, previous)
    q_bot = grid_tracers(:, :, num_levels, previous, nhum)
    p_surf = psg(:, :, previous)
    call compute_pressures_and_heights(tg(:, :, :, previous), psg(:, :, previous), &
                                       z_full_prev, z_half_prev, p_full_prev, p_half_prev, &
                                       grid_tracers(:, :, :, previous, nhum))
    p_bot = p_full_prev(:, :, num_levels)
    call compute_z_bot(psg(:, :, previous), tg(:, :, num_levels, previous), z_bot_out, &
                       grid_tracers(:, :, num_levels, previous, nhum))

    return
  end subroutine get_bottom_mass
!####################################################################################################################

  !> Returns the wind at the lowest model level, at the previous time level.
  subroutine get_bottom_wind(u_bot, v_bot)

    real, intent(out), dimension(:, :) :: u_bot, v_bot  !! zonal and meridional wind [m/s]

    if (.not. module_is_initialized) then
      call error_mesg('get_bottom_wind', 'atmosphere module has not been initialized.', FATAL)
    end if

    u_bot = ug(:, :, num_levels, previous)
    v_bot = vg(:, :, num_levels, previous)

    return
  end subroutine get_bottom_wind
!####################################################################################################################

  !> Returns the number of longitudes and latitudes of the global grid or of this PE's
  !> subdomain.
  subroutine atmosphere_resolution(num_lon_out, num_lat_out, global)

    integer, intent(out)          :: num_lon_out, num_lat_out  !! number of longitudes and latitudes
    logical, intent(in), optional :: global  !! `.true.`: global grid; `.false.` (default): this PE's subdomain
    logical :: global_tmp

    if (.not. module_is_initialized) then
      call error_mesg('atmosphere_resolution', 'atmosphere module has not been initialized.', FATAL)
    end if

    if (present(global)) then
      global_tmp = global
    else
      global_tmp = .false.
    end if

    if (global_tmp) then
      call get_lon_max(num_lon_out)
      call get_lat_max(num_lat_out)
    else
      num_lon_out = ie - is + 1
      num_lat_out = je - js + 1
    end if

    return
  end subroutine atmosphere_resolution
!####################################################################################################################

  !> Returns the diagnostic axis ids of the atmosphere grid.
  subroutine get_atmosphere_axes(axes_out)
    integer, intent(out), dimension(:) :: axes_out  !! axis ids (lon, lat, pfull, phalf)

    if (.not. module_is_initialized) then
      call error_mesg('get_atmosphere_axes', 'atmosphere module has not been initialized.', FATAL)
    end if

    axes_out = get_axis_id()

    return
  end subroutine get_atmosphere_axes
!####################################################################################################################

  !> Returns the longitudes and latitudes of the grid-box boundaries.
  subroutine atmosphere_boundary(lon_boundaries, lat_boundaries, global)

    real, intent(out), dimension(:) :: lon_boundaries, lat_boundaries
    !! longitudes and latitudes of the grid-box boundaries [rad]
    logical, intent(in), optional     :: global  !! `.true.`: global grid; `.false.` (default): this PE's subdomain

    logical :: global_tmp

    if (.not. module_is_initialized) then
      call error_mesg('atmosphere_boundary', 'atmosphere module has not been initialized.', FATAL)
    end if

    if (present(global)) then
      global_tmp = global
    else
      global_tmp = .false.
    end if
    call get_grid_boundaries(lon_boundaries, lat_boundaries, global_tmp)

    return
  end subroutine atmosphere_boundary
!####################################################################################################################
  !> Returns the domain decomposition of the atmosphere grid.
  subroutine atmosphere_domain(domain)
    type(domain2d), intent(out) :: domain  !! domain decomposition of the grid

    if (.not. atmos_domain_is_computed) then
      call error_mesg('atmosphere_domain', 'spec_mpp has not been initialized.', FATAL)
    end if

    domain = grid_domain

  end subroutine atmosphere_domain
!####################################################################################################################

  !> Writes the restart file `RESTART/atmosphere.res.nc` and terminates the physics and the
  !> dynamics.
  subroutine atmosphere_end(Time)
    type(time_type), intent(in) :: Time  !! current time
    integer :: ntr, nt
    character(len=64) :: tr_name
    type(restart_file_type) :: rst

    if (.not. module_is_initialized) return

    call open_restart_write(rst, 'RESTART/atmosphere.res.nc', grid_domain)
    call write_restart_field(rst, 'previous', real(previous))
    call write_restart_field(rst, 'current', real(current))
    do nt = 1, num_time_levels
      call write_restart_field(rst, 'ug', ug(:, :, :, nt), nt)
      call write_restart_field(rst, 'vg', vg(:, :, :, nt), nt)
      call write_restart_field(rst, 'tg', tg(:, :, :, nt), nt)
      call write_restart_field(rst, 'psg', psg(:, :, nt), nt)
      do ntr = 1, num_tracers
        tr_name = trim(tracer_attributes(ntr)%name)
        call write_restart_field(rst, trim(tr_name), grid_tracers(:, :, :, nt, ntr), nt)
      end do
    end do
    call write_restart_field(rst, 'wg_full', wg_full)
    call close_restart(rst)

    deallocate (dt_psg, dt_ug, dt_vg, dt_tg, dt_tracers)

    call spectral_physics_end(Time)
    call spectral_dynamics_end(tracer_attributes, Time)

    module_is_initialized = .false.

    return
  end subroutine atmosphere_end
!####################################################################################################################

end module atmosphere_mod
