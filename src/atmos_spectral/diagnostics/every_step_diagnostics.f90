!> Diagnostics of the dynamical core sent at every dynamics (sub)step (module name
!> `dynamics_every`).
!>
!> Besides the state (`ps_every`, `u_every`, ...) these include fluxes weighted by `ps/p00`,
!> static energies, the tendencies from the hyperdiffusion and the energy and water
!> corrections, the global advection tendencies of the tracers, eddy fluxes and two residual
!> mean streamfunctions (`psi_star`, `psi_dwc`; zonal means, independent of the
!> decomposition). The amplitude of the 2*dt computational
!> mode (the mean of `(-1)**n` times the field over all steps n) is sent as a static field
!> (`2dt_*`) at the end of the run.
module every_step_diagnostics_mod

  use fms_mod, only: error_mesg, FATAL

  use mpp_mod, only: mpp_sum

  use time_manager_mod, only: time_type

  use field_manager_mod, only: MODEL_ATMOS

  use tracer_manager_mod, only: get_number_tracers, get_tracer_names, get_tracer_index, NO_TRACER

  use diag_manager_mod, only: diag_axis_init, register_diag_field, register_static_field, send_data

  use press_and_geopot_mod, only: pressure_variables

  use transforms_mod, only: get_deg_lon, get_deg_lat, get_grid_domain
  use spec_mpp_mod, only: grid_domain

  use constants_mod, only: hlv, cp_air, grav, kappa, rvgas, rdgas, radius, radian, omega

!===============================================================================================
  implicit none
  private
!===============================================================================================

  public :: every_step_diagnostics_init, every_step_diagnostics, every_step_diagnostics_end

  ! set by every_step_diagnostics_init: whether the diagnostics need the damping tendencies of
  ! u and v or of T, the kegen fields, and the global advective tendencies of each tracer
  logical, public, protected :: need_damp_uv = .false., need_damp_t = .false., need_kegen = .false.
  logical, public, protected, allocatable, dimension(:) :: need_hadv, need_vadv

!===============================================================================================

  real, parameter :: p00 = 1.e5, virtual_factor = (rvgas/rdgas) - 1.0

  integer :: id_ps, id_u, id_v, id_t, num_levels, num_tracers

  integer :: id_drystaten, id_moiststaten, id_uv, id_vq, id_vdse, &
             id_vp, id_wu, id_wv, id_wq, id_wdse, id_wup, id_wvp, id_wqp, id_wdsep, id_wp, &
             id_wsubt, id_wsubtv, id_vqint, id_udt_damp, id_vdt_damp, id_tdt_damp, &
             id_entrop_dampuv, id_entrop_dampt, id_tdt_dampuv, id_tdt_tempcor, id_qdt_watercor, &
             id_entrop_tempcor, id_kegen, id_kegenq, id_kegenqtinv
  integer :: id_upvp, id_upwp, id_vptp, id_psi_star, id_psi_dwc

  integer, allocatable, dimension(:) :: id_tr, id_dt_hadv, id_dt_vadv

  logical :: module_is_initialized = .false.

  integer :: two_dt_id_ps, two_dt_id_u, two_dt_id_v, two_dt_id_t
  integer, allocatable, dimension(:) :: two_dt_id_tr
  integer :: iwt, num_time_steps

  real, allocatable, dimension(:, :)     :: two_dt_ps
  real, allocatable, dimension(:, :, :)   :: two_dt_u, two_dt_v, two_dt_t
  real, allocatable, dimension(:, :, :, :) :: two_dt_tr

  real, allocatable, dimension(:) :: lon, lat   ! all longitudes and latitudes [deg]

  integer :: is, ie, js, je
  ! tracer index of specific humidity. Initialized to NO_TRACER, then set in subroutine every_step_diagnostics_init.
  integer :: nsphum = NO_TRACER

  type(time_type) :: Time_save ! every_step_diagnostics_end needs time because it is a required argument of send_data,
  ! even though the fields are static. When Time is not available to every_step_diagnostics_end
  ! it uses Time_save instead.
!===============================================================================================

contains

!===============================================================================================

  !> Defines the `dynamics_every` axes and registers the diagnostics.
  subroutine every_step_diagnostics_init(Time, lon_max, lat_max, num_levels_in, reference_sea_level_press)

    type(time_type), intent(in) :: Time  !! current time
    integer, intent(in) :: lon_max, lat_max, num_levels_in
    !! `lon_max`, `lat_max`, `num_levels_in`: numbers of longitudes, latitudes and levels
    real, intent(in) :: reference_sea_level_press  !! surface pressure of the `pfull_every` axis values [Pa]

    integer, dimension(2) :: axes_2d
    integer, dimension(2) :: axes_zm
    integer, dimension(3) :: axes_3d
    integer :: id_lon, id_lat, id_pfull, ntr
    real, dimension(2) :: vrange, trange
    character(len=128) :: tname, longname, units
    character(len=16) :: mod_name = 'dynamics_every'
    real, dimension(num_levels_in)   :: p_full, ln_p_full
    real, dimension(num_levels_in + 1) :: p_half, ln_p_half

    allocate (lon(lon_max), lat(lat_max))
    call get_deg_lon(lon)
    call get_deg_lat(lat)
    id_lon = diag_axis_init('lon_every', lon, 'degrees_E', 'x', 'longitude', set_name=mod_name, Domain2=grid_domain)
    id_lat = diag_axis_init('lat_every', lat, 'degrees_N', 'y', 'latitude', set_name=mod_name, Domain2=grid_domain)
    call pressure_variables(p_half, ln_p_half, p_full, ln_p_full, reference_sea_level_press)
    p_full = .01*p_full
    id_pfull = diag_axis_init('pfull_every', p_full, 'hPa', 'z', 'approx full pressure level', direction=-1, set_name=mod_name)

    vrange = (/-400., 400./)
    trange = (/100., 400./)
    axes_2d = (/id_lon, id_lat/)
    axes_3d = (/id_lon, id_lat, id_pfull/)
    axes_zm = (/id_lat, id_pfull/)

    id_ps = register_diag_field(mod_name, 'ps_every', axes_2d, Time, 'surface pressure', 'Pa')
    id_u = register_diag_field(mod_name, 'u_every', axes_3d, Time, 'zonal wind component', 'm/s', range=vrange)
    id_v = register_diag_field(mod_name, 'v_every', axes_3d, Time, 'meridional wind component', 'm/s', range=vrange)
    id_t = register_diag_field(mod_name, 't_every', axes_3d, Time, 'temperature', 'K', range=trange)

    id_vp = register_diag_field(mod_name, 'vp', axes_3d, time, 'meridional wind, weighted by ps/p00', 'm/s', range=vrange)
    id_uv = register_diag_field(mod_name, 'uv', axes_3d, time, &
                                'zonal wind times meridional wind, weighted by ps/p00', 'm2/s2')
    id_vq = register_diag_field(mod_name, 'vq', axes_3d, time, &
                                'meridional wind times specific humidity, weighted by ps/p00', 'm/s kg/kg')
    id_vqint = register_diag_field(mod_name, 'vqint', (/id_lon, id_lat/), time, &
                                   'vertically integrated meridional moisture flux', 'kg/m/s')
    id_vdse = register_diag_field(mod_name, 'vdse', axes_3d, time, &
                                  'meridional wind times dry static energy, weighted by ps/p00', 'm/s J/kg')
    id_wu = register_diag_field(mod_name, 'wu', axes_3d, time, 'omega u', 'Pa/s m/s')
    id_wv = register_diag_field(mod_name, 'wv', axes_3d, time, 'omega v', 'Pa/s m/s')
    id_wq = register_diag_field(mod_name, 'wq', axes_3d, time, 'omega q', 'Pa/s kg/kg')
    id_wdse = register_diag_field(mod_name, 'wdse', axes_3d, time, 'omega DSE', 'Pa/s J/kg')
    id_wup = register_diag_field(mod_name, 'wup', axes_3d, time, 'omega u, weighted by ps/p00', 'Pa/s m/s')
    id_wvp = register_diag_field(mod_name, 'wvp', axes_3d, time, 'omega v, weighted by ps/p00', 'Pa/s m/s')
    id_wqp = register_diag_field(mod_name, 'wqp', axes_3d, time, 'omega q, weighted by ps/p00', 'Pa/s kg/kg')
    id_wdsep = register_diag_field(mod_name, 'wdsep', axes_3d, time, 'omega DSE, weighted by ps/p00', 'Pa/s J/kg')
    id_wp = register_diag_field(mod_name, 'omegap', axes_3d, time, 'dp/dt vertical velocity, weighted by ps/p00', 'Pa/s')
    id_wsubt = register_diag_field(mod_name, 'wsubt', axes_3d, time, 'kappa omega T/p, weighted by ps/p00', 'K/s')
    id_wsubtv = register_diag_field(mod_name, 'wsubtv', axes_3d, time, 'kappa omega Tv/p, weighted by ps/p00', 'K/s')
    id_kegen = register_diag_field(mod_name, 'kegen', axes_3d, time, &
                                   'kappa omega T/p from the dynamical core (Tv if use_virtual_temperature), '// &
                                   'weighted by ps/p00', 'K/s')
    id_kegenq = register_diag_field(mod_name, 'kegenq', axes_3d, time, &
                                    'kappa omega T q/p from the dynamical core (Tv if use_virtual_temperature), '// &
                                    'weighted by ps/p00', 'K/s kg/kg')
    id_kegenqtinv = register_diag_field(mod_name, 'kegenqtinv', axes_3d, time, &
                                        'kappa omega q/p from the dynamical core, weighted by ps/p00', 'kg/kg/s')
    id_drystaten = register_diag_field(mod_name, 'dry_stat_en', axes_3d, time, 'dry static energy', 'J/kg')
    id_moiststaten = register_diag_field(mod_name, 'moist_stat_en', axes_3d, time, 'moist static energy', 'J/kg')
    id_udt_damp = register_diag_field(mod_name, 'udt_damp', axes_3d, time, 'Zonal wind tend from horiz diff', 'm/s2')
    id_vdt_damp = register_diag_field(mod_name, 'vdt_damp', axes_3d, time, 'Merid wind tend from horiz diff', 'm/s2')
    id_tdt_damp = register_diag_field(mod_name, 'tdt_damp', axes_3d, time, 'Temperature tend from horiz diff', 'K/s')
    id_entrop_dampuv = register_diag_field(mod_name, 'entrop_dampuv', &
                                           axes_3d, time, 'Entropy prod from horiz diff of vel', '1/s')
    id_entrop_dampt = register_diag_field(mod_name, 'entrop_dampt', &
                                          axes_3d, time, 'Entropy prod from horiz diff of temp', '1/s')
    id_tdt_dampuv = register_diag_field(mod_name, 'tdt_dampuv', &
                                        axes_3d, time, 'Temp tend (not applied) due to diff of winds', 'K/s')
    id_tdt_tempcor = register_diag_field(mod_name, 'tdt_tempcor', &
                                         axes_3d, time, 'Temp tend due to dynamics energy correction', 'K/s')
    id_qdt_watercor = register_diag_field(mod_name, 'qdt_watercor', &
                                          axes_3d, time, 'Humidity tend due to dynamics water correction', 'kg/kg/s')
    id_entrop_tempcor = register_diag_field(mod_name, 'entrop_tempcor', axes_3d, time, &
                                            'Entropy tend due to dynamics energy correction, weighted by ps/p00', '1/s')

    id_upvp = register_diag_field(mod_name, 'upvp', axes_3d, time, 'meridional eddy momentum flux', 'm2/s2')
    id_upwp = register_diag_field(mod_name, 'upwp', axes_3d, time, 'vertical eddy momentum flux', 'm/s Pa/s')
    id_vptp = register_diag_field(mod_name, 'vpTp', axes_3d, time, 'meridional eddy heat flux', 'm/s K')
    id_psi_dwc = register_diag_field(mod_name, 'psi_dwc', axes_zm, time, &
                                     'residual mean streamfunction from downward control (EP flux divergence)', 'm/s Pa')
    id_psi_star = register_diag_field(mod_name, 'psi_star', axes_zm, time, &
                                      'residual mean streamfunction (Eulerian mean minus eddy heat flux term)', 'm/s Pa')

    call get_number_tracers(MODEL_ATMOS, num_prog=num_tracers)
    allocate (id_tr(num_tracers), id_dt_hadv(num_tracers), id_dt_vadv(num_tracers))
    do ntr = 1, num_tracers
      call get_tracer_names(MODEL_ATMOS, ntr, tname, longname, units)
      id_tr(ntr) = register_diag_field(mod_name, trim(tname)//'_every', axes_3d, Time, longname, units)
      id_dt_hadv(ntr) = register_diag_field(mod_name, trim(tname)//'_hadv', axes_3d, Time, &
                                            'Global mean column integrated '//trim(longname)// &
                                            ' tendency due to horizontal advection', &
                                            trim(units)//' kg/m2/s')
      id_dt_vadv(ntr) = register_diag_field(mod_name, trim(tname)//'_vadv', axes_3d, Time, &
                                            'Global mean column integrated '//trim(longname)// &
                                            ' tendency due to vertical advection', &
                                            trim(units)//' kg/m2/s')
    end do
    nsphum = get_tracer_index(MODEL_ATMOS, 'sphum')

    need_damp_uv = id_udt_damp > 0 .or. id_vdt_damp > 0 .or. id_entrop_dampuv > 0 .or. id_tdt_dampuv > 0
    need_damp_t = id_tdt_damp > 0 .or. id_entrop_dampt > 0
    need_kegen = id_kegen > 0 .or. id_kegenq > 0 .or. id_kegenqtinv > 0
    allocate (need_hadv(num_tracers), need_vadv(num_tracers))
    need_hadv = id_dt_hadv > 0
    need_vadv = id_dt_vadv > 0

    call get_grid_domain(is, ie, js, je)
    num_levels = num_levels_in
    allocate (two_dt_ps(is:ie, js:je))
    allocate (two_dt_u(is:ie, js:je, num_levels))
    allocate (two_dt_v(is:ie, js:je, num_levels))
    allocate (two_dt_t(is:ie, js:je, num_levels))
    allocate (two_dt_tr(is:ie, js:je, num_levels, num_tracers))
    allocate (two_dt_id_tr(num_tracers))

    two_dt_ps = 0.
    two_dt_u = 0.
    two_dt_v = 0.
    two_dt_t = 0.
    two_dt_tr = 0.

    two_dt_id_ps = register_static_field(mod_name, '2dt_ps', axes_2d, 'Amplitude of 2*dt wave in surface pressure', 'Pa')
    two_dt_id_u = register_static_field(mod_name, '2dt_u', axes_3d, 'Amplitude of 2*dt wave in zonal wind', 'm/s')
    two_dt_id_v = register_static_field(mod_name, '2dt_v', axes_3d, 'Amplitude of 2*dt wave in meridional wind', 'm/s')
    two_dt_id_t = register_static_field(mod_name, '2dt_t', axes_3d, 'Amplitude of 2*dt wave in temperature', 'K')
    do ntr = 1, num_tracers
      call get_tracer_names(MODEL_ATMOS, ntr, tname, longname, units)
      two_dt_id_tr(ntr) = &
        register_static_field(mod_name, '2dt_'//trim(tname), axes_3d, 'Amplitude of 2*dt wave in '//trim(longname), units)
    end do
    iwt = 1
    num_time_steps = 0

    module_is_initialized = .true.

    return
  end subroutine every_step_diagnostics_init
!===================================================================================
  !> Sends the `dynamics_every` diagnostics and accumulates the 2*dt wave amplitudes.
  !>
  !> Called at the end of each dynamics (sub)step. `Time` is the end of the (sub)step, where
  !> `p_surf`, `u_grid`, `v_grid`, `t_grid` and `tr_grid` (the new time level) are valid.
  !> `wg_full`, `p_full`, `p_half`, `z_full`, `kegen*` and the damping tendencies were
  !> evaluated at the start of the (sub)step; like the physics tendencies, which are sent at
  !> `Time_next`, they are attributed to the end of the step over which they act.
  subroutine every_step_diagnostics(Time, p_surf, u_grid, v_grid, t_grid, tr_grid, &
                                    wg_full, p_full, p_half, z_full, dt_ug_damp, dt_vg_damp, dt_tg_damp, &
                                    tdt_tempcor, water_correction, dt_hadv, dt_vadv, kegen, kegenq, kegenqtinv)

    type(time_type), intent(in) :: Time  !! end of the (sub)step
    real, intent(in), dimension(is:ie, js:je)                         :: p_surf  !! surface pressure [Pa]
    real, intent(in), dimension(is:ie, js:je, num_levels)             :: u_grid, v_grid, t_grid
    !! `u_grid`, `v_grid`: zonal and meridional wind [m/s]; `t_grid`: temperature [K]
    real, intent(in), dimension(is:ie, js:je, num_levels, num_tracers) :: tr_grid  !! grid-point tracers
    real, intent(in), dimension(is:ie, js:je, num_levels)             :: wg_full, p_full, z_full
    !! `wg_full`: vertical pressure velocity [Pa/s]; `p_full`: pressure [Pa]; `z_full`: height [m];
    !! all at full levels
    real, intent(in), dimension(is:ie, js:je, num_levels)             :: dt_ug_damp, dt_vg_damp, dt_tg_damp
    !! `dt_ug_damp`, `dt_vg_damp`, `dt_tg_damp`: tendencies of zonal and meridional wind [m/s2] and
    !! temperature [K/s] from the hyperdiffusion
    real, intent(in), dimension(is:ie, js:je, num_levels + 1)           :: p_half  !! pressure at half levels [Pa]
    real, intent(in)                                                  :: tdt_tempcor
    !! temperature tendency of the energy correction, the same at all points [K/s]
    real, intent(in), dimension(num_tracers)                          :: dt_hadv, dt_vadv
    !! `dt_hadv`, `dt_vadv`: global mean column-integrated tracer tendencies due to horizontal and
    !! vertical advection [(tracer units) kg/m2/s]
    real, intent(in), dimension(is:ie, js:je, num_levels)             :: water_correction, kegen, kegenq, kegenqtinv
    !! `water_correction`: humidity tendency of the water correction [kg/kg/s]; `kegen`:
    !! kappa omega T/p [K/s]; `kegenq`: `kegen` times specific humidity [K/s kg/kg]; `kegenqtinv`:
    !! kappa omega q/p [kg/kg/s] (T is the virtual temperature with `use_virtual_temperature`)

    real, dimension(is:ie, js:je, num_levels) :: tempdiag_3d, dse, q_grid
    real, dimension(is:ie, js:je)             :: tempdiag_2d

    logical :: used
    integer :: i, j, k, ntr
    character(len=8) :: err_msg_1, err_msg_2

    ! zonal means [.] (latitude, level) for the eddy fluxes and the streamfunctions
    real, dimension(is:ie, js:je, num_levels) :: theta
    real, dimension(js:je, num_levels) :: tempdiag_zm, u_zm, v_zm, w_zm, th_zm, p_zm, vpthp, dthdp, &
                                          upvp, upwp, dudp, vth_strat, fhat, f_p, df_p, div_f, vstar
    real, dimension(size(lat), num_levels) :: f_global, df_global
    real, dimension(js:je) :: coslat

    if (.not. module_is_initialized) then
      call error_mesg('every_step_diagnostics', 'module is not initialized', FATAL)
    end if

    Time_save = Time

    if (id_ps > 0) used = send_data(id_ps, p_surf, Time)
    if (id_u > 0) used = send_data(id_u, u_grid, Time)
    if (id_v > 0) used = send_data(id_v, v_grid, Time)
    if (id_t > 0) used = send_data(id_t, t_grid, Time)

    if (id_drystaten > 0 .or. id_vdse > 0 .or. id_wdse > 0 .or. id_wdsep > 0) then
      dse = cp_air*t_grid + grav*z_full ! Dry static energy
    end if

    if (id_moiststaten > 0 .or. id_vq > 0 .or. id_vqint > 0 .or. id_wq > 0 .or. id_wqp > 0 .or. id_wsubtv > 0) then
      q_grid = tr_grid(:, :, :, nsphum)
    end if

    if (id_drystaten > 0) then
      used = send_data(id_drystaten, dse, Time)
    end if
    if (id_moiststaten > 0) then
      tempdiag_3d = dse + hlv*q_grid ! Moist static energy
      used = send_data(id_moiststaten, tempdiag_3d, Time)
    end if
    if (id_uv > 0) then
      do k = 1, num_levels
        tempdiag_3d(:, :, k) = u_grid(:, :, k)*v_grid(:, :, k)*p_surf/p00
      end do
      used = send_data(id_uv, tempdiag_3d, Time)
    end if
    if (id_vq > 0) then
      do k = 1, num_levels
        tempdiag_3d(:, :, k) = v_grid(:, :, k)*q_grid(:, :, k)*p_surf/p00
      end do
      used = send_data(id_vq, tempdiag_3d, Time)
    end if
    if (id_vqint > 0) then
      tempdiag_2d(:, :) = 0.
      do k = 1, num_levels
        tempdiag_2d(:, :) = tempdiag_2d(:, :) + &
                            v_grid(:, :, k)*q_grid(:, :, k)*(p_half(:, :, k + 1) - p_half(:, :, k))/grav
      end do
      used = send_data(id_vqint, tempdiag_2d, Time)
    end if
    if (id_vdse > 0) then
      do k = 1, num_levels
        tempdiag_3d(:, :, k) = v_grid(:, :, k)*dse(:, :, k)*p_surf/p00
      end do
      used = send_data(id_vdse, tempdiag_3d, Time)
    end if
    if (id_vp > 0) then
      do k = 1, num_levels
        tempdiag_3d(:, :, k) = v_grid(:, :, k)*p_surf/p00
      end do
      used = send_data(id_vp, tempdiag_3d, Time)
    end if
    if (id_wu > 0) then
      tempdiag_3d = wg_full*u_grid
      used = send_data(id_wu, tempdiag_3d, Time)
    end if
    if (id_wv > 0) then
      tempdiag_3d = wg_full*v_grid
      used = send_data(id_wv, tempdiag_3d, Time)
    end if
    if (id_wq > 0) then
      tempdiag_3d = wg_full*q_grid
      used = send_data(id_wq, tempdiag_3d, Time)
    end if
    if (id_wdse > 0) then
      tempdiag_3d = wg_full*dse
      used = send_data(id_wdse, tempdiag_3d, Time)
    end if
    if (id_wup > 0) then
      do k = 1, num_levels
        tempdiag_3d(:, :, k) = wg_full(:, :, k)*u_grid(:, :, k)*p_surf/p00
      end do
      used = send_data(id_wup, tempdiag_3d, Time)
    end if
    if (id_wvp > 0) then
      do k = 1, num_levels
        tempdiag_3d(:, :, k) = wg_full(:, :, k)*v_grid(:, :, k)*p_surf/p00
      end do
      used = send_data(id_wvp, tempdiag_3d, Time)
    end if
    if (id_wqp > 0) then
      do k = 1, num_levels
        tempdiag_3d(:, :, k) = wg_full(:, :, k)*q_grid(:, :, k)*p_surf/p00
      end do
      used = send_data(id_wqp, tempdiag_3d, Time)
    end if
    if (id_wdsep > 0) then
      do k = 1, num_levels
        tempdiag_3d(:, :, k) = wg_full(:, :, k)*dse(:, :, k)*p_surf/p00
      end do
      used = send_data(id_wdsep, tempdiag_3d, Time)
    end if
    if (id_wp > 0) then
      do k = 1, num_levels
        tempdiag_3d(:, :, k) = wg_full(:, :, k)*p_surf/p00
      end do
      used = send_data(id_wp, tempdiag_3d, Time)
    end if
    if (id_wsubt > 0) then
      do k = 1, num_levels
        tempdiag_3d(:, :, k) = kappa*wg_full(:, :, k)*t_grid(:, :, k)/p_full(:, :, k)*p_surf/p00
      end do
      used = send_data(id_wsubt, tempdiag_3d, Time)
    end if
    if (id_wsubtv > 0) then
      do k = 1, num_levels
        tempdiag_3d(:, :, k) = kappa*wg_full(:, :, k)*t_grid(:, :, k)*(1. + virtual_factor*q_grid(:, :, k))/ &
                               p_full(:, :, k)*p_surf/p00
      end do
      used = send_data(id_wsubtv, tempdiag_3d, Time)
    end if
    if (id_kegen > 0) then
      do k = 1, num_levels
        tempdiag_3d(:, :, k) = kegen(:, :, k)*p_surf/p00
      end do
      used = send_data(id_kegen, tempdiag_3d, Time)
    end if
    if (id_kegenq > 0) then
      do k = 1, num_levels
        tempdiag_3d(:, :, k) = kegenq(:, :, k)*p_surf/p00
      end do
      used = send_data(id_kegenq, tempdiag_3d, Time)
    end if
    if (id_kegenqtinv > 0) then
      do k = 1, num_levels
        tempdiag_3d(:, :, k) = kegenqtinv(:, :, k)*p_surf/p00
      end do
      used = send_data(id_kegenqtinv, tempdiag_3d, Time)
    end if
    if (id_udt_damp > 0) then
      used = send_data(id_udt_damp, dt_ug_damp, Time)
    end if
    if (id_vdt_damp > 0) then
      used = send_data(id_vdt_damp, dt_vg_damp, Time)
    end if
    if (id_tdt_damp > 0) then
      used = send_data(id_tdt_damp, dt_tg_damp, Time)
    end if
    if (id_entrop_dampuv > 0) then
      do k = 1, num_levels
        tempdiag_3d(:, :, k) = -1./cp_air*(u_grid(:, :, k)*dt_ug_damp(:, :, k) + &
                                           v_grid(:, :, k)*dt_vg_damp(:, :, k))*p_surf/p00/t_grid(:, :, k)
      end do
      used = send_data(id_entrop_dampuv, tempdiag_3d, Time)
    end if
    if (id_entrop_dampt > 0) then
      do k = 1, num_levels
        tempdiag_3d(:, :, k) = dt_tg_damp(:, :, k)/t_grid(:, :, k)*p_surf/p00
      end do
      used = send_data(id_entrop_dampt, tempdiag_3d, Time)
    end if
    if (id_tdt_dampuv > 0) then
      tempdiag_3d = -1./cp_air*(u_grid*dt_ug_damp + v_grid*dt_vg_damp)
      used = send_data(id_tdt_dampuv, tempdiag_3d, Time)
    end if
    if (id_tdt_tempcor > 0) then
      tempdiag_3d = tdt_tempcor
      used = send_data(id_tdt_tempcor, tempdiag_3d, Time)
    end if
    if (id_qdt_watercor > 0) then
      used = send_data(id_qdt_watercor, water_correction, Time)
    end if
    if (id_entrop_tempcor > 0) then
      do k = 1, num_levels
        tempdiag_3d(:, :, k) = tdt_tempcor/t_grid(:, :, k)*p_surf/p00
      end do
      used = send_data(id_entrop_tempcor, tempdiag_3d, Time)
    end if
    ! zonal means: each is summed over the longitudes once (in the same order as a sum per point)
    if (id_upvp > 0 .or. id_upwp > 0 .or. id_psi_dwc > 0 .or. id_psi_star > 0) then
      u_zm = sum(u_grid, 1)/max(1, size(u_grid, 1))
    end if
    if (id_upvp > 0 .or. id_vptp > 0 .or. id_psi_dwc > 0 .or. id_psi_star > 0) then
      v_zm = sum(v_grid, 1)/max(1, size(v_grid, 1))
    end if
    if (id_upwp > 0 .or. id_psi_dwc > 0) then
      w_zm = sum(wg_full, 1)/max(1, size(wg_full, 1))
    end if

    if (id_upvp > 0) then
      do k = 1, num_levels
        do j = js, je
          tempdiag_3d(:, j, k) = (u_grid(:, j, k) - u_zm(j, k))*(v_grid(:, j, k) - v_zm(j, k))
        end do
      end do
      used = send_data(id_upvp, tempdiag_3d, Time)
    end if
    if (id_upwp > 0) then
      do k = 1, num_levels
        do j = js, je
          tempdiag_3d(:, j, k) = (u_grid(:, j, k) - u_zm(j, k))*(wg_full(:, j, k) - w_zm(j, k))
        end do
      end do
      used = send_data(id_upwp, tempdiag_3d, Time)
    end if
    if (id_vptp > 0) then
      tempdiag_zm = sum(t_grid, 1)/max(1, size(t_grid, 1))
      do k = 1, num_levels
        do j = js, je
          tempdiag_3d(:, j, k) = (v_grid(:, j, k) - v_zm(j, k))*(t_grid(:, j, k) - tempdiag_zm(j, k))
        end do
      end do
      used = send_data(id_vptp, tempdiag_3d, Time)
    end if

    if (id_psi_dwc > 0 .or. id_psi_star > 0) then
      ! potential temperature, [v'theta'] and the stratification d[theta]/dp
      theta = t_grid*(p00/p_full)**kappa
      th_zm = sum(theta, 1)/max(1, size(theta, 1))
      p_zm = sum(p_full, 1)/max(1, size(p_full, 1))
      do k = 1, num_levels
        do j = js, je
          vpthp(j, k) = sum((v_grid(:, j, k) - v_zm(j, k))*(theta(:, j, k) - th_zm(j, k)))/max(1, size(v_grid, 1))
        end do
      end do
      call p_derivative(th_zm, p_zm, dthdp)
      vth_strat = vpthp/dthdp   ! [v'theta']/[theta]_p [m/s Pa]
    end if

    if (id_psi_dwc > 0) then
      ! Residual mean meridional velocity from the steady TEM zonal momentum balance (downward
      ! control; Andrews et al. 1987, eqs. 3.5.2-3.5.3, in pressure coordinates):
      !   -fhat v* = d(cos**2 F_phi)/dphi/(a cos**2) + d(F_p)/dp,
      !   F_phi = [u]_p [v'theta']/[theta]_p - [u'v'],   F_p = fhat [v'theta']/[theta]_p - [u'omega'],
      !   fhat = f - d([u] cos)/dphi/(a cos);
      ! psi_dwc is v* integrated in pressure from the top level. The latitude derivatives use the
      ! zonal means of all latitudes, gathered from all PEs.
      do k = 1, num_levels
        do j = js, je
          upvp(j, k) = sum((u_grid(:, j, k) - u_zm(j, k))*(v_grid(:, j, k) - v_zm(j, k)))/max(1, size(u_grid, 1))
          upwp(j, k) = sum((u_grid(:, j, k) - u_zm(j, k))*(wg_full(:, j, k) - w_zm(j, k)))/max(1, size(u_grid, 1))
        end do
      end do
      call p_derivative(u_zm, p_zm, dudp)
      do j = js, je
        coslat(j) = cos(lat(j)/radian)
      end do

      ! fhat
      do k = 1, num_levels
        tempdiag_zm(:, k) = u_zm(:, k)*coslat
      end do
      call gather_latitudes(tempdiag_zm, f_global)
      call lat_derivative(f_global, df_global)
      do j = js, je
        fhat(j, :) = 2*omega*sin(lat(j)/radian) - df_global(j, :)/(radius*coslat(j))
      end do

      ! divergence of the EP flux
      do k = 1, num_levels
        tempdiag_zm(:, k) = (dudp(:, k)*vth_strat(:, k) - upvp(:, k))*coslat**2
      end do
      call gather_latitudes(tempdiag_zm, f_global)
      call lat_derivative(f_global, df_global)
      do j = js, je
        div_f(j, :) = df_global(j, :)/(radius*coslat(j)**2)
      end do
      f_p = fhat*vth_strat - upwp
      call p_derivative(f_p, p_zm, df_p)
      div_f = div_f + df_p

      vstar = -div_f/fhat
      tempdiag_zm(:, 1) = 0.
      do k = 2, num_levels
        tempdiag_zm(:, k) = tempdiag_zm(:, k - 1) + 0.5*(p_zm(:, k) - p_zm(:, k - 1))*(vstar(:, k) + vstar(:, k - 1))
      end do
      used = send_data(id_psi_dwc, tempdiag_zm, Time)
    end if

    if (id_psi_star > 0) then
      ! psi* = [int v dp] - [v'theta']/[theta]_p, the column integral from the top level
      tempdiag_3d(:, :, 1) = 0.
      do k = 2, num_levels
        tempdiag_3d(:, :, k) = tempdiag_3d(:, :, k - 1) + 0.5*(p_full(:, :, k) - p_full(:, :, k - 1))* &
                               (v_grid(:, :, k) + v_grid(:, :, k - 1))
      end do
      tempdiag_zm = sum(tempdiag_3d, 1)/max(1, size(tempdiag_3d, 1)) - vth_strat
      used = send_data(id_psi_star, tempdiag_zm, Time)
    end if

    if (size(tr_grid, 4) /= num_tracers) then
      write (err_msg_1, '(i8)') size(tr_grid, 4)
      write (err_msg_2, '(i8)') num_tracers
      call error_mesg('every_step_diagnostics', 'size(tracers)='//err_msg_1//' Should be='//err_msg_2, FATAL)
    end if
    do ntr = 1, num_tracers
      if (id_tr(ntr) > 0) used = send_data(id_tr(ntr), tr_grid(:, :, :, ntr), Time)
      ! global means, sent as 3-d fields of one value
      if (id_dt_hadv(ntr) > 0) then
        tempdiag_3d = dt_hadv(ntr)
        used = send_data(id_dt_hadv(ntr), tempdiag_3d, Time)
      end if
      if (id_dt_vadv(ntr) > 0) then
        tempdiag_3d = dt_vadv(ntr)
        used = send_data(id_dt_vadv(ntr), tempdiag_3d, Time)
      end if
    end do

    if (two_dt_id_ps > 0) then
      do j = js, je
        do i = is, ie
          two_dt_ps(i, j) = two_dt_ps(i, j) + iwt*p_surf(i, j)
        end do
      end do
    end if
    if (two_dt_id_u > 0) then
      do k = 1, num_levels
        do j = js, je
          do i = is, ie
            two_dt_u(i, j, k) = two_dt_u(i, j, k) + iwt*u_grid(i, j, k)
          end do
        end do
      end do
    end if
    if (two_dt_id_v > 0) then
      do k = 1, num_levels
        do j = js, je
          do i = is, ie
            two_dt_v(i, j, k) = two_dt_v(i, j, k) + iwt*v_grid(i, j, k)
          end do
        end do
      end do
    end if
    if (two_dt_id_t > 0) then
      do k = 1, num_levels
        do j = js, je
          do i = is, ie
            two_dt_t(i, j, k) = two_dt_t(i, j, k) + iwt*t_grid(i, j, k)
          end do
        end do
      end do
    end if
    do ntr = 1, num_tracers
      if (two_dt_id_tr(ntr) > 0) then
        do k = 1, num_levels
          do j = js, je
            do i = is, ie
              two_dt_tr(i, j, k, ntr) = two_dt_tr(i, j, k, ntr) + iwt*tr_grid(i, j, k, ntr)
            end do
          end do
        end do
      end if
    end do

    iwt = -1*iwt
    num_time_steps = num_time_steps + 1

    return
  end subroutine every_step_diagnostics
!===================================================================================
  !> Sends the 2*dt wave amplitudes and deallocates the module arrays.
  subroutine every_step_diagnostics_end(Time_in)
    type(time_type), intent(in), optional :: Time_in
    !! time for `send_data` (default: the time of the last call of `every_step_diagnostics`)
    logical :: used
    integer :: ntr
    type(time_type) :: Time

    if (.not. module_is_initialized) return

    if (present(Time_in)) then
      Time = Time_in
    else
      Time = Time_save
    end if

    if (num_time_steps > 0) then
      two_dt_ps = two_dt_ps/num_time_steps
      two_dt_u = two_dt_u/num_time_steps
      two_dt_v = two_dt_v/num_time_steps
      two_dt_t = two_dt_t/num_time_steps
      two_dt_tr = two_dt_tr/num_time_steps
    end if

    if (two_dt_id_ps > 0) used = send_data(two_dt_id_ps, two_dt_ps)
    if (two_dt_id_u > 0) used = send_data(two_dt_id_u, two_dt_u)
    if (two_dt_id_v > 0) used = send_data(two_dt_id_v, two_dt_v)
    if (two_dt_id_t > 0) used = send_data(two_dt_id_t, two_dt_t)
    do ntr = 1, num_tracers
      if (two_dt_id_tr(ntr) > 0) used = send_data(two_dt_id_tr(ntr), two_dt_tr(:, :, :, ntr))
    end do

    deallocate (id_tr, id_dt_hadv, id_dt_vadv, two_dt_id_tr, need_hadv, need_vadv)
    deallocate (two_dt_ps, two_dt_u, two_dt_v, two_dt_t, two_dt_tr)
    deallocate (lon, lat)
    module_is_initialized = .false.

    return
  end subroutine every_step_diagnostics_end
!===================================================================================

  ! Vertical derivative d(f)/dp of zonal-mean fields (latitude, level): centred inside,
  ! one-sided at the top and bottom levels.
  subroutine p_derivative(f, p, dfdp)
    real, intent(in), dimension(:, :) :: f, p
    real, intent(out), dimension(:, :) :: dfdp
    integer :: k, n

    n = size(f, 2)
    if (n < 2) then
      dfdp = 0.
      return
    end if
    dfdp(:, 1) = (f(:, 2) - f(:, 1))/(p(:, 2) - p(:, 1))
    do k = 2, n - 1
      dfdp(:, k) = (f(:, k + 1) - f(:, k - 1))/(p(:, k + 1) - p(:, k - 1))
    end do
    dfdp(:, n) = (f(:, n) - f(:, n - 1))/(p(:, n) - p(:, n - 1))

  end subroutine p_derivative
!===================================================================================

  ! Meridional derivative d(f)/d(latitude) [per radian] of zonal-mean fields on all
  ! latitudes (lat_max, level): centred inside, one-sided at the first and last latitude.
  subroutine lat_derivative(f, dfdphi)
    real, intent(in), dimension(:, :) :: f
    real, intent(out), dimension(:, :) :: dfdphi
    integer :: j, n

    n = size(f, 1)
    dfdphi(1, :) = (f(2, :) - f(1, :))/((lat(2) - lat(1))/radian)
    do j = 2, n - 1
      dfdphi(j, :) = (f(j + 1, :) - f(j - 1, :))/((lat(j + 1) - lat(j - 1))/radian)
    end do
    dfdphi(n, :) = (f(n, :) - f(n - 1, :))/((lat(n) - lat(n - 1))/radian)

  end subroutine lat_derivative
!===================================================================================

  ! Returns the zonal-mean field f of this PE's latitudes js:je on all latitudes. Every PE
  ! must call it. Each latitude is on one PE only (the grid is decomposed in latitude), so
  ! the sum over the PEs adds exact zeros and does not depend on the decomposition.
  subroutine gather_latitudes(f, f_global)
    real, intent(in), dimension(js:, :) :: f
    real, intent(out), dimension(:, :) :: f_global
    real, dimension(size(f_global, 1)*size(f_global, 2)) :: buffer

    f_global = 0.
    f_global(js:je, :) = f(js:je, :)
    buffer = reshape(f_global, (/size(buffer)/))
    call mpp_sum(buffer, size(buffer))
    f_global = reshape(buffer, shape(f_global))

  end subroutine gather_latitudes
!===================================================================================

end module every_step_diagnostics_mod
