
!> High-level interface to the atmospheric physics.
!>
!> Calls the physics modules and returns the tendencies and the boundary fluxes that drive
!> the atmosphere and force the surface model. It is designed around the implicit
!> vertical diffusion scheme and advances the model one time step in two passes, which
!> correspond to the down and up sweeps of the tridiagonal solver:
!>
!> * `physics_driver_down`: radiation, the Held-Suarez forcing (`do_held_suarez`), local
!>   heating (`do_local_heating`), damping and gravity-wave drag (`do_damping`),
!>   boundary-layer turbulence, tracers, and the downward pass of the vertical diffusion
!>   (`do_boundary_layer`).
!> * `physics_driver_up`: the upward pass of the vertical diffusion, then convection and
!>   large-scale condensation (`do_moist_physics`).
!>
!> The diffusion coefficients, optionally smoothed in time (`diffusion_smooth`), are kept
!> in the restart file `RESTART/physics_driver.res.nc`.
!>
!> Namelist: `physics_driver_nml`
!> ([namelist reference](https://eddy-stanford.github.io/MiMA/Parameters/#physics_driver_nml)).
!>
!> Original authors: Fei Liu.
module physics_driver_mod
!   shared modules:

  use time_manager_mod, only: time_type, get_time, operator(-), &
                              time_manager_init
  use field_manager_mod, only: field_manager_init, MODEL_ATMOS
  use tracer_manager_mod, only: tracer_manager_init, &
                                get_number_tracers

  use atmos_tracer_driver_mod, only: atmos_tracer_driver_init, &
                                     atmos_tracer_driver, &
                                     atmos_tracer_driver_end
  use fms_mod, only: mpp_clock_id, mpp_clock_begin, &
                     mpp_clock_end, CLOCK_MODULE_DRIVER, &
                     MPP_CLOCK_SYNC, fms_init, &
                     input_nml_file, stdlog, &
                     write_version_number, &
                     error_mesg, FATAL, &
                     WARNING, NOTE, check_nml_error, &
                     mpp_pe, mpp_root_pe, &
                     mpp_error, mpp_chksum
  use mpp_domains_mod, only: domain2d
  use restart_file_mod, only: restart_file_type, open_restart_read, &
                              open_restart_write, close_restart, &
                              read_restart_field, write_restart_field
  use sat_vapor_pres_mod, only: sat_vapor_pres_init, lookup_es
  use constants_mod, only: ES0, HLV, RVGAS, TFREEZE, RADIUS, OMEGA, GRAV, &
                           RDGAS, KAPPA, HLF, STEFAN

!    component modules:

  use moist_processes_mod, only: moist_processes, &
                                 moist_processes_init, &
                                 moist_processes_end

  use vert_turb_driver_mod, only: vert_turb_driver, &
                                  vert_turb_driver_init, &
                                  vert_turb_driver_end

  use vert_diff_driver_mod, only: vert_diff_driver_down, &
                                  vert_diff_driver_up, &
                                  vert_diff_driver_init, &
                                  vert_diff_driver_end, &
                                  surf_diff_type

  use damping_driver_mod, only: damping_driver, &
                                damping_driver_init, &
                                damping_driver_end

  use radiation_mod, only: radiation_init, radiation_down, radiation_end

  use held_suarez_mod, only: held_suarez_init, held_suarez_forcing, held_suarez_end

  use local_heating_mod, only: local_heating_init, local_heating

!-----------------------------------------------------------------

  implicit none
  private

!---------------------------------------------------------------------
!----------- version number for this module -------------------

  character(len=128) :: version = '$Id: physics_driver.f90,v 12.0.6.2 2005/05/16 13:56:54 pjp Exp $'
  character(len=128) :: tagname = '$Name:  $'

!---------------------------------------------------------------------
!-------  interfaces --------

  public physics_driver_init, physics_driver_down, &
    physics_driver_up, physics_driver_end, &
    do_local_heating, surface_is_coupled

  private &
    !  called from physics_driver_init:
    read_restart_nc, check_constants, check_sat_vapor_pres, &
    !  called from physics_driver_down:
    check_args, &
    !  called from check_args:
    check_dim

  interface check_dim
    module procedure check_dim_2d, check_dim_3d, check_dim_4d
  end interface

!---------------------------------------------------------------------
!------- namelist ------

  real    :: tau_diff = 3600.    !! [s] time scale for smoothing the diffusion coefficients in time

  logical :: do_damping = .true.  !! Rayleigh sponge and gravity-wave drag (`damping_driver_nml`)

  logical :: do_local_heating = .false.  !! add prescribed local heating (`local_heating_nml`)

  logical :: do_held_suarez = .false.    !! add the Held-Suarez (1994) forcing (`held_suarez_nml`)
  logical :: do_boundary_layer = .true.
  !! boundary-layer turbulence, vertical diffusion and coupling to the surface fluxes. With
  !! `.false.` the surface state is not updated.
  logical :: do_moist_physics = .true.   !! convection and large-scale condensation (`moist_processes_nml`)

  real    :: diff_min = 1.e-3    !! [m2/s] diffusion coefficients below this are set to zero
  logical :: diffusion_smooth = .true.  !! smooth the diffusion coefficients in time
  namelist /physics_driver_nml/ tau_diff, &
    diff_min, diffusion_smooth, &
    do_damping, do_local_heating, &
    do_held_suarez, do_boundary_layer, do_moist_physics

!---------------------------------------------------------------------
!------- public data ------

  public surf_diff_type   ! defined in  vert_diff_driver_mod, republished
  ! here

!---------------------------------------------------------------------
!------- private data ------

!--------------------------------------------------------------------
! list of restart versions readable by this module:
!
! version 1: initial implementation 1/2003, contains diffusion coef-
!            ficient contribution from cu_mo_trans_mod. This variable
!            is generated in physics_driver_up (moist_processes) and
!            used on the next step in vert_diff_down, necessitating
!            its storage.
!
! version 2: adds pbltop as generated in vert_turb_driver_mod. This
!            variable is then used on the next timestep by topo_drag
!            (called from damping_driver_mod), necessitating its
!            storage.
!
! version 3: adds the diffusion coefficients which are passed to
!            vert_diff_driver.  These diffusion are saved should
!            smoothing of vertical diffusion coefficients be turned
!            on.
!
! version 4: adds a logical variable, convect, which indicates whether
!            or not the grid column is convecting. This diagnostic is
!            needed by the entrain_module in vert_turb_driver.
!
! version 5: adds radturbten when strat_cloud_mod is active, adds
!            lw_tendency when edt_mod or entrain_mod is active.
!
!---------------------------------------------------------------------
  integer, dimension(5) :: restart_versions = (/1, 2, 3, 4, 5/)

!--------------------------------------------------------------------
!    the following allocatable arrays are either used to hold physics
!    data between timesteps when required, or hold physics data between
!    physics_down and physics_up.
!
!    diff_t         vertical diffusion coefficient for temperature
!                   which optionally may be time smoothed, meaning
!                   values must be saved between steps
!    diff_m         vertical diffusion coefficient for momentum
!                   which optionally may be time smoothed, meaning
!                   values must be saved between steps
!    lw_tendency    longwave heating rate, generated in radiation and
!                   needed in vert_turb_driver when either edt_mod
!                   or entrain_mod is active. must be saved because
!                   radiation is not calculated on each step.
!    pbltop         top of boundary layer obtained from vert_turb_driver
!                   and then used on the next timestep in topo_drag_mod
!                   called from damping_driver_down
!    convect        flag indicating whether convection is occurring in
!                   a grid column. generated in physics_driver_up and
!                   then used in vert_turb_driver called from
!                   physics_driver_down on the next step.
!----------------------------------------------------------------------
  real, dimension(:, :, :), allocatable :: diff_t, diff_m

  type(domain2d) :: domain  ! grid domain, for the restart file

!---------------------------------------------------------------------
!    internal timing clock variables:
!---------------------------------------------------------------------
  integer :: radiation_clock, damping_clock, turb_clock, &
             tracer_clock, diff_up_clock, diff_down_clock, &
             moist_processes_clock

!--------------------------------------------------------------------
!    miscellaneous control variables:
!---------------------------------------------------------------------
  logical   :: do_check_args = .true.   ! argument dimensions should
  ! be checked ?
  logical   :: module_is_initialized = .false.
  ! module has been initialized ?
  integer   :: nt                       ! total no. of tracers
  integer   :: ntp                      ! total no. of prognostic tracers
!---------------------------------------------------------------------
!---------------------------------------------------------------------

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################
  !> Initializes the module: reads `physics_driver_nml`, checks the FMS physical constants
  !> and saturation vapour pressure, initializes the physics modules and reads
  !> `INPUT/physics_driver.res.nc` if it exists.
  subroutine physics_driver_init(Time, lonb, latb, domain_in, axes, pref, &
                                 trs, Surf_diff, phalf, mask, kbot, &
                                 diffm, difft)

    type(time_type), intent(in)              :: Time  !! current time
    real, dimension(:), intent(in)              :: lonb, latb
    !! longitudes and latitudes of the grid box edges [rad]
    type(domain2d), intent(in)              :: domain_in  !! domain decomposition of the model grid
    integer, dimension(4), intent(in)              :: axes  !! diagnostic axes (lon, lat, pfull, phalf)
    real, dimension(:, :), intent(in)              :: pref
    !! two reference pressure profiles at nlev+1 levels, with surface pressures
    !! `pref(nlev+1,1)` = 101325 and `pref(nlev+1,2)` = 81060 [Pa]
    real, dimension(:, :, :, :), intent(inout)           :: trs  !! atmospheric tracer fields
    type(surf_diff_type), intent(inout)           :: Surf_diff  !! surface data of the implicit vertical diffusion
    real, dimension(:, :, :), intent(in)              :: phalf  !! pressure at half levels [Pa]
    real, dimension(:, :, :), intent(in), optional  :: mask
    !! mask to remove points below ground (eta vertical coordinate)
    integer, dimension(:, :), intent(in), optional  :: kbot
    !! index of the lowest model level above ground (eta vertical coordinate)
    real, dimension(:, :, :), intent(out), optional  :: diffm, difft
    !! diffusion coefficients for momentum and temperature from the restart file (zero without
    !! one) [m2/s]

!---------------------------------------------------------------------
!  local variables:

    real, dimension(size(lonb(:)) - 1, size(latb(:)) - 1) :: sgsmtn
    integer          ::  id, jd, kd
    integer          ::  ierr, io, unit
    integer          ::  ndum

!---------------------------------------------------------------------
!  local variables:
!
!       sgsmtn        sgs orography obtained from mg_drag_mod;
!                     appears to not be currently used
!       id,jd,kd      model dimensions on the processor
!       ierr          error code
!       io            io status returned from an io call
!       unit          unit number used for an i/ operation

!---------------------------------------------------------------------
!    if routine has already been executed, return.
!---------------------------------------------------------------------
    if (module_is_initialized) return

!---------------------------------------------------------------------
!    verify that the modules used by this module that are not called
!    later in this subroutine have already been initialized.
!---------------------------------------------------------------------
    call fms_init

!--------------------------------------------------------------------
!    read namelist.
!--------------------------------------------------------------------
    read (input_nml_file, nml=physics_driver_nml, iostat=io)
    ierr = check_nml_error(io, 'physics_driver_nml')

    call time_manager_init
    call tracer_manager_init
    call field_manager_init(ndum)
    call check_constants
    call check_sat_vapor_pres

!--------------------------------------------------------------------
!    write version number and namelist to log file.
!--------------------------------------------------------------------
    call write_version_number(version, tagname)
    if (mpp_pe() == mpp_root_pe()) &
      write (stdlog(), nml=physics_driver_nml)

!---------------------------------------------------------------------
!    define the model dimensions on the local processor.
!---------------------------------------------------------------------
    id = size(lonb(:)) - 1
    jd = size(latb(:)) - 1
    kd = size(trs, 3)
    call get_number_tracers(MODEL_ATMOS, num_tracers=nt, &
                            num_prog=ntp)

!-----------------------------------------------------------------------
    call moist_processes_init(id, jd, kd, lonb, latb, pref(:, 1), &
                              axes, Time)

!-----------------------------------------------------------------------
!    initialize damping_driver_mod.
!-----------------------------------------------------------------------
    if (do_damping) &
      call damping_driver_init(lonb, latb, domain_in, pref(:, 1), axes, Time, &
                               sgsmtn)

!-----------------------------------------------------------------------
!    initialize vert_turb_driver_mod.
!-----------------------------------------------------------------------
    call vert_turb_driver_init(axes, Time)

!-----------------------------------------------------------------------
!    initialize vert_diff_driver_mod.
!-----------------------------------------------------------------------
    call vert_diff_driver_init(Surf_diff, id, jd, kd, axes, Time)

    call radiation_init(axes, Time, id, jd, kd, lonb, latb, domain_in)

    if (do_held_suarez) call held_suarez_init(axes, Time)

    if (do_local_heating) call local_heating_init(axes, Time)

!-----------------------------------------------------------------------
!    initialize atmos_tracer_driver_mod.
!-----------------------------------------------------------------------
    call atmos_tracer_driver_init(lonb, latb, trs, axes, time, &
                                  phalf, mask)

!---------------------------------------------------------------------
!    initialize  various clocks used to time the physics components.
!---------------------------------------------------------------------
    radiation_clock = &
      mpp_clock_id('   Physics_down: Radiation', &
                   grain=CLOCK_MODULE_DRIVER, flags=MPP_CLOCK_SYNC)
    damping_clock = &
      mpp_clock_id('   Physics_down: Damping', &
                   grain=CLOCK_MODULE_DRIVER, flags=MPP_CLOCK_SYNC)
    turb_clock = &
      mpp_clock_id('   Physics_down: Vert. Turb.', &
                   grain=CLOCK_MODULE_DRIVER, flags=MPP_CLOCK_SYNC)
    tracer_clock = &
      mpp_clock_id('   Physics_down: Tracer', &
                   grain=CLOCK_MODULE_DRIVER, flags=MPP_CLOCK_SYNC)
    diff_down_clock = &
      mpp_clock_id('   Physics_down: Vert. Diff.', &
                   grain=CLOCK_MODULE_DRIVER, flags=MPP_CLOCK_SYNC)
    diff_up_clock = &
      mpp_clock_id('   Physics_up: Vert. Diff.', &
                   grain=CLOCK_MODULE_DRIVER, flags=MPP_CLOCK_SYNC)
    moist_processes_clock = &
      mpp_clock_id('   Physics_up: Moist Processes', &
                   grain=CLOCK_MODULE_DRIVER, flags=MPP_CLOCK_SYNC)

!---------------------------------------------------------------------
!    allocate space for the module variables.
!---------------------------------------------------------------------
    allocate (diff_t(id, jd, kd))
    allocate (diff_m(id, jd, kd))

!--------------------------------------------------------------------
!    obtain initial values for the module variables from the restart
!    file, or initialize them if there is none.
!--------------------------------------------------------------------
    domain = domain_in
    call read_restart_nc

!---------------------------------------------------------------------
!    if desired, define variables to return diff_m and diff_t.
!---------------------------------------------------------------------
    if (present(difft)) then
      difft = diff_t
    end if
    if (present(diffm)) then
      diffm = diff_m
    end if

!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
    module_is_initialized = .true.

!-----------------------------------------------------------------------

  end subroutine physics_driver_init

!######################################################################
  !> Computes the first-pass physics tendencies (radiation, prescribed forcings, damping
  !> and turbulence) and the downward pass of the implicit vertical diffusion, whose
  !> surface terms are passed to the surface model.
  subroutine physics_driver_down(is, ie, js, je, &
                                 Time_prev, Time, Time_next, &
                                 lat, lon, area, &
                                 p_half, p_full, z_half, z_full, &
                                 u, v, t, q, r, um, vm, tm, qm, rm, &
                                 frac_land, rough_mom, &
                                 albedo, t_surf_rad, &
                                 u_star, b_star, q_star, &
                                 dtau_du, dtau_dv, tau_x, tau_y, &
                                 udt, vdt, tdt, qdt, rdt, &
                                 flux_sw, flux_lw, gust, &
                                 Surf_diff, &
                                 mask, kbot, &
                                 diffm, difft)

    integer, intent(in)             :: is, ie, js, je
    !! starting and ending subdomain i, j indices of the physics window
    type(time_type), intent(in)             :: Time_prev, Time, &
                                               Time_next
    !! times of the previous level (of `um`, `vm`, `tm`, `qm`, `rm`), the current level (of
    !! `u`, `v`, `t`, `q`, `r`) and the next level (used for the diagnostics)
    real, dimension(:, :), intent(in)             :: lat, lon, area
    !! latitudes and longitudes of the model points [rad]; grid box area [m2] (not used)
    real, dimension(:, :, :), intent(in)             :: p_half, p_full, &
                                                        z_half, z_full, &
                                                        u, v, t, q, &
                                                        um, vm, tm, qm
    !! `p_half`, `p_full`: pressure at half and full levels [Pa]; `z_half`, `z_full`: height
    !! at half and full levels [m]; `u`, `v`, `t`, `q`: zonal and meridional wind [m/s],
    !! temperature [K] and specific humidity [kg/kg] at the current time level; `um`, `vm`,
    !! `tm`, `qm`: the same at the previous time level
    real, dimension(:, :, :, :), intent(inout)          :: r  !! tracers at the current time level
    real, dimension(:, :, :, :), intent(inout)          :: rm  !! tracers at the previous time level
    real, dimension(:, :), intent(in)             :: frac_land, &
                                                     rough_mom, &
                                                     albedo, t_surf_rad, &
                                                     u_star, b_star, &
                                                     q_star, dtau_du, dtau_dv
    !! surface fields: land fraction; roughness length for momentum [m]; albedo; radiative
    !! surface temperature [K]; friction velocity [m/s]; buoyancy scale [m/s2]; moisture
    !! scale [kg/kg]; derivatives of the zonal and meridional surface stress with respect to
    !! the lowest-level wind [kg/m2/s]
    real, dimension(:, :), intent(inout)          :: tau_x, tau_y  !! zonal and meridional surface stress [N/m2]
    real, dimension(:, :, :), intent(inout)          :: udt, vdt, tdt, qdt
    !! tendencies of the zonal and meridional wind [m/s2], temperature [K/s] and specific
    !! humidity [kg/kg/s]
    real, dimension(:, :, :, :), intent(inout)          :: rdt  !! tracer tendencies
    real, dimension(:, :), intent(out)            :: flux_sw, flux_lw, gust
    !! `flux_sw`: net downward shortwave flux at the surface [W/m2]; `flux_lw`: downward
    !! longwave flux at the surface [W/m2]; `gust`: gustiness [m/s]
    type(surf_diff_type), intent(inout)          :: Surf_diff  !! surface data of the implicit vertical diffusion
    real, dimension(:, :, :), intent(in), optional :: mask
    !! mask of the levels with data (0: below ground, 1: data); present together with `kbot`
    integer, dimension(:, :), intent(in), optional :: kbot  !! lowest level with data
    real, dimension(:, :, :), intent(out), optional :: diffm, difft
    !! diffusion coefficients for momentum and temperature [m2/s]

!---------------------------------------------------------------------
!    local variables:

    real, dimension(size(u, 1), size(u, 2), size(u, 3)) :: diff_t_vert, &
                                                           diff_m_vert
    real, dimension(size(u, 1), size(u, 2))           :: z_pbl
    integer          ::    sec, day
    real             ::    dt, alpha, dt2

!---------------------------------------------------------------------
!   local variables:
!
!      diff_t_vert     vertical diffusion coefficient for temperature
!                      calculated on the current step
!      diff_m_vert     vertical diffusion coefficient for momentum
!                      calculated on the current step
!      z_pbl           height of planetary boundary layer
!      sec, day        second and day components of the time_type
!                      variable
!      dt              model physics time step [ seconds ]
!      alpha           ratio of physics time step to diffusion-smoothing
!                      time scale
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    verify that the module is initialized.
!---------------------------------------------------------------------
    if (.not. module_is_initialized) then
      call error_mesg('physics_driver_mod', &
                      'module has not been initialized', FATAL)
    end if

!---------------------------------------------------------------------
!    check the size of the input arguments. this is only done on the
!    first call to physics_driver_down.
!---------------------------------------------------------------------
    if (do_check_args) call check_args &
      (lat, lon, area, p_half, p_full, z_half, z_full, &
       u, v, t, q, r, um, vm, tm, qm, rm, &
       udt, vdt, tdt, qdt, rdt)

!---------------------------------------------------------------------
!    compute the physics time step (from tau-1 to tau+1).
!---------------------------------------------------------------------
    call get_time(Time_next - Time_prev, sec, day)
    dt = real(sec + day*86400)

    flux_sw = 0.0
    flux_lw = 0.0

    call mpp_clock_begin(radiation_clock)
    call radiation_down(is, js, Time, Time_next, lat, lon, p_full, p_half, z_full, z_half, &
                        t, q, t_surf_rad, albedo, tdt, flux_sw, flux_lw)
    call mpp_clock_end(radiation_clock)

!----------------------------------------------------------------------
!    Held-Suarez forcing, computed from the previous time level
!----------------------------------------------------------------------
    if (do_held_suarez) then
      call held_suarez_forcing(is, js, Time_next, lat, p_full, p_half, &
                               um, vm, tm, udt, vdt, tdt)
    end if
!----------------------------------------------------------------------
!    artificial local heating if required
!----------------------------------------------------------------------
    if (do_local_heating) then
      call local_heating(is, js, Time, lon, lat, p_full, tdt)
    end if

!----------------------------------------------------------------------
!    call damping_driver to calculate the various model dampings that
!    are desired.
!----------------------------------------------------------------------
    if (do_damping) then
      call mpp_clock_begin(damping_clock)
      call damping_driver(is, js, lat, Time_next, dt, &
                          p_full, p_half, z_full, z_half, &
                          um, vm, tm, qm, rm(:, :, :, 1:ntp), &
                          udt, vdt, tdt, qdt, rdt, &
                          mask=mask, kbot=kbot)
      call mpp_clock_end(damping_clock)
    end if

!---------------------------------------------------------------------
!    call vert_turb_driver to calculate diffusion coefficients. save
!    the planetary boundary layer height on return.
!---------------------------------------------------------------------
    if (do_boundary_layer) then
      call mpp_clock_begin(turb_clock)
      call vert_turb_driver(is, js, Time_next, dt, &
                            p_half, p_full, z_half, z_full, u_star, &
                            b_star, u, v, t, q, um, vm, tm, qm, &
                            udt, vdt, tdt, qdt, &
                            diff_t_vert, diff_m_vert, gust, z_pbl, &
                            mask=mask, kbot=kbot)
      call mpp_clock_end(turb_clock)
    else
      gust = 0.0
    end if

!-----------------------------------------------------------------------
!    process any tracer fields.
!-----------------------------------------------------------------------
    call mpp_clock_begin(tracer_clock)
    call atmos_tracer_driver(is, ie, js, je, Time, lon, lat, &
                             frac_land, p_half, p_full, r, u, v, t, &
                             q, u_star, rdt, rm, dt, z_half, &
                             z_full, t_surf_rad, albedo, &
                             Time_next, kbot)
    call mpp_clock_end(tracer_clock)

!-----------------------------------------------------------------------
!    optionally use an implicit calculation of the vertical diffusion
!    coefficients.
!
!    the vertical diffusion coefficients are solved using an implicit
!    solution to the following equation:
!
!    dK/dt   = - ( K - K_cur) / tau_diff
!
!    where K         = diffusion coefficient
!          K_cur     = diffusion coefficient diagnosed from current
!                      time steps' state
!          tau_diff  = time scale for adjustment
!
!    in the code below alpha = dt / tau_diff
!---------------------------------------------------------------------
    if (do_boundary_layer) then
    if (diffusion_smooth) then
      call get_time(Time_next - Time, sec, day)
      dt2 = real(sec + day*86400)
      alpha = dt2/tau_diff
      diff_m(is:ie, js:je, :) = (diff_m(is:ie, js:je, :) + &
                                 alpha*diff_m_vert(:, :, :))/ &
                                (1.+alpha)
      where (diff_m(is:ie, js:je, :) < diff_min)
        diff_m(is:ie, js:je, :) = 0.0
      end where
      diff_t(is:ie, js:je, :) = (diff_t(is:ie, js:je, :) + &
                                 alpha*diff_t_vert(:, :, :))/ &
                                (1.+alpha)
      where (diff_t(is:ie, js:je, :) < diff_min)
        diff_t(is:ie, js:je, :) = 0.0
      end where
    else
      diff_t(is:ie, js:je, :) = diff_t_vert
      diff_m(is:ie, js:je, :) = diff_m_vert
    end if

!-----------------------------------------------------------------------
!    call vert_diff_driver_down to calculate the first pass atmos-
!    pheric vertical diffusion.
!-----------------------------------------------------------------------
    call mpp_clock_begin(diff_down_clock)
    call vert_diff_driver_down(is, js, Time_next, dt, p_half, &
                               p_full, z_full, &
                               diff_m(is:ie, js:je, :), &
                               diff_t(is:ie, js:je, :), &
                               um, vm, tm, qm, rm(:, :, :, 1:ntp), &
                               dtau_du, dtau_dv, tau_x, tau_y, &
                               udt, vdt, tdt, qdt, rdt, &
                               Surf_diff, &
                               mask=mask, kbot=kbot)
    end if ! do_boundary_layer

!---------------------------------------------------------------------
!    if desired, return diff_m and diff_t to calling routine.
!-----------------------------------------------------------------------
    if (present(difft)) then
      difft = diff_t(is:ie, js:je, :)
    end if
    if (present(diffm)) then
      diffm = diff_m(is:ie, js:je, :)
    end if

    call mpp_clock_end(diff_down_clock)

  end subroutine physics_driver_down

!#######################################################################
  !> Completes the vertical diffusion and computes the moist processes and the convective
  !> gustiness.
  subroutine physics_driver_up(is, ie, js, je, &
                               Time_prev, Time, Time_next, &
                               lat, lon, area, &
                               p_half, p_full, z_half, z_full, &
                               omega, &
                               u, v, t, q, r, um, vm, tm, qm, rm, &
                               frac_land, &
                               udt, vdt, tdt, qdt, rdt, &
                               Surf_diff, &
                               lprec, fprec, gust, &
                               mask, kbot)

    integer, intent(in)             :: is, ie, js, je
    !! starting and ending subdomain i, j indices of the physics window
    type(time_type), intent(in)             :: Time_prev, Time, &
                                               Time_next
    !! times of the previous level (of `um`, `vm`, `tm`, `qm`, `rm`), the current level (of
    !! `u`, `v`, `t`, `q`, `r`) and the next level (used for the diagnostics)
    real, dimension(:, :), intent(in)             :: lat, lon, area
    !! latitudes and longitudes of the model points [rad]; grid box area [m2] (not used)
    real, dimension(:, :, :), intent(in)             :: p_half, p_full, &
                                                        omega, &
                                                        z_half, z_full, &
                                                        u, v, t, q, &
                                                        um, vm, tm, qm
    !! `p_half`, `p_full`: pressure at half and full levels [Pa]; `omega`: vertical pressure
    !! velocity [Pa/s]; `z_half`, `z_full`: height at half and full levels [m]; `u`, `v`,
    !! `t`, `q`: zonal and meridional wind [m/s], temperature [K] and specific humidity
    !! [kg/kg] at the current time level; `um`, `vm`, `tm`, `qm`: the same at the previous
    !! time level
    real, dimension(:, :, :, :), intent(in)             :: r, rm
    !! tracers at the current and previous time levels
    real, dimension(:, :), intent(in)             :: frac_land  !! land fraction
    real, dimension(:, :, :), intent(inout)          :: udt, vdt, tdt, qdt
    !! tendencies of the zonal and meridional wind [m/s2], temperature [K/s] and specific
    !! humidity [kg/kg/s]
    real, dimension(:, :, :, :), intent(inout)          :: rdt  !! tracer tendencies
    type(surf_diff_type), intent(inout)          :: Surf_diff  !! surface data of the implicit vertical diffusion
    real, dimension(:, :), intent(out)            :: lprec, fprec
    !! liquid and frozen precipitation rates [kg/m2/s]
    real, dimension(:, :), intent(inout)          :: gust
    !! gustiness [m/s]; the convective gustiness is added on return
    real, dimension(:, :, :), intent(in), optional :: mask
    !! mask of the levels with data (0: below ground, 1: data); present together with `kbot`
    integer, dimension(:, :), intent(in), optional :: kbot  !! lowest level with data

!--------------------------------------------------------------------
!   local variables:

    real, dimension(size(u, 1), size(u, 2))            :: gust_cv
    integer :: sec, day
    real    :: dt

!---------------------------------------------------------------------
!   local variables:
!
!        gust_cv
!        sec, day         second and day components of the time_type
!                         variable
!        dt               physics time step [ seconds ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    verify that the module is initialized.
!---------------------------------------------------------------------
    if (.not. module_is_initialized) then
      call error_mesg('physics_driver_mod', &
                      'module has not been initialized', FATAL)
    end if

!---------------------------------------------------------------------
!    compute the physics time step (from tau-1 to tau+1).
!---------------------------------------------------------------------
    call get_time(Time_next - Time_prev, sec, day)
    dt = real(sec + day*86400)

!------------------------------------------------------------------
!    call vert_diff_driver_up to complete the vertical diffusion
!    calculation.
!------------------------------------------------------------------
    if (do_boundary_layer) then
      call mpp_clock_begin(diff_up_clock)
! XXX df's version of vert_diff_driver_up requires t for one of his new diagnostic fields
      call vert_diff_driver_up(is, js, Time_next, dt, p_half, &
                               Surf_diff, tdt, qdt, mask=mask, &
                               kbot=kbot, t=t)
      call mpp_clock_end(diff_up_clock)
    end if

!-----------------------------------------------------------------------
!    if the fms integration path is being followed, call moist processes
!    to compute moist physics, including convection and processes
!    involving condenstion.
!-----------------------------------------------------------------------
    if (do_moist_physics) then
      call mpp_clock_begin(moist_processes_clock)
      call moist_processes(is, ie, js, je, Time_next, dt, frac_land, &
                           p_half, p_full, z_half, z_full, omega, &
                           diff_t(is:ie, js:je, :), &
                           t, q, r, u, v, tm, qm, rm, um, vm, &
                           tdt, qdt, rdt, udt, vdt, &
                           lprec, fprec, &
                           gust_cv, area, lat, mask=mask, kbot=kbot)
      call mpp_clock_end(moist_processes_clock)

!---------------------------------------------------------------------
!    add the convective gustiness effect to that previously obtained
!    from non-convective parameterizations.
!---------------------------------------------------------------------
      gust = sqrt(gust*gust + gust_cv*gust_cv)
    else
      lprec = 0.0
      fprec = 0.0
    end if

!-----------------------------------------------------------------------

  end subroutine physics_driver_up

!#######################################################################
  !> Writes the restart file `RESTART/physics_driver.res.nc` and finalizes the physics
  !> modules.
  subroutine physics_driver_end(Time)

    type(time_type), intent(in) :: Time  !! current time

!---------------------------------------------------------------------
!   local variable:

    type(restart_file_type) :: rst
!---------------------------------------------------------------------
!    verify that the module is initialized.
!---------------------------------------------------------------------
    if (.not. module_is_initialized) then
      call error_mesg('physics_driver_mod', &
                      'module has not been initialized', FATAL)
    end if

    if (mpp_pe() == mpp_root_pe()) then
      call error_mesg('physics_driver_mod', 'Writing netCDF formatted restart file: RESTART/physics_driver.res.nc', NOTE)
    end if
    call open_restart_write(rst, 'RESTART/physics_driver.res.nc', domain)
    call write_restart_field(rst, 'vers', real(restart_versions(size(restart_versions(:)))))
    !--------------------------------------------------------------------
    !    write out the data fields that are relevant for this experiment.
    !--------------------------------------------------------------------
    call write_restart_field(rst, 'diff_t', diff_t)
    call write_restart_field(rst, 'diff_m', diff_m)
    call close_restart(rst)
!--------------------------------------------------------------------
!    call the destructor routines for those modules who were initial-
!    ized from this module.
!--------------------------------------------------------------------
    call vert_turb_driver_end
    call vert_diff_driver_end
    call radiation_end
    if (do_held_suarez) call held_suarez_end
    call moist_processes_end
    call atmos_tracer_driver_end
    if (do_damping) call damping_driver_end

!---------------------------------------------------------------------
!    deallocate the module variables.
!---------------------------------------------------------------------
    deallocate (diff_t, diff_m)

!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
    module_is_initialized = .false.

!-----------------------------------------------------------------------

  end subroutine physics_driver_end

!#####################################################################

!#######################################################################

  !> Returns `.true.` if the atmosphere exchanges heat, moisture and momentum with the
  !> surface (`do_boundary_layer`); otherwise the coupler does not update the surface state.
  logical function surface_is_coupled()

    surface_is_coupled = do_boundary_layer

  end function surface_is_coupled

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################

  !> Checks that FMS was built with its GFDL physical constants (the default, CMake option
  !> `CONSTANTS=GFDL`), which MiMA's configurations and tuning assume. The GFS and GEOS
  !> sets differ, e.g. in `RADIUS` and `GRAV`.
  subroutine check_constants

    real, parameter :: tol = 1.e-12
    real, dimension(9) :: fms, gfdl

    fms = (/RADIUS, OMEGA, GRAV, RDGAS, RVGAS, KAPPA, HLV, HLF, STEFAN/)
    gfdl = (/6371.0e3, 7.292e-5, 9.80, 287.04, 461.50, 2./7., 2.500e6, 3.34e5, 5.6734e-8/)
    if (any(abs(fms - gfdl) > tol*abs(gfdl))) then
      call error_mesg('physics_driver_init', &
                      'FMS was not built with the GFDL physical constants MiMA '// &
                      'assumes: rebuild FMS with CONSTANTS=GFDL', FATAL)
    end if

  end subroutine check_constants

!#####################################################################

  !> Initializes the FMS saturation vapour pressure table and checks that it is the simple
  !> Clausius-Clapeyron form MiMA uses.
  !>
  !> The simple form (constant latent heat of vaporization, no ice) is computed by FMS only
  !> with `do_simple = .true.` in `sat_vapor_pres_nml`. The table is interpolated, so it is
  !> compared with the formula to a relative tolerance.
  subroutine check_sat_vapor_pres

    real, parameter :: tol = 1.e-5
    real, dimension(5), parameter :: temp = &
                                     (/180., 230., 273.16, 300., 330./)
    real, dimension(5) :: es_table, es_simple

    call sat_vapor_pres_init
    call lookup_es(temp, es_table)
    es_simple = ES0*610.78*exp(-HLV/RVGAS*(1./temp - 1./TFREEZE))
    if (any(abs(es_table - es_simple) > tol*es_simple)) then
      call error_mesg('physics_driver_init', &
                      'the saturation vapour pressure is not the simple form '// &
                      'MiMA uses: set do_simple = .true. in sat_vapor_pres_nml', FATAL)
    end if

  end subroutine check_sat_vapor_pres

!#####################################################################
  !> Reads the diffusion coefficients from `INPUT/physics_driver.res.nc`, or sets them to
  !> zero if there is no restart file.
  subroutine read_restart_nc

    type(restart_file_type) :: rst

    if (open_restart_read(rst, 'INPUT/physics_driver.res.nc', domain)) then
      if (mpp_pe() == mpp_root_pe()) call mpp_error('physics_driver_mod', &
                                                    'Reading NetCDF formatted restart file: INPUT/physics_driver.res.nc', NOTE)
      call read_restart_field(rst, 'diff_t', diff_t)
      call read_restart_field(rst, 'diff_m', diff_m)
      call close_restart(rst)
    else
      diff_t = 0.0
      diff_m = 0.0
    end if

  end subroutine read_restart_nc

!#####################################################################
  !> Checks that the input arrays of `physics_driver_down` have consistent sizes.
  subroutine check_args(lat, lon, area, p_half, p_full, z_half, z_full, &
                        u, v, t, q, r, um, vm, tm, qm, rm, &
                        udt, vdt, tdt, qdt, rdt, mask, kbot)

!----------------------------------------------------------------------
!    check_args determines if the input arrays to physics_driver_down
!    are of a consistent size.
!-----------------------------------------------------------------------

    real, dimension(:, :), intent(in)          :: lat, lon, area
    real, dimension(:, :, :), intent(in)          :: p_half, p_full, &
                                                     z_half, z_full, &
                                                     u, v, t, q, um, vm, &
                                                     tm, qm
    real, dimension(:, :, :, :), intent(in)          :: r, rm
    real, dimension(:, :, :), intent(in)          :: udt, vdt, tdt, qdt
    real, dimension(:, :, :, :), intent(in)          :: rdt
    real, dimension(:, :, :), intent(in), optional :: mask
    integer, dimension(:, :), intent(in), optional :: kbot

!-----------------------------------------------------------------------
!   intent(in) variables:
!
!      lat            latitude of model points [ radians ]
!      lon            longitude of model points [ radians ]
!      area           grid box area - currently not used [ m**2 ]
!      p_half         pressure at half levels (offset from t,q,u,v,r)
!                     [ Pa ]
!      p_full         pressure at full levels [ Pa }
!      z_half         height at half levels [ m ]
!      z_full         height at full levels [ m ]
!      u              zonal wind at current time step [ m / s ]
!      v              meridional wind at current time step [ m / s ]
!      t              temperature at current time step [ deg k ]
!      q              specific humidity at current time step  kg / kg ]
!      r              multiple 3d tracer fields at current time step
!      um,vm          zonal and meridional wind at previous time step
!      tm,qm          temperature and specific humidity at previous
!                     time step
!      rm             multiple 3d tracer fields at previous time step
!      udt            zonal wind tendency [ m / s**2 ]
!      vdt            meridional wind tendency [ m / s**2 ]
!      tdt            temperature tendency [ deg k / sec ]
!      qdt            specific humidity tendency
!                     [  kg vapor / kg air / sec ]
!      rdt            multiple tracer tendencies [ unit / unit / sec ]
!
!   intent(in), optional:
!
!       mask        mask that designates which levels do not have data
!                   present (i.e., below ground); 0.=no data, 1.=data
!       kbot        lowest level which has data
!                   note:  both mask and kbot must be present together.
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:

    integer ::  id, jd, kd  ! model dimensions on the processor
    integer ::  ierr        ! error flag

!--------------------------------------------------------------------
!    define the sizes that the arrays should be.
!--------------------------------------------------------------------
    id = size(u, 1)
    jd = size(u, 2)
    kd = size(u, 3)

!--------------------------------------------------------------------
!    check the dimensions of each input array. if they are incompat-
!    ible in size with the standard, the error flag is set to so
!    indicate.
!--------------------------------------------------------------------
    ierr = 0
    ierr = ierr + check_dim(lat, 'lat', id, jd)
    ierr = ierr + check_dim(lon, 'lon', id, jd)
    ierr = ierr + check_dim(area, 'area', id, jd)

    ierr = ierr + check_dim(p_half, 'p_half', id, jd, kd + 1)
    ierr = ierr + check_dim(p_full, 'p_full', id, jd, kd)
    ierr = ierr + check_dim(z_half, 'z_half', id, jd, kd + 1)
    ierr = ierr + check_dim(z_full, 'z_full', id, jd, kd)

    ierr = ierr + check_dim(u, 'u', id, jd, kd)
    ierr = ierr + check_dim(v, 'v', id, jd, kd)
    ierr = ierr + check_dim(t, 't', id, jd, kd)
    ierr = ierr + check_dim(q, 'q', id, jd, kd)
    ierr = ierr + check_dim(um, 'um', id, jd, kd)
    ierr = ierr + check_dim(vm, 'vm', id, jd, kd)
    ierr = ierr + check_dim(tm, 'tm', id, jd, kd)
    ierr = ierr + check_dim(qm, 'qm', id, jd, kd)

    ierr = ierr + check_dim(udt, 'udt', id, jd, kd)
    ierr = ierr + check_dim(vdt, 'vdt', id, jd, kd)
    ierr = ierr + check_dim(tdt, 'tdt', id, jd, kd)
    ierr = ierr + check_dim(qdt, 'qdt', id, jd, kd)

    if (nt > 0) then
      ierr = ierr + check_dim(r, 'r', id, jd, kd, nt)
      ierr = ierr + check_dim(rm, 'rm', id, jd, kd, nt)
    end if
    if (ntp > 0) then
      ierr = ierr + check_dim(rdt, 'rdt', id, jd, kd, ntp)
    end if

!--------------------------------------------------------------------
!    if any problems were detected, exit with an error message.
!--------------------------------------------------------------------
    if (ierr > 0) then
      call error_mesg('physics_driver_mod', 'bad dimensions', FATAL)
    end if

!--------------------------------------------------------------------
!    set a flag to indicate that this check was done and need not be
!    done again.
!--------------------------------------------------------------------
    do_check_args = .false.

!-----------------------------------------------------------------------

  end subroutine check_args

!#######################################################################
  function check_dim_2d(data, name, id, jd) result(ierr)

!--------------------------------------------------------------------
!    check_dim_2d compares the size of two-dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!--------------------------------------------------------------------

    real, intent(in), dimension(:, :) :: data
    character(len=*), intent(in)        :: name
    integer, intent(in)                 :: id, jd
    integer                             :: ierr

!---------------------------------------------------------------------
!  intent(in) variables:
!
!     data        array to be checked
!     name        name associated with array to be checked
!     id, jd      expected i and j dimensions
!
!  result variable:
!
!     ierr        set to 0 if ok, otherwise is a count of the number
!                 of incompatible dimensions
!
!--------------------------------------------------------------------

    ierr = 0
    if (size(data, 1) /= id) then
      call error_mesg('physics_driver_mod', &
                      'dimension 1 of argument '// &
                      name(1:len_trim(name))//' has wrong size.', NOTE)
      ierr = ierr + 1
    end if
    if (size(data, 2) /= jd) then
      call error_mesg('physics_driver_mod', &
                      'dimension 2 of argument '// &
                      name(1:len_trim(name))//' has wrong size.', NOTE)
      ierr = ierr + 1
    end if

!----------------------------------------------------------------------

  end function check_dim_2d

!#######################################################################
  function check_dim_3d(data, name, id, jd, kd) result(ierr)

!--------------------------------------------------------------------
!    check_dim_3d compares the size of thr1eedimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!--------------------------------------------------------------------

    real, intent(in), dimension(:, :, :) :: data
    character(len=*), intent(in)          :: name
    integer, intent(in)                   :: id, jd, kd
    integer ierr

!---------------------------------------------------------------------
!  intent(in) variables:
!
!     data        array to be checked
!     name        name associated with array to be checked
!     id, jd,kd   expected i, j and k dimensions
!
!  result variable:
!
!     ierr        set to 0 if ok, otherwise is a count of the number
!                 of incompatible dimensions
!
!--------------------------------------------------------------------

    ierr = 0
    if (size(data, 1) /= id) then
      call error_mesg('physics_driver_mod', &
                      'dimension 1 of argument '// &
                      name(1:len_trim(name))//' has wrong size.', NOTE)
      ierr = ierr + 1
    end if
    if (size(data, 2) /= jd) then
      call error_mesg('physics_driver_mod', &
                      'dimension 2 of argument '// &
                      name(1:len_trim(name))//' has wrong size.', NOTE)
      ierr = ierr + 1
    end if
    if (size(data, 3) /= kd) then
      call error_mesg('physics_driver_mod', &
                      'dimension 3 of argument '// &
                      name(1:len_trim(name))//' has wrong size.', NOTE)
      ierr = ierr + 1
    end if

!---------------------------------------------------------------------

  end function check_dim_3d

!#######################################################################
  function check_dim_4d(data, name, id, jd, kd, nt) result(ierr)

!--------------------------------------------------------------------
!    check_dim_4d compares the size of four dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!--------------------------------------------------------------------
    real, intent(in), dimension(:, :, :, :) :: data
    character(len=*), intent(in)            :: name
    integer, intent(in)                     :: id, jd, kd, nt
    integer                                 :: ierr

!---------------------------------------------------------------------
!  intent(in) variables:
!
!     data          array to be checked
!     name          name associated with array to be checked
!     id,jd,kd,nt   expected i, j and k dimensions
!
!  result variable:
!
!     ierr          set to 0 if ok, otherwise is a count of the number
!                   of incompatible dimensions
!
!--------------------------------------------------------------------

    ierr = 0
    if (size(data, 1) /= id) then
      call error_mesg('physics_driver_mod', &
                      'dimension 1 of argument '// &
                      name(1:len_trim(name))//' has wrong size.', NOTE)
      ierr = ierr + 1
    end if
    if (size(data, 2) /= jd) then
      call error_mesg('physics_driver_mod', &
                      'dimension 2 of argument '// &
                      name(1:len_trim(name))//' has wrong size.', NOTE)
      ierr = ierr + 1
    end if
    if (size(data, 3) /= kd) then
      call error_mesg('physics_driver_mod', &
                      'dimension 3 of argument '// &
                      name(1:len_trim(name))//' has wrong size.', NOTE)
      ierr = ierr + 1
    end if
    if (size(data, 4) /= nt) then
      call error_mesg('physics_driver_mod', &
                      'dimension 4 of argument '// &
                      name(1:len_trim(name))//' has wrong size.', NOTE)
      ierr = ierr + 1
    end if

!---------------------------------------------------------------------

  end function check_dim_4d

!#######################################################################

end module physics_driver_mod
