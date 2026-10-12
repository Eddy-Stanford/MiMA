!> The exchange of fluxes between the atmosphere and the surface components (ocean and land).
!>
!> MiMA's atmosphere, ocean and land all use the model grid, so the exchange is a set of
!> point-by-point operations, organised as in the GFDL coupler:
!>
!> 1. `sfc_boundary_layer`: the surface fluxes and their derivatives at the current surface
!>    temperature (`surface_flux_mod`), with the roughness lengths and the background albedo
!>    (`surface_albedo_mod`) set at initialization;
!> 2. `flux_down_from_atmos`: after the down sweep of the atmosphere's vertical diffusion,
!>    completes the implicit coupling for the lowest level and hands the fluxes and
!>    derivatives to the ocean and land models (`surface_exchange_type`), together with the
!>    local surface heating of `local_heating_nml`;
!> 3. `flux_up_to_atmos`: after the ocean and land have stepped their surface temperature,
!>    updates the fluxes to the new surface temperature and returns the changes of the
!>    lowest-level temperature and humidity for the up sweep.
!>
!> With a prescribed or fixed SST (`ocean_model_nml`), `flux_up_to_atmos` also sends the
!> implied ocean heat flux `ocean/qflux_implied`: the Q-flux with which a slab ocean of the
!> same heat capacity would reproduce the SST change of the step. Its time mean (or monthly
!> means) can be read back as `qflux_file` (`qflux_nml`).
!>
!> The surface temperature seen by the atmosphere is the land temperature on land points
!> and the ocean temperature elsewhere. The module also sends the merged surface
!> diagnostics (module `flux`), including the reference-height diagnostics `t_ref`,
!> `rh_ref`, `u_ref` and `v_ref`.
!>
!> Namelist: `flux_exchange_nml`
!> ([namelist reference](https://eddy-stanford.github.io/MiMA/Parameters/#flux_exchange_nml)).
module flux_exchange_mod

  use atmos_model_mod, only: atmos_data_type, land_ice_atmos_boundary_type

  use surface_flux_mod, only: surface_flux, surface_flux_init
  use mima_monin_obukhov_mod, only: mo_profile
  use sat_vapor_pres_mod, only: compute_qs

  use mima_diag_integral_mod, only: diag_integral_field_init, &
                                    sum_diag_integral_field

  use fms_mod, only: input_nml_file, check_nml_error, &
                     error_mesg, FATAL, mpp_pe, mpp_root_pe, stdlog
  use fms2_io_mod, only: file_exists
  use restart_file_mod, only: restart_file_type, open_restart_read, close_restart, &
                              read_restart_field, restart_field_exists

  use diag_manager_mod, only: register_diag_field, send_data

  use time_manager_mod, only: time_type, operator(-)
  use field_manager_mod, only: MODEL_ATMOS
  use tracer_manager_mod, only: get_number_tracers, query_method

  use constants_mod, only: cp_air, hlv, hlf

!mj local surface heating
  use physics_driver_mod, only: do_local_heating, surface_is_coupled
  use local_heating_mod, only: horizontal_heating, ngauss, hamp, pcenter

  use surface_types_mod, only: surface_grid_type, surface_exchange_type
  use namelist_check_mod, only: nml_has_group
  use surface_albedo_mod, only: surface_albedo_init
  use ocean_model_mod, only: ocean_data_type, ocean_model_init, ocean_model_end
  use land_model_mod, only: land_data_type, land_model_read_nml_and_mask, &
                            land_model_init_state, land_model_read_restart, &
                            land_model_roughness, land_model_boundary, land_model_end

  implicit none
  private

  public :: flux_exchange_init, sfc_boundary_layer, flux_down_from_atmos, &
            flux_up_to_atmos, flux_exchange_end

!-----------------------------------------------------------------------
  character(len=4), parameter :: mod_name = 'flux'

  integer :: id_drag_moist, id_drag_heat, id_drag_mom, &
             id_rough_moist, id_rough_heat, id_rough_mom, &
             id_u_star, id_b_star, id_u_flux, id_v_flux, &
             id_t_atm, id_u_atm, id_v_atm, id_wind, &
             id_t_ref, id_rh_ref, id_u_ref, id_v_ref, &
             id_del_h, id_del_m, id_del_q, id_albedo
! Sent by flux_up_to_atmos, which runs only with a coupled surface: registered only then
  integer :: id_t_surf = 0, id_t_flux = 0, id_q_flux = 0, id_r_flux = 0, id_heat = 0, &
             id_entrop_evap = 0, id_entrop_shflx = 0, id_entrop_lwflx = 0
  integer :: id_qflux_implied = 0

  logical :: do_init = .true.
  logical :: do_surface_heating = .false.
  logical :: need_fluxes = .true.  ! sfc_boundary_layer computes the fluxes (see there)

!-----------------------------------------------------------------------

  real :: z_ref_heat = 2.  !! [m] reference height of the diagnostics `t_ref` and `rh_ref`
  real :: z_ref_mom = 10.  !! [m] reference height of the diagnostics `u_ref` and `v_ref`

  namelist /flux_exchange_nml/ z_ref_heat, z_ref_mom

!-----------------------------------------------------------------------
!---- module storage ------

  type(surface_grid_type), save :: grid
  real :: dt_step

  real, allocatable, dimension(:, :) :: albedo_bg, rough_mom, rough_heat, rough_moist

! Static fields for surface_flux and the local surface heating, set by flux_exchange_init
  real, allocatable, dimension(:, :) :: u_surf, v_surf  ! surface velocity: 0 (no currents)
  logical, allocatable, dimension(:, :) :: all_points, no_points  ! .true. and .false. everywhere
  real, allocatable, dimension(:, :) :: lon2d, lat2d  ! grid-box centres [rad] (surface heating only)

  real, allocatable, dimension(:, :) :: e_t_n, f_t_delt_n, &
                                        e_q_n, f_q_delt_n

  real, allocatable, dimension(:, :) :: dhdt_surf, dedt_surf, &
                                        drdt_surf, dhdt_atm, dedq_atm, &
                                        flux_t, flux_q, flux_lw, flux_u, flux_v

contains

!#######################################################################

  !> Computes the surface fluxes and their derivatives at the current surface temperature
  !> (`surface_flux`) and passes the surface state to the atmosphere; sends the flux
  !> diagnostics.
  !>
  !> The fluxes and derivatives needed by `flux_down_from_atmos` and `flux_up_to_atmos` are
  !> kept in module storage. With the bucket hydrology, the land model first sets the
  !> evaporation efficiency `beta` and the snow albedo (`land_model_boundary`), and on the
  !> land points the evaporation, its derivatives and the moisture scale `q_star` are then
  !> multiplied by `beta`; the surface humidity of the diagnostics becomes
  !> `q_atm + beta*(q_surf - q_atm)`. The buoyancy in `surface_flux` still uses the saturated
  !> surface humidity (an approximation).
  !>
  !> Without a coupled surface (`surface_is_coupled`, e.g. Held-Suarez) the atmosphere uses
  !> only the surface temperature and albedo, and the fluxes are skipped unless a diagnostic
  !> of the module `flux` is requested or a tracer has a `dry_deposition` method (which uses
  !> `u_star`); `Boundary%u_star` and the other flux fields then stay 0.
  subroutine sfc_boundary_layer(dt, Time, Atm, Ocean, Land, Boundary)

    real, intent(in)  :: dt  !! time step [s]
    type(time_type), intent(in)  :: Time  !! current time
    type(atmos_data_type), intent(in)  :: Atm  !! atmospheric state at the lowest level
    type(ocean_data_type), intent(in) :: Ocean  !! the ocean state
    type(land_data_type), intent(inout) :: Land  !! the land state (on return, with the bucket, `beta`)
    type(land_ice_atmos_boundary_type), intent(inout) :: Boundary
    !! surface fields for the atmosphere: on return the surface temperature, albedo and
    !! roughness, the surface stress and its derivatives, `u_star`, `b_star` and the land
    !! fraction (0)

    real, dimension(size(Atm%t_bot, 1), size(Atm%t_bot, 2)) :: &
      q_star, q_surf, cd_q, cd_t, cd_m, wind, &
      dedq_surf  ! derivative of the evaporation with respect to q_surf (0, and not used)

    logical :: used

!-----------------------------------------------------------------------

    if (do_init) call error_mesg('sfc_boundary_layer', &
                                 'must call flux_exchange_init first', FATAL)

    dt_step = dt

! Boundary%rough_mom and, without the bucket, Boundary%albedo are static (flux_exchange_init)
    if (Land%do_bucket) then
      Boundary%albedo = albedo_bg
      call land_model_boundary(Land, Boundary%albedo)
    end if

    where (Land%mask)
      Boundary%t = Land%t_surf
    elsewhere
      Boundary%t = Ocean%t_surf
    end where

    if (.not. need_fluxes) return

! surface_flux reads q_surf only at its land points (there are none here, see below)
    q_surf = 0.0

! Every point, land included, is a saturated surface at Boundary%t for surface_flux: no
! point is passed as land, so t_ca and the input q_surf are not used and dedq_surf is 0
    call surface_flux(Atm%t_bot, Atm%q_bot, Atm%u_bot, Atm%v_bot, &
                      Atm%p_bot, Atm%z_bot, Atm%p_surf, &
                      Boundary%t, &  ! t_surf
                      Boundary%t, &  ! t_ca (not used)
                      q_surf, u_surf, v_surf, &
                      Boundary%rough_mom, rough_heat, rough_moist, Atm%gust, &
                      flux_t, flux_q, flux_lw, flux_u, flux_v, &
                      cd_m, cd_t, cd_q, wind, Boundary%u_star, Boundary%b_star, q_star, &
                      dhdt_surf, dedt_surf, dedq_surf, drdt_surf, &
                      dhdt_atm, dedq_atm, Boundary%dtaudu, Boundary%dtaudv, &
                      no_points, &   ! land: none
                      Ocean%mask, &  ! seawater: the ocean points (used only with ncar_ocean_flux)
                      all_points)    ! avail: all points

! Bucket: the land evaporates beta times the potential (saturated) rate
    if (Land%do_bucket) then
      where (Land%mask)
        flux_q = flux_q*Land%beta
        dedt_surf = dedt_surf*Land%beta
        dedq_atm = dedq_atm*Land%beta
        q_star = q_star*Land%beta
        q_surf = Atm%q_bot + Land%beta*(q_surf - Atm%q_bot)
      end where
    end if

    Boundary%u_flux = flux_u
    Boundary%v_flux = flux_v

!=======================================================================
!-------------------- diagnostics section ------------------------------

    if (id_wind > 0) used = send_data(id_wind, wind, Time)
    if (id_drag_moist > 0) used = send_data(id_drag_moist, cd_q, Time)
    if (id_drag_heat > 0) used = send_data(id_drag_heat, cd_t, Time)
    if (id_drag_mom > 0) used = send_data(id_drag_mom, cd_m, Time)
    if (id_rough_heat > 0) used = send_data(id_rough_heat, rough_heat, Time)
    if (id_rough_mom > 0) used = send_data(id_rough_mom, Boundary%rough_mom, Time)
    if (id_rough_moist > 0) used = send_data(id_rough_moist, rough_moist, Time)
    if (id_u_star > 0) used = send_data(id_u_star, Boundary%u_star, Time)
    if (id_b_star > 0) used = send_data(id_b_star, Boundary%b_star, Time)
    if (id_t_atm > 0) used = send_data(id_t_atm, Atm%t_bot, Time)
    if (id_u_atm > 0) used = send_data(id_u_atm, Atm%u_bot, Time)
    if (id_v_atm > 0) used = send_data(id_v_atm, Atm%v_bot, Time)
    if (id_albedo > 0) used = send_data(id_albedo, Boundary%albedo, Time)
    if (id_u_flux > 0) used = send_data(id_u_flux, flux_u, Time)
    if (id_v_flux > 0) used = send_data(id_v_flux, flux_v, Time)

    if (id_t_ref > 0 .or. id_rh_ref > 0 .or. id_u_ref > 0 .or. id_v_ref > 0 .or. &
        id_del_h > 0 .or. id_del_m > 0 .or. id_del_q > 0) &
      call reference_height_diagnostics(Time, Atm, Boundary%t, q_surf, &
                                        u_surf, v_surf, Boundary%rough_mom, rough_heat, &
                                        rough_moist, Boundary%u_star, Boundary%b_star, q_star)

!=======================================================================

  end subroutine sfc_boundary_layer

!#######################################################################

  !> Sends the diagnostics at the reference heights `z_ref_heat` and `z_ref_mom`.
  subroutine reference_height_diagnostics(Time, Atm, t_surf, q_surf, &
                                          u_surf, v_surf, rough_mom, &
                                          rough_heat, rough_moist, &
                                          u_star, b_star, q_star)

! Diagnostics at the reference heights z_ref_heat (t_ref, rh_ref) and
! z_ref_mom (u_ref, v_ref) above the surface, computed as in the FMS/AM2
! flux_exchange. Between the surface and the lowest model level (height
! Atm%z_bot) a quantity f is given by the Monin-Obukhov profile, so that
!
!     f(z_ref) = f_surf + (f_atm - f_surf) * del_f,
!     del_f    = (f(z_ref) - f_surf) / (f(z_atm) - f_surf),
!
! where del_m (winds), del_h (temperature) and del_q (specific humidity)
! come from mo_profile, for the same roughness lengths and u_star, b_star
! as the surface fluxes. The surface values are t_surf, the surface
! specific humidity q_surf from surface_flux, and zero wind. Then
!
!     rh_ref = 100 * q_ref / q_sat(t_ref, p_surf),
!     q_sat  = eps*e_s / (p_surf - (1-eps)*e_s),   eps = rdgas/rvgas,
!
! with e_s(t_ref) from the sat_vapor_pres_mod table (compute_qs; with
! do_simple = .true. the table is the simple Clausius-Clapeyron fit over
! liquid water). These fields do not affect the model state.

    type(time_type), intent(in) :: Time
    type(atmos_data_type), intent(in) :: Atm
    real, dimension(:, :), intent(in) :: t_surf, q_surf, u_surf, v_surf, &
                                         rough_mom, rough_heat, rough_moist, &
                                         u_star, b_star, q_star

    real, dimension(size(t_surf, 1), size(t_surf, 2)) :: del_m, del_h, del_q, &
                                                         del_unused, t_ref, q_ref, qs_ref
    logical :: used

! mo_profile returns all three factors at its first height (its second height
! argument is not used), so del_h and del_q are recomputed at z_ref_heat.
    call mo_profile(z_ref_mom, z_ref_heat, Atm%z_bot, &
                    rough_mom, rough_heat, rough_moist, &
                    u_star, b_star, q_star, del_m, del_h, del_q)
    if (z_ref_heat /= z_ref_mom) &
      call mo_profile(z_ref_heat, z_ref_heat, Atm%z_bot, &
                      rough_mom, rough_heat, rough_moist, &
                      u_star, b_star, q_star, del_unused, del_h, del_q)

    t_ref = t_surf + (Atm%t_bot - t_surf)*del_h
    if (id_t_ref > 0) used = send_data(id_t_ref, t_ref, Time)

    if (id_rh_ref > 0) then
      q_ref = q_surf + (Atm%q_bot - q_surf)*del_q
      call compute_qs(t_ref, Atm%p_surf, qs_ref)
      used = send_data(id_rh_ref, 100.*q_ref/qs_ref, Time)
    end if

    if (id_u_ref > 0) used = send_data(id_u_ref, u_surf + (Atm%u_bot - u_surf)*del_m, Time)
    if (id_v_ref > 0) used = send_data(id_v_ref, v_surf + (Atm%v_bot - v_surf)*del_m, Time)

    if (id_del_h > 0) used = send_data(id_del_h, del_h, Time)
    if (id_del_m > 0) used = send_data(id_del_m, del_m, Time)
    if (id_del_q > 0) used = send_data(id_del_q, del_q, Time)

  end subroutine reference_height_diagnostics

!#######################################################################

  !> Completes the implicit coupling with the atmosphere's vertical diffusion for the
  !> lowest level, and fills `Exch` for the ocean and land models.
  !>
  !> Called after the down sweep of the vertical diffusion. The sensible heat and moisture
  !> fluxes and their derivatives with respect to the surface temperature are updated with
  !> the lowest-level changes of the down sweep. The local surface heating of
  !> `local_heating_nml` is added to `Exch` (zero if there is none); it is evaluated at the
  !> start of the step (`Time - Atm%Time_step`), the time at which the physics evaluates
  !> the atmospheric part of the local heating.
  subroutine flux_down_from_atmos(Time, Atm, Exch)

    type(time_type), intent(in)  :: Time  !! time at the end of the step (for `Exch%Time`)
    type(atmos_data_type), intent(in)  :: Atm
    !! atmospheric state, including the surface radiation and precipitation and the data of
    !! the down sweep of the vertical diffusion (`Atm%Surf_diff`)
    type(surface_exchange_type), intent(inout) :: Exch  !! fields for the ocean and land models

    real, dimension(size(Atm%t_bot, 1), size(Atm%t_bot, 2)) :: &
      gamma, dtmass, delta_t, delta_q, dflux_t, dflux_q

    real    :: cp_inv

    flux_lw = Atm%flux_lw - flux_lw

    dtmass = Atm%Surf_Diff%dtmass
    delta_t = Atm%Surf_Diff%delta_t
    delta_q = Atm%Surf_Diff%delta_q
    dflux_t = Atm%Surf_Diff%dflux_t
    dflux_q = Atm%Surf_Diff%dflux_q

    cp_inv = 1.0/cp_air

    ! temperature

    gamma = 1./(1.0 - dtmass*(dflux_t + dhdt_atm*cp_inv))
    e_t_n = dtmass*dhdt_surf*cp_inv*gamma
    f_t_delt_n = (delta_t + dtmass*flux_t*cp_inv)*gamma

    flux_t = flux_t + dhdt_atm*f_t_delt_n
    dhdt_surf = dhdt_surf + dhdt_atm*e_t_n

! moisture

    gamma = 1./(1.0 - dtmass*(dflux_q + dedq_atm))
    e_q_n = dtmass*dedt_surf*gamma
    f_q_delt_n = (delta_q + dtmass*flux_q)*gamma

    flux_q = flux_q + dedq_atm*f_q_delt_n
    dedt_surf = dedt_surf + dedq_atm*e_q_n

! The local surface heating (Exch%extra_heating stays 0 without it), at the start of the
! step like the atmospheric part in physics_driver_down
    if (do_surface_heating) &
      call horizontal_heating(Time - Atm%Time_step, lon2d, lat2d, Exch%extra_heating)

    Exch%flux_lw = flux_lw
    Exch%flux_sw = Atm%flux_sw
    Exch%lprec = Atm%lprec
    Exch%fprec = Atm%fprec
    Exch%flux_t = flux_t
    Exch%flux_q = flux_q
    Exch%dhdt_surf = dhdt_surf
    Exch%dedt_surf = dedt_surf
    Exch%drdt_surf = drdt_surf
    Exch%dt = dt_step
    Exch%Time = Time

  end subroutine flux_down_from_atmos

!#######################################################################

  !> Updates the surface fluxes to the new surface temperature, returns the changes of the
  !> lowest-level temperature and humidity for the up sweep of the vertical diffusion, and
  !> sends the surface diagnostics.
  subroutine flux_up_to_atmos(Time, Atm, Ocean, Land, Exch, Boundary)

    type(time_type), intent(in)  :: Time  !! current time
    type(atmos_data_type), intent(in)  :: Atm  !! atmospheric state
    type(ocean_data_type), intent(in) :: Ocean  !! the ocean state, after its step
    type(land_data_type), intent(in) :: Land  !! the land state, after its step
    type(surface_exchange_type), intent(in) :: Exch  !! fields handed to the ocean and land models
    type(land_ice_atmos_boundary_type), intent(inout) :: Boundary
    !! on return, `dt_t` and `dt_q`: changes of the lowest-level temperature and specific
    !! humidity for the up sweep of the vertical diffusion

    real, dimension(size(Atm%t_bot, 1), size(Atm%t_bot, 2)) :: &
      dt_t_surf, t_surf, heat_capacity, &
      entrop_evap, entrop_shflx, entrop_lwflx, qflux_implied

    logical :: used

    where (Land%mask)
      dt_t_surf = Land%dt_t_surf
    elsewhere
      dt_t_surf = Ocean%dt_t_surf
    end where

    flux_t = flux_t + dt_t_surf*dhdt_surf
    flux_q = flux_q + dt_t_surf*dedt_surf
    flux_lw = flux_lw - dt_t_surf*drdt_surf
    Boundary%dt_t = f_t_delt_n + dt_t_surf*e_t_n
    Boundary%dt_q = f_q_delt_n + dt_t_surf*e_q_n

!=======================================================================
!-------------------- diagnostics section ------------------------------

! Implied ocean heat flux (prescribed or fixed SST): with the final fluxes of the step,
! cap*dT/dt = F_net + Q, where dT excludes the local surface heating
    if (id_qflux_implied > 0) then
      where (Ocean%mask)
        qflux_implied = Ocean%cap*(Ocean%dt_t_surf/Exch%dt - Exch%extra_heating) &
                        - (flux_lw + Exch%flux_sw - hlf*Exch%fprec - (flux_t + hlv*flux_q))
      elsewhere
        qflux_implied = 0.0
      end where
      used = send_data(id_qflux_implied, qflux_implied, Time)
    end if

    if (id_t_surf > 0 .or. id_entrop_evap > 0 .or. id_entrop_shflx > 0 .or. id_entrop_lwflx > 0) then
      where (Land%mask)
        t_surf = Land%t_surf
      elsewhere
        t_surf = Ocean%t_surf
      end where
    end if

    if (id_t_surf > 0) used = send_data(id_t_surf, t_surf, Time)
    if (id_t_flux > 0) used = send_data(id_t_flux, flux_t, Time)
    if (id_r_flux > 0) used = send_data(id_r_flux, flux_lw, Time)
    if (id_q_flux > 0) used = send_data(id_q_flux, flux_q, Time)
    if (id_heat > 0) then
      where (Land%mask)
        heat_capacity = Land%cap
      elsewhere
        heat_capacity = Ocean%cap
      end where
      used = send_data(id_heat, heat_capacity, Time)
    end if
    if (id_entrop_evap > 0) then
      entrop_evap = flux_q/t_surf
      used = send_data(id_entrop_evap, entrop_evap, Time)
    end if
    if (id_entrop_shflx > 0) then
      entrop_shflx = flux_t/t_surf
      used = send_data(id_entrop_shflx, entrop_shflx, Time)
    end if
    if (id_entrop_lwflx > 0) then
      entrop_lwflx = flux_lw/t_surf
      used = send_data(id_entrop_lwflx, entrop_lwflx, Time)
    end if

    call sum_diag_integral_field('evap', flux_q*86400.)

!-----------------------------------------------------------------------

  end subroutine flux_up_to_atmos

!#######################################################################

  !> Initializes the flux exchange and the surface components.
  !>
  !> Stops if `input.nml` still holds the v1 group `simple_surface_nml`, or `diag_table` the
  !> v1 module `simple_surface`, or if the only surface restart in `INPUT/` is the v1
  !> `simple_surface.res.nc`. Reads `flux_exchange_nml`, sets the land mask, initializes
  !> the ocean and then the land state (each from its restart file if it exists), sets the
  !> background albedo and the roughness lengths, allocates `Boundary` and `Exch`, and
  !> registers the diagnostics (`ocean/qflux_implied` only with a prescribed or fixed SST).
  !> Also stops if an initial surface temperature, from a cold start or a restart, is outside
  !> 150-400 K on its component's points. The diagnostics that are sent only with a coupled
  !> surface (`physics_driver_mod`: `surface_is_coupled`) are registered only then.
  subroutine flux_exchange_init(Time, Atm, Ocean, Land, Boundary, Exch)

    type(time_type), intent(in)  :: Time  !! current time
    type(atmos_data_type), intent(in)  :: Atm  !! atmospheric grid and domain
    type(ocean_data_type), intent(inout) :: Ocean  !! the ocean state
    type(land_data_type), intent(inout) :: Land  !! the land state
    type(land_ice_atmos_boundary_type), intent(inout) :: Boundary  !! surface fields for the atmosphere
    type(surface_exchange_type), intent(inout) :: Exch  !! fields for the ocean and land models

    integer :: ierr, io, i, j, id, jd
    real, allocatable, dimension(:, :) :: t_all
    real :: cap_default
    logical :: ocean_rst, land_rst

    call check_v1_input

!------ read namelist ------

    read (input_nml_file, nml=flux_exchange_nml, iostat=io)
    ierr = check_nml_error(io, 'flux_exchange_nml')

    if (mpp_pe() == mpp_root_pe()) then
      write (stdlog(), nml=flux_exchange_nml)
    end if
! Also when sfc_boundary_layer skips the fluxes, so that surface_flux_nml is always checked
    call surface_flux_init

!------ the surface grid ------

    id = size(Atm%t_bot, 1)
    jd = size(Atm%t_bot, 2)
    grid%domain = Atm%domain
    grid%axes = Atm%axes(1:2)
    grid%lon_bnd => Atm%lon_bnd
    grid%lat_bnd => Atm%lat_bnd
    grid%nx = id
    grid%ny = jd
    allocate (grid%area(id, jd))
    do j = 1, jd
      do i = 1, id
        grid%area(i, j) = (Atm%lon_bnd(i + 1) - Atm%lon_bnd(i))*(sin(Atm%lat_bnd(j + 1)) - sin(Atm%lat_bnd(j)))
      end do
    end do

    call diag_integral_field_init('evap', 'f6.3')
    call diag_field_init(Time, grid%axes)
    need_fluxes = surface_is_coupled() .or. any_flux_diagnostic() .or. any_dry_deposition()

!------ components: land mask, then ocean, then land state ------

    ocean_rst = file_exists('INPUT/ocean_model.res.nc') .or. file_exists('INPUT/ocean_model.res.nc.0000')
    land_rst = file_exists('INPUT/land_model.res.nc') .or. file_exists('INPUT/land_model.res.nc.0000')
    if (ocean_rst .neqv. land_rst) &
      call error_mesg('flux_exchange_init', 'only one of INPUT/ocean_model.res.nc and '// &
                      'INPUT/land_model.res.nc exists; a restart needs both', FATAL)
! A v1 surface restart alone would otherwise give a silent cold start of the surface
    if (.not. ocean_rst .and. .not. land_rst .and. &
        (file_exists('INPUT/simple_surface.res.nc') .or. file_exists('INPUT/simple_surface.res.nc.0000'))) &
      call error_mesg('flux_exchange_init', 'INPUT/simple_surface.res.nc is a v1 surface restart, '// &
                      'which is no longer read: convert it to INPUT/ocean_model.res.nc and '// &
                      'INPUT/land_model.res.nc with the recipe in Migration_v2.md, or delete it '// &
                      'to cold-start the surface', FATAL)

    call land_model_read_nml_and_mask(Land, grid)

    allocate (t_all(id, jd))
    call ocean_model_init(Ocean, Time, grid, .not. Land%mask, t_all, cap_default)
    call land_model_init_state(Land, Time, grid, t_all, cap_default)
    deallocate (t_all)
    if (ocean_rst) then
      call read_surface_restarts(Ocean, Land)
      call land_model_read_restart(Land, grid)
    end if
    call check_t_surf_range(Ocean, Land)

! The implied ocean heat flux: only where the SST is not stepped by the slab
    if ((trim(Ocean%model) == 'prescribed' .or. trim(Ocean%model) == 'fixed') .and. surface_is_coupled()) &
      id_qflux_implied = register_diag_field('ocean', 'qflux_implied', grid%axes, Time, &
                                             'implied ocean heat flux convergence: the Q-flux a slab '// &
                                             'ocean needs to reproduce the SST change (0 over land)', 'W/m2')

!------ albedo and roughness ------

    allocate (albedo_bg(id, jd), rough_mom(id, jd), rough_heat(id, jd), rough_moist(id, jd))
    call surface_albedo_init(grid, albedo_bg)
    rough_mom = Ocean%rough_ocean
    rough_heat = Ocean%rough_ocean
    rough_moist = Ocean%rough_ocean
    call land_model_roughness(Land, grid, Ocean%rough_ocean, rough_mom, rough_heat, rough_moist)

!------ module storage ------

    allocate (e_t_n(id, jd), e_q_n(id, jd), f_t_delt_n(id, jd), f_q_delt_n(id, jd), &
              dhdt_surf(id, jd), dedt_surf(id, jd), drdt_surf(id, jd), &
              dhdt_atm(id, jd), dedq_atm(id, jd), flux_t(id, jd), flux_q(id, jd), &
              flux_lw(id, jd), flux_u(id, jd), flux_v(id, jd))

!------ the boundary fields for the atmosphere ------

    allocate (Boundary%t(id, jd), &
              Boundary%albedo(id, jd), &
              Boundary%land_frac(id, jd), &
              Boundary%dtaudu(id, jd), &
              Boundary%dtaudv(id, jd), &
              Boundary%dt_t(id, jd), &
              Boundary%dt_q(id, jd), &
              Boundary%u_flux(id, jd), &
              Boundary%v_flux(id, jd), &
              Boundary%u_star(id, jd), &
              Boundary%b_star(id, jd), &
              Boundary%q_star(id, jd), &
              Boundary%rough_mom(id, jd))

    Boundary%t = 0.0
    Boundary%albedo = albedo_bg
    Boundary%land_frac = 0.0
    Boundary%dtaudu = 0.0
    Boundary%dtaudv = 0.0
    Boundary%dt_t = 0.0
    Boundary%dt_q = 0.0
    Boundary%u_flux = 0.0
    Boundary%v_flux = 0.0
    Boundary%u_star = 0.0
    Boundary%b_star = 0.0
    Boundary%q_star = 0.0
    Boundary%rough_mom = rough_mom

!------ the exchange fields for the ocean and land models ------

    allocate (Exch%flux_lw(id, jd), Exch%flux_sw(id, jd), Exch%lprec(id, jd), Exch%fprec(id, jd), &
              Exch%flux_t(id, jd), Exch%flux_q(id, jd), Exch%dhdt_surf(id, jd), Exch%dedt_surf(id, jd), &
              Exch%drdt_surf(id, jd), Exch%extra_heating(id, jd), Exch%qflux_land(id, jd))
    Exch%flux_lw = 0.0
    Exch%flux_sw = 0.0
    Exch%lprec = 0.0
    Exch%fprec = 0.0
    Exch%flux_t = 0.0
    Exch%flux_q = 0.0
    Exch%dhdt_surf = 0.0
    Exch%dedt_surf = 0.0
    Exch%drdt_surf = 0.0
    Exch%extra_heating = 0.0
! The Q-flux over land: Ocean%qflux holds the Q-flux applied at every point, which is 0
! over land unless qflux_ocean_only = .false. (then the land gets it as in v1)
    Exch%qflux_land = 0.0
    where (Land%mask) Exch%qflux_land = Ocean%qflux
    Exch%dt = 0.0
    Exch%Time = Time

!------ static fields of the step ------

    allocate (u_surf(id, jd), v_surf(id, jd), all_points(id, jd), no_points(id, jd))
    u_surf = 0.0
    v_surf = 0.0
    all_points = .true.
    no_points = .false.

!mj adding local surface heating
    if (do_local_heating) then
      do j = 1, ngauss
        if (hamp(j) .ne. 0. .and. pcenter(j) .lt. 0.) then
          do_surface_heating = .true.
        end if
      end do
    end if
    if (do_surface_heating) then
      allocate (lon2d(id, jd), lat2d(id, jd))
      do j = 1, jd
        do i = 1, id
          lat2d(i, j) = 0.5*(Atm%lat_bnd(j + 1) + Atm%lat_bnd(j))
          lon2d(i, j) = 0.5*(Atm%lon_bnd(i + 1) + Atm%lon_bnd(i))
        end do
      end do
    end if

    do_init = .false.

!-----------------------------------------------------------------------

  end subroutine flux_exchange_init

!#######################################################################

  ! Sets the ocean and land surface temperatures from INPUT/ocean_model.res.nc and
  ! INPUT/land_model.res.nc, which take precedence over the cold-start fields. Each file
  ! holds t_surf on the points of the component that wrote it (0 elsewhere), so the two are
  ! combined (ocean file where it is positive, land file elsewhere) and each component then
  ! takes the combined field on its current mask. Files written by the model cover every
  ! point, so a change of the land mask is accepted: a point that is now land takes the
  ! SST it had, and a point that is now ocean takes the land temperature it had. This also
  ! accepts files converted from a v1 restart (the same SST on all points in both). Stops
  ! if a point of the current mask has no temperature in either file.
  subroutine read_surface_restarts(Ocean, Land)
    type(ocean_data_type), intent(inout) :: Ocean
    type(land_data_type), intent(inout) :: Land

    real, dimension(grid%nx, grid%ny) :: t_ocean, t_land, t_rst

    call read_t_surf('INPUT/ocean_model.res.nc', t_ocean)
    call read_t_surf('INPUT/land_model.res.nc', t_land)
    where (t_ocean > 0.)
      t_rst = t_ocean
    elsewhere
      t_rst = t_land
    end where
    if (any((Ocean%mask .or. Land%mask) .and. t_rst .le. 0.)) &
      call error_mesg('flux_exchange_init', 'INPUT/ocean_model.res.nc and INPUT/land_model.res.nc '// &
                      'have no surface temperature (t_surf <= 0 in both) at some points. Files '// &
                      'written by the model cover every point, even after a change of the land mask '// &
                      '(a point that is now land takes its old SST, one that is now ocean its old '// &
                      'land temperature), so the files were probably edited or converted by hand. '// &
                      'Fill t_surf at those points, or cold-start the surface (remove both files)', FATAL)
    where (Ocean%mask) Ocean%t_surf = t_rst
    where (Land%mask) Land%t_surf = t_rst

  end subroutine read_surface_restarts

  ! Stops if the initial surface temperature, from a cold start or a restart, is outside
  ! 150-400 K at a point of the component's mask. The usual cause is an SST file with fill
  ! values over land, which the land takes through t_init = 'ocean'.
  subroutine check_t_surf_range(Ocean, Land)
    type(ocean_data_type), intent(in) :: Ocean
    type(land_data_type), intent(in) :: Land

    if (any(Ocean%mask .and. .not. (Ocean%t_surf >= 150. .and. Ocean%t_surf <= 400.))) &
      call error_mesg('flux_exchange_init', 'the initial ocean surface temperature is outside '// &
                      '150-400 K at some ocean points. Check the cold-start SST (sst_init, Tm, '// &
                      'deltaT, and the values of sst_file) or, on a restart, INPUT/ocean_model.res.nc', FATAL)
    if (any(Land%mask .and. .not. (Land%t_surf >= 150. .and. Land%t_surf <= 400.))) &
      call error_mesg('flux_exchange_init', 'the initial land surface temperature is outside '// &
                      '150-400 K at some land points. A likely cause: an SST file with fill '// &
                      "values over land, which the land takes on a cold start with t_init = 'ocean' "// &
                      "(sst_init = 'file' or ocean_model = 'prescribed'). Fill the land points of the "// &
                      "SST file, or set t_init = 'file' with a t_surf_file. Otherwise check t_surf_file "// &
                      'or, on a restart, INPUT/land_model.res.nc', FATAL)

  end subroutine check_t_surf_range

  ! Reads t_surf from a surface restart file.
  subroutine read_t_surf(path, t_surf)
    character(len=*), intent(in) :: path
    real, intent(out) :: t_surf(:, :)
    type(restart_file_type) :: rst

    if (.not. open_restart_read(rst, path, grid%domain)) &
      call error_mesg('flux_exchange_init', 'cannot open '//path, FATAL)
    if (.not. restart_field_exists(rst, 't_surf')) &
      call error_mesg('flux_exchange_init', path//' has no field t_surf', FATAL)
    call read_restart_field(rst, 't_surf', t_surf)
    call close_restart(rst)

  end subroutine read_t_surf

!#######################################################################

  ! Stops if the input files are still in the v1 format of the surface: unknown namelist
  ! groups and diagnostic modules are otherwise ignored without a message.
  subroutine check_v1_input

    integer :: unit, io
    character(len=1024) :: line

    if (nml_has_group('simple_surface_nml')) &
      call error_mesg('flux_exchange_init', 'input.nml contains the v1 namelist simple_surface_nml, '// &
                      'which is no longer read; convert it to the v2 surface namelists '// &
                      '(see Migration_v2.md)', FATAL)

    if (mpp_pe() == mpp_root_pe() .and. file_exists('diag_table')) then
      open (newunit=unit, file='diag_table', form='formatted', action='read', status='old')
      do
        read (unit, '(a)', iostat=io) line
        if (io /= 0) exit
        line = adjustl(line)
        if (line(1:1) == '#') cycle
        if (index(line, '"simple_surface"') > 0 .or. index(line, "'simple_surface'") > 0) then
          close (unit)
          call error_mesg('flux_exchange_init', 'diag_table requests fields of the v1 module '// &
                          '"simple_surface", which no longer exists; use the modules "flux", "ocean" '// &
                          'and "land" (see Migration_v2.md)', FATAL)
        end if
      end do
      close (unit)
    end if

  end subroutine check_v1_input

!#######################################################################

  ! Whether a diagnostic of the module flux is requested (called after diag_field_init)
  logical function any_flux_diagnostic()

    any_flux_diagnostic = any([id_drag_moist, id_drag_heat, id_drag_mom, &
                               id_rough_moist, id_rough_heat, id_rough_mom, &
                               id_u_star, id_b_star, id_u_flux, id_v_flux, &
                               id_t_atm, id_u_atm, id_v_atm, id_wind, &
                               id_t_ref, id_rh_ref, id_u_ref, id_v_ref, &
                               id_del_h, id_del_m, id_del_q, id_albedo, &
                               id_t_surf, id_t_flux, id_q_flux, id_r_flux, id_heat, &
                               id_entrop_evap, id_entrop_shflx, id_entrop_lwflx] > 0)

  end function any_flux_diagnostic

  ! Whether a prognostic tracer has a dry_deposition method, which uses u_star
  ! (atmos_tracer_driver)
  logical function any_dry_deposition()

    integer :: n, ntprog
    character(len=128) :: name, control

    any_dry_deposition = .false.
    call get_number_tracers(MODEL_ATMOS, num_prog=ntprog)
    do n = 1, ntprog
      if (query_method('dry_deposition', MODEL_ATMOS, n, name, control)) any_dry_deposition = .true.
    end do

  end function any_dry_deposition

!#######################################################################

  subroutine diag_field_init(Time, atmos_axes)

    type(time_type), intent(in) :: Time
    integer, intent(in) :: atmos_axes(2)

    integer :: iref
    character(len=6) :: label_zm, label_zh
    real, dimension(2) :: trange = (/100., 400./), &
                          vrange = (/-400., 400./)
!-----------------------------------------------------------------------
!  initializes diagnostic fields that may be output from this module
!  (the id numbers may be referenced anywhere in this module)
!-----------------------------------------------------------------------

!------ labels for diagnostics -------
!  (z_ref_mom, z_ref_heat are namelist variables)

    iref = int(z_ref_mom + 0.5)
    if (real(iref) == z_ref_mom) then
      write (label_zm, 105) iref
      if (iref < 10) write (label_zm, 100) iref
    else
      write (label_zm, 110) z_ref_mom
    end if

    iref = int(z_ref_heat + 0.5)
    if (real(iref) == z_ref_heat) then
      write (label_zh, 105) iref
      if (iref < 10) write (label_zh, 100) iref
    else
      write (label_zh, 110) z_ref_heat
    end if

100 format(i1, ' m', 3x)
105 format(i2, ' m', 2x)
110 format(f4.1, ' m')

    id_wind = &
      register_diag_field(mod_name, 'wind', atmos_axes, Time, &
                          'wind speed for flux calculations', 'm/s', &
                          range=(/0., vrange(2)/))

    id_drag_moist = &
      register_diag_field(mod_name, 'drag_moist', atmos_axes, Time, &
                          'drag coeff for moisture', '1')

    id_drag_heat = &
      register_diag_field(mod_name, 'drag_heat', atmos_axes, Time, &
                          'drag coeff for heat', '1')

    id_drag_mom = &
      register_diag_field(mod_name, 'drag_mom', atmos_axes, Time, &
                          'drag coeff for momentum', '1')

    id_rough_moist = &
      register_diag_field(mod_name, 'rough_moist', atmos_axes, Time, &
                          'surface roughness for moisture', 'm')

    id_rough_heat = &
      register_diag_field(mod_name, 'rough_heat', atmos_axes, Time, &
                          'surface roughness for heat', 'm')

    id_rough_mom = &
      register_diag_field(mod_name, 'rough_mom', atmos_axes, Time, &
                          'surface roughness for momentum', 'm')

    id_u_star = &
      register_diag_field(mod_name, 'u_star', atmos_axes, Time, &
                          'friction velocity', 'm/s')

    id_b_star = &
      register_diag_field(mod_name, 'b_star', atmos_axes, Time, &
                          'buoyancy scale', 'm/s2')

    id_u_flux = &
      register_diag_field(mod_name, 'tau_x', atmos_axes, Time, &
                          'zonal surface stress on the atmosphere (positive eastward)', 'N/m2')

    id_v_flux = &
      register_diag_field(mod_name, 'tau_y', atmos_axes, Time, &
                          'meridional surface stress on the atmosphere (positive northward)', 'N/m2')

    if (surface_is_coupled()) then
      id_t_surf = &
        register_diag_field(mod_name, 't_surf', atmos_axes, Time, &
                            'surface temperature', 'K', &
                            range=trange)

      id_t_flux = &
        register_diag_field(mod_name, 'shflx', atmos_axes, Time, &
                            'sensible heat flux', 'W/m2')

      id_q_flux = &
        register_diag_field(mod_name, 'evap', atmos_axes, Time, &
                            'evaporation rate', 'kg/m2/s')

      id_r_flux = &
        register_diag_field(mod_name, 'lwflx', atmos_axes, Time, &
                            'net (down-up) longwave flux', 'W/m2')
    end if

    id_t_atm = &
      register_diag_field(mod_name, 't_atm', atmos_axes, Time, &
                          'temperature at btm level', 'K', &
                          range=trange)

    id_u_atm = &
      register_diag_field(mod_name, 'u_atm', atmos_axes, Time, &
                          'u wind component at btm level', 'm/s', &
                          range=vrange)

    id_v_atm = &
      register_diag_field(mod_name, 'v_atm', atmos_axes, Time, &
                          'v wind component at btm level', 'm/s', &
                          range=vrange)

    id_t_ref = &
      register_diag_field(mod_name, 't_ref', atmos_axes, Time, &
                          'air temperature at '//trim(label_zh), 'K', &
                          range=trange)

    id_rh_ref = &
      register_diag_field(mod_name, 'rh_ref', atmos_axes, Time, &
                          'relative humidity at '//trim(label_zh)//' (100 q/q_sat)', 'percent')

    id_u_ref = &
      register_diag_field(mod_name, 'u_ref', atmos_axes, Time, &
                          'zonal wind at '//trim(label_zm), 'm/s', &
                          range=vrange)

    id_v_ref = &
      register_diag_field(mod_name, 'v_ref', atmos_axes, Time, &
                          'meridional wind at '//trim(label_zm), 'm/s', &
                          range=vrange)

    id_del_h = &
      register_diag_field(mod_name, 'del_h', atmos_axes, Time, &
                          'Monin-Obukhov profile factor (T('//trim(label_zh)//')-T_surf)/(T_atm-T_surf)', '1')
    id_del_m = &
      register_diag_field(mod_name, 'del_m', atmos_axes, Time, &
                          'Monin-Obukhov profile factor u('//trim(label_zm)//')/u_atm', '1')
    id_del_q = &
      register_diag_field(mod_name, 'del_q', atmos_axes, Time, &
                          'Monin-Obukhov profile factor (q('//trim(label_zh)//')-q_surf)/(q_atm-q_surf)', '1')
    id_albedo = &
      register_diag_field(mod_name, 'albedo', atmos_axes, Time, &
                          'surface albedo', '1')
    if (surface_is_coupled()) then
      id_heat = & !mj
        register_diag_field(mod_name, 'heat_capacity', atmos_axes, Time, &
                            'mixed layer heat capacity', 'J/m2/K')
      id_entrop_evap = &
        register_diag_field(mod_name, 'entrop_evap', atmos_axes, Time, &
                            'entropy source from evap', 'kg/m2/s/K')

      id_entrop_shflx = &
        register_diag_field(mod_name, 'entrop_shflx', atmos_axes, Time, &
                            'entropy source from SH flux', 'W/m2/K')

      id_entrop_lwflx = &
        register_diag_field(mod_name, 'entrop_lwflx', atmos_axes, Time, &
                            'entropy source from LW flux', 'W/m2/K')
    end if

!-----------------------------------------------------------------------

  end subroutine diag_field_init

!########################################################################

  !> Terminates the surface components, which write their restart files
  !> (`RESTART/ocean_model.res.nc` and `RESTART/land_model.res.nc`).
  subroutine flux_exchange_end(Ocean, Land)

    type(ocean_data_type), intent(in) :: Ocean  !! the ocean state
    type(land_data_type), intent(in) :: Land  !! the land state

    call ocean_model_end(Ocean, grid)
    call land_model_end(Land, grid)
! The boundaries belong to the atmosphere (Atm%lon_bnd, Atm%lat_bnd)
    nullify (grid%lon_bnd, grid%lat_bnd)

  end subroutine flux_exchange_end

!#######################################################################

end module flux_exchange_mod
