
!> Tendencies due to vertical diffusion, from an implicit tridiagonal solver.
!>
!> Momentum, temperature (diffused as `t + grav*z/cp_air`), specific humidity and the
!> other prognostic tracers are diffused implicitly in time. The solution is split into a
!> downward sweep (`gcm_vert_diff_down`) and an upward sweep (`gcm_vert_diff_up`), so
!> that the surface model can compute the surface fluxes of heat and moisture implicitly
!> in between, from the data stored in a `surf_diff_type`. The air density in the
!> diffusion uses the temperature, or the virtual temperature if
!> `use_virtual_temp_vert_diff`. With `do_conserve_energy`, the kinetic energy removed by
!> the momentum diffusion heats the air. Tracers whose `diff_vert` method in the field
!> table is `'none'` are not diffused.
!>
!> The switches are set by `vert_diff_driver_mod` from `vert_diff_driver_nml`. A
!> technical description is in `vert_diff.tech.ps` in the source directory.
module vert_diff_mod

  use constants_mod, only: GRAV, RDGAS, RVGAS, CP_AIR

  use fms_mod, only: error_mesg, FATAL, uppercase, &
                     write_version_number, stdlog, &
                     mpp_pe, mpp_root_pe

  use field_manager_mod, only: MODEL_ATMOS
  use tracer_manager_mod, only: query_method, get_tracer_index, NO_TRACER

  implicit none
  private

! public interfaces
!=======================================================================
  public :: vert_diff_init, &
            vert_diff_end, &
            gcm_vert_diff, &
            gcm_vert_diff_down, &
            gcm_vert_diff_up, &
            vert_diff, &
            surf_diff_type

!=======================================================================

! form of interfaces
!=======================================================================

  !> Data at the lowest model level passed between the downward sweep of the vertical
  !> diffusion, the surface model and the upward sweep.
  type surf_diff_type

    real, pointer, dimension(:, :) :: dtmass => null(), &
                                      dflux_t => null(), &
                                      dflux_q => null(), &
                                      delta_t => null(), &
                                      delta_q => null(), &
                                      delta_u => null(), &
                                      delta_v => null()
    !! `dtmass`: `grav*delt` over the pressure thickness of the lowest level [m2 s/kg];
    !! `dflux_t`, `dflux_q`: derivatives of the diffusive fluxes into the lowest level with
    !! respect to its temperature and specific humidity [kg/m2/s]; `delta_t`, `delta_q`:
    !! changes of the lowest-level temperature [K] and specific humidity [kg/kg] over the
    !! step (without the surface flux; the surface model replaces them with the values that
    !! include it); `delta_u`, `delta_v`: changes of the lowest-level u and v over the step [m/s]

  end type surf_diff_type

  real, allocatable, dimension(:, :, :) :: e_global, f_t_global, f_q_global

  logical :: do_conserve_energy = .true.
  logical :: use_virtual_temp_vert_diff
  integer :: sphum, mix_rat

!--------------------- version number ---------------------------------

  character(len=128) :: version = '$Id: vert_diff.f90,v 12.0 2005/04/14 15:50:38 fms Exp $'
  character(len=128) :: tagname = '$Name: lima $'
  logical            :: module_is_initialized = .false.

  real, parameter :: d608 = (RVGAS - RDGAS)/RDGAS

contains

!#######################################################################

  !> Initializes the module: finds the specific-humidity tracer, sets the switches and
  !> allocates the arrays kept between the downward and upward sweeps and those of
  !> `Tri_surf`.
  subroutine vert_diff_init(Tri_surf, idim, jdim, kdim, &
                            do_conserve_energy_in, &
                            use_virtual_temp_vert_diff_in)

    type(surf_diff_type), intent(inout) :: Tri_surf  !! surface diffusion data (allocated here)
    integer, intent(in)    :: idim, jdim, kdim  !! sizes of the local grid (longitude, latitude, levels)
    logical, intent(in)    :: do_conserve_energy_in
    !! heat the air by the kinetic energy dissipated by the momentum diffusion
    logical, optional, intent(in)    :: use_virtual_temp_vert_diff_in
    !! use the virtual temperature for the air density (default `.false.`)

    call write_version_number(version, tagname)

! get the tracer number for specific humidity
    sphum = get_tracer_index(MODEL_ATMOS, 'sphum')
    mix_rat = get_tracer_index(MODEL_ATMOS, 'mix_rat')
    if (sphum /= NO_TRACER .and. mix_rat /= NO_TRACER) then
      call error_mesg('gcm_vert_diff_init', 'sphum and mix_rat cannot both'// &
                      'be present in the field_table at the same time', FATAL)
    end if

    if (mpp_pe() == mpp_root_pe()) then
      write (stdlog(), '(a,i4)') 'Tracer number for specific humidity =', sphum
      write (stdlog(), '(a,i4)') 'Tracer number for mixing ratio      =', mix_rat
    end if

    if (present(use_virtual_temp_vert_diff_in)) then
      use_virtual_temp_vert_diff = use_virtual_temp_vert_diff_in
    else
      use_virtual_temp_vert_diff = .false.
    end if

    if (.not. module_is_initialized) then

      if (allocated(e_global)) deallocate (e_global)
      if (allocated(f_t_global)) deallocate (f_t_global)
      if (allocated(f_q_global)) deallocate (f_q_global)

      allocate (e_global(idim, jdim, kdim - 1)); e_global = 0.0
      allocate (f_t_global(idim, jdim, kdim - 1)); f_t_global = 0.0
      allocate (f_q_global(idim, jdim, kdim - 1)); f_q_global = 0.0

      module_is_initialized = .true.

    end if

    call alloc_surf_diff_type(Tri_surf, idim, jdim)

    do_conserve_energy = do_conserve_energy_in

  end subroutine vert_diff_init

!#######################################################################

  subroutine alloc_surf_diff_type(Tri_surf, idim, jdim)

    type(surf_diff_type), intent(inout) :: Tri_surf
    integer, intent(in)    :: idim, jdim

    allocate (Tri_surf%dtmass(idim, jdim)); Tri_surf%dtmass = 0.0
    allocate (Tri_surf%dflux_t(idim, jdim)); Tri_surf%dflux_t = 0.0
    allocate (Tri_surf%dflux_q(idim, jdim)); Tri_surf%dflux_q = 0.0
    allocate (Tri_surf%delta_t(idim, jdim)); Tri_surf%delta_t = 0.0
    allocate (Tri_surf%delta_q(idim, jdim)); Tri_surf%delta_q = 0.0
    allocate (Tri_surf%delta_u(idim, jdim)); Tri_surf%delta_u = 0.0
    allocate (Tri_surf%delta_v(idim, jdim)); Tri_surf%delta_v = 0.0

  end subroutine alloc_surf_diff_type

!#######################################################################

  subroutine dealloc_surf_diff_type(Tri_surf)

    type(surf_diff_type), intent(inout) :: Tri_surf

    deallocate (Tri_surf%dtmass)
    deallocate (Tri_surf%dflux_t)
    deallocate (Tri_surf%dflux_q)
    deallocate (Tri_surf%delta_t)
    deallocate (Tri_surf%delta_q)
    deallocate (Tri_surf%delta_u)
    deallocate (Tri_surf%delta_v)

  end subroutine dealloc_surf_diff_type

!#######################################################################

  !> Deallocates the arrays of the module.
  subroutine vert_diff_end

    if (module_is_initialized) then

      if (allocated(e_global)) deallocate (e_global)
      if (allocated(f_t_global)) deallocate (f_t_global)
      if (allocated(f_q_global)) deallocate (f_q_global)

    end if
    module_is_initialized = .false.

  end subroutine vert_diff_end

!#######################################################################

  !> Diffuses momentum and the tracers and does the downward sweep for temperature and
  !> specific humidity.
  !>
  !> The tendencies of u, v and the tracers are complete on return; those of temperature
  !> and specific humidity are completed by `gcm_vert_diff_up` after the surface model has
  !> updated `Tri_surf`.
  subroutine gcm_vert_diff_down(is, js, delt, &
                                u, v, t, q, tr, &
                                diff_m, diff_t, p_half, p_full, &
                                z_full, tau_u, tau_v, &
                                dtau_du, dtau_dv, &
                                dt_u, dt_v, dt_t, dt_q, dt_tr, &
                                dissipative_heat, Tri_surf, &
                                kbot)

    integer, intent(in)                        :: is, js  !! starting i,j indices of the physics window
    real, intent(in)                        :: delt  !! time step [s]
    real, intent(in), dimension(:, :, :)   :: u, v, t, q, &
                                              diff_m, diff_t, &
                                              p_half, p_full, &
                                              z_full
    !! `u`, `v`: zonal and meridional wind [m/s]; `t`: temperature [K]; `q`: specific humidity
    !! [kg/kg]; `diff_m`, `diff_t`: diffusivities for momentum and for heat and tracers at
    !! half levels [m2/s]; `p_half`, `p_full`: pressure at half and full levels [Pa];
    !! `z_full`: height of full levels [m]
    real, intent(in), dimension(:, :, :, :) :: tr  !! prognostic tracers
    real, intent(in), dimension(:, :)     :: dtau_du, dtau_dv
    !! derivatives of the surface stresses with respect to the lowest-level u and v [kg/m2/s]
    real, intent(inout), dimension(:, :)     :: tau_u, tau_v
    !! zonal and meridional surface stress, updated for the implicit change of the
    !! lowest-level wind [N/m2]
    real, intent(inout), dimension(:, :, :)   :: dt_u, dt_v, dt_t
    !! tendencies of u and v [m/s2], including the diffusion on return, and of temperature
    !! [K/s], with the dissipative heating added on return
    real, intent(in), dimension(:, :, :)   :: dt_q  !! specific humidity tendency [kg/kg/s]
    real, intent(inout), dimension(:, :, :, :) :: dt_tr  !! tracer tendencies, including the diffusion on return
    real, intent(out), dimension(:, :, :)   :: dissipative_heat
    !! heating by the kinetic energy dissipated by the momentum diffusion (0 unless
    !! `do_conserve_energy`) [K/s]
    type(surf_diff_type), intent(inout)        :: Tri_surf  !! surface diffusion data for the surface model

    integer, intent(in), dimension(:, :), optional :: kbot  !! index of the lowest model level

    real, dimension(size(u, 1), size(u, 2), size(u, 3)) :: tt, mu, nu

    real, dimension(size(u, 1), size(u, 2)) :: f_t_delt_n1, f_q_delt_n1, &
                                               mu_delt_n, nu_n, e_n1, delta_t_n, delta_q_n, &
                                               delta_u_n, delta_v_n

    real    :: gcp
    integer :: i, j, kb, ie, je

!-----------------------------------------------------------------------

    if (.not. module_is_initialized) call error_mesg('gcm_vert_diff_down in vert_diff_mod', &
                                                     'the initialization routine gcm_vert_diff_init has not been called', &
                                                     FATAL)

    ie = is + size(t, 1) - 1
    je = js + size(t, 2) - 1

    gcp = GRAV/CP_AIR
    tt = t + z_full*gcp   ! the vertical gradient of tt determines the
    ! diffusive flux of temperature

    call compute_mu(p_half, mu)
    call compute_nu(diff_m, p_half, p_full, z_full, t, q, nu)

!  diffuse u-momentum and v_momentum

    call uv_vert_diff(delt, mu, nu, u, v, dtau_du, dtau_dv, tau_u, tau_v, &
                      dt_u, dt_v, dt_t, delta_u_n, delta_v_n, &
                      dissipative_heat, kbot)

!  recompute nu for a different diffusivity
    call compute_nu(diff_t, p_half, p_full, z_full, t, q, nu)

!  diffuse tracers
    call tr_vert_diff(delt, mu, nu, tr, dt_tr, kbot)

!  downward sweep of tridiagonal solver for temperature and specific humidity
    call vert_diff_down_2 &
      (delt, mu, nu, tt, q, dt_t, dt_q, &
       e_global(is:ie, js:je, :), &
       f_t_global(is:ie, js:je, :), &
       f_q_global(is:ie, js:je, :), &
       mu_delt_n, nu_n, e_n1, f_t_delt_n1, f_q_delt_n1, &
       delta_t_n, delta_q_n, kbot)

! store information needed by flux_exchange module

    Tri_surf%delta_t(is:ie, js:je) = delta_t_n + mu_delt_n*nu_n*f_t_delt_n1
    Tri_surf%delta_q(is:ie, js:je) = delta_q_n + mu_delt_n*nu_n*f_q_delt_n1
    Tri_surf%dflux_t(is:ie, js:je) = -nu_n*(1.0 - e_n1)
    Tri_surf%dflux_q(is:ie, js:je) = -nu_n*(1.0 - e_n1)
    Tri_surf%dtmass(is:ie, js:je) = mu_delt_n
    Tri_surf%delta_u(is:ie, js:je) = delta_u_n
    Tri_surf%delta_v(is:ie, js:je) = delta_v_n

!-----------------------------------------------------------------------

  end subroutine gcm_vert_diff_down

!#######################################################################

  !> Does the upward sweep for temperature and specific humidity, from the lowest-level
  !> changes `Tri_surf%delta_t` and `Tri_surf%delta_q` computed by the surface model.
  subroutine gcm_vert_diff_up(is, js, delt, Tri_surf, dt_t, dt_q, kbot)

    integer, intent(in)                      :: is, js  !! starting i,j indices of the physics window
    real, intent(in)                      :: delt  !! time step [s]
    type(surf_diff_type), intent(in)         :: Tri_surf  !! surface diffusion data from the surface model
    real, intent(out), dimension(:, :, :) :: dt_t, dt_q
    !! tendencies of temperature [K/s] and specific humidity [kg/kg/s], including the
    !! diffusion
    integer, intent(in), dimension(:, :), optional :: kbot  !! index of the lowest model level

    integer :: ie, je

    ie = is + size(dt_t, 1) - 1
    je = js + size(dt_t, 2) - 1

    call vert_diff_up(delt, &
                      e_global(is:ie, js:je, :), &
                      f_t_global(is:ie, js:je, :), &
                      Tri_surf%delta_t(is:ie, js:je), &
                      dt_t, kbot)

    call vert_diff_up(delt, &
                      e_global(is:ie, js:je, :), &
                      f_q_global(is:ie, js:je, :), &
                      Tri_surf%delta_q(is:ie, js:je), &
                      dt_q, kbot)

  end subroutine gcm_vert_diff_up

!#######################################################################

  !> Diffuses momentum, temperature, specific humidity and the tracers in one step, for a
  !> model in which the surface fluxes do not depend implicitly on the surface
  !> temperature.
  subroutine gcm_vert_diff(delt, u, v, t, q, tr, &
                           diff_m, diff_t, p_half, p_full, z_full, &
                           dtau_du, dtau_dv, dsens_datmos, devap_datmos, &
                           sens, evap, tau_u, tau_v, &
                           dt_u, dt_v, dt_t, dt_q, dt_tr, &
                           dissipative_heat, kbot)

    real, intent(in)                          :: delt  !! time step [s]
    real, intent(in), dimension(:, :, :)     :: u, v, t, q, p_half, p_full, &
                                                z_full, diff_m, diff_t
    !! `u`, `v`: zonal and meridional wind [m/s]; `t`: temperature [K]; `q`: specific humidity
    !! [kg/kg]; `p_half`, `p_full`: pressure at half and full levels [Pa]; `z_full`: height of
    !! full levels [m]; `diff_m`, `diff_t`: diffusivities for momentum and for heat and
    !! tracers at half levels [m2/s]
    real, intent(in), dimension(:, :, :, :)   :: tr  !! prognostic tracers
    real, intent(in), dimension(:, :)       :: dtau_du, dtau_dv, dsens_datmos, &
                                               devap_datmos
    !! derivatives of the surface stresses, sensible heat flux and evaporation with respect
    !! to the lowest-level u, v, temperature and specific humidity
    real, intent(inout), dimension(:, :)       :: tau_u, tau_v, sens, evap
    !! surface stresses [N/m2], sensible heat flux [W/m2] and evaporation [kg/m2/s], updated
    !! for the implicit change of the lowest-level values
    real, intent(inout), dimension(:, :, :)     :: dt_u, dt_v, dt_t, dt_q
    !! tendencies of u, v [m/s2], temperature [K/s] and specific humidity [kg/kg/s],
    !! including the diffusion on return
    real, intent(inout), dimension(:, :, :, :)   :: dt_tr  !! tracer tendencies, including the diffusion on return
    real, intent(out), dimension(:, :, :)     :: dissipative_heat
    !! heating by the kinetic energy dissipated by the momentum diffusion (0 unless
    !! `do_conserve_energy`) [K/s]

    integer, intent(in), dimension(:, :), optional :: kbot  !! index of the lowest model level

    real, dimension(size(u, 1), size(u, 2), size(u, 3)) :: mu, nu
    real, dimension(size(u, 1), size(u, 2))           :: delta_u_n, delta_v_n

!-----------------------------------------------------------------------

    call compute_mu(p_half, mu)

    call compute_nu(diff_m, p_half, p_full, z_full, t, q, nu)

    call uv_vert_diff(delt, mu, nu, u, v, dtau_du, dtau_dv, tau_u, tau_v, &
                      dt_u, dt_v, dt_t, delta_u_n, delta_v_n, &
                      dissipative_heat, kbot)

    call compute_nu(diff_t, p_half, p_full, z_full, t, q, nu)

    call tq_vert_diff(delt, mu, nu, t, q, z_full, &
                      dsens_datmos, devap_datmos, &
                      sens, evap, dt_t, dt_q, kbot)

    call tr_vert_diff(delt, mu, nu, tr, dt_tr, kbot)

  end subroutine gcm_vert_diff

!#######################################################################

  !> Diffuses a single field in one step, with a surface flux.
  subroutine vert_diff(delt, xi, t, q, diff, p_half, p_full, z_full, &
                       flux, dflux_datmos, factor, dt_xi, kbot)

    real, intent(in)                          :: delt  !! time step [s]
    real, intent(in), dimension(:, :, :)     :: xi, t, q, diff, p_half, p_full, z_full
    !! `xi`: the field to diffuse; `t`: temperature [K]; `q`: specific humidity [kg/kg];
    !! `diff`: diffusivity at half levels [m2/s]; `p_half`, `p_full`: pressure at half and
    !! full levels [Pa]; `z_full`: height of full levels [m]
    real, intent(inout), dimension(:, :)       :: flux
    !! surface flux, updated for the implicit change of the lowest-level value
    real, intent(in), dimension(:, :)       :: dflux_datmos
    !! derivative of `flux` with respect to the lowest-level value
    real, intent(in)                          :: factor  !! `flux/factor` is the flux of `xi`
    real, intent(inout), dimension(:, :, :)     :: dt_xi  !! tendency of `xi`, including the diffusion on return

    integer, intent(in), dimension(:, :), optional :: kbot  !! index of the lowest model level

    real, dimension(size(xi, 1), size(xi, 2), size(xi, 3)) :: mu, nu
    real, dimension(size(xi, 1), size(xi, 2), size(xi, 3) - 1) :: e, f

    real, dimension(size(xi, 1), size(xi, 2))  :: mu_delt_n, nu_n, e_n1, &
                                                  f_delt_n1, delta_xi_n

!-----------------------------------------------------------------------

    call compute_mu(p_half, mu)

    call compute_nu(diff, p_half, p_full, z_full, t, q, nu)

    call vert_diff_down &
      (delt, mu, nu, xi, dt_xi, e, f, mu_delt_n, nu_n, e_n1, &
       f_delt_n1, delta_xi_n, kbot)

    call diff_surface(mu_delt_n, nu_n, e_n1, f_delt_n1, &
                      dflux_datmos, flux, factor, delta_xi_n)

    call vert_diff_up(delt, e, f, delta_xi_n, dt_xi, kbot)

  end subroutine vert_diff

!#######################################################################

  subroutine uv_vert_diff(delt, mu, nu, u, v, &
                          dtau_du, dtau_dv, tau_u, tau_v, dt_u, dt_v, dt_t, &
                          delta_u_n, delta_v_n, dissipative_heat, kbot)

    real, intent(in)                        :: delt
    real, intent(in), dimension(:, :, :)   :: u, v, mu, nu
    real, intent(in), dimension(:, :)     :: dtau_du, dtau_dv
    real, intent(inout), dimension(:, :)     :: tau_u, tau_v
    real, intent(inout), dimension(:, :, :)   :: dt_u, dt_v, dt_t
    real, intent(out), dimension(:, :, :)   :: dissipative_heat
    real, intent(out), dimension(:, :)     :: delta_u_n, delta_v_n

! Note (IH)
!   delta_u_n = dt_u/delt at lowest model level, and similarly
!   for delta_v_n  -- it is convenient to output them separately

    integer, intent(in), dimension(:, :), optional :: kbot

    real, dimension(size(u, 1), size(u, 2)) :: mu_delt_n, nu_n, e_n1, &
                                               f_u_delt_n1, f_v_delt_n1

    real, dimension(size(u, 1), size(u, 2), size(u, 3)) :: dt_u_temp, dt_v_temp

    real, dimension(size(u, 1), size(u, 2), size(u, 3) - 1) :: e, f_u, f_v
    integer :: i, j, kb

    real    :: half_delt, cp_inv

!-----------------------------------------------------------------------

    half_delt = 0.5*delt
    cp_inv = 1.0/CP_AIR

    if (do_conserve_energy) then
      dt_u_temp = dt_u
      dt_v_temp = dt_v
    end if

    call vert_diff_down_2 &
      (delt, mu, nu, u, v, dt_u, dt_v, e, f_u, f_v, &
       mu_delt_n, nu_n, e_n1, f_u_delt_n1, f_v_delt_n1, &
       delta_u_n, delta_v_n, kbot)

    call diff_surface(mu_delt_n, nu_n, e_n1, f_u_delt_n1, &
                      dtau_du, tau_u, 1.0, delta_u_n)
    call diff_surface(mu_delt_n, nu_n, e_n1, f_v_delt_n1, &
                      dtau_dv, tau_v, 1.0, delta_v_n)

    call vert_diff_up(delt, e, f_u, delta_u_n, dt_u, kbot)
    call vert_diff_up(delt, e, f_v, delta_v_n, dt_v, kbot)

    if (do_conserve_energy) then
      dt_u_temp = dt_u - dt_u_temp
      dt_v_temp = dt_v - dt_v_temp
      dissipative_heat = -cp_inv*((u + half_delt*dt_u_temp)*dt_u_temp &
                                  + (v + half_delt*dt_v_temp)*dt_v_temp)
      dt_t = dt_t + dissipative_heat
    else
      dissipative_heat = 0.0
    end if

!-----------------------------------------------------------------------

  end subroutine uv_vert_diff

!#######################################################################

  subroutine tq_vert_diff(delt, mu, nu, t, q, z_full, &
                          dsens_datmos, devap_datmos, sens, evap, &
                          dt_t, dt_q, kbot)

    real, intent(in)                        :: delt
    real, intent(in), dimension(:, :, :)   :: t, q, z_full, mu, nu
    real, intent(in), dimension(:, :)     :: dsens_datmos, devap_datmos
    real, intent(inout), dimension(:, :)     :: sens, evap
    real, intent(inout), dimension(:, :, :)   :: dt_t, dt_q

    integer, intent(in), dimension(:, :), optional :: kbot

    real, dimension(size(t, 1), size(t, 2)) :: mu_delt_n, nu_n, &
                                               e_n1, f_t_delt_n1, f_q_delt_n1, &
                                               delta_t_n, delta_q_n, &
                                               flux1, dflux1_dsurf

    real, dimension(size(t, 1), size(t, 2), size(t, 3) - 1) :: e, f_t, f_q
    real, dimension(size(t, 1), size(t, 2), size(t, 3)) :: tt

    integer :: i, j, kb
    real    :: gcp
!-----------------------------------------------------------------------

    gcp = GRAV/CP_AIR
    tt = t + z_full*gcp

    call vert_diff_down_2 &
      (delt, mu, nu, tt, q, dt_t, dt_q, e, f_t, f_q, &
       mu_delt_n, nu_n, e_n1, f_t_delt_n1, f_q_delt_n1, &
       delta_t_n, delta_q_n, kbot)

    call diff_surface(mu_delt_n, nu_n, e_n1, f_t_delt_n1, &
                      dsens_datmos, sens, CP_AIR, delta_t_n)

    call diff_surface(mu_delt_n, nu_n, e_n1, f_q_delt_n1, &
                      devap_datmos, evap, 1.0, delta_q_n)

    call vert_diff_up(delt, e, f_t, delta_t_n, dt_t, kbot)
    call vert_diff_up(delt, e, f_q, delta_q_n, dt_q, kbot)

!-----------------------------------------------------------------------

  end subroutine tq_vert_diff

!#######################################################################

  subroutine tr_vert_diff(delt, mu, nu, tr, dt_tr, kbot)

    real, intent(in)                        :: delt
    real, intent(in), dimension(:, :, :)   :: mu, nu
    real, intent(in), dimension(:, :, :, :) :: tr
    real, intent(inout), dimension(:, :, :, :) :: dt_tr

    integer, intent(in), dimension(:, :), optional :: kbot

    real, dimension(size(tr, 1), size(tr, 2)) :: mu_delt_n, nu_n, e_n1

    real, dimension(size(tr, 1), size(tr, 2)) :: f_delt_n1, delta_tr_n
    real, dimension(size(tr, 1), size(tr, 2)) :: dflux_dtr, flux
    real, dimension(size(tr, 1), size(tr, 2), size(tr, 3) - 1) :: ftr
    real, dimension(size(tr, 1), size(tr, 2), size(tr, 3) - 1) :: etr
    real, dimension(size(tr, 1), size(tr, 2), size(tr, 3)) :: a, b, c, g
    integer :: i, j, k, kb, n, ntr, nlev
    character(len=128) :: scheme
!-----------------------------------------------------------------------

    ntr = size(tr, 4) ! number of prognostic tracers

    dflux_dtr = 0.0
    call compute_e(delt, mu, nu, etr, a, b, c, g)
    if (present(kbot)) then
      do j = 1, size(tr, 2)
      do i = 1, size(tr, 1)
        kb = kbot(i, j)
        mu_delt_n(i, j) = mu(i, j, kb)*delt
        nu_n(i, j) = nu(i, j, kb)
        e_n1(i, j) = etr(i, j, kb - 1)
      end do
      end do
    else
      nlev = size(mu, 3)
      mu_delt_n(:, :) = mu(:, :, nlev)*delt
      nu_n(:, :) = nu(:, :, nlev)
      e_n1(:, :) = etr(:, :, nlev - 1)
    end if

    do n = 1, ntr
      if (n == sphum .or. n == mix_rat) cycle
      if (query_method('diff_vert', MODEL_ATMOS, n, scheme)) then
        if (uppercase(trim(scheme)) == 'NONE') cycle
      end if
      call explicit_tend(mu, nu, tr(:, :, :, n), dt_tr(:, :, :, n))
      call compute_f(dt_tr(:, :, :, n), b, c, g, ftr)
      if (present(kbot)) then
        do j = 1, size(tr, 2)
        do i = 1, size(tr, 1)
          kb = kbot(i, j)
          f_delt_n1(i, j) = ftr(i, j, kb - 1)*delt
          delta_tr_n(i, j) = dt_tr(i, j, kb, n)*delt
        end do
        end do
      else
        f_delt_n1(:, :) = ftr(:, :, nlev - 1)*delt
        delta_tr_n(:, :) = dt_tr(:, :, nlev, n)*delt
      end if
      flux = 0.0
      call diff_surface(mu_delt_n, nu_n, e_n1, f_delt_n1, dflux_dtr, flux, 1.0, delta_tr_n)

! If flux needs to be saved then it should be made a module variable.
! vert_diff_init must allocate it and then call assign_tracer_field
! to set a pointer in tracer_manager_mod. It can be allocated as a
! 3 dimensional array with the 3'd index for tracer number.

      call vert_diff_up(delt, etr, ftr, delta_tr_n, dt_tr(:, :, :, n), kbot)
    end do

!-----------------------------------------------------------------------

  end subroutine tr_vert_diff

!#######################################################################

  !> Does the downward sweep of the tridiagonal solver for one field.
  subroutine vert_diff_down &
    (delt, mu, nu, tr, dt_tr, e, f, mu_delt_n, nu_n, &
     e_n1, f_delt_n1, delta_tr_n, kbot)

!-----------------------------------------------------------------------

    real, intent(in)                         :: delt
    real, intent(in), dimension(:, :, :)   :: mu, nu
    real, intent(in), dimension(:, :, :)   :: tr
    real, intent(inout), dimension(:, :, :)   :: dt_tr
    real, intent(out), dimension(:, :, :)   :: e
    real, intent(out), dimension(:, :, :)   :: f
    real, intent(out), dimension(:, :)     :: mu_delt_n, nu_n, e_n1
    real, intent(out), dimension(:, :)     :: f_delt_n1, delta_tr_n

    integer, intent(in), dimension(:, :), optional :: kbot

    real, dimension(size(tr, 1), size(tr, 2), size(tr, 3)) :: a, b, c, g

    integer :: i, j, k, kb, nlev

!-----------------------------------------------------------------------

    call explicit_tend(mu, nu, tr, dt_tr)

    call compute_e(delt, mu, nu, e, a, b, c, g)

    call compute_f(dt_tr, b, c, g, f)

    if (present(kbot)) then
      do j = 1, size(tr, 2)
      do i = 1, size(tr, 1)
        kb = kbot(i, j)
        mu_delt_n(i, j) = mu(i, j, kb)*delt
        nu_n(i, j) = nu(i, j, kb)
        e_n1(i, j) = e(i, j, kb - 1)
      end do
      end do
      do j = 1, size(tr, 2)
      do i = 1, size(tr, 1)
        kb = kbot(i, j)
        f_delt_n1(i, j) = f(i, j, kb - 1)*delt
        delta_tr_n(i, j) = dt_tr(i, j, kb)*delt
      end do
      end do
    else
      nlev = size(mu, 3)
      mu_delt_n(:, :) = mu(:, :, nlev)*delt
      nu_n(:, :) = nu(:, :, nlev)
      e_n1(:, :) = e(:, :, nlev - 1)
      f_delt_n1(:, :) = f(:, :, nlev - 1)*delt
      delta_tr_n(:, :) = dt_tr(:, :, nlev)*delt
    end if

!-----------------------------------------------------------------------

  end subroutine vert_diff_down

!#######################################################################

  !> Does the downward sweep of the tridiagonal solver for two fields.
  subroutine vert_diff_down_2 &
    (delt, mu, nu, xi_1, xi_2, dt_xi_1, dt_xi_2, e, f_1, f_2, &
     mu_delt_n, nu_n, e_n1, f_1_delt_n1, f_2_delt_n1, &
     delta_1_n, delta_2_n, kbot)

!-----------------------------------------------------------------------

    real, intent(in)                       :: delt
    real, intent(in), dimension(:, :, :) :: mu, nu, xi_1, xi_2
    real, intent(in), dimension(:, :, :) :: dt_xi_1, dt_xi_2
    real, intent(out), dimension(:, :, :) :: e, f_1, f_2
    real, intent(out), dimension(:, :)   :: mu_delt_n, nu_n, e_n1, &
                                            f_1_delt_n1, f_2_delt_n1, &
                                            delta_1_n, delta_2_n

    integer, intent(in), dimension(:, :), optional :: kbot

    real, dimension(size(xi_1, 1), size(xi_1, 2), size(xi_1, 3)) :: a, b, c, g, &
                                                                    dt_xi_11, dt_xi_22

    integer :: i, j, k, kb, nlev

!-----------------------------------------------------------------------

! local copy of input
    dt_xi_11 = dt_xi_1
    dt_xi_22 = dt_xi_2

    call explicit_tend(mu, nu, xi_1, dt_xi_11)
    call explicit_tend(mu, nu, xi_2, dt_xi_22)

    call compute_e(delt, mu, nu, e, a, b, c, g)

    call compute_f(dt_xi_11, b, c, g, f_1)
    call compute_f(dt_xi_22, b, c, g, f_2)

    if (present(kbot)) then
      do j = 1, size(xi_1, 2)
      do i = 1, size(xi_1, 1)
        kb = kbot(i, j)
        mu_delt_n(i, j) = mu(i, j, kb)*delt
        nu_n(i, j) = nu(i, j, kb)
        e_n1(i, j) = e(i, j, kb - 1)
        f_1_delt_n1(i, j) = f_1(i, j, kb - 1)*delt
        f_2_delt_n1(i, j) = f_2(i, j, kb - 1)*delt
        delta_1_n(i, j) = dt_xi_11(i, j, kb)*delt
        delta_2_n(i, j) = dt_xi_22(i, j, kb)*delt
      end do
      end do
    else
      nlev = size(mu, 3)
      mu_delt_n(:, :) = mu(:, :, nlev)*delt
      nu_n(:, :) = nu(:, :, nlev)
      e_n1(:, :) = e(:, :, nlev - 1)
      f_1_delt_n1(:, :) = f_1(:, :, nlev - 1)*delt
      f_2_delt_n1(:, :) = f_2(:, :, nlev - 1)*delt
      delta_1_n(:, :) = dt_xi_11(:, :, nlev)*delt
      delta_2_n(:, :) = dt_xi_22(:, :, nlev)*delt
    end if

!-----------------------------------------------------------------------

  end subroutine vert_diff_down_2

!#######################################################################

  !> Solves for the change of the lowest-level value including the surface flux, and
  !> updates the flux.
  subroutine diff_surface(mu_delt, nu, e_n1, f_delt_n1, &
                          dflux_datmos, flux, factor, delta_xi)

!-----------------------------------------------------------------------

    real, intent(in), dimension(:, :) :: mu_delt, nu, e_n1, f_delt_n1, &
                                         dflux_datmos
    real, intent(inout), dimension(:, :) :: flux, delta_xi
    real, intent(in) :: factor

!-----------------------------------------------------------------------

    real, dimension(size(flux, 1), size(flux, 2)) :: dflux
    real :: fff

    fff = 1.0/factor

    dflux = -nu*(1.0 - e_n1)
    delta_xi = delta_xi + mu_delt*nu*f_delt_n1

    delta_xi = (delta_xi + mu_delt*flux*fff)/ &
               (1.0 - mu_delt*(dflux + dflux_datmos*fff))

    flux = flux + dflux_datmos*delta_xi

!-----------------------------------------------------------------------

  end subroutine diff_surface

!#######################################################################

  !> Does the upward sweep of the tridiagonal solver.
  subroutine vert_diff_up(delt, e, f, delta_xi_n, dt_xi, kbot)

!-----------------------------------------------------------------------

    real, intent(in)                      :: delt
    real, intent(in), dimension(:, :, :) :: e, f
    real, intent(in), dimension(:, :)   :: delta_xi_n
    real, intent(out), dimension(:, :, :) :: dt_xi
    integer, intent(in), dimension(:, :), optional :: kbot

    integer :: i, j, k, kb, nlev
!-----------------------------------------------------------------------

    if (present(kbot)) then
      do j = 1, size(dt_xi, 2)
      do i = 1, size(dt_xi, 1)
        kb = kbot(i, j)
        dt_xi(i, j, kb) = delta_xi_n(i, j)/delt
        do k = kb - 1, 1, -1
          dt_xi(i, j, k) = e(i, j, k)*dt_xi(i, j, k + 1) + f(i, j, k)
        end do
      end do
      end do
    else
      nlev = size(dt_xi, 3)
      dt_xi(:, :, nlev) = delta_xi_n/delt
      do k = size(dt_xi, 3) - 1, 1, -1
        dt_xi(:, :, k) = e(:, :, k)*dt_xi(:, :, k + 1) + f(:, :, k)
      end do
    end if

!-----------------------------------------------------------------------

  end subroutine vert_diff_up

!#######################################################################

  !> Computes the tridiagonal coefficients and the `e` coefficients of the solver.
  subroutine compute_e(delt, mu, nu, e, a, b, c, g)

!-----------------------------------------------------------------------

    real, intent(in)                       :: delt
    real, intent(in), dimension(:, :, :) :: mu, nu
    real, intent(out), dimension(:, :, :) :: e, a, b, c, g

    integer :: k, nlev

!-----------------------------------------------------------------------

    nlev = size(mu, 3)

    a(:, :, 1:nlev - 1) = -mu(:, :, 1:nlev - 1)*nu(:, :, 2:nlev)*delt
    a(:, :, nlev) = 0.0
    c(:, :, 2:nlev) = -mu(:, :, 2:nlev)*nu(:, :, 2:nlev)*delt
    c(:, :, 1) = 0.0

    b = 1.0 - a - c

    e(:, :, 1) = -a(:, :, 1)/b(:, :, 1)
    do k = 2, nlev - 1
      g(:, :, k) = 1.0/(b(:, :, k) + c(:, :, k)*e(:, :, k - 1))
      e(:, :, k) = -a(:, :, k)*g(:, :, k)
    end do

!-----------------------------------------------------------------------

  end subroutine compute_e

!#######################################################################

  !> Computes the `f` coefficients of the solver for one field.
  subroutine compute_f(dt_xi, b, c, g, f)

!-----------------------------------------------------------------------
    real, intent(in), dimension(:, :, :) :: dt_xi, b, c, g
    real, intent(out), dimension(:, :, :) :: f

    integer :: k
!-----------------------------------------------------------------------

    f(:, :, 1) = dt_xi(:, :, 1)/b(:, :, 1)

    do k = 2, size(b, 3) - 1
      f(:, :, k) = (dt_xi(:, :, k) - c(:, :, k)*f(:, :, k - 1))*g(:, :, k)
    end do

!-----------------------------------------------------------------------

  end subroutine compute_f

!#######################################################################

  !> Adds the explicit diffusion tendency to `dt_xi`.
  subroutine explicit_tend(mu, nu, xi, dt_xi)

!-----------------------------------------------------------------------

    real, intent(in), dimension(:, :, :) :: mu, nu, xi
    real, intent(inout), dimension(:, :, :) :: dt_xi

    real, dimension(size(xi, 1), size(xi, 2), size(xi, 3)) :: fluxx

    integer :: k, nlev

!-----------------------------------------------------------------------

    nlev = size(mu, 3)

    fluxx(:, :, 1) = 0.0
    fluxx(:, :, 2:nlev) = nu(:, :, 2:nlev)*(xi(:, :, 2:nlev) - xi(:, :, 1:nlev - 1))

    dt_xi(:, :, 1:nlev - 1) = dt_xi(:, :, 1:nlev - 1) + &
                              mu(:, :, 1:nlev - 1)*(fluxx(:, :, 2:nlev) - fluxx(:, :, 1:nlev - 1))
    dt_xi(:, :, nlev) = dt_xi(:, :, nlev) - mu(:, :, nlev)*fluxx(:, :, nlev)

!-----------------------------------------------------------------------

  end subroutine explicit_tend

!#######################################################################

  !> Computes `grav` divided by the pressure thickness of each level.
  subroutine compute_mu(p_half, mu)

!-----------------------------------------------------------------------
    real, intent(in), dimension(:, :, :) :: p_half
    real, intent(out), dimension(:, :, :) :: mu

    integer :: nlev
!-----------------------------------------------------------------------

    nlev = size(mu, 3)

    mu(:, :, 1:nlev) = GRAV/(p_half(:, :, 2:nlev + 1) - p_half(:, :, 1:nlev))

!-----------------------------------------------------------------------

  end subroutine compute_mu

!#######################################################################

  !> Computes the air density at half levels times the diffusivity, divided by the level
  !> spacing.
  subroutine compute_nu(diff, p_half, p_full, z_full, t, q, nu)

!-----------------------------------------------------------------------
    real, intent(in), dimension(:, :, :) :: diff, p_half, p_full, &
                                            z_full, t, q
    real, intent(out), dimension(:, :, :) :: nu

    real, dimension(size(t, 1), size(t, 2), size(t, 3)) :: rho_half, tt
    integer :: nlev
!-----------------------------------------------------------------------

    nlev = size(nu, 3)

    if (use_virtual_temp_vert_diff) then
      tt = t*(1.0 + d608*q)           ! virtual temperature
    else
      tt = t ! Take out virtual temperature effect here to mimic supersource
    end if

    rho_half(:, :, 2:nlev) = &         ! density at half levels
      2.0*p_half(:, :, 2:nlev)/(RDGAS*(tt(:, :, 2:nlev) + tt(:, :, 1:nlev - 1)))

    nu(:, :, 2:nlev) = rho_half(:, :, 2:nlev)*diff(:, :, 2:nlev)/ &
                       (z_full(:, :, 1:nlev - 1) - z_full(:, :, 2:nlev))
!-----------------------------------------------------------------------

  end subroutine compute_nu

!#######################################################################

end module vert_diff_mod

