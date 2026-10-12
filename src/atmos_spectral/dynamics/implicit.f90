!> Semi-implicit correction of the gravity-wave terms of the spectral tendencies.
!>
!> The divergence, temperature and log surface pressure equations are linearized about a
!> reference state (a temperature profile, 300 K at all levels as set by `spectral_dynamics`,
!> and a surface pressure, `reference_sea_level_press`). The gravity-wave terms are treated
!> implicitly with weight `alpha` (`alpha_implicit` in `spectral_dynamics_nml`): for each
!> total wavenumber L the divergence tendency is multiplied by the inverse of
!> `I + (alpha*dt)**2 * L(L+1)/a**2 * D`, where D couples the levels, and the temperature and
!> surface-pressure tendencies are corrected to match.
!> `implicit_init` works with global fields, the rest on the spectral domain.
module implicit_mod
!||zn, pjk 5/13/99
!following vb's changes to old implicit_mod.

  use fms_mod, only: error_mesg, FATAL, write_version_number

  use constants_mod, only: rdgas, radius, cp_air

  use press_and_geopot_mod, only: pressure_variables

  use matrix_invert_mod, only: invert

  use transforms_mod, only: get_spec_domain

  implicit none

  private

  public :: implicit_init, implicit_correction, implicit_end

  interface linear_geopotential
    module procedure linear_geopotential_1d, linear_geopotential_3d
  end interface

  interface linear_tp_tendency
    module procedure linear_tp_tendency_1d, linear_tp_tendency_3d
  end interface

  logical :: module_is_initialized = .false.

  character(len=128), parameter :: version = '$Id: implicit.f90,v 11.0 2004/09/28 19:29:40 fms Exp $'
  character(len=128), parameter :: tagname = '$Name: lima $'

  real, allocatable, dimension(:)   :: ref_temperature_implicit
  real, allocatable, dimension(:)   :: ref_ln_p_half, ref_ln_p_full, del_ln_p_half, del_ln_p_full
  real, allocatable, dimension(:, :) :: eigen
  integer, allocatable, dimension(:, :) :: wavenumber

  real :: alpha, ref_surf_p_implicit

  real, allocatable, dimension(:, :)   :: div_mat
  real, allocatable, dimension(:)     :: h
  real, allocatable, dimension(:, :, :) :: wave_matrix
  real, allocatable, dimension(:)     :: pk, bk, dpk, dbk

  real :: dt = 0.0
  real :: xi

  integer :: num_levels, num_total_wavenumbers
  character(len=64) :: vert_difference_option
  integer :: ms, me, ns, ne

contains

!------------------------------------------------------------------------

  !> Stores the vertical coordinate and the reference state, and builds the matrix that
  !> couples the levels in the linearized divergence equation.
  subroutine implicit_init(pk_in, bk_in, ref_temperature_implicit_in, ref_surf_p_implicit_in, &
                           num_total_wavenumbers_in, eigen_in, wavenumber_in, alpha_in, &
                           vert_difference_option_in)

    real, intent(in), dimension(:)     :: pk_in, bk_in, ref_temperature_implicit_in
    !! `pk_in` [Pa], `bk_in`: coefficients of the half-level pressures `pk + bk*ps`;
    !! `ref_temperature_implicit_in`: reference temperature at each level [K]
    real, intent(in), dimension(0:, 0:) :: eigen_in  !! eigenvalues L(L+1)/a**2 of minus the Laplacian [1/m2]
    integer, intent(in), dimension(0:, 0:) :: wavenumber_in  !! total wavenumber L of each spectral coefficient

    real, intent(in) :: ref_surf_p_implicit_in, alpha_in
    !! `ref_surf_p_implicit_in`: reference surface pressure [Pa]; `alpha_in`: implicitness
    !! (0.5: centred, 1: backward)
    integer, intent(in) :: num_total_wavenumbers_in  !! largest total wavenumber
    character(len=*), intent(in) :: vert_difference_option_in  !! vertical differencing (`'simmons_and_burridge'`)

    real, dimension(size(ref_temperature_implicit_in, 1)) :: p_full_work_1, p_full_work_2, ln_p_full_work_1, ln_p_full_work_2
    real, dimension(size(ref_temperature_implicit_in, 1) + 1) :: p_half_work_1, p_half_work_2, ln_p_half_work_1, ln_p_half_work_2

    real :: surface_p_1, eps = 1.e-5
    integer :: k

    if (module_is_initialized) return

    call write_version_number(version, tagname)

    call get_spec_domain(ms, me, ns, ne)

    num_levels = size(ref_temperature_implicit_in, 1)
    num_total_wavenumbers = num_total_wavenumbers_in

    allocate (ref_temperature_implicit(num_levels))
    allocate (ref_ln_p_half(num_levels + 1))
    allocate (ref_ln_p_full(num_levels))
    allocate (del_ln_p_half(num_levels + 1))
    allocate (del_ln_p_full(num_levels))
    allocate (div_mat(num_levels, num_levels))
    allocate (wave_matrix(num_levels, num_levels, 0:num_total_wavenumbers))
    allocate (h(num_levels))
    allocate (eigen(0:size(eigen_in, 1) - 1, 0:size(eigen_in, 2) - 1))
    allocate (wavenumber(0:size(eigen_in, 1) - 1, 0:size(eigen_in, 2) - 1))

    ref_temperature_implicit = ref_temperature_implicit_in
    ref_surf_p_implicit = ref_surf_p_implicit_in
    eigen = eigen_in
    wavenumber = wavenumber_in
    alpha = alpha_in
    vert_difference_option = vert_difference_option_in

    allocate (pk(size(pk_in(:))), bk(size(bk_in, 1)))
    allocate (dpk(num_levels), dbk(num_levels))

    pk = pk_in
    bk = bk_in

    do k = 1, num_levels
      dpk(k) = pk(k + 1) - pk(k)
      dbk(k) = bk(k + 1) - bk(k)
    end do

    call pressure_variables(p_half_work_1, ref_ln_p_half, p_full_work_1, ref_ln_p_full, ref_surf_p_implicit)

!del_ln_p_half = dln_p_half_dps(ref_surf_p_implicit)
!del_ln_p_full = dln_p_full_dps(ref_surf_p_implicit)

!  functions dln_p_half_dps and del_ln_p_full are used nowhere else and
!  they perform very simple calculations, so it only makes the code more
!  confusing to do these calculations in functions.

    do k = 2, num_levels + 1
      del_ln_p_half(k) = bk(k)/(pk(k) + bk(k)*ref_surf_p_implicit)
    end do

    if (pk(1) .eq. 0.0) then
      del_ln_p_half(1) = 1.0/ref_surf_p_implicit
    else
      del_ln_p_half(1) = bk(1)/(pk(1) + bk(1)*ref_surf_p_implicit)
    end if

    surface_p_1 = ref_surf_p_implicit*(1.0 - 0.5*eps)
    call pressure_variables(p_half_work_1, ln_p_half_work_1, p_full_work_1, ln_p_full_work_1, surface_p_1)

    surface_p_1 = ref_surf_p_implicit*(1.0 + 0.5*eps)
    call pressure_variables(p_half_work_2, ln_p_half_work_2, p_full_work_2, ln_p_full_work_2, surface_p_1)

    del_ln_p_full = (ln_p_full_work_2 - ln_p_full_work_1)/(eps*ref_surf_p_implicit)

    call build_matrix

    module_is_initialized = .true.

    return
  end subroutine implicit_init

!------------------------------------------------------------------------

  subroutine build_matrix

    real, dimension(num_levels) :: input, zero, dt_t
    real, dimension(num_levels + 1) :: zero_1
    real :: dt_p

    real, dimension(num_levels, num_levels) :: gamma, tau
    real, dimension(num_levels) :: nu, h1, h2

    integer :: k, kk, kkk

    zero = 0.0
    zero_1 = 0.0

    tau = 0.0
    nu = 0.0
    do k = 1, num_levels
      input = 0.0
      input(k) = 1.0

      call linear_tp_tendency(input, ref_temperature_implicit, ref_surf_p_implicit, &
                              ref_ln_p_half, ref_ln_p_full, dt_p, dt_t)
      nu(k) = -dt_p
      tau(:, k) = -dt_t

      gamma(:, k) = linear_geopotential(input, zero_1, zero, &
                                        ref_temperature_implicit, ref_ln_p_half, ref_ln_p_full)
    end do

    h1 = pres_grad_funct(ref_temperature_implicit, ref_ln_p_half, ref_ln_p_full, ref_surf_p_implicit)

    h2 = linear_geopotential(zero, del_ln_p_half, del_ln_p_full, &
                             ref_temperature_implicit, ref_ln_p_half, ref_ln_p_full)

    h = h1 + h2

    do k = 1, num_levels
      do kk = 1, num_levels
        div_mat(k, kk) = h(k)*nu(kk)
        do kkk = 1, num_levels
          div_mat(k, kk) = div_mat(k, kk) + gamma(k, kkk)*tau(kkk, kk)
        end do
      end do
    end do

    return
  end subroutine build_matrix

!-------------------------------------------------------------------------

  !> Builds and inverts the implicit matrix for each total wavenumber.
  subroutine build_wave_matrices

    real :: factor, det
    integer :: k, L

    wave_matrix = 0.0
    do L = 0, num_total_wavenumbers
      factor = xi*xi*L*(L + 1)/radius**2
      do k = 1, num_levels
        wave_matrix(k, k, L) = 1.0
      end do
      wave_matrix(:, :, L) = wave_matrix(:, :, L) + factor*div_mat(:, :)
      call invert(wave_matrix(:, :, L), det)
    end do

    return
  end subroutine build_wave_matrices

!---------------------------------------------------------------------

  !> Replaces the explicit tendencies of divergence, temperature and log surface pressure by
  !> the semi-implicit ones.
  !>
  !> The matrices for each total wavenumber are rebuilt when the time step changes.
  subroutine implicit_correction(dt_divs, dt_ts, dt_ln_ps, divs, ts, ln_ps, dt_in, previous, current)
!in parallel, wavenumber (0,0) actually means (ms,ns)
!it is necessary to use the (0:,0:) dimensioning here because of the m-n loop below
    complex, intent(inout), dimension(0:, 0:, :) :: dt_divs, dt_ts
    !! `dt_divs`, `dt_ts`: spectral tendencies of divergence [1/s2] and temperature [K/s]
    complex, intent(inout), dimension(0:, 0:) :: dt_ln_ps  !! spectral tendency of log surface pressure [1/s]
    complex, intent(in), dimension(0:, 0:, :, :) :: divs, ts
    !! `divs`, `ts`: spectral divergence [1/s] and temperature [K] (last dimension: time level)
    complex, intent(in), dimension(0:, 0:, :) :: ln_ps  !! spectral log surface pressure (last dimension: time level)
    real, intent(in) :: dt_in  !! time step of the leapfrog step (from `previous` to the future level) [s]
    integer, intent(in) :: previous, current  !! `previous`, `current`: indices of the previous and current time levels

    complex, dimension(0:size(divs, 1) - 1, 0:size(divs, 2) - 1, size(divs, 3)) :: dt_ts_temp
    complex, dimension(num_levels) ::  work
    complex, dimension(0:size(divs, 1) - 1, 0:size(divs, 2) - 1) :: dt_ps_temp
    integer :: m, n, L, kk

    if (.not. module_is_initialized) then
      call error_mesg('implicit_correction', 'failed to initialize implicit', FATAL)
    end if

    if (dt_in .ne. dt) then
      dt = dt_in
      xi = dt*alpha
      call build_wave_matrices
    end if

    call adjust_dt_divs(dt_divs, dt_ts, dt_ln_ps, divs, ts, ln_ps, previous, current)

    do n = 0, size(divs, 2) - 1
      do m = 0, size(divs, 1) - 1
!||zn, pjk 5/13/99, following vb's dmsm/implicit.f90
        L = wavenumber(m + ms, n + ns)
        if (L .le. num_total_wavenumbers) then
          work = matmul(wave_matrix(:, :, L), dt_divs(m, n, :))
          dt_divs(m, n, :) = work
        end if
      end do
    end do

    call linear_tp_tendency(dt_divs, ref_temperature_implicit, ref_surf_p_implicit, &
                            ref_ln_p_half, ref_ln_p_full, dt_ps_temp, dt_ts_temp)

    dt_ts = dt_ts + xi*dt_ts_temp
    dt_ln_ps = dt_ln_ps + xi*dt_ps_temp/ref_surf_p_implicit

    return
  end subroutine implicit_correction

!---------------------------------------------------------------------
  !> Adds the implicit gravity-wave terms of the previous and current levels to the
  !> divergence tendency.
  subroutine adjust_dt_divs(dt_divs, dt_ts, dt_ln_ps, divs, ts, ln_ps, previous, current)
    !||zn, pjk 5/13/99, following vb's dmsm/implicit.f90. (0:,0:) -> (:,:)
    complex, intent(inout), dimension(:, :, :) :: dt_divs
    complex, intent(inout), dimension(:, :, :) :: dt_ts
    complex, intent(inout), dimension(:, :) :: dt_ln_ps
    complex, intent(in), dimension(:, :, :, :) :: divs, ts
    complex, intent(in), dimension(:, :, :) :: ln_ps
    integer, intent(in) :: previous, current

    complex, dimension(size(ts, 1), size(ts, 2), size(ts, 3)) :: divs_temp, ts_temp, geopot, dt_ts_temp, ln_p_full_temp
    complex, dimension(size(ts, 1), size(ts, 2), size(ts, 3) + 1) :: ln_p_half_temp
    complex, dimension(size(ts, 1), size(ts, 2)) :: ps_temp, dt_ps_temp

    integer :: k

    divs_temp = divs(:, :, :, previous) - divs(:, :, :, current)

    call linear_tp_tendency(divs_temp, ref_temperature_implicit, ref_surf_p_implicit, &
                            ref_ln_p_half, ref_ln_p_full, dt_ps_temp, dt_ts_temp)

    dt_ts = dt_ts + dt_ts_temp
    dt_ln_ps = dt_ln_ps + dt_ps_temp/ref_surf_p_implicit

    ts_temp = ts(:, :, :, previous) - ts(:, :, :, current) + xi*dt_ts
    ps_temp = ln_ps(:, :, previous) - ln_ps(:, :, current) + xi*dt_ln_ps

    ln_p_full_temp = (0., 0.)
    ln_p_half_temp = (0., 0.)

    geopot = linear_geopotential(ts_temp, ln_p_half_temp, ln_p_full_temp, ref_temperature_implicit, ref_ln_p_half, ref_ln_p_full)

    do k = 1, num_levels
      dt_divs(:, :, k) = dt_divs(:, :, k) + eigen(ms:me, ns:ne)*(geopot(:, :, k) + h(k)*ps_temp(:, :)*ref_surf_p_implicit)
    end do

    return
  end subroutine adjust_dt_divs

!-----------------------------------------------------------------------

  function linear_geopotential_3d(del_t, del_ln_p_half, del_ln_p_full, t, ln_p_half, ln_p_full) result(geopot)
!||zn, pjk 5/13/99, following vb's dmsm/implicit.f90. (0:,0:) -> (:,:)

    complex, intent(in), dimension(:, :, :) :: del_t
    complex, intent(in), dimension(:, :, :) :: del_ln_p_half
    complex, intent(in), dimension(:, :, :) :: del_ln_p_full
    real, intent(in), dimension(:)       :: t
    real, intent(in), dimension(:)       :: ln_p_half
    real, intent(in), dimension(:)       :: ln_p_full

    complex, dimension(size(del_t, 1), size(del_t, 2), size(del_t, 3))   :: geopot
    complex, dimension(size(del_t, 1), size(del_t, 2), size(del_t, 3) + 1) :: geopot_half

    integer :: k

    geopot_half(:, :, num_levels + 1) = (0., 0.)

    do k = num_levels, 2, -1
      geopot_half(:, :, k) = geopot_half(:, :, k + 1) &
                             + rdgas*(del_t(:, :, k)*(ln_p_half(k + 1) - ln_p_half(k)) &
                                      + t(k)*(del_ln_p_half(:, :, k + 1) - del_ln_p_half(:, :, k)))
    end do

    do k = 1, num_levels
      geopot(:, :, k) = geopot_half(:, :, k + 1) &
                        + rdgas*(del_t(:, :, k)*(ln_p_half(k + 1) - ln_p_full(k)) &
                                 + t(k)*(del_ln_p_half(:, :, k + 1) - del_ln_p_full(:, :, k)))
    end do

    return
  end function linear_geopotential_3d
!-----------------------------------------------------------------------

  function linear_geopotential_1d(del_t, del_ln_p_half, del_ln_p_full, t, ln_p_half, ln_p_full) result(geopot)

    real, intent(in), dimension(:) :: del_t
    real, intent(in), dimension(:) :: del_ln_p_half
    real, intent(in), dimension(:) :: del_ln_p_full
    real, intent(in), dimension(:) :: t
    real, intent(in), dimension(:) :: ln_p_half
    real, intent(in), dimension(:) :: ln_p_full

    real, dimension(size(del_t, 1)) :: geopot

    complex, dimension(0:0, 0:0, size(del_t, 1)) :: del_t_3d, geopot_3d
    complex, dimension(0:0, 0:0, size(del_t, 1) + 1) :: del_ln_p_half_3d
    complex, dimension(0:0, 0:0, size(del_t, 1)) :: del_ln_p_full_3d

    del_t_3d(0, 0, :) = cmplx(del_t, 0.0)
    del_ln_p_half_3d(0, 0, :) = cmplx(del_ln_p_half, 0.0)
    del_ln_p_full_3d(0, 0, :) = cmplx(del_ln_p_full, 0.0)

    geopot_3d = linear_geopotential(del_t_3d, del_ln_p_half_3d, del_ln_p_full_3d, t, ln_p_half, ln_p_full)

    geopot = real(geopot_3d(0, 0, :))

    return
  end function linear_geopotential_1d
!-----------------------------------------------------------------------

  function pres_grad_funct(tg, ln_p_half, ln_p_full, p_surf) result(x)

    real, intent(in), dimension(:) :: tg, ln_p_half, ln_p_full
    real, intent(in)                :: p_surf
    real, dimension(size(tg, 1)) :: x

    real :: dlog_1, dlog_2
    integer :: k

    if (trim(vert_difference_option) == 'simmons_and_burridge') then
      do k = 1, size(ln_p_full, 1)
        dlog_1 = ln_p_half(k + 1) - ln_p_full(k)
        dlog_2 = ln_p_full(k) - ln_p_half(k)
        x(k) = rdgas*tg(k)*(bk(k + 1)*dlog_1 + bk(k)*dlog_2)/(dpk(k) + dbk(k)*p_surf)
      end do
    end if

    return
  end function pres_grad_funct
!-----------------------------------------------------------------------

  subroutine linear_tp_tendency_3d(div, t_ref, p_surf_ref, ln_p_half_ref, ln_p_full_ref, dt_p_surf, dt_t)

    complex, intent(in), dimension(0:, 0:, :) :: div
    real, intent(in), dimension(:) :: t_ref, ln_p_half_ref, ln_p_full_ref
    real, intent(in)  :: p_surf_ref

    complex, intent(out), dimension(0:, 0:)   :: dt_p_surf
    complex, intent(out), dimension(0:, 0:, :) :: dt_t

    real ::  dp, dp_inv, dlog_1, dlog_3, p_full_ref
    complex, dimension(0:size(div, 1) - 1, 0:size(div, 2) - 1) :: dmean, dmean_tot
    complex, dimension(0:size(div, 1) - 1, 0:size(div, 2) - 1, size(div, 3) + 1) :: vert_vel, temp

    real :: kappa
    integer :: k

    kappa = rdgas/cp_air

    dmean_tot = 0.0

    if (vert_difference_option == 'simmons_and_burridge') then

      do k = 1, num_levels
        dp = dpk(k) + dbk(k)*p_surf_ref
        dp_inv = 1/dp
        dlog_1 = ln_p_half_ref(k + 1) - ln_p_full_ref(k)
        dlog_3 = ln_p_half_ref(k + 1) - ln_p_half_ref(k)
        dmean = div(:, :, k)*dp
        dt_t(:, :, k) = -kappa*t_ref(k)*(dmean_tot*dlog_3 + dmean*dlog_1)*dp_inv
        dmean_tot = dmean_tot + dmean
        vert_vel(:, :, k + 1) = -dmean_tot
      end do

    end if

    dt_p_surf = -dmean_tot

    do k = 2, num_levels
      vert_vel(:, :, k) = vert_vel(:, :, k) + dmean_tot*bk(k)
    end do

    do k = 2, num_levels
      temp(:, :, k) = -vert_vel(:, :, k)*(t_ref(k) - t_ref(k - 1))
    end do

    temp(:, :, 1) = 0.0
    temp(:, :, num_levels + 1) = 0.0

    do k = 1, num_levels
      dp = dpk(k) + dbk(k)*p_surf_ref
      dp_inv = 1/dp
      dt_t(:, :, k) = dt_t(:, :, k) + .5*dp_inv*(temp(:, :, k + 1) + temp(:, :, k))
    end do

    return
  end subroutine linear_tp_tendency_3d
!-----------------------------------------------------------------------

  subroutine linear_tp_tendency_1d(div, t_ref, p_surf_ref, ln_p_half_ref, ln_p_full_ref, dt_p_surf, dt_t)

    real, intent(in), dimension(:) :: div
    real, intent(in), dimension(:) :: t_ref, ln_p_half_ref, ln_p_full_ref
    real, intent(in)  :: p_surf_ref

    real, intent(out) :: dt_p_surf
    real, intent(out), dimension(:) :: dt_t

    complex, dimension(0:0, 0:0, size(div, 1)) :: div_3d, dt_t_3d
    complex, dimension(0:0, 0:0)           :: dt_p_surf_3d

    div_3d(0, 0, :) = cmplx(div, 0.)

    call linear_tp_tendency(div_3d, t_ref, p_surf_ref, ln_p_half_ref, ln_p_full_ref, dt_p_surf_3d, dt_t_3d)

    dt_t = real(dt_t_3d(0, 0, :))
    dt_p_surf = real(dt_p_surf_3d(0, 0))

    return
  end subroutine linear_tp_tendency_1d
!-----------------------------------------------------------------------
  !> Deallocates the module arrays.
  subroutine implicit_end

    if (.not. module_is_initialized) return

    deallocate (ref_temperature_implicit)
    deallocate (ref_ln_p_half, ref_ln_p_full, del_ln_p_half, del_ln_p_full)
    deallocate (div_mat, wave_matrix, h, eigen, wavenumber, pk, bk, dpk, dbk)
    module_is_initialized = .false.

    return
  end subroutine implicit_end
!-----------------------------------------------------------------------

end module implicit_mod
