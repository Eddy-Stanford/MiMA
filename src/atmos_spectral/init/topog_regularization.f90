!> Regularization (smoothing) of the topography over the ocean, after Lindberg and
!> Broccoli (1996).
!>
!> The spectrally truncated topography `b` has Gibbs ripples over the ocean. The
!> regularized topography `a` minimizes
!>
!>     J(a) = |a - b|**2 + lambda * |M del2(a)|**2
!>
!> where `|.|` is the global area-weighted norm, `M` is the ocean mask and `del2` the
!> Laplacian on the unit sphere, so that it stays close to `b` everywhere but is smooth over
!> the ocean. `lambda` is chosen so that a fraction `fraction_target` of the ocean roughness,
!> `|M del2(b)|**2`, is removed.
!>
!> The minimum solves the normal equations `(1 + lambda L K L) a = b`, where `L` multiplies
!> each coefficient by `n(n+1)` and `K` applies the ocean mask in grid space. They are
!> solved by conjugate gradients with a diagonal preconditioner, and `lambda` by a search in
!> `log(lambda)` against `log(1 - fraction)` that starts from a guess scaled with the
!> truncation, so the method works at any resolution. The results do not depend on the
!> domain decomposition.
!>
!> The original GFDL code (MiMA 1.x) instead used a preconditioned fixed-point iteration
!> that also multiplied every zonal wavenumber `m > 0` by `sin(x)/x`, `x = pi m/(2
!> num_fourier)`, on each iteration, and a secant search from `lambda` = 1e-7. It did not
!> converge to the minimum of `J` (its mountains were lower) and failed at high resolution.
!>
!> Namelist: `topog_regularization_nml`
!> ([namelist reference](https://eddy-stanford.github.io/MiMA/Parameters/#topog_regularization_nml)).
!>
!> References:
!>
!> * Lindberg, C., and A. J. Broccoli, 1996: Representation of topography in spectral
!>   climate models and its effect on simulated precipitation. J. Climate, 9, 2641-2659.
!>
!> Original authors: Charles Jackson, Peter Phillipps.
module topog_regularization_mod

  use fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL, write_version_number, &
                     input_nml_file, check_nml_error, stdlog

  use mpp_mod, only: mpp_sum

  use transforms_mod, only: compute_legendre, trans_grid_to_spherical, trans_spherical_to_grid, &
                            get_sin_lat, get_wts_lat, transforms_are_initialized, get_lon_max, &
                            get_lat_max, get_num_fourier, get_fourier_inc, get_num_spherical, &
                            get_grid_domain, get_spec_domain, get_spherical_wave, &
                            area_weighted_global_mean

  implicit none
  private

  character(len=128), parameter :: version = '$Id: topog_regularization.f90 MiMA v2 $'
  character(len=128), parameter :: tagname = '$Name: MiMA v2 $'

  public :: regularize_topography

  real :: fraction_tolerance = 1.e-4
  !! the search for `lambda` stops when the fraction of the ocean roughness removed is within
  !! this of `ocean_topog_smoothing`
  real :: solver_tolerance = 1.e-10
  !! the conjugate-gradient solver stops when the norm of the residual is below this times the
  !! norm of the unsmoothed field
  integer :: max_solver_iterations = 2000  !! largest number of conjugate-gradient iterations
  integer :: max_lambda_iterations = 30  !! largest number of solves in the search for `lambda`
  logical :: verbose = .false.  !! print the progress of the solver and of the `lambda` search

  namelist /topog_regularization_nml/ fraction_tolerance, solver_tolerance, max_solver_iterations, &
    max_lambda_iterations, verbose

  integer :: is, ie, js, je, ms, me, ns, ne
  integer :: lon_max, lat_max, num_fourier, num_spherical, fourier_inc

  ! L = n(n+1), the diagonal D of the masked operator K, the weights of the coefficients in
  ! the area-weighted inner product (1 for m = 0, 2 for m > 0), and the ocean mask
  real, allocatable, dimension(:, :) :: eig, diag_mask
  real, allocatable, dimension(:) :: wt_m
  logical, allocatable, dimension(:, :) :: mask

contains

!============================================================================================

  !> Regularizes `field` over the ocean so that the fraction `fraction_target` of its ocean
  !> roughness is removed.
  subroutine regularize_topography(fraction_target, ocean_mask, field, smoothed_field, lambda, fraction_smoothed)

    real, intent(in) :: fraction_target  !! fraction of the ocean roughness to remove, 0 < `fraction_target` < 1
    logical, intent(in), dimension(:, :) :: ocean_mask  !! `.true.` at ocean points
    real, intent(in), dimension(:, :) :: field  !! field to smooth (the surface geopotential)
    real, intent(out), dimension(:, :) :: smoothed_field  !! regularized, spectrally truncated `field`
    real, intent(out) :: lambda  !! smoothing parameter found
    real, intent(out) :: fraction_smoothed  !! fraction of the ocean roughness removed

    integer :: io, ierr

    read (input_nml_file, nml=topog_regularization_nml, iostat=io)
    ierr = check_nml_error(io, 'topog_regularization_nml')
    call write_version_number(version, tagname)
    if (mpp_pe() == mpp_root_pe()) write (stdlog(), nml=topog_regularization_nml)

    if (.not. transforms_are_initialized()) &
      call error_mesg('regularize_topography', 'transforms are not initialized', FATAL)
    if (fraction_target <= 0. .or. fraction_target >= 1.) &
      call error_mesg('regularize_topography', 'ocean_topog_smoothing must be between 0 and 1', FATAL)

    call get_grid_domain(is, ie, js, je)
    call get_spec_domain(ms, me, ns, ne)
    call get_lon_max(lon_max)
    call get_lat_max(lat_max)
    call get_num_fourier(num_fourier)
    call get_num_spherical(num_spherical)
    call get_fourier_inc(fourier_inc)
    call check_shape(ocean_mask, field, smoothed_field)

    call regularization_init(ocean_mask)
    call find_lambda(fraction_target, field, smoothed_field, lambda, fraction_smoothed)
    call regularization_end

    if (mpp_pe() == mpp_root_pe()) &
      print '("Topography regularization: lambda=",1pe16.8,"  fraction_smoothed=",1pe16.8)', &
      lambda, fraction_smoothed

  end subroutine regularize_topography
!============================================================================================

  subroutine check_shape(ocean_mask, field, smoothed_field)
    logical, intent(in), dimension(:, :) :: ocean_mask
    real, intent(in), dimension(:, :) :: field, smoothed_field
    character(len=16) :: got, expected

    write (expected, '(2i6)') ie - is + 1, je - js + 1
    if (any(shape(ocean_mask) /= (/ie - is + 1, je - js + 1/)) .or. &
        any(shape(field) /= (/ie - is + 1, je - js + 1/)) .or. &
        any(shape(smoothed_field) /= (/ie - is + 1, je - js + 1/))) then
      write (got, '(2i6)') shape(field)
      call error_mesg('regularize_topography', 'arguments must have the shape of the grid domain,'// &
                      expected//'; field has'//got, FATAL)
    end if
  end subroutine check_shape

!============================================================================================

  !> Sets up the Laplacian eigenvalues, the inner-product weights and the diagonal of the
  !> masked operator, `D(m,n) = sum_j w_j f_j P_mn(mu_j)**2`, where `f_j` is the ocean
  !> fraction of latitude j.
  subroutine regularization_init(ocean_mask)
    logical, intent(in), dimension(:, :) :: ocean_mask

    integer, dimension(ms:me, ns:ne) :: wave
    real, dimension(lat_max) :: sin_lat, wts_lat, ocean_count
    real, dimension(0:num_fourier, 0:num_spherical, 1) :: legendre
    integer :: j, m

    allocate (eig(ms:me, ns:ne), diag_mask(ms:me, ns:ne), wt_m(ms:me), mask(is:ie, js:je))
    mask = ocean_mask

    call get_spherical_wave(wave)
    eig = real(wave)*real(wave + 1)
    do m = ms, me
      wt_m(m) = merge(1., 2., m == 0)
    end do

    ! number of ocean points on each latitude (integers, so the sum is exact)
    ocean_count = 0.
    do j = js, je
      ocean_count(j) = count(ocean_mask(:, j - js + 1))
    end do
    call mpp_sum(ocean_count, lat_max)

    call get_sin_lat(sin_lat)
    call get_wts_lat(wts_lat)
    diag_mask = 0.
    do j = 1, lat_max
      if (ocean_count(j) == 0.) cycle
      call compute_legendre(legendre, num_fourier, fourier_inc, num_spherical, sin_lat(j:j), 1)
      diag_mask = diag_mask + (wts_lat(j)*ocean_count(j)/lon_max)*legendre(ms:me, ns:ne, 1)**2
    end do

  end subroutine regularization_init
!--------------------------------------------------------------------------------------------

  subroutine regularization_end
    deallocate (eig, diag_mask, wt_m, mask)
  end subroutine regularization_end
!--------------------------------------------------------------------------------------------

  !> Finds `lambda` such that the fraction of the ocean roughness removed is within
  !> `fraction_tolerance` of `fraction_target`.
  !>
  !> `g(x) = log(1 - fraction) - log(1 - fraction_target)`, `x = log(lambda)`, decreases
  !> with `x` and is close to linear. Until `g` changes sign the search takes secant steps
  !> (from the first guess, a step that assumes a slope of -1); then it uses false position
  !> with the Illinois modification. Each solve starts from the solution of the previous one.
  subroutine find_lambda(fraction_target, field, smoothed_field, lambda, fraction_smoothed)
    real, intent(in) :: fraction_target
    real, intent(in), dimension(:, :) :: field
    real, intent(out), dimension(:, :) :: smoothed_field
    real, intent(out) :: lambda, fraction_smoothed

    complex, dimension(ms:me, ns:ne) :: field_spec, smooth_spec
    real :: rough_b, x, x_new, g, x_lo, g_lo, x_hi, g_hi, x_prev, g_prev
    integer :: it, n_eval, side
    logical :: have_lo, have_hi
    character(len=16) :: chtmp

    call trans_grid_to_spherical(field, field_spec)
    rough_b = roughness(field_spec)
    if (rough_b <= 0.) call error_mesg('regularize_topography', 'the topography is flat over the ocean, '// &
                                       'or there is no ocean; set ocean_topog_smoothing = 0', FATAL)

    ! within a factor 2 of the lambda that removes 99.5% at T42 to T341
    lambda = 800./(real(num_fourier*fourier_inc)*real(num_fourier*fourier_inc + 1))**2
    smooth_spec = field_spec/(1.+lambda*eig**2)

    have_lo = .false.
    have_hi = .false.
    side = 0
    x = log(lambda)
    x_prev = x
    g_prev = 0.
    do n_eval = 1, max_lambda_iterations
      lambda = exp(x)
      call cg_solve(lambda, field_spec, smooth_spec, it)
      fraction_smoothed = 1.-roughness(smooth_spec)/rough_b
      if (verbose .and. mpp_pe() == mpp_root_pe()) &
        print '("regularize_topography: lambda=",1pe16.8,"  fraction_smoothed=",1pe16.8,"  cg iterations=",i6)', &
        lambda, fraction_smoothed, it
      if (abs(fraction_smoothed - fraction_target) <= fraction_tolerance) then
        call trans_spherical_to_grid(smooth_spec, smoothed_field)
        return
      end if
      g = log(max(1.-fraction_smoothed, epsilon(1.))/(1.-fraction_target))

      if (g > 0.) then           ! too little smoothing: lambda must grow
        x_lo = x
        g_lo = g
        have_lo = .true.
      else                       ! too much smoothing
        x_hi = x
        g_hi = g
        have_hi = .true.
      end if

      if (have_lo .and. have_hi) then
        ! false position; halve the function value at the end that is kept twice (Illinois)
        if (g > 0.) then
          if (side == 1) g_hi = 0.5*g_hi
          side = 1
        else
          if (side == -1) g_lo = 0.5*g_lo
          side = -1
        end if
        x = x_lo + g_lo*(x_hi - x_lo)/(g_lo - g_hi)
      else if (n_eval > 1 .and. (g - g_prev)*(x - x_prev) < 0.) then
        ! no bracket yet: secant step through the last two points, at most a factor 1000
        x_new = x - g*(x - x_prev)/(g - g_prev)
        x_prev = x
        g_prev = g
        x = x + sign(min(abs(x_new - x), log(1000.)), g)
      else
        ! no bracket yet: step assuming a slope of -1, by a factor 1.1 to 1000
        x_prev = x
        g_prev = g
        x = x + sign(min(max(abs(g), 0.1), log(1000.)), g)
      end if
    end do

    write (chtmp, '(i8)') max_lambda_iterations
    call error_mesg('regularize_topography', 'no lambda found that removes the requested fraction of the roughness in'// &
                    trim(chtmp)//' solves; increase max_lambda_iterations or fraction_tolerance', FATAL)

  end subroutine find_lambda
!--------------------------------------------------------------------------------------------

  !> Solves `(1 + lambda L K L) a = b` by conjugate gradients with the diagonal (Jacobi)
  !> preconditioner `1/(1 + lambda L**2 D)`, starting from `a`.
  subroutine cg_solve(lambda, field_spec, smooth_spec, it)
    real, intent(in) :: lambda
    complex, intent(in), dimension(ms:me, ns:ne) :: field_spec
    complex, intent(inout), dimension(ms:me, ns:ne) :: smooth_spec
    integer, intent(out) :: it

    complex, dimension(ms:me, ns:ne) :: r, z, p, q
    real, dimension(ms:me, ns:ne) :: precond
    real :: rz, rz_old, alpha, rr, rr_stop
    character(len=16) :: chtmp

    precond = 1./(1.+lambda*eig**2*diag_mask)
    call apply_operator(lambda, smooth_spec, q)
    r = field_spec - q
    rr_stop = solver_tolerance**2*dot(field_spec, field_spec)
    rr = dot(r, r)
    z = precond*r
    p = z
    rz = dot(r, z)
    it = 0
    do while (rr > rr_stop)
      if (it == max_solver_iterations) then
        write (chtmp, '(i8)') max_solver_iterations
        call error_mesg('regularize_topography', 'conjugate gradients did not converge in'//trim(chtmp)// &
                        ' iterations; increase max_solver_iterations or solver_tolerance', FATAL)
      end if
      it = it + 1
      call apply_operator(lambda, p, q)
      alpha = rz/dot(p, q)
      smooth_spec = smooth_spec + alpha*p
      r = r - alpha*q
      rr = dot(r, r)
      z = precond*r
      rz_old = rz
      rz = dot(r, z)
      p = z + (rz/rz_old)*p
    end do
  end subroutine cg_solve
!--------------------------------------------------------------------------------------------

  !> `y = (1 + lambda L K L) x`.
  subroutine apply_operator(lambda, x, y)
    real, intent(in) :: lambda
    complex, intent(in), dimension(ms:me, ns:ne) :: x
    complex, intent(out), dimension(ms:me, ns:ne) :: y

    real, dimension(is:ie, js:je) :: grid

    call trans_spherical_to_grid(eig*x, grid)
    where (.not. mask) grid = 0.
    call trans_grid_to_spherical(grid, y)
    y = x + lambda*eig*y
  end subroutine apply_operator
!--------------------------------------------------------------------------------------------

  !> Area mean over the globe of the squared Laplacian of `x` (on the unit sphere) at the
  !> ocean points.
  real function roughness(x)
    complex, intent(in), dimension(ms:me, ns:ne) :: x

    real, dimension(is:ie, js:je) :: grid

    call trans_spherical_to_grid(eig*x, grid)
    where (mask)
      grid = grid**2
    elsewhere
      grid = 0.
    end where
    roughness = area_weighted_global_mean(grid)
  end function roughness
!--------------------------------------------------------------------------------------------

  !> Area-weighted inner product of two real fields given by their coefficients. The sum
  !> over n is taken on each process, and the sum over m in a fixed order, so the result does
  !> not depend on the decomposition.
  real function dot(x, y)
    complex, intent(in), dimension(ms:me, ns:ne) :: x, y

    real, dimension(0:num_fourier) :: partial
    integer :: m

    partial = 0.
    do m = ms, me
      partial(m) = wt_m(m)*sum(real(x(m, :))*real(y(m, :)) + aimag(x(m, :))*aimag(y(m, :)))
    end do
    call mpp_sum(partial, num_fourier + 1)
    dot = sum(partial)
  end function dot

end module topog_regularization_mod
