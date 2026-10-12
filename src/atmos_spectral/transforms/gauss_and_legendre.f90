!> Gaussian latitudes and weights, and the associated Legendre functions.
module gauss_and_legendre_mod

  use fms_mod, only: error_mesg, FATAL

  use constants_mod, only: pi

  implicit none
  private

  public :: compute_legendre, compute_gaussian

contains

!-----------------------------------------------------------------------
  !> Computes the normalized associated Legendre functions at the given latitudes by recursion.
  !>
  !> `legendre(m, n, j)` is the function of zonal wavenumber `M = m*fourier_inc` and degree
  !> `M + n` at latitude j, normalized so that its square integrates to 1 over
  !> sin(latitude) from -1 to 1.
  subroutine compute_legendre(legendre, num_fourier, fourier_inc, &
                              num_spherical, sin_lat, n_lat)
!-----------------------------------------------------------------------

    integer, intent(in) :: num_fourier, fourier_inc, num_spherical, n_lat
    !! `num_fourier`, `num_spherical`: largest indices m and n; `fourier_inc`: zonal wavenumber
    !! increment; `n_lat`: number of latitudes
    real, intent(in), dimension(n_lat) :: sin_lat  !! sines of the latitudes

    real, intent(out), dimension(0:num_fourier, 0:num_spherical, n_lat)  :: legendre  !! associated Legendre functions

    integer :: j, m, n, fourier_max
    real, dimension(0:num_fourier*fourier_inc, 0:num_spherical) :: poly, eps, m2, l2
    real, dimension(0:num_fourier*fourier_inc) :: b
    real, dimension(n_lat) :: cos_lat

    fourier_max = num_fourier*fourier_inc

    do j = 1, n_lat
      cos_lat(j) = sqrt(1 - sin_lat(j)*sin_lat(j))
    end do

    do n = 0, num_spherical
      do m = 0, fourier_max
        m2(m, n) = m**2
        l2(m, n) = (m + n)**2
      end do
    end do

    eps = sqrt((l2 - m2)/(4.0*l2 - 1.0))

    do m = 1, fourier_max
      b(m) = sqrt(0.5*(2.0*real(m) + 1.0)/real(m))
    end do

    poly(0, 0) = sqrt(0.5)
    do j = 1, n_lat
      do m = 1, fourier_max
        poly(m, 0) = b(m)*cos_lat(j)*poly(m - 1, 0)
      end do
      do m = 0, fourier_max
        poly(m, 1) = sin_lat(j)*poly(m, 0)/eps(m, 1)
      end do

      do n = 2, num_spherical
        do m = 0, fourier_max
          poly(m, n) = (sin_lat(j)*poly(m, n - 1) - eps(m, n - 1)*poly(m, n - 2))/eps(m, n)
        end do
      end do

      do n = 0, num_spherical
        do m = 0, num_fourier
          legendre(m, n, j) = poly(m*fourier_inc, n)
        end do
      end do
    end do

    return
  end subroutine compute_legendre

!----------------------------------------------------------------------
  !> Computes the Gaussian latitudes and weights of one hemisphere by Newton iteration
  !> (Press et al., 1990, Numerical Recipes).
  !>
  !> The latitudes are the roots of the Legendre polynomial of degree `2*n_hem`, ordered from
  !> the pole towards the equator; the weights of the whole globe sum to 2.
  subroutine compute_gaussian(sin_hem, wts_hem, n_hem)
!----------------------------------------------------------------------

    integer, intent(in) :: n_hem  !! number of latitudes in one hemisphere
    real, intent(out), dimension(n_hem) :: sin_hem, wts_hem
    !! `sin_hem`: sines of the Gaussian latitudes; `wts_hem`: Gaussian weights

    real :: converg
    integer :: itermax
    integer :: i, iter, j, n, nprec
    real :: pp, p1, p2, p3, z, z1
    logical :: converged

    ! convergence criterion from the decimal precision of the reals (1e-15 for 64-bit reals)
    nprec = precision(converg)
    converg = .1**nprec

    itermax = 10

    n = 2*n_hem
    do i = 1, n_hem
      z = cos(pi*(i - 0.25)/(n + 0.5))
      converged = .false.
      do iter = 1, itermax
        p1 = 1.0
        p2 = 0.0

        do j = 1, n
          p3 = p2
          p2 = p1
          p1 = ((2.0*j - 1.0)*z*p2 - (j - 1.0)*p3)/j
        end do

        pp = n*(z*p1 - p2)/(z*z - 1.0e+00)
        z1 = z
        z = z1 - p1/pp
        if (abs(z - z1) < converg) then
          converged = .true.
          exit
        end if
      end do
      if (.not. converged) &
        call error_mesg('compute_gaussian', 'abscissas failed to converge in itermax iterations', FATAL)

      sin_hem(i) = z
      wts_hem(i) = 2.0/((1.0 - z*z)*pp*pp)

    end do

    return
  end subroutine compute_gaussian
!----------------------------------------------------------------------

end module gauss_and_legendre_mod
