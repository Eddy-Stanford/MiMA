
!> Fast Fourier transforms between real grid space and complex Fourier space, for many
!> sequences at once.
!>
!> Computes multiple 1-D FFTs and inverse FFTs of 2-D and 3-D arrays, between real
!> gridpoint values and complex Fourier coefficients, with Temperton's FFT of `fft99_mod`.
!> `fft_init` sets the length `n` of the transforms. The complex Fourier components are
!> stored as
!>
!> ```text
!> fourier(1)     = cmplx(a(0), b(0))
!> fourier(2)     = cmplx(a(1), b(1))
!>     ...
!> fourier(n/2+1) = cmplx(a(n/2), b(n/2))
!> ```
!>
!> The 32-bit and 64-bit interfaces both do the transform at the default real precision
!> (64-bit in MiMA).
!>
!> Original authors: Bruce Wyman.
module fft_mod

  use platform_mod, only: R8_KIND, R4_KIND
  use fms_mod, only: error_mesg, FATAL
  use fft99_mod, only: fft991, set99

  implicit none
  private

  public :: fft_init, fft_end, fft_grid_to_fourier, fft_fourier_to_grid

  !> Transforms multiple sequences of real gridpoint values to complex Fourier
  !> coefficients.
  !>
  !> `fourier = fft_grid_to_fourier(grid)`, for 2-D or 3-D arrays of 32- or 64-bit reals.
  !> Each sequence along the first dimension of `grid` has the length `n` set by
  !> `fft_init`; that dimension must be at least `n` (only the first `n` values are used).
  !> The first dimension of `fourier` is n/2+1 and the other dimensions are those of
  !> `grid`. Stops if `fft_init` has not been called.
  interface fft_grid_to_fourier
    module procedure fft_grid_to_fourier_float_2d, fft_grid_to_fourier_double_2d, &
      fft_grid_to_fourier_float_3d, fft_grid_to_fourier_double_3d
  end interface

  !> Transforms multiple sequences of complex Fourier coefficients to real gridpoint
  !> values (the inverse of `fft_grid_to_fourier`).
  !>
  !> `grid = fft_fourier_to_grid(fourier)`, for 2-D or 3-D arrays of 32- or 64-bit
  !> complex values. The first dimension of `fourier` must be at least n/2+1 (only the
  !> first n/2+1 values are used), where `n` is the length set by `fft_init`. The first
  !> dimension of `grid` is `n`; the other dimensions are those of `fourier`. Stops if
  !> `fft_init` has not been called.
  interface fft_fourier_to_grid
    module procedure fft_fourier_to_grid_float_2d, fft_fourier_to_grid_double_2d, &
      fft_fourier_to_grid_float_3d, fft_fourier_to_grid_double_3d
  end interface

  real, allocatable, dimension(:) :: table99   ! trigonometric values for fft991, from set99
  integer, allocatable, dimension(:) :: ifax   ! factors of n/2, from set99

  integer :: leng, leng1, leng2, lenc    ! n, n+1, n+2 and n/2+1

  logical :: module_is_initialized = .false.

contains

!#######################################################################

  !> Transforms 2-D 32-bit gridpoint data to Fourier space (see `fft_grid_to_fourier`).
  function fft_grid_to_fourier_float_2d(grid) result(fourier)

    real(R4_KIND), intent(in), dimension(:, :)  :: grid
    complex(R4_KIND), dimension(lenc, size(grid, 2)) :: fourier

    real, dimension(leng2, size(grid, 2)) :: data
    real, dimension(leng1, size(grid, 2)) :: work
    integer :: j, k, num

    if (.not. module_is_initialized) &
      call error_handler('fft_grid_to_fourier', 'fft_init must be called.')

    if (size(grid, 1) < leng) &
      call error_handler('fft_grid_to_fourier', 'length of input data too small.')

    num = size(grid, 2)    ! number of transforms

    do j = 1, num
      data(1:leng, j) = grid(1:leng, j)
    end do
    call fft991(data, work, table99, ifax, 1, leng2, leng, num, -1)
    do j = 1, size(grid, 2)
    do k = 1, lenc
      fourier(k, j) = cmplx(data(2*k - 1, j), data(2*k, j))
    end do
    end do

  end function fft_grid_to_fourier_float_2d

!#######################################################################

  !> Transforms 2-D 32-bit Fourier data to grid space (see `fft_fourier_to_grid`).
  function fft_fourier_to_grid_float_2d(fourier) result(grid)

    complex(R4_KIND), intent(in), dimension(:, :)     :: fourier
    real(R4_KIND), dimension(leng, size(fourier, 2)) :: grid

    real, dimension(leng2, size(fourier, 2)) :: data
    real, dimension(leng1, size(fourier, 2)) :: work
    integer :: j, k, num

    if (.not. module_is_initialized) &
      call error_handler('fft_fourier_to_grid', 'fft_init must be called.')

    num = size(fourier, 2)    ! number of transforms

    if (size(fourier, 1) < lenc) &
      call error_handler('fft_fourier_to_grid', 'length of input data too small.')

    do j = 1, num
    do k = 1, lenc
      data(2*k - 1, j) = real(fourier(k, j))
      data(2*k, j) = aimag(fourier(k, j))
    end do
    end do
    call fft991(data, work, table99, ifax, 1, leng2, leng, num, +1)
    do j = 1, num
      grid(1:leng, j) = data(1:leng, j)
    end do

  end function fft_fourier_to_grid_float_2d

!#######################################################################

  !> Transforms 2-D 64-bit gridpoint data to Fourier space (see `fft_grid_to_fourier`).
  function fft_grid_to_fourier_double_2d(grid) result(fourier)

    real(R8_KIND), intent(in), dimension(:, :)  :: grid
    complex(R8_KIND), dimension(lenc, size(grid, 2)) :: fourier

    real, dimension(leng2, size(grid, 2)) :: data
    real, dimension(leng1, size(grid, 2)) :: work
    integer :: j, k, num

    if (.not. module_is_initialized) &
      call error_handler('fft_grid_to_fourier', 'fft_init must be called.')

    if (size(grid, 1) < leng) &
      call error_handler('fft_grid_to_fourier', 'length of input data too small.')

    num = size(grid, 2)    ! number of transforms

    do j = 1, num
      data(1:leng, j) = grid(1:leng, j)
    end do
    call fft991(data, work, table99, ifax, 1, leng2, leng, num, -1)
    do j = 1, size(grid, 2)
    do k = 1, lenc
      fourier(k, j) = cmplx(data(2*k - 1, j), data(2*k, j))
    end do
    end do

  end function fft_grid_to_fourier_double_2d

!#######################################################################

  !> Transforms 2-D 64-bit Fourier data to grid space (see `fft_fourier_to_grid`).
  function fft_fourier_to_grid_double_2d(fourier) result(grid)

    complex(R8_KIND), intent(in), dimension(:, :)     :: fourier
    real(R8_KIND), dimension(leng, size(fourier, 2)) :: grid

    real, dimension(leng2, size(fourier, 2)) :: data
    real, dimension(leng1, size(fourier, 2)) :: work
    integer :: j, k, num

    if (.not. module_is_initialized) &
      call error_handler('fft_fourier_to_grid', 'fft_init must be called.')

    num = size(fourier, 2)    ! number of transforms

    if (size(fourier, 1) < lenc) &
      call error_handler('fft_fourier_to_grid', 'length of input data too small.')

    do j = 1, num
    do k = 1, lenc
      data(2*k - 1, j) = real(fourier(k, j))
      data(2*k, j) = aimag(fourier(k, j))
    end do
    end do
    call fft991(data, work, table99, ifax, 1, leng2, leng, num, +1)
    do j = 1, num
      grid(1:leng, j) = data(1:leng, j)
    end do

  end function fft_fourier_to_grid_double_2d

!#######################################################################
!                   interface overloads
!#######################################################################

  !> Transforms 3-D 32-bit gridpoint data to Fourier space (see `fft_grid_to_fourier`).
  function fft_grid_to_fourier_float_3d(grid) result(fourier)

    real(R4_KIND), intent(in), dimension(:, :, :) :: grid
    complex(R4_KIND), dimension(lenc, size(grid, 2), size(grid, 3)) :: fourier
    integer :: n

!$omp parallel do if(size(grid, 3) > 1)
    do n = 1, size(grid, 3)
      fourier(:, :, n) = fft_grid_to_fourier_float_2d(grid(:, :, n))
    end do
!$omp end parallel do

  end function fft_grid_to_fourier_float_3d

!#######################################################################

  !> Transforms 3-D 32-bit Fourier data to grid space (see `fft_fourier_to_grid`).
  function fft_fourier_to_grid_float_3d(fourier) result(grid)

    complex(R4_KIND), intent(in), dimension(:, :, :) :: fourier
    real(R4_KIND), dimension(leng, size(fourier, 2), size(fourier, 3)) :: grid
    integer :: n

!$omp parallel do if(size(fourier, 3) > 1)
    do n = 1, size(fourier, 3)
      grid(:, :, n) = fft_fourier_to_grid_float_2d(fourier(:, :, n))
    end do
!$omp end parallel do

  end function fft_fourier_to_grid_float_3d

!#######################################################################

  !> Transforms 3-D 64-bit gridpoint data to Fourier space (see `fft_grid_to_fourier`).
  function fft_grid_to_fourier_double_3d(grid) result(fourier)

    real(R8_KIND), intent(in), dimension(:, :, :) :: grid
    complex(R8_KIND), dimension(lenc, size(grid, 2), size(grid, 3)) :: fourier
    integer :: n

!$omp parallel do if(size(grid, 3) > 1)
    do n = 1, size(grid, 3)
      fourier(:, :, n) = fft_grid_to_fourier_double_2d(grid(:, :, n))
    end do
!$omp end parallel do

  end function fft_grid_to_fourier_double_3d

!#######################################################################

  !> Transforms 3-D 64-bit Fourier data to grid space (see `fft_fourier_to_grid`).
  function fft_fourier_to_grid_double_3d(fourier) result(grid)

    complex(R8_KIND), intent(in), dimension(:, :, :) :: fourier
    real(R8_KIND), dimension(leng, size(fourier, 2), size(fourier, 3)) :: grid
    integer :: n

!$omp parallel do if(size(fourier, 3) > 1)
    do n = 1, size(fourier, 3)
      grid(:, :, n) = fft_fourier_to_grid_double_2d(fourier(:, :, n))
    end do
!$omp end parallel do

  end function fft_fourier_to_grid_double_3d

!#######################################################################

  !> Sets the length of the transforms and sets up the trigonometric tables.
  !>
  !> It must be called once before the transforms. To change the length, call `fft_end`
  !> first: calling `fft_init` again without `fft_end` stops the model.
  subroutine fft_init(n)

    integer, intent(in) :: n  !! number of real values in a single sequence; the transformed data have
                              !! n/2+1 complex values

    if (module_is_initialized) &
      call error_handler('fft_init', 'attempted to reinitialize fft')

    leng = n; leng1 = n + 1; leng2 = n + 2; lenc = n/2 + 1

    allocate (table99(3*leng/2 + 1))
    allocate (ifax(13))
    call set99(table99, ifax, leng)

    module_is_initialized = .true.

  end subroutine fft_init

!#######################################################################

  !> Unsets the transform length and frees the tables; stops if `fft_init` has not been
  !> called.
  subroutine fft_end

    if (.not. module_is_initialized) &
      call error_handler('fft_end', 'attempt to un-initialize fft that has not been initialized')

    leng = 0; leng1 = 0; leng2 = 0; lenc = 0

    if (allocated(table99)) deallocate (table99)
    if (allocated(ifax)) deallocate (ifax)

    module_is_initialized = .false.

  end subroutine fft_end

!#######################################################################

  subroutine error_handler(routine, message)
    character(len=*), intent(in) :: routine, message

    call error_mesg(routine, message, FATAL)

  end subroutine error_handler

!#######################################################################

end module fft_mod
