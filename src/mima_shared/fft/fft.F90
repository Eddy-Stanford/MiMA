
!> Fast Fourier transforms between real grid space and complex Fourier space, for many
!> sequences at once.
!>
!> Computes multiple 1-D FFTs and inverse FFTs of 2-D and 3-D arrays, between real
!> gridpoint values and complex Fourier coefficients, in single (32-bit) and double
!> (64-bit) precision. `fft_init` sets the length `n` of the transforms. The complex
!> Fourier components are stored as
!>
!> ```text
!> fourier(1)     = cmplx(a(0), b(0))
!> fourier(2)     = cmplx(a(1), b(1))
!>     ...
!> fourier(n/2+1) = cmplx(a(n/2), b(n/2))
!> ```
!>
!> By default the stand-alone Temperton FFT of `fft99_mod` is used, at the real
!> precision chosen at compile time (with 32-bit reals the transforms cannot be done at
!> 64-bit precision). Compiled with `-D NAGFFT`, the NAG library routines C06FPF, C06FQF
!> and C06GQF are used instead (64-bit data only); on Cray and SGI systems the vendor
!> scientific library routines SCFFTM, CSFFTM, DZFFTM and ZDFFTM are used. Compiled with
!> `-D test_fft`, the file also contains a test program that transforms random data to
!> Fourier space and back and prints it with the original.
!>
!> Original authors: Bruce Wyman.
module fft_mod

!-----------------------------------------------------------------------
!these are used to determine hardware/OS/compiler

#ifdef __sgi
#  ifdef _COMPILER_VERSION
!the MIPSPro compiler defines _COMPILER_VERSION
#    define sgi_mipspro
#  else
#    define sgi_generic
#  endif
#endif

!fft uses the SCILIB on SGICRAY, and the NAG library otherwise
#if defined(_CRAY) || defined(sgi_mipspro)
#  define SGICRAY
#endif

  use platform_mod, only: R8_KIND, R4_KIND
  use fms_mod, only: write_version_number, &
                     error_mesg, FATAL
#if !defined(SGICRAY) && !defined(NAGFFT)
  use fft99_mod, only: fft991, set99
#endif

  implicit none
  private

!----------------- interfaces --------------------

  public :: fft_init, fft_end, fft_grid_to_fourier, fft_fourier_to_grid

  !> Transforms multiple sequences of real gridpoint values to complex Fourier
  !> coefficients.
  !>
  !> `fourier = fft_grid_to_fourier(grid)`, for 2-D or 3-D arrays of 32- or 64-bit reals.
  !> Each sequence along the first dimension of `grid` has the length `n` set by
  !> `fft_init`; that dimension must be at least `n` (exactly n+1 with the Cray/SGI
  !> library). The first dimension of `fourier` is n/2+1 and the other dimensions are
  !> those of `grid`. Stops if `fft_init` has not been called, or for 32-bit data with the
  !> NAG library.
  interface fft_grid_to_fourier
    module procedure fft_grid_to_fourier_float_2d, fft_grid_to_fourier_double_2d, &
      fft_grid_to_fourier_float_3d, fft_grid_to_fourier_double_3d
  end interface

  !> Transforms multiple sequences of complex Fourier coefficients to real gridpoint
  !> values (the inverse of `fft_grid_to_fourier`).
  !>
  !> `grid = fft_fourier_to_grid(fourier)`, for 2-D or 3-D arrays of 32- or 64-bit
  !> complex values. The first dimension of `fourier` must be at least n/2+1 (exactly
  !> n/2+1 with the Cray/SGI library), where `n` is the length set by `fft_init`. The first
  !> dimension of `grid` is n+1, with the values in the first `n`; the other dimensions
  !> are those of `fourier`. Stops if `fft_init` has not been called, or for 32-bit data
  !> with the NAG library.
  interface fft_fourier_to_grid
    module procedure fft_fourier_to_grid_float_2d, fft_fourier_to_grid_double_2d, &
      fft_fourier_to_grid_float_3d, fft_fourier_to_grid_double_3d
  end interface

!---------------------- private data -----------------------------------

! tables for trigonometric constants and factors
! (not all will be used)
  real(R8_KIND), allocatable, dimension(:) :: table8
  real(R4_KIND), allocatable, dimension(:) :: table4
  real, allocatable, dimension(:) :: table99
  integer, allocatable, dimension(:) :: ifax

  logical :: do_log = .true.
  integer :: leng, leng1, leng2, lenc    ! related to transform size

  logical :: module_is_initialized = .false.

!  cvs version and tag name
  character(len=128) :: version = '$Id: fft.F90,v 10.0 2003/10/24 22:01:29 fms Exp $'
  character(len=128) :: tagname = '$Name: lima $'

contains

!#######################################################################

  !> Transforms 2-D 32-bit gridpoint data to Fourier space (see `fft_grid_to_fourier`).
  function fft_grid_to_fourier_float_2d(grid) result(fourier)

!-----------------------------------------------------------------------

    real(R4_KIND), intent(in), dimension(:, :)  :: grid
    complex(R4_KIND), dimension(lenc, size(grid, 2)) :: fourier

!-----------------------------------------------------------------------
!
!  input
!  -----
!   grid = Multiple transforms in grid point space, the first dimension
!          must be n+1 (where n is the size of each real transform).
!
!  returns
!  -------
!    Multiple transforms in complex fourier space, the first dimension
!    must equal n/2+1 (where n is the size of each real transform).
!    The remaining dimensions must be the same size as the input
!    argument "grid".
!
!-----------------------------------------------------------------------
#ifdef SGICRAY
#  ifdef _CRAY
!  local storage for cray fft
    real(R4_KIND), dimension((2*leng + 4)*size(grid, 2)) :: work
#  else
!  local storage for sgi fft
    real(R4_KIND), dimension(leng2) :: work
#  endif
#else
#  ifdef NAGFFT
!  local storage for nag fft
    real(R4_KIND), dimension(size(grid, 2), leng) :: data, work
#  else
!  local storage for temperton fft
    real, dimension(leng2, size(grid, 2)) :: data
    real, dimension(leng1, size(grid, 2)) :: work
#  endif
#endif

    real(R4_KIND) :: scale
    integer :: j, k, num, len_grid, ifail

!-----------------------------------------------------------------------

    if (.not. module_is_initialized) &
      call error_handler('fft_grid_to_fourier', &
                         'fft_init must be called.')

!-----------------------------------------------------------------------

    len_grid = size(grid, 1)
#ifdef SGICRAY
    if (len_grid /= leng1) call error_handler('fft_grid_to_fourier', &
                                              'size of first dimension of input data is wrong')
#else
    if (len_grid < leng) call error_handler('fft_grid_to_fourier', &
                                            'length of input data too small.')
#endif
!-----------------------------------------------------------------------
!----------------transform to fourier coefficients (+1)-----------------

    num = size(grid, 2)    ! number of transforms

#ifdef SGICRAY
!  Cray/SGI fft
    scale = 1./real(leng)
#  ifdef _CRAY
    call scfftm(-1, leng, num, scale, grid, leng1, fourier, lenc, &
                table4, work, 0)
#  else
    call scfftm(-1, leng, num, scale, grid, leng1, fourier, lenc, &
                table4, work, 0)
#  endif
#else
#  ifdef NAGFFT
!  NAG fft
!  will not allow float kind for NAG
    call error_handler('fft_grid_to_fourier', &
                       'float kind not supported for nag fft')
    do j = 1, size(grid, 2)
      data(j, 1:leng) = grid(1:leng, j)
    end do
! call c06fpe ( num, leng, data, 's', table4, work, ifail )
    scale = 1./sqrt(float(leng))
    data = data*scale
    fourier(1, :) = cmplx(data(:, 1), 0.)
    do k = 2, lenc - 1
      fourier(k, :) = cmplx(data(:, k), data(:, leng - k + 2))
    end do
    fourier(lenc, :) = cmplx(data(:, lenc), 0.)
#  else
!  Temperton fft
    do j = 1, num
      data(1:leng, j) = grid(1:leng, j)
    end do
    call fft991(data, work, table99, ifax, 1, leng2, leng, num, -1)
    do j = 1, size(grid, 2)
    do k = 1, lenc
      fourier(k, j) = cmplx(data(2*k - 1, j), data(2*k, j))
    end do
    end do
#  endif
#endif
!-----------------------------------------------------------------------

  end function fft_grid_to_fourier_float_2d

!#######################################################################

  !> Transforms 2-D 32-bit Fourier data to grid space (see `fft_fourier_to_grid`).
  function fft_fourier_to_grid_float_2d(fourier) result(grid)

!-----------------------------------------------------------------------

    complex(R4_KIND), intent(in), dimension(:, :)     :: fourier
    real(R4_KIND), dimension(leng1, size(fourier, 2)) :: grid

!-----------------------------------------------------------------------
!
!  input
!  -----
!  fourier = Multiple transforms in complex fourier space, the first
!            dimension must equal n/2+1 (where n is the size of each
!            real transform).
!
!  returns
!  -------
!    Multiple transforms in grid point space, the first dimension
!    must be n+1 (where n is the size of each real transform).
!    The remaining dimensions must be the same size as the input
!    argument "fourier".
!
!-----------------------------------------------------------------------
#ifdef SGICRAY
#  ifdef _CRAY
!  local storage for cray fft
    real(R4_KIND), dimension((2*leng + 4)*size(fourier, 2)) :: work
#  else
!  local storage for sgi fft
    real(R4_KIND), dimension(leng2) :: work
#  endif
#else
#  ifdef NAGFFT
!  local storage for nag fft
    real(R4_KIND), dimension(size(fourier, 2), leng) :: data, work
#  else
!  local storage for temperton fft
    real, dimension(leng2, size(fourier, 2)) :: data
    real, dimension(leng1, size(fourier, 2)) :: work
#  endif
#endif

    real(R4_KIND) :: scale
    integer :: j, k, num, len_fourier, ifail

!-----------------------------------------------------------------------

    if (.not. module_is_initialized) &
      call error_handler('fft_grid_to_fourier', &
                         'fft_init must be called.')

!-----------------------------------------------------------------------

    len_fourier = size(fourier, 1)
    num = size(fourier, 2)    ! number of transforms

#ifdef SGICRAY
    if (len_fourier /= lenc) call error_handler('fft_fourier_to_grid', &
                                                'size of first dimension of input data is wrong')
#else
    if (len_fourier < lenc) call error_handler('fft_fourier_to_grid', &
                                               'length of input data too small.')
#endif
!-----------------------------------------------------------------------
!----------------inverse transform to real space (-1)-------------------

#ifdef SGICRAY
!  Cray/SGI fft
    scale = 1.0
#  ifdef _CRAY
    call csfftm(+1, leng, num, scale, fourier, len_fourier, &
                grid, leng1, table4, work, 0)
#  else
    call csfftm(+1, leng, num, scale, fourier, len_fourier, &
                grid, leng1, table4, work, 0)
#  endif
#else
#  ifdef NAGFFT
!  NAG fft
!  will not allow float kind for nag
    call error_handler('fft_fourier_to_grid', &
                       'float kind not supported for nag fft')

    ! save input complex array in real format (herm.)
    do k = 1, lenc
      data(:, k) = real(fourier(k, :))
    end do
    do k = 2, lenc - 1
      data(:, leng - k + 2) = aimag(fourier(k, :))
    end do

! call c06gqe ( num, leng, data, ifail )
! call c06fqe ( num, leng, data, 's', table4, work, ifail )

    ! scale and transpose data
    scale = sqrt(real(leng))
    do j = 1, num
      grid(1:leng, j) = data(j, 1:leng)*scale
    end do
#  else
!  Temperton fft
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
#  endif
#endif

!-----------------------------------------------------------------------

  end function fft_fourier_to_grid_float_2d

!#######################################################################
  !> Transforms 2-D 64-bit gridpoint data to Fourier space (see `fft_grid_to_fourier`).
  function fft_grid_to_fourier_double_2d(grid) result(fourier)

!-----------------------------------------------------------------------

    real(R8_KIND), intent(in), dimension(:, :)  :: grid
    complex(R8_KIND), dimension(lenc, size(grid, 2)) :: fourier

!-----------------------------------------------------------------------
!
!  input
!  -----
!   grid = Multiple transforms in grid point space, the first dimension
!          must be n+1 (where n is the size of each real transform).
!
!  returns
!  -------
!    Multiple transforms in complex fourier space, the first dimension
!    must equal n/2+1 (where n is the size of each real transform).
!    The remaining dimensions must be the same size as the input
!    argument "grid".
!
!-----------------------------------------------------------------------
#ifdef SGICRAY
#  ifdef _CRAY
!  local storage for cray fft
    real(R8_KIND), dimension((2*leng + 4)*size(grid, 2)) :: work
#  else
!  local storage for sgi fft
    real(R8_KIND), dimension(leng2) :: work
#  endif
#else
#  ifdef NAGFFT
!  local storage for nag fft
    real(R8_KIND), dimension(size(grid, 2), leng) :: data, work
#  else
!  local storage for temperton fft
    real, dimension(leng2, size(grid, 2)) :: data
    real, dimension(leng1, size(grid, 2)) :: work
#  endif
#endif

    real(R8_KIND) :: scale
    integer :: j, k, num, len_grid, ifail

!-----------------------------------------------------------------------

    if (.not. module_is_initialized) &
      call error_handler('fft_grid_to_fourier', &
                         'fft_init must be called.')

!-----------------------------------------------------------------------

    len_grid = size(grid, 1)
#ifdef SGICRAY
    if (len_grid /= leng1) call error_handler('fft_grid_to_fourier', &
                                              'size of first dimension of input data is wrong')
#else
    if (len_grid < leng) call error_handler('fft_grid_to_fourier', &
                                            'length of input data too small.')
#endif
!-----------------------------------------------------------------------
!----------------transform to fourier coefficients (+1)-----------------

    num = size(grid, 2)    ! number of transforms
#ifdef SGICRAY
!  Cray/SGI fft
    scale = 1./float(leng)
#  ifdef _CRAY
    call scfftm(-1, leng, num, scale, grid, leng1, fourier, lenc, &
                table8, work, 0)
#  else
    call dzfftm(-1, leng, num, scale, grid, leng1, fourier, lenc, &
                table8, work, 0)
#  endif
#else
#  ifdef NAGFFT
!  NAG fft
    do j = 1, size(grid, 2)
      data(j, 1:leng) = grid(1:leng, j)
    end do
    call c06fpf(num, leng, data, 's', table8, work, ifail)
    scale = 1./sqrt(float(leng))
    data = data*scale
    fourier(1, :) = cmplx(data(:, 1), 0.)
    do k = 2, lenc - 1
      fourier(k, :) = cmplx(data(:, k), data(:, leng - k + 2))
    end do
    fourier(lenc, :) = cmplx(data(:, lenc), 0.)
#  else
!  Temperton fft
    do j = 1, num
      data(1:leng, j) = grid(1:leng, j)
    end do
    call fft991(data, work, table99, ifax, 1, leng2, leng, num, -1)
    do j = 1, size(grid, 2)
    do k = 1, lenc
      fourier(k, j) = cmplx(data(2*k - 1, j), data(2*k, j))
    end do
    end do
#  endif
#endif
!-----------------------------------------------------------------------

  end function fft_grid_to_fourier_double_2d

!#######################################################################

  !> Transforms 2-D 64-bit Fourier data to grid space (see `fft_fourier_to_grid`).
  function fft_fourier_to_grid_double_2d(fourier) result(grid)

!-----------------------------------------------------------------------

    complex(R8_KIND), intent(in), dimension(:, :)     :: fourier
    real(R8_KIND), dimension(leng1, size(fourier, 2)) :: grid

!-----------------------------------------------------------------------
!
!  input
!  -----
!  fourier = Multiple transforms in complex fourier space, the first
!            dimension must equal n/2+1 (where n is the size of each
!            real transform).
!
!  returns
!  -------
!    Multiple transforms in grid point space, the first dimension
!    must be n+1 (where n is the size of each real transform).
!    The remaining dimensions must be the same size as the input
!    argument "fourier".
!
!-----------------------------------------------------------------------
#ifdef SGICRAY
#  ifdef _CRAY
!  local storage for cray fft
    real(R8_KIND), dimension((2*leng + 4)*size(fourier, 2)) :: work
#  else
!  local storage for sgi fft
    real(R8_KIND), dimension(leng2) :: work
#  endif
#else
#  ifdef NAGFFT
!  local storage for nag fft
    real(R8_KIND), dimension(size(fourier, 2), leng) :: data, work
#  else
!  local storage for temperton fft
    real, dimension(leng2, size(fourier, 2)) :: data
    real, dimension(leng1, size(fourier, 2)) :: work
#  endif
#endif

    real(R8_KIND) :: scale
    integer :: j, k, num, len_fourier, ifail

!-----------------------------------------------------------------------

    if (.not. module_is_initialized) &
      call error_handler('fft_grid_to_fourier', &
                         'fft_init must be called.')

!-----------------------------------------------------------------------

    len_fourier = size(fourier, 1)
    num = size(fourier, 2)    ! number of transforms

#ifdef SGICRAY
    if (len_fourier /= lenc) call error_handler('fft_fourier_to_grid', &
                                                'size of first dimension of input data is wrong')
#else
    if (len_fourier < lenc) call error_handler('fft_fourier_to_grid', &
                                               'length of input data too small.')
#endif
!-----------------------------------------------------------------------
!----------------inverse transform to real space (-1)-------------------

#ifdef SGICRAY
!  Cray/SGI fft
    scale = 1.0
#  ifdef _CRAY
    call csfftm(+1, leng, num, scale, fourier, len_fourier, &
                grid, leng1, table8, work, 0)
#  else
    call zdfftm(+1, leng, num, scale, fourier, len_fourier, &
                grid, leng1, table8, work, 0)
#  endif
#else
#  ifdef NAGFFT
!  NAG fft

    ! save input complex array in real format (herm.)
    do k = 1, lenc
      data(:, k) = real(fourier(k, :))
    end do
    do k = 2, lenc - 1
      data(:, leng - k + 2) = aimag(fourier(k, :))
    end do

    call c06gqf(num, leng, data, ifail)
    call c06fqf(num, leng, data, 's', table8, work, ifail)

    ! scale and transpose data
    scale = sqrt(real(leng))
    do j = 1, num
      grid(1:leng, j) = data(j, 1:leng)*scale
    end do
#  else
!  Temperton fft
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
#  endif
#endif

!-----------------------------------------------------------------------

  end function fft_fourier_to_grid_double_2d

!#######################################################################
!                   interface overloads
!#######################################################################
  !> Transforms 3-D 32-bit gridpoint data to Fourier space (see `fft_grid_to_fourier`).
  function fft_grid_to_fourier_float_3d(grid) result(fourier)

!-----------------------------------------------------------------------
    real(R4_KIND), intent(in), dimension(:, :, :) :: grid
    complex(R4_KIND), dimension(lenc, size(grid, 2), size(grid, 3)) :: fourier
    integer :: n
!-----------------------------------------------------------------------

    do n = 1, size(grid, 3)
      fourier(:, :, n) = fft_grid_to_fourier_float_2d(grid(:, :, n))
    end do

!-----------------------------------------------------------------------

  end function fft_grid_to_fourier_float_3d

!#######################################################################

  !> Transforms 3-D 32-bit Fourier data to grid space (see `fft_fourier_to_grid`).
  function fft_fourier_to_grid_float_3d(fourier) result(grid)

!-----------------------------------------------------------------------
    complex(R4_KIND), intent(in), dimension(:, :, :) :: fourier
    real(R4_KIND), dimension(leng1, size(fourier, 2), size(fourier, 3)) :: grid
    integer :: n
!-----------------------------------------------------------------------

    do n = 1, size(fourier, 3)
      grid(:, :, n) = fft_fourier_to_grid_float_2d(fourier(:, :, n))
    end do

!-----------------------------------------------------------------------

  end function fft_fourier_to_grid_float_3d

!#######################################################################

  !> Transforms 3-D 64-bit gridpoint data to Fourier space (see `fft_grid_to_fourier`).
  function fft_grid_to_fourier_double_3d(grid) result(fourier)

!-----------------------------------------------------------------------
    real(R8_KIND), intent(in), dimension(:, :, :) :: grid
    complex(R8_KIND), dimension(lenc, size(grid, 2), size(grid, 3)) :: fourier
    integer :: n
!-----------------------------------------------------------------------

    do n = 1, size(grid, 3)
      fourier(:, :, n) = fft_grid_to_fourier_double_2d(grid(:, :, n))
    end do

!-----------------------------------------------------------------------

  end function fft_grid_to_fourier_double_3d

!#######################################################################

  !> Transforms 3-D 64-bit Fourier data to grid space (see `fft_fourier_to_grid`).
  function fft_fourier_to_grid_double_3d(fourier) result(grid)

!-----------------------------------------------------------------------
    complex(R8_KIND), intent(in), dimension(:, :, :) :: fourier
    real(R8_KIND), dimension(leng1, size(fourier, 2), size(fourier, 3)) :: grid
    integer :: n
!-----------------------------------------------------------------------

    do n = 1, size(fourier, 3)
      grid(:, :, n) = fft_fourier_to_grid_double_2d(fourier(:, :, n))
    end do

!-----------------------------------------------------------------------

  end function fft_fourier_to_grid_double_3d

!#######################################################################

  !> Sets the length of the transforms and sets up the trigonometric tables.
  !>
  !> It must be called once before the transforms. To change the length, call `fft_end`
  !> first: calling `fft_init` again without `fft_end` stops the model.
  subroutine fft_init(n)

!-----------------------------------------------------------------------
    integer, intent(in) :: n  !! number of real values in a single sequence; the transformed data have
                              !! n/2+1 complex values
!-----------------------------------------------------------------------
#ifdef SGICRAY
    real(R4_KIND) ::  dummy4(1)
    complex(R4_KIND) :: cdummy4(1)
    real(R8_KIND) ::  dummy8(1)
    complex(R8_KIND) :: cdummy8(1)
    integer :: isys(0:1)
#else
#  ifdef NAGFFT
    real(R8_KIND) :: data8(n), work8(n)
    real(R4_KIND) :: data4(n), work4(n)
    integer       :: ifail4, ifail8
#  endif
#endif
!-----------------------------------------------------------------------
!   --- fourier transform initialization ----

    if (module_is_initialized) &
      call error_handler('fft_init', 'attempted to reinitialize fft')

!  write version and tag name to log file
    if (do_log) then
      call write_version_number(version, tagname)
      do_log = .false.
    end if

!  variables that save length of transform
    leng = n; leng1 = n + 1; leng2 = n + 2; lenc = n/2 + 1

#ifdef SGICRAY
#  ifdef _CRAY
!  initialization for cray
!  float kind may not apply for cray
    allocate (table4(100 + 2*leng), table8(100 + 2*leng))   ! size may be too large?
    call scfftm(0, leng, 1, 0.0, dummy4, 1, cdummy4, 1, table4, dummy4, 0)
    call scfftm(0, leng, 1, 0.0, dummy8, 1, cdummy8, 1, table8, dummy8, 0)
#  else
!  initialization for sgi
    allocate (table4(leng + 256), table8(leng + 256))
    isys(0) = 1
    call scfftm(0, leng, 1, 0.0, dummy4, 1, cdummy4, 1, table4, dummy8, isys)
    call dzfftm(0, leng, 1, 0.0, dummy8, 1, cdummy8, 1, table8, dummy8, isys)
#  endif
#else
#  ifdef NAGFFT
!  initialization for nag fft
    ifail8 = 0
    allocate (table8(100 + 2*leng))   ! size may be too large?
    call c06fpf(1, leng, data8, 'i', table8, work8, ifail8)

!  will not allow float kind for nag
    ifail4 = 0
! allocate (table4(100+2*leng))
! call c06fpe ( 1, leng, data4, 'i', table4, work4, ifail4 )

    if (ifail4 /= 0 .or. ifail8 /= 0) then
      call error_handler('fft_init', 'nag fft initialization error')
    end if
#  else
!  initialization for Temperton fft
    allocate (table99(3*leng/2 + 1))
    allocate (ifax(10))
    call set99(table99, ifax, leng)
#  endif
#endif

    module_is_initialized = .true.

!-----------------------------------------------------------------------

  end subroutine fft_init

!#######################################################################
  !> Unsets the transform length and frees the tables; stops if `fft_init` has not been
  !> called.
  subroutine fft_end

!-----------------------------------------------------------------------
!
!   unsets transform size and deallocates memory
!
!-----------------------------------------------------------------------
!   --- fourier transform un-initialization ----

    if (.not. module_is_initialized) &
      call error_handler('fft_end', &
                         'attempt to un-initialize fft that has not been initialized')

    leng = 0; leng1 = 0; leng2 = 0; lenc = 0

    if (allocated(table4)) deallocate (table4)
    if (allocated(table8)) deallocate (table8)
    if (allocated(table99)) deallocate (table99)

    module_is_initialized = .false.

!-----------------------------------------------------------------------

  end subroutine fft_end

!#######################################################################
! wrapper for handling errors

  subroutine error_handler(routine, message)
    character(len=*), intent(in) :: routine, message

    call error_mesg(routine, message, FATAL)

!  print *, 'ERROR: ',trim(routine)
!  print *, 'ERROR: ',trim(message)
!  stop 111

  end subroutine error_handler

!#######################################################################

end module fft_mod

#ifdef test_fft
program test
  use fft_mod
  integer, parameter :: lot = 2
  real, allocatable :: ain(:, :), aout(:, :)
  complex, allocatable :: four(:, :)
  integer :: i, j, m, n
  integer :: ntrans(2) = (/60, 90/)

! test multiple transform lengths
  do m = 1, 2

    ! set up input data
    n = ntrans(m)
    allocate (ain(n + 1, lot), aout(n + 1, lot), four(n/2 + 1, lot))
    call random_number(ain(1:n, :))
    aout(1:n, :) = ain(1:n, :)

    call fft_init(n)
    ! transform grid to fourier and back
    four = fft_grid_to_fourier(aout)
    aout = fft_fourier_to_grid(four)

    ! print original and transformed
    do j = 1, lot
    do i = 1, n
      write (*, '(2i4,3(2x,f15.9))') j, i, ain(i, j), aout(i, j), aout(i, j) - ain(i, j)
    end do
    end do

    call fft_end
    deallocate (ain, aout, four)
  end do

end program test
#endif
