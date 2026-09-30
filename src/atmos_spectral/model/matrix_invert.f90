!> Inversion of a small nonsingular matrix, with its determinant.
!>
!> Original authors: Triveni N. Upadhyay.
module matrix_invert_mod

  use fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL, &
                     write_version_number

  implicit none

  public  :: invert
  integer, private :: maxmag

  character(len=128), parameter :: version = '$Id matrix_invert.f90 $'  !! version string
  character(len=128), parameter :: tagname = '$Name: lima $'  !! tag name
  logical :: entry_to_logfile_done = .false.  !! whether the version has been written to the log file

contains

  !> Inverts an n by n nonsingular matrix in place and returns its determinant.
  !>
  !> Uses elementary transformations with the pivotal-element method (column interchanges
  !> to put the largest element of the row on the diagonal). If the magnitude of the
  !> determinant falls below 1e-30 the matrix is taken to be singular and the model stops.
  subroutine invert(matrix, det)

    real, intent(inout), dimension(:, :) :: matrix  !! the matrix on input, its inverse on output
    real, intent(out) :: det  !! determinant of the input matrix

    real, dimension(2*size(matrix, 1)) :: dd, h
    real, dimension(2*size(matrix, 1), size(matrix, 1)) :: ac, temp
    real :: min_det = 1.0e-30
    character(len=24) :: chtmp
    integer :: n, i, j, L, m, k

    if (.not. entry_to_logfile_done) then
      call write_version_number(version, tagname)
      entry_to_logfile_done = .true.
    end if

    n = size(matrix, 1)

!   INITIALIZE

    det = 1.0
    ac(1:n, :) = matrix(:, :)
    do j = 1, n
      ac(n + 1:2*n, j) = 0.0
      ac(n + j, j) = 1.0
    end do

    do k = 1, n

! FIND  LARGEST ELEMENT IN THE ROW

      h(k:n) = ac(k, k:n)
      m = n - k + 1
      L = max_mag(h(k:n), M) + k

! INTERCHANGE COLUMNS IF THE LARGEST ELEMENT IS NOT THE DIAGONAL ELEMENT.

      if (k - L < 0) then
        do i = k, 2*n
          dd(i) = ac(i, k)
          ac(i, k) = ac(i, L)
          ac(i, L) = dd(i)
        end do
        det = -det
      end if

! DIVIDE THE COLUMN BY THE LARGEST ELEMENT

      det = det*ac(k, k)
      if (abs(det) < min_det) then
        write (chtmp, '(1pe24.16)') det
        call error_mesg('invert', 'DETERMINANT OF MATRIX ='//chtmp// &
        & ' THE MAGNITUDE OF THE DETERMINANT IS LESS THAN THE MINIMUM ALLOWED. &
        &  THE INPUT MATRIX APPEARS TO BE SINGULAR.', FATAL)
      end if
      h(k:2*n) = ac(k:2*n, k)/ac(k, k)
      do j = 1, n
        temp(k:2*n, j) = h(k:2*n)*ac(k, j)
      end do
      ac(k:2*n, :) = ac(k:2*n, :) - temp(k:2*n, :)
      ac(k:2*n, k) = h(k:2*n)
    end do

    matrix(1:n, :) = ac(n + 1:2*n, :)

    return
  end subroutine invert

  !> Returns the offset from the start of `h` (0 for the first element) of the element with
  !> the largest magnitude.
  function max_mag(h, m) result(max)

    integer, intent(in) :: m  !! length of `h`
    real, intent(in) :: h(m)  !! values to search
    integer :: max, i
    real :: rmax

    max = 0
    rmax = abs(h(1))
    do i = 1, m
      if (abs(h(i)) > rmax) then
        rmax = abs(h(i))
        max = i - 1
      end if
    end do
    return
  end function max_mag

end module matrix_invert_mod
