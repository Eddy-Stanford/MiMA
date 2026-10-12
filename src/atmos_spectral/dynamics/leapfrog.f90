!> Leapfrog time stepping with a Robert filter.
!>
!> The time levels are held in the last dimension of the field and are selected by the
!> indices `previous`, `current` and `future`. When `previous == current` (the first step)
!> the step is a forward step. The Robert filter of the current level is
!> `a(current) + robert_coeff*(a(previous) - 2*a(current) + a(future))`.
!> `leapfrog` does the step and the whole filter; `leapfrog_2level_A` does the step and the
!> part of the filter without `a(future)`, and `leapfrog_2level_B` adds that last part once
!> the new level is final.
module leapfrog_mod

  use fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL, write_version_number

!===================================================================================================
  implicit none
  private
!===================================================================================================

  !> Steps a spectral field from the previous to the future level and Robert-filters the
  !> current level.
  interface leapfrog
    module procedure leapfrog_2d_complex, leapfrog_3d_complex
  end interface

  !> Steps a field from the previous to the future level and applies the first part of the
  !> Robert filter to the current level (the part that does not involve the future level).
  interface leapfrog_2level_A
    module procedure leapfrog_2level_A_2d_complex, &
      leapfrog_2level_A_3d_complex, &
      leapfrog_2level_A_3d_real
  end interface

  !> Completes the Robert filter started by `leapfrog_2level_A`: adds `robert_coeff` times the
  !> new level to the filtered level.
  interface leapfrog_2level_B
    module procedure leapfrog_2level_B_2d_complex, &
      leapfrog_2level_B_3d_complex, &
      leapfrog_2level_B_3d_real
  end interface

  character(len=128), parameter :: version = '$Id leapfrog.f90 $'
  character(len=128), parameter :: tagname = '$Name: lima $'

  public :: leapfrog, leapfrog_2level_A, leapfrog_2level_B

  logical :: entry_to_logfile_done = .false.

contains

!================================================================================

  subroutine leapfrog_2level_A_3d_complex(a, dt_a, previous, current, future, delta_t, robert_coeff)

    complex, intent(inout), dimension(:, :, :, :) :: a
    complex, intent(in), dimension(:, :, :) :: dt_a
    integer, intent(in) :: previous, current, future
    real, intent(in) :: delta_t, robert_coeff

    if (.not. entry_to_logfile_done) then
      call write_version_number(version, tagname)
      entry_to_logfile_done = .true.
    end if

    if (previous == current) then
      a(:, :, :, future) = a(:, :, :, previous) + delta_t*dt_a
      a(:, :, :, current) = a(:, :, :, current) + robert_coeff*(a(:, :, :, previous) - 2.0*a(:, :, :, current))
    else
      a(:, :, :, current) = a(:, :, :, current) + robert_coeff*(a(:, :, :, previous) - 2.0*a(:, :, :, current))
      a(:, :, :, future) = a(:, :, :, previous) + delta_t*dt_a
    end if

    return
  end subroutine leapfrog_2level_A_3d_complex

!================================================================================

  subroutine leapfrog_2level_B_3d_complex(a, current, future, robert_coeff)

    complex, intent(inout), dimension(:, :, :, :) :: a
    integer, intent(in) :: current, future
    real, intent(in) :: robert_coeff

    if (.not. entry_to_logfile_done) then
      call write_version_number(version, tagname)
      entry_to_logfile_done = .true.
    end if

    a(:, :, :, current) = a(:, :, :, current) + robert_coeff*a(:, :, :, future)

    return
  end subroutine leapfrog_2level_B_3d_complex

!================================================================================

  subroutine leapfrog_2level_A_3d_real(a, dt_a, previous, current, future, delta_t, robert_coeff)

    real, intent(inout), dimension(:, :, :, :) :: a
    real, intent(in), dimension(:, :, :) :: dt_a
    integer, intent(in) :: previous, current, future
    real, intent(in) :: delta_t, robert_coeff

    if (.not. entry_to_logfile_done) then
      call write_version_number(version, tagname)
      entry_to_logfile_done = .true.
    end if

    if (previous == current) then
      a(:, :, :, future) = a(:, :, :, previous) + delta_t*dt_a
      a(:, :, :, current) = a(:, :, :, current) + robert_coeff*(a(:, :, :, previous) - 2.0*a(:, :, :, current))
    else
      a(:, :, :, current) = a(:, :, :, current) + robert_coeff*(a(:, :, :, previous) - 2.0*a(:, :, :, current))
      a(:, :, :, future) = a(:, :, :, previous) + delta_t*dt_a
    end if

    return
  end subroutine leapfrog_2level_A_3d_real

!================================================================================

  subroutine leapfrog_2level_B_3d_real(a, current, future, robert_coeff)

    real, intent(inout), dimension(:, :, :, :) :: a
    integer, intent(in) :: current, future
    real, intent(in) :: robert_coeff

    if (.not. entry_to_logfile_done) then
      call write_version_number(version, tagname)
      entry_to_logfile_done = .true.
    end if

    a(:, :, :, current) = a(:, :, :, current) + robert_coeff*a(:, :, :, future)

    return
  end subroutine leapfrog_2level_B_3d_real

!================================================================================

  subroutine leapfrog_2level_A_2d_complex(a, dt_a, previous, current, future, delta_t, robert_coeff)

    complex, intent(inout), dimension(:, :, :) :: a
    complex, intent(in), dimension(:, :) :: dt_a
    integer, intent(in) :: previous, current, future
    real, intent(in) :: delta_t, robert_coeff

    complex, dimension(size(a, 1), size(a, 2), 1, size(a, 3)) :: a_3d
    complex, dimension(size(a, 1), size(a, 2), 1)           :: dt_a_3d

    if (.not. entry_to_logfile_done) then
      call write_version_number(version, tagname)
      entry_to_logfile_done = .true.
    end if

    a_3d(:, :, 1, :) = a
    dt_a_3d(:, :, 1) = dt_a
    call leapfrog_2level_A_3d_complex(a_3d, dt_a_3d, previous, current, future, delta_t, robert_coeff)
    a = a_3d(:, :, 1, :)

  end subroutine leapfrog_2level_A_2d_complex
!================================================================================

  subroutine leapfrog_2level_B_2d_complex(a, current, future, robert_coeff)

    complex, intent(inout), dimension(:, :, :) :: a
    integer, intent(in) :: current, future
    real, intent(in) :: robert_coeff
    complex, dimension(size(a, 1), size(a, 2), 1, size(a, 3)) :: a_3d

    if (.not. entry_to_logfile_done) then
      call write_version_number(version, tagname)
      entry_to_logfile_done = .true.
    end if

    a_3d(:, :, 1, :) = a
    call leapfrog_2level_B_3d_complex(a_3d, current, future, robert_coeff)
    a = a_3d(:, :, 1, :)

  end subroutine leapfrog_2level_B_2d_complex
!================================================================================

  subroutine leapfrog_3d_complex(a, dt_a, previous, current, future, delta_t, robert_coeff)

    complex, intent(inout), dimension(:, :, :, :) :: a
    complex, intent(in), dimension(:, :, :) :: dt_a
    integer, intent(in) :: previous, current, future
    real, intent(in) :: delta_t, robert_coeff

    if (.not. entry_to_logfile_done) then
      call write_version_number(version, tagname)
      entry_to_logfile_done = .true.
    end if

    if (previous == current) then
      a(:, :, :, future) = a(:, :, :, previous) + delta_t*dt_a
      a(:, :, :, current) = a(:, :, :, current) + robert_coeff*(a(:, :, :, previous) - 2.0*a(:, :, :, current) + a(:, :, :, future))
    else
      a(:, :, :, current) = a(:, :, :, current) + robert_coeff*(a(:, :, :, previous) - 2.0*a(:, :, :, current))
      a(:, :, :, future) = a(:, :, :, previous) + delta_t*dt_a
      a(:, :, :, current) = a(:, :, :, current) + robert_coeff*a(:, :, :, future)
    end if

    return
  end subroutine leapfrog_3d_complex

!================================================================================

  subroutine leapfrog_2d_complex(a, dt_a, previous, current, future, delta_t, robert_coeff)

    complex, intent(inout), dimension(:, :, :) :: a
    complex, intent(in), dimension(:, :) :: dt_a
    integer, intent(in) :: previous, current, future
    real, intent(in) :: delta_t, robert_coeff

    complex, dimension(size(a, 1), size(a, 2), 1, size(a, 3)) :: a_3d
    complex, dimension(size(a, 1), size(a, 2), 1)           :: dt_a_3d

    if (.not. entry_to_logfile_done) then
      call write_version_number(version, tagname)
      entry_to_logfile_done = .true.
    end if

    a_3d(:, :, 1, :) = a
    dt_a_3d(:, :, 1) = dt_a
    call leapfrog_3d_complex(a_3d, dt_a_3d, previous, current, future, delta_t, robert_coeff)
    a = a_3d(:, :, 1, :)

    return
  end subroutine leapfrog_2d_complex

!================================================================================

end module leapfrog_mod
