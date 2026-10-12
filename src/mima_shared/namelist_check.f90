!> Scan of `input.nml` for namelist groups and variables that were removed in v2.
!>
!> The namelist read reports a variable it does not know only as a generic error, and a
!> group that no module reads is ignored silently. The components use these functions
!> before their namelist read to stop with a migration message instead. Comments (from
!> the first `!` of a line) are ignored, so a mention of a removed name in a comment is
!> not a match.
module namelist_check_mod

  use fms_mod, only: input_nml_file, lowercase

  implicit none
  private

  public :: nml_has_group, nml_sets

contains

  !> Whether a line of `input.nml` opens the namelist group `group` (lower case, without
  !> the `&`).
  logical function nml_has_group(group)
    character(len=*), intent(in) :: group  !! name of the group, in lower case

    integer :: i

    nml_has_group = .false.
    do i = 1, size(input_nml_file)
      if (index(nml_line(i), '&'//group) == 1) then
        nml_has_group = .true.
        return
      end if
    end do

  end function nml_has_group

  !> Whether a line of `input.nml` sets the namelist variable `key` (lower case): `key` as a
  !> whole word, followed by optional blanks and `=`, before the comment of the line (from
  !> its first `!`). A mention in a comment or inside a longer name is not a match.
  logical function nml_sets(key)
    character(len=*), intent(in) :: key  !! name of the variable, in lower case

    character(len=:), allocatable :: line
    integer :: i, k, pos
    logical :: word_start

    nml_sets = .false.
    do i = 1, size(input_nml_file)
      line = nml_line(i)
      pos = 0
      do
        k = index(line(pos + 1:), key)
        if (k == 0) exit
        pos = pos + k
        word_start = pos == 1
        if (.not. word_start) word_start = verify(line(pos - 1:pos - 1), &
                                                  'abcdefghijklmnopqrstuvwxyz0123456789_%') /= 0
        if (word_start .and. index(adjustl(line(pos + len(key):)), '=') == 1) then
          nml_sets = .true.
          return
        end if
      end do
    end do

  end function nml_sets

  ! Line i of input.nml, left-adjusted, in lower case and without its comment (from the
  ! first '!'; a '!' inside a quoted string is not expected in the keys scanned for)
  function nml_line(i) result(line)
    integer, intent(in) :: i
    character(len=:), allocatable :: line

    integer :: k

    line = adjustl(lowercase(input_nml_file(i)))
    k = index(line, '!')
    if (k > 0) line = line(1:k - 1)

  end function nml_line

end module namelist_check_mod
