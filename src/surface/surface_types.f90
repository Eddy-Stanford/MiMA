!> Derived types shared by the surface components: the grid and the per-step exchange fields.
!>
!> `surface_grid_type` describes this PE's part of the model grid, on which the ocean, the
!> land and the flux exchange all work (MiMA has one grid, so there is no exchange grid).
!> `surface_exchange_type` holds the fields that `flux_down_from_atmos` hands to the
!> ocean and land models at each step. The module also holds the scan of `input.nml` with
!> which the components stop with a migration message when a namelist group or variable
!> that was removed in v2 is still set (unknown groups are otherwise ignored silently).
module surface_types_mod

  use mpp_domains_mod, only: domain2d
  use time_manager_mod, only: time_type
  use fms_mod, only: input_nml_file, lowercase

  implicit none
  private

  public :: surface_grid_type, surface_exchange_type, nml_has_group, nml_sets

  !> This PE's part of the model grid, as seen by the surface components.
  type surface_grid_type
    type(domain2d) :: domain  !! domain decomposition (that of the atmosphere)
    integer :: axes(2)  !! diagnostic axis ids (longitude, latitude)
    real, pointer :: lon_bnd(:) => null()  !! longitudes of the grid-box boundaries on this PE [rad]
    real, pointer :: lat_bnd(:) => null()  !! latitudes of the grid-box boundaries on this PE [rad]
    real, allocatable :: area(:, :)
    !! relative grid-box area `dlon*(sin(lat_n) - sin(lat_s))` (the Earth's radius squared omitted)
    integer :: nx  !! number of longitudes on this PE
    integer :: ny  !! number of latitudes on this PE
  end type surface_grid_type

  !> Fields handed from `flux_down_from_atmos` to the ocean and land models at each step.
  !>
  !> All arrays are on this PE's grid (`nx`, `ny`). The fluxes and derivatives are the
  !> values after the implicit coupling with the atmosphere's vertical diffusion has been
  !> completed for the lowest level, at the surface temperature of the start of the step.
  type surface_exchange_type
    real, allocatable :: flux_lw(:, :)  !! net downward longwave flux at the surface [W/m2]
    real, allocatable :: flux_sw(:, :)  !! net downward shortwave flux at the surface [W/m2]
    real, allocatable :: lprec(:, :)  !! liquid precipitation rate [kg/m2/s]
    real, allocatable :: fprec(:, :)  !! frozen precipitation rate [kg/m2/s]
    real, allocatable :: flux_t(:, :)  !! sensible heat flux (upward) [W/m2]
    real, allocatable :: flux_q(:, :)  !! evaporation rate (upward) [kg/m2/s]
    real, allocatable :: dhdt_surf(:, :)
    !! derivative of `flux_t` with respect to the surface temperature [W/m2/K]
    real, allocatable :: dedt_surf(:, :)
    !! derivative of `flux_q` with respect to the surface temperature [kg/m2/s/K]
    real, allocatable :: drdt_surf(:, :)
    !! derivative of the upward longwave flux with respect to the surface temperature [W/m2/K]
    real, allocatable :: extra_heating(:, :)  !! local surface heating (`local_heating_nml`; 0 if off) [K/s]
    real, allocatable :: qflux_land(:, :)
    !! Q-flux applied over land: 0 unless `qflux_ocean_only = .false.` (`qflux_nml`), when it
    !! is the Q-flux on the land points, as in v1; set at initialization and, for a Q-flux
    !! read from a file, at every step by `update_ocean_model_fast` [W/m2]
    real :: dt  !! time step [s]
    type(time_type) :: Time  !! time at the end of the step
  end type surface_exchange_type

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

end module surface_types_mod
