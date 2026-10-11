!> Derived types shared by the surface components: the grid and the per-step exchange fields.
!>
!> `surface_grid_type` describes this PE's part of the model grid, on which the ocean, the
!> land and the flux exchange all work (MiMA has one grid, so there is no exchange grid).
!> `surface_exchange_type` holds the fields that `flux_down_from_atmos` hands to the
!> ocean and land models at each step. The module holds no logic.
module surface_types_mod

  use mpp_domains_mod, only: domain2d
  use time_manager_mod, only: time_type

  implicit none
  private

  public :: surface_grid_type, surface_exchange_type

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
    real, allocatable :: qflux_legacy(:, :)
    !! prescribed ocean heat flux added at all points, ocean and land, as in v1 [W/m2]
    real :: dt  !! time step [s]
    type(time_type) :: Time  !! time at the end of the step
  end type surface_exchange_type

end module surface_types_mod
