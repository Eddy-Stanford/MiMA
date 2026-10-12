!> Mass-weighted global integral of a grid-point field.
module global_integral_mod

  use press_and_geopot_mod, only: half_level_pressures

  use transforms_mod, only: area_weighted_global_mean

  use constants_mod, only: grav

  implicit none
  private

  public :: mass_weighted_global_integral

contains

!---------------------------------------------------------------------------------------------

  !> Returns the mass-weighted vertical integral of `field`, averaged over the globe.
  !>
  !> The units of the result are (units of `field`) * kg/m2.
  function mass_weighted_global_integral(field, surf_press)

    real :: mass_weighted_global_integral  !! global mean of the vertical integral [(units of `field`) kg/m2]
    real, intent(in), dimension(:, :, :) :: field  !! field on the model levels
    real, intent(in), dimension(:, :)   :: surf_press  !! surface pressure [Pa]
    real, dimension(size(field, 1), size(field, 2), size(field, 3)) :: dp
    real, dimension(size(field, 1), size(field, 2), size(field, 3) + 1) :: p_half

    real, dimension(size(field, 1), size(field, 2)) :: vert_integral

    integer :: k, num_levels

    num_levels = size(field, 3)
    p_half = half_level_pressures(surf_press)
    dp = p_half(:, :, 2:num_levels + 1) - p_half(:, :, 1:num_levels)
    vert_integral = 0.
    do k = 1, num_levels
      vert_integral = vert_integral + field(:, :, k)*dp(:, :, k)
    end do
    mass_weighted_global_integral = area_weighted_global_mean(vert_integral)/grav

    return
  end function mass_weighted_global_integral
!---------------------------------------------------------------------------------------------
end module global_integral_mod
