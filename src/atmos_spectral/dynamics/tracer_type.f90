!> The derived type that holds the attributes of a prognostic tracer.
module tracer_type_mod

  implicit none
  private

  public :: tracer_type

  !> Attributes of a prognostic tracer, set from the field table by `spectral_dynamics_init`.
  type tracer_type
    character(len=32) :: name, numerical_representation, advect_horiz, advect_vert, hole_filling
    !! `name`: tracer name (lower case); `numerical_representation`: `'spectral'` or `'grid'`;
    !! `advect_horiz`: horizontal advection (`'spectral'` or `'van_leer'`); `advect_vert`: vertical
    !! advection scheme; `hole_filling`: `'on'` or `'off'` (spectral tracers only)
    real :: robert_coeff  !! Robert filter coefficient of the tracer
  end type

end module tracer_type_mod
