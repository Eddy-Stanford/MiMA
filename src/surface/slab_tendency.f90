!> The implicit step of a slab (mixed-layer) surface temperature, shared by the ocean and
!> land models.
!>
!> The change of the slab temperature over one step follows from the surface energy budget
!> `cap*dT/dt = F + dF/dT*dT`, where `F` is the net downward flux (radiation, sensible and
!> latent heat, the melting of snowfall and any extra heat flux, such as a Q-flux or the
!> local surface heating) and `dF/dT` its derivative with respect to the surface temperature,
!> so that the fluxes are taken at the end of the step (implicit in time).
module slab_tendency_mod

  use constants_mod, only: hlv

  implicit none
  private

  public :: slab_dt

contains

  !> Returns the implicit change of a slab temperature over one step [K].
  !>
  !> `dt_t = flux/(1 - deriv)`, with `flux = (flux_lw + flux_sw - melt - (flux_t + hlv*flux_q)
  !> + q_extra)*dt/cap` and `deriv = -(dhdt + hlv*dedt + drdt)*dt/cap`. Call it only where
  !> `cap` > 0 (for example under `where (mask)`).
  elemental real function slab_dt(flux_lw, flux_sw, melt, flux_t, flux_q, q_extra, &
                                  dhdt, dedt, drdt, dt, cap) result(dt_t)
    real, intent(in) :: flux_lw  !! net downward longwave flux [W/m2]
    real, intent(in) :: flux_sw  !! net downward shortwave flux [W/m2]
    real, intent(in) :: melt  !! heat used to melt the snowfall (`hlf*fprec`) [W/m2]
    real, intent(in) :: flux_t  !! upward sensible heat flux [W/m2]
    real, intent(in) :: flux_q  !! evaporation rate [kg/m2/s]
    real, intent(in) :: q_extra
    !! extra downward heat flux, such as a Q-flux and the local surface heating [W/m2]
    real, intent(in) :: dhdt  !! derivative of `flux_t` with respect to the surface temperature [W/m2/K]
    real, intent(in) :: dedt  !! derivative of `flux_q` with respect to the surface temperature [kg/m2/s/K]
    real, intent(in) :: drdt
    !! derivative of the upward longwave flux with respect to the surface temperature [W/m2/K]
    real, intent(in) :: dt  !! time step [s]
    real, intent(in) :: cap  !! heat capacity [J/m2/K]
    real :: flux, deriv

    flux = (flux_lw + flux_sw - melt - (flux_t + hlv*flux_q) + q_extra)*dt/cap
    deriv = -(dhdt + hlv*dedt + drdt)*dt/cap
    dt_t = flux/(1.0 - deriv)

  end function slab_dt

end module slab_tendency_mod
