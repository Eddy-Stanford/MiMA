! ============================================================================
!> Surface fluxes of sensible heat, water vapour, momentum and upward longwave
!> radiation, and their derivatives for the implicit coupling with the atmosphere.
!>
!> Uses bulk formulae with drag coefficients from Monin-Obukhov similarity
!> (`mima_monin_obukhov_mod`), or over sea water optionally the NCAR formulation of
!> Large and Yeager (`ncar_ocean_flux`). The surface humidity is saturated except over
!> land. Negative evaporation is set to zero. The momentum, moisture, sensible heat and
!> radiative fluxes can each be switched off.
!>
!> Namelist: `surface_flux_nml`
!> ([namelist reference](https://eddy-stanford.github.io/MiMA/Parameters/#surface_flux_nml)).
!>
!> References:
!>
!> * Large, W. G., and S. G. Yeager, 2004: Diurnal to decadal global forcing for ocean and
!>   sea-ice models: the data sets and flux climatologies. NCAR Technical Note
!>   NCAR/TN-460+STR.
!>
!> Original authors: Steve Klein, Isaac Held, Bruce Wyman; NCAR ocean fluxes: Mike Winton.
module surface_flux_mod
!-----------------------------------------------------------------------
!                   GNU General Public License
!
! This program is free software; you can redistribute it and/or modify it and
! are expected to follow the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2 of
! the License, or (at your option) any later version.
!
! MOM is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
! License for more details.
!
! For the full text of the GNU General Public License,
! write to: Free Software Foundation, Inc.,
!           675 Mass Ave, Cambridge, MA 02139, USA.
! or see:   http://www.gnu.org/licenses/gpl.html
!-----------------------------------------------------------------------

  use fms_mod, only: FATAL, mpp_pe, mpp_root_pe, write_version_number
  use fms_mod, only: check_nml_error, input_nml_file, stdlog
  use mima_monin_obukhov_mod, only: mo_drag, mo_profile
  use sat_vapor_pres_mod, only: lookup_es
  use constants_mod, only: cp_air, hlv, stefan, rdgas, rvgas, grav, vonkarm

  implicit none
  private

! ==== public interface ======================================================
  public surface_flux
! ==== end of public interface ===============================================

  !> Computes the surface fluxes and their derivatives, for 1-D or 2-D arrays of points.
  !>
  !> The arguments are described in `surface_flux_1d`.
  interface surface_flux
!    module procedure surface_flux_0d
    module procedure surface_flux_1d
    module procedure surface_flux_2d
  end interface

!-----------------------------------------------------------------------

  character(len=*), parameter :: version = '$Id: surface_flux.f90,v 12.0.4.1 2005/05/13 18:16:39 pjp Exp $'
  character(len=*), parameter :: tagname = '$Name:  $'

  logical :: do_init = .true.

  real, parameter :: d622 = rdgas/rvgas
  real, parameter :: d378 = 1.-d622
  real, parameter :: hlars = hlv/rvgas
  real, parameter :: gcp = grav/cp_air
  real, parameter :: kappa = rdgas/cp_air
  real            :: d608 = d378/d622
  ! d608 set to zero at initialization if the use of
  ! virtual temperatures is turned off in namelist

! ---- namelist with default values ------------------------------------------
  logical :: no_neg_q = .false.  !! set negative lowest-level humidity to zero
  ! for backwards compatibility
  logical :: use_virtual_temp = .false.  !! use virtual potential temperature for the surface-layer stability
  logical :: alt_gustiness = .false.  !! alternative gustiness: a lower bound `gust_const` on the wind speed
  logical :: old_dtaudv = .true.  !! use the same d(stress)/d(wind) for both wind components
  logical :: use_mixing_ratio = .false.  !! Manabe Climate Model form of the moisture flux (legacy)
  real    :: gust_const = 1.0  !! [m/s] see `alt_gustiness`
  logical :: ncar_ocean_flux = .false.  !! NCAR (Large and Yeager) ocean flux formulation
  logical :: no_surface_momentum_flux = .false. !! switch off the surface momentum flux
  logical :: no_surface_moisture_flux = .false. !! switch off the surface moisture flux
  logical :: no_surface_heat_flux = .false. !! switch off the surface sensible heat flux
  logical :: no_surface_radiative_flux = .false. !! switch off the surface radiative flux
  ! epg: the no_surface_* options

  namelist /surface_flux_nml/ no_neg_q, &
    use_virtual_temp, &
    alt_gustiness, &
    gust_const, &
    old_dtaudv, &
    use_mixing_ratio, &
    ncar_ocean_flux, &
    no_surface_momentum_flux, &
    no_surface_moisture_flux, &
    no_surface_heat_flux, &
    no_surface_radiative_flux

contains

! ============================================================================
  !> Computes the surface fluxes and their derivatives for a 1-D array of points.
  subroutine surface_flux_1d( &
    t_atm, q_atm_in, u_atm, v_atm, p_atm, z_atm, &
    p_surf, t_surf, t_ca, q_surf, &
    u_surf, v_surf, &
    rough_mom, rough_heat, rough_moist, rough_scale, gust, &
    flux_t, flux_q, flux_r, flux_u, flux_v, &
    cd_m, cd_t, cd_q, &
    w_atm, u_star, b_star, q_star, &
    dhdt_surf, dedt_surf, dedq_surf, drdt_surf, &
    dhdt_atm, dedq_atm, dtaudu_atm, dtaudv_atm, &
    dt, land, seawater, avail)
!  slm Mar 28 2002 -- remove agument drag_q since it is just cd_q*wind
! ============================================================================
    ! ---- arguments -----------------------------------------------------------
    logical, intent(in), dimension(:) :: land, seawater, avail
    !! `land`: land points (surface temperature `t_ca`, surface humidity `q_surf` given);
    !! `seawater`: points where the NCAR ocean fluxes are used (with `ncar_ocean_flux`);
    !! `avail`: points where the fluxes are computed (zero elsewhere)
    real, intent(in), dimension(:) :: &
      t_atm, q_atm_in, u_atm, v_atm, &
      p_atm, z_atm, t_ca, &
      p_surf, t_surf, u_surf, v_surf, &
      rough_mom, rough_heat, rough_moist, rough_scale, gust
    !! `t_atm` [K], `q_atm_in` [kg/kg], `u_atm`, `v_atm` [m/s], `p_atm` [Pa]: temperature,
    !! specific humidity, wind and pressure at the lowest model level; `z_atm`: height of the
    !! lowest level above the surface [m]; `t_ca`: canopy air temperature, used over land [K];
    !! `p_surf`: surface pressure [Pa]; `t_surf`: surface temperature [K]; `u_surf`, `v_surf`:
    !! surface velocity [m/s]; `rough_mom`, `rough_heat`, `rough_moist`: roughness lengths for
    !! momentum, heat and moisture [m]; `rough_scale`: roughness length to which the momentum
    !! drag is scaled (`rough_mom`: no scaling) [m]; `gust`: gustiness [m/s]
    real, intent(out), dimension(:) :: &
      flux_t, flux_q, flux_r, flux_u, flux_v, &
      dhdt_surf, dedt_surf, dedq_surf, drdt_surf, &
      dhdt_atm, dedq_atm, dtaudu_atm, dtaudv_atm, &
      w_atm, u_star, b_star, q_star, &
      cd_m, cd_t, cd_q
    !! `flux_t`: sensible heat flux [W/m2]; `flux_q`: evaporation [kg/m2/s]; `flux_r`: upward
    !! longwave flux of the surface [W/m2]; `flux_u`, `flux_v`: zonal and meridional surface
    !! stress on the atmosphere [Pa]; `dhdt_surf`, `dedt_surf`, `drdt_surf`: derivatives of
    !! `flux_t` [W/m2/K], `flux_q` [kg/m2/s/K] and `flux_r` [W/m2/K] with respect to the
    !! surface temperature; `dedq_surf`: derivative of `flux_q` with respect to the surface
    !! humidity, over land [kg/m2/s]; `dhdt_atm`, `dedq_atm`: derivatives of `flux_t` and
    !! `flux_q` with respect to the lowest-level temperature [W/m2/K] and humidity
    !! [kg/m2/s]; `dtaudu_atm`, `dtaudv_atm`: derivatives of `flux_u` and `flux_v` with
    !! respect to the lowest-level wind [kg/m2/s]; `w_atm`: wind speed used for the fluxes,
    !! including gustiness [m/s]; `u_star`: friction velocity [m/s]; `b_star`: buoyancy scale
    !! [m/s2]; `q_star`: moisture scale [kg/kg]; `cd_m`, `cd_t`, `cd_q`: drag coefficients for
    !! momentum, heat and moisture
    real, intent(inout), dimension(:) :: q_surf
    !! surface specific humidity: given over land, diagnosed from the flux on output [kg/kg]
    real, intent(in) :: dt  !! time step [s] (not used)

    ! ---- local constants -----------------------------------------------------
    ! temperature increment and its reciprocal value for comp. of derivatives
    real, parameter:: del_temp = 0.1, del_temp_inv = 1.0/del_temp

    ! ---- local vars ----------------------------------------------------------
    real, dimension(size(t_atm(:))) :: &
      thv_atm, th_atm, tv_atm, thv_surf, &
      e_sat, e_sat1, q_sat, q_sat1, p_ratio, &
      t_surf0, t_surf1, u_dif, v_dif, &
      rho_drag, drag_t, drag_m, drag_q, rho, &
      q_atm, q_surf0, dw_atmdu, dw_atmdv
    logical, dimension(size(t_atm(:))) :: evap_clipped

    integer :: i, nbad

    !epg: these coefficients allow one to zero out the surface fluxes
    !     by default, nothing happens
    real :: cfm = 1.0
    real :: cfq = 1.0
    real :: cft = 1.0
    real :: cfr = 1.0

    if (do_init) call surface_flux_init

    ! epg: this allows you to turn off surface fluxes, which can be useful in
    ! running eddy life cycle experiments
    if (no_surface_momentum_flux) cfm = 0.0
    if (no_surface_moisture_flux) cfq = 0.0
    if (no_surface_heat_flux) cft = 0.0
    if (no_surface_radiative_flux) cfr = 0.0

    !---- use local value of surf temp ----

    t_surf0 = 200.   !  avoids out-of-bounds in es lookup
    where (avail)
      where (land)
        t_surf0 = t_ca
      elsewhere
        t_surf0 = t_surf
      end where
    end where

    t_surf1 = t_surf0 + del_temp

    call lookup_es(t_surf0, e_sat)  ! saturation vapor pressure
    call lookup_es(t_surf1, e_sat1)  ! perturbed  vapor pressure

    if (use_mixing_ratio) then
      ! surface mixing ratio at saturation
      q_sat = d622*e_sat/(p_surf - e_sat)
      q_sat1 = d622*e_sat1/(p_surf - e_sat1)
    else
      q_sat = d622*e_sat/p_surf
      q_sat1 = d622*e_sat1/p_surf
    end if

    ! initilaize surface air humidity according to surface type
    where (land)
      q_surf0 = q_surf ! land calculates it
    elsewhere
      q_surf0 = q_sat  ! everything else assumes saturated sfc humidity
    end where

    ! check for negative atmospheric humidities
    where (avail) q_atm = q_atm_in
    if (no_neg_q) then
      where (avail .and. q_atm_in < 0.0) q_atm = 0.0
    end if

    ! generate information needed by monin_obukhov
    where (avail)
      p_ratio = (p_surf/p_atm)**kappa

      tv_atm = t_atm*(1.0 + d608*q_atm)     ! virtual temperature
      th_atm = t_atm*p_ratio                ! potential T, using p_surf as refernce
      thv_atm = tv_atm*p_ratio                ! virt. potential T, using p_surf as reference
      thv_surf = t_surf0*(1.0 + d608*q_surf0) ! surface virtual (potential) T
!     thv_surf= t_surf0                        ! surface virtual (potential) T -- just for testing tun off the q_surf

      u_dif = u_surf - u_atm                    ! velocity components relative to surface
      v_dif = v_surf - v_atm
    end where

    if (alt_gustiness) then
      do i = 1, size(avail)
        if (.not. avail(i)) cycle
        w_atm(i) = max(sqrt(u_dif(i)**2 + v_dif(i)**2), gust_const)
        ! derivatives of surface wind w.r.t. atm. wind components
        if (w_atm(i) > gust_const) then
          dw_atmdu(i) = u_dif(i)/w_atm(i)
          dw_atmdv(i) = v_dif(i)/w_atm(i)
        else
          dw_atmdu(i) = 0.0
          dw_atmdv(i) = 0.0
        end if
      end do
    else
      where (avail)
        w_atm = sqrt(u_dif*u_dif + v_dif*v_dif + gust*gust)
        ! derivatives of surface wind w.r.t. atm. wind components
        dw_atmdu = u_dif/w_atm
        dw_atmdv = v_dif/w_atm
      end where
    end if

    !  monin-obukhov similarity theory
    call mo_drag(thv_atm, thv_surf, z_atm, &
                 rough_mom, rough_heat, rough_moist, w_atm, &
                 cd_m, cd_t, cd_q, u_star, b_star, avail)

    ! override with ocean fluxes from NCAR calculation
    if (ncar_ocean_flux) then
      call ncar_ocean_fluxes(w_atm, th_atm, t_surf0, q_atm, q_surf0, z_atm, &
                             seawater, cd_m, cd_t, cd_q, u_star, b_star)
    end if

    where (avail)
      ! scale momentum drag coefficient on orographic roughness
      cd_m = cd_m*(log(z_atm/rough_mom + 1)/log(z_atm/rough_scale + 1))**2
      ! surface layer drag coefficients

      ! epg: the coeffieciets cft, cfq, and cfm default to 1.0, but are set to
      !      0.0 if no_surface_fluxes has turned on
      drag_t = cft*cd_t*w_atm
      drag_q = cfq*cd_q*w_atm
      drag_m = cfm*cd_m*w_atm

      ! density
      rho = p_atm/(rdgas*tv_atm)

      ! sensible heat flux
      rho_drag = cp_air*drag_t*rho
      flux_t = rho_drag*(t_surf0 - th_atm)  ! flux of sensible heat (W/m**2)
      dhdt_surf = rho_drag                   ! d(sensible heat flux)/d(surface temperature)
      dhdt_atm = -rho_drag*p_ratio           ! d(sensible heat flux)/d(atmos temperature)

      ! evaporation
      rho_drag = drag_q*rho
      flux_q = rho_drag*(q_surf0 - q_atm) ! flux of water vapor  (Kg/(m**2 s))
      evap_clipped = flux_q < 0.0
      where (evap_clipped) !added by CIG on May 31 2018; never should have negative evaporation
        flux_q = 0.0
      end where

      where (land)
        dedq_surf = rho_drag
        dedt_surf = 0
      elsewhere
        dedq_surf = 0
        dedt_surf = rho_drag*(q_sat1 - q_sat)*del_temp_inv
      end where

      dedq_atm = -rho_drag   ! d(latent heat flux)/d(atmospheric mixing ratio)

      ! where negative evaporation was clipped to zero, the flux does not
      ! depend on the surface or atmospheric state
      where (evap_clipped)
        dedq_surf = 0.0
        dedt_surf = 0.0
        dedq_atm = 0.0
      end where

      q_star = flux_q/(u_star*rho)             ! moisture scale
      ! ask Chris and Steve K if we still want to keep this for diagnostics
      q_surf = q_atm + flux_q/(rho*cd_q*w_atm)   ! surface specific humidity

      ! epg: cfr is 1.0 by default, but set to 0.0 when no_surface_radiative_flux is activated
      flux_r = cfr*stefan*t_surf**4               ! (W/m**2)
      drdt_surf = 4*stefan*t_surf**3               ! d(upward longwave)/d(surface temperature)

      ! stresses
      rho_drag = drag_m*rho
      flux_u = rho_drag*u_dif   ! zonal      component of stress (Nt/m**2)
      flux_v = rho_drag*v_dif   ! meridional component of stress

    elsewhere
      ! zero-out un-available data in output only fields
      flux_t = 0.0
      flux_q = 0.0
      flux_r = 0.0
      flux_u = 0.0
      flux_v = 0.0
      dhdt_surf = 0.0
      dedt_surf = 0.0
      dedq_surf = 0.0
      drdt_surf = 0.0
      dhdt_atm = 0.0
      dedq_atm = 0.0
      u_star = 0.0
      b_star = 0.0
      q_star = 0.0
      q_surf = 0.0
      w_atm = 0.0
    end where

    !CIG - diagnose negative evaporation  - may 31 2018
    !    write (*,*) "delq",   MINVAL( (flux_q) )

    ! calculate d(stress component)/d(atmos wind component)
    dtaudu_atm = 0.0
    dtaudv_atm = 0.0
    if (old_dtaudv) then
      where (avail)
        dtaudv_atm = -rho_drag
        dtaudu_atm = -rho_drag
      end where
    else
      where (avail)
        dtaudu_atm = -cd_m*rho*(dw_atmdu*u_dif + w_atm)
        dtaudv_atm = -cd_m*rho*(dw_atmdv*v_dif + w_atm)
      end where
    end if

  end subroutine surface_flux_1d

!#######################################################################

  subroutine surface_flux_0d( &
    t_atm_0, q_atm_0, u_atm_0, v_atm_0, p_atm_0, z_atm_0, &
    p_surf_0, t_surf_0, t_ca_0, q_surf_0, &
    u_surf_0, v_surf_0, &
    rough_mom_0, rough_heat_0, rough_moist_0, rough_scale_0, gust_0, &
    flux_t_0, flux_q_0, flux_r_0, flux_u_0, flux_v_0, &
    cd_m_0, cd_t_0, cd_q_0, &
    w_atm_0, u_star_0, b_star_0, q_star_0, &
    dhdt_surf_0, dedt_surf_0, dedq_surf_0, drdt_surf_0, &
    dhdt_atm_0, dedq_atm_0, dtaudu_atm_0, dtaudv_atm_0, &
    dt, land_0, seawater_0, avail_0)

    ! ---- arguments -----------------------------------------------------------
    logical, intent(in) :: land_0, seawater_0, avail_0
    real, intent(in) :: &
      t_atm_0, q_atm_0, u_atm_0, v_atm_0, &
      p_atm_0, z_atm_0, t_ca_0, &
      p_surf_0, t_surf_0, u_surf_0, v_surf_0, &
      rough_mom_0, rough_heat_0, rough_moist_0, rough_scale_0, gust_0
    real, intent(out) :: &
      flux_t_0, flux_q_0, flux_r_0, flux_u_0, flux_v_0, &
      dhdt_surf_0, dedt_surf_0, dedq_surf_0, drdt_surf_0, &
      dhdt_atm_0, dedq_atm_0, dtaudu_atm_0, dtaudv_atm_0, &
      w_atm_0, u_star_0, b_star_0, q_star_0, &
      cd_m_0, cd_t_0, cd_q_0
    real, intent(inout) :: q_surf_0
    real, intent(in)    :: dt

    ! ---- local vars ----------------------------------------------------------
    logical, dimension(1) :: land, seawater, avail
    real, dimension(1) :: &
      t_atm, q_atm, u_atm, v_atm, &
      p_atm, z_atm, t_ca, &
      p_surf, t_surf, u_surf, v_surf, &
      rough_mom, rough_heat, rough_moist, rough_scale, gust
    real, dimension(1) :: &
      flux_t, flux_q, flux_r, flux_u, flux_v, &
      dhdt_surf, dedt_surf, dedq_surf, drdt_surf, &
      dhdt_atm, dedq_atm, dtaudu_atm, dtaudv_atm, &
      w_atm, u_star, b_star, q_star, &
      cd_m, cd_t, cd_q
    real, dimension(1) :: q_surf

    avail = .true.

    t_atm(1) = t_atm_0
    q_atm(1) = q_atm_0
    u_atm(1) = u_atm_0
    v_atm(1) = v_atm_0
    p_atm(1) = p_atm_0
    z_atm(1) = z_atm_0
    t_ca(1) = t_ca_0
    p_surf(1) = p_surf_0
    t_surf(1) = t_surf_0
    u_surf(1) = u_surf_0
    v_surf(1) = v_surf_0
    rough_mom(1) = rough_mom_0
    rough_heat(1) = rough_heat_0
    rough_moist(1) = rough_moist_0
    rough_scale(1) = rough_scale_0
    gust(1) = gust_0
    q_surf(1) = q_surf_0
    land(1) = land_0
    seawater(1) = seawater_0
    avail(1) = avail_0

    call surface_flux_1d( &
      t_atm, q_atm, u_atm, v_atm, p_atm, z_atm, &
      p_surf, t_surf, t_ca, q_surf, &
      u_surf, v_surf, &
      rough_mom, rough_heat, rough_moist, rough_scale, gust, &
      flux_t, flux_q, flux_r, flux_u, flux_v, &
      cd_m, cd_t, cd_q, &
      w_atm, u_star, b_star, q_star, &
      dhdt_surf, dedt_surf, dedq_surf, drdt_surf, &
      dhdt_atm, dedq_atm, dtaudu_atm, dtaudv_atm, &
      dt, land, seawater, avail)

    flux_t_0 = flux_t(1)
    flux_q_0 = flux_q(1)
    flux_r_0 = flux_r(1)
    flux_u_0 = flux_u(1)
    flux_v_0 = flux_v(1)
    dhdt_surf_0 = dhdt_surf(1)
    dedt_surf_0 = dedt_surf(1)
    dedq_surf_0 = dedq_surf(1)
    drdt_surf_0 = drdt_surf(1)
    dhdt_atm_0 = dhdt_atm(1)
    dedq_atm_0 = dedq_atm(1)
    dtaudu_atm_0 = dtaudu_atm(1)
    dtaudv_atm_0 = dtaudv_atm(1)
    w_atm_0 = w_atm(1)
    u_star_0 = u_star(1)
    b_star_0 = b_star(1)
    q_star_0 = q_star(1)
    q_surf_0 = q_surf(1)
    cd_m_0 = cd_m(1)
    cd_t_0 = cd_t(1)
    cd_q_0 = cd_q(1)

  end subroutine surface_flux_0d

  !> Computes the surface fluxes for a 2-D array of points, calling `surface_flux_1d` for
  !> each row; the arguments are as in `surface_flux_1d`.
  subroutine surface_flux_2d( &
    t_atm, q_atm_in, u_atm, v_atm, p_atm, z_atm, &
    p_surf, t_surf, t_ca, q_surf, &
    u_surf, v_surf, &
    rough_mom, rough_heat, rough_moist, rough_scale, gust, &
    flux_t, flux_q, flux_r, flux_u, flux_v, &
    cd_m, cd_t, cd_q, &
    w_atm, u_star, b_star, q_star, &
    dhdt_surf, dedt_surf, dedq_surf, drdt_surf, &
    dhdt_atm, dedq_atm, dtaudu_atm, dtaudv_atm, &
    dt, land, seawater, avail)

    ! ---- arguments -----------------------------------------------------------
    logical, intent(in), dimension(:, :) :: land, seawater, avail
    real, intent(in), dimension(:, :) :: &
      t_atm, q_atm_in, u_atm, v_atm, &
      p_atm, z_atm, t_ca, &
      p_surf, t_surf, u_surf, v_surf, &
      rough_mom, rough_heat, rough_moist, rough_scale, gust
    real, intent(out), dimension(:, :) :: &
      flux_t, flux_q, flux_r, flux_u, flux_v, &
      dhdt_surf, dedt_surf, dedq_surf, drdt_surf, &
      dhdt_atm, dedq_atm, dtaudv_atm, dtaudu_atm, &
      w_atm, u_star, b_star, q_star, &
      cd_m, cd_t, cd_q
    real, intent(inout), dimension(:, :) :: q_surf
    real, intent(in) :: dt

    ! ---- local vars -----------------------------------------------------------
    integer :: j

    do j = 1, size(t_atm, 2)
      call surface_flux_1d( &
        t_atm(:, j), q_atm_in(:, j), u_atm(:, j), v_atm(:, j), p_atm(:, j), z_atm(:, j), &
        p_surf(:, j), t_surf(:, j), t_ca(:, j), q_surf(:, j), &
        u_surf(:, j), v_surf(:, j), &
        rough_mom(:, j), rough_heat(:, j), rough_moist(:, j), rough_scale(:, j), gust(:, j), &
        flux_t(:, j), flux_q(:, j), flux_r(:, j), flux_u(:, j), flux_v(:, j), &
        cd_m(:, j), cd_t(:, j), cd_q(:, j), &
        w_atm(:, j), u_star(:, j), b_star(:, j), q_star(:, j), &
        dhdt_surf(:, j), dedt_surf(:, j), dedq_surf(:, j), drdt_surf(:, j), &
        dhdt_atm(:, j), dedq_atm(:, j), dtaudu_atm(:, j), dtaudv_atm(:, j), &
        dt, land(:, j), seawater(:, j), avail(:, j))
    end do
  end subroutine surface_flux_2d

! ============================================================================
!  Initialization of the surface flux module--reads the nml.
!
  subroutine surface_flux_init

! ---- local vars ----------------------------------------------------------
    integer :: unit, ierr, io

    ! read namelist
    read (input_nml_file, nml=surface_flux_nml, iostat=io)
    ierr = check_nml_error(io, 'surface_flux_nml')

    ! write version number
    call write_version_number(version, tagname)

    if (mpp_pe() == mpp_root_pe()) write (stdlog(), nml=surface_flux_nml)

    if (.not. use_virtual_temp) d608 = 0.0

    do_init = .false.

  end subroutine surface_flux_init

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Over-ocean fluxes following Large and Yeager (used in NCAR models)           !
! Coded by Mike Winton (Michael.Winton@noaa.gov)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!
  subroutine ncar_ocean_fluxes(u_del, t, ts, q, qs, z, avail, &
                               cd, ch, ce, ustar, bstar)
    real, intent(in), dimension(:) :: u_del, t, ts, q, qs, z
    logical, intent(in), dimension(:) :: avail
    real, intent(inout), dimension(:) :: cd, ch, ce, ustar, bstar

    real :: cd_n10, ce_n10, ch_n10, cd_n10_rt    ! neutral 10m drag coefficients
    real :: cd_rt                                ! full drag coefficients @ z
    real :: zeta, x2, x, psi_m, psi_h            ! stability parameters
    real :: u, u10, tv, tstar, qstar, z0, xx, stab
    integer, parameter :: n_itts = 2
    integer i, j

    do i = 1, size(u_del(:))
      if (avail(i)) then
        tv = t(i)*(1 + 0.608*q(i)); 
        u = max(u_del(i), 0.5); ! 0.5 m/s floor on wind (undocumented NCAR)
        u10 = u; ! first guess 10m wind

        cd_n10 = (2.7/u10 + 0.142 + 0.0764*u10)/1e3; ! L-Y eqn. 6a
        cd_n10_rt = sqrt(cd_n10); 
        ce_n10 = 34.6*cd_n10_rt/1e3; ! L-Y eqn. 6b
        stab = 0.5 + sign(0.5, t(i) - ts(i))
        ch_n10 = (18.0*stab + 32.7*(1 - stab))*cd_n10_rt/1e3; ! L-Y eqn. 6c

        cd(i) = cd_n10; ! first guess for exchange coeff's at z
        ch(i) = ch_n10; 
        ce(i) = ce_n10; 
        do j = 1, n_itts                                           ! Monin-Obukhov iteration
          cd_rt = sqrt(cd(i)); 
          ustar(i) = cd_rt*u; ! L-Y eqn. 7a
          tstar = (ch(i)/cd_rt)*(t(i) - ts(i)); ! L-Y eqn. 7b
          qstar = (ce(i)/cd_rt)*(q(i) - qs(i)); ! L-Y eqn. 7c
          bstar(i) = grav*(tstar/tv + qstar/(q(i) + 1/0.608)); 
          zeta = vonkarm*bstar(i)*z(i)/(ustar(i)*ustar(i)); ! L-Y eqn. 8a
          zeta = sign(min(abs(zeta), 10.0), zeta); ! undocumented NCAR
          x2 = sqrt(abs(1 - 16*zeta)); ! L-Y eqn. 8b
          x2 = max(x2, 1.0); ! undocumented NCAR
          x = sqrt(x2); 
          if (zeta > 0) then
            psi_m = -5*zeta; ! L-Y eqn. 8c
            psi_h = -5*zeta; ! L-Y eqn. 8c
          else
            psi_m = log((1 + 2*x + x2)*(1 + x2)/8) - 2*(atan(x) - atan(1.0)); ! L-Y eqn. 8d
            psi_h = 2*log((1 + x2)/2); ! L-Y eqn. 8e
          end if

          u10 = u/(1 + cd_n10_rt*(log(z(i)/10) - psi_m)/vonkarm); ! L-Y eqn. 9
          cd_n10 = (2.7/u10 + 0.142 + 0.0764*u10)/1e3; ! L-Y eqn. 6a again
          cd_n10_rt = sqrt(cd_n10); 
          ce_n10 = 34.6*cd_n10_rt/1e3; ! L-Y eqn. 6b again
          stab = 0.5 + sign(0.5, zeta)
          ch_n10 = (18.0*stab + 32.7*(1 - stab))*cd_n10_rt/1e3; ! L-Y eqn. 6c again
          z0 = 10*exp(-vonkarm/cd_n10_rt); ! diagnostic

          xx = (log(z(i)/10) - psi_m)/vonkarm; 
          cd(i) = cd_n10/(1 + cd_n10_rt*xx)**2; ! L-Y 10a
          xx = (log(z(i)/10) - psi_h)/vonkarm; 
          ch(i) = ch_n10/(1 + ch_n10*xx/cd_n10_rt)**2; !       b
          ce(i) = ce_n10/(1 + ce_n10*xx/cd_n10_rt)**2; !       c
        end do
      end if
    end do

  end subroutine ncar_ocean_fluxes

end module surface_flux_mod

