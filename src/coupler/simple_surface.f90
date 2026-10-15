
!> Mixed-layer (slab) ocean and surface properties.
!>
!> Sets the surface albedo, the roughness lengths and the land-sea contrast of the heat
!> capacity, computes the surface fluxes (`surface_flux_mod`) at the current surface
!> temperature, and steps the surface temperature with the implicit coupling to the
!> atmosphere's vertical diffusion. The surface energy budget includes the net radiation,
!> the sensible and latent heat fluxes, the melting of snowfall, prescribed ocean heat fluxes
!> (`qflux_mod`) and optional local surface heating (`local_heating_mod`). The surface
!> temperature can instead be held fixed or prescribed from a file. It also sends the
!> surface diagnostics, including the reference-height diagnostics `t_ref`, `rh_ref`,
!> `u_ref` and `v_ref`, and saves the state in `RESTART/simple_surface.res.nc`.
!>
!> On a cold start the initial SST is `Tm - deltaT*(3 sin^2(lat) - 1)/3` (uniform `-Tm` if
!> `Tm` <= 0), or read from `sst_file`; values of `Tm` of 399 and above select the
!> Aqua-Planet Experiment and other idealized profiles coded in `simple_surface_init`.
!>
!> Namelist: `simple_surface_nml`
!> ([namelist reference](https://eddy-stanford.github.io/MiMA/Parameters/#simple_surface_nml)).
module simple_surface_mod

!use atmos_coupled_mod, only: atmos_boundary_data_type
  use atmos_model_mod, only: atmos_data_type

  use surface_flux_mod, only: surface_flux
  use mima_monin_obukhov_mod, only: mo_profile
  use sat_vapor_pres_mod, only: compute_qs

  use mima_diag_integral_mod, only: diag_integral_field_init, &
                                    sum_diag_integral_field

  use fms_mod, only: input_nml_file, check_nml_error, &
                     error_mesg, FATAL, mpp_pe, mpp_root_pe, &
                     write_version_number, stdlog
  use fms2_io_mod, only: FmsNetcdfFile_t, open_file, close_file, read_data
  use restart_file_mod, only: restart_file_type, open_restart_read, open_restart_write, &
                              close_restart, read_restart_field, write_restart_field, &
                              check_field_size

  use diag_manager_mod, only: register_diag_field, &
                              register_static_field, send_data

  use time_manager_mod, only: time_type, get_time

  use constants_mod, only: rdgas, rvgas, cp_air, hlv, hlf

! mj know about surface topography
  use spectral_dynamics_mod, only: get_surf_geopotential
  use topography_mod, only: get_ocean_mask
! mj read SSTs
  use mima_interpolator_mod, only: interpolate_type, interpolator_init&
       &, CONSTANT, interpolator
!mj q-flux
  use qflux_mod, only: qflux_init, qflux, warmpool
!mj local surface heating
  use physics_driver_mod, only: do_local_heating
  use local_heating_mod, only: horizontal_heating, ngauss, hamp, pcenter

  use mpp_domains_mod, only: mpp_global_field, mpp_get_compute_domain, mpp_get_global_domain

  implicit none
  private

  public :: simple_surface_init, &
            compute_flux, &
            update_simple_surface, &
            simple_surface_end

!-----------------------------------------------------------------------
  character(len=128) :: version = '$Id: simple_surface.f90,v 1.1.2.3 2005/05/21 02:02:04 pjp Exp $'
  character(len=128) :: tagname = '$Name:  $'

!-----------------------------------------------------------------------
!-------- namelist (for diagnostics) ------

  character(len=14), parameter :: mod_name = 'simple_surface'

  integer :: id_drag_moist, id_drag_heat, id_drag_mom, &
             id_rough_moist, id_rough_heat, id_rough_mom, &
             id_u_star, id_b_star, id_u_flux, id_v_flux, id_t_surf, &
             id_t_flux, id_q_flux, id_o_flux, id_r_flux, &
             id_t_atm, id_u_atm, id_v_atm, id_wind, &
             id_t_ref, id_rh_ref, id_u_ref, id_v_ref, &
             id_del_h, id_del_m, id_del_q, id_albedo, id_entrop_evap, &
             id_entrop_shflx, id_entrop_lwflx, &
             id_heat, id_tdt_horiz !mj

  logical :: first_static = .true.
  logical :: do_init = .true.
  logical :: do_surface_heating = .false.

!-----------------------------------------------------------------------

  real ::   z_ref_heat = 2., & !! [m] reference height of the diagnostics `t_ref` and `rh_ref`
          z_ref_mom = 10., & !! [m] reference height of the diagnostics `u_ref` and `v_ref`
          heat_capacity = 3.e08, & !! [J/m2/K] mixed-layer heat capacity (poleward of `heat_cap_limit`)
          land_capacity = 1.e07, & !! [J/m2/K] heat capacity over land (as set by `land_option`); <= 0: `heat_capacity`
          trop_capacity = 1.e08, &
          !! [J/m2/K] heat capacity equatorward of `trop_cap_limit`, varying linearly to
          !! `heat_capacity` at `heat_cap_limit`; <= 0: `heat_capacity`
          trop_cap_limit = 20., & !! [deg] latitude of the tropical heat capacity `trop_capacity`
          heat_cap_limit = 60., & !! [deg] latitude of the extratropical heat capacity `heat_capacity`
          zsurf_cap_limit = 10., & !! [m] with `land_option = 'zsurf'`, points higher than this are land
          np_cap_factor = 1., & !! factor on `heat_capacity` in the Northern Hemisphere
          const_roughness = 3.21e-05, & !! [m] roughness length
          const_albedo = 0.23, & !! surface albedo (low-latitude value for choices 2-7)
          albedo_exp = 2., & !! exponent for choice 4
          albedo_cntrSH = 64., & !! [deg] centre latitude of the Southern Hemisphere albedo increase for choices 5 and 7
          albedo_cntrNH = 68., & !! [deg] centre latitude of the Northern Hemisphere albedo increase for choices 5 and 7
          albedo_desert = 0.20, & !! albedo added over the deserts for choice 7
          albedo_wdth = 5., & !! [deg] width of the albedo increase for choices 5 and 7
          higher_albedo = 0.80, & !! high-latitude albedo for choices 2-7
          lat_glacier = -70., & !! [deg] latitude of the albedo step for choices 2 and 3
          Tm = 285., &
          !! [K] initial SST profile `Tm - deltaT*(3 sin^2(lat) - 1)/3`; if `Tm` <= 0, a uniform SST
          !! of `-Tm`
          deltaT = 40., & !! [K] equator-to-pole difference of the initial SST
          qflux_amp = 30., & !mj
          qflux_width = 16.        !mj
  ! mj: land_capacity, trop_capacity, trop_cap_limit, heat_cap_limit, zsurf_cap_limit,
  ! np_cap_factor, albedo_exp, albedo_cntrSH, albedo_wdth; cig: albedo_cntrNH, albedo_desert
  ! (qflux_amp and qflux_width are not used; the Q-flux parameters are in qflux_nml)
!cig
  real ::   mom_roughness_land = 5.e3, & !! factor for the land momentum roughness with `roughness_choice` = 3, 4
          q_roughness_land = 1.e-12 !! factor for the land moisture roughness with `roughness_choice` = 3, 4

  integer :: surface_choice = 1  !! 1: slab mixed layer (interactive SST); 2: SST fixed at its initial value
  integer :: roughness_choice = 4
  !! `1`: `const_roughness` everywhere; `3`: over land, momentum and moisture roughness
  !! multiplied by `mom_roughness_land` and `q_roughness_land`; `4`: as 3, with larger moisture
  !! roughness over tropical and midlatitude land than over subtropical land. 3 and 4 need
  !! `land_option = 'interpolated'` or `'oceanmaskpole'`.
  integer :: albedo_choice = 7
  !! `1`: `const_albedo`; `2`: `higher_albedo` poleward of `lat_glacier` in one hemisphere (NH if
  !! `lat_glacier` > 0); `3`: `higher_albedo` poleward of `lat_glacier` in both hemispheres; `4`:
  !! increase as `(lat/90)^albedo_exp`; `5`: tanh increase centred at `albedo_cntrNH`,
  !! `albedo_cntrSH` with width `albedo_wdth`; `6`: sin^2 increase from equator to pole; `7`: as
  !! 5, plus `albedo_desert` over the Sahara, Gobi and Australian deserts
  logical :: do_qflux = .true. !! add the meridional ocean heat flux of `qflux_nml`
  logical :: do_warmpool = .true. !! add the zonally asymmetric ocean heat fluxes of `qflux_nml`
  logical :: do_read_sst = .false. !! take the initial SST from `sst_file` (cold start)
  logical :: do_sc_sst = .false. !! prescribe the SST from `sst_file` at every step (implies `do_read_sst`)
  ! mj: do_qflux, do_warmpool, do_read_sst, do_sc_sst
  character(len=256) :: sst_file  !! SST file name, without `.nc`, in `INPUT/`
  character(len=256) :: land_option = 'interpolated'
  !! where the land is: `'none'`; `'interpolated'`: Navy land-sea mask (the `water_file` of
  !! `topography_nml`); `'oceanmaskpole'`: as `'interpolated'`, with the latitude-dependent
  !! ocean heat capacity; `'zsurf'`: surface height above `zsurf_cap_limit`; `'lonlat'`: the
  !! boxes `slandlon`..`elandlon`, `slandlat`..`elandlat`; `'input'`: land-sea mask file
  !! `INPUT/lmask.nc`
  character(len=256) :: land_sea_mask_file = 'lmask'
  real, dimension(10) :: slandlon = 0, slandlat = 0, elandlon = -1, elandlat = -1
  !! with `land_option = 'lonlat'`, start and end longitude and latitude [deg] of up to 10
  !! land boxes

  namelist /simple_surface_nml/ z_ref_heat, z_ref_mom, &
    surface_choice, heat_capacity, &
    land_capacity, trop_capacity, & !mj
    trop_cap_limit, heat_cap_limit, & !mj
    np_cap_factor, zsurf_cap_limit, & !mj
    roughness_choice, const_roughness, &
    albedo_choice, const_albedo, &
    higher_albedo, lat_glacier, &
    Tm, &
    deltaT, mom_roughness_land, q_roughness_land, &  !cig
    do_qflux, do_warmpool, &  !mj
    do_read_sst, do_sc_sst, sst_file, &  !mj
    land_option, slandlon, slandlat, &  !mj
    elandlon, elandlat, &
    albedo_exp, albedo_cntrSH, albedo_cntrNH, albedo_wdth, albedo_desert     !mj

!-----------------------------------------------------------------------

!---- allocatable module storage ------

  real, allocatable, dimension(:, :) :: e_t_n, f_t_delt_n, &
                                        e_q_n, f_q_delt_n

  real, allocatable, dimension(:, :) :: dhdt_surf, dedt_surf, dedq_surf, &
                                        drdt_surf, dhdt_atm, dedq_atm, &
                                        flux_t, flux_q, flux_lw

  real, allocatable, dimension(:, :) :: sst, flux_u, flux_v, flux_o

! mj know about topography
  real, allocatable, dimension(:, :) :: zsurf, land_sea_heat_capacity
!mj read sst and land sea mask from input file
  real, allocatable, dimension(:, :) :: land_sea_mask
  logical, allocatable, dimension(:, :):: lmask_navy
  type(interpolate_type), save :: sst_interp, lmask_interp

contains

!#######################################################################

  !> Sets the roughness lengths and the albedo, and computes the surface fluxes and their
  !> derivatives at the current surface temperature (`surface_flux`); sends the flux
  !> diagnostics.
  !>
  !> The fluxes and derivatives needed by `update_simple_surface` are kept in module storage.
  subroutine compute_flux(dt, Time, Atm, land_frac, &
                          t_surf_atm, albedo, rough_mom, &
                          flux_u_atm, flux_v_atm, dtaudu_atm, &
                          dtaudv_atm, u_star, b_star)

    real, intent(in)  :: dt  !! time step [s]
    type(time_type), intent(in)  :: Time  !! current time
    type(atmos_data_type), intent(in)  :: Atm  !! atmospheric state at the lowest level
    real, dimension(:, :), intent(out) :: albedo, rough_mom, &
                                          land_frac, dtaudu_atm, &
                                          dtaudv_atm, &
                                          flux_u_atm, flux_v_atm, &
                                          u_star, b_star
    !! `albedo`: surface albedo; `rough_mom`: roughness length for momentum [m]; `land_frac`:
    !! land fraction (set to 0); `dtaudu_atm`, `dtaudv_atm`: derivatives of the zonal and
    !! meridional surface stress with respect to the lowest-level wind [kg/m2/s];
    !! `flux_u_atm`, `flux_v_atm`: zonal and meridional surface stress on the atmosphere
    !! [Pa]; `u_star`: friction velocity [m/s]; `b_star`: buoyancy scale [m/s2]

    real, dimension(:, :), intent(out) :: t_surf_atm  !! surface temperature [K]

    real, dimension(size(Atm%t_bot, 1), size(Atm%t_bot, 2)) :: &
      u_surf, v_surf, rough_heat, rough_moist, &
      stomatal, snow, water, max_water, &
      q_star, q_surf, cd_q, cd_t, cd_m, wind

    logical, dimension(size(Atm%t_bot, 1), size(Atm%t_bot, 2)) :: &
      mask, glacier, seawater

    logical :: used
    logical :: ocean_mask_worked

    integer :: j, i
    real :: lat, pi, lon

    pi = 4.0*atan(1.)

!-----------------------------------------------------------------------

    if (do_init) call error_mesg('compute_flux', &
                                 'must call simple_surface_init first', FATAL)

!-----------------------------------------------------------------------
!------ allocate storage also needed in flux_up_to_atmos -----

    ! (still allocated if the previous step did not update the surface)
    if (.not. allocated(e_t_n)) &
      allocate (e_t_n(size(Atm%t_bot, 1), size(Atm%t_bot, 2)), &
                e_q_n(size(Atm%t_bot, 1), size(Atm%t_bot, 2)), &
                f_t_delt_n(size(Atm%t_bot, 1), size(Atm%t_bot, 2)), &
                f_q_delt_n(size(Atm%t_bot, 1), size(Atm%t_bot, 2)), &
                dhdt_surf(size(Atm%t_bot, 1), size(Atm%t_bot, 2)), &
                dedt_surf(size(Atm%t_bot, 1), size(Atm%t_bot, 2)), &
                dedq_surf(size(Atm%t_bot, 1), size(Atm%t_bot, 2)), &
                drdt_surf(size(Atm%t_bot, 1), size(Atm%t_bot, 2)), &
                dhdt_atm(size(Atm%t_bot, 1), size(Atm%t_bot, 2)), &
                dedq_atm(size(Atm%t_bot, 1), size(Atm%t_bot, 2)), &
                flux_t(size(Atm%t_bot, 1), size(Atm%t_bot, 2)), &
                flux_q(size(Atm%t_bot, 1), size(Atm%t_bot, 2)), &
                flux_lw(size(Atm%t_bot, 1), size(Atm%t_bot, 2)))

    u_surf = 0.0
    v_surf = 0.0
    stomatal = 0.0
    snow = 0.0
    water = 1.0
    max_water = 1.0

    mask = .true.
    glacier = .false.
    seawater = .false.

    if (roughness_choice == 1) then
      rough_mom = const_roughness
      rough_heat = const_roughness
      rough_moist = const_roughness
    elseif (roughness_choice == 3) then   !cig: set higher roughness values over land as compared to ocean
      rough_mom = const_roughness
      rough_heat = const_roughness
      rough_moist = const_roughness

      if (trim(land_option) .eq. 'interpolated' .or. trim(land_option) .eq. 'oceanmaskpole') then
        where (.not. lmask_navy) rough_mom = const_roughness*mom_roughness_land
        where (.not. lmask_navy) rough_moist = const_roughness*q_roughness_land
      end if
    elseif (roughness_choice == 4) then
      !cig: set higher roughness values over land as compared to ocean, and more evaporation over tropics and
      !midlatitudes as compared to subtropics
      rough_mom = const_roughness
      rough_heat = const_roughness
      rough_moist = const_roughness
      if (trim(land_option) .eq. 'interpolated' .or. trim(land_option) .eq. 'oceanmaskpole') then
        where (.not. lmask_navy) rough_mom = const_roughness*mom_roughness_land

        do j = 1, size(Atm%t_bot, 2)
          lat = 0.5*(Atm%lat_bnd(j + 1) + Atm%lat_bnd(j))*180./pi
          where (.not. lmask_navy(:, j)) rough_moist(:, j) = const_roughness*q_roughness_land + &
                      &         +(1.e-7)*exp(-abs(lat - 0.)**3./(2*15.)) &
                      &         + (1.e-25)*exp(-abs(lat - 45.)**3./(2*30.)) &
                      &         + (1.e-25)*exp(-abs(lat + 45.)**3./(2*30.))
        end do
      end if

    end if

    if (albedo_choice == 1) then
      albedo = const_albedo
    elseif (albedo_choice == 2) then
      do j = 1, size(Atm%t_bot, 2)
        lat = 0.5*(Atm%lat_bnd(j + 1) + Atm%lat_bnd(j))*180/pi
        ! mj SH or NH only
        if (lat_glacier .ge. 0.) then
          if (lat > lat_glacier) then

            albedo(:, j) = higher_albedo

          else

            albedo(:, j) = const_albedo

          end if
        else
          if (lat < lat_glacier) then

            albedo(:, j) = higher_albedo

          else

            albedo(:, j) = const_albedo

          end if
        end if
      end do
    elseif (albedo_choice == 3) then ! NH and SH albedo step
      do j = 1, size(Atm%t_bot, 2)
        lat = 0.5*(Atm%lat_bnd(j + 1) + Atm%lat_bnd(j))*180/pi

        if (abs(lat) > lat_glacier) then

          albedo(:, j) = higher_albedo

        else

          albedo(:, j) = const_albedo

        end if

      end do
!mj add symmetric higher_albedo - exponential increase from equator to pole
    elseif (albedo_choice == 4) then
      do j = 1, size(Atm%t_bot, 2)
        lat = 0.5*(Atm%lat_bnd(j + 1) + Atm%lat_bnd(j))*180/pi
        lat = abs(lat)
        albedo(:, j) = const_albedo + (higher_albedo - const_albedo)*(lat/90.)**albedo_exp
      end do
!mj add symmetric higher_albedo - tanh increase around albedo_cntr with width
! albedo_wdth. albedo_cntr can differ between the hemispheres
    elseif (albedo_choice .eq. 5) then
      do j = 1, size(Atm%t_bot, 2)
        lat = 0.5*(Atm%lat_bnd(j + 1) + Atm%lat_bnd(j))*180/pi
        lat = abs(lat)
        albedo(:, j) = const_albedo + &
                       (higher_albedo - const_albedo)*0.5*(1 + tanh((lat - albedo_cntrNH)/albedo_wdth)) + &
                       (higher_albedo - const_albedo)*0.5*(1 - tanh((lat + albedo_cntrSH)/albedo_wdth))

      end do
!mj add symmetric higher albedo - sin2 increase from equator to pole
    elseif (albedo_choice .eq. 6) then
      do j = 1, size(Atm%t_bot, 2)
        lat = 0.5*(Atm%lat_bnd(j + 1) + Atm%lat_bnd(j))
        albedo(:, j) = const_albedo + (higher_albedo - const_albedo)* &
                       sin(lat)**2
      end do

    elseif (albedo_choice .eq. 7) then  !cig: as in option 5, but Gobi, Sahara, and Australian deserts have increased albedo
      do j = 1, size(Atm%t_bot, 2)
        lat = 0.5*(Atm%lat_bnd(j + 1) + Atm%lat_bnd(j))*180/pi

        albedo(:, j) = const_albedo + &
                       (higher_albedo - const_albedo)*0.5*(1 + tanh((lat - albedo_cntrNH)/albedo_wdth)) + &
                       (higher_albedo - const_albedo)*0.5*(1 - tanh((lat + albedo_cntrSH)/albedo_wdth))

        do i = 1, size(Atm%t_bot, 1)
          lon = 0.5*(Atm%lon_bnd(i + 1) + Atm%lon_bnd(i))*180/pi

          if ((lon .gt. 118. .and. lon .lt. 145. .and. lat .gt. -30. .and. lat .lt. -19.) .or. &
              (lon .gt. 80. .and. lon .lt. 105. .and. lat .gt. 32. .and. lat .lt. 40.) .or. &
              (lon .gt. 80. .and. lon .lt. 115. .and. lat .gt. 40. .and. lat .lt. 52.) .or. &
              ((lon .gt. 345. .or. lon .lt. 50.) .and. lat .gt. 13. .and. lat .lt. 30.)) then
            albedo(i, j) = const_albedo + albedo_desert

          end if
        end do
      end do

    end if

    cd_t = 0.0
    cd_m = 0.0
    cd_q = 0.0

    t_surf_atm = sst

!  call surface_flux (Atm%t_bot, Atm%q_bot, Atm%u_bot, Atm%v_bot,        & ! Fez
!                     Atm%p_bot, Atm%z_bot,                              &
!                     Atm%p_surf, t_surf_atm, u_surf, v_surf,            &
!                     rough_mom, rough_heat, rough_moist,                &
!                     Atm%gust, stomatal,                                &
!                     snow, water,  max_water,                           &
!                     flux_t, flux_q, flux_lw, flux_u, flux_v,           &
!                     cd_m,   cd_t, cd_q, wind,                          &
!                     u_star, b_star, q_star, q_surf,                    &
!                     dhdt_surf, dedt_surf,  drdt_surf,                  &
!                     dhdt_atm,  dedq_atm,   dtaudv_atm,                 &
!                     dt, mask, glacier)

    call surface_flux(Atm%t_bot, Atm%q_bot, Atm%u_bot, Atm%v_bot, & ! Lima
                      Atm%p_bot, Atm%z_bot, Atm%p_surf, t_surf_atm, &
                      t_surf_atm, & ! Required argument, intent(in). t_surf_atm instead of Land%t_ca
                      q_surf, u_surf, v_surf, &
                      rough_mom, rough_heat, rough_moist, &
                      rough_mom, & ! Required argument, intent(in). rough_mom instead of Land%rough_scale
                      Atm%gust, flux_t, flux_q, flux_lw, flux_u, flux_v, &
                      cd_m, cd_t, cd_q, wind, u_star, b_star, q_star, &
                      dhdt_surf, dedt_surf, &
                      dedq_surf, & ! Required argument, intent(out), but not needed by this model.
                      drdt_surf, dhdt_atm, dedq_atm, &
                      dtaudu_atm, & ! returned to the coupler for the implicit zonal stress
                      dtaudv_atm, dt, & ! Required argument, intent(in). Looks like it should be .false. everywhere.
                      .not. mask, &
                      seawater, & ! Required argument, intent(in). Looks like fudgefactor for salt water. Use .false.
                      mask) ! Required argument, intent(in). Looks like it should be .true. everywhere.

! intent(out):: flux_t, flux_q, flux_lw, flux_u, flux_v, cd_m, cd_t, cd_q, wind, u_star, b_star, q_star,
! intent(out):: dhdt_surf, dedt_surf, dedq_surf, drdt_surf, dhdt_atm, dedq_atm, dtaudu_atm, dtaudv_atm
! intent(inout) :: q_surf
! All others intent(in)

    flux_u_atm = flux_u
    flux_v_atm = flux_v

    land_frac = 0.0

!=======================================================================
!-------------------- diagnostics section ------------------------------

    if (id_wind > 0) used = send_data(id_wind, wind, Time)
    if (id_drag_moist > 0) used = send_data(id_drag_moist, cd_q, Time)
    if (id_drag_heat > 0) used = send_data(id_drag_heat, cd_t, Time)
    if (id_drag_mom > 0) used = send_data(id_drag_mom, cd_m, Time)
    if (id_rough_heat > 0) used = send_data(id_rough_heat, rough_heat, Time)
    if (id_rough_mom > 0) used = send_data(id_rough_mom, rough_mom, Time)
    if (id_rough_moist > 0) used = send_data(id_rough_moist, rough_moist, Time)
    if (id_u_star > 0) used = send_data(id_u_star, u_star, Time)
    if (id_b_star > 0) used = send_data(id_b_star, b_star, Time)
    if (id_t_atm > 0) used = send_data(id_t_atm, Atm%t_bot, Time)
    if (id_u_atm > 0) used = send_data(id_u_atm, Atm%u_bot, Time)
    if (id_v_atm > 0) used = send_data(id_v_atm, Atm%v_bot, Time)
    if (id_albedo > 0) used = send_data(id_albedo, albedo, Time)
    if (id_u_flux > 0) used = send_data(id_u_flux, flux_u, Time)
    if (id_v_flux > 0) used = send_data(id_v_flux, flux_v, Time)

    if (id_t_ref > 0 .or. id_rh_ref > 0 .or. id_u_ref > 0 .or. id_v_ref > 0 .or. &
        id_del_h > 0 .or. id_del_m > 0 .or. id_del_q > 0) &
      call reference_height_diagnostics(Time, Atm, t_surf_atm, q_surf, &
                                        u_surf, v_surf, rough_mom, rough_heat, &
                                        rough_moist, u_star, b_star, q_star)

!=======================================================================

  end subroutine compute_flux

!#######################################################################

  !> Sends the diagnostics at the reference heights `z_ref_heat` and `z_ref_mom`.
  subroutine reference_height_diagnostics(Time, Atm, t_surf, q_surf, &
                                          u_surf, v_surf, rough_mom, &
                                          rough_heat, rough_moist, &
                                          u_star, b_star, q_star)

! Diagnostics at the reference heights z_ref_heat (t_ref, rh_ref) and
! z_ref_mom (u_ref, v_ref) above the surface, computed as in the FMS/AM2
! flux_exchange. Between the surface and the lowest model level (height
! Atm%z_bot) a quantity f is given by the Monin-Obukhov profile, so that
!
!     f(z_ref) = f_surf + (f_atm - f_surf) * del_f,
!     del_f    = (f(z_ref) - f_surf) / (f(z_atm) - f_surf),
!
! where del_m (winds), del_h (temperature) and del_q (specific humidity)
! come from mo_profile, for the same roughness lengths and u_star, b_star
! as the surface fluxes. The surface values are t_surf, the surface
! specific humidity q_surf from surface_flux, and zero wind. Then
!
!     rh_ref = 100 * q_ref / q_sat(t_ref, p_surf),
!     q_sat  = eps*e_s / (p_surf - (1-eps)*e_s),   eps = rdgas/rvgas,
!
! with e_s(t_ref) from the sat_vapor_pres_mod table (compute_qs; with
! do_simple = .true. the table is the simple Clausius-Clapeyron fit over
! liquid water). These fields do not affect the model state.

    type(time_type), intent(in) :: Time
    type(atmos_data_type), intent(in) :: Atm
    real, dimension(:, :), intent(in) :: t_surf, q_surf, u_surf, v_surf, &
                                         rough_mom, rough_heat, rough_moist, &
                                         u_star, b_star, q_star

    real, dimension(size(t_surf, 1), size(t_surf, 2)) :: del_m, del_h, del_q, &
                                                         del_unused, t_ref, q_ref, qs_ref
    logical :: used

! mo_profile returns all three factors at its first height (its second height
! argument is not used), so del_h and del_q are recomputed at z_ref_heat.
    call mo_profile(z_ref_mom, z_ref_heat, Atm%z_bot, &
                    rough_mom, rough_heat, rough_moist, &
                    u_star, b_star, q_star, del_m, del_h, del_q)
    if (z_ref_heat /= z_ref_mom) &
      call mo_profile(z_ref_heat, z_ref_heat, Atm%z_bot, &
                      rough_mom, rough_heat, rough_moist, &
                      u_star, b_star, q_star, del_unused, del_h, del_q)

    t_ref = t_surf + (Atm%t_bot - t_surf)*del_h
    if (id_t_ref > 0) used = send_data(id_t_ref, t_ref, Time)

    if (id_rh_ref > 0) then
      q_ref = q_surf + (Atm%q_bot - q_surf)*del_q
      call compute_qs(t_ref, Atm%p_surf, qs_ref)
      used = send_data(id_rh_ref, 100.*q_ref/qs_ref, Time)
    end if

    if (id_u_ref > 0) used = send_data(id_u_ref, u_surf + (Atm%u_bot - u_surf)*del_m, Time)
    if (id_v_ref > 0) used = send_data(id_v_ref, v_surf + (Atm%v_bot - v_surf)*del_m, Time)

    if (id_del_h > 0) used = send_data(id_del_h, del_h, Time)
    if (id_del_m > 0) used = send_data(id_del_m, del_m, Time)
    if (id_del_q > 0) used = send_data(id_del_q, del_q, Time)

  end subroutine reference_height_diagnostics

!#######################################################################

  !> Steps the surface temperature and completes the implicit coupling with the vertical
  !> diffusion of the atmosphere; sends the surface diagnostics.
  !>
  !> Called after the down sweep of the vertical diffusion. With `surface_choice = 1`
  !> the change of the mixed-layer temperature follows from the implicit surface energy
  !> budget (or from `sst_file` with `do_sc_sst`); with `surface_choice = 2` it is zero. The
  !> local surface heating of `local_heating_nml` is then added. The fluxes are updated
  !> to the new surface temperature.
  subroutine update_simple_surface(dt, Time, Atm, dt_t_atm, dt_q_atm)

    real, intent(in) :: dt  !! time step [s]
    type(time_type), intent(in)  :: Time  !! current time
    type(atmos_data_type), intent(in)  :: Atm
    !! atmospheric state, including the surface radiation and precipitation and the data of
    !! the down sweep of the vertical diffusion (`Atm%Surf_diff`)

    real, dimension(:, :), intent(out) :: dt_t_atm, dt_q_atm
    !! changes of the lowest-level temperature and specific humidity for the up sweep of the
    !! vertical diffusion (passed to `Surf_diff%delta_t`, `Surf_diff%delta_q`)

    real, dimension(size(Atm%t_bot, 1), size(Atm%t_bot, 2)) :: &
      gamma, dtmass, delta_t, delta_q, dflux_t, dflux_q, &
      flux, deriv, dt_t_surf, &
      entrop_evap, entrop_shflx, entrop_lwflx
    real, dimension(size(Atm%t_bot, 1), size(Atm%t_bot, 2)) :: &
      lon2d, lat2d, horiz_heat

    real    :: cp_inv
    logical :: used

! mj input SST
    real, dimension(size(Atm%t_bot, 1), size(Atm%t_bot, 2)) :: sst_new
! mj shallower ocean in tropics, land-sea contrast
    real ::  pi
    integer :: i, j

    pi = 4.*atan(1.)

    flux_lw = Atm%flux_lw - flux_lw

    dtmass = Atm%Surf_Diff%dtmass
    delta_t = Atm%Surf_Diff%delta_t
    delta_q = Atm%Surf_Diff%delta_q
    dflux_t = Atm%Surf_Diff%dflux_t
    dflux_q = Atm%Surf_Diff%dflux_q

    cp_inv = 1.0/cp_air

    ! temperature

    gamma = 1./(1.0 - dtmass*(dflux_t + dhdt_atm*cp_inv))
    e_t_n = dtmass*dhdt_surf*cp_inv*gamma
    f_t_delt_n = (delta_t + dtmass*flux_t*cp_inv)*gamma

    flux_t = flux_t + dhdt_atm*f_t_delt_n
    dhdt_surf = dhdt_surf + dhdt_atm*e_t_n

! moisture

    gamma = 1./(1.0 - dtmass*(dflux_q + dedq_atm))
    e_q_n = dtmass*dedt_surf*gamma
    f_q_delt_n = (delta_q + dtmass*flux_q)*gamma

    flux_q = flux_q + dedq_atm*f_q_delt_n
    dedt_surf = dedt_surf + dedq_atm*e_q_n

    if (surface_choice == 1) then
      if (do_sc_sst) then !mj sst read from input file
        call interpolator(sst_interp, Time, sst_new, trim(sst_file))
        dt_t_surf = sst_new - sst

      else
        flux = (flux_lw + Atm%flux_sw - hlf*Atm%fprec &
                - (flux_t + hlv*flux_q) + flux_o)*dt/land_sea_heat_capacity

        deriv = -(dhdt_surf + hlv*dedt_surf + drdt_surf)*dt/land_sea_heat_capacity
        !flux    = (flux_lw + Atm%flux_sw - hlf*Atm%fprec &
        !        - (flux_t + hlv*flux_q) + flux_o)*dt/heat_capacity

        !   deriv   = - (dhdt_surf + hlv*dedt_surf + drdt_surf)*dt/heat_capacity
! mj end

        dt_t_surf = flux/(1.0 - deriv)

      end if

    elseif (surface_choice == 2) then

      dt_t_surf = 0.0

    end if

    !###############################
    ! additional heating
    !
    if (do_surface_heating) then
      do j = 1, size(Atm%t_bot, 2)
        do i = 1, size(Atm%t_bot, 1)
          lat2d(i, j) = 0.5*(Atm%lat_bnd(j + 1) + Atm%lat_bnd(j))
          lon2d(i, j) = 0.5*(Atm%lon_bnd(i + 1) + Atm%lon_bnd(i))
        end do
      end do
      call horizontal_heating(Time, lon2d, lat2d, horiz_heat)
    else
      horiz_heat = 0.0 ! for diagnostics output
    end if
    dt_t_surf = dt_t_surf + horiz_heat*dt
    !###############################
    ! apply all heating
    sst = sst + dt_t_surf

    flux_t = flux_t + dt_t_surf*dhdt_surf
    flux_q = flux_q + dt_t_surf*dedt_surf
    flux_lw = flux_lw - dt_t_surf*drdt_surf
    dt_t_atm = f_t_delt_n + dt_t_surf*e_t_n
    dt_q_atm = f_q_delt_n + dt_t_surf*e_q_n

!=======================================================================
!-------------------- diagnostics section ------------------------------

    if (id_t_surf > 0) used = send_data(id_t_surf, sst, Time)
    if (id_t_flux > 0) used = send_data(id_t_flux, flux_t, Time)
    if (id_r_flux > 0) used = send_data(id_r_flux, flux_lw, Time)
    if (id_q_flux > 0) used = send_data(id_q_flux, flux_q, Time)
    if (id_o_flux > 0) used = send_data(id_o_flux, flux_o, Time)
    if (id_heat > 0) used = send_data(id_heat, land_sea_heat_capacity, Time)
    if (id_entrop_evap > 0) then
      entrop_evap = flux_q/sst
      used = send_data(id_entrop_evap, entrop_evap, Time)
    end if
    if (id_entrop_shflx > 0) then
      entrop_shflx = flux_t/sst
      used = send_data(id_entrop_shflx, entrop_shflx, Time)
    end if
    if (id_entrop_lwflx > 0) then
      entrop_lwflx = flux_lw/sst
      used = send_data(id_entrop_lwflx, entrop_lwflx, Time)
    end if

    call sum_diag_integral_field('evap', flux_q*86400.)

!=======================================================================
!---- deallocate module storage ----

    deallocate (f_t_delt_n, f_q_delt_n, e_t_n, e_q_n)
    deallocate (dhdt_surf, dedt_surf, dedq_surf, drdt_surf, dhdt_atm, dedq_atm, &
                flux_t, flux_q, flux_lw)

!-----------------------------------------------------------------------

  end subroutine update_simple_surface

!#######################################################################

  !> Initializes the module: reads `simple_surface_nml`, sets up the heat capacity and the
  !> land-sea mask, the initial SST (from `INPUT/simple_surface.res.nc` if it exists) and the
  !> Q-fluxes, and registers the diagnostics.
  subroutine simple_surface_init(Time, Atm)

    type(time_type), intent(in)  :: Time  !! current time
    type(atmos_data_type), intent(in)  :: Atm  !! atmospheric grid and domain

    integer :: ierr, io

    integer :: i, j, k
    real :: xx, xx2, lat, lon, pi, y0
    real :: coslat !mj
    real, dimension(32) :: ssttabl
    ! mj shallower ocean in tropics, land-sea contrast
    real :: loc_cap
    logical :: ocean_mask_worked
    type(FmsNetcdfFile_t) :: mask_file
    real, allocatable, dimension(:, :) :: global_mask
    type(restart_file_type) :: rst
    integer :: is, ie, js, je, nlon, nlat

    pi = 4.0*atan(1.)

    !-----------------------------------------------------------------------
!------ read namelist ------

    read (input_nml_file, nml=simple_surface_nml, iostat=io)
    ierr = check_nml_error(io, 'simple_surface_nml')

!mj make choices compatible
    !if(do_read_sst .or. do_sc_sst) call error_mesg ('simple_surface',  &
    !              'THERE IS A BUG WITH DO_READ_SST, SO I AM STOPPING', FATAL)
    if (do_sc_sst) do_read_sst = .true.
    if (trop_capacity .le. 0.) trop_capacity = heat_capacity
    if (land_capacity .le. 0.) land_capacity = heat_capacity
    if (roughness_choice /= 1 .and. roughness_choice /= 3 .and. roughness_choice /= 4) &
      call error_mesg('simple_surface', 'roughness_choice must be 1, 3 or 4', FATAL)

!--------- write version number and namelist ------------------

    call write_version_number(version, tagname)
    if (mpp_pe() == mpp_root_pe()) then
      write (stdlog(), nml=simple_surface_nml)
    end if

    call diag_integral_field_init('evap', 'f6.3')
    call diag_field_init(Time, Atm%axes(1:2))

    allocate (sst(size(Atm%t_bot, 1), size(Atm%t_bot, 2)))
    allocate (flux_u(size(Atm%t_bot, 1), size(Atm%t_bot, 2)))
    allocate (flux_v(size(Atm%t_bot, 1), size(Atm%t_bot, 2)))
    allocate (flux_o(size(Atm%t_bot, 1), size(Atm%t_bot, 2)))

!mj read fixed SSTs
    if (do_read_sst) then
      call interpolator_init(sst_interp, trim(sst_file)//'.nc', Atm%lon_bnd, Atm%lat_bnd, data_out_of_bounds=(/CONSTANT/))
    end if

! Set up the heat capacity and the land-sea mask whatever the SST choice: the
! mask is also used by roughness_choice = 3, 4 and the heat capacity by the
! heat_capacity diagnostic.
    allocate (land_sea_heat_capacity(size(Atm%t_bot, 1), size(Atm%t_bot, 2)))
    land_sea_heat_capacity = heat_capacity
    !mj ocean depth function of latitude
    if (trop_capacity .ne. heat_capacity .or. np_cap_factor .ne. 1.0) then
      do j = 1, size(Atm%t_bot, 2)
        lat = 0.5*180/pi*(Atm%lat_bnd(j + 1) + Atm%lat_bnd(j))
        if (lat > 0.) then
          loc_cap = heat_capacity*np_cap_factor
        else
          loc_cap = heat_capacity
        end if
        if (abs(lat) < trop_cap_limit) then
          land_sea_heat_capacity(:, j) = trop_capacity
        elseif (abs(lat) < heat_cap_limit) then
          land_sea_heat_capacity(:, j) = trop_capacity*(1.-(abs(lat) - trop_cap_limit)/(heat_cap_limit - trop_cap_limit)) &
                                         + (abs(lat) - trop_cap_limit)/(heat_cap_limit - trop_cap_limit)*loc_cap
        elseif (lat > heat_cap_limit) then
          land_sea_heat_capacity(:, j) = loc_cap
        end if
      end do
    end if

    if (trim(land_option) .eq. 'input') then
      allocate (land_sea_mask(size(Atm%t_bot, 1), size(Atm%t_bot, 2)))
      if (mpp_pe() .eq. mpp_root_pe()) write (*, '(a)') 'Reading land-sea mask from file INPUT/'//trim(land_sea_mask_file)//'.nc'
      if (.not. open_file(mask_file, 'INPUT/'//trim(land_sea_mask_file)//'.nc', 'read')) &
        call error_mesg('simple_surface_init', &
                        'cannot open INPUT/'//trim(land_sea_mask_file)//'.nc', FATAL)
      call mpp_get_compute_domain(Atm%domain, is, ie, js, je)
      call mpp_get_global_domain(Atm%domain, xsize=nlon, ysize=nlat)
      call check_field_size(mask_file, trim(land_sea_mask_file), (/nlon, nlat/))
      allocate (global_mask(nlon, nlat))
      call read_data(mask_file, trim(land_sea_mask_file), global_mask)
      call close_file(mask_file)
      land_sea_mask = global_mask(is:ie, js:je)
      deallocate (global_mask)
      where (land_sea_mask .gt. 0) land_sea_heat_capacity = land_capacity
! mj use navy land-sea mask
    else if (trim(land_option) .eq. 'interpolated') then
      allocate (lmask_navy(size(land_sea_heat_capacity, 1), size(land_sea_heat_capacity, 2)))
      ocean_mask_worked = get_ocean_mask(Atm%lon_bnd, Atm%lat_bnd, lmask_navy)
      if (.not. ocean_mask_worked) then
        call error_mesg('get_ocean_mask', 'land_option="'//trim(land_option)//'"'// &
                        ' and ocean_mask is not present but water data file does not exist', FATAL)
      end if
      where (.not. lmask_navy) land_sea_heat_capacity = land_capacity
! mj land heat capacity function of surface topography
    else if (trim(land_option) .eq. 'zsurf') then
      allocate (zsurf(size(Atm%t_bot, 1), size(Atm%t_bot, 2)))
      call get_surf_geopotential(zsurf)
      where (zsurf > zsurf_cap_limit) land_sea_heat_capacity = land_capacity
      ! mj land heat capacity given in inputfile
! mj land heat capacity given through ?landlon, ?landlat
    else if (trim(land_option) .eq. 'lonlat') then
      do j = 1, size(Atm%t_bot, 2)
        lat = 0.5*180/pi*(Atm%lat_bnd(j + 1) + Atm%lat_bnd(j))
        do i = 1, size(Atm%t_bot, 1)
          lon = 0.5*180/pi*(Atm%lon_bnd(i + 1) + Atm%lon_bnd(i))
          do k = 1, size(slandlat)
            if (lon >= slandlon(k) .and. lon <= elandlon(k) &
                 &.and. lat >= slandlat(k) .and. lat <= elandlat(k)) then
              land_sea_heat_capacity(i, j) = land_capacity
            end if
          end do
        end do
      end do
      ! cig land heat capacity function of ocean_mask (if ocean mask exists), and use MJ's algorithm for deeper ocean mixed layer
      ! depth for poles vs tropics
    else if (trim(land_option) .eq. 'oceanmaskpole' .and. (trop_capacity .ne. heat_capacity .or. np_cap_factor .ne. 1.0)) then

      allocate (lmask_navy(size(land_sea_heat_capacity, 1), size(land_sea_heat_capacity, 2)))
      ocean_mask_worked = get_ocean_mask(Atm%lon_bnd, Atm%lat_bnd, lmask_navy)

      if (.not. ocean_mask_worked) then
        call error_mesg('get_ocean_mask', 'land_option="'//trim(land_option)//'"'// &
                        ' and ocean_mask is not present but water data file does not exist', FATAL)
      end if

      do j = 1, size(Atm%t_bot, 2)
        lat = 0.5*180/pi*(Atm%lat_bnd(j + 1) + Atm%lat_bnd(j))
        if (lat > 0.) then
          loc_cap = heat_capacity*np_cap_factor
        else
          loc_cap = heat_capacity
        end if
        if (abs(lat) < trop_cap_limit) then
          land_sea_heat_capacity(:, j) = trop_capacity
        elseif (abs(lat) < heat_cap_limit) then
          land_sea_heat_capacity(:, j) = trop_capacity*(1.-(abs(lat) - trop_cap_limit)/(heat_cap_limit - trop_cap_limit)) &
                                         + (abs(lat) - trop_cap_limit)/(heat_cap_limit - trop_cap_limit)*loc_cap
        elseif (lat > heat_cap_limit) then
          land_sea_heat_capacity(:, j) = loc_cap
        end if
      end do
      where (.not. lmask_navy) land_sea_heat_capacity = land_capacity
    end if

    if (open_restart_read(rst, 'INPUT/simple_surface.res.nc', Atm%domain)) then
      call read_restart_field(rst, 'sst', sst)
      call read_restart_field(rst, 'flux_u', flux_u)
      call read_restart_field(rst, 'flux_v', flux_v)
      call close_restart(rst)
!mj read fixed SSTs
    else if (do_read_sst) then
      call interpolator(sst_interp, Time, sst, trim(sst_file))
    else
      do j = 1, size(Atm%t_bot, 2)
        lat = 0.5*(Atm%lat_bnd(j + 1) + Atm%lat_bnd(j))
!    xx = 1. - sin(1.5*lat)*sin(1.5*lat)
!    if(abs(lat) .gt. atan(1.0)*60.0/45.0) xx = 0.0
!    sst(:,j) = 273.15 + 27.0*xx
        xx = sin(lat)*sin(lat)
        xx2 = xx*xx
!    sst(:,j) = 305.0 - 10.0*xx2 - 30.0*xx
!    sst(:,j) = 290.0
!    sst(:,j) = 300.0 - 35.0*xx2
! CL functional form
        sst(:, j) = Tm - deltaT*(3.*xx - 1.)/3.
! mj equal temperature
        if (Tm .le. 0.) then
          sst(:, j) = -Tm
        end if
! APE "control+2" experiment
        if (Tm .ge. 399.) then
          sst(:, j) = 273.15 + 29.0*(1.0 - sin(1.40369*lat)*sin(1.40369*lat))
          if (abs(lat) .ge. pi/3.) sst(:, j) = 273.15
        end if
! APE "control" experiment
        if (Tm .ge. 499.) then
          sst(:, j) = 273.15 + 27.0*(1.0 - sin(1.5*lat)*sin(1.5*lat))
          if (abs(lat) .ge. pi/3.) sst(:, j) = 273.15
        end if
! APE "qobs" experiment
        if (Tm .ge. 599.) then
          sst(:, j) = 273.15 + .5*27.*(1.-sin(1.5*lat)*sin(1.5*lat)) &
                      + .5*27.*(1.-sin(1.5*lat)*sin(1.5*lat)*sin(1.5*lat)*sin(1.5*lat))
          if (abs(lat) .ge. pi/3.) sst(:, j) = 273.15
        end if
! APE "qobs+5" experiment
        if (Tm .ge. 699.) then
          sst(:, j) = sst(:, j) + 5.
        end if
! Jian's qobs + El Nino forcing
        if (Tm .ge. 704.) then
          sst(:, j) = sst(:, j) - 5.+2.*(1.-sin(6.*max(min(abs(lat), pi/12.), 0.))**4)
        end if
! Jian's qobs + expanding 1
        if (Tm .ge. 709.) then
          sst(:, j) = sst(:, j) - 2.*(1.-sin(6.*max(min(abs(lat), pi/12.), 0.))**4)
          y0 = 0.
          sst(:, j) = sst(:, j) + 2.*(1.-sin(3.*max(min((abs(lat) - y0), pi/6.), 0.))**4)
        end if
! Jian's qobs + expanding 2
        if (Tm .ge. 714.) then
          sst(:, j) = sst(:, j) - 2.*(1.-sin(3.*max(min((abs(lat) - y0), pi/6.), 0.))**4)
          y0 = 20.*pi/180.
          sst(:, j) = sst(:, j) + 2.*(1.-sin(3.*max(min((abs(lat) - y0), pi/6.), 0.))**4)
        end if
      end do
      if (Tm .ge. 800.) then
        data ssttabl/244.2464, 244.9859, 246.1535, 247.7335, 249.7301, 252.0308, 254.5561, &
          257.2648, 260.1125, 263.0303, 265.9484, 268.7581, 271.4256, 273.9761, &
          276.4128, 278.7544, 281.0249, 283.2028, 285.2909, 287.3179, 289.2411, &
          291.0211, 292.6213, 293.9812, 295.1642, 296.2877, 297.2654, 298.1137, &
          298.8906, 299.5288, 300.2687, 301.1037/
        do j = 1, size(Atm%t_bot, 2)
          lat = 0.5*180./pi*(Atm%lat_bnd(j + 1) + Atm%lat_bnd(j))
          sst(:, j) = ssttabl(32)
          if (abs(lat) .ge. 2.) then
            sst(:, j) = ssttabl(31)
          end if
          if (abs(lat) .ge. 5.) then
            sst(:, j) = ssttabl(30)
          end if
          if (abs(lat) .ge. 8.) then
            sst(:, j) = ssttabl(29)
          end if
          if (abs(lat) .ge. 11.) then
            sst(:, j) = ssttabl(28)
          end if
          if (abs(lat) .ge. 14.) then
            sst(:, j) = ssttabl(27)
          end if
          if (abs(lat) .ge. 17.) then
            sst(:, j) = ssttabl(26)
          end if
          if (abs(lat) .ge. 19.) then
            sst(:, j) = ssttabl(25)
          end if
          if (abs(lat) .ge. 22.) then
            sst(:, j) = ssttabl(24)
          end if
          if (abs(lat) .ge. 25.) then
            sst(:, j) = ssttabl(23)
          end if
          if (abs(lat) .ge. 28.) then
            sst(:, j) = ssttabl(22)
          end if
          if (abs(lat) .ge. 31.) then
            sst(:, j) = ssttabl(21)
          end if
          if (abs(lat) .ge. 33.) then
            sst(:, j) = ssttabl(20)
          end if
          if (abs(lat) .ge. 36.) then
            sst(:, j) = ssttabl(19)
          end if
          if (abs(lat) .ge. 39.) then
            sst(:, j) = ssttabl(18)
          end if
          if (abs(lat) .ge. 42.) then
            sst(:, j) = ssttabl(17)
          end if
          if (abs(lat) .ge. 45.) then
            sst(:, j) = ssttabl(16)
          end if
          if (abs(lat) .ge. 47.) then
            sst(:, j) = ssttabl(15)
          end if
          if (abs(lat) .ge. 50.) then
            sst(:, j) = ssttabl(14)
          end if
          if (abs(lat) .ge. 53.) then
            sst(:, j) = ssttabl(13)
          end if
          if (abs(lat) .ge. 56.) then
            sst(:, j) = ssttabl(12)
          end if
          if (abs(lat) .ge. 58.) then
            sst(:, j) = ssttabl(11)
          end if
          if (abs(lat) .ge. 61.) then
            sst(:, j) = ssttabl(10)
          end if
          if (abs(lat) .ge. 64.) then
            sst(:, j) = ssttabl(9)
          end if
          if (abs(lat) .ge. 67.) then
            sst(:, j) = ssttabl(8)
          end if
          if (abs(lat) .ge. 70.) then
            sst(:, j) = ssttabl(7)
          end if
          if (abs(lat) .ge. 72.) then
            sst(:, j) = ssttabl(6)
          end if
          if (abs(lat) .ge. 75.) then
            sst(:, j) = ssttabl(5)
          end if
          if (abs(lat) .ge. 78.) then
            sst(:, j) = ssttabl(4)
          end if
          if (abs(lat) .ge. 81.) then
            sst(:, j) = ssttabl(3)
          end if
          if (abs(lat) .ge. 84.) then
            sst(:, j) = ssttabl(2)
          end if
          if (abs(lat) .ge. 86.) then
            sst(:, j) = ssttabl(1)
          end if
        end do
      end if
      flux_u = 0.0
      flux_v = 0.0
    end if

    flux_o = 0.

    if (do_qflux .or. do_warmpool) then
      call qflux_init
!mj q-flux as in Merlis et al (2013) [Part II]
      if (do_qflux) call qflux(Atm%lat_bnd, flux_o)
!mj q-flux to create a tropical temperature perturbation
      if (do_warmpool) call warmpool(Atm%lon_bnd, Atm%lat_bnd, flux_o)
    end if

!mj adding local surface heating
    if (do_local_heating) then
      do j = 1, ngauss
        if (hamp(j) .ne. 0. .and. pcenter(j) .lt. 0.) then
          do_surface_heating = .true.
        end if
      end do
    end if

    do_init = .false.

!-----------------------------------------------------------------------

  end subroutine simple_surface_init

!#######################################################################

  subroutine diag_field_init(Time, atmos_axes)

    type(time_type), intent(in) :: Time
    integer, intent(in) :: atmos_axes(2)

    integer :: iref
    character(len=6) :: label_zm, label_zh
    real, dimension(2) :: trange = (/100., 400./), &
                          vrange = (/-400., 400./), &
                          frange = (/-0.01, 1.01/)
!-----------------------------------------------------------------------
!  initializes diagnostic fields that may be output from this module
!  (the id numbers may be referenced anywhere in this module)
!-----------------------------------------------------------------------

!------ labels for diagnostics -------
!  (z_ref_mom, z_ref_heat are namelist variables)

    iref = int(z_ref_mom + 0.5)
    if (real(iref) == z_ref_mom) then
      write (label_zm, 105) iref
      if (iref < 10) write (label_zm, 100) iref
    else
      write (label_zm, 110) z_ref_mom
    end if

    iref = int(z_ref_heat + 0.5)
    if (real(iref) == z_ref_heat) then
      write (label_zh, 105) iref
      if (iref < 10) write (label_zh, 100) iref
    else
      write (label_zh, 110) z_ref_heat
    end if

100 format(i1, ' m', 3x)
105 format(i2, ' m', 2x)
110 format(f4.1, ' m')

    id_wind = &
      register_diag_field(mod_name, 'wind', atmos_axes, Time, &
                          'wind speed for flux calculations', 'm/s', &
                          range=(/0., vrange(2)/))

    id_drag_moist = &
      register_diag_field(mod_name, 'drag_moist', atmos_axes, Time, &
                          'drag coeff for moisture', '1')

    id_drag_heat = &
      register_diag_field(mod_name, 'drag_heat', atmos_axes, Time, &
                          'drag coeff for heat', '1')

    id_drag_mom = &
      register_diag_field(mod_name, 'drag_mom', atmos_axes, Time, &
                          'drag coeff for momentum', '1')

    id_rough_moist = &
      register_diag_field(mod_name, 'rough_moist', atmos_axes, Time, &
                          'surface roughness for moisture', 'm')

    id_rough_heat = &
      register_diag_field(mod_name, 'rough_heat', atmos_axes, Time, &
                          'surface roughness for heat', 'm')

    id_rough_mom = &
      register_diag_field(mod_name, 'rough_mom', atmos_axes, Time, &
                          'surface roughness for momentum', 'm')

    id_u_star = &
      register_diag_field(mod_name, 'u_star', atmos_axes, Time, &
                          'friction velocity', 'm/s')

    id_b_star = &
      register_diag_field(mod_name, 'b_star', atmos_axes, Time, &
                          'buoyancy scale', 'm/s2')

    id_u_flux = &
      register_diag_field(mod_name, 'tau_x', atmos_axes, Time, &
                          'zonal surface stress on the atmosphere (positive eastward)', 'N/m2')

    id_v_flux = &
      register_diag_field(mod_name, 'tau_y', atmos_axes, Time, &
                          'meridional surface stress on the atmosphere (positive northward)', 'N/m2')

    id_t_surf = &
      register_diag_field(mod_name, 't_surf', atmos_axes, Time, &
                          'surface temperature', 'K', &
                          range=trange)

    id_t_flux = &
      register_diag_field(mod_name, 'shflx', atmos_axes, Time, &
                          'sensible heat flux', 'W/m2')

    id_q_flux = &
      register_diag_field(mod_name, 'evap', atmos_axes, Time, &
                          'evaporation rate', 'kg/m2/s')

    id_o_flux = &
      register_diag_field(mod_name, 'oflx', atmos_axes, Time, &
                          'prescribed ocean heat divergence', 'W/m2')

    id_r_flux = &
      register_diag_field(mod_name, 'lwflx', atmos_axes, Time, &
                          'net (down-up) longwave flux', 'W/m2')

    id_t_atm = &
      register_diag_field(mod_name, 't_atm', atmos_axes, Time, &
                          'temperature at btm level', 'K', &
                          range=trange)

    id_u_atm = &
      register_diag_field(mod_name, 'u_atm', atmos_axes, Time, &
                          'u wind component at btm level', 'm/s', &
                          range=vrange)

    id_v_atm = &
      register_diag_field(mod_name, 'v_atm', atmos_axes, Time, &
                          'v wind component at btm level', 'm/s', &
                          range=vrange)

    id_t_ref = &
      register_diag_field(mod_name, 't_ref', atmos_axes, Time, &
                          'air temperature at '//trim(label_zh), 'K', &
                          range=trange)

    id_rh_ref = &
      register_diag_field(mod_name, 'rh_ref', atmos_axes, Time, &
                          'relative humidity at '//trim(label_zh)//' (100 q/q_sat)', 'percent')

    id_u_ref = &
      register_diag_field(mod_name, 'u_ref', atmos_axes, Time, &
                          'zonal wind at '//trim(label_zm), 'm/s', &
                          range=vrange)

    id_v_ref = &
      register_diag_field(mod_name, 'v_ref', atmos_axes, Time, &
                          'meridional wind at '//trim(label_zm), 'm/s', &
                          range=vrange)

    id_del_h = &
      register_diag_field(mod_name, 'del_h', atmos_axes, Time, &
                          'Monin-Obukhov profile factor (T('//trim(label_zh)//')-T_surf)/(T_atm-T_surf)', '1')
    id_del_m = &
      register_diag_field(mod_name, 'del_m', atmos_axes, Time, &
                          'Monin-Obukhov profile factor u('//trim(label_zm)//')/u_atm', '1')
    id_del_q = &
      register_diag_field(mod_name, 'del_q', atmos_axes, Time, &
                          'Monin-Obukhov profile factor (q('//trim(label_zh)//')-q_surf)/(q_atm-q_surf)', '1')
    id_albedo = &
      register_diag_field(mod_name, 'albedo', atmos_axes, Time, &
                          'surface albedo', '1')
    id_heat = & !mj
      register_diag_field(mod_name, 'heat_capacity', atmos_axes, Time, &
                          'mixed layer heat capacity', 'J/m2/K')
    id_entrop_evap = &
      register_diag_field(mod_name, 'entrop_evap', atmos_axes, Time, &
                          'entropy source from evap', 'kg/m2/s/K')

    id_entrop_shflx = &
      register_diag_field(mod_name, 'entrop_shflx', atmos_axes, Time, &
                          'entropy source from SH flux', 'W/m2/K')

    id_entrop_lwflx = &
      register_diag_field(mod_name, 'entrop_lwflx', atmos_axes, Time, &
                          'entropy source from LW flux', 'W/m2/K')

!-----------------------------------------------------------------------

  end subroutine diag_field_init

!########################################################################

  !> Writes the restart file `RESTART/simple_surface.res.nc` (SST and surface stress).
  subroutine simple_surface_end(Atm)

    type(atmos_data_type), intent(in)  :: Atm  !! atmospheric grid and domain
    type(restart_file_type) :: rst

    call open_restart_write(rst, 'RESTART/simple_surface.res.nc', Atm%domain)
    call write_restart_field(rst, 'sst', sst)
    call write_restart_field(rst, 'flux_u', flux_u)
    call write_restart_field(rst, 'flux_v', flux_v)
    call close_restart(rst)

  end subroutine simple_surface_end

!#######################################################################

end module simple_surface_mod
