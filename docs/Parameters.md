# Parameter settings

All run-time parameters are set in namelists in `input.nml`. This page lists every namelist group that MiMA reads and every variable in it, with its default value in the source code. For the physical meaning of the parameters, see the MiMA reference papers and the comments in the source code.

**The defaults are the standard test case.** Since v2.0 every default in MiMA's own namelists equals the value in the shipped [`input/input.nml`](https://github.com/Eddy-Stanford/MiMA/blob/main/input/input.nml): the RRTM, stationary-wave configuration of [Garfinkel et al. (2020)](https://doi.org/10.1175/JCLI-D-19-0181.1) with Q-fluxes, `cg_drag` and the Navy land-sea mask. An `input.nml` therefore only needs the settings that differ from it, plus the FMS settings listed under [FMS namelists](#fms-namelists). If you are moving from v1, where many defaults were different, see the [migration guide](Migration_v2.md#changed-defaults).

A few rules:

* A variable that is not in a namelist group is a fatal error at start-up ("Unknown namelist, or mistyped namelist variable"). A namelist group that MiMA does not read is ignored.
* A group that is missing from `input.nml` takes all its defaults.
* Every namelist is written, with the values actually used, to `logfile.000000.out`.
* Array variables (e.g. `slandlon`, the `local_heating_nml` variables) are given as comma-separated lists.

Contents:

* [General](#general): `coupler_nml`, `atmos_model_nml`, `spec_mpp_nml`, `diag_integral_nml`
* [Dynamics](#dynamics): `spectral_dynamics_nml`, `vert_coordinate_nml`, `spectral_init_cond_nml`, `transforms_nml`
* [Physics driver](#physics-driver): `physics_driver_nml`
* [Radiation](#radiation): `radiation_nml`, `rrtm_radiation_nml`, `astro_nml`, `gray_radiation_nml`
* [Surface](#surface): `simple_surface_nml`, `qflux_nml`, `surface_flux_nml`, `monin_obukhov_nml`
* [Boundary layer](#boundary-layer): `vert_turb_driver_nml`, `diffusivity_nml`, `vert_diff_driver_nml`
* [Moisture](#moisture): `moist_processes_nml`, `betts_miller_nml`, `lscale_cond_nml`, `moist_conv_nml`
* [Damping and gravity-wave drag](#damping-and-gravity-wave-drag): `damping_driver_nml`, `cg_drag_nml`, `mg_drag_nml`
* [Held-Suarez forcing](#held-suarez-forcing): `held_suarez_nml`
* [Local heating](#local-heating): `local_heating_nml`
* [Tracers and input data](#tracers-and-input-data): `atmos_radon_nml`, `atmos_convection_tracer_nml`, `interpolator_nml`
* [FMS namelists](#fms-namelists): `sat_vapor_pres_nml`, `topography_nml`, `gaussian_topog_nml`, `fms_nml`, `diag_manager_nml`

## General

### `coupler_nml`

Run length, time step and calendar (`coupler/coupler_main.f90`).

 Variable | Default | Meaning
 :--- | :---: | :---
 `current_date` | `1,1,1,0,0,0` | start date (year, month, day, hour, minute, second); used on a cold start, or on a restart if `force_date_from_namelist = .true.`
 `calendar` | `'thirty_day'` | `'thirty_day'` (12 months of 30 days), `'julian'`, `'noleap'` or `'no_calendar'`
 `force_date_from_namelist` | `.false.` | take the date from `current_date` even when `INPUT/coupler.res` exists
 `months`, `days`, `hours`, `minutes`, `seconds` | `0, 360, 0, 0, 0` | run length; it should be a whole multiple of `dt_atmos`
 `dt_atmos` | `500` | model time step [s]
 `atmos_npes` | `0` | number of processes for the atmosphere (0: all)
 `do_atmos` | `.true.` | run the atmosphere

### `atmos_model_nml`

 Variable | Default | Meaning
 :--- | :---: | :---
 `restart_tbot_qbot` | `.false.` | also store the lowest-level temperature and humidity (`t_bot`, `q_bot`) in `atmos_coupled.res.nc`

### `spec_mpp_nml`

 Variable | Default | Meaning
 :--- | :---: | :---
 `io_layout` | `1,1` | I/O layout of the grid domain: each I/O domain (group of processes) writes one file of the diagnostics and restarts. `1,1` gives single files. Each entry must divide the processor layout, which is `1, npes`; with more than one I/O domain the files are split (`atmos_daily.nc.0000`, ...) and must be joined with `mppnccombine`. Split restart files can only be read with the same `io_layout`.

### `diag_integral_nml`

Global integrals printed during the run (`atmos_param/diag_integral/diag_integral.f90`).

 Variable | Default | Meaning
 :--- | :---: | :---
 `output_interval` | `-1.0` | interval at which the integrals are written, in `time_units`; negative: once, at the end of the run (averaged over the whole run)
 `time_units` | `'hours'` | units of `output_interval`
 `file_name` | `' '` | if not blank, write the integrals to this file instead of standard output
 `print_header` | `.true.` | print a header line
 `fields_per_print_line` | `4` | number of fields per line

## Dynamics

### `spectral_dynamics_nml`

The spectral dynamical core (`atmos_spectral/model/spectral_dynamics.f90`).

 Variable | Default | Meaning
 :--- | :---: | :---
 `lon_max`, `lat_max` | `128`, `64` | number of longitudes and latitudes of the Gaussian grid (T42)
 `num_fourier`, `num_spherical` | `42`, `43` | number of zonal and meridional waves retained (T42)
 `fourier_inc` | `1` | if > 1, a sector model with `fourier_inc`-fold symmetry in longitude
 `triang_trunc` | `.true.` | triangular (`.true.`) or rhomboidal truncation
 `num_levels` | `40` | number of vertical levels
 `vert_coord_option` | `'uneven_sigma'` | `'even_sigma'`: equally spaced sigma levels; `'uneven_sigma'`: sigma levels set by `surf_res`, `scale_heights` and `exponent`; `'hybrid'`: as `'uneven_sigma'` with a transition to pressure levels set by `p_sigma` and `p_press`; `'input'`: `pk` and `bk` from `vert_coordinate_nml`
 `surf_res` | `0.1` | parameter 1 of the `uneven_sigma` and `hybrid` levels (resolution near the surface)
 `scale_heights` | `7.9` | parameter 2 of the `uneven_sigma` and `hybrid` levels (model top in scale heights)
 `exponent` | `1.4` | parameter 3 of the `uneven_sigma` and `hybrid` levels
 `p_press`, `p_sigma` | `0.1`, `0.3` | `hybrid` levels: transition between pressure and sigma levels (`p_sigma` > `p_press`)
 `vert_difference_option` | `'simmons_and_burridge'` | vertical differencing; the only option
 `reference_sea_level_press` | `1.e5` | [Pa] cold-start surface pressure and reference pressure of the implicit scheme
 `topography_option` | `'interpolated'` | `'interpolated'`: realistic topography interpolated from the file in `topography_nml`; `'flat'`; `'gaussian'`: idealized mountains from `gaussian_topog_nml`; `'input'`: `zsurf` from `INPUT/topography.data.nc` on the model grid
 `ocean_topog_smoothing` | `0.995` | fractional smoothing of the topography over the ocean (with `'interpolated'`); 0: spectrally truncated but not regularized
 `damping_option` | `'resolution_dependent'` | hyperdiffusion: `'resolution_dependent'` (`damping_coeff` in 1/s, the damping rate of the smallest wave) or `'resolution_independent'`
 `damping_order` | `4` | order of the hyperdiffusion (Laplacian to this power; 4 is del<sup>8</sup>)
 `damping_coeff` | `1.15740741e-4` | [1/s] hyperdiffusion coefficient (one tenth of a day)
 `damping_order_vor`, `damping_coeff_vor` | `-1`, `-1.` | separate values for vorticity (negative: as for the other fields)
 `damping_order_div`, `damping_coeff_div` | `-1`, `-1.` | separate values for divergence (negative: as for the other fields)
 `eddy_sponge_coeff`, `zmu_sponge_coeff`, `zmv_sponge_coeff` | `0.`, `0.`, `0.` | del<sup>2</sup> sponge at the top level for the eddy, zonal-mean zonal and zonal-mean meridional winds (0: off)
 `do_mass_correction` | `.true.` | keep the global mean surface pressure fixed
 `do_water_correction` | `.true.` | keep the global water vapour fixed in the dynamics (below `water_correction_limit`); set `.false.` for dry runs
 `water_correction_limit` | `200.e2` | [Pa] correct water only below this pressure; correcting in the stratosphere introduces an artificial sink there
 `do_energy_correction` | `.true.` | keep the global mean energy fixed
 `use_virtual_temperature` | `.false.` | use virtual temperature in the geopotential
 `use_implicit` | `.true.` | semi-implicit time stepping
 `alpha_implicit` | `0.5` | implicitness of the gravity-wave terms (0.5: centred, 1: backward)
 `robert_coeff` | `0.03` | Robert filter coefficient
 `vert_advect_uv`, `vert_advect_t` | `'second_centered'` | vertical advection of wind and temperature: `'second_centered'`, `'fourth_centered'`, `'van_leer_linear'` or `'finite_volume_parabolic'`
 `longitude_origin` | `0.` | longitude of the first grid point [deg]
 `initial_sphum` | `2.e-6` | [kg/kg] cold-start specific humidity (so the stratosphere does not start dry)
 `valid_range_t` | `100., 500.` | [K] the model stops if the temperature leaves this range
 `num_steps` | `1` | number of dynamics substeps per time step
 `print_interval` | `1, 0` | interval (days, seconds) for printing global integrals of the dynamics
 `specify_initial_conditions` | `.false.` | read the initial state from `INPUT/initial_conditions.nc` on a cold start (see [specified initial conditions](Configurations.md#specified-initial-conditions))
 `add_noise` | `-1.` | if > 0, amplitude [K] of random noise added to the spectral temperature at start-up (see [noise](Configurations.md#adding-noise-to-the-initial-conditions))
 `add_noise_seed` | `-1` | random seed for `add_noise`; set >= 0 for reproducible noise
 `noise_spectral_cutoff_minimum`, `noise_spectral_cutoff_maximum` | `1`, `20` | range of spectral indices that receive noise

### `vert_coordinate_nml`

Used only with `vert_coord_option = 'input'` (`atmos_spectral/init/vert_coordinate.f90`).

 Variable | Default | Meaning
 :--- | :---: | :---
 `pk` | (unset) | [Pa] `num_levels`+1 values; the interface pressures are `pk + bk*ps`
 `bk` | (unset) | `num_levels`+1 sigma values between 0 and 1; `bk(num_levels+1)` must be 1

### `spectral_init_cond_nml`

 Variable | Default | Meaning
 :--- | :---: | :---
 `initial_temperature` | `264.` | [K] temperature of the isothermal atmosphere on a cold start

### `transforms_nml`

 Variable | Default | Meaning
 :--- | :---: | :---
 `check_fourier_imag` | `.false.` | debugging check in the Fourier transforms

## Physics driver

### `physics_driver_nml`

Which physics components are used (`atmos_param/physics_driver/physics_driver.f90`). The radiation scheme is chosen in [`radiation_nml`](#radiation_nml).

 Variable | Default | Meaning
 :--- | :---: | :---
 `do_moist_physics` | `.true.` | convection and large-scale condensation (`moist_processes_nml`)
 `do_boundary_layer` | `.true.` | boundary-layer turbulence, vertical diffusion and coupling to the surface fluxes. With `.false.` the surface state is not updated.
 `do_damping` | `.true.` | Rayleigh sponge and gravity-wave drag (`damping_driver_nml`)
 `do_held_suarez` | `.false.` | add the Held-Suarez (1994) forcing (`held_suarez_nml`)
 `do_local_heating` | `.false.` | add prescribed local heating (`local_heating_nml`)
 `tau_diff` | `3600.` | [s] time scale for smoothing the diffusion coefficients in time
 `diffusion_smooth` | `.true.` | smooth the diffusion coefficients in time
 `diff_min` | `1.e-3` | [m2/s] diffusion coefficients below this are set to zero

## Radiation

### `radiation_nml`

 Variable | Default | Meaning
 :--- | :---: | :---
 `radiation_scheme` | `'rrtm'` | `'rrtm'`: RRTMG clear-sky radiation (`rrtm_radiation_nml`, `astro_nml`); `'gray'`: gray radiation (`gray_radiation_nml`); `'none'`: no radiative heating and no radiative surface fluxes. See [radiation options](Configurations.md#radiation-options).

### `rrtm_radiation_nml`

The RRTM wrapper (`atmos_param/radiation/rrtm/rrtm_radiation.f90`). File names are given without `.nc` and are read from `INPUT/`; the field in the file must have the same name as the file.

 Variable | Default | Meaning
 :--- | :---: | :---
 `do_read_ozone` | `.true.` | read ozone from `ozone_file` (the only way to have non-constant ozone)
 `ozone_file` | `'ozone_1990'` | ozone file (in `input/INPUT/`)
 `scale_ozone` | `1.0` | factor applied to the ozone from the file
 `o3_val` | `0.0` | constant ozone used if `do_read_ozone = .false.`
 `co2ppmv` | `390.` | [ppmv] CO<sub>2</sub> concentration
 `include_secondary_gases` | `.false.` | use the following values for CH<sub>4</sub>, N<sub>2</sub>O, O<sub>2</sub>, CFC-11, CFC-12, CFC-22 and CCl<sub>4</sub> (otherwise they are zero)
 `ch4_val`, `n2o_val`, `o2_val` | `0.` | volume mixing ratios if `include_secondary_gases`
 `cfc11_val`, `cfc12_val`, `cfc22_val`, `ccl4_val` | `0.` | volume mixing ratios if `include_secondary_gases`
 `h2o_lower_limit` | `2.e-7` | [kg/kg] smallest specific humidity passed to RRTM
 `temp_lower_limit`, `temp_upper_limit` | `100.`, `370.` | [K] temperature limits applied before calling RRTM
 `do_zm_tracers` | `.false.` | pass only the zonal mean of the absorbers to RRTM
 `do_zm_rad` | `.false.` | compute radiation for the zonal mean only
 `dt_rad` | `4500` | [s] radiation time step; radiation is computed every step if `dt_rad` < `dt_atmos`
 `store_intermediate_rad` | `.true.` | keep the radiative heating constant between radiation steps (`.false.`: heat only on radiation steps)
 `do_rad_time_avg` | `.true.` | average the solar zenith angle over `dt_rad_avg`
 `dt_rad_avg` | `4500` | [s] averaging interval for the zenith angle; 0: no averaging; < 0: `dt_rad`. 86400 removes the diurnal cycle.
 `lonstep` | `4` | compute radiation only at every `lonstep`-th longitude and interpolate
 `slowdown_rad` | `1.0` | factor on the speed of the seasonal cycle (> 1 faster, < 1 slower)
 `do_precip_albedo` | `.false.` | increase the surface albedo where it rains (a crude cloud effect)
 `precip_albedo` | `0.35` | if so, the albedo of a fully precipitating grid box
 `precip_lat` | `0.0` | [deg] if so, apply only poleward of this latitude
 `precip_albedo_mode` | `'full'` | if so, use total (`'full'`), large-scale (`'lscale'`) or convective (`'conv'`) precipitation

The other RRTMG inputs are not namelist variables. MiMA calls RRTMG without clouds and aerosols:

 RRTMG input | Value
 :--- | :---
 `icld`, `iaer` | 0: no clouds, no aerosols (cloud and aerosol optical properties are zero)
 `h2ovmr` | from the specific humidity (at least `h2o_lower_limit`), every radiation step
 `o3vmr` | from `ozone_file`, or `o3_val`
 `co2vmr` | `co2ppmv`
 `ch4vmr`, `n2ovmr`, `o2vmr`, `cfc11vmr`, `cfc12vmr`, `cfc22vmr`, `ccl4vmr` | the `*_val` values if `include_secondary_gases`, else 0
 `asdir`, `asdif`, `aldir`, `aldif` | the surface albedo from `simple_surface` (modified by `do_precip_albedo`), the same in all four
 `emis` | 1 (black-body surface)
 `coszen` | computed from `astro_nml`, every radiation step
 `adjes` | `solrad`
 `dyofyr` | day of the year if `use_dyofyr`, else 0
 `scon` | `solr_cnst`

### `astro_nml`

The orbit and solar constant for RRTM (`atmos_param/radiation/rrtm/astro.f90`).

 Variable | Default | Meaning
 :--- | :---: | :---
 `solr_cnst` | `1370.` | [W/m2] solar constant
 `obliq` | `23.439` | [deg] obliquity
 `use_dyofyr` | `.false.` | let RRTM compute the Earth-Sun distance from the day of the year (assumes 365 days per year)
 `solrad` | `1.0` | Earth-Sun distance factor if `use_dyofyr = .false.`
 `solday` | `0` | if > 0, perpetual run at this day of the year
 `equinox_day` | `0.25` | fraction of the year at which the March equinox occurs

### `gray_radiation_nml`

The gray radiation of [Frierson, Held and Zurita-Gotor (2006)](https://doi.org/10.1175/JAS3753.1) (`atmos_param/radiation/gray_radiation.f90`). The insolation is annual-mean and zonally symmetric, `solar_constant/4 * (1 + del_sol*P2(lat) + del_sw*sin(lat))`.

 Variable | Default | Meaning
 :--- | :---: | :---
 `solar_constant` | `1360.0` | [W/m2] solar constant
 `del_sol` | `0.0` | amplitude of the second Legendre polynomial in the insolation (equator-to-pole contrast)
 `del_sw` | `0.0` | amplitude of a hemispherically asymmetric insolation term (winter/summer hemisphere)
 `ir_tau_eq`, `ir_tau_pole` | `4.0`, `4.0` | longwave optical depth at the surface at the equator and at the poles
 `linear_tau` | `0.1` | fraction of the longwave optical depth that is linear in pressure (the rest goes as p<sup>4</sup>)
 `atm_abs` | `0.2` | shortwave optical depth of the atmosphere
 `sw_diff` | `0.0` | equator-to-pole reduction of the shortwave optical depth
 `size_pert`, `long_pert`, `del_long` | `0.`, `180.`, `30.` | amplitude [W/m2], longitude [deg] and width [deg] of a zonally localized insolation perturbation
 `fcng_pert`, `lat_pert`, `lon_pert`, `del_lat`, `del_lon` | `0.`, `0.`, `180.`, `30.`, `90.` | amplitude [W/m2], centre [deg] and half-widths [deg] of a localized (Walker-type) insolation perturbation

## Surface

### `simple_surface_nml`

The mixed-layer ocean and surface properties (`coupler/simple_surface.f90`).

 Variable | Default | Meaning
 :--- | :---: | :---
 `surface_choice` | `1` | 1: slab mixed layer (interactive SST); 2: SST fixed at its initial value
 `Tm` | `285.` | [K] initial SST profile `Tm - deltaT*(3 sin²(lat) - 1)/3`; if `Tm` <= 0, a uniform SST of `-Tm`
 `deltaT` | `40.` | [K] equator-to-pole difference of the initial SST
 `do_read_sst` | `.false.` | take the initial SST from `sst_file` (cold start)
 `do_sc_sst` | `.false.` | prescribe the SST from `sst_file` at every step (implies `do_read_sst`)
 `sst_file` | (blank) | SST file name, without `.nc`, in `INPUT/`
 `heat_capacity` | `3.e8` | [J/m2/K] mixed-layer heat capacity (poleward of `heat_cap_limit`)
 `trop_capacity` | `1.e8` | [J/m2/K] heat capacity equatorward of `trop_cap_limit`, varying linearly to `heat_capacity` at `heat_cap_limit`; <= 0: `heat_capacity`
 `trop_cap_limit`, `heat_cap_limit` | `20.`, `60.` | [deg] latitudes of the tropical and extratropical heat capacities
 `np_cap_factor` | `1.` | factor on `heat_capacity` in the Northern Hemisphere
 `land_capacity` | `1.e7` | [J/m2/K] heat capacity over land (as set by `land_option`); <= 0: `heat_capacity`
 `land_option` | `'interpolated'` | where the land is: `'none'`; `'interpolated'`: Navy land-sea mask (the `water_file` of `topography_nml`); `'oceanmaskpole'`: as `'interpolated'`, with the latitude-dependent ocean heat capacity; `'zsurf'`: surface height above `zsurf_cap_limit`; `'lonlat'`: the boxes `slandlon`..`elandlon`, `slandlat`..`elandlat`; `'input'`: land-sea mask file `INPUT/lmask.nc`
 `zsurf_cap_limit` | `10.` | [m] with `land_option = 'zsurf'`, points higher than this are land
 `slandlon`, `elandlon`, `slandlat`, `elandlat` | `0`, `-1`, `0`, `-1` (10 values each) | with `land_option = 'lonlat'`, start and end longitude and latitude [deg] of up to 10 land boxes
 `roughness_choice` | `4` | 1: `const_roughness` everywhere; 3: over land, momentum and moisture roughness multiplied by `mom_roughness_land` and `q_roughness_land`; 4: as 3, with larger moisture roughness over tropical and midlatitude land than over subtropical land. 3 and 4 need `land_option = 'interpolated'` or `'oceanmaskpole'`.
 `const_roughness` | `3.21e-5` | [m] roughness length
 `mom_roughness_land`, `q_roughness_land` | `5.e3`, `1.e-12` | factors for the land roughness with `roughness_choice` = 3, 4
 `albedo_choice` | `7` | 1: `const_albedo`; 2: `higher_albedo` poleward of `lat_glacier` in one hemisphere (NH if `lat_glacier` > 0); 3: `higher_albedo` poleward of `lat_glacier` in both hemispheres; 4: increase as `(lat/90)^albedo_exp`; 5: tanh increase centred at `albedo_cntrNH`, `albedo_cntrSH` with width `albedo_wdth`; 6: sin² increase from equator to pole; 7: as 5, plus `albedo_desert` over the Sahara, Gobi and Australian deserts
 `const_albedo` | `0.23` | surface albedo (low-latitude value for choices 2-7)
 `higher_albedo` | `0.80` | high-latitude albedo for choices 2-7
 `lat_glacier` | `-70.` | [deg] latitude of the albedo step for choices 2 and 3
 `albedo_exp` | `2.` | exponent for choice 4
 `albedo_cntrNH`, `albedo_cntrSH` | `68.`, `64.` | [deg] centre latitudes of the albedo increase for choices 5 and 7
 `albedo_wdth` | `5.` | [deg] width of the albedo increase for choices 5 and 7
 `albedo_desert` | `0.20` | albedo added over the deserts for choice 7
 `do_qflux` | `.true.` | add the meridional ocean heat flux of `qflux_nml`
 `do_warmpool` | `.true.` | add the zonally asymmetric ocean heat fluxes of `qflux_nml`
 `z_ref_heat`, `z_ref_mom` | `2.`, `10.` | [m] reference heights of the diagnostics `t_ref`, `rh_ref` and `u_ref`, `v_ref`

### `qflux_nml`

Prescribed ocean heat fluxes (Q-fluxes) (`atmos_param/qflux/qflux.f90`), used with `do_qflux` and `do_warmpool` in `simple_surface_nml`. The zonally asymmetric fluxes are those of [Garfinkel et al. (2020)](https://doi.org/10.1175/JCLI-D-19-0181.1).

 Variable | Default | Meaning
 :--- | :---: | :---
 `qflux_amp` | `26.` | [W/m2] amplitude of the meridional Q-flux
 `qflux_width` | `16.` | [deg] half-width of the meridional Q-flux
 `warmpool_amp` | `18.` | [W/m2] amplitude of the warm pool
 `warmpool_width` | `35.` | [deg] latitudinal width of the warm pool
 `warmpool_centr` | `0.` | [deg] central latitude of the warm pool
 `warmpool_k` | `1.66666` | zonal wave number of the warm pool
 `warmpool_phase` | `140.` | [deg] longitude phase of the warm pool
 `warmpool_localization_choice` | `3` | 1: cosine in longitude; 2: cosine restricted to the Indo-Pacific, plus Gulf Stream, Kuroshio and tropical Atlantic terms; 3: the localized patterns of Garfinkel et al. (2020). Which of the regional amplitudes below are used depends on this choice (see `qflux.f90`).
 `gulf_k` | `4` | zonal wave number of the Gulf Stream perturbation (choice 2 only)
 `gulf_phase` | `310.` | [deg] longitude phase of the Gulf Stream perturbation (choice 2 only)
 `gulf_amp` | `70.` | [W/m2] Gulf Stream amplitude (choices 2 and 3; with choice 3 it scales a fixed, localized Gulf Stream pattern, and the tropical Atlantic term is only applied if `gulf_amp` > 0)
 `kuroshio_amp` | `40.` | [W/m2] Kuroshio amplitude (choices 2 and 3)
 `trop_atlantic_amp` | `50.` | [W/m2] tropical Atlantic amplitude (choices 2 and 3)
 `Hawaiiextra` | `30.0` | [W/m2] extra flux near Hawaii
 `north_sea_heat` | `0.` | [1] factor on `gulf_amp` for moving heat from Canada to the North Sea (choice 2 only)
 `Pac_ITCZextra` | `0.` | [W/m2] extra flux in the tropical South Pacific (strengthens the local ITCZ)
 `Pac_SPCZextra` | `0.` | [W/m2] extra flux in the subtropical Pacific (modulates the SPCZ)
 `Africaextra` | `0.` | [W/m2] extra flux near the Agulhas current
 `Sampeextra` | `0.` | [W/m2] extra flux off South America

### `surface_flux_nml`

Surface fluxes (`coupler/surface_flux.f90`).

 Variable | Default | Meaning
 :--- | :---: | :---
 `use_virtual_temp` | `.false.` | use virtual potential temperature for the surface-layer stability
 `old_dtaudv` | `.true.` | use the same d(stress)/d(wind) for both wind components
 `alt_gustiness` | `.false.` | alternative gustiness: a lower bound `gust_const` on the wind speed
 `gust_const` | `1.0` | [m/s] see `alt_gustiness`
 `no_neg_q` | `.false.` | set negative lowest-level humidity to zero
 `use_mixing_ratio` | `.false.` | Manabe Climate Model form of the moisture flux (legacy)
 `ncar_ocean_flux` | `.false.` | NCAR (Large and Yeager) ocean flux formulation
 `no_surface_momentum_flux` | `.false.` | switch off the surface momentum flux
 `no_surface_moisture_flux` | `.false.` | switch off the surface moisture flux
 `no_surface_heat_flux` | `.false.` | switch off the surface sensible heat flux
 `no_surface_radiative_flux` | `.false.` | switch off the surface radiative flux

### `monin_obukhov_nml`

Monin-Obukhov similarity for the surface layer (`atmos_param/monin_obukhov/monin_obukhov.f90`).

 Variable | Default | Meaning
 :--- | :---: | :---
 `rich_crit` | `2.0` | critical Richardson number (must be > 0.25)
 `drag_min` | `4.e-5` | minimum drag coefficient
 `neutral` | `.false.` | neutral stability everywhere
 `stable_option` | `1` | stability function for stable conditions (1 or 2)
 `zeta_trans` | `0.5` | transition value of z/L for `stable_option = 2`

## Boundary layer

### `vert_turb_driver_nml`

Boundary-layer turbulence (`atmos_param/vert_turb_driver/vert_turb_driver.f90`). The only scheme is the non-local K-profile scheme of `diffusivity_nml`.

 Variable | Default | Meaning
 :--- | :---: | :---
 `do_diffusivity` | `.true.` | compute diffusion coefficients with the non-local K scheme (`.false.`: no boundary-layer diffusion)
 `do_molecular_diffusion` | `.false.` | add molecular diffusion
 `use_tau` | `.false.` | use the current time level (`.true.`) or the updated values (`.false.`)
 `gust_scheme` | `'constant'` | surface gustiness: `'constant'` (`constant_gust`) or `'beljaars'` (from u* and b*)
 `constant_gust` | `0.` | [m/s] constant gustiness
 `gust_factor` | `1.0` | factor for the `'beljaars'` gustiness

### `diffusivity_nml`

The non-local K scheme (`atmos_param/diffusivity/diffusivity.f90`).

 Variable | Default | Meaning
 :--- | :---: | :---
 `fixed_depth` | `.false.` | use a fixed boundary-layer depth `depth_0`
 `depth_0` | `5000.0` | [m] boundary-layer depth if `fixed_depth`
 `frac_inner` | `0.1` | fraction of the boundary layer that is the surface layer (between 0 and 1)
 `rich_crit_pbl` | `1.0` | critical bulk Richardson number defining the boundary-layer top
 `background_m`, `background_t` | `0.0`, `0.0` | [m2/s] minimum diffusivities for momentum and heat

### `vert_diff_driver_nml`

 Variable | Default | Meaning
 :--- | :---: | :---
 `do_conserve_energy` | `.true.` | heat the air by the dissipation of kinetic energy in the vertical diffusion
 `use_virtual_temp_vert_diff` | `.false.` | use virtual temperature in the vertical diffusion

## Moisture

Following [Frierson (2007)](https://doi.org/10.1175/JAS3935.1), MiMA uses large-scale condensation with the Betts-Miller convection scheme and its "shallower" shallow convection.

### `moist_processes_nml`

 Variable | Default | Meaning
 :--- | :---: | :---
 `do_bm` | `.true.` | Betts-Miller convection (`betts_miller_nml`)
 `do_mca` | `.false.` | moist convective adjustment (`moist_conv_nml`)
 `do_lsc` | `.true.` | large-scale condensation (`lscale_cond_nml`)
 `use_tau` | `.false.` | use the current time level (`.true.`) or the updated values (`.false.`)
 `pdepth` | `150.e2` | [Pa] boundary-layer depth used to decide between rain and snow
 `tfreeze` | `273.16` | [K] freezing temperature for that decision
 `do_gust_cv` | `.false.` | convective gustiness
 `gustmax` | `3.` | [m/s] maximum convective gustiness
 `gustconst` | `10./86400.` | [kg/m2/s] precipitation rate at which convective gustiness starts to matter

### `betts_miller_nml`

 Variable | Default | Meaning
 :--- | :---: | :---
 `tau_bm` | `7200.` | [s] relaxation time
 `rhbm` | `0.7` | relative humidity of the reference profile
 `do_simp` | `.false.` | adjust the time scales so that precipitation is always continuous
 `do_shallower` | `.true.` | shallow convection: choose a smaller depth so that precipitation is zero
 `do_changeqref` | `.false.` | shallow convection: change both q and T so that precipitation is zero
 `do_envsat` | `.false.` | reference humidity relative to the environment (`.true.`) or the parcel (`.false.`)
 `do_taucape` | `.false.` | make `tau_bm` proportional to CAPE<sup>-1/2</sup>
 `capetaubm` | `900.` | [J/kg] CAPE at which the relaxation time is `tau_bm` (with `do_taucape`)
 `tau_min` | `2400.` | [s] minimum relaxation time (with `do_taucape`)
 `buoyancy_kick` | `0.` | [K] temperature added to the parcel at the lowest level

### `lscale_cond_nml`

 Variable | Default | Meaning
 :--- | :---: | :---
 `hc` | `1.0` | relative humidity at which condensation occurs
 `do_evap` | `.true.` | re-evaporate falling precipitation in sub-saturated layers below

### `moist_conv_nml`

Moist convective adjustment (used with `do_mca = .true.`).

 Variable | Default | Meaning
 :--- | :---: | :---
 `HC` | `1.0` | relative humidity at which the adjustment occurs
 `TOLmin`, `TOLmax` | `0.02`, `0.10` | initial and maximum tolerance of the iterative adjustment
 `ITSMOD` | `30` | number of iterations at each tolerance

## Damping and gravity-wave drag

### `damping_driver_nml`

Upper boundary and gravity-wave drag (`atmos_param/damping_driver/damping_driver.f90`), used with `do_damping = .true.` in `physics_driver_nml`.

 Variable | Default | Meaning
 :--- | :---: | :---
 `do_cg_drag` | `.true.` | non-orographic (convective) gravity-wave drag (`cg_drag_nml`)
 `do_rayleigh` | `.false.` | Rayleigh friction (sponge) at the top of the model
 `trayfric` | `-0.5` | Rayleigh friction time scale: [s] if > 0, [days] if < 0
 `sponge_pbottom` | `50.` | [Pa] bottom of the Rayleigh sponge
 `do_mg_drag` | `.false.` | orographic gravity-wave drag (`mg_drag_nml`)
 `do_const_drag` | `.false.` | idealized seasonal "gravity-wave" drag in the stratosphere
 `const_drag_amp` | `3.e-4` | [m/s2] amplitude of the constant drag
 `const_drag_off` | `0.` | offset of its latitudinal profile
 `do_conserve_energy` | `.true.` | heat the air by the momentum lost to the damping

### `cg_drag_nml`

The Alexander and Dunkerton (1999) non-orographic gravity-wave scheme, with the changes of [Cohen et al. (2013)](https://doi.org/10.1175/JAS-D-12-0240.1) and [Garfinkel et al. (2020)](https://doi.org/10.1175/JCLI-D-19-0181.1) (`atmos_param/cg_drag/cg_drag.f90`).

 Variable | Default | Meaning
 :--- | :---: | :---
 `cg_drag_freq` | `21600` | [s] interval between drag calculations; the drag is held fixed in between (and saved in `cg_drag.res.nc`)
 `cg_drag_offset` | `0` | [s] on a cold start, time to the first calculation (0: `cg_drag_freq`)
 `source_level_pressure` | `315.e2` | [Pa] the source level is the highest level with a pressure greater than this at the equator
 `damp_level_pressure` | `0.85e2` | [Pa] momentum flux reaching the model top is deposited from the top down to this level
 `Bt_0` | `0.0043` | [Pa] total source momentum flux poleward of `phi0n`/`phi0s`
 `Bt_eq` | `0.0043` | [Pa] total source momentum flux between `dphis` and `dphin`; it varies linearly to `Bt_0` at `phi0n`/`phi0s`
 `Bt_nh`, `Bt_sh` | `0.00`, `0.00` | [Pa] additional flux in the Northern and Southern Hemisphere extratropics (tanh transition of width `dphin`/`dphis` at `phi0n`/`phi0s`)
 `phi0n`, `phi0s` | `15.`, `-15.` | [deg] latitudes where the flux reaches `Bt_0`
 `dphin`, `dphis` | `10.`, `-10.` | [deg] edges of the tropical band (flux `Bt_eq`, spectrum width `cwtropics`), and widths of the `Bt_nh`/`Bt_sh` transitions
 `Bw`, `Bn` | `0.4`, `0.0` | [m2/s2] amplitudes of the wide and narrow parts of the phase-speed spectrum (`Bn` is 0 in the tropical band)
 `cw`, `cwtropics` | `35.0`, `35.0` | [m/s] half-width of the wide spectrum outside and inside the tropical band
 `cn` | `2.0` | [m/s] half-width of the narrow spectrum
 `flag` | `0` | 1: spectrum peaks at c = 0; 0: at c - u = 0 (always 0 in the tropical band)
 `cmax` | `99.6` | [m/s] maximum phase speed
 `dc` | `1.2` | [m/s] phase-speed resolution
 `nk` | `1` | number of wavelengths in the spectrum
 `kelvin_kludge` | `1.` | factor on the source flux of waves with c - u < 0 in the tropical band

### `mg_drag_nml`

Orographic gravity-wave drag (`atmos_param/mg_drag/mg_drag.f90`), used with `do_mg_drag = .true.`.

 Variable | Default | Meaning
 :--- | :---: | :---
 `xl_mtn` | `1.0e5` | [m] effective mountain length
 `gmax` | `2.0` | order-one tuning parameter (larger: more drag)
 `acoef` | `1.0` | order-one tuning parameter
 `rho` | `1.13` | [kg/m3] standard sea-level air density
 `low_lev_frac` | `0.23` | fraction of the atmosphere (from the bottom) used for the base flux, where no wave breaking is allowed
 `flux_cut_level` | `0.0` | [Pa] above this level the flux divergence is set to zero
 `source_of_sgsmtn` | `'input'` | sub-grid orography: `'input'` (read from `INPUT/mg_drag.res.nc`) or `'computed'` (from the high-resolution topography)
 `do_conserve_energy` | `.false.` | heat the air by the dissipated kinetic energy

## Held-Suarez forcing

### `held_suarez_nml`

The [Held and Suarez (1994)](https://doi.org/10.1175/1520-0477(1994)075<1825:APFTIO>2.0.CO;2) forcing (`atmos_param/held_suarez/held_suarez.f90`), used with `do_held_suarez = .true.`; the defaults are the HS94 values. See [Held-Suarez forcing](Configurations.md#held-suarez-forcing).

 Variable | Default | Meaning
 :--- | :---: | :---
 `t_zero` | `315.` | [K] surface equilibrium temperature at the equator
 `t_strat` | `200.` | [K] minimum (stratospheric) equilibrium temperature
 `delh` | `60.` | [K] equator-to-pole temperature difference
 `delv` | `10.` | [K] vertical potential temperature difference
 `p_ref` | `1.e5` | [Pa] reference pressure
 `sigma_b` | `0.7` | top of the frictional boundary layer (sigma)
 `ka` | `40.` | [days] free-atmosphere relaxation time
 `ks` | `4.` | [days] surface relaxation time at the equator
 `kf` | `1.` | [days] boundary-layer Rayleigh friction time
 `do_rayleigh_friction` | `.true.` | apply the boundary-layer friction
 `do_conserve_energy` | `.false.` | heat the air by the frictional dissipation

## Local heating

### `local_heating_nml`

Prescribed Gaussian heating (`atmos_param/local_heating/local_heating.f90`), used with `do_local_heating = .true.`. Every variable is an array of up to 10 entries, one per heat source.

 Variable | Default | Meaning
 :--- | :---: | :---
 `hamp` | `0.` | [K/day] amplitude of the heating
 `loncenter` | `-1.` | [deg] longitude of the centre; zonally symmetric if < 0
 `lonwidth` | `-1.` | [deg] zonal width, if `loncenter` >= 0
 `lonmove` | `0.` | [deg/day] zonal speed of the source
 `latcenter` | `0.` | [deg] latitude of the centre
 `latwidth` | `15.` | [deg] meridional width
 `latmove` | `0.` | [deg/day] meridional speed of the source
 `pcenter` | `-1.` | [hPa] pressure of the centre; surface heating if < 0
 `pwidth` | `1.` | [log10(hPa)] vertical width; constant in the vertical if < 0
 `pmove` | `0.` | [hPa/day] vertical speed of the source
 `is_periodic` | `.false.` | reset the position periodically (with `tphase` and `tperiod`): periodic in longitude and pressure, back and forth in latitude
 `twidth` | `-1.` | [days] temporal width; constant in time if < 0
 `tphase` | `0.` | [days] temporal phase
 `tperiod` | `-1.` | temporal period: [fraction of a year] if < 0, [days] if > 0

## Tracers and input data

### `atmos_radon_nml`, `atmos_convection_tracer_nml`

Passive tracers of the tracer driver, used only if they are in the `field_table`.

 Variable | Default | Meaning
 :--- | :---: | :---
 `ncopies_radon` (`atmos_radon_nml`) | `9` | number of copies of the radon tracer (`radon`, `radon_2`, ...) looked for in the field table; at most 9
 `ncopies_cnvct_trcr` (`atmos_convection_tracer_nml`) | `9` | number of copies of the convection tracer looked for in the field table; at most 9

### `interpolator_nml`

MiMA's interpolator of climatology files, e.g. ozone and SST (`atmos_shared/interpolator/interpolator.F90`).

 Variable | Default | Meaning
 :--- | :---: | :---
 `read_all_on_init` | `.false.` | read all time levels of a file at initialization
 `verbose` | `0` | amount of diagnostic printout

## FMS namelists

These groups belong to the [FMS library](https://github.com/NOAA-GFDL/FMS) and keep FMS's defaults; see the FMS documentation for their other variables. The shipped `input.nml` sets:

 Group | Variable | Value in `input.nml` | Meaning
 :--- | :--- | :---: | :---
 `sat_vapor_pres_nml` | `do_simple` | `.true.` | **required**: the simple saturation vapour pressure table (Clausius-Clapeyron with constant latent heat, no ice) that MiMA uses. The FMS default (`.false.`, Goff-Gratch) is rejected at start-up.
 `topography_nml` | `topog_file`, `water_file` | `'INPUT/navy_topography.data.nc'`, `'INPUT/navy_pctwater.data.nc'` | high-resolution topography and water fraction, used by `topography_option = 'interpolated'` and `land_option = 'interpolated'`
 `fms_nml` | `domains_stack_size` | `600000` | stack size of the domain communication buffer; may need to be larger at higher resolution
 `diag_manager_nml` | `do_diag_field_log` | not set | `.true.` writes the list of registered diagnostic fields to `diag_field_log.out.0`

Namelist `gaussian_topog_nml` sets idealized mountains for `topography_option = 'gaussian'`. For example, 3 km wave-two mountains in midlatitudes:

```fortran
&gaussian_topog_nml
    height = 3000., 3000.,
    olat   =   45.,   45.,
    olon   =   90.,  270.,
    wlat   =   20.,   20.,
    wlon   =   20.,   20.,
    rlat   =    0.,    0.,
    rlon   =    0.,    0. /
```

`height` is the height [m], `olon`, `olat` the centre [deg], `wlon`, `wlat` the half-widths [deg], and `rlon`, `rlat` the ridge lengths [deg] of each mountain.
