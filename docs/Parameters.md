# Parameter settings

All run-time parameters are set in namelists in `input.nml`. This page lists every namelist group that MiMA reads and every variable in it, with its default value in the source code. The tables are generated from the doc comments in the source code by `tools/mimadoc` (see [Fortran API reference](FortranAPI.md#writing-doc-comments)); the text around them is written by hand. For the physics behind the parameters, see the MiMA reference papers.

**The defaults are the standard test case.** Since v2.0 every default in MiMA's own namelists equals the value in the shipped [`input/input.nml`](https://github.com/Eddy-Stanford/MiMA/blob/main/input/input.nml): the RRTM, stationary-wave configuration of [Garfinkel et al. (2020)](https://doi.org/10.1175/JCLI-D-19-0181.1) with Q-fluxes, `cg_drag` and the Navy land-sea mask. An `input.nml` therefore only needs the settings that differ from it, plus the FMS settings listed under [FMS namelists](#fms-namelists). If you are moving from v1, where many defaults were different, see the [migration guide](Migration_v2.md#changed-defaults).

A few rules:

* A variable that is not in a namelist group is a fatal error at start-up ("Unknown namelist, or mistyped namelist variable"). A namelist group that MiMA does not read is ignored, **except** `&simple_surface_nml` of v1, which stops the model with a message (see [Migration to v2.0](Migration_v2.md#checks-for-old-settings)).
* A group that is missing from `input.nml` takes all its defaults.
* Every namelist is written, with the values actually used, to `logfile.000000.out`.
* Array variables (e.g. `slandlon`, the [`local_heating_nml`](#local_heating_nml) variables) are given as comma-separated lists.

Contents:

* [General](#general): [`coupler_nml`](#coupler_nml), [`atmos_model_nml`](#atmos_model_nml), [`spec_mpp_nml`](#spec_mpp_nml), [`diag_integral_nml`](#diag_integral_nml)
* [Dynamics](#dynamics): [`spectral_dynamics_nml`](#spectral_dynamics_nml), [`vert_coordinate_nml`](#vert_coordinate_nml), [`spectral_init_cond_nml`](#spectral_init_cond_nml), [`topog_regularization_nml`](#topog_regularization_nml), [`transforms_nml`](#transforms_nml)
* [Physics driver](#physics-driver): [`physics_driver_nml`](#physics_driver_nml)
* [Radiation](#radiation): [`radiation_nml`](#radiation_nml), [`rrtm_radiation_nml`](#rrtm_radiation_nml), [`astro_nml`](#astro_nml), [`gray_radiation_nml`](#gray_radiation_nml)
* [Surface](#surface): [`flux_exchange_nml`](#flux_exchange_nml), [`surface_albedo_nml`](#surface_albedo_nml), [`ocean_model_nml`](#ocean_model_nml), [`qflux_nml`](#qflux_nml), [`land_model_nml`](#land_model_nml), [`surface_flux_nml`](#surface_flux_nml), [`monin_obukhov_nml`](#monin_obukhov_nml)
* [Boundary layer](#boundary-layer): [`vert_turb_driver_nml`](#vert_turb_driver_nml), [`diffusivity_nml`](#diffusivity_nml), [`vert_diff_driver_nml`](#vert_diff_driver_nml)
* [Moisture](#moisture): [`moist_processes_nml`](#moist_processes_nml), [`betts_miller_nml`](#betts_miller_nml), [`lscale_cond_nml`](#lscale_cond_nml), [`moist_conv_nml`](#moist_conv_nml)
* [Damping and gravity-wave drag](#damping-and-gravity-wave-drag): [`damping_driver_nml`](#damping_driver_nml), [`cg_drag_nml`](#cg_drag_nml), [`mg_drag_nml`](#mg_drag_nml)
* [Held-Suarez forcing](#held-suarez-forcing): [`held_suarez_nml`](#held_suarez_nml)
* [Local heating](#local-heating): [`local_heating_nml`](#local_heating_nml)
* [Tracers and input data](#tracers-and-input-data): [`atmos_radon_nml`](#atmos_radon_nml), [`interpolator_nml`](#interpolator_nml)
* [FMS namelists](#fms-namelists): `sat_vapor_pres_nml`, `topography_nml`, `gaussian_topog_nml`, `fms_nml`, `diag_manager_nml`

## General

### `coupler_nml`

Run length, time step and calendar.

<!-- mimadoc:namelist coupler_nml -->

Declared in [`coupler_main`](https://eddy-stanford.github.io/MiMA/api/coupler_main/) (`src/coupler/coupler_main.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `current_date` | integer, dimension(6) | `(/1, 1, 1, 0, 0, 0/)` | start date (year, month, day, hour, minute, second); used on a cold start, or on a restart if `force_date_from_namelist = .true.` |
| `calendar` | character(len=17) | `'thirty_day '` | `'thirty_day'` (12 months of 30 days), `'julian'`, `'noleap'` or `'no_calendar'` |
| `force_date_from_namelist` | logical | `.false.` | take the date from `current_date` even when `INPUT/coupler.res` exists |
| `months`, `days`, `hours`, `minutes`, `seconds` | integer | `0`, `360`, `0`, `0`, `0` | `months`, `days`, `hours`, `minutes`, `seconds`: run length; it should be a whole multiple of `dt_atmos` |
| `dt_atmos` | integer | `500` | [s] model time step |
| `do_atmos` | logical | `.true.` | run the atmosphere |
| `atmos_npes` | integer | `0` | number of processes for the atmosphere (0: all) |

<!-- mimadoc:end -->

### `atmos_model_nml`

<!-- mimadoc:namelist atmos_model_nml -->

Declared in [`atmos_model_mod`](https://eddy-stanford.github.io/MiMA/api/atmos_model_mod/) (`src/atmos_coupled/atmos_model.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `restart_tbot_qbot` | logical | `.false.` | also store the lowest-level temperature and humidity (`t_bot`, `q_bot`) in `atmos_coupled.res.nc` |

<!-- mimadoc:end -->

### `spec_mpp_nml`

<!-- mimadoc:namelist spec_mpp_nml -->

Declared in [`spec_mpp_mod`](https://eddy-stanford.github.io/MiMA/api/spec_mpp_mod/) (`src/atmos_spectral/tools/spec_mpp.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `io_layout` | integer, dimension(2) | `(/1, 1/)` | I/O layout of the grid domain: each I/O domain (group of processes) writes one file of the diagnostics and restarts. `1,1` gives single files. Each entry must divide the processor layout, which is `1, npes`; with more than one I/O domain the files are split (`atmos_daily.nc.0000`, ...) and must be joined with `mppnccombine`. Split restart files can only be read with the same `io_layout`. |
| `spectral_decomposition` | character(len=16) | `'balanced'` | how the zonal indices m are divided among the processes, in contiguous blocks: `'equal'`: blocks of equal numbers of m; `'balanced'`: blocks of equal estimated work. With triangular truncation the Legendre transforms of zonal index m cost about `num_spherical - m` (the number of its retained coefficients), so equal blocks leave the first process with the most work. `'balanced'` gives each m the cost `num_spherical - m + spectral_cost_per_m` and minimizes the largest block cost. It applies only to triangular truncation with `fourier_inc = 1` (otherwise `'equal'` is used). The results and restart files are the same as with `'equal'`, bit for bit, and a run may restart with either. |
| `spectral_cost_per_m` | real | `20.` | for `spectral_decomposition = 'balanced'`: the cost of each zonal index m that does not depend on its number of retained coefficients (the time step and damping work on all `num_spherical + 1` coefficients, the transposes on all latitudes), in units of the cost of one retained coefficient. Measured at T42 on 4 processes: about 16 with 1 OpenMP thread and 28 with 2, as the threads speed up the Legendre transforms more than the rest. |

<!-- mimadoc:end -->

### `diag_integral_nml`

Global integrals printed during the run.

<!-- mimadoc:namelist diag_integral_nml -->

Declared in [`mima_diag_integral_mod`](https://eddy-stanford.github.io/MiMA/api/mima_diag_integral_mod/) (`src/atmos_param/diag_integral/diag_integral.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `output_interval` | real | `-1.0` | interval at which the integrals are written, in `time_units`; negative: once, at the end of the run (averaged over the whole run) |
| `time_units` | character(len=8) | `'hours'` | units of `output_interval`: `'seconds'`, `'minutes'`, `'hours'` or `'days'` |
| `file_name` | character(len=mxch) | `' '` | if not blank, write the integrals to this file instead of standard output |
| `print_header` | logical | `.true.` | print a header line |
| `fields_per_print_line` | integer | `4` | number of fields per line |

<!-- mimadoc:end -->

## Dynamics

### `spectral_dynamics_nml`

The spectral dynamical core.

<!-- mimadoc:namelist spectral_dynamics_nml -->

Declared in [`spectral_dynamics_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_dynamics_mod/) (`src/atmos_spectral/model/spectral_dynamics.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `use_virtual_temperature` | logical | `.false.` | use virtual temperature in the geopotential |
| `damping_option` | character(len=64) | `'resolution_dependent'` | hyperdiffusion: `'resolution_dependent'` (`damping_coeff` in 1/s, the damping rate of the smallest wave) or `'resolution_independent'` |
| `damping_order` | integer | `4` | order of the hyperdiffusion (Laplacian to this power; 4 is `del^8`) |
| `damping_coeff` | real | `1.15740741e-4` | [1/s] hyperdiffusion coefficient (one tenth of a day) |
| `damping_order_vor` | integer | `-1` | separate order for vorticity (negative: as for the other fields) |
| `damping_coeff_vor` | real | `-1.` | separate coefficient for vorticity (negative: as for the other fields) |
| `damping_order_div` | integer | `-1` | separate order for divergence (negative: as for the other fields) |
| `damping_coeff_div` | real | `-1.` | separate coefficient for divergence (negative: as for the other fields) |
| `do_mass_correction` | logical | `.true.` | keep the global mean surface pressure fixed |
| `do_water_correction` | logical | `.true.` | keep the global water vapour fixed in the dynamics (below `water_correction_limit`); set `.false.` for dry runs |
| `do_energy_correction` | logical | `.true.` | keep the global mean energy fixed |
| `vert_advect_uv` | character(len=64) | `'second_centered'` | vertical advection of wind: `'second_centered'`, `'fourth_centered'`, `'van_leer_linear'` or `'finite_volume_parabolic'` |
| `vert_advect_t` | character(len=64) | `'second_centered'` | vertical advection of temperature: `'second_centered'`, `'fourth_centered'`, `'van_leer_linear'` or `'finite_volume_parabolic'` |
| `use_implicit` | logical | `.true.` | semi-implicit time stepping |
| `longitude_origin` | real | `0.` | [rad] longitude of the first grid point |
| `robert_coeff` | real | `.03` | Robert filter coefficient |
| `alpha_implicit` | real | `.5` | implicitness of the gravity-wave terms (0.5: centred, 1: backward) |
| `vert_difference_option` | character(len=64) | `'simmons_and_burridge'` | vertical differencing; the only option |
| `reference_sea_level_press` | real | `1.e5` | [Pa] cold-start surface pressure and reference pressure of the implicit scheme |
| `lon_max` | integer | `128` | number of longitudes of the Gaussian grid (T42) |
| `lat_max` | integer | `64` | number of latitudes of the Gaussian grid (T42) |
| `num_levels` | integer | `40` | number of vertical levels |
| `num_fourier` | integer | `42` | number of zonal waves retained (T42) |
| `num_spherical` | integer | `43` | number of meridional waves retained (T42) |
| `fourier_inc` | integer | `1` | if > 1, a sector model with `fourier_inc`-fold symmetry in longitude |
| `triang_trunc` | logical | `.true.` | triangular (`.true.`) or rhomboidal truncation |
| `topography_option` | character(len=64) | `'interpolated'` | `'interpolated'`: realistic topography interpolated from the file in `topography_nml`; `'flat'`; `'gaussian'`: idealized mountains from `gaussian_topog_nml`; `'input'`: `zsurf` from `INPUT/topography.data.nc` on the model grid |
| `vert_coord_option` | character(len=64) | `'uneven_sigma'` | `'even_sigma'`: equally spaced sigma levels; `'uneven_sigma'`: sigma levels set by `surf_res`, `scale_heights` and `exponent`; `'hybrid'`: as `'uneven_sigma'` with a transition to pressure levels set by `p_sigma` and `p_press`; `'input'`: `pk` and `bk` from `vert_coordinate_nml` |
| `scale_heights` | real | `7.9` | parameter 2 of the `uneven_sigma` and `hybrid` levels (model top in scale heights) |
| `surf_res` | real | `.1` | parameter 1 of the `uneven_sigma` and `hybrid` levels (resolution near the surface) |
| `p_press` | real | `.1` | `hybrid` levels: transition between pressure and sigma levels (`p_sigma` > `p_press`) |
| `p_sigma` | real | `.3` | `hybrid` levels: transition between pressure and sigma levels (`p_sigma` > `p_press`) |
| `exponent` | real | `1.4` | parameter 3 of the `uneven_sigma` and `hybrid` levels |
| `ocean_topog_smoothing` | real | `0.995` | with `'interpolated'`: fraction of the roughness (squared Laplacian) of the topography over the ocean removed by the regularization (`topog_regularization_nml`), 0 to below 1; 0: spectrally truncated but not regularized |
| `initial_sphum` | real | `2.e-06` | [kg/kg] cold-start specific humidity (so the stratosphere does not start dry) |
| `valid_range_t` | real, dimension(2) | `(/100., 500./)` | [K] the model stops if the temperature leaves this range |
| `eddy_sponge_coeff` | real | `0.` | `del^2` sponge at the top level for the eddy winds (0: off) |
| `zmu_sponge_coeff` | real | `0.` | `del^2` sponge at the top level for the zonal-mean zonal wind (0: off) |
| `zmv_sponge_coeff` | real | `0.` | `del^2` sponge at the top level for the zonal-mean meridional wind (0: off) |
| `print_interval` | integer, dimension(2) | `(/1, 0/)` | interval (days, seconds) for printing global integrals of the dynamics |
| `num_steps` | integer | `1` | number of dynamics substeps per time step |
| `noise_spectral_cutoff_minimum` | integer | `1` | lowest spectral index that receives noise |
| `noise_spectral_cutoff_maximum` | integer | `20` | highest spectral index that receives noise |
| `water_correction_limit` | real | `200.e2` | [Pa] correct water only below this pressure; correcting in the stratosphere introduces an artificial sink there |
| `specify_initial_conditions` | logical | `.false.` | read the initial state from `INPUT/initial_conditions.nc` on a cold start (see [specified initial conditions](Configurations.md#specified-initial-conditions)) |
| `add_noise` | real | `-1.` | if > 0, amplitude [K] of random noise added to the spectral temperature at start-up (see [noise](Configurations.md#adding-noise-to-the-initial-conditions)) |
| `add_noise_seed` | integer | `-1` | random seed for `add_noise`; set >= 0 for reproducible noise |

<!-- mimadoc:end -->

### `vert_coordinate_nml`

Used only with `vert_coord_option = 'input'`.

<!-- mimadoc:namelist vert_coordinate_nml -->

Declared in [`vert_coordinate_mod`](https://eddy-stanford.github.io/MiMA/api/vert_coordinate_mod/) (`src/atmos_spectral/init/vert_coordinate.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `pk`, `bk` | real, dimension(max_levels + 1) | unset, unset | `pk`: [Pa] `num_levels`+1 values; the interface pressures are `pk + bk*ps`. `bk`: `num_levels`+1 sigma values between 0 and 1; `bk(num_levels+1)` must be 1 |

<!-- mimadoc:end -->

### `spectral_init_cond_nml`

<!-- mimadoc:namelist spectral_init_cond_nml -->

Declared in [`spectral_init_cond_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_init_cond_mod/) (`src/atmos_spectral/init/spectral_init_cond.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `initial_temperature` | real | `264.` | [K] temperature of the isothermal atmosphere on a cold start |

<!-- mimadoc:end -->

### `topog_regularization_nml`

The regularization of the topography over the ocean after Lindberg and Broccoli (1996), used with `topography_option = 'interpolated'` and `ocean_topog_smoothing` > 0 in [`spectral_dynamics_nml`](#spectral_dynamics_nml). It removes the Gibbs ripples of the spectrally truncated topography over the ocean while keeping it close to the truncated topography elsewhere; `ocean_topog_smoothing` is the fraction of the ocean roughness (the squared Laplacian) that it removes. The defaults work at any resolution, so this group is normally not needed; `verbose = .true.` shows the progress of the search for the smoothing parameter. See [`topog_regularization_mod`](https://eddy-stanford.github.io/MiMA/api/topog_regularization_mod/) for the method.

<!-- mimadoc:namelist topog_regularization_nml -->

Declared in [`topog_regularization_mod`](https://eddy-stanford.github.io/MiMA/api/topog_regularization_mod/) (`src/atmos_spectral/init/topog_regularization.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `fraction_tolerance` | real | `1.e-4` | the search for `lambda` stops when the fraction of the ocean roughness removed is within this of `ocean_topog_smoothing` |
| `solver_tolerance` | real | `1.e-10` | the conjugate-gradient solver stops when the norm of the residual is below this times the norm of the unsmoothed field |
| `max_solver_iterations` | integer | `2000` | largest number of conjugate-gradient iterations |
| `max_lambda_iterations` | integer | `30` | largest number of solves in the search for `lambda` |
| `verbose` | logical | `.false.` | print the progress of the solver and of the `lambda` search |

<!-- mimadoc:end -->

### `transforms_nml`

<!-- mimadoc:namelist transforms_nml -->

Declared in [`transforms_mod`](https://eddy-stanford.github.io/MiMA/api/transforms_mod/) (`src/atmos_spectral/tools/transforms.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `check_fourier_imag` | logical | `.false.` | debugging check in the Fourier transforms: stop if the imaginary part of the m = 0 or m = `num_lon`/2 Fourier coefficient is not zero |

<!-- mimadoc:end -->

## Physics driver

### `physics_driver_nml`

Which physics components are used. The radiation scheme is chosen in [`radiation_nml`](#radiation_nml).

<!-- mimadoc:namelist physics_driver_nml -->

Declared in [`physics_driver_mod`](https://eddy-stanford.github.io/MiMA/api/physics_driver_mod/) (`src/atmos_param/physics_driver/physics_driver.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `tau_diff` | real | `3600.` | [s] time scale for smoothing the diffusion coefficients in time |
| `diff_min` | real | `1.e-3` | [m2/s] diffusion coefficients below this are set to zero |
| `diffusion_smooth` | logical | `.true.` | smooth the diffusion coefficients in time |
| `do_damping` | logical | `.true.` | Rayleigh sponge and gravity-wave drag (`damping_driver_nml`) |
| `do_local_heating` | logical | `.false.` | add prescribed local heating (`local_heating_nml`) |
| `do_held_suarez` | logical | `.false.` | add the Held-Suarez (1994) forcing (`held_suarez_nml`) |
| `do_boundary_layer` | logical | `.true.` | boundary-layer turbulence, vertical diffusion and coupling to the surface fluxes. With `.false.` the surface state is not updated. |
| `do_moist_physics` | logical | `.true.` | convection and large-scale condensation (`moist_processes_nml`) |

<!-- mimadoc:end -->

## Radiation

### `radiation_nml`

<!-- mimadoc:namelist radiation_nml -->

Declared in [`radiation_mod`](https://eddy-stanford.github.io/MiMA/api/radiation_mod/) (`src/atmos_param/radiation/radiation.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `radiation_scheme` | character(len=16) | `'rrtm'` | `'rrtm'`: RRTMG clear-sky radiation (`rrtm_radiation_nml`, `astro_nml`); `'gray'`: gray radiation (`gray_radiation_nml`); `'none'`: no radiative heating and no radiative surface fluxes. See [radiation options](Configurations.md#radiation-options). |

<!-- mimadoc:end -->

### `rrtm_radiation_nml`

The RRTM wrapper. File names are given without `.nc` and are read from `INPUT/`; the field in the file must have the same name as the file.

<!-- mimadoc:namelist rrtm_radiation_nml -->

Declared in [`rrtm_radiation`](https://eddy-stanford.github.io/MiMA/api/rrtm_radiation/) (`src/atmos_param/radiation/rrtm/rrtm_radiation.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `include_secondary_gases` | logical | `.false.` | use the following values for CH<sub>4</sub>, N<sub>2</sub>O, O<sub>2</sub>, CFC-11, CFC-12, CFC-22 and CCl<sub>4</sub> (otherwise they are zero) |
| `do_read_ozone` | logical | `.true.` | read ozone from `ozone_file` (the only way to have non-constant ozone) |
| `ozone_file` | character(len=256) | `'ozone_1990'` | ozone file (in `input/INPUT/`) |
| `scale_ozone` | real(kind=rb) | `1.0` | factor applied to the ozone from the file |
| `o3_val` | real(kind=rb) | `0.0` | constant ozone used if `do_read_ozone = .false.` |
| `ch4_val` | real(kind=rb) | `0.` | CH<sub>4</sub> volume mixing ratio if `include_secondary_gases` |
| `n2o_val` | real(kind=rb) | `0.` | N<sub>2</sub>O volume mixing ratio if `include_secondary_gases` |
| `o2_val` | real(kind=rb) | `0.` | O<sub>2</sub> volume mixing ratio if `include_secondary_gases` |
| `cfc11_val` | real(kind=rb) | `0.` | CFC-11 volume mixing ratio if `include_secondary_gases` |
| `cfc12_val` | real(kind=rb) | `0.` | CFC-12 volume mixing ratio if `include_secondary_gases` |
| `cfc22_val` | real(kind=rb) | `0.` | CFC-22 volume mixing ratio if `include_secondary_gases` |
| `ccl4_val` | real(kind=rb) | `0.` | CCl<sub>4</sub> volume mixing ratio if `include_secondary_gases` |
| `h2o_lower_limit` | real(kind=rb) | `2.e-7` | [kg/kg] smallest specific humidity passed to RRTM |
| `temp_lower_limit` | real(kind=rb) | `100.` | [K] lower temperature limit applied before calling RRTM |
| `temp_upper_limit` | real(kind=rb) | `370.` | [K] upper temperature limit applied before calling RRTM |
| `co2ppmv` | real(kind=rb) | `390.` | [ppmv] CO<sub>2</sub> concentration |
| `slowdown_rad` | real(kind=rb) | `1.0` | factor on the speed of the seasonal cycle (> 1 faster, < 1 slower) |
| `store_intermediate_rad` | logical | `.true.` | keep the radiative heating constant between radiation steps (`.false.`: heat only on radiation steps) |
| `do_rad_time_avg` | logical | `.true.` | average the solar zenith angle over `dt_rad_avg` |
| `dt_rad` | integer(kind=im) | `4500` | [s] radiation time step; radiation is computed every step if `dt_rad` < `dt_atmos` |
| `dt_rad_avg` | integer(kind=im) | `4500` | [s] averaging interval for the zenith angle; 0: no averaging; < 0: `dt_rad`. 86400 removes the diurnal cycle. |
| `lonstep` | integer(kind=im) | `4` | compute radiation only at every `lonstep`-th longitude and interpolate |
| `balance_sw_columns` | logical | `.true.` | share the sunlit columns of the shortwave calculation evenly among the processes. RRTMG skips the columns where the sun is down, so the processes whose latitude bands are more in daylight have more work and the others wait for them. With `.true.` processes with more sunlit columns than the mean send the extra columns to processes with fewer, which compute them and send the fluxes back. The results are the same as with `.false.`, bit for bit. |
| `do_zm_tracers` | logical | `.false.` | pass only the zonal mean of the absorbers to RRTM |
| `do_zm_rad` | logical | `.false.` | compute radiation for the zonal mean only |
| `do_precip_albedo` | logical | `.false.` | increase the surface albedo where it rains (a crude cloud effect) |
| `precip_albedo_mode` | character(len=14) | `'full'` | if so, use total (`'full'`), large-scale (`'lscale'`) or convective (`'conv'`) precipitation |
| `precip_albedo` | real(kind=rb) | `0.35` | if so, the albedo of a fully precipitating grid box |
| `precip_lat` | real(kind=rb) | `0.0` | [deg] if so, apply only poleward of this latitude |

<!-- mimadoc:end -->

The other RRTMG inputs are not namelist variables. MiMA calls RRTMG without clouds and aerosols:

 RRTMG input | Value
 :--- | :---
 `icld`, `iaer` | 0: no clouds, no aerosols (cloud and aerosol optical properties are zero)
 `h2ovmr` | from the specific humidity (at least `h2o_lower_limit`), every radiation step
 `o3vmr` | from `ozone_file`, or `o3_val`
 `co2vmr` | `co2ppmv`
 `ch4vmr`, `n2ovmr`, `o2vmr`, `cfc11vmr`, `cfc12vmr`, `cfc22vmr`, `ccl4vmr` | the `*_val` values if `include_secondary_gases`, else 0
 `asdir`, `asdif`, `aldir`, `aldif` | the surface albedo from [`surface_albedo_mod`](https://eddy-stanford.github.io/MiMA/api/surface_albedo_mod/), with the snow on bucket land (modified by `do_precip_albedo`), the same in all four
 `emis` | 1 (black-body surface)
 `coszen` | computed from [`astro_nml`](#astro_nml), every radiation step
 `adjes` | `solrad`
 `dyofyr` | day of the year if `use_dyofyr`, else 0
 `scon` | `solr_cnst`

### `astro_nml`

The orbit and solar constant for RRTM.

<!-- mimadoc:namelist astro_nml -->

Declared in [`rrtm_astro`](https://eddy-stanford.github.io/MiMA/api/rrtm_astro/) (`src/atmos_param/radiation/rrtm/astro.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `obliq` | real(kind=rb) | `23.439` | [deg] obliquity |
| `use_dyofyr` | logical | `.false.` | let RRTM compute the Earth-Sun distance from the day of the year (assumes 365 days per year) |
| `solr_cnst` | real(kind=rb) | `1370.` | [W/m2] solar constant |
| `solrad` | real(kind=rb) | `1.0` | Earth-Sun distance factor if `use_dyofyr = .false.` |
| `solday` | integer(kind=im) | `0` | if > 0, perpetual run at this day of the year |
| `equinox_day` | real(kind=rb) | `0.25` | fraction of the year at which the March equinox occurs |

<!-- mimadoc:end -->

### `gray_radiation_nml`

The gray radiation of [Frierson, Held and Zurita-Gotor (2006)](https://doi.org/10.1175/JAS3753.1). The insolation is annual-mean and zonally symmetric, `solar_constant/4 * (1 + del_sol*P2(lat) + del_sw*sin(lat))`.

<!-- mimadoc:namelist gray_radiation_nml -->

Declared in [`gray_radiation_mod`](https://eddy-stanford.github.io/MiMA/api/gray_radiation_mod/) (`src/atmos_param/radiation/gray_radiation.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `solar_constant` | real | `1360.0` | [W/m2] solar constant |
| `del_sol` | real | `0.0` | amplitude of the second Legendre polynomial in the insolation (equator-to-pole contrast) |
| `ir_tau_eq` | real | `4.0` | longwave optical depth at the surface at the equator |
| `ir_tau_pole` | real | `4.0` | longwave optical depth at the surface at the poles |
| `atm_abs` | real | `0.2` | shortwave optical depth of the atmosphere |
| `sw_diff` | real | `0.0` | equator-to-pole reduction of the shortwave optical depth |
| `long_pert` | real | `180.` | [deg] longitude of the zonally localized insolation perturbation |
| `del_long` | real | `30.` | [deg] width of the zonally localized insolation perturbation |
| `size_pert` | real | `0.` | [W/m2] amplitude of a zonally localized insolation perturbation |
| `linear_tau` | real | `0.1` | fraction of the longwave optical depth that is linear in pressure (the rest goes as p<sup>4</sup>) |
| `del_sw` | real | `0.0` | amplitude of a hemispherically asymmetric insolation term (winter/summer hemisphere) |
| `lat_pert` | real | `0.0` | [deg] latitude of the centre of the localized (Walker-type) insolation perturbation |
| `lon_pert` | real | `180.0` | [deg] longitude of the centre of the localized (Walker-type) insolation perturbation |
| `del_lat` | real | `30.0` | [deg] latitudinal half-width of the localized (Walker-type) insolation perturbation |
| `del_lon` | real | `90.0` | [deg] longitudinal half-width of the localized (Walker-type) insolation perturbation |
| `fcng_pert` | real | `0.0` | [W/m2] amplitude of a localized (Walker-type) insolation perturbation |

<!-- mimadoc:end -->

## Surface

The surface consists of an ocean model, a land model, a surface albedo and a flux exchange, each with its own namelist. See [Surface options](Configurations.md#surface-options) for how to combine them, and the [migration table](Migration_v2.md#old-to-new-namelist-mapping) for the old `&simple_surface_nml` variables.

### `flux_exchange_nml`

The surface flux computation and the reference heights of the diagnostics.

<!-- mimadoc:namelist flux_exchange_nml -->

Declared in [`flux_exchange_mod`](https://eddy-stanford.github.io/MiMA/api/flux_exchange_mod/) (`src/surface/flux_exchange.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `z_ref_heat` | real | `2.` | [m] reference height of the diagnostics `t_ref` and `rh_ref` |
| `z_ref_mom` | real | `10.` | [m] reference height of the diagnostics `u_ref` and `v_ref` |

<!-- mimadoc:end -->

### `surface_albedo_nml`

The background surface albedo, the same on land and ocean. With `albedo_profile = 'step'` the albedo is high where the latitude is poleward of `step_lat_nh` or `step_lat_sh` (strict comparison, degrees). The desert boxes (Sahara, Gobi, Australia) exist only with `'tanh_desert'`. Snow on bucket land brightens the albedo of the land points, and the radiation sees it (see [`land_model_nml`](#land_model_nml)).

<!-- mimadoc:namelist surface_albedo_nml -->

Declared in [`surface_albedo_mod`](https://eddy-stanford.github.io/MiMA/api/surface_albedo_mod/) (`src/surface/surface_albedo.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `albedo_profile` | character(len=32) | `'tanh_desert'` | latitude profile of the albedo: `'constant'`: `const_albedo`; `'step'`: `higher_albedo` where the latitude is above `step_lat_nh` or below `step_lat_sh`; `'power'`: increase as `(\|lat\|/90)^albedo_exp`; `'tanh'`: tanh increase centred at `albedo_cntrNH`, `-albedo_cntrSH` with width `albedo_wdth`; `'sin2'`: `sin^2(lat)` increase from equator to pole; `'tanh_desert'`: as `'tanh'`, plus `albedo_desert` over the Sahara, Gobi and Australian deserts |
| `const_albedo` | real | `0.23` | surface albedo (the low-latitude value of the profiles) |
| `higher_albedo` | real | `0.80` | high-latitude albedo of the profiles |
| `albedo_exp` | real | `2.` | exponent of the `'power'` profile |
| `albedo_cntrNH` | real | `68.` | [deg] centre latitude of the Northern Hemisphere albedo increase (`'tanh'`, `'tanh_desert'`) |
| `albedo_cntrSH` | real | `64.` | [deg] centre latitude (positive) of the Southern Hemisphere albedo increase (`'tanh'`, `'tanh_desert'`) |
| `albedo_wdth` | real | `5.` | [deg] width of the albedo increase (`'tanh'`, `'tanh_desert'`) |
| `step_lat_nh` | real | `90.` | [deg] with `'step'`, `higher_albedo` where the latitude is above this value (90: none) |
| `step_lat_sh` | real | `-90.` | [deg] with `'step'`, `higher_albedo` where the latitude is below this value (-90: none) |
| `albedo_desert` | real | `0.20` | with `'tanh_desert'`, the albedo over the desert boxes is `const_albedo + albedo_desert` |

<!-- mimadoc:end -->

### `ocean_model_nml`

The ocean: a slab (`'slab'`), a prescribed SST (`'prescribed'`) or a fixed temperature (`'fixed'`), its initial SST, its heat capacity and its roughness length. `sst_file` is read from `INPUT/<sst_file>.nc`; see [Prescribed SST with interactive land](Configurations.md#prescribed-sst-with-interactive-land) for the file conventions.

<!-- mimadoc:namelist ocean_model_nml -->

Declared in [`ocean_model_mod`](https://eddy-stanford.github.io/MiMA/api/ocean_model_mod/) (`src/surface/ocean_model.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `ocean_model` | character(len=32) | `'slab'` | `'slab'`: mixed layer with an interactive SST; `'prescribed'`: SST read from `sst_file` at every step (implies `sst_init = 'file'`); `'fixed'`: SST held at its initial value |
| `sst_init` | character(len=32) | `'profile'` | initial SST on a cold start: `'profile'`: `Tm - deltaT*(3 sin^2(lat) - 1)/3`; `'uniform'`: `Tm` everywhere; `'ape_control'`, `'ape_qobs'`: the Aqua-Planet Experiment profiles "control" and "qobs" (Neale and Hoskins, 2000); `'file'`: from `sst_file` |
| `Tm` | real | `285.` | [K] mean initial SST of `'profile'`; the SST of `'uniform'` |
| `sst_offset` | real | `0.` | [K] added to the cold-start SST of every `sst_init` (including `'file'`), and at every step to the SST of a prescribed ocean (`ocean_model = 'prescribed'`; e.g. for AMIP+4K experiments); a fixed ocean (`ocean_model = 'fixed'`) gets it through its initial SST. For the APE profiles it is added after the profile: `'ape_qobs'` with `sst_offset = 5.` is the "qobs+5" experiment |
| `deltaT` | real | `40.` | [K] equator-to-pole difference of the initial SST of `'profile'` |
| `sst_file` | character(len=256) | `''` | SST file name, without `.nc`, in `INPUT/` |
| `sst_file_var` | character(len=256) | `''` | name of the SST variable in `sst_file` (empty: `sst_file`) |
| `rough_ocean` | real | `3.21e-05` | [m] roughness length of the ocean (momentum, heat and moisture) |
| `heat_capacity` | real | `3.e08` | [J/m2/K] mixed-layer heat capacity (poleward of `heat_cap_limit`); must be positive |
| `trop_capacity` | real | `1.e08` | [J/m2/K] heat capacity equatorward of `trop_cap_limit`, varying linearly to `heat_capacity` at `heat_cap_limit`; <= 0: `heat_capacity` |
| `trop_cap_limit` | real | `20.` | [deg] latitude of the tropical heat capacity `trop_capacity` |
| `heat_cap_limit` | real | `60.` | [deg] latitude of the extratropical heat capacity `heat_capacity` |
| `np_cap_factor` | real | `1.` | factor on `heat_capacity` in the Northern Hemisphere |

<!-- mimadoc:end -->

### `qflux_nml`

Prescribed ocean heat flux convergence (Q-flux) of the slab ocean: analytic, with the zonally asymmetric fluxes of [Garfinkel et al. (2020)](https://doi.org/10.1175/JCLI-D-19-0181.1), or read from a file ([Q-flux from a file](Configurations.md#q-flux-from-a-file)). With the defaults the Q-flux is zero over land and its area-weighted ocean mean is removed; both switches `.false.` give the Q-flux of MiMA v1.

<!-- mimadoc:namelist qflux_nml -->

Declared in [`ocean_qflux_mod`](https://eddy-stanford.github.io/MiMA/api/ocean_qflux_mod/) (`src/surface/ocean_qflux.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `qflux_source` | character(len=32) | `'analytic'` | `'analytic'`: the meridional Q-flux of Merlis et al. (2013) (`qflux_amp`, `qflux_width`); `'none'`: no meridional Q-flux; `'file'`: read from `qflux_file` and interpolated in time at every step (1 record: constant; 2 to 12 records: a climatology repeated every year; more: a dated time series) |
| `do_warmpool` | logical | `.true.` | add the zonally asymmetric Q-flux `warmpool_pattern`, whatever `qflux_source` |
| `qflux_file` | character(len=256) | `''` | with `qflux_source = 'file'`, Q-flux file name, without `.nc`, in `INPUT/` |
| `qflux_file_var` | character(len=256) | `'qflux'` | name of the Q-flux variable [W/m2] in `qflux_file` (empty: `qflux_file`); for the output of a prescribed-SST run, `'qflux_implied'` |
| `qflux_ocean_only` | logical | `.true.` | apply the Q-flux over the ocean only (the default); `.false.`: also over land, as in v1 (set both `qflux_ocean_only` and `qflux_remove_mean` to `.false.` to reproduce v1) |
| `qflux_remove_mean` | logical | `.true.` | subtract a constant from the Q-flux at the ocean points so that its global area integral is zero (the default): the constant is the area integral of the applied Q-flux over all points (land included with `qflux_ocean_only = .false.`) divided by the ocean area; with `qflux_ocean_only` it is the ocean mean (the global mean without land). `.false.`: no mean removed, as in v1 (set both `qflux_ocean_only` and `qflux_remove_mean` to `.false.` to reproduce v1) |
| `qflux_amp` | real | `26.` | [W/m2] amplitude of the meridional Q-flux; with `warmpool_pattern = 'garfinkel2020'` it also scales several regional terms (Caribbean and South America, the West African coast, the Barents Sea and, with `warmpool_amp`, Australia and the Indian Ocean), even with `qflux_source = 'none'` |
| `qflux_width` | real | `16.` | [deg] half-width of the meridional Q-flux |
| `warmpool_pattern` | character(len=32) | `'garfinkel2020'` | the zonally asymmetric Q-flux added with `do_warmpool`: `'garfinkel2020'`: the warm pool and the regional patterns of Garfinkel et al. (2020) (v1's `warmpool_localization_choice = 3`); `'warm_pool'`: the warm pool alone, one wavelength in longitude (v1's choice 1, which went round the whole latitude band). v1's choice 2 was removed. See `ocean_qflux.f90`. |
| `warmpool_amp` | real | `18.` | [W/m2] amplitude of the warm pool; with `'garfinkel2020'`, also of its wave-8 eastern edge, and the West Pacific, Pacific, Hawaii and Australian terms are only applied if `warmpool_amp` > 0.001 (the Australian ones scale with `qflux_amp + warmpool_amp`) |
| `warmpool_width` | real | `35.` | [deg] latitudinal half-width of the warm pool (and, with `'garfinkel2020'`, of its eastern edge and of the Africa wave) |
| `warmpool_centr` | real | `0.` | [deg] central latitude of the warm pool (and, with `'garfinkel2020'`, of its eastern edge and of the Africa wave) |
| `warmpool_k` | real | `1.66666` | zonal wave number of the warm pool; with `'warm_pool'` it sets the window (one wavelength); `'garfinkel2020'` uses the fixed window [`warmpool_phase` - 54, `warmpool_phase` + 162], which assumes `warmpool_k` = 5/3 |
| `warmpool_phase` | real | `140.` | [deg] longitude of the warm-pool maximum; with `'garfinkel2020'` it also moves the longitude windows (not the centres) of the warm-pool edge, Africa, Kuroshio and Gulf Stream terms |
| `gulf_k` | integer | `4` | with `'garfinkel2020'`, the zonal wave number plus 1 of the heating near 67N between 318E and 18E (the default gives wave number 3) |
| `gulf_phase` | real | `310.` | [deg] with `'garfinkel2020'`, the phase of the heating near 67N: its maximum is at `gulf_phase` + 38 |
| `gulf_amp` | real | `70.` | [W/m2] with `'garfinkel2020'`, the Gulf Stream amplitude; the Atlantic terms (tropical Atlantic, Caribbean and South America, near 67N, West African coast, Barents Sea) are only applied if `gulf_amp` > 0.001 and `trop_atlantic_amp` > 0.001 |
| `kuroshio_amp` | real | `40.` | [W/m2] with `'garfinkel2020'`, the Kuroshio amplitude; with `kuroshio_amp` <= 0.001, a reduced set of West Pacific terms is used |
| `trop_atlantic_amp` | real | `50.` | [W/m2] with `'garfinkel2020'`, the tropical Atlantic amplitude, also of the heating near 67N; half of it scales the Africa wave and switches the Agulhas terms on (if > 0.001); the Arctic terms are only applied if it is > 0.001 |
| `Pac_ITCZextra` | real | `0.` | [W/m2] extra cooling in the subtropical South-East Pacific and zonal heating near 55S (`'garfinkel2020'` only) |
| `Sampeextra` | real | `0.` | [W/m2] zonal dipole near 40S and 48S (`'garfinkel2020'` only) |
| `Pac_SPCZextra` | real | `0.` | [W/m2] extra flux in the SPCZ and South Pacific, which modulates the SPCZ (`'garfinkel2020'` only) |
| `Africaextra` | real | `0.` | [W/m2] extra flux near the Agulhas current (`'garfinkel2020'` only) |
| `Hawaiiextra` | real | `30.0` | [W/m2] heating in the western North Pacific near 16N, which also changes many Pacific and Indian Ocean terms (`'garfinkel2020'` only; see the table) |

<!-- mimadoc:end -->

### `land_model_nml`

The land: where it is (`land_mask`), what it is (`land_model`: slab, prescribed or fixed temperature), its heat capacity, roughness and, optionally, a bucket hydrology and snow. See [Land-sea mask](Configurations.md#land-sea-mask) and [Bucket land and snow](Configurations.md#bucket-land-and-snow). The bucket and snow are verified for conservation but not climate-validated. `hydrology = 'bucket'` requires `rough_moist_factor = 1.` and `rough_moist_profile = 'uniform'`, and `do_snow` requires the bucket and `land_model = 'simple'`.

<!-- mimadoc:namelist land_model_nml -->

Declared in [`land_model_mod`](https://eddy-stanford.github.io/MiMA/api/land_model_mod/) (`src/surface/land_model.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `land_model` | character(len=32) | `'simple'` | `'simple'`: slab land surface with an interactive temperature; `'prescribed'`: temperature read from `t_surf_file` at every step; `'fixed'`: temperature held at its initial value |
| `land_mask` | character(len=32) | `'navy'` | where the land is: `'none'`; `'navy'`: Navy land-sea mask (the `water_file` of `topography_nml`); `'file'`: land where the variable `land_mask_file` of `INPUT/<land_mask_file>.nc` (on the model grid) is positive; `'boxes'`: the boxes `slandlon`..`elandlon`, `slandlat`..`elandlat`; `'height'`: surface height above `land_height_limit` |
| `land_mask_file` | character(len=256) | `'lmask'` | with `land_mask = 'file'`, name of the land-sea mask file (without `.nc`, in `INPUT/`) and of its variable |
| `slandlon`, `slandlat`, `elandlon`, `elandlat` | real, dimension(10) | `0`, `0`, `-1`, `-1` | with `land_mask = 'boxes'`, start and end longitude and latitude [deg] of up to 10 land boxes |
| `land_height_limit` | real | `10.` | [m] with `land_mask = 'height'`, points higher than this are land |
| `heat_capacity` | real | `1.e07` | [J/m2/K] heat capacity of the land; <= 0: the ocean `heat_capacity` |
| `t_init` | character(len=32) | `'ocean'` | initial land temperature on a cold start: `'ocean'`: the cold-start SST of `ocean_model_nml` (`sst_init`); `'file'`: from `t_surf_file`. `land_model = 'prescribed'` always starts from `t_surf_file` |
| `t_surf_file` | character(len=256) | `''` | land temperature file name, without `.nc`, in `INPUT/` |
| `t_surf_file_var` | character(len=256) | `''` | name of the temperature variable in `t_surf_file` (empty: `t_surf_file`) |
| `rough_mom_factor` | real | `5.e3` | factor on the ocean roughness `rough_ocean` for the land momentum roughness |
| `rough_moist_factor` | real | `1.e-12` | factor on `rough_ocean` for the land moisture roughness |
| `rough_moist_profile` | character(len=32) | `'legacy_tropical'` | land moisture roughness: `'uniform'`: `rough_ocean*rough_moist_factor`; `'legacy_tropical'`: as `'uniform'`, plus a larger value over tropical land (legacy) |
| `hydrology` | character(len=32) | `'saturated'` | land hydrology: `'saturated'`: the land surface is saturated, as the ocean; `'bucket'`: Manabe bucket, evaporation limited by the soil water (needs `rough_moist_factor = 1.` and `rough_moist_profile = 'uniform'`) |
| `water_capacity` | real | `150.` | [kg/m2] with `hydrology = 'bucket'`, field capacity of the bucket; more water runs off |
| `beta_frac` | real | `0.75` | with `hydrology = 'bucket'`, fraction of `water_capacity` above which the land evaporates at the potential rate (below it, in proportion to the water) |
| `water_init_frac` | real | `1.0` | with `hydrology = 'bucket'`, initial soil water on a cold start, as a fraction of `water_capacity` |
| `do_snow` | logical | `.false.` | snow on land: stored snowfall, sublimation, frost, melt and snow albedo (needs `hydrology = 'bucket'` and `land_model = 'simple'`) |
| `albedo_snow` | real | `0.60` | with `do_snow`, albedo of deep snow (snow never lowers the albedo below the background) |
| `snow_mask_depth` | real | `10.` | [kg/m2] with `do_snow`, snow at which the albedo is halfway between the background and `albedo_snow` |
| `snow_max` | real | `1000.` | [kg/m2] with `do_snow`, maximum snow; more runs off (a crude glacier) |

<!-- mimadoc:end -->

### `surface_flux_nml`

Surface fluxes.

<!-- mimadoc:namelist surface_flux_nml -->

Declared in [`surface_flux_mod`](https://eddy-stanford.github.io/MiMA/api/surface_flux_mod/) (`src/surface/surface_flux.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `no_neg_q` | logical | `.false.` | set negative lowest-level humidity to zero |
| `use_virtual_temp` | logical | `.false.` | use virtual potential temperature for the surface-layer stability |
| `alt_gustiness` | logical | `.false.` | alternative gustiness: a lower bound `gust_const` on the wind speed |
| `gust_const` | real | `1.0` | [m/s] see `alt_gustiness` |
| `old_dtaudv` | logical | `.true.` | use the same d(stress)/d(wind) for both wind components |
| `use_mixing_ratio` | logical | `.false.` | Manabe Climate Model form of the moisture flux (legacy) |
| `ncar_ocean_flux` | logical | `.false.` | NCAR (Large and Yeager) ocean flux formulation |
| `no_surface_momentum_flux` | logical | `.false.` | switch off the surface momentum flux |
| `no_surface_moisture_flux` | logical | `.false.` | switch off the surface moisture flux |
| `no_surface_heat_flux` | logical | `.false.` | switch off the surface sensible heat flux |
| `no_surface_radiative_flux` | logical | `.false.` | switch off the surface radiative flux |

<!-- mimadoc:end -->

### `monin_obukhov_nml`

Monin-Obukhov similarity for the surface layer.

<!-- mimadoc:namelist monin_obukhov_nml -->

Declared in [`mima_monin_obukhov_mod`](https://eddy-stanford.github.io/MiMA/api/mima_monin_obukhov_mod/) (`src/atmos_param/monin_obukhov/monin_obukhov.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `rich_crit` | real | `2.0` | critical Richardson number (must be > 0.25) |
| `neutral` | logical | `.false.` | neutral stability everywhere |
| `drag_min` | real | `4.e-05` | minimum drag coefficient |
| `stable_option` | integer | `1` | stability function for stable conditions (1 or 2); used only by `stable_mix`, which is not called in MiMA |
| `zeta_trans` | real | `0.5` | transition value of z/L for `stable_option = 2` |

<!-- mimadoc:end -->

## Boundary layer

### `vert_turb_driver_nml`

Boundary-layer turbulence. The only scheme is the non-local K-profile scheme of [`diffusivity_nml`](#diffusivity_nml).

<!-- mimadoc:namelist vert_turb_driver_nml -->

Declared in [`vert_turb_driver_mod`](https://eddy-stanford.github.io/MiMA/api/vert_turb_driver_mod/) (`src/atmos_param/vert_turb_driver/vert_turb_driver.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `gust_scheme` | character(len=24) | `'constant'` | surface gustiness: `'constant'` (`constant_gust`) or `'beljaars'` (from u* and b*) |
| `constant_gust` | real | `0.` | [m/s] constant gustiness |
| `use_tau` | logical | `.false.` | use the current time level (`.true.`) or the updated values (`.false.`) |
| `do_molecular_diffusion` | logical | `.false.` | add molecular diffusion |
| `do_diffusivity` | logical | `.true.` | compute diffusion coefficients with the non-local K scheme (`.false.`: no boundary-layer diffusion) |
| `gust_factor` | real | `1.0` | factor for the `'beljaars'` gustiness |

<!-- mimadoc:end -->

### `diffusivity_nml`

The non-local K scheme.

<!-- mimadoc:namelist diffusivity_nml -->

Declared in [`diffusivity_mod`](https://eddy-stanford.github.io/MiMA/api/diffusivity_mod/) (`src/atmos_param/diffusivity/diffusivity.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `fixed_depth` | logical | `.false.` | use a fixed boundary-layer depth `depth_0` |
| `depth_0` | real | `5000.0` | [m] boundary-layer depth if `fixed_depth` |
| `frac_inner` | real | `0.1` | fraction of the boundary layer that is the surface layer (between 0 and 1) |
| `rich_crit_pbl` | real | `1.0` | critical bulk Richardson number defining the boundary-layer top |
| `background_m` | real | `0.0` | [m2/s] minimum diffusivity for momentum |
| `background_t` | real | `0.0` | [m2/s] minimum diffusivity for heat |

<!-- mimadoc:end -->

### `vert_diff_driver_nml`

<!-- mimadoc:namelist vert_diff_driver_nml -->

Declared in [`vert_diff_driver_mod`](https://eddy-stanford.github.io/MiMA/api/vert_diff_driver_mod/) (`src/atmos_param/vert_diff_driver/vert_diff_driver.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `do_conserve_energy` | logical | `.true.` | heat the air by the dissipation of kinetic energy in the vertical diffusion |
| `use_virtual_temp_vert_diff` | logical | `.false.` | use virtual temperature in the vertical diffusion |

<!-- mimadoc:end -->

## Moisture

Following [Frierson (2007)](https://doi.org/10.1175/JAS3935.1), MiMA uses large-scale condensation with the Betts-Miller convection scheme and its "shallower" shallow convection.

### `moist_processes_nml`

**Precipitation phase.** With `precip_phase = 'temperature'` (the default) the precipitation of a column falls as snow where the mean temperature of its lowest `pdepth` (150 hPa) is at or below `tfreeze`, and as rain elsewhere (the GFDL test, reused unchanged). The schemes then release the latent heat of sublimation, `hls = hlv + hlf`, for the snow, so that the energy balances when the surface takes `hlf` to melt it: `lscale_cond` and `moist_conv` use `hls` in the cold columns (condensation, and the re-evaporation of falling snow in `lscale_cond`), and Betts-Miller adds `hlf` times the precipitation as a uniform temperature change over its convecting layers, on top of its enthalpy correction. `precip_phase = 'liquid'` makes all precipitation rain, as in v1 (bit for bit). The frozen precipitation is the `moist/snow_conv` and `moist/snow_ls` diagnostics, and it feeds the snow of [`land_model_nml`](#land_model_nml). Known limits: the saturation vapour pressure is still over liquid water everywhere (the simple table of `sat_vapor_pres_nml`), and the Betts-Miller parcel and CAPE calculation use `hlv`. See [Changes to the answers](Migration_v2.md#changes-to-the-answers) for the impact.

<!-- mimadoc:namelist moist_processes_nml -->

Declared in [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) (`src/atmos_param/moist_processes/moist_processes.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `do_mca`, `do_lsc` | logical | `.false.`, `.true.` | `do_mca`: moist convective adjustment (`moist_conv_nml`); `do_lsc`: large-scale condensation (`lscale_cond_nml`) |
| `precip_phase` | character(len=32) | `'temperature'` | phase of the precipitation: `'temperature'`: snow in the columns whose mean temperature over the lowest `pdepth` is at or below `tfreeze`, rain elsewhere; `'liquid'`: all rain (v1) |
| `pdepth` | real | `150.e2` | [Pa] with `precip_phase = 'temperature'`, depth of the lowest layer whose mean temperature decides between rain and snow |
| `tfreeze` | real | `273.16` | [K] with `precip_phase = 'temperature'`, freezing temperature for that decision |
| `use_tau`, `do_gust_cv` | logical | `.false.`, `.false.` | `use_tau`: use the current time level (`.true.`) or the updated values (`.false.`); `do_gust_cv`: convective gustiness |
| `gustmax` | real | `3.` | [m/s] maximum convective gustiness |
| `gustconst` | real | `10./86400.` | [kg/m2/s] precipitation rate at which convective gustiness starts to matter |
| `do_bm` | logical | `.true.` | Betts-Miller convection (`betts_miller_nml`); not with `do_mca` |

<!-- mimadoc:end -->

### `betts_miller_nml`

<!-- mimadoc:namelist betts_miller_nml -->

Declared in [`betts_miller_mod`](https://eddy-stanford.github.io/MiMA/api/betts_miller_mod/) (`src/atmos_param/betts_miller/betts_miller.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `tau_bm` | real | `7200.` | [s] relaxation time |
| `rhbm` | real | `.7` | relative humidity of the reference profile |
| `do_simp` | logical | `.false.` | adjust the time scales so that precipitation is always continuous |
| `do_shallower` | logical | `.true.` | shallow convection: choose a smaller depth so that precipitation is zero |
| `do_changeqref` | logical | `.false.` | shallow convection: change both q and T so that precipitation is zero |
| `do_envsat` | logical | `.false.` | reference humidity relative to the environment (`.true.`) or the parcel (`.false.`) |
| `do_taucape` | logical | `.false.` | make `tau_bm` proportional to `CAPE**(-1/2)` |
| `capetaubm` | real | `900.` | [J/kg] CAPE at which the relaxation time is `tau_bm` (with `do_taucape`) |
| `tau_min` | real | `2400.` | [s] minimum relaxation time (with `do_taucape`) |
| `buoyancy_kick` | real | `0.` | [K] temperature added to the parcel at the lowest level |

<!-- mimadoc:end -->

### `lscale_cond_nml`

<!-- mimadoc:namelist lscale_cond_nml -->

Declared in [`lscale_cond_mod`](https://eddy-stanford.github.io/MiMA/api/lscale_cond_mod/) (`src/atmos_param/lscale_cond/lscale_cond.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `hc` | real | `1.00` | relative humidity at which condensation occurs (0 <= `hc` <= 1) |
| `do_evap` | logical | `.true.` | re-evaporate falling precipitation in sub-saturated layers below |

<!-- mimadoc:end -->

### `moist_conv_nml`

Moist convective adjustment (used with `do_mca = .true.`).

<!-- mimadoc:namelist moist_conv_nml -->

Declared in [`moist_conv_mod`](https://eddy-stanford.github.io/MiMA/api/moist_conv_mod/) (`src/atmos_param/moist_conv/moist_conv.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `HC` | real | `1.00` | relative humidity at which the adjustment occurs |
| `TOLmin`, `TOLmax` | real | `.02`, `.10` | [K] `TOLmin`, `TOLmax`: initial and maximum tolerance of the iterative adjustment |
| `ITSMOD` | integer | `30` | number of iterations at each tolerance |

<!-- mimadoc:end -->

## Damping and gravity-wave drag

### `damping_driver_nml`

Upper boundary and gravity-wave drag, used with `do_damping = .true.` in [`physics_driver_nml`](#physics_driver_nml).

<!-- mimadoc:namelist damping_driver_nml -->

Declared in [`damping_driver_mod`](https://eddy-stanford.github.io/MiMA/api/damping_driver_mod/) (`src/atmos_param/damping_driver/damping_driver.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `trayfric` | real | `-0.5` | Rayleigh friction time scale: [s] if > 0, [days] if < 0 |
| `do_rayleigh` | logical | `.false.` | Rayleigh friction (sponge) at the top of the model |
| `sponge_pbottom` | real | `50.` | [Pa] bottom of the Rayleigh sponge |
| `do_cg_drag` | logical | `.true.` | non-orographic (convective) gravity-wave drag (`cg_drag_nml`) |
| `do_mg_drag` | logical | `.false.` | orographic gravity-wave drag (`mg_drag_nml`) |
| `do_conserve_energy` | logical | `.true.` | heat the air by the momentum lost to the damping |
| `do_const_drag` | logical | `.false.` | idealized seasonal "gravity-wave" drag in the stratosphere |
| `const_drag_amp` | real | `3.e-04` | [m/s2] amplitude of the constant drag |
| `const_drag_off` | real | `0.` | offset of its latitudinal profile |

<!-- mimadoc:end -->

### `cg_drag_nml`

The Alexander and Dunkerton (1999) non-orographic gravity-wave scheme, with the changes of [Cohen et al. (2013)](https://doi.org/10.1175/JAS-D-12-0240.1) and [Garfinkel et al. (2020)](https://doi.org/10.1175/JCLI-D-19-0181.1).

<!-- mimadoc:namelist cg_drag_nml -->

Declared in [`cg_drag_mod`](https://eddy-stanford.github.io/MiMA/api/cg_drag_mod/) (`src/atmos_param/cg_drag/cg_drag.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `cg_drag_freq` | integer | `21600` | [s] interval between drag calculations; the drag is held fixed in between (and saved in `cg_drag.res.nc`) |
| `cg_drag_offset` | integer | `0` | [s] on a cold start, time to the first calculation (0: `cg_drag_freq`) |
| `source_level_pressure` | real | `315.e+02` | [Pa] the source level is the highest level with a pressure greater than this at the equator |
| `damp_level_pressure` | real | `0.85e+02` | [Pa] momentum flux reaching the model top is deposited from the top down to this level |
| `nk` | integer | `1` | number of wavelengths in the spectrum |
| `cmax` | real | `99.6` | [m/s] maximum phase speed |
| `dc` | real | `1.2` | [m/s] phase-speed resolution |
| `Bt_0` | real | `0.0043` | [Pa] total source momentum flux poleward of `phi0n`/`phi0s` |
| `Bt_sh` | real | `0.00` | [Pa] additional flux in the Southern Hemisphere extratropics (tanh transition of width `dphis` at `phi0s`) |
| `Bt_nh` | real | `0.00` | [Pa] additional flux in the Northern Hemisphere extratropics (tanh transition of width `dphin` at `phi0n`) |
| `Bt_eq` | real | `0.0043` | [Pa] total source momentum flux between `dphis` and `dphin`; it varies linearly to `Bt_0` at `phi0n`/`phi0s` |
| `phi0n`, `phi0s`, `dphin`, `dphis` | real | `15.`, `-15.`, `10.`, `-10.` | [deg] `phi0n`, `phi0s`: latitudes where the flux reaches `Bt_0`. `dphin`, `dphis`: edges of the tropical band (flux `Bt_eq`, spectrum width `cwtropics`), and widths of the `Bt_nh`/`Bt_sh` transitions |
| `Bw` | real | `0.4` | [m2/s2] amplitude of the wide part of the phase-speed spectrum |
| `Bn` | real | `0.0` | [m2/s2] amplitude of the narrow part of the phase-speed spectrum (0 in the tropical band) |
| `cw` | real | `35.0` | [m/s] half-width of the wide spectrum outside the tropical band |
| `cwtropics` | real | `35.0` | [m/s] half-width of the wide spectrum inside the tropical band |
| `cn` | real | `2.0` | [m/s] half-width of the narrow spectrum |
| `flag` | integer | `0` | `1`: spectrum peaks at c = 0; `0`: at c - u = 0 (always 0 in the tropical band) |
| `kelvin_kludge` | real | `1.` | factor on the source flux of waves with c - u < 0 in the tropical band |

<!-- mimadoc:end -->

### `mg_drag_nml`

Orographic gravity-wave drag, used with `do_mg_drag = .true.`.

<!-- mimadoc:namelist mg_drag_nml -->

Declared in [`mg_drag_mod`](https://eddy-stanford.github.io/MiMA/api/mg_drag_mod/) (`src/atmos_param/mg_drag/mg_drag.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `xl_mtn` | real | `1.0e5` | [m] effective mountain length |
| `gmax`, `acoef`, `rho` | real | `2.0`, `1.0`, `1.13` | `gmax`: order-one tuning parameter (larger: more drag); `acoef`: order-one tuning parameter; `rho` [kg/m3]: standard sea-level air density |
| `low_lev_frac` | real | `.23` | fraction of the atmosphere (from the bottom) used for the base flux, where no wave breaking is allowed |
| `do_conserve_energy` | logical | `.false.` | heat the air by the dissipated kinetic energy |
| `source_of_sgsmtn` | character(len=128) | `'input'` | sub-grid orography: `'input'` (read from `INPUT/mg_drag.res.nc`) or `'computed'` (from the high-resolution topography) |
| `flux_cut_level` | real | `0.0` | [Pa] above this level the flux divergence is set to zero |

<!-- mimadoc:end -->

## Held-Suarez forcing

### `held_suarez_nml`

The [Held and Suarez (1994)](https://doi.org/10.1175/1520-0477(1994)075<1825:APFTIO>2.0.CO;2) forcing, used with `do_held_suarez = .true.`; the defaults are the HS94 values. See [Held-Suarez forcing](Configurations.md#held-suarez-forcing).

<!-- mimadoc:namelist held_suarez_nml -->

Declared in [`held_suarez_mod`](https://eddy-stanford.github.io/MiMA/api/held_suarez_mod/) (`src/atmos_param/held_suarez/held_suarez.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `t_zero` | real | `315.` | [K] surface equilibrium temperature at the equator |
| `t_strat` | real | `200.` | [K] minimum (stratospheric) equilibrium temperature |
| `delh` | real | `60.` | [K] equator-to-pole temperature difference |
| `delv` | real | `10.` | [K] vertical potential temperature difference |
| `p_ref` | real | `1.e5` | [Pa] reference pressure |
| `sigma_b` | real | `0.7` | top of the frictional boundary layer (sigma) |
| `ka` | real | `40.` | [days] free-atmosphere relaxation time |
| `ks` | real | `4.` | [days] surface relaxation time at the equator |
| `kf` | real | `1.` | [days] boundary-layer Rayleigh friction time |
| `do_rayleigh_friction` | logical | `.true.` | apply the boundary-layer friction |
| `do_conserve_energy` | logical | `.false.` | heat the air by the frictional dissipation |

<!-- mimadoc:end -->

## Local heating

### `local_heating_nml`

Prescribed Gaussian heating, used with `do_local_heating = .true.`. Every variable is an array of up to 10 entries, one per heat source.

<!-- mimadoc:namelist local_heating_nml -->

Declared in [`local_heating_mod`](https://eddy-stanford.github.io/MiMA/api/local_heating_mod/) (`src/atmos_param/local_heating/local_heating.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `hamp` | real, dimension(ngauss) | `0.` | [K/day] amplitude of the heating |
| `lonwidth` | real, dimension(ngauss) | `-1.` | [deg] zonal width, if `loncenter` >= 0 |
| `loncenter` | real, dimension(ngauss) | `-1.` | [deg] longitude of the centre; zonally symmetric if < 0 |
| `lonmove` | real, dimension(ngauss) | `0.` | [deg/day] zonal speed of the source |
| `latwidth` | real, dimension(ngauss) | `15.` | [deg] meridional width |
| `latcenter` | real, dimension(ngauss) | `0.` | [deg] latitude of the centre |
| `latmove` | real, dimension(ngauss) | `0.` | [deg/day] meridional speed of the source |
| `pwidth` | real, dimension(ngauss) | `1.` | [log10(hPa)] vertical width; constant in the vertical if < 0 |
| `pcenter` | real, dimension(ngauss) | `-1.` | [hPa] pressure of the centre; surface heating if < 0 |
| `pmove` | real, dimension(ngauss) | `0.` | [hPa/day] vertical speed of the source |
| `is_periodic` | logical, dimension(ngauss) | `.false.` | reset the position periodically (with `tphase` and `tperiod`): periodic in longitude and pressure, back and forth in latitude |
| `twidth` | real, dimension(ngauss) | `-1.` | [days] temporal width; constant in time if < 0 |
| `tphase` | real, dimension(ngauss) | `0.` | [days] temporal phase |
| `tperiod` | real, dimension(ngauss) | `-1.` | temporal period: [fraction of a year] if < 0, [days] if > 0 |

<!-- mimadoc:end -->

## Tracers and input data

### `atmos_radon_nml`

Passive tracers of the tracer driver, used only if they are in the `field_table`: the radon tracer.

<!-- mimadoc:namelist atmos_radon_nml -->

Declared in [`atmos_radon_mod`](https://eddy-stanford.github.io/MiMA/api/atmos_radon_mod/) (`src/atmos_shared/tracer_driver/atmos_radon.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `ncopies_radon` | integer | `9` | number of copies of the radon tracer (`radon`, `radon_2`, ...) looked for in the field table; at most 9 |

<!-- mimadoc:end -->

### `atmos_convection_tracer_nml`

The convection tracer, used only if it is in the `field_table`.

<!-- mimadoc:namelist atmos_convection_tracer_nml -->

Declared in [`atmos_convection_tracer_mod`](https://eddy-stanford.github.io/MiMA/api/atmos_convection_tracer_mod/) (`src/atmos_shared/tracer_driver/atmos_convection_tracer.f90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `ncopies_cnvct_trcr` | integer | `9` | number of copies of the convection tracer looked for in the field table; at most 9 |

<!-- mimadoc:end -->

### `interpolator_nml`

MiMA's interpolator of climatology files, e.g. ozone and SST.

<!-- mimadoc:namelist interpolator_nml -->

Declared in [`mima_interpolator_mod`](https://eddy-stanford.github.io/MiMA/api/mima_interpolator_mod/) (`src/atmos_shared/interpolator/interpolator.F90`).

| Variable | Type | Default | Description |
|---|---|---|---|
| `read_all_on_init` | logical | `.false.` | read all time levels of a file at initialization |
| `verbose` | integer | `0` | amount of diagnostic printout |

<!-- mimadoc:end -->

## FMS namelists

These groups belong to the [FMS library](https://github.com/NOAA-GFDL/FMS) and keep FMS's defaults; see the FMS documentation for their other variables. The shipped `input.nml` sets:

 Group | Variable | Value in `input.nml` | Meaning
 :--- | :--- | :---: | :---
 `sat_vapor_pres_nml` | `do_simple` | `.true.` | **required**: the simple saturation vapour pressure table (Clausius-Clapeyron with constant latent heat, no ice) that MiMA uses. The FMS default (`.false.`, Goff-Gratch) is rejected at start-up.
 `topography_nml` | `topog_file`, `water_file` | `'INPUT/navy_topography.data.nc'`, `'INPUT/navy_pctwater.data.nc'` | high-resolution topography and water fraction, used by `topography_option = 'interpolated'` and `land_mask = 'navy'` of [`land_model_nml`](#land_model_nml)
 `fms_nml` | `domains_stack_size` | `600000` | stack size of the domain communication buffer; must be larger at higher resolution (the error message gives the size needed, e.g. about 1.2e6 at T255)
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
