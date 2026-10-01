# Migrating from MiMA v1 to v2.0

This page is for users who have a working MiMA v1 setup (v1.2.x or older) and want to run it with v2.0. It lists everything that changed and what you have to do about it.

* [Overview](#overview)
* [Quick checklist](#quick-checklist)
* [Building](#building)
* [Namelist changes](#namelist-changes)
* [field_table changes](#field_table-changes)
* [diag_table changes](#diag_table-changes)
* [Output files](#output-files)
* [Restart files](#restart-files)
* [Removed components](#removed-components)
* [Changes to the answers](#changes-to-the-answers)
* [Other bug fixes](#other-bug-fixes)
* [New features](#new-features)

## Overview

MiMA v2.0 is a deliberate clean break from v1. It:

* **is trimmed to the idealized configurations MiMA is used for**: RRTM radiation with the simple mixed-layer surface (the default test case), gray radiation, and the Held-Suarez (1994) benchmark. Physics packages that no configuration used (the AM2 radiation, Donner, RAS, stratiform clouds, Mellor-Yamada and others; see [Removed components](#removed-components)) are gone, together with their namelist variables and diagnostics.
* **builds against the external [FMS](https://github.com/NOAA-GFDL/FMS) library**, release 2026.02 or later, instead of the copy of FMS from about 2005 that was bundled with v1.
* **changes the answers.** Several bugs that affected the shipped configurations were fixed, the Earth radius and the saturation vapour pressure table now follow FMS, and restarted runs now reproduce continuous runs. v2.0 results are not bit-for-bit identical to v1. The measured impact is listed in [Changes to the answers](#changes-to-the-answers).
* **makes the code defaults equal to the shipped `input/input.nml`**, so an `input.nml` that sets only the run length and the FMS settings (`&topography_nml`, `&fms_nml`, `&sat_vapor_pres_nml`) gives the standard MiMA setup. If your `input.nml` relied on the old defaults, you get different values now; see [Changed defaults](#changed-defaults).

The model itself (spectral dynamical core, RRTM, gray radiation, Betts-Miller convection, large-scale condensation, moist convective adjustment, the non-local boundary layer, the mixed-layer ocean with Q-fluxes, [`cg_drag`](https://eddy-stanford.github.io/MiMA/api/cg_drag_mod/), [`mg_drag`](https://eddy-stanford.github.io/MiMA/api/mg_drag_mod/), local heating, passive tracers) is unchanged apart from the bug fixes listed below. The [Fortran API reference](FortranAPI.md) describes every module of v2.0.

## Quick checklist

To move a v1 run directory to v2.0:

1. Build MiMA v2.0 with CMake 3.22 or later; FMS is found or downloaded automatically ([Building](#building)).
2. Convert `input.nml` ([Converting an old input.nml](#converting-an-old-inputnml)). The easiest way is often to start from the new `input/input.nml` and re-apply your own changes. At the very least:
    * replace `do_grey_radiation`/`do_rrtm_radiation` with `&radiation_nml radiation_scheme`;
    * rename `&grey_radiation_nml` to `&gray_radiation_nml`;
    * delete `&fms_io_nml` and add `&sat_vapor_pres_nml do_simple = .true. /`;
    * delete every variable that no longer exists (an unknown variable stops the model);
    * check the [changed defaults](#changed-defaults) for every variable you did not set.
3. Make sure `field_table` has a `sphum` tracer ([field_table changes](#field_table-changes)).
4. Update module and field names in `diag_table` (mainly the radiation diagnostics), and check it with `tools/mimadoc validate` ([diag_table changes](#diag_table-changes)).
5. Expect single netCDF output files, no need for `mppnccombine`, and a few format changes ([Output files](#output-files)).
6. v1 restart files can be read, but a new baseline is needed because the answers change ([Restart files](#restart-files)).

## Building

See [Getting started](GettingStarted.md#dependencies) for the full instructions.

| | v1 | v2.0 |
|---|---|---|
| FMS | bundled copy in `src/shared` (about 2005) | external [FMS](https://github.com/NOAA-GFDL/FMS) 2026.02 or later, with 8-byte reals (CMake target `FMS::fms_r8`) |
| How FMS is provided | built with MiMA | an installed FMS is used if CMake finds it (`CMAKE_PREFIX_PATH` or `FMS_ROOT`); otherwise FMS 2026.02 is downloaded (pinned by SHA256) and built with MiMA; offline: `-DFETCHCONTENT_SOURCE_DIR_FMS=/path/to/FMS-2026.02` |
| CMake | 3.16 or later | 3.22 or later |
| MPI | any MPI library | must provide the Fortran `mpi_f08` module (FMS 2026.02 uses it) |
| New CMake option | | `MIMA_OPENMP` (default `OFF` on macOS, `ON` elsewhere). MiMA has no OpenMP code; this only affects a downloaded FMS. |
| `mppnccombine` | needed to join per-processor output | not needed by default ([Output files](#output-files)) |

Notes:

* An installed FMS must use its default GFDL physical constants (`-DCONSTANTS=GFDL`; Spack `constants=GFDL`). MiMA checks this at start-up and stops if, for example, the GFS constants were chosen.
* MiMA also checks at start-up that FMS uses the simple saturation vapour pressure table, which needs `do_simple = .true.` in `&sat_vapor_pres_nml` (see below).
* FMS keeps its own compiler flags. Any `-ffp-contract` flag given for MiMA is passed on to a downloaded FMS, so that fused multiply-add contraction is the same in both.
* With the Intel compilers, a downloaded FMS is compiled with FMS's own Intel flags, so answers will not match builds that used v1's bundled FMS even apart from the changes below.
* Four MiMA modules were renamed so that they can be linked with FMS, which has modules of the same names: `diag_integral_mod`, `interpolator_mod` and `monin_obukhov_mod` are now [`mima_diag_integral_mod`](https://eddy-stanford.github.io/MiMA/api/mima_diag_integral_mod/), [`mima_interpolator_mod`](https://eddy-stanford.github.io/MiMA/api/mima_interpolator_mod/) and [`mima_monin_obukhov_mod`](https://eddy-stanford.github.io/MiMA/api/mima_monin_obukhov_mod/); [`fft`](https://eddy-stanford.github.io/MiMA/api/fft_mod/) and [`fft99`](https://eddy-stanford.github.io/MiMA/api/fft99_mod/) moved to `src/mima_shared`. This matters only if you have your own Fortran code that `use`s them.
* The radiation code moved to `src/atmos_param/radiation/` (`radiation.f90`, `gray_radiation.f90` and `rrtm/`). The gray module is now [`gray_radiation_mod`](https://eddy-stanford.github.io/MiMA/api/gray_radiation_mod/) (was `grey_radiation_mod`), and [`radiation_mod`](https://eddy-stanford.github.io/MiMA/api/radiation_mod/) chooses the scheme.

## Namelist changes

### How errors show up

v2.0 reads namelists with FMS's `check_nml_error`:

* A variable that no longer exists in a group that still exists is a **fatal error** ("Unknown namelist, or mistyped namelist variable in namelist ..."). You must delete it.
* A whole group that no longer exists (e.g. `&ocean_rough_nml`, `&fms_io_nml`) is **silently ignored**. Delete it anyway, so that your `input.nml` does not suggest settings that have no effect.
* A MiMA group that is missing takes the code defaults, which are now the values of the shipped `input/input.nml`. (The FMS groups keep FMS's defaults; e.g. `&topography_nml` must still point to the files in `INPUT/`.)

### Replaced settings

| v1 | v2.0 |
|---|---|
| `&physics_driver_nml do_rrtm_radiation = .true.` (the v1 default) | `&radiation_nml radiation_scheme = 'rrtm' /` (the default) |
| `&physics_driver_nml do_grey_radiation = .true., do_rrtm_radiation = .false.` | `&radiation_nml radiation_scheme = 'gray' /` |
| both `.false.` | `&radiation_nml radiation_scheme = 'none' /` (no radiative heating, zero radiative surface fluxes) |
| `&grey_radiation_nml` | `&gray_radiation_nml` (same variables) |
| `&fms_io_nml threading_write = 'single', fileset_write = 'single'` | removed. Restarts **and** diagnostics are single files by default; see `&spec_mpp_nml io_layout` |
| (nothing) | `&sat_vapor_pres_nml do_simple = .true. /` is **required**. Without it FMS uses its Goff-Gratch table and MiMA stops at start-up. |
| `&physics_driver_nml do_moist_processes = .false.` (moist physics moved to after the dynamics) | removed; moist physics always runs in the physics step. **Do not** translate it to `do_moist_physics = .false.`, which switches convection and condensation off. |
| `&vert_turb_driver_nml do_mellor_yamada = .true.` | removed; only the non-local K scheme (`do_diffusivity`, default now `.true.`) remains |
| `use_df_stuff = .true.` (in five groups) | removed; the `.true.` formulation is the only one |
| `&*_nml do_netcdf_restart` ([`physics_driver_nml`](Parameters.md#physics_driver_nml), [`atmos_model_nml`](Parameters.md#atmos_model_nml), [`mg_drag_nml`](Parameters.md#mg_drag_nml)) | removed; restarts are always netCDF |

### Removed variables in groups that still exist

Delete these from your `input.nml`. Variables marked with * were set in the v1 `input/input.nml`.

| Group | Removed variables | Why |
|---|---|---|
| [`physics_driver_nml`](Parameters.md#physics_driver_nml) | `do_grey_radiation`*, `do_rrtm_radiation`*, `do_radiation`, `do_moist_processes`, `do_netcdf_restart` | see [Replaced settings](#replaced-settings); `do_radiation` was the non-functional AM2 radiation |
| [`moist_processes_nml`](Parameters.md#moist_processes_nml) | `do_strat`*, `do_ras`*, `do_rh_clouds`*, `do_diag_clouds`*, `do_bmmass`*, `do_bmomp`*, `use_df_stuff`*, `do_donner_deep`, `do_cmt`, `do_dryadj`, `do_correct_q`, `qsrc` | schemes removed |
| [`moist_conv_nml`](Parameters.md#moist_conv_nml) | `beta`*, `use_df_stuff`* | detrainment into stratiform cloud removed |
| [`lscale_cond_nml`](Parameters.md#lscale_cond_nml) | `use_df_stuff`* | |
| [`diffusivity_nml`](Parameters.md#diffusivity_nml) | `do_entrain`*, `use_df_stuff`*, `entr_ratio`, `parcel_buoy`, `znom`, `free_atm_diff`, `free_atm_skyhi_diff`, `pbl_mcm`, `rich_crit_diff`, `mix_len`, `rich_prandtl`, `ampns`, `ampns_max` | free-atmosphere diffusion, PBL-top entrainment, parcel PBL depth and MCM option removed. The group is no longer in the shipped `input.nml`. |
| [`surface_flux_nml`](Parameters.md#surface_flux_nml) | `use_df_stuff`*, `raoult_sat_vap` | `raoult_sat_vap` had no effect (the surface never flags sea water) |
| [`vert_turb_driver_nml`](Parameters.md#vert_turb_driver_nml) | `do_mellor_yamada`*, `do_shallow_conv`*, `use_df_stuff`*, `do_edt`, `do_entrain`, `do_stable_bl` | only the non-local K scheme remains |
| [`vert_diff_driver_nml`](Parameters.md#vert_diff_driver_nml) | `do_mcm_no_neg_q`, `do_mcm_plev`, `do_mcm_vert_diff_tq` | MCM options removed |
| [`damping_driver_nml`](Parameters.md#damping_driver_nml) | `do_topo_drag` | `topo_drag` was a fatal stub |
| [`mg_drag_nml`](Parameters.md#mg_drag_nml) | `do_mcm_mg_drag`, `do_netcdf_restart` | |
| [`cg_drag_nml`](Parameters.md#cg_drag_nml) | `weighttop`*, `weightminus1`*, `weightminus2`*, `Bt_aug`, `Bt_eq_width`, `calculate_ked`, `num_diag_pts_ij`, `num_diag_pts_latlon`, `i_coords_gl`, `j_coords_gl`, `lat_coords_gl`, `lon_coords_gl` | read but never used |
| [`coupler_nml`](Parameters.md#coupler_nml) | `do_flux` | not used |
| [`gray_radiation_nml`](Parameters.md#gray_radiation_nml) | `wave_amp`, `wave_lon`, `wave_lat`, `wave_del_lon`, `wave_del_lat`, `wave_period`, `wave_env`, `wave_source` | the travelling-wave forcing they configured was disabled in the code |
| [`monin_obukhov_nml`](Parameters.md#monin_obukhov_nml) | `relax_time` | any value other than 0 was a fatal error |
| [`rrtm_radiation_nml`](Parameters.md#rrtm_radiation_nml) | `do_read_radiation`, `radiation_file`, `do_read_sw_flux`, `sw_flux_file`, `do_read_lw_flux`, `lw_flux_file`, `do_read_h2o`, `h2o_file`, `do_fixed_water`, `fixed_water`, `fixed_water_pres`, `fixed_water_lat`, `rad_missing_value` | file-driven radiation and water vapour removed (reading ozone, `do_read_ozone`, stays) |
| [`simple_surface_nml`](Parameters.md#simple_surface_nml) | `do_oflx`, `max_of`, `lonmax_of`, `latmax_of`, `latwidth_of`, `lonwidth_of`, `do_oflxmerid`, `maxofmerid`, `latmaxofmerid` | superseded by `&qflux_nml` |
| [`atmos_model_nml`](Parameters.md#atmos_model_nml) | `do_netcdf_restart` | |

### Removed option values

| Variable | Removed values | What happens now |
|---|---|---|
| `spectral_dynamics_nml vert_coord_option` | `'mcm'`, `'v197'` | fatal error; use `'even_sigma'`, `'uneven_sigma'`, `'hybrid'` or `'input'` |
| `spectral_dynamics_nml vert_difference_option` | `'mcm'` | fatal error; `'simmons_and_burridge'` is the only option |
| `simple_surface_nml roughness_choice` | `2` (`ocean_rough`) | fatal error; use 1, 3 or 4 (v1 left the roughness lengths unset for other values) |
| field_table `convection` method | `"ras"`, `"donner_and_ras"` | the tracer is not transported by convection (`"mca"`, `"mca_and_ras"`, `"donner_and_mca"`, `"all"` still mean MCA) |

### Removed namelist groups

These groups belong to removed code. They are ignored if present; delete them.

* Radiation (AM2): `radiation_driver_nml`, `sea_esf_rad_nml`, `astronomy_nml`, `aerosol_nml`, `aerosolrad_package_nml`, `radiative_gases_nml`, `ozone_nml`, `cloudrad_package_nml`, `cloud_spec_nml`, `isccp_clouds_nml`, `longwave_*_nml`, `shortwave_driver_nml`, `esfsw_*_nml`, `sealw99_nml`, `lhsw_driver_nml`, `gas_tf_nml`, `microphys_*_nml`, `*_clouds_w_nml`, and the other `sea_esf_rad` groups
* Clouds: `strat_cloud_nml`, `cloud_rad_nml`, `diag_cloud*`/`rh_based_clouds_nml`
* Convection: `donner_deep_nml`, `ras_nml`, `cu_mo_trans_nml`, `bm_massflux_nml`, `bm_omp_nml`, `dry_adj_nml`, `shallow_conv_nml`
* Turbulence: `my25_turb_nml`, `edt_nml`, `entrain_nml`, `stable_bl_turb_nml`
* Other: `ocean_rough_nml`, `atmosphere_nml` (`do_mcm_moist_processes`), `mcm_moist_processes_nml`, `mcm_mca_lsc_nml`, `flux_exchange_nml`, `fms_io_nml`, `mpp_io_nml`

The FMS groups that MiMA's input files use, `&fms_nml`, `&topography_nml` and `&gaussian_topog_nml`, still exist (now in the FMS library). `&diag_manager_nml` also exists, but modern FMS no longer has v1's `init_verbose` and `iospec`, and `&fms_nml` no longer has `iospec_ieee32`: setting them is a fatal error.

### New namelist variables

| Group | Variable | Default | Meaning |
|---|---|---|---|
| [`radiation_nml`](Parameters.md#radiation_nml) (new) | `radiation_scheme` | `'rrtm'` | `'rrtm'`, `'gray'` or `'none'` |
| [`physics_driver_nml`](Parameters.md#physics_driver_nml) | `do_held_suarez` | `.false.` | add the Held-Suarez (1994) forcing |
| | `do_boundary_layer` | `.true.` | boundary-layer turbulence, vertical diffusion and coupling to the surface fluxes. With `.false.` the surface state is left unchanged. |
| | `do_moist_physics` | `.true.` | convection and large-scale condensation |
| [`held_suarez_nml`](Parameters.md#held_suarez_nml) (new) | `t_zero`, `t_strat`, `delh`, `delv`, `p_ref`, `sigma_b`, `ka`, `ks`, `kf`, `do_rayleigh_friction`, `do_conserve_energy` | HS94 values | see [Configurations](Configurations.md#held-suarez-forcing) |
| [`spec_mpp_nml`](Parameters.md#spec_mpp_nml) (new) | `io_layout` | `1,1` | I/O layout for diagnostics and restarts; `1,1` writes one file each |
| `sat_vapor_pres_nml` (FMS) | `do_simple` | FMS default `.false.`; **must be set to `.true.`** | simple Clausius-Clapeyron table, as MiMA has always used |

See [Parameter settings](Parameters.md) for every variable and its default.

### Changed defaults

The code defaults were changed to the values in the shipped `input/input.nml` (BUG-08). An `input.nml` that sets a variable is not affected; an `input.nml` that relied on a v1 default now gets the v2.0 value. If your `input.nml` was derived from the shipped v1 file, almost all of these were already set in it and nothing changes.

| Group | Variable | v1 default | v2.0 default |
|---|---|---|---|
| [`coupler_nml`](Parameters.md#coupler_nml) | `dt_atmos` | 0 | 500 |
| | `days` | 0 | 360 |
| [`spectral_dynamics_nml`](Parameters.md#spectral_dynamics_nml) | `damping_order` | 2 | 4 |
| | `num_levels` | 18 | 40 |
| | `vert_coord_option` | `'even_sigma'` | `'uneven_sigma'` |
| | `scale_heights` | 4.0 | 7.9 |
| | `exponent` | 2.5 | 1.4 |
| | `ocean_topog_smoothing` | 0.93 | 0.995 |
| | `initial_sphum` | 0.0 | 2.e-6 |
| | `reference_sea_level_press` | 101325. | 1.e5 |
| | `water_correction_limit` | 0. | 200.e2 |
| [`rrtm_radiation_nml`](Parameters.md#rrtm_radiation_nml) | `do_read_ozone` | `.false.` | `.true.` |
| | `ozone_file` | `'ozone'` | `'ozone_1990'` |
| | `co2ppmv` | 300. | 390. |
| | `dt_rad` | 0 (every step) | 4500 |
| | `dt_rad_avg` | 86400 | 4500 |
| | `lonstep` | 1 | 4 |
| [`astro_nml`](Parameters.md#astro_nml) | `solr_cnst` | 1368.22 | 1370. |
| [`simple_surface_nml`](Parameters.md#simple_surface_nml) | `heat_capacity` | 4.e8 | 3.e8 |
| | `land_capacity` | -1 (= `heat_capacity`) | 1.e7 |
| | `trop_capacity` | -1 (= `heat_capacity`) | 1.e8 |
| | `trop_cap_limit` | 15. | 20. |
| | `const_albedo` | 0.30 | 0.23 |
| | `albedo_choice` | 1 | 7 |
| | `albedo_cntrSH`, `albedo_cntrNH` | 45., 65. | 64., 68. |
| | `albedo_wdth` | 10. | 5. |
| | `higher_albedo` | 0.38 | 0.80 |
| | `lat_glacier` | 45. | -70. |
| | `Tm` | 305. | 285. |
| | `roughness_choice` | 1 | 4 |
| | `mom_roughness_land`, `q_roughness_land` | 1., 1. | 5.e3, 1.e-12 |
| | `do_qflux`, `do_warmpool` | `.false.` | `.true.` |
| | `land_option` | `'none'` | `'interpolated'` |
| [`qflux_nml`](Parameters.md#qflux_nml) | `qflux_amp` | 30. | 26. |
| | `warmpool_amp`, `warmpool_width` | 5., 20. | 18., 35. |
| | `warmpool_k`, `warmpool_phase` | 1, 0. | 1.66666, 140. |
| | `warmpool_localization_choice` | 1 | 3 |
| | `gulf_phase`, `gulf_amp` | 140., 0. | 310., 70. |
| | `kuroshio_amp`, `trop_atlantic_amp`, `Hawaiiextra` | 0., 0., 0. | 40., 50., 30. |
| [`betts_miller_nml`](Parameters.md#betts_miller_nml) | `rhbm` | 0.8 | 0.7 |
| | `do_simp` | `.true.` | `.false.` |
| [`monin_obukhov_nml`](Parameters.md#monin_obukhov_nml) | `drag_min` | 1.e-5 | 4.e-5 |
| [`surface_flux_nml`](Parameters.md#surface_flux_nml) | `use_virtual_temp` | `.true.` | `.false.` |
| | `old_dtaudv` | `.false.` | `.true.` |
| [`vert_turb_driver_nml`](Parameters.md#vert_turb_driver_nml) | `do_diffusivity` | `.false.` | `.true.` |
| | `use_tau` | `.true.` | `.false.` |
| | `constant_gust` | 1.0 | 0. |
| [`vert_diff_driver_nml`](Parameters.md#vert_diff_driver_nml) | `do_conserve_energy` | `.false.` | `.true.` |
| | `use_virtual_temp_vert_diff` | `.true.` | `.false.` |
| [`damping_driver_nml`](Parameters.md#damping_driver_nml) | `do_cg_drag` | `.false.` | `.true.` |
| | `trayfric` | 0. | -0.5 |
| | `do_conserve_energy` | `.false.` | `.true.` |
| [`cg_drag_nml`](Parameters.md#cg_drag_nml) | `cg_drag_freq` | 0 | 21600 |
| | `damp_level_pressure` | 80. | 85. |
| | `Bt_0`, `Bt_eq` | 0.004, 0. | 0.0043, 0.0043 |
| | `Bt_nh`, `Bt_sh` | 0.001, -0.001 | 0., 0. |
| | `phi0n`, `phi0s` | 30., -30. | 15., -15. |
| | `dphin`, `dphis` | 5., -5. | 10., -10. |
| | `cw`, `cwtropics` | 40., 40. | 35., 35. |
| | `flag` | 1 | 0 |

Behaviour that depended on defaults of **removed** switches also changes:

* `vert_turb_driver_nml`: v1 ran Mellor-Yamada 2.5 unless `do_mellor_yamada = .false.` was set. v2.0 has no MY2.5; with the new default `do_diffusivity = .true.` it runs the non-local K scheme (as the v1 `input.nml` did).
* `diffusivity_nml do_entrain` defaulted to `.true.` in v1; there is no PBL-top entrainment in v2.0 (the v1 `input.nml` switched it off).
* `use_df_stuff` defaulted to `.false.` in `surface_flux_nml` and [`diffusivity_nml`](Parameters.md#diffusivity_nml); v2.0 always uses the `.true.` formulation (saturation specific humidity `d622*es/p`, latent heat of vaporization only, dry Richardson-number PBL depth), as every shipped configuration did.
* Radiation: v1 defaulted to RRTM (`do_rrtm_radiation = .true.`), so switching on only gray radiation was fatal; v2.0 still defaults to RRTM, selected with `radiation_scheme`.
* With `do_damping = .true.` (the default), [`cg_drag`](https://eddy-stanford.github.io/MiMA/api/cg_drag_mod/) is now on unless you set `do_cg_drag = .false.` in [`damping_driver_nml`](Parameters.md#damping_driver_nml).

### Converting an old input.nml

1. `&physics_driver_nml`: delete `do_grey_radiation`, `do_rrtm_radiation`, `do_radiation`, `do_moist_processes` and `do_netcdf_restart`. Add `&radiation_nml radiation_scheme = 'rrtm' /` (or `'gray'`, `'none'`).
2. Rename `&grey_radiation_nml` to `&gray_radiation_nml`.
3. Delete `&fms_io_nml`. Add

   ```fortran
   &sat_vapor_pres_nml
       do_simple = .true. /
   ```

4. Delete the variables in [Removed variables](#removed-variables-in-groups-that-still-exist) from the groups that still exist. In a namelist derived from the v1 `input/input.nml` these are: `do_bmmass`, `do_bmomp`, `do_strat`, `do_ras`, `do_diag_clouds`, `do_rh_clouds`, `use_df_stuff` ([`moist_processes_nml`](Parameters.md#moist_processes_nml)); `beta`, `use_df_stuff` ([`moist_conv_nml`](Parameters.md#moist_conv_nml)); `use_df_stuff` ([`lscale_cond_nml`](Parameters.md#lscale_cond_nml), [`surface_flux_nml`](Parameters.md#surface_flux_nml)); `do_entrain`, `use_df_stuff` ([`diffusivity_nml`](Parameters.md#diffusivity_nml)); `do_mellor_yamada`, `do_shallow_conv`, `use_df_stuff` ([`vert_turb_driver_nml`](Parameters.md#vert_turb_driver_nml)); `weighttop`, `weightminus1`, `weightminus2` ([`cg_drag_nml`](Parameters.md#cg_drag_nml)).
5. Delete the [groups of removed schemes](#removed-namelist-groups), e.g. `&ocean_rough_nml`.
6. Replace any [removed option values](#removed-option-values) (`vert_coord_option = 'mcm'`/`'v197'`, `roughness_choice = 2`, ...).
7. Go through the [changed defaults](#changed-defaults). For each variable your `input.nml` does not set, either accept the new value or set the v1 value explicitly.
8. Optionally set `&spec_mpp_nml io_layout` if you want per-processor files (e.g. for very large runs).
9. Do a short test run (e.g. `days = 1`) and check `logfile.000000.out`: every namelist is written there with the values actually used.

The shipped `input/input.nml` (RRTM), `input/examples/gray/input.nml` and `input/examples/held_suarez/input.nml` are complete v2.0 examples.

## field_table changes

* A humidity tracer (`sphum`, or `mix_rat`) is now **required**. v1 ran without one but then read memory out of bounds; v2.0 stops at initialization. Dry runs (e.g. Held-Suarez) keep a `sphum` tracer that stays zero and set `do_moist_physics = .false.`. The `dry_model` flag is removed.
* Tracer convection methods `"ras"` and `"donner_and_ras"` now mean no convective transport ([Removed option values](#removed-option-values)).

## diag_table changes

The `diag_table` format is unchanged. What changed are some module and field names, and the removed schemes' diagnostics. The default `input/diag_table` works unchanged.

A field that is requested but not registered is not an error: FMS prints a warning (`module/field_name (...) NOT registered`) and leaves the field out. So an old `diag_table` runs, but silently loses fields. Check it before you run:

```bash
python3 /path/to/MiMA/tools/mimadoc validate diag_table --nml input.nml
```

It reports fields that do not exist (with a hint when the field exists under another module, e.g. `rrtm_radiation/olr`), fields that are not available with your namelist settings, and format errors. [Diagnostics](Diagnostics.md) lists every field MiMA can write, with units, long names and the namelist settings it needs.

### Radiation diagnostics

Both schemes now register their diagnostics under module **`radiation`**, with one name per quantity. Change the module column (and the field name where it changed):

| v1 module / field | v2.0 `radiation` field | Note |
|---|---|---|
| `rrtm_radiation` / `tdt_rad`, `tdt_sw`, `tdt_lw` | `tdt_rad`, `tdt_sw`, `tdt_lw` | |
| `rrtm_radiation` / `olr` | `olr` | |
| `rrtm_radiation` / `flux_sw` | `swnet_sfc` | net SW at the surface (positive down) |
| `rrtm_radiation` / `flux_lw` | `lwdn_sfc` | downward LW at the surface |
| `rrtm_radiation` / `isr` | `swnet_toa` | it was net, not incoming, SW at TOA |
| `rrtm_radiation` / `rrtm_albedo` | `albedo_rad` | |
| `rrtm_radiation` / `coszen`, `ozone` | `coszen`, `ozone` | `ozone` now reports the constant `o3_val` when ozone is not read from a file (was uninitialized) |
| `rrtm_radiation` / `thalf` | `thalf` | now on the half-level axis `phalf`, including the top level (was on the full-level axis, BUG-09) |
| `grey_radiation` / `olr`, `tdt_rad`, `swdn_toa`, `lwdn_sfc`, `lwup_sfc`, `entrop_rad` | same names | |
| `grey_radiation` / `swdn_sfc` | `swnet_sfc` | it was net SW at the surface |
| `grey_radiation` / `flux_sw` | `swnet_half` | net SW on half levels (positive up) |
| `grey_radiation` / `flux_lw` | `lwnet_half` | net LW on half levels (positive up) |
| `grey_radiation` / `flux_rad` | `netrad_half` | net radiative flux on half levels (positive up) |
| `grey_radiation` / `tau` | `tau_lw` | LW optical depth on half levels |
| `grey_radiation` / `tau_rad` | `tau_sw` | SW optical depth on half levels |
| `grey_radiation` / `phalf_rad` | removed | use `pk + bk*ps` |

New shared fields: RRTM now also provides `swdn_toa` and `lwup_sfc`, and gray provides `swnet_toa` and `albedo_rad`, with the same names, units and long names. The v1 module `radiation` belonged to the AM2 radiation, which never ran in MiMA; none of its fields exist any more.

### Removed diagnostics

All diagnostics of removed schemes are gone:

| v1 module | Removed fields |
|---|---|
| `radiation` (AM2), `cloudrad`, `isccp` | all |
| `donner_deep`, `ras`, `cu_mo_trans`, `strat` | all |
| `edt`, `entrain`, `stable_bl_turb` | all |
| `vert_turb` | `tke`, `lscale`, `lscale_0`, `diff_t_stab`, `diff_m_stab`, `diff_t_entr`, `diff_m_entr`, `diff_sc` (the Mellor-Yamada, stable-BL and entrainment fields; the other `vert_turb` fields remain) |
| `moist` | the stratiform-cloud fields (`qldt_*`, `qidt_*`, `qadt_*`, `*_col` of liquid, ice and cloud fraction, `LWP`, `IWP`, `AWP`, `mc_full`), the Donner fields (`*_donner`, `<tracer>_donmca*`), `tdt_dadj` (dry adjustment), `massflux` (BM mass-flux), `capeflag` (never set) |
| `damping` | `udt_topo`, `vdt_topo` |
| `tracers` | `rnemiss` (never sent), `hook_no` |

### Diagnostics that now contain data

These were registered in v1 but never sent, so requesting them gave empty fields (and could corrupt the averaging bounds of the last record in the file):

* `simple_surface`: `tau_x`, `tau_y` (surface stress on the atmosphere, N/m2, positive eastward/northward), `t_ref`, `rh_ref` (at `z_ref_heat` = 2 m), `u_ref`, `v_ref` (at `z_ref_mom` = 10 m), and the Monin-Obukhov profile factors `del_h`, `del_m`, `del_q`. `rh_ref` can exceed 100% where evaporation is clipped (e.g. over cold land).
* `moist`: `<tracer>dt_conv` and `<tracer>dt_conv_col` for tracers transported by MCA.

### Units and long names

Only the `units` and `long_name` attributes changed (no field was renamed and no data changed). Scripts that match on the old strings need updating.

* Units are written one way throughout, parseable by UDUNITS: `K` (was `deg_k`, `deg_K`), `K/s` (was `deg_k/s`, `deg_K/s`), `Pa` (was `pascals`), `Pa/s` (was `Pa/sec`), `m/s` (was `m/sec`, `meters/second`), `m` (was `meters`), `1/s` (was `sec**-1`, ` /s`), `m/s2` (was `m/s**2`, `m/s^2`), `m2/s2` (was `(m/sec)**2`, `m**2/s**2`, `m^2/s^2`), `m2/s` (was `m^2/s`), `W/m2` and `W/m2/K` (was `w/m2`, `watts/m2`, `W/m^2`, `w/m2/K`), `J/kg` (was `J/Kg`), `N/m2` for the `mg_drag` stresses `taubx`, `tauby`, `taus` (was `kg/m/s2`), `1` for dimensionless fields (was `none`, `no units`, `?`), `kg/kg` for `ozone` (was `mmr`).
* Wrong units corrected (module `dynamics_every`): `tdt_damp` K/s (was m/s**2); `wsubt`, `wsubtv`, `kegen` K/s (were K); `kegenq` K/s kg/kg (was K); `kegenqtinv` kg/kg/s (was K); `qdt_watercor` kg/kg/s (was 1/s); `vq` m/s kg/kg (was m/s); `wq`, `wqp` Pa/s kg/kg (were Pa/s); `vqint` kg/m/s (was m^2/s); `<tracer>_hadv`, `<tracer>_vadv` `<tracer units> kg/m2/s`.
* Other wrong units: `moist/<tracer>_col` `<tracer units> kg/m2` and `moist/<tracer>dt_conv_col` `<tracer units> kg/m2/s` (both were the tracer's own units); `mca/<tracer>dt_MCA_col` `<tracer units> kg/m2/s`; `simple_surface/heat_capacity` J/m2/K (was none); tracer dry/wet deposition `<units> kg/m2/s` (was `kg/(m2 s)`).
* Corrected long names, mostly in `dynamics_every` (`uv`, `vq`, `vqint`, `vdse`, `vp`, `psi_dwc` and `psi_star`, which had the same name, `kegen*`, `<tracer>_hadv/_vadv`), and `radiation/ozone`, `moist/klzbs`, `moist/rhsurf`, `cg_drag/kedx_cgwd`, `kedy_cgwd`.

## Output files

| | v1 | v2.0 |
|---|---|---|
| Diagnostic files | one per MPI process (`atmos_daily.nc.0000`, ...), combined with `mppnccombine` | **one file** (`atmos_daily.nc`) by default. Set `&spec_mpp_nml io_layout` for per-I/O-domain files (each entry must divide the processor layout, which is `1, npes`); these still need `mppnccombine`. |
| Axis variables (`lon`, `lat`, `pfull`, `phalf`, `lonb`, `latb`, `nv`) | float | double |
| Time bounds | `time_bounds`, units `days` | `time_bnds`, units `days since 0001-01-01 00:00:00` (the same as `time`) |
| Fill value | `missing_value` only | `_FillValue` and `missing_value` on every field (`1.e+20` for most fields) |
| Axis attributes | `cartesian_axis = "X"` etc. | CF `axis = "X"` etc.; `standard_name = "none"` is no longer written |
| `time:calendar_type` | `THIRTY_DAY_MONTHS` | `360_DAY` (`time:calendar = "360_day"` in both) |
| Log files | `logfile.0000.out` | `logfile.000000.out`, plus `warnfile.000000.out` for warnings |
| Field log (with `&diag_manager_nml do_diag_field_log = .true.`) | `diag_field_log.out` | `diag_field_log.out.0` |

### Changes to diagnostic values (the model state is unchanged)

* **RRTM diagnostics** (`radiation/tdt_rad`, `olr`, ... with `radiation_scheme = 'rrtm'`) are sent at the end of the time step (`Time_next`), like all other physics, instead of at its start (BUG-69). They shift by one model step relative to v1; the first daily mean after a restart now matches a continuous run, and the last instantaneous record of a run is no longer empty.
* **Every-step dynamics diagnostics** (module `dynamics_every`) are stamped at the time they are valid, `Time + step*dt/num_steps`, instead of `Time + step*int(dt/2)` (BUG-20). With the default `num_steps = 1`, `t_every` and `ps_every` now equal `dynamics/temp` and `ps` exactly.
* Betts-Miller `invtaubmt`/`invtaubmq` were undefined in several branches, including the default one; they are now set.
* `mca/<tracer>dt_conv`, the `simple_surface` fields above and `radiation/ozone` contain data where they were empty or undefined.
* The global integrals ([`diag_integral`](https://eddy-stanford.github.io/MiMA/api/mima_diag_integral_mod/)) use 64-bit counters: with the default `output_interval = -1` (one print at the end of the run) the counter overflowed after about four model years at T42. A field that is never sent (e.g. `prec` in a dry run) is written as zero instead of stopping the model.

## Restart files

Restart files are still written to `RESTART/` and read from `INPUT/` with the same names, layout and variable names, so **v1 restart files can be read by v2.0**, apart from the fields that no longer exist. Because the answers change, a run continued from v1 restarts is not a continuation of the v1 run.

| File | Change |
|---|---|
| `cg_drag.res.nc` | **new**: `gwd_u`, `gwd_v` and the time to the next [`cg_drag`](https://eddy-stanford.github.io/MiMA/api/cg_drag_mod/) calculation (BUG-02). Without it (e.g. from v1 restarts) `cg_drag` cold-starts: no drag until `cg_drag_freq` has elapsed, as v1 did after every restart. |
| `rrtm_radiation.res.nc` | **new**: the time of the last radiation call, the stored heating rates and fluxes, and the precipitation-albedo accumulators when used (BUG-03). Without it RRTM recomputes radiation on the first step, as v1 did. |
| `physics_driver.res.nc` | now only `vers`, `diff_t`, `diff_m`. Dropped: `diff_cu_mo`, `pbltop`, `convect`, `doing_strat`, `doing_edt`, `doing_entrain`, `radturbten`, `lw_tendency`. |
| `spectral_physics.res.nc` | no longer written or read (it held state for the removed moist-after-dynamics path) |
| `atmos_coupled.res.nc` | the unused grid bounds `glon_bnd`, `glat_bnd` are no longer written (BUG-21). `dt` is still written but no longer read: `lprec` and `fprec` are rates and are not rescaled when `dt_atmos` changes (BUG-67). |
| `atmosphere.res.nc`, `spectral_dynamics.res.nc`, `simple_surface.res.nc`, `mg_drag.res.nc`, `coupler.res` | unchanged |

Other changes:

* **A run split into segments now ends bit-for-bit identical to a continuous run**, also when a segment ends between radiation calls. In v1, each segment started without gravity-wave drag for up to `cg_drag_freq` and recomputed radiation on its first step.
* Restart files are single files by default, like the output (v1 wrote single restart files only with `&fms_io_nml fileset_write = 'single'`, which the shipped `input.nml` set). Restart files split into per-processor pieces (`*.res.nc.0000`, ...) are read only with the `io_layout` that wrote them; otherwise the model stops and asks you to combine them with `mppnccombine`.
* The native (unformatted) restart format is no longer supported (`do_netcdf_restart` is removed). A native-format restart file in `INPUT/` without a netCDF one is a fatal error.
* A restart field whose size does not match the model grid is a fatal error.
* `cg_drag_freq` no longer has to divide a day.

## Removed components

| Component | Namelist switch in v1 | Why it was removed |
|---|---|---|
| AM2 radiation (`sea_esf_rad`, `radiation_driver`, `astronomy`, clouds, aerosols, ISCCP) | `physics_driver_nml do_radiation` | never functional in MiMA: the call to `radiation_driver` was commented out |
| MCM moist processes (`mcm_moist_processes`, `mcm_mca_lsc`) | `atmosphere_nml do_mcm_moist_processes` | could never run (its init called a fatal stub) |
| Diagnostic and RH cloud schemes (`diag_cloud`, `diag_cloud_rad`, `cloud_zonal`, `rh_clouds`) | `do_diag_clouds`, `do_rh_clouds` | null stubs |
| Donner deep convection | `do_donner_deep` | not used by any configuration |
| Relaxed Arakawa-Schubert and cumulus momentum transport | `do_ras`, `do_cmt` | not used by any configuration |
| Moist processes after the dynamics | `physics_driver_nml do_moist_processes = .false.` | not used; moist physics always runs in the physics step |
| Betts-Miller mass-flux and OMP variants | `do_bmmass`, `do_bmomp` | never initialized, so their namelists were ignored; standard Betts-Miller (`do_bm`) is unchanged |
| Dry convective adjustment | `do_dryadj` | not used by any configuration |
| Stratiform clouds (`strat_cloud`, `cloud_rad`, `cloud_generator`) | `do_strat` | not used by any configuration |
| Mellor-Yamada 2.5, EDT, entrainment, stable-BL and shallow-convection turbulence schemes | `do_mellor_yamada`, `do_edt`, `do_entrain`, `do_stable_bl`, `do_shallow_conv` | not used; the non-local K scheme (`do_diffusivity`) remains |
| `topo_drag` | `do_topo_drag` | null stub, fatal when called; `mg_drag` and `cg_drag` remain |
| `ocean_rough` roughness | `roughness_choice = 2` | not used by any configuration |
| The `use_df_stuff = .false.` moisture formulation | `use_df_stuff` | every configuration used `.true.` |
| Manabe Climate Model options | `vert_coord_option = 'mcm'`, `'v197'`, `vert_difference_option = 'mcm'`, `pbl_mcm`, `do_mcm_*` | not used by any configuration |
| Prescribed ocean heat fluxes in [`simple_surface`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) | `do_oflx`, `do_oflxmerid` | superseded by [`qflux_nml`](Parameters.md#qflux_nml) ([`qflux_mod`](https://eddy-stanford.github.io/MiMA/api/qflux_mod/)); `do_oflxmerid` read outside its table (BUG-46) |
| Ad-hoc humidity source | `do_correct_q`, `qsrc` | not used (BUG-53) |
| Free-atmosphere diffusion and PBL-top entrainment in [`diffusivity`](https://eddy-stanford.github.io/MiMA/api/diffusivity_mod/) | `free_atm_diff`, `do_entrain`, ... | not used by the shipped configurations |
| RRTM radiation and water vapour read from files | `do_read_radiation`, `do_read_sw_flux`, `do_read_lw_flux`, `do_read_h2o`, `do_fixed_water` | not used; `do_fixed_water` also used the wrong units (BUG-40) |
| Inert `cg_drag` parameters and column diagnostics | `weighttop`, ..., `num_diag_pts_*` | read but never used; the column-diagnostic code was never compiled |
| Native-format restarts | `do_netcdf_restart` | netCDF restarts only |
| Bundled FMS (`src/shared`) | | replaced by the external FMS library |
| Uncompiled and unused code (old `flux_exchange` coupler, solo driver, RRTMG McICA, `beta_dist`, test programs, CVS files) | | never built or never called |

## Changes to the answers

v2.0 does not reproduce v1. The fixes below change the results of the shipped configurations. Unless stated otherwise the impact was measured with 30-day runs of the default configuration from a common spun-up state (one year), each fix compared with the code just before it. For reference, a perturbation of 1e-6 in the CO2 concentration gives, after 30 days, a zonal-mean *u* rms difference of 0.25 m/s (max 1.4 m/s) and a zonal-mean *T* rms difference of 0.08 K (max 0.55 K): differences of that size are chaotic noise.

| Change | Affects | Zonal-mean *u*: rms (max) diff. [m/s] | Zonal-mean *T*: rms (max) diff. [K] | Notes |
|---|---|---|---|---|
| **BUG-01**: `cg_drag` is recomputed every `cg_drag_freq` (6 h by default); it was recomputed about every `cg_drag_freq/2` because the alarm was decremented by the leapfrog step 2Δt | `do_cg_drag = .true.` | 0.50 (2.9) | 0.15 (0.92) | after 1 day `udt_cgwd` differs by 0.36 m/s/day rms |
| **BUG-04**: the zonal surface stress is treated implicitly, like the meridional one (`dtaudu` was never returned, so it was explicit) | all | 0.47 (2.7) | 0.15 (0.81) | global precipitation +0.016 mm/day |
| **BUG-15**: the evaporation derivatives are zero where negative evaporation is clipped, so the implicit surface step cannot produce evaporation there | where the surface would take up water vapour | 0.36 (1.6) | 0.12 (1.2) | |
| **BUG-16**: Betts-Miller `do_shallower` conserves heat when shallow convection is confined to the lowest level | `do_bm` with `do_shallower = .true.` | 0 | 0 | not triggered in the 30-day test |
| **BUG-13**: momentum that `cg_drag` deposits above the model top is spread over the levels down to `damp_level_pressure` so that it is conserved (mass-weighted); v1 added the same acceleration at each level, about 1.8 times the escaping momentum | `do_cg_drag = .true.` | 2.7 (27) | 1.2 (13) | see below |
| **BUG-14**: all bottom-level fields passed to the surface fluxes (t, q, ps, pressure and height of the lowest level) are at the same time level | all | 0.25 (1.4) | 0.07 (0.66) | global precipitation -0.016 mm/day, sensible heat flux -0.11 W/m2 |
| All six fixes above together | | 2.8 (27) | 1.2 (14) | |
| **Earth radius** `RADIUS` = 6371 km, the FMS value (was 6376 km, although its comment said 6371 km) | all | 0.38 (2.4) | 0.09 (0.66) | |
| **Saturation vapour pressure table**: FMS's `do_simple` table, -173 to 350 °C at 0.1 K (was -200 to 250 °C at 0.2 K); same Clausius-Clapeyron formula. es changes by up to 1.5e-4 (relative), des/dT by up to 4.3e-3. | all moist runs | 0.25 (1.9) | 0.06 (0.70) | global precipitation -0.016 mm/day |
| Radius and SVP table together | | 0.44 (2.5) | 0.10 (0.67) | |

The radius and SVP changes, and BUG-01, 04, 14 and 15 individually, change the 30-day zonal means by one to two times the noise level: they are systematic, but small. **BUG-13 is not small**: it weakens the gravity-wave drag above `damp_level_pressure` (0.85 hPa) by about a third and changes the winds above about 1 hPa substantially. After 30 days:

| Pressure [hPa] | 0.18 | 0.56 | 0.94 | 1.6 | 3.3 | 11 | 34 | 92 |
|---|---|---|---|---|---|---|---|---|
| Zonal-mean *u* rms (max) diff. [m/s] | 11 (27) | 6.7 (17) | 4.9 (14) | 3.3 (9.1) | 1.7 (4.5) | 0.8 (2.0) | 0.4 (1.6) | 0.3 (1.3) |
| GWD rms, before → after the fix [m/s/day] | 11.2 → 7.3 | 11.5 → 7.5 | 1.6 → 1.5 | 0.8 → 0.7 | 0.6 → 0.6 | 0.2 → 0.2 | 0.1 → 0.1 | 0.1 → 0.1 |

If your work focuses on the upper stratosphere or mesosphere, the `cg_drag` tuning (`Bt_0`, `Bt_eq`, `cw`, ... in [`cg_drag_nml`](Parameters.md#cg_drag_nml)) may need revisiting.

Changes that affect only some runs:

* **Restarts** (BUG-02, BUG-03): a run in segments now equals a continuous run. In v1 each segment started without gravity-wave drag until the first `cg_drag` calculation and recomputed radiation at its first step; in the 5-day regression test the second segment changes by up to 2.6 m/s in *u* and 0.64 K in *T* after 5 days. Runs without restart files are not affected.
* **BUG-67**: after a restart with a different `dt_atmos`, the precipitation rates `lprec` and `fprec` read from `atmos_coupled.res.nc` are no longer rescaled by the ratio of the time steps.
* **Removed switches and changed defaults**: an `input.nml` that relied on v1 defaults now runs a different model; see [Changed defaults](#changed-defaults).
* **Saturation vapour pressure limits**: temperatures below about 100 K or above about 623 K in a saturation vapour pressure lookup are now fatal errors (v1: below 73 K or above 523 K).
* **Intel builds** with a downloaded FMS use FMS's own compiler flags for the FMS part.

Diagnostic output also changes where the model state does not: see [Changes to diagnostic values](#changes-to-diagnostic-values-the-model-state-is-unchanged).

## Other bug fixes

These do not change the answers of the shipped configurations.

* Betts-Miller: guarded two out-of-bounds accesses (a parcel buoyant up to the model top, BUG-05; the LCL table read past its end, BUG-10); also in the CAPE/CIN diagnostics. Removed a stray `hi lcl` print (BUG-27).
* [`cg_drag`](https://eddy-stanford.github.io/MiMA/api/cg_drag_mod/): latitudes indexed with global offsets in a local array (BUG-32); undefined damping level if no level lies above `damp_level_pressure` (BUG-33).
* [`spectral_init_cond`](https://eddy-stanford.github.io/MiMA/api/spectral_init_cond_mod/): an input topography file is accepted only if both dimensions match the grid (BUG-48).
* The interpolator ([`mima_interpolator_mod`](https://eddy-stanford.github.io/MiMA/api/mima_interpolator_mod/)) accepts the calendar name `'360'` as well as `'360_day'` (BUG-52). It no longer writes out of bounds for files whose record count is not a multiple of 12 (BUG-62) or reads past its input in linear-in-pressure interpolation (BUG-63); seasonal files with 1-4 records now work (BUG-64); files without a time dimension, or with a fixed-size `time` dimension as xarray writes, give one field valid at all times (BUG-65). Monthly climatologies such as the ozone file give identical results.
* [`simple_surface`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) `do_sc_sst = .true.` (SSTs read from a file) no longer crashes (BUG-66).
* `specify_initial_conditions`: `initial_conditions.nc` is read with the netCDF-Fortran 90 interface ([`spectral_initialize_fields`](https://eddy-stanford.github.io/MiMA/api/spectral_initialize_fields_mod/#spectral_initialize_fields)). Each variable must have the model grid's shape (a larger file used to be read as a corner sub-block without warning), float variables are converted correctly, a trailing time dimension of length 1 is accepted, and errors name the file and variable.
* RRTM: [`rrtm_radiation_end`](https://eddy-stanford.github.io/MiMA/api/rrtm_radiation/#rrtm_radiation_end) deallocates its arrays (BUG-30); the allocation guard for file-driven radiation (BUG-38) is gone with that option.
* NULL pointer arguments of the removed AM2 surface fields are gone (BUG-18); `radturbten` and `lw_tendency`, previously unset, are gone.
* A missing humidity tracer is fatal instead of an out-of-bounds read (BUG-06).
* [`local_heating`](https://eddy-stanford.github.io/MiMA/api/local_heating_mod/) no longer depends on the RRTM code (BUG-23) or on the C preprocessor for its namelist (BUG-22).
* [`vert_diff`](https://eddy-stanford.github.io/MiMA/api/vert_diff_mod/) writes its log line only on the root PE (BUG-25); `time_stamp.out` and `coupler.res` are written by the root PE only.

## New features

* **Radiation scheme selection** with `&radiation_nml radiation_scheme = 'rrtm' | 'gray' | 'none'`, one module (`radiation`) and one set of names for the radiation diagnostics.
* **Held-Suarez (1994) forcing** (`do_held_suarez`, `&held_suarez_nml`), with new switches `do_boundary_layer` and `do_moist_physics`. It works at any resolution and with any vertical levels, and can be combined with the Rayleigh sponge, `cg_drag`, topography, and MiMA's moist physics and surface (moist variants). Example: `input/examples/held_suarez/` (dry HS94, T42 L20). See [Held-Suarez forcing](Configurations.md#held-suarez-forcing).
* **Gray radiation example** `input/examples/gray/`: the default test case with gray radiation, and a `diag_table` with the gray fluxes.
* **Diagnostics reference** [Diagnostics](Diagnostics.md), generated from the source by `tools/mimadoc`, which can also **validate a diag_table** against the source and an `input.nml`.
* **Surface stress and reference-height diagnostics** in `simple_surface`: `tau_x`, `tau_y`, `t_ref`, `rh_ref`, `u_ref`, `v_ref`, `del_h`, `del_m`, `del_q`.
* **More shared radiation fields**: RRTM `swdn_toa`, `lwup_sfc`; gray `swnet_toa`, `albedo_rad`.
* **Single-file output and restarts** by default, with `&spec_mpp_nml io_layout` for split files.
* **Exact restarts**: `cg_drag` and RRTM state are saved, so a run in segments equals a continuous run.
* **Consistent defaults**: the code defaults are the standard test case, so an `input.nml` needs only the settings that differ from it (plus the FMS settings `&topography_nml`, `&fms_nml` and `&sat_vapor_pres_nml`).
* **Start-up checks** that FMS uses the GFDL physical constants and the simple saturation vapour pressure table.
* **Modern FMS** (2026.02) and build: CMake 3.22, `find_package(FMS)` with a pinned download as fallback, optional OpenMP.
* **CF/UDUNITS-style units** for all diagnostics.
