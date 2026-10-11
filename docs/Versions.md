# Version history

This file is to be updated whenever major changes are implemented in MiMA.
This is an important part of the documentation, as it might be quite important to know which version
includes which physical feature, bugfix, etc.

Any contributor to the code should change this file when merging into the main branch.
The newest changes should be on the top of this list.

## v2.0

A clean break from v1: results are not bit-for-bit identical to v1, and `input.nml`, `diag_table` and output-processing scripts need updating. See the [migration guide](Migration_v2.md).

* Unreleased:
    * The model is trimmed to its idealized configurations: RRTM radiation with the mixed-layer surface (the default), gray radiation, and Held-Suarez. Unused physics (the AM2 radiation, Donner, RAS, stratiform clouds, Mellor-Yamada, EDT, dry adjustment, `topo_drag` and others) and their namelist variables and diagnostics are removed.
    * Built against the external [FMS](https://github.com/NOAA-GFDL/FMS) library, release 2026.01.01 or later, instead of the bundled copy; CMake finds an installed FMS or downloads it. CMake 3.22 or later is required. `&sat_vapor_pres_nml do_simple = .true.` is required and `&fms_io_nml` is gone.
    * The radiation scheme is chosen with [`radiation_nml`](Parameters.md#radiation_nml) `radiation_scheme` (`'rrtm'`, `'gray'` or `'none'`), replacing `do_rrtm_radiation`/`do_grey_radiation`; the radiation diagnostics of both schemes are under module `radiation` with unified names.
    * The code defaults equal the shipped `input/input.nml`.
    * The surface is split into components in a new `src/surface/` directory, following the names and the call sequence of the GFDL AM4 coupler but running in the same process on the atmospheric grid: [`flux_exchange_mod`](https://eddy-stanford.github.io/MiMA/api/flux_exchange_mod/), [`ocean_model_mod`](https://eddy-stanford.github.io/MiMA/api/ocean_model_mod/), [`land_model_mod`](https://eddy-stanford.github.io/MiMA/api/land_model_mod/), [`surface_albedo_mod`](https://eddy-stanford.github.io/MiMA/api/surface_albedo_mod/) and [`ocean_qflux_mod`](https://eddy-stanford.github.io/MiMA/api/ocean_qflux_mod/). `&simple_surface_nml` is replaced by `&ocean_model_nml`, `&land_model_nml`, `&surface_albedo_nml` and `&flux_exchange_nml` (the model stops if it is still present; [mapping](Migration_v2.md#old-to-new-namelist-mapping)), the diagnostics are in the modules `flux`, `ocean` and `land`, and the surface restart files are `ocean_model.res.nc` and `land_model.res.nc`.
    * New surface options: prescribed SST over the ocean with interactive land, fixed ocean and land, Q-flux read from a file with the new diagnostic `ocean/qflux_implied`, the land-sea mask from a file, boxes, height or the Navy data, and an optional bucket land with snow (verified for conservation, not climate-validated). See [Surface options](Configurations.md#surface-options).
    * Precipitation can now fall as snow: `&moist_processes_nml precip_phase` (`'temperature'` by default, `'liquid'` gives v1) partitions it by the temperature of the lowest 150 hPa, with the latent heat of sublimation released for the snow (BUG-121; known limits: saturation over liquid, Betts-Miller parcel with `hlv`).
    * Answer-changing surface fixes: the Q-flux is applied to the ocean only and its ocean mean is removed (`qflux_ocean_only`, `qflux_remove_mean`; both `.false.` restores v1), precipitation falls as snow in cold columns (`precip_phase = 'liquid'` restores v1), `albedo_choice = 5` uses the signed latitude (BUG-95), `zsurf_cap_limit` is in metres (BUG-96), the land roughness applies to every land-mask source, the derivatives of the radiative and momentum fluxes respect `no_surface_radiative_flux` and `no_surface_momentum_flux` (BUG-98), `surface_flux` receives the ocean mask as `seawater` (BUG-97, only with `ncar_ocean_flux`), and the local surface heating is coupled consistently with the atmosphere (BUG-120, only with `pcenter < 0`). See [Changes to the answers](Migration_v2.md#changes-to-the-answers).
    * Held-Suarez (1994) forcing (`do_held_suarez`, [`held_suarez_nml`](Parameters.md#held_suarez_nml), [`held_suarez_mod`](https://eddy-stanford.github.io/MiMA/api/held_suarez_mod/)) with the new switches `do_boundary_layer` and `do_moist_physics`; example configurations `input/examples/held_suarez` and `input/examples/gray`.
    * Answer-changing fixes: [`cg_drag`](https://eddy-stanford.github.io/MiMA/api/cg_drag_mod/) recomputed at the intended interval and conserving the momentum deposited above the model top (this changes the winds above about 1 hPa), implicit zonal surface stress, consistent bottom-level time levels, evaporation derivatives, Betts-Miller shallow-convection energy, Earth radius 6371 km and the FMS saturation vapour pressure table.
    * `cg_drag` and RRTM state are saved in new restart files, so runs in segments reproduce continuous runs. v1 restart files can still be read, except for the surface state (`simple_surface.res.nc`, which has to be converted to the new surface restart files, see [Restart files](Migration_v2.md#restart-files)).
    * Output and restarts are single files by default (`spec_mpp_nml io_layout`), so `mppnccombine` is no longer needed. Axes are written in double precision, time bounds are named `time_bnds`, fields have a `_FillValue`, and units follow UDUNITS.
    * New [diagnostics reference](Diagnostics.md), generated from the source by `tools/mimadoc`, which can also check a `diag_table`. The surface stress and reference-height diagnostics (now in module `flux`) are now written.
    * New documentation site built with MkDocs, with a [Fortran API reference](FortranAPI.md) generated from doc comments added throughout the source. The [namelist reference](Parameters.md) is now generated from the source too.

## v1.2.X

Maintained at [Eddy-Stanford/MiMA](https://github.com/Eddy-Stanford/MiMA).

* Unreleased: Documentation reorganised and updated for the CMake build. Outdated Apptainer definition `mima.def` removed. Python is no longer a build dependency. The `postprocessing/` directory (`mppnccombine`, `plevel_interpolation`, `output_to_input.py`) and the `BUILD_COMBINE` CMake option were removed. Use `mppnccombine` and `plevel.sh` from [FRE-NCtools](https://github.com/NOAA-GFDL/FRE-NCtools) instead.
* v1.2.5 (March 2026): Improved initial-condition noise. `add_noise_seed`, `noise_spectral_cutoff_minimum` and `noise_spectral_cutoff_maximum` added to [`spectral_dynamics_nml`](Parameters.md#spectral_dynamics_nml). See [Adding noise to the initial conditions](Configurations.md#adding-noise-to-the-initial-conditions).
* v1.2.3 – v1.2.4 (February 2026): `add_noise` option in `spectral_dynamics_nml` adds random thermal noise to the temperature field at start-up. Development container fixed.
* v1.2.2 (January 2026): `diag_integral` flushes its output buffer and gives more helpful error messages. `CITATION.cff` added.
* v1.2.1 (October 2025): The old `mkmf` build system (`bin/`) was removed, and CMake is now the only build system. New CMake options: `BUILD_COMBINE` builds `mppnccombine`, and `INSTALL_EXEC` creates a ready-to-run `exec/` directory.
* v1.2 (January 2025): CMake build supports GNU and Intel (including `ifx`) compilers on both x86-64 and arm64. Apptainer container definition (`mima.def`) added. Output uses the CF-compliant calendar name `360_day`.

## Unreleased upstream changes (2020–2023)

These changes were made after v1.1 and are included from v1.2 onwards (the `legacy` tag marks the state of the code just before the v1.2 changes).

* CMake build system added alongside `mkmf`.
* Stationary-wave configuration of [Garfinkel et al. (2020)](https://doi.org/10.1175/JCLI-D-19-0181.1), previously referred to as v2.0: zonally asymmetric Q fluxes ([`qflux_nml`](Parameters.md#qflux_nml)), desert albedo (`albedo_choice = 7`), land surface changes in `simple_surface`, and fixes to the gravity-wave drag (`cg_drag`). Navy topography and land-sea mask files added to `input/INPUT/`, and the sample `input.nml` updated to this configuration.

## v1.X

* v1.1 (February 2020): Restart from arbitrary initial conditions (`specify_initial_conditions`), switches to turn off surface fluxes, local (and moving) Gaussian heating ([`local_heating_nml`](Parameters.md#local_heating_nml)), option to use the Navy high-resolution land-sea mask, and `postprocessing/output_to_input.py`.
* v1.0.1: Patch for v1.0. Addresses incoming solar radiation issues when diurnal cycle averaging is used.
* v1.0: initial published version, as described in [Jucker and Gerber, J Clim (2017)](https://doi.org/10.1175/JCLI-D-17-0127.1).
