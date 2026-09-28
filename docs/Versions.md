[back to contents](README.md)

# Version history

This file is to be updated whenever major changes are implemented in MiMA.
This is an important part of the documentation, as it might be quite important to know which version
includes which physical feature, bugfix, etc.

Any contributor to the code should change this file when merging into the master branch.
The newest changes should be on the top of this list.

## v1.2.X

Maintained at [Eddy-Stanford/MiMA](https://github.com/Eddy-Stanford/MiMA).

* Unreleased: Documentation reorganised and updated for the CMake build. Outdated Apptainer definition `mima.def` removed. Python is no longer a build dependency. The `postprocessing/` directory (`mppnccombine`, `plevel_interpolation`, `output_to_input.py`) and the `BUILD_COMBINE` CMake option were removed. Use `mppnccombine` and `plevel.sh` from [FRE-NCtools](https://github.com/NOAA-GFDL/FRE-NCtools) instead.
* v1.2.5 (March 2026): Improved initial-condition noise. `add_noise_seed`, `noise_spectral_cutoff_minimum` and `noise_spectral_cutoff_maximum` added to `spectral_dynamics_nml`. See [Adding noise to the initial conditions](Configurations.md#adding-noise-to-the-initial-conditions).
* v1.2.3 – v1.2.4 (February 2026): `add_noise` option in `spectral_dynamics_nml` adds random thermal noise to the temperature field at start-up. Development container fixed.
* v1.2.2 (January 2026): `diag_integral` flushes its output buffer and gives more helpful error messages. `CITATION.cff` added.
* v1.2.1 (October 2025): The old `mkmf` build system (`bin/`) was removed, and CMake is now the only build system. New CMake options: `BUILD_COMBINE` builds `mppnccombine`, and `INSTALL_EXEC` creates a ready-to-run `exec/` directory.
* v1.2 (January 2025): CMake build supports GNU and Intel (including `ifx`) compilers on both x86-64 and arm64. Apptainer container definition (`mima.def`) added. Output uses the CF-compliant calendar name `360_day`.

## Unreleased upstream changes (2020–2023)

These changes were made after v1.1 and are included from v1.2 onwards (the `legacy` tag marks the state of the code just before the v1.2 changes).

* CMake build system added alongside `mkmf`.
* Stationary-wave configuration of [Garfinkel et al. (2020)](https://doi.org/10.1175/JCLI-D-19-0181.1), previously referred to as v2.0: zonally asymmetric Q fluxes (`qflux_nml`), desert albedo (`albedo_choice = 7`), land surface changes in `simple_surface`, and fixes to the gravity-wave drag (`cg_drag`). Navy topography and land-sea mask files added to `input/INPUT/`, and the sample `input.nml` updated to this configuration.

## v1.X

* v1.1 (February 2020): Restart from arbitrary initial conditions (`specify_initial_conditions`, as used in [Yamada and Pauluis (2017)](https://doi.org/10.1175/JAS-D-16-0329.1)), switches to turn off surface fluxes for life cycle experiments, local (and moving) Gaussian heating (`local_heating_nml`), option to use the Navy high-resolution land-sea mask, and `postprocessing/output_to_input.py`.
* v1.0.1: Patch for v1.0. Addresses incoming solar radiation issues when diurnal cycle averaging is used.
* v1.0: initial published version, as described in [Jucker and Gerber, J Clim (2017)](https://doi.org/10.1175/JCLI-D-17-0127.1).
