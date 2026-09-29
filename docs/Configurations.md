[back to contents](README.md)

# Model configurations

This page describes some common ways of changing the model setup beyond the default [test case](GettingStarted.md#the-test-case). All options are set in `input.nml`; see [Parameter settings](Parameters.md) for the full list of namelist variables.

* [Radiation options](#radiation-options)
* [Specified initial conditions](#specified-initial-conditions)
* [Adding noise to the initial conditions](#adding-noise-to-the-initial-conditions)

## Radiation options

The radiation scheme is chosen with `radiation_scheme` in `radiation_nml`:

```fortran
&radiation_nml
    radiation_scheme = 'rrtm' /
```

* `'rrtm'` (default): RRTMG clear-sky radiation, configured with `rrtm_radiation_nml` and `astro_nml`.
* `'gray'`: the gray radiation scheme of Dargan Frierson ([Frierson, Held, Zurita-Gotor, JAS (2006)](https://doi.org/10.1175/JAS3753.1)), configured with `grey_radiation_nml`.
* `'none'`: no radiative heating and no radiative surface fluxes.

## Specified initial conditions

Without restart files, MiMA starts from an isothermal atmosphere at rest with a small vorticity perturbation. To start from your own initial state instead, set this flag in `spectral_dynamics_nml`:

```fortran
specify_initial_conditions = .true.
```

Then provide a netCDF file named `initial_conditions.nc` in the `INPUT/` directory where the model runs. It must contain zonal wind, meridional wind, temperature, specific humidity and surface pressure (`ucomp`, `vcomp`, `temp`, `sphum` and `ps`) at the model resolution. The file is only read on a cold start: if restart files are present in `INPUT/`, the model restarts from them instead. A zonally symmetric initial state stays zonally symmetric unless it is perturbed, e.g. with [noise](#adding-noise-to-the-initial-conditions).

## Adding noise to the initial conditions

MiMA can add random noise to the temperature field when the model state is loaded. This is useful for breaking symmetry in idealized experiments, or for generating ensembles of runs from the same starting point. It is controlled from `spectral_dynamics_nml`:

```fortran
add_noise                     = 0.1,
add_noise_seed                = 42,
noise_spectral_cutoff_minimum = 1,
noise_spectral_cutoff_maximum = 20,
```

Each spectral temperature coefficient between `noise_spectral_cutoff_minimum` and `noise_spectral_cutoff_maximum` (applied to both spectral indices) receives a perturbation drawn uniformly from [-`add_noise`, `add_noise`] K, on every model level. Noise is never added at the highest wavenumbers, since this leads to numerical instability. The default minimum of 1 leaves the zonal mean untouched.

* Noise is only added if `add_noise` is greater than zero (the default is `-1`, i.e. no noise).
* Set `add_noise_seed` to a non-negative integer to make the noise reproducible. Use different seeds for different ensemble members. With the default (`-1`) the compiler's default random seed is used, so runs may not be reproducible.
* The noise is added every time the model starts, whether from a cold start or from a restart file. If you run a long simulation in several segments, remove `add_noise` from `input.nml` after the first segment unless you want the model to be perturbed again at every restart.
