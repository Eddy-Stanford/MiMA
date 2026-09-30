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
* `'gray'`: the gray radiation scheme of Dargan Frierson ([Frierson, Held, Zurita-Gotor, JAS (2006)](https://doi.org/10.1175/JAS3753.1)), configured with `gray_radiation_nml`.
* `'none'`: no radiative heating and no radiative surface fluxes.

The two schemes provide different diagnostics (the shared ones, such as `olr` and `tdt_rad`, have the same names); see the `radiation` module in [Diagnostics](Diagnostics.md#module-radiation). A complete gray-radiation setup is in `input/examples/gray/`: its `input.nml` is the default test case with `radiation_scheme = 'gray'` and `gray_radiation_nml` in place of `rrtm_radiation_nml` and `astro_nml`, and its `diag_table` writes daily and 30-day means including the gray radiative fluxes. It uses the same `INPUT/` files as the test case.

## Held-Suarez forcing

MiMA can run the [Held and Suarez (1994)](https://doi.org/10.1175/1520-0477(1994)075<1825:APFTIO>2.0.CO;2) idealized forcing: Newtonian relaxation of temperature towards a zonally symmetric equilibrium profile, and Rayleigh friction of the winds in the boundary layer. It is switched on in `physics_driver_nml`, and its parameters are set in `held_suarez_nml` (the defaults are the HS94 values):

 Variable | Default | Meaning
 :--- | :---: | :---
 `t_zero`, `t_strat` | 315, 200 K | surface equilibrium temperature at the equator, minimum equilibrium temperature
 `delh`, `delv` | 60, 10 K | equator-to-pole temperature difference, vertical potential temperature difference
 `sigma_b` | 0.7 | top of the frictional boundary layer (sigma)
 `ka`, `ks`, `kf` | 40, 4, 1 days | free-atmosphere and surface relaxation times, boundary-layer friction time
 `do_rayleigh_friction` | `.true.` | apply the boundary-layer friction
 `do_conserve_energy` | `.false.` | heat the air by the frictional dissipation

The forcing depends only on latitude and sigma, so it works at any horizontal resolution and with any vertical levels.

For the standard **dry** benchmark, switch off radiation, moist physics and the boundary layer, which also leaves the surface state unchanged. A complete example (T42, 20 evenly spaced sigma levels, flat topography) is in `input/examples/held_suarez/`; it needs no input data files:

```fortran
&radiation_nml
    radiation_scheme = 'none' /

&physics_driver_nml
    do_held_suarez    = .true.,
    do_boundary_layer = .false.,
    do_moist_physics  = .false.,
    do_damping        = .false. /
```

The model always carries a humidity tracer (`sphum` in the field table); in the dry setup it stays zero. Also set `use_virtual_temperature = .false.` and `do_water_correction = .false.` in `spectral_dynamics_nml`.

The HS forcing can be combined with other parts of the model:
* `do_damping = .true.` with `damping_driver_nml` enables the Rayleigh sponge (`do_rayleigh`) and/or the convective gravity-wave drag (`do_cg_drag`). Note that `do_cg_drag` defaults to `.true.`, so set `do_cg_drag = .false.` if you want only the sponge.
* Non-flat topography through `topography_option` in `spectral_dynamics_nml`.
* **Moist variants:** with `do_moist_physics = .true.` and `do_boundary_layer = .true.` (and `do_rayleigh_friction = .false.`), the HS temperature relaxation replaces radiation while MiMA's moist physics, boundary layer and surface fluxes stay active. This is similar in spirit to the moist Held-Suarez test of [Thatcher and Jablonowski (2016)](https://doi.org/10.5194/gmd-9-1263-2016), but uses MiMA's own boundary-layer and surface schemes.
  With `radiation_scheme = 'none'` the surface receives no radiation, so hold the SST fixed with `surface_choice = 2` in `simple_surface_nml` (its initial profile is set by `Tm` and `deltaT`); a slab ocean would otherwise cool without limit. The HS equilibrium temperature near the equatorial surface (315 K) is warmer than typical SSTs, so the lowest layers are stably stratified over the ocean and the hydrological cycle is weak: with the default SST profile (about 298 K at the equator) precipitation takes about three weeks to start and settles near 0.5 mm/day in the global mean.

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
