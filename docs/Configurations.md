[back to contents](README.md)

# Model configurations

This page describes some common ways of changing the model setup beyond the default [test case](GettingStarted.md#the-test-case). All options are set in `input.nml`; see [Parameter settings](Parameters.md) for the full list of namelist variables.

* [Radiation options](#radiation-options)
* [Life cycle calculations](#life-cycle-calculations)
* [Adding noise to the initial conditions](#adding-noise-to-the-initial-conditions)

## Radiation options

By default, MiMA uses the RRTM radiation code. This is set by `do_rrtm_radiation = .true.` (default). There are, however, two more options for radiation, described below.

MiMA includes the gray radiation scheme developed by Dargan Frierson ([Frierson, Held, Zurita-Gotor, JAS (2006)](https://doi.org/10.1175/JAS3753.1)). To switch between the radiation schemes, the flags `do_grey_radiation` and `do_rrtm_radiation` in the namelist `physics_driver_nml` can be set accordingly (only one of them should be `.true.` of course).

Theoretically, there is also the possibility of running the full AM2 radiation scheme, with the flag `do_radiation` in `physics_driver_nml`. However, this option will need a lot of input files for tracer concentration, which are not part of the MiMA repository. This option, although all the relevant files are present and being compiled, has never been tested, and should only be used with great caution.

## Life cycle calculations

MiMA can be used to run life cycle experiments, as explored in Yamada and Pauluis (2017). To specify the initial conditions, activate this flag in `spectral_dynamics_nml`:

```fortran
specify_initial_conditions = .true.
```

Then a netCDF file containing the initial conditions for zonal wind, meridional wind, temperature, specific humidity, and surface pressure (`ucomp`, `vcomp`, `temp`, `sphum`, and `ps`, respectively) must be provided. It should be at the resolution of the model. It should be named `initial_conditions.nc` and placed in the `INPUT/` directory where the model is executed. Note that if you do not include a slight zonal perturbation, the model will maintain a zonally symmetric state, stuck to the unstable fixed point. There are different strategies for exciting zonal asymmetries. You can add random noise (see [below](#adding-noise-to-the-initial-conditions)), or focus in on a particular wavenumber, as detailed below.

A traditional life cycle is run with no forcing. To shut off all diabatic processes, you must make these adjustments to the namelist. To turn off radiation and damping (except for hyperdiffusion), add these options to `physics_driver_nml`:

```fortran
do_grey_radiation = .false.,
do_rrtm_radiation = .false.,
do_damping = .false.
```

Then, to turn off any diabatic forcings at the lower boundary, in `surface_flux_nml` add these options:

```fortran
no_surface_momentum_flux  = .true.,
no_surface_moisture_flux  = .true.,
no_surface_heat_flux      = .true.,
no_surface_radiative_flux = .true.
```

The spectral dynamical core allows one to focus in on a particular wavenumber, as was done by Yamada and Pauluis (2017). For example, to run a T170 resolution model, but enforce 6-fold symmetry (i.e., only capture instabilities at wave 6 and harmonics), use these options in `spectral_dynamics_nml`:

* `lon_max                 = 128,`     [an ideal number for the Fourier transforms, close to 512/6]
* `lat_max                 = 256,`     [this grid corresponds to T170 resolution]
* `num_fourier             = 29,`      [this is approximately 170/6]
* `num_spherical           = 171,`     [this is always the T-resolution + 1]
* `fourier_inc             = 6,`       [this allows zonal waves 0, 6, 12, ...]

This trick allows you to run a higher resolution integration about 6 times faster.

Lastly, note that hyperdiffusion is still required for stability. For the Yamada and Pauluis life cycles, these options were selected in `spectral_dynamics_nml`:

```fortran
damping_option          = 'resolution_dependent',
damping_order           = 3,
damping_coeff           = 6.94444444e-5,
damping_order_vor       = 3,
damping_order_div       = 3,
damping_coeff_vor       = 6.94444444e-5,
damping_coeff_div       = 6.94444444e-5,
```

Reference: [Yamada, R., and O. Pauluis, 2017: Wave-mean-flow interactions in moist baroclinic lifecycles. J. Atmos. Sci., 74, 2143-2162, doi:10.1175/JAS-D-16-0329.1](https://doi.org/10.1175/JAS-D-16-0329.1).

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
