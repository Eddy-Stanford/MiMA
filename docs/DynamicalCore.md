# The dynamical core

MiMA's dynamical core is the GFDL spectral dynamical core: it solves the hydrostatic primitive equations on the sphere in vorticity-divergence form, with spherical harmonics in the horizontal and the Simmons and Burridge (1981) finite differences in the vertical. The numerics are described in the GFDL note [The Spectral Dynamical Core](references/gfdl_spectral_core.pdf); this page summarizes them in MiMA's notation and lists [what MiMA does differently](#what-mima-changes-relative-to-the-gfdl-note). The options are set in [`spectral_dynamics_nml`](Parameters.md#spectral_dynamics_nml), and the code is in `src/atmos_spectral/` (`atmosphere/`, `dynamics/`, `diagnostics/`, `init/` and `transforms/`), with [`spectral_dynamics_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_dynamics_mod/) as the main module.

* [Equations](#equations)
* [Spectral layout](#spectral-layout)
* [Transforms and decomposition](#transforms-and-decomposition)
* [Time step](#time-step)
* [Damping and sponge](#damping-and-sponge)
* [Mass, energy and water corrections](#mass-energy-and-water-corrections)
* [Initial state](#initial-state)
* [What MiMA changes relative to the GFDL note](#what-mima-changes-relative-to-the-gfdl-note)

## Equations

The prognostic variables are the relative vorticity ζ, the divergence D, the temperature T and ln p<sub>s</sub> as spherical harmonics, and the tracers of the field table (at least the specific humidity `sphum`). Water vapour enters the dynamics only through the virtual temperature T<sub>v</sub> = T(1 + 0.608 q) in the geopotential and the pressure gradient, and only with `use_virtual_temperature = .true.` (the default is `.false.`).

The levels are bounded by the half-level pressures p<sub>k+1/2</sub> = `pk`(k) + `bk`(k) p<sub>s</sub>, from the top down ([`vert_coordinate_mod`](https://eddy-stanford.github.io/MiMA/api/vert_coordinate_mod/)). `vert_coord_option` chooses equally spaced sigma levels (`'even_sigma'`), sigma levels set by `surf_res`, `scale_heights` and `exponent` (`'uneven_sigma'`, the default), a transition from these to pressure levels between `p_sigma` and `p_press` (`'hybrid'`), or `pk` and `bk` from [`vert_coordinate_nml`](Parameters.md#vert_coordinate_nml) (`'input'`, at most 100 levels). The full-level pressures follow Simmons and Burridge, ln p<sub>k</sub> = ln p<sub>k+1/2</sub> − α<sub>k</sub>, and with a top at zero pressure ln p<sub>1</sub> = ln p<sub>3/2</sub> − 1 ([`press_and_geopot_mod`](https://eddy-stanford.github.io/MiMA/api/press_and_geopot_mod/)). The geopotential is summed up from the surface geopotential, and the vertical mass flux down from the top, where it is zero (equations 25-46 of the note).

In each layer the momentum equations are written as ∂u/∂t = A − (a cos θ)<sup>-1</sup> ∂E/∂λ and ∂v/∂t = B − a<sup>-1</sup> ∂E/∂θ, with E = Φ + (u² + v²)/2 and A, B the rotational, pressure-gradient (R T<sub>v</sub> ∇ln p) and vertical-advection terms; then ∂ζ/∂t = ∇·(B, −A) and ∂D/∂t = ∇·(A, B) − ∇²E (equations 55-66). The energy conversion κT ω/p in the temperature equation is differenced consistently with the pressure gradient, so that the scheme conserves energy and angular momentum.

## Spectral layout

A spectral field is a complex array `c(m, n)` with the zonal index m = 0, ..., `num_fourier` and the "meridional" index n = 0, ..., `num_spherical`. The zonal wavenumber is M = m × `fourier_inc` (`fourier_inc` > 1 gives a sector model with that symmetry in longitude) and the total wavenumber is L = M + n. Each processor holds the m of its block and all n, as `c(ms:me, 0:num_spherical)`.

**Truncation.** With `triang_trunc = .true.` (the default) the coefficients with L > `num_spherical` − 1 are zero, so T42 is `num_fourier = 42`, `num_spherical = 43`. The extra row is needed for the meridional derivatives, which couple n − 1 and n + 1. With `triang_trunc = .false.` (rhomboidal truncation) the row n = `num_spherical` is zero.

**Grid.** `lon_max` longitudes and `lat_max` Gaussian latitudes, ordered from south to north. Quadratic terms are alias-free when `lon_max` ≥ 3 `num_fourier` + 1 and 2 `lat_max` ≥ 3 (`num_spherical` − 1) + 1 (5 instead of 3 for rhomboidal truncation); the model stops otherwise. `lon_max` must be even with no prime factors other than 2, 3 and 5 (for the FFT). The usual triangular truncations are:

| | `num_fourier` | `num_spherical` | `lon_max` | `lat_max` |
|---|:---:|:---:|:---:|:---:|
| T21 | 21 | 22 | 64 | 32 |
| T42 (default) | 42 | 43 | 128 | 64 |
| T85 | 85 | 86 | 256 | 128 |
| T170 | 170 | 171 | 512 | 256 |
| T255 | 255 | 256 | 768 | 384 |

**Normalisation.** The associated Legendre functions satisfy ∫ P<sub>L</sub><sup>M</sup>(μ)² dμ = 1 over μ = sin θ from −1 to 1, so P<sub>0</sub><sup>0</sup> = 1/√2 ([`gauss_and_legendre_mod`](https://eddy-stanford.github.io/MiMA/api/gauss_and_legendre_mod/)), and a grid field is

x(λ, μ) = Σ<sub>n</sub> c(0, n) P(μ) + 2 Re Σ<sub>m>0</sub> Σ<sub>n</sub> c(m, n) P(μ) e<sup>iMλ</sup>,

with P = P<sub>L</sub><sup>M</sup> for each (m, n). The global mean of x is therefore c(0, 0)/√2 (adding a constant Δ to the grid field adds √2 Δ to `c(0, 0)`, which is how the corrections below change the spectral fields), and its global-mean square is ½ Σ<sub>n</sub> c(0, n)² + Σ<sub>m>0, n</sub> |c(m, n)|².

## Transforms and decomposition

A field goes from the grid to spherical harmonics by an FFT in longitude ([`grid_fourier_mod`](https://eddy-stanford.github.io/MiMA/api/grid_fourier_mod/), with Temperton's FFT in [`fft_mod`](https://eddy-stanford.github.io/MiMA/api/fft_mod/)), a transpose of the Fourier coefficients, and a Legendre transform by Gaussian quadrature ([`spherical_fourier_mod`](https://eddy-stanford.github.io/MiMA/api/spherical_fourier_mod/)), which uses the symmetry of the Legendre functions about the equator to treat each latitude with its mirror latitude. The inverse runs the other way. [`transforms_mod`](https://eddy-stanford.github.io/MiMA/api/transforms_mod/) combines them with the operators of [`spherical_mod`](https://eddy-stanford.github.io/MiMA/api/spherical_mod/): `uv_grid_from_vor_div` (u cos θ and v cos θ in spectral space, transformed and divided by cos θ), `vor_div_from_uv_grid` and `horizontal_advection` (−**v**·∇ξ on the grid for a spectral ξ).

The decomposition is set by [`spec_mpp_mod`](https://eddy-stanford.github.io/MiMA/api/spec_mpp_mod/) ([`spec_mpp_nml`](Parameters.md#spec_mpp_nml)):

* the grid is divided in latitude only (layout 1 × npes), so `lat_max` must be a multiple of the number of processes, and each process holds whole latitude circles;
* the spectral coefficients are divided in m only (layout npes × 1), in contiguous blocks. `spectral_decomposition = 'balanced'` (the default) gives each block about the same work, with a cost of `num_spherical` − m + `spectral_cost_per_m` for each m; `'equal'` gives each block the same number of m, as in v1. The results are the same either way, and a run may restart with either;
* the transposes between the two exchange the Fourier coefficients in a ring of point-to-point messages.

With `MIMA_OPENMP` the Legendre transforms, the FFTs and the grid-point loops of the dynamics are also split over the levels among the threads of each process ([OpenMP](GettingStarted.md#openmp)); the answers do not depend on the number of threads.

## Time step

Each atmospheric time step first computes the physics tendencies (`atmosphere_down` and `atmosphere_up` in [`atmosphere_mod`](https://eddy-stanford.github.io/MiMA/api/atmosphere_mod/)), then steps the dynamics with them in `num_steps` substeps (default 1). Each (sub)step of `spectral_dynamics` follows the outline on pages 16-17 of the note:

1. store the global mass, energy and water for the corrections;
2. compute the pressures, the geopotential and the gradient of p<sub>s</sub> (from ln p<sub>s</sub> in spectral space);
3. compute the pressure-gradient and energy-conversion terms, the vertical mass flux and the ln p<sub>s</sub> tendency (`four_in_one`);
4. add the vertical advection of u, v and T (`vert_advect_uv`, `vert_advect_t`), the horizontal advection of T and (f + ζ) × **v**, and form the spectral vorticity and divergence tendencies, with −∇²E;
5. apply the semi-implicit correction and the damping;
6. take the leapfrog step and transform the new ζ, D, u, v, T and p<sub>s</sub> to the grid;
7. step the tracers, then apply the corrections.

* **Leapfrog and Robert filter** ([`leapfrog_mod`](https://eddy-stanford.github.io/MiMA/api/leapfrog_mod/)): X(t+Δt) = X(t−Δt) + 2Δt G, and X(t) is filtered to X(t) + r [X(t−Δt) − 2X(t) + X(t+Δt)] with r = `robert_coeff` (default 0.03). Only two time levels are stored: the new level overwrites the old one. On the last substep the filter is applied in two parts, the part without X(t+Δt) before the step (`leapfrog_2level_A`) and r X(t+Δt) once the new level is final (`complete_robert_filter`, called by `atmosphere_up`). The first step of a cold start is a forward step of Δt.
* **Semi-implicit scheme** ([`implicit_mod`](https://eddy-stanford.github.io/MiMA/api/implicit_mod/), `use_implicit = .true.`): the gravity-wave terms of the divergence, temperature and ln p<sub>s</sub> equations are linearised about a state at rest with T = 300 K at all levels and p<sub>s</sub> = `reference_sea_level_press`, and treated with the weight `alpha_implicit` (0.5: centred; 1: backward). This gives one N × N matrix per total wavenumber L (N levels), inverted when the time step changes (equations 68-78).
* The physics receives `atmosphere_mod`'s own copy of the grid state at the two levels; the Robert filter changes only the dynamics' copy, so the physics steps from the unfiltered previous level.

## Damping and sponge

[`spectral_damping_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_damping_mod/) damps vorticity, divergence, temperature and the spectral tracers with a hyperdiffusion of rate

* ν(L) = `damping_coeff` × [L(L+1)/(L<sub>0</sub>(L<sub>0</sub>+1))]<sup>`damping_order`</sup> with L<sub>0</sub> = `num_spherical` − 1 (`damping_option = 'resolution_dependent'`, the default), so the smallest retained wave is damped at `damping_coeff` (default 1.157 × 10<sup>-4</sup> s<sup>-1</sup>, a tenth of a day) and `damping_order = 4` is a ∇<sup>8</sup> diffusion; or
* ν(L) = `damping_coeff` × [L(L+1)/a²]<sup>`damping_order`</sup> (`'resolution_independent'`).

Vorticity and divergence can have their own `damping_coeff_vor`, `damping_order_vor`, `damping_coeff_div` and `damping_order_div`. The damping is implicit over the leapfrog step: the tendency G becomes [G − ν X(t−Δt)]/(1 + 2Δt ν).

On the top level only, a ∇² sponge with rate `coeff` × L(L+1)/a² acts on the eddy (m > 0) vorticity and divergence (`eddy_sponge_coeff`), the zonal-mean vorticity, i.e. the zonal-mean zonal wind (`zmu_sponge_coeff`), and the zonal-mean divergence, i.e. the meridional wind (`zmv_sponge_coeff`). All three are 0 (off) by default; the Rayleigh sponge of the upper levels is in the physics ([`damping_driver_nml`](Parameters.md#damping_driver_nml)).

## Mass, energy and water corrections

ln p<sub>s</sub> as the prognostic variable cannot conserve mass exactly, so after each step [`spectral_dynamics_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_dynamics_mod/) can restore three global integrals (all on by default; the energy and water corrections need the mass correction):

* `do_mass_correction`: p<sub>s</sub> is multiplied by the ratio of its global mean at the previous level to that at the new level.
* `do_energy_correction`: a uniform ΔT is added so that the mass-weighted global mean of c<sub>p</sub>T + (u² + v²)/2 equals that of the previous level stepped with the physics tendencies alone. `tdt_tempcor` in the `dynamics_every` [diagnostics](Diagnostics.md#module-dynamics_every) is ΔT divided by the leapfrog step.
* `do_water_correction`: the specific humidity is multiplied by one factor where the pressure is at least `water_correction_limit` (default 200 hPa), chosen so that the global water equals that of the previous level stepped with the physics tendencies alone. Above the limit, in the stratosphere, it is not changed: correcting there would act as an artificial sink. This needs a grid-point humidity (`numerical_representation = 'grid'`), unless `water_correction_limit` ≤ 0.

The model stops if the temperature leaves `valid_range_t` (default 100 to 500 K), printing where.

Tracers are spectral or grid-point, as set by `numerical_representation` in the field table. Grid-point tracers (`sphum` in all the shipped configurations) are advected horizontally by a finite-volume van Leer scheme ([`fv_advection_mod`](https://eddy-stanford.github.io/MiMA/api/fv_advection_mod/), after Lin and Rood 1996) and are not diffused; spectral tracers are advected spectrally, diffused, and can have their negative values filled from their neighbours ([`water_borrowing_mod`](https://eddy-stanford.github.io/MiMA/api/water_borrowing_mod/), `hole_filling = 'on'`).

## Initial state

If the restart files `INPUT/spectral_dynamics.res.nc` and `INPUT/atmosphere.res.nc` exist, the model restarts from them. Otherwise [`spectral_init_cond`](https://eddy-stanford.github.io/MiMA/api/spectral_init_cond_mod/#spectral_init_cond) sets up a cold start:

* the vertical coordinate, as above;
* the surface geopotential (`topography_option`): `'interpolated'` (the default) from the file of `topography_nml` ([FMS namelists](Parameters.md#fms-namelists)), spectrally truncated and, with `ocean_topog_smoothing` > 0 (default 0.995), regularized over the ocean so that this fraction of its roughness is removed ([`topog_regularization_mod`](https://eddy-stanford.github.io/MiMA/api/topog_regularization_mod/), [`topog_regularization_nml`](Parameters.md#topog_regularization_nml)); `'flat'`; `'gaussian'` mountains from `gaussian_topog_nml`; or `'input'`, `zsurf` from `INPUT/topography.data.nc`;
* an isothermal atmosphere at rest at `initial_temperature` ([`spectral_init_cond_nml`](Parameters.md#spectral_init_cond_nml), default 264 K), with the surface pressure in hydrostatic balance with the topography (`reference_sea_level_press` at height 0), a vorticity perturbation of 10<sup>-7</sup> s<sup>-1</sup> in four coefficients on the lowest three levels, and `sphum = initial_sphum` (default 2 × 10<sup>-6</sup> kg/kg) ([`spectral_initialize_fields_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_initialize_fields_mod/));
* or the state of `INPUT/initial_conditions.nc` ([Specified initial conditions](Configurations.md#specified-initial-conditions)).

Noise: see [Adding noise to the initial conditions](Configurations.md#adding-noise-to-the-initial-conditions).

## What MiMA changes relative to the GFDL note

* **No adjustment step.** The note's option of passing the new temperature and humidity to a second physics routine after the dynamics has been removed (with `spectral_physics.res.nc`, see [Migration to v2](Migration_v2.md)). The Robert filter is still applied in two parts, but nothing happens between them, so it is the standard filter.
* **Water correction.** The note's correction rescales the humidity everywhere; MiMA's rescales it only below `water_correction_limit`.
* **Damping.** v2 uses `damping_order = 4` (∇<sup>8</sup>; v1 used 2) and has the separate vorticity and divergence parameters and the top-level sponges described above, which the note does not mention.
* **Hole filling** takes the tracer from the neighbours in longitude and in the vertical (the note says latitude and height).
* **Topography** over the ocean is regularized by a conjugate-gradient solver that was rewritten for v2, so that it also works at high resolution (see [Migration to v2](Migration_v2.md)).
* **Parallel layout.** The `'balanced'` division of the m among the processes and the OpenMP threading are v2 additions; the FFT is always Temperton's (the vendor-library and NAG versions are gone).
* **Initial states** can be read from a file (`specify_initial_conditions`) or perturbed with noise (`add_noise`).
* **Vertical coordinates**: `'mcm'` and `'v197'` have been removed, and `'simmons_and_burridge'` is the only vertical differencing.
* **Diagnostics**: `dynamics_every` adds eddy fluxes (`upvp`, `upwp`, `vpTp`) and two residual-mean streamfunctions: `psi_star` = [∫v dp] − [v′θ′]/∂<sub>p</sub>[θ], with θ here the potential temperature and [ ] the zonal mean, and `psi_dwc`, the pressure integral of v* from the downward-control balance −f̂ v* = (a cos²φ)<sup>-1</sup> ∂(cos²φ F<sub>φ</sub>)/∂φ + ∂F<sub>p</sub>/∂p of the TEM equations (Andrews, Holton and Leovy 1987; φ is the latitude).

The note is the GFDL documentation of the FMS spectral core, which MiMA v1 kept as `src/atmos_spectral/documentation/spectral_core.pdf`. It refers to companion notes on the barotropic vorticity and shallow-water models, which are not included.

References:

* Simmons, A. J., and D. M. Burridge, 1981: An energy and angular-momentum conserving vertical finite-difference scheme and hybrid vertical coordinates. [Mon. Wea. Rev., 109, 758-766](https://doi.org/10.1175/1520-0493(1981)109<0758:AEAAMC>2.0.CO;2).
* Lin, S.-J., and R. B. Rood, 1996: Multidimensional flux-form semi-Lagrangian transport schemes. [Mon. Wea. Rev., 124, 2046-2069](https://doi.org/10.1175/1520-0493(1996)124<2046:MFFSLT>2.0.CO;2).
* Andrews, D. G., J. R. Holton, and C. B. Leovy, 1987: Middle Atmosphere Dynamics. Academic Press.
