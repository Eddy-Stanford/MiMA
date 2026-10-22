# Model configurations

This page describes some common ways of changing the model setup beyond the default [test case](GettingStarted.md#the-test-case). All options are set in `input.nml`; see [Parameter settings](Parameters.md) for the full list of namelist variables.

* [Radiation options](#radiation-options)
* [Held-Suarez forcing](#held-suarez-forcing)
* [Surface options](#surface-options)
* [Specified initial conditions](#specified-initial-conditions)
* [Adding noise to the initial conditions](#adding-noise-to-the-initial-conditions)

## Radiation options

The radiation scheme is chosen with `radiation_scheme` in [`radiation_nml`](Parameters.md#radiation_nml):

```fortran
&radiation_nml
    radiation_scheme = 'rrtm' /
```

* `'rrtm'` (default): RRTMG clear-sky radiation, configured with [`rrtm_radiation_nml`](Parameters.md#rrtm_radiation_nml) and [`astro_nml`](Parameters.md#astro_nml).
* `'gray'`: the gray radiation scheme of Dargan Frierson ([Frierson, Held, Zurita-Gotor, JAS (2006)](https://doi.org/10.1175/JAS3753.1)), configured with [`gray_radiation_nml`](Parameters.md#gray_radiation_nml).
* `'none'`: no radiative heating and no radiative surface fluxes.

In the code, [`radiation_mod`](https://eddy-stanford.github.io/MiMA/api/radiation_mod/) reads `radiation_nml` and passes the physics driver's calls on to the chosen scheme: [`rrtm_radiation`](https://eddy-stanford.github.io/MiMA/api/rrtm_radiation/) (which calls RRTMG, with the solar zenith angle from [`rrtm_astro`](https://eddy-stanford.github.io/MiMA/api/rrtm_astro/)) or [`gray_radiation_mod`](https://eddy-stanford.github.io/MiMA/api/gray_radiation_mod/).
The two schemes provide different diagnostics (the shared ones, such as `olr` and `tdt_rad`, have the same names); see the `radiation` module in [Diagnostics](Diagnostics.md#module-radiation). A complete gray-radiation setup is in `examples/gray/`: its `input.nml` is the default test case with `radiation_scheme = 'gray'` and `gray_radiation_nml` in place of `rrtm_radiation_nml` and `astro_nml`, and its `diag_table` writes daily and 30-day means including the gray radiative fluxes. It needs `navy_topography.data.nc` and `navy_pctwater.data.nc` from `input/INPUT/` (see [Running the examples](GettingStarted.md#running-the-examples)).

## Held-Suarez forcing

MiMA can run the [Held and Suarez (1994)](https://doi.org/10.1175/1520-0477(1994)075<1825:APFTIO>2.0.CO;2) idealized forcing: Newtonian relaxation of temperature towards a zonally symmetric equilibrium profile, and Rayleigh friction of the winds in the boundary layer. It is switched on in [`physics_driver_nml`](Parameters.md#physics_driver_nml), and its parameters are set in [`held_suarez_nml`](Parameters.md#held_suarez_nml) (the defaults are the HS94 values). The forcing is computed by [`held_suarez_forcing`](https://eddy-stanford.github.io/MiMA/api/held_suarez_mod/#held_suarez_forcing) in [`held_suarez_mod`](https://eddy-stanford.github.io/MiMA/api/held_suarez_mod/), which the physics driver calls when `do_held_suarez = .true.`:

 Variable | Default | Meaning
 :--- | :---: | :---
 `t_zero`, `t_strat` | 315, 200 K | surface equilibrium temperature at the equator, minimum equilibrium temperature
 `delh`, `delv` | 60, 10 K | equator-to-pole temperature difference, vertical potential temperature difference
 `sigma_b` | 0.7 | top of the frictional boundary layer (sigma)
 `ka`, `ks`, `kf` | 40, 4, 1 days | free-atmosphere and surface relaxation times, boundary-layer friction time
 `do_rayleigh_friction` | `.true.` | apply the boundary-layer friction
 `do_conserve_energy` | `.false.` | heat the air by the frictional dissipation

The forcing depends only on latitude and sigma, so it works at any horizontal resolution and with any vertical levels.

For the standard **dry** benchmark, switch off radiation, moist physics and the boundary layer, which also leaves the surface state unchanged. A complete example (T42, 20 evenly spaced sigma levels, flat topography) is in `examples/held_suarez/`; it needs a `RESTART/` directory but no `INPUT/` data files (`land_model_nml land_mask = 'none'`, so no land-sea mask file; see [Running the examples](GettingStarted.md#running-the-examples)):

```fortran
&radiation_nml
    radiation_scheme = 'none' /

&physics_driver_nml
    do_held_suarez    = .true.,
    do_boundary_layer = .false.,
    do_moist_physics  = .false.,
    do_damping        = .false. /
```

The model always carries a humidity tracer (`sphum` in the field table); in the dry setup it stays zero. Also set `use_virtual_temperature = .false.` and `do_water_correction = .false.` in [`spectral_dynamics_nml`](Parameters.md#spectral_dynamics_nml).

The HS forcing can be combined with other parts of the model:

* `do_damping = .true.` with [`damping_driver_nml`](Parameters.md#damping_driver_nml) enables the Rayleigh sponge (`do_rayleigh`) and/or the convective gravity-wave drag (`do_cg_drag`, [`cg_drag_mod`](https://eddy-stanford.github.io/MiMA/api/cg_drag_mod/)); see [`damping_driver_mod`](https://eddy-stanford.github.io/MiMA/api/damping_driver_mod/). Note that `do_cg_drag` defaults to `.true.`, so set `do_cg_drag = .false.` if you want only the sponge.
* Non-flat topography through `topography_option` in `spectral_dynamics_nml`.
* **Moist variants:** with `do_moist_physics = .true.` and `do_boundary_layer = .true.` (and `do_rayleigh_friction = .false.`), the HS temperature relaxation replaces radiation while MiMA's moist physics, boundary layer and surface fluxes stay active. This is similar in spirit to the moist Held-Suarez test of [Thatcher and Jablonowski (2016)](https://doi.org/10.5194/gmd-9-1263-2016), but uses MiMA's own boundary-layer and surface schemes.
  With `radiation_scheme = 'none'` the surface receives no radiation, so hold the surface temperature fixed with `ocean_model = 'fixed'` in [`ocean_model_nml`](Parameters.md#ocean_model_nml) and `land_model = 'fixed'` in [`land_model_nml`](Parameters.md#land_model_nml) (see [Fixed surface](#fixed-surface); the initial profile is set by `Tm` and `deltaT`); a slab ocean would otherwise cool without limit. The HS equilibrium temperature near the equatorial surface (315 K) is warmer than typical SSTs, so the lowest layers are stably stratified over the ocean and the hydrological cycle is weak: with the default SST profile (about 298 K at the equator) precipitation takes about three weeks to start and settles near 0.5 mm/day in the global mean.

## Surface options

The surface is made of an ocean and a land component, set in [`ocean_model_nml`](Parameters.md#ocean_model_nml) and [`land_model_nml`](Parameters.md#land_model_nml). The surface albedo is set in [`surface_albedo_nml`](Parameters.md#surface_albedo_nml), the prescribed ocean heat flux convergence in [`qflux_nml`](Parameters.md#qflux_nml) and the reference heights of the 2 m and 10 m diagnostics in [`flux_exchange_nml`](Parameters.md#flux_exchange_nml). The defaults are the standard test case, a slab ocean and a saturated slab land with the Navy land-sea mask. Coming from v1, see the [mapping of the old `&simple_surface_nml`](Migration_v2.md#old-to-new-namelist-mapping).

* [Slab ocean](#slab-ocean)
* [Prescribed SST with interactive land](#prescribed-sst-with-interactive-land)
* [Fixed surface](#fixed-surface)
* [Q-flux from a file](#q-flux-from-a-file)
* [Bucket land and snow](#bucket-land-and-snow)
* [Land-sea mask](#land-sea-mask)

| | `'slab'` (ocean), `'simple'` (land) | `'prescribed'` | `'fixed'` |
|---|---|---|---|
| `ocean_model` | temperature stepped by the surface energy budget (plus the Q-flux over the ocean) | temperature set to the SST file at each step | initial temperature kept (`sst_init`) |
| `land_model` | the slab (optionally with bucket and snow) | temperature set to `t_surf_file` | initial temperature kept (`t_init`) |

### Slab ocean

The default. Each ocean point is a slab of heat capacity `heat_capacity` [J/m2/K] that is heated by the net surface flux and by the Q-flux. The initial temperature is set by `sst_init` (`'profile'` by default: `Tm - deltaT*(3 sin^2(lat) - 1)/3`, which has the area mean `Tm` and an equator-to-pole difference `deltaT`), and a restart file always takes precedence. `sst_init` also takes `'uniform'` (`Tm` everywhere), `'file'` (from `sst_file`) and the Aqua-Planet Experiment profiles `'ape_control'` and `'ape_qobs'` ([Neale and Hoskins, 2000](https://doi.org/10.1006/asle.2000.0019)). `sst_offset` [K] (default 0) is added to the cold-start SST of every `sst_init`, for example `sst_init = 'ape_qobs'`, `sst_offset = 5.` is the "qobs+5" experiment. The heat capacity is `trop_capacity` equatorward of `trop_cap_limit`, `heat_capacity` poleward of `heat_cap_limit` and linear in latitude in between; in the Northern Hemisphere the `heat_capacity` end is multiplied by `np_cap_factor`. With `trop_capacity = heat_capacity` and `np_cap_factor = 1.` it is uniform. As a guide, `3.e8` J/m2/K is a water layer about 72 m deep and `1.e8` about 24 m. The land has its own `heat_capacity` in `&land_model_nml`.

```fortran
&ocean_model_nml
    ocean_model   = 'slab',
    sst_init      = 'profile',
    Tm            = 285.,
    sst_offset    = 0.,
    deltaT        = 40.,
    heat_capacity = 3.e08,
    trop_capacity = 1.e08,
    trop_cap_limit = 20.,
    heat_cap_limit = 60.,
    np_cap_factor = 1. /

&qflux_nml
    qflux_source      = 'analytic',   ! 'analytic', 'file' or 'none'
    do_warmpool       = .true.,       ! the zonally asymmetric warm-pool and Gulf Stream terms
    qflux_ocean_only  = .true.,       ! no Q-flux over land
    qflux_remove_mean = .true. /      ! subtract the area-weighted ocean mean
```

The Q-flux options:

* `qflux_source = 'analytic'`: the zonally symmetric meridional Q-flux (`qflux_amp`, `qflux_width`) and, with `do_warmpool = .true.`, the zonally asymmetric terms of [Garfinkel et al. (2020)](https://doi.org/10.1175/JCLI-D-19-0181.1) (`warmpool_*`, `gulf_*`, `kuroshio_amp`, ...).
* The zonally asymmetric terms are set by `warmpool_pattern`: `'garfinkel2020'` (default) is the published pattern of Garfinkel et al. (2020), a sum of 79 labelled components (waves and Gaussian blobs: warm pool, equatorial cold tongue, SPCZ, ITCZ, Kuroshio, Gulf Stream, Agulhas, Peru upwelling, Southern Ocean, Norwegian Sea, ...) listed with their centres, widths and scaling knobs in [`ocean_qflux_mod`](https://eddy-stanford.github.io/MiMA/api/ocean_qflux_mod/); `'warm_pool'` is only the tropical warm pool `warmpool_amp*(1 - lat^4)*cos(warmpool_k*(lon - warmpool_phase))`, restricted to one wavelength. The sensitivity knobs (`Hawaiiextra`, `Pac_SPCZextra`, `Pac_ITCZextra`, `Africaextra`, `Sampeextra`) mostly move heat between regions: their net global effect is below 0.001 W/m2 per unit. Several regional terms of `'garfinkel2020'` scale with `qflux_amp` (and `qflux_amp + warmpool_amp`), also when `qflux_source = 'none'`.
* `qflux_source = 'file'`: read from a file, see [Q-flux from a file](#q-flux-from-a-file).
* `qflux_source = 'none'`: no base Q-flux. `do_warmpool = .true.` adds the warm-pool and related terms to **any** `qflux_source` (`'none'`, `'analytic'`, `'file'`), as in v1; masking and mean removal apply to the sum.
* `qflux_ocean_only` and `qflux_remove_mean` (both `.true.` by default) set the Q-flux to zero over land and subtract its area-weighted ocean mean (over the whole sphere with `land_mask = 'none'`), so that it only redistributes heat. Setting both to `.false.` adds the Q-flux at every point exactly as MiMA v1 did, which is needed to reproduce setups published with v1 (see [Reproducing the v1 surface](Migration_v2.md#reproducing-the-v1-surface) and [Changes to the answers](Migration_v2.md#changes-to-the-answers)). With `qflux_ocean_only = .false.` the land heat budget receives the Q-flux as well; by default it never does. The statistics of the raw and the applied Q-flux are written to the log file at start-up.

The Q-flux is written as `ocean/qflux` (v1's output name `oflx` can be kept: `"ocean", "qflux", "oflx", ...`).

### Prescribed SST with interactive land

An AMIP-style run: the SST over the ocean follows a file, while the land keeps its own slab energy budget (or the bucket below). Put the file in `INPUT/` and name it in `sst_file`, without the `.nc`:

```fortran
&ocean_model_nml
    ocean_model  = 'prescribed',
    sst_file     = 'sst_amip',      ! INPUT/sst_amip.nc
    sst_file_var = 'sst' /          ! name of the variable in the file; default = the file name

&land_model_nml
    land_model = 'simple' /         ! the default: interactive land
```

* **Fill the file over land.** The file is regridded to the model grid with the conservative interpolator. If the file grid is not the model grid, an ocean cell next to the coast mixes in the file's land values (usually missing values or zeros), which gives wrong SSTs along every coast. Extrapolate the SST over land (nearest ocean value is enough) before using the file. The ocean does not use the values on land and, with `land_model = 'simple'`, the land steps its own temperature. But on a cold start with `t_init = 'ocean'` (the default) the land **starts** from the file's land values, so a fill value there stops the model (the initial surface temperature must be within 150-400 K); fill the land, or set `t_init = 'file'` with a `t_surf_file`.
* **Units and variable.** Kelvin, a 2-D (lon, lat) field with a time axis. `sst_file_var` names the variable; the default is the file name, which is v1's convention. MiMA reads only that variable, so extra variables in the file (e.g. `time_bnds`) do no harm.
* **Records.** One record is a constant SST. Two to twelve records are a climatology that repeats every year, whatever the year of the time axis. More than twelve records are a dated series, and the run must lie inside its time span. Time is read from the `time` axis (`units` as `days since 0001-01-01 00:00:00`).
* **Calendars.** The model calendar is set by `calendar` in [`coupler_nml`](Parameters.md#coupler_nml) (`'thirty_day'` by default, `'julian'` or `'noleap'`), and the `calendar` attribute of the file's time axis must match it: `360_day` (or `360`) for `thirty_day`, `julian` for `julian`, `noleap` for `noleap`. A `julian` file is also converted for a `noleap` model, and a `noleap` file for a `julian` model; any other combination is a fatal error. A time axis with the base year 0000 (`days since 0000-01-01`) is a climatology for all years and has no calendar restriction. MiMA's own output files are `360_day` and can be used as they are in a `thirty_day` run.
* `sst_offset` [K] is added to the file's SST at every step (and to the cold-start SST), for example `sst_offset = 4.` for an AMIP+4K experiment. The default 0 leaves the SST as in the file. The diagnostic `ocean/sst_interp` is the file's SST without the offset.
* The temperature at a time step is set to the file's value; at start-up `sst_init` is forced to `'file'`, so a cold start begins from the file.
* With the default `qflux_ocean_only = .true.` the Q-flux is not used in this configuration: the ocean temperature is imposed and the land receives none. (With `qflux_ocean_only = .false.` the land would be heated by it.)
* The surface component writes `flux/t_surf` and, because the SST is imposed, [`ocean/qflux_implied`](#q-flux-from-a-file). No special treatment of sea ice: the SST file and the background albedo apply everywhere.

To reproduce the v1 behaviour, where the file SST was also imposed on land points, set `land_model = 'prescribed'` and `t_surf_file` to the same file (`t_surf_file_var` defaults to the file name).

### Fixed surface

The surface temperature stays at its initial value. This is the surface of the moist Held-Suarez configuration ([Held-Suarez forcing](#held-suarez-forcing)), and of any run in which the surface should not respond:

```fortran
&ocean_model_nml
    ocean_model = 'fixed',
    sst_init    = 'profile',       ! or 'uniform' (Tm), 'file', 'ape_control', ...
    Tm          = 285.,
    deltaT      = 40. /

&land_model_nml
    land_model = 'fixed' /
```

The initial temperature of the land is the ocean's (`t_init = 'ocean'`), so there is no land-sea contrast unless you set `t_init` or `sst_init` otherwise. The surface fluxes are still computed; only the temperature is held. `sst_init` also takes the APE profiles (`'ape_control'`, `'ape_qobs'`), `'uniform'` (temperature `Tm`) and `'file'`, and `sst_offset` [K] is added to the initial temperature, which a fixed ocean keeps.

A complete aquaplanet is in `examples/aquaplanet/`: the Aqua-Planet Experiment "qobs" case with `land_mask = 'none'`, flat topography, `ocean_model = 'fixed'` and `sst_init = 'ape_qobs'`, perpetual equinox (`obliq = 0.` in `astro_nml`), `solr_cnst = 1365.` and 348 ppmv CO<sub>2</sub>, and otherwise the physics of the test case. Unlike the APE protocol it has no clouds and uses the seasonal `ozone_1990`. `sst_offset = 4.` gives the qobs+4K experiment and `sst_init = 'ape_control'` the control SST.

### Q-flux from a file

A slab ocean can read its Q-flux from a file instead of the analytic formula. The usual use is to obtain the Q-flux that keeps a slab ocean at a given SST climatology, so that a slab-ocean run reproduces the SST of a prescribed-SST run:

1. **Run with prescribed SST.** Use `ocean_model = 'prescribed'` (see above) with the SST climatology and `heat_capacity`, `trop_capacity`, ... equal to the values of the slab run you want (the diagnostic uses the heat capacity profile of `ocean_model_nml`). Request the diagnostic in the `diag_table`, as a monthly or a time mean:

    ```
    "atmos_monthly", 30, "days", 1, "days", "time",
    "ocean", "qflux_implied", "qflux_implied", "atmos_monthly", "all", .true., "none", 2,
    ```

2. **Average the diagnostic** over the years of the run (after the spin-up of the land), for example with `ncra` or `cdo ymonmean`: one record for a time-mean Q-flux, or 12 monthly records for a seasonal cycle. The file may be a MiMA output file as it is.
3. **Use it in the slab run.**

    ```fortran
    &ocean_model_nml
        ocean_model = 'slab' /

    &qflux_nml
        qflux_source      = 'file',
        qflux_file        = 'qflux_implied',   ! INPUT/qflux_implied.nc
        qflux_file_var    = 'qflux_implied',   ! variable name in the file (default 'qflux')
        do_warmpool       = .false.,           ! .true. would add the warm-pool terms to the file
        qflux_ocean_only  = .true.,
        qflux_remove_mean = .true. /
    ```

What to know:

* `qflux_implied` is the ocean heat flux convergence `Q = cap*(dT/dt - extra_heating) - F_net`, where `F_net` is the net surface flux after the implicit step (longwave, shortwave, snow melt, sensible and latent heat). It is **zero over land**, and the **local surface heating** of [`local_heating_nml`](Parameters.md#local_heating_nml) is excluded, so a later slab run that has the same local heating does not count it twice. A slab with `Q = qflux_implied` reproduces the SST change that the prescribed run imposed.
* When the file is read, the Q-flux is **masked to the ocean points and its area-weighted ocean mean is removed**, whatever its source (the two logicals of `&qflux_nml`). The mean of `qflux_implied` over the ocean is not zero in general (it is the heat the prescribed SST takes up or loses), and removing it is what keeps the slab ocean from drifting in global mean. Set both logicals `.false.` to use the file as it is. With the means removed, the slab run reproduces the SST pattern, not the global-mean SST of the prescribed run, which is set by the energy balance and the initial condition.
* **For an exact round trip of the SST** (a slab run that follows the prescribed SST), derive the Q-flux from a **spun-up** prescribed run and set `qflux_remove_mean = .false.` (and `do_warmpool = .false.`). With the default `qflux_remove_mean = .true.` the slab drifts by the prescribed run's net surface imbalance (its global mean), which is intended: a Q-flux that only redistributes heat cannot hold a global mean. In the 10-day acceptance test, the applied Q-flux equalled `qflux_implied` minus its ocean mean to 1e-11 W/m2, the SST drift of the slab run differed from that of the prescribed run by `Q*t/cap` with a correlation of 0.999, and MiMA diagnostic output was read directly as `qflux_file`. The test's drift ratio (0.84) was limited by the prescribed run starting cold, whose spin-up transient does not apply to the following days.
* `do_warmpool = .true.` adds the warm-pool terms to a file Q-flux as it does to any other `qflux_source`.
* The file conventions are those of the SST file: one record is constant, 2-12 records are cyclic over the year, more than 12 are a dated series. The variable name is `qflux_file_var` (default `'qflux'`), the grid may differ from the model grid (conservative regridding; the file is zero over land, so coastal cells are not corrupted by missing values).
* A prescribed-SST run is not in equilibrium: the SST file fixes T, so `qflux_implied` includes the heat that a free slab would exchange with the atmosphere, and the first days of a run that starts cold are a spin-up transient (a mean of 30 W/m2 over days 0 to 5 in the test). Average over several years of the prescribed run, and judge the resulting slab run by its SST drift.

### Bucket land and snow

By default the land is saturated: evaporation proceeds as over the ocean, subject to the roughness settings. With `hydrology = 'bucket'` the land has a water store (a Manabe bucket), and with `do_snow = .true.` also a snow store with melt, sublimation, frost and a snow albedo. The settings:

| Variable | Meaning | Default |
|---|---|---|
| `hydrology` | `'saturated'` or `'bucket'` | `'saturated'` |
| `water_capacity` | bucket capacity [kg/m2] | 150. |
| `beta_frac` | fraction of the capacity above which evaporation is not limited | 0.75 |
| `water_init_frac` | initial water as a fraction of the capacity | 1.0 |
| `do_snow` | store snow, melt it, sublimate it, brighten the albedo (requires the bucket and `land_model = 'simple'`) | `.false.` |
| `albedo_snow` | albedo of deep snow | 0.60 |
| `snow_mask_depth` | snow mass at which the albedo is half-way to `albedo_snow` [kg/m2] | 10. |
| `snow_max` | snow mass above which the excess runs off [kg/m2] | 1000. |

Rain fills the bucket and evaporation empties it; what exceeds `water_capacity` runs off. The evaporation is multiplied by `beta`, which is 1 with snow on the ground and otherwise `min(W/(beta_frac*water_capacity), 1)` for the soil water `W`; the factor applies to the evaporation, its derivatives and `q_star`, so the implicit coupling and the reference-height diagnostics see it. Without `do_snow`, snowfall melts on arrival (the surface takes `hlf*fprec`) and joins the bucket. With `do_snow`, snowfall is stored as snow of mass `S`. Each step the snow sublimates first (at the evaporation rate, taking the latent heat of fusion from the surface), takes up frost where the evaporation is negative (releasing it again), and melts into the bucket where the surface is above freezing, limited by the heat above freezing and by the snow. Snow above `snow_max` runs off (a crude glacier). The snow brightens the albedo as `S/(S + snow_mask_depth)` toward `albedo_snow`, never below the background albedo, and **the radiation sees it**: the albedo passed to RRTM and the gray scheme, and the diagnostic `flux/albedo`, include the snow.

**Where the snow comes from.** The snowfall is that of the moist physics, set by `precip_phase` in [`moist_processes_nml`](Parameters.md#moist_processes_nml): with the default `'temperature'` it snows in the columns whose lowest 150 hPa are at or below 273.16 K ([Changes to the answers](Migration_v2.md#changes-to-the-answers)). The phase aloft uses `tfreeze` of `&moist_processes_nml`, whereas the melt at the surface uses the constant `TFREEZE` of `constants_mod`. With `precip_phase = 'liquid'` nothing is frozen, and snow appears only from a restart or from frost. The same frozen precipitation also falls on the ocean and on the default saturated land, where it melts on arrival and cools the surface by `hlf` per kg.

**Required settings.** The legacy moisture-roughness throttle of the default land (`rough_moist_factor = 1.e-12`, `rough_moist_profile = 'legacy_tropical'`) suppresses evaporation over land by design, and would suppress the bucket's evaporation twice. The model therefore stops if the bucket is switched on together with it. Set `rough_moist_factor = 1.` and `rough_moist_profile = 'uniform'`:

```fortran
&land_model_nml
    land_model            = 'simple',
    hydrology             = 'bucket',
    rough_moist_factor    = 1.,
    rough_moist_profile   = 'uniform',
    water_capacity        = 150.,
    beta_frac             = 0.75,
    water_init_frac       = 1.0,
    do_snow               = .true.,
    albedo_snow           = 0.60,
    snow_mask_depth       = 10.,
    snow_max              = 1000. /
```

The model also stops for `do_snow` without `hydrology = 'bucket'` or with `land_model /= 'simple'` (the melt, sublimation and frost change the surface temperature).

Diagnostics: bucket: `land/water`, `land/beta`, `land/runoff`, `land/water_clip`; snow: `land/snow`, `land/melt`, `land/sublimation`, `land/frost`; the rates are in kg/m2/s. `water_clip` is the water that had to be added to keep the bucket from going negative; it is the only place where water is not conserved. `water` and `snow` are saved in `RESTART/land_model.res.nc`.

**Status: verified for conservation, not climate-validated.** In short runs (days to weeks) the water budget closes to 1e-13 kg/m2 and the land energy budget to 1e-9 W/m2, and the snow albedo and the restart round trip were checked. The slow processes, such as the spin-up of the soil moisture over years or the snow accumulation on the polar ice sheets (bounded here by `snow_max`), have not been examined. `surface_flux` still uses the saturated humidity at the surface for the buoyancy, an approximation over dry land. Check the water budget and the snow before relying on it for a long run.

A complete setup is in `examples/bucket_snow/`: the default test case with the bucket and snow switched on as above, and a `diag_table` with the soil water, the snow and their budget terms (see [Running the examples](GettingStarted.md#running-the-examples)).

### Land-sea mask

`land_mask` in [`land_model_nml`](Parameters.md#land_model_nml) chooses where the land is:

| `land_mask` | Land is | Needs |
|---|---|---|
| `'navy'` (default) | the Navy land-sea mask, interpolated to the model grid: land where at least 50 % of the cell is land | `INPUT/navy_pctwater.data.nc` |
| `'file'` | where the variable `land_mask_file` of the file `INPUT/<land_mask_file>.nc` (default `'lmask'`, so `INPUT/lmask.nc`; the variable has the same name as the file) is positive, on the model grid | the file; the shipped `input/INPUT/lmask.nc` is an Earth land-sea mask (0 or 1) for T42 (128 × 64), so `land_mask = 'file'` works at T42 without further files |
| `'boxes'` | inside the longitude-latitude boxes `slandlon`, `slandlat`, `elandlon`, `elandlat` (degrees; up to 10 boxes) | |
| `'height'` | where the surface height (from the topography) exceeds `land_height_limit` (default 10, **metres**) | |
| `'none'` | nowhere: an aquaplanet. The restart files are still written. | |

The mask decides which component owns each point, and the land roughness, land heat capacity and the bucket apply only on land. The ocean `heat_capacity` profile is used on the ocean points, and `&land_model_nml heat_capacity` on the land. With `land_mask = 'none'` the `land_model` is irrelevant (the model prints a note if it is not `'simple'`).

## Specified initial conditions

Without restart files, MiMA starts from an isothermal atmosphere at rest with a small vorticity perturbation. To start from your own initial state instead, set this flag in [`spectral_dynamics_nml`](Parameters.md#spectral_dynamics_nml):

```fortran
specify_initial_conditions = .true.
```

Then provide a netCDF file named `initial_conditions.nc` in the `INPUT/` directory where the model runs. It must contain zonal wind, meridional wind, temperature, specific humidity and surface pressure (`ucomp`, `vcomp`, `temp`, `sphum` and `ps`) at the model resolution. The file is only read on a cold start: if restart files are present in `INPUT/`, the model restarts from them instead. A zonally symmetric initial state stays zonally symmetric unless it is perturbed, e.g. with [noise](#adding-noise-to-the-initial-conditions). The cold-start fields are set up by [`spectral_init_cond`](https://eddy-stanford.github.io/MiMA/api/spectral_init_cond_mod/#spectral_init_cond), which calls [`spectral_initialize_fields`](https://eddy-stanford.github.io/MiMA/api/spectral_initialize_fields_mod/#spectral_initialize_fields) (`choice_of_init = 3` reads the file).

## Adding noise to the initial conditions

MiMA can add a random perturbation to the temperature, and optionally to the vorticity, when the model starts. This breaks the symmetry of an idealized initial state (a zonally symmetric state stays zonally symmetric otherwise) and makes the members of an ensemble differ. It is set in [`initial_noise_nml`](Parameters.md#initial_noise_nml) and applied by [`initial_noise_mod`](https://eddy-stanford.github.io/MiMA/api/initial_noise_mod/) after the restart files have been read or the cold start has been set up:

```fortran
&initial_noise_nml
    noise_type  = 'white',
    amplitude_t = 0.5,
    seed        = 1 /
```

With the default `noise_type = 'none'` nothing is added.

**What is added.** The perturbation is built in spherical-harmonic space and limited to the total wavenumbers L from `l_min` to `l_max` (default 1 to 21) and to the zonal wavenumbers from `m_min` (default 1, which leaves the zonal mean untouched; 0 also perturbs it). The global mean (L = 0) and the last total wavenumber of the truncation (L = 42 at T42), where noise would be unstable, are never perturbed, and the truncation is respected.

* `'white'`: every spherical harmonic in the band gets the same variance, with random phases. The field is isotropic: its statistics do not depend on the position or the direction on the sphere. Its variance per total wavenumber grows like 2L + 1, so most of it is in the smallest scales of the band.
* `'spectral'`: an isotropic field with a chosen variance per total wavenumber V(L). With `spectral_shape = 'power'` (default), V(L) ~ L<sup>-`spectral_slope`</sup>: `spectral_slope = 3` (default) is red noise like the enstrophy-cascade range, 0 is a flat spectrum in L. With `spectral_shape = 'gaussian'`, V(L) ~ exp(-(L - `l_peak`)<sup>2</sup>/(2 `l_width`<sup>2</sup>)) puts the perturbation at one scale.
* `'smooth_grid'`: smooth gradient noise (Perlin noise) on the grid, summed over `octaves` (default 4) whose amplitude falls by `gain` (default 0.5) from one to the next while the scale halves; `base_wavenumber` (default 4) is the approximate total wavenumber of the first octave. The noise is evaluated in three dimensions at the position of each grid point on the unit sphere, so it has no seam at the date line and no artefacts at the poles. It is then transformed to spectral space and limited to the band like the other types. With `radius` > 0 [deg] it is multiplied by a Gaussian blob of that e-folding great-circle distance around `lon0`, `lat0`, to perturb one region only (with `m_min = 1` the zonal mean of the blob is removed, which leaves a weak opposite anomaly along its latitude band; use `m_min = 0` to keep it local). It needs the whole sphere (`fourier_inc = 1`).

**Amplitude.** `amplitude_t` [K] is the global grid-point rms of the temperature perturbation on each level, exactly (to rounding): the perturbation of each level is rescaled after it is drawn. Peak values are typically 3 to 4 times larger. The model prints the rms it obtained on each level, top to bottom, in a line `initial noise: grid rms per level of temperature [K]:`. With `perturb_vorticity = .true.` the vorticity is perturbed in the same way, with the rms `amplitude_vor` [1/s] (default 10<sup>-6</sup> s<sup>-1</sup>), and the winds are recomputed from it (the divergence is not perturbed). The surface pressure and the humidity are not perturbed.

**Vertical structure.** With `vertical_structure = 'correlated'` (default), the patterns of two levels have the correlation exp(-|ln p<sub>k</sub> - ln p<sub>j</sub>|/`vertical_scale`), with `vertical_scale` in scale heights (default 1) and p the pressures of the levels for a surface pressure of `reference_sea_level_press`. This avoids putting much of the noise into grid-scale vertical structures, which excite fast gravity waves. `'independent'` draws an independent pattern on each level, and `'column'` puts the same pattern on every level. In all cases the amplitude is multiplied by a taper that is 1 between `p_top` and `p_bottom` [Pa] (default 0 and 10<sup>5</sup> Pa: every level) and 0 outside, with a sin<sup>2</sup> ramp of depth `taper_depth` [Pa] inside both ends (default 0: no ramp). For example, `p_top = 10000., p_bottom = 90000., taper_depth = 20000.` leaves the stratosphere and the boundary layer unperturbed.

**Reproducibility.** The random numbers are not taken from the compiler's random number generator but computed from a hash of `seed`, `ensemble_member`, the field, the level and the global indices of each spherical harmonic (or, for `'smooth_grid'`, of each lattice point of the noise). The perturbation is therefore the same for any number of MPI ranks or OpenMP threads, and the integer random numbers are the same with any compiler (the conversion to Gaussian numbers uses `log`, `cos` and `sqrt`, whose last bit can differ between math libraries). Widening the band only adds harmonics: those already in it keep their random numbers (the whole perturbation is then rescaled to `amplitude_t`). So for `'white'` and `'spectral'`, runs at T42 and T85 with the same levels and a band inside both truncations get the same perturbation.

**When it is applied.** `apply` says at which starts the perturbation is added:

* `'cold_start'` (default): only when there are no restart files. Later segments of the run, which start from its restart files, are not perturbed again.
* `'once'`: also at a restart whose restart files do not record a perturbation. This is the way to branch an ensemble from a spun-up state; the later segments of each member are not perturbed again.
* `'every_start'`: at every start, restarts included (the behaviour of v1's `add_noise`).

The record is a scalar `initial_noise_applied` in `spectral_dynamics.res.nc`, which is written once a perturbation has been added and kept by the later segments; the restart files of a run that was never perturbed do not have it. At a restart the same increments are added to the state in both `spectral_dynamics.res.nc` and `atmosphere.res.nc`.

**An ensemble.** To make members that differ only by their initial perturbation, give each member the same namelist except `ensemble_member` (0, 1, 2, ...); keep `seed` fixed, or change it for a new set of members. To branch the members from a spun-up control run, copy the control's restart files into the `INPUT/` directory of each member and use `apply = 'once'`:

```fortran
&initial_noise_nml
    noise_type      = 'spectral',
    spectral_slope  = 3.,
    amplitude_t     = 0.1,
    seed            = 2026,
    ensemble_member = 3,
    apply           = 'once' /
```

If the control itself started from a perturbed cold start, its restart files record that perturbation and `'once'` adds nothing (the log then says `no perturbation added at this restart`). In that case use `apply = 'every_start'` for the first segment of each member and `'once'` for the later ones.

The removed v1 keys `add_noise`, `add_noise_seed`, `noise_spectral_cutoff_minimum` and `noise_spectral_cutoff_maximum` of `spectral_dynamics_nml` stop the model with a message; see [Migration to v2](Migration_v2.md) for their replacements.
