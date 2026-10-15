# Diagnostics

MiMA writes only the diagnostics you ask for in the `diag_table` file in the run
directory. This page explains how to write a `diag_table` and lists every field the
model can provide, with the namelist settings it needs.

* [How to write a diag_table](#how-to-write-a-diag_table)
* [Checking a diag_table](#checking-a-diag_table)
* [Reading the field tables](#reading-the-field-tables)
* [Fields by module](#fields-by-module)

Ready-made tables are in the repository: `input/diag_table` (the default RRTM test
case), `input/examples/gray/diag_table` (gray radiation) and
`input/examples/held_suarez/diag_table` (dry Held-Suarez).

## How to write a diag_table

A `diag_table` is a plain text file of comma-separated values. Anything after a `#` is a
comment. It has three kinds of lines: the two global lines, which must be lines 1 and 2
of the file, then file lines and field lines in any order (define a file before the
fields that go into it).

```
"MiMA experiment"
0001 1 1 0 0 0
#  file name,   output frequency, units, format, time units, time axis name
"atmos_daily",    1, "days", 1, "days", "time",
"atmos_monthly", 30, "days", 1, "days", "time",
#  module,     field,    output name, file,           sampling, reduction, region, packing
 "dynamics", "ps",     "ps",        "atmos_daily",   "all",    .true.,    "none", 2,
 "dynamics", "ucomp",  "ucomp",     "atmos_monthly", "all",    .true.,    "none", 2,
 "moist",    "precip", "precip",    "atmos_daily",   "all",    .true.,    "none", 2,
 "moist",    "precip", "precip",    "atmos_daily",   "all",    "max",     "none", 2,
```

This writes daily means of surface pressure and precipitation plus the daily maximum
precipitation (`precip_max`) to `atmos_daily`, and 30-day means of the zonal wind to
`atmos_monthly`.

**Global lines.** Line 1 is a title (written to the file metadata). Line 2 is the base
date, `year month day hour minute second`: the reference time of the time axis. Use the
model start date (`current_date` in [`coupler_nml`](Parameters.md#coupler_nml)).

**File lines** define an output file:

| column | example | meaning |
|---|---|---|
| file name | `"atmos_daily"` | Each MPI process writes `atmos_daily.nc.NNNN`. Combine them with `mppnccombine` (see [Getting started](GettingStarted.md#output)). |
| output frequency | `1` | Write every N units (`> 0`), every time step (`0`), or once at the end of the run (`-1`). |
| units | `"days"` | Units of the output frequency: `"seconds"`, `"minutes"`, `"hours"`, `"days"`, `"months"` or `"years"`. |
| format | `1` | Always `1` (netCDF). |
| time units | `"days"` | Units of the time axis in the file. |
| time axis name | `"time"` | Name of the time axis. It must contain `time`. |

With the 30-day calendar of the test case, `30, "days"` and `1, "months"` are the same.

**Field lines** send one diagnostic field to a file:

| column | example | meaning |
|---|---|---|
| module | `"dynamics"` | Module name, from the tables below. It is case sensitive. |
| field | `"ucomp"` | Field name, from the tables below. It is not case sensitive. |
| output name | `"ucomp"` | Name of the variable in the output file. |
| file | `"atmos_daily"` | A file defined by a file line (without `.nc`). |
| sampling | `"all"` | Not used. Always `"all"`. |
| reduction | `.true.` | What to write at each output time. `.true.`: the time mean since the last output. `.false.`: the instantaneous value. `"max"` / `"min"`: the maximum / minimum since the last output (`_max` / `_min` is appended to the output name unless it already ends that way). `"rms"`, `"sum"`, `"pow2"` (mean of the square) and `"diurnal24"` (24 diurnal-cycle means) also work. Use `.false.` for static fields. |
| region | `"none"` | `"none"` writes the whole globe. A region is given as `"lon_min lon_max lat_min lat_max k_min k_max"` in degrees and level indices (`-1 -1` for all levels). |
| packing | `2` | Precision of the output: `1` = double (64-bit), `2` = single (32-bit float). `4` and `8` pack to 16-bit and 8-bit integers. |

A line may be at most 256 characters long. Two fields in the same file must have
different output names; to write the same field twice (for example its mean and its
maximum), give it a different output name or reduction.

**Listing the fields of a run.** The tables below are generated from the source code.
To see which fields are actually registered in a given configuration, add this to
`input.nml`:

```fortran
&diag_manager_nml
    do_diag_field_log = .true. /
```

The root process then writes `diag_field_log.out.0` in the run directory, with one
line per registered field, whether or not it is in `diag_table`:

```
Module|Field|Long Name|Units|Number of Axis|Time Axis|Missing Value|Min Value|Max Value|AXES LIST
dynamics|ucomp|zonal wind component|m/s|3|T||  -400.00000000000000|   400.00000000000000|lon,lat,pfull
```

**Fields that are not registered.** If a field line names a module/field that is not
registered in the run (a typo, or a field whose scheme or switch is off, such as an
RRTM field in a gray run), the model does not stop. When it opens the file, diag_manager
prints a warning, for example

```
WARNING from PE 0: diag_util_mod::opening_file: module/field_name (radiation/tdt_sw) NOT registered
```

and leaves the field out of the file. Errors in the table itself (a wrong number of
columns, an undefined file, a packing value outside 1-8, a duplicated output name) are
fatal.

## Checking a diag_table

`tools/mimadoc` checks a `diag_table` against the fields in the source code before you
run. It needs only Python 3.7 or later (no packages). For example, in a run directory:

```bash
python3 /path/to/MiMA/tools/mimadoc validate diag_table --nml input.nml
```

It reports field lines whose module/field is not registered anywhere (with a hint when
the field exists under another module), fields that are not available with the settings
in `input.nml` (`--nml` is optional; namelist variables not set there take their default
values from the source), and format errors: undefined files, unknown reduction methods,
packing values and duplicated output names. It exits with status 1 if it finds a
problem.

## Reading the field tables

* **field**: `<tracer>` stands for the name of any tracer in `field_table` (for example
  `sphum`). `{expr}` means the name is built at run time from `expr`.
* **dims**: `lon, lat` are the horizontal grid; `pfull` are the full model levels and
  `phalf` the half levels between them (one more than `pfull`); `lat, pfull` is a zonal
  mean. *static* fields have no time axis: they are written once, so use `.false.`
  for them.
* **available with**: the settings for which the field is registered. Namelist
  variables are followed by their namelist in parentheses: `do_bm` (moist_processes)
  is `do_bm` in `&moist_processes_nml`. Variables without a namelist are set inside
  the model (see the source). An empty cell means that the field is always registered
  when its module is active. Conditions that apply to every field of a module are given
  once, above the module's table.
* **source**: the Fortran module that registers the field (with `register_diag_field` or
  `register_static_field`), linked to its page in the [Fortran API reference](FortranAPI.md). Fields registered in several places, for example by
  different radiation schemes, have one row with all the places.
* A field whose units, long name or dims differ between the places where it is
  registered is marked **differs**. A field marked *never sent* is registered, but
  the model never passes it any data, so its output contains no model values.

## Fields by module

This list is generated from the Fortran sources by `tools/mimadoc`. After adding or
changing a diagnostic, regenerate it with `python3 tools/mimadoc generate`;
`python3 tools/mimadoc check` fails if it is out of date. Source files that CMake does
not compile are left out.

<!-- mimadoc:diagnostics -->

| module | fields | registered when | source |
|---|---:|---|---|
| [`cg_drag`](#module-cg_drag) | 5 | `do_cg_drag` (damping_driver) and `do_damping` (physics_driver) | [`cg_drag_mod`](https://eddy-stanford.github.io/MiMA/api/cg_drag_mod/) |
| [`climo`](#module-climo) | 1 |  | [`mima_interpolator_mod`](https://eddy-stanford.github.io/MiMA/api/mima_interpolator_mod/) |
| [`damping`](#module-damping) | 14 | `do_damping` (physics_driver) | [`damping_driver_mod`](https://eddy-stanford.github.io/MiMA/api/damping_driver_mod/) |
| [`dynamics`](#module-dynamics) | 24 |  | [`spectral_dynamics_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_dynamics_mod/) |
| [`dynamics_every`](#module-dynamics_every) | 47 |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| [`held_suarez`](#module-held_suarez) | 5 | `do_held_suarez` (physics_driver) | [`held_suarez_mod`](https://eddy-stanford.github.io/MiMA/api/held_suarez_mod/) |
| [`hinterp`](#module-hinterp) | 1 |  | [`mima_interpolator_mod`](https://eddy-stanford.github.io/MiMA/api/mima_interpolator_mod/) |
| [`local_heating`](#module-local_heating) | 1 | `do_3d_heating` and `do_local_heating` (physics_driver) | [`local_heating_mod`](https://eddy-stanford.github.io/MiMA/api/local_heating_mod/) |
| [`mca`](#module-mca) | 8 | `do_mca` (moist_processes) | [`moist_conv_mod`](https://eddy-stanford.github.io/MiMA/api/moist_conv_mod/) |
| [`moist`](#module-moist) | 30 |  | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |
| [`radiation`](#module-radiation) | 19 | `radiation_scheme` (radiation): gray, rrtm | [`gray_radiation_mod`](https://eddy-stanford.github.io/MiMA/api/gray_radiation_mod/), [`rrtm_radiation`](https://eddy-stanford.github.io/MiMA/api/rrtm_radiation/) |
| [`simple_surface`](#module-simple_surface) | 31 |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| [`tracers`](#module-tracers) | 6 |  | [`atmos_carbon_aerosol_mod`](https://eddy-stanford.github.io/MiMA/api/atmos_carbon_aerosol_mod/), [`atmos_sulfur_hex_mod`](https://eddy-stanford.github.io/MiMA/api/atmos_sulfur_hex_mod/), [`atmos_tracer_utilities_mod`](https://eddy-stanford.github.io/MiMA/api/atmos_tracer_utilities_mod/) |
| [`vert_diff`](#module-vert_diff) | 10 |  | [`vert_diff_driver_mod`](https://eddy-stanford.github.io/MiMA/api/vert_diff_driver_mod/) |
| [`vert_turb`](#module-vert_turb) | 8 |  | [`vert_turb_driver_mod`](https://eddy-stanford.github.io/MiMA/api/vert_turb_driver_mod/) |

<a id="module-cg_drag"></a>

### `cg_drag`

Registered only when `do_cg_drag` (damping_driver) and `do_damping` (physics_driver).

| field | dims | units | long_name | available with | source |
|---|---|---|---|---|---|
| `bf_cgwd` | lon, lat, pfull | 1/s | buoyancy frequency from cg_drag |  | [`cg_drag_mod`](https://eddy-stanford.github.io/MiMA/api/cg_drag_mod/) |
| `gwfu_cgwd` | lon, lat, pfull | m/s2 | gravity wave forcing on mean zonal flow |  | [`cg_drag_mod`](https://eddy-stanford.github.io/MiMA/api/cg_drag_mod/) |
| `gwfv_cgwd` | lon, lat, pfull | m/s2 | gravity wave forcing on mean meridional flow |  | [`cg_drag_mod`](https://eddy-stanford.github.io/MiMA/api/cg_drag_mod/) |
| `kedx_cgwd` | lon, lat, pfull | m2/s | effective eddy viscosity from cg_drag (zonal) |  | [`cg_drag_mod`](https://eddy-stanford.github.io/MiMA/api/cg_drag_mod/) |
| `kedy_cgwd` | lon, lat, pfull | m2/s | effective eddy viscosity from cg_drag (meridional) |  | [`cg_drag_mod`](https://eddy-stanford.github.io/MiMA/api/cg_drag_mod/) |

<a id="module-climo"></a>

### `climo`

| field | dims | units | long_name | available with | source |
|---|---|---|---|---|---|
| `{clim_type%field_name(i)}` | lon, lat | kg/m2 | column integral of {clim_type%field_name(i)} (climatology grid) |  | [`mima_interpolator_mod`](https://eddy-stanford.github.io/MiMA/api/mima_interpolator_mod/) |

<a id="module-damping"></a>

### `damping`

Registered only when `do_damping` (physics_driver).

| field | dims | units | long_name | available with | source |
|---|---|---|---|---|---|
| `diss_heat_gwd` | lon, lat | W/m2 | Integrated dissipative heating from gravity wave drag | `do_mg_drag` (damping_driver) | [`damping_driver_mod`](https://eddy-stanford.github.io/MiMA/api/damping_driver_mod/) |
| `diss_heat_rdamp` | lon, lat | W/m2 | Integrated dissipative heating from Rayleigh damping | `do_rayleigh` (damping_driver) | [`damping_driver_mod`](https://eddy-stanford.github.io/MiMA/api/damping_driver_mod/) |
| `sgsmtn` | lon, lat; static | m | sub-grid scale topography for gravity wave drag | `do_mg_drag` (damping_driver) | [`damping_driver_mod`](https://eddy-stanford.github.io/MiMA/api/damping_driver_mod/) |
| `taubx` | lon, lat | N/m2 | x base flux for grav wave drag | `do_mg_drag` (damping_driver) | [`damping_driver_mod`](https://eddy-stanford.github.io/MiMA/api/damping_driver_mod/) |
| `tauby` | lon, lat | N/m2 | y base flux for grav wave drag | `do_mg_drag` (damping_driver) | [`damping_driver_mod`](https://eddy-stanford.github.io/MiMA/api/damping_driver_mod/) |
| `taus` | lon, lat, pfull | N/m2 | saturation flux for gravity wave drag | `do_mg_drag` (damping_driver) | [`damping_driver_mod`](https://eddy-stanford.github.io/MiMA/api/damping_driver_mod/) |
| `tdt_diss_gwd` | lon, lat, pfull | K/s | Dissipative heating from gravity wave drag | `do_mg_drag` (damping_driver) | [`damping_driver_mod`](https://eddy-stanford.github.io/MiMA/api/damping_driver_mod/) |
| `tdt_diss_rdamp` | lon, lat, pfull | K/s | Dissipative heating from Rayleigh damping | `do_rayleigh` (damping_driver) | [`damping_driver_mod`](https://eddy-stanford.github.io/MiMA/api/damping_driver_mod/) |
| `udt_cgwd` | lon, lat, pfull | m/s2 | u wind tendency for cg gravity wave drag | `do_cg_drag` (damping_driver) | [`damping_driver_mod`](https://eddy-stanford.github.io/MiMA/api/damping_driver_mod/) |
| `udt_cnstd` | lon, lat, pfull | m/s2 | u wind tendency for constant drag | `do_const_drag` (damping_driver) | [`damping_driver_mod`](https://eddy-stanford.github.io/MiMA/api/damping_driver_mod/) |
| `udt_gwd` | lon, lat, pfull | m/s2 | u wind tendency for gravity wave drag | `do_mg_drag` (damping_driver) | [`damping_driver_mod`](https://eddy-stanford.github.io/MiMA/api/damping_driver_mod/) |
| `udt_rdamp` | lon, lat, pfull | m/s2 | u wind tendency for Rayleigh damping | `do_rayleigh` (damping_driver) | [`damping_driver_mod`](https://eddy-stanford.github.io/MiMA/api/damping_driver_mod/) |
| `vdt_gwd` | lon, lat, pfull | m/s2 | v wind tendency for gravity wave drag | `do_mg_drag` (damping_driver) | [`damping_driver_mod`](https://eddy-stanford.github.io/MiMA/api/damping_driver_mod/) |
| `vdt_rdamp` | lon, lat, pfull | m/s2 | v wind tendency for Rayleigh damping | `do_rayleigh` (damping_driver) | [`damping_driver_mod`](https://eddy-stanford.github.io/MiMA/api/damping_driver_mod/) |

<a id="module-dynamics"></a>

### `dynamics`

| field | dims | units | long_name | available with | source |
|---|---|---|---|---|---|
| `<tracer>` | lon, lat, pfull | &lt;tracer_units&gt; | &lt;tracer_longname&gt; |  | [`spectral_dynamics_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_dynamics_mod/) |
| `bk` | phalf; static | 1 | vertical coordinate sigma values |  | [`spectral_dynamics_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_dynamics_mod/) |
| `div` | lon, lat, pfull | 1/s | divergence |  | [`spectral_dynamics_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_dynamics_mod/) |
| `div217` | lon, lat | 1/s | divergence at model level 9 |  | [`spectral_dynamics_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_dynamics_mod/) |
| `height` | lon, lat, pfull | m | geopotential height at full model levels |  | [`spectral_dynamics_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_dynamics_mod/) |
| `height_half` | lon, lat, phalf | m | geopotential height at half model levels |  | [`spectral_dynamics_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_dynamics_mod/) |
| `omega` | lon, lat, pfull | Pa/s | dp/dt vertical velocity |  | [`spectral_dynamics_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_dynamics_mod/) |
| `omega_sq` | lon, lat, pfull | Pa2/s2 | omega squared |  | [`spectral_dynamics_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_dynamics_mod/) |
| `omega_temp` | lon, lat, pfull | Pa K/s | dp/dt * temperature |  | [`spectral_dynamics_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_dynamics_mod/) |
| `pk` | phalf; static | Pa | vertical coordinate pressure values |  | [`spectral_dynamics_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_dynamics_mod/) |
| `pres_full` | lon, lat, pfull | Pa | pressure at full model levels |  | [`spectral_dynamics_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_dynamics_mod/) |
| `pres_half` | lon, lat, phalf | Pa | pressure at half model levels |  | [`spectral_dynamics_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_dynamics_mod/) |
| `ps` | lon, lat | Pa | surface pressure |  | [`spectral_dynamics_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_dynamics_mod/) |
| `slp` | lon, lat | Pa | sea level pressure |  | [`spectral_dynamics_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_dynamics_mod/) |
| `temp` | lon, lat, pfull | K | temperature |  | [`spectral_dynamics_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_dynamics_mod/) |
| `temp_sq` | lon, lat, pfull | K2 | temperature squared |  | [`spectral_dynamics_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_dynamics_mod/) |
| `ucomp` | lon, lat, pfull | m/s | zonal wind component |  | [`spectral_dynamics_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_dynamics_mod/) |
| `ucomp_sq` | lon, lat, pfull | m2/s2 | zonal wind squared |  | [`spectral_dynamics_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_dynamics_mod/) |
| `ucomp_vcomp` | lon, lat, pfull | m2/s2 | zonal wind * meridional wind |  | [`spectral_dynamics_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_dynamics_mod/) |
| `vcomp` | lon, lat, pfull | m/s | meridional wind component |  | [`spectral_dynamics_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_dynamics_mod/) |
| `vcomp_sq` | lon, lat, pfull | m2/s2 | meridional wind squared |  | [`spectral_dynamics_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_dynamics_mod/) |
| `vor` | lon, lat, pfull | 1/s | vorticity |  | [`spectral_dynamics_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_dynamics_mod/) |
| `wspd` | lon, lat, pfull | m/s | wind speed |  | [`spectral_dynamics_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_dynamics_mod/) |
| `zsurf` | lon, lat; static | m | geopotential height at the surface |  | [`spectral_dynamics_mod`](https://eddy-stanford.github.io/MiMA/api/spectral_dynamics_mod/) |

<a id="module-dynamics_every"></a>

### `dynamics_every`

| field | dims | units | long_name | available with | source |
|---|---|---|---|---|---|
| `2dt_<tracer>` | lon, lat, pfull; static | &lt;tracer_units&gt; | Amplitude of 2*dt wave in &lt;tracer_longname&gt; |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `2dt_ps` | lon, lat; static | Pa | Amplitude of 2*dt wave in surface pressure |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `2dt_t` | lon, lat, pfull; static | K | Amplitude of 2*dt wave in temperature |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `2dt_u` | lon, lat, pfull; static | m/s | Amplitude of 2*dt wave in zonal wind |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `2dt_v` | lon, lat, pfull; static | m/s | Amplitude of 2*dt wave in meridional wind |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `<tracer>_every` | lon, lat, pfull | &lt;tracer_units&gt; | &lt;tracer_longname&gt; |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `<tracer>_hadv` | lon, lat, pfull | &lt;tracer_units&gt; kg/m2/s | Global mean column integrated &lt;tracer_longname&gt; tendency due to horizontal advection |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `<tracer>_vadv` | lon, lat, pfull | &lt;tracer_units&gt; kg/m2/s | Global mean column integrated &lt;tracer_longname&gt; tendency due to vertical advection |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `dry_stat_en` | lon, lat, pfull | J/kg | dry static energy |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `entrop_dampt` | lon, lat, pfull | 1/s | Entropy prod from horiz diff of temp |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `entrop_dampuv` | lon, lat, pfull | 1/s | Entropy prod from horiz diff of vel |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `entrop_tempcor` | lon, lat, pfull | 1/s | Entropy tend due to hor diff temp corr |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `kegen` | lon, lat, pfull | K/s | kappa omega T/p from the dynamical core (Tv if use_virtual_temperature), weighted by ps/p00 |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `kegenq` | lon, lat, pfull | K/s kg/kg | kappa omega T q/p from the dynamical core (Tv if use_virtual_temperature), weighted by ps/p00 |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `kegenqtinv` | lon, lat, pfull | kg/kg/s | kappa omega q/p from the dynamical core, weighted by ps/p00 |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `moist_stat_en` | lon, lat, pfull | J/kg | moist static energy |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `omegap` | lon, lat, pfull | Pa/s | dp/dt vertical velocity, weighted by ps/p00 |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `ps_every` | lon, lat | Pa | surface pressure |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `psi_dwc` | lat, pfull | m/s Pa | residual mean streamfunction from downward control (EP flux divergence) |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `psi_star` | lat, pfull | m/s Pa | residual mean streamfunction (Eulerian mean minus eddy heat flux term) |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `qdt_watercor` | lon, lat, pfull | kg/kg/s | Humidity tend due to dynamics water correction |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `t_every` | lon, lat, pfull | K | temperature |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `tdt_damp` | lon, lat, pfull | K/s | Temperature tend from horiz diff |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `tdt_dampuv` | lon, lat, pfull | K/s | Temp tend (not applied) due to diff of winds |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `tdt_tempcor` | lon, lat, pfull | K/s | Temp tend due to dynamics temp correction |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `u_every` | lon, lat, pfull | m/s | zonal wind component |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `udt_damp` | lon, lat, pfull | m/s2 | Zonal wind tend from horiz diff |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `upvp` | lon, lat, pfull | m2/s2 | meridional eddy momentum flux |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `upwp` | lon, lat, pfull | m/s Pa/s | vertical eddy momentum flux |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `uv` | lon, lat, pfull | m2/s2 | zonal wind times meridional wind, weighted by ps/p00 |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `v_every` | lon, lat, pfull | m/s | meridional wind component |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `vdse` | lon, lat, pfull | m/s J/kg | meridional wind times dry static energy, weighted by ps/p00 |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `vdt_damp` | lon, lat, pfull | m/s2 | Merid wind tend from horiz diff |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `vp` | lon, lat, pfull | m/s | meridional wind, weighted by ps/p00 |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `vpTp` | lon, lat, pfull | m/s K | meridional eddy heat flux |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `vq` | lon, lat, pfull | m/s kg/kg | meridional wind times specific humidity, weighted by ps/p00 |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `vqint` | lon, lat | kg/m/s | vertically integrated meridional moisture flux |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `wdse` | lon, lat, pfull | Pa/s J/kg | omega DSE |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `wdsep` | lon, lat, pfull | Pa/s J/kg | omega DSE, weighted by ps/p00 |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `wq` | lon, lat, pfull | Pa/s kg/kg | omega q |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `wqp` | lon, lat, pfull | Pa/s kg/kg | omega q, weighted by ps/p00 |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `wsubt` | lon, lat, pfull | K/s | kappa omega T/p, weighted by ps/p00 |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `wsubtv` | lon, lat, pfull | K/s | kappa omega Tv/p, weighted by ps/p00 |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `wu` | lon, lat, pfull | Pa/s m/s | omega u |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `wup` | lon, lat, pfull | Pa/s m/s | omega u, weighted by ps/p00 |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `wv` | lon, lat, pfull | Pa/s m/s | omega v |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |
| `wvp` | lon, lat, pfull | Pa/s m/s | omega v, weighted by ps/p00 |  | [`every_step_diagnostics_mod`](https://eddy-stanford.github.io/MiMA/api/every_step_diagnostics_mod/) |

<a id="module-held_suarez"></a>

### `held_suarez`

Registered only when `do_held_suarez` (physics_driver).

| field | dims | units | long_name | available with | source |
|---|---|---|---|---|---|
| `diss_heat_hs` | lon, lat, pfull | K/s | Heating from Held-Suarez frictional dissipation |  | [`held_suarez_mod`](https://eddy-stanford.github.io/MiMA/api/held_suarez_mod/) |
| `tdt_hs` | lon, lat, pfull | K/s | Temperature tendency from Held-Suarez relaxation |  | [`held_suarez_mod`](https://eddy-stanford.github.io/MiMA/api/held_suarez_mod/) |
| `teq` | lon, lat, pfull | K | Held-Suarez equilibrium temperature |  | [`held_suarez_mod`](https://eddy-stanford.github.io/MiMA/api/held_suarez_mod/) |
| `udt_hs` | lon, lat, pfull | m/s2 | Zonal wind tendency from Held-Suarez friction |  | [`held_suarez_mod`](https://eddy-stanford.github.io/MiMA/api/held_suarez_mod/) |
| `vdt_hs` | lon, lat, pfull | m/s2 | Meridional wind tendency from Held-Suarez friction |  | [`held_suarez_mod`](https://eddy-stanford.github.io/MiMA/api/held_suarez_mod/) |

<a id="module-hinterp"></a>

### `hinterp`

| field | dims | units | long_name | available with | source |
|---|---|---|---|---|---|
| `{clim_type%field_name(i)}` | lon, lat | kg/m2 | column integral of {clim_type%field_name(i)} (interpolated to the model grid) |  | [`mima_interpolator_mod`](https://eddy-stanford.github.io/MiMA/api/mima_interpolator_mod/) |

<a id="module-local_heating"></a>

### `local_heating`

Registered only when `do_3d_heating` and `do_local_heating` (physics_driver).

| field | dims | units | long_name | available with | source |
|---|---|---|---|---|---|
| `tdt_lheat` | lon, lat, pfull | K/s | Temperature tendency due to local heating |  | [`local_heating_mod`](https://eddy-stanford.github.io/MiMA/api/local_heating_mod/) |

<a id="module-mca"></a>

### `mca`

Registered only when `do_mca` (moist_processes).

| field | dims | units | long_name | available with | source |
|---|---|---|---|---|---|
| `<tracer>dt_MCA` | lon, lat, pfull | &lt;tracer_units&gt;/s | &lt;tracer&gt; tendency from MCA | `tracers_in_mca(tr)` | [`moist_conv_mod`](https://eddy-stanford.github.io/MiMA/api/moist_conv_mod/) |
| `<tracer>dt_MCA_col` | lon, lat | &lt;tracer_units&gt; kg/m2/s | &lt;tracer&gt; path tendency from MCA | `tracers_in_mca(tr)` | [`moist_conv_mod`](https://eddy-stanford.github.io/MiMA/api/moist_conv_mod/) |
| `prec_conv` | lon, lat | kg/m2/s | Precipitation rate from moist conv adj |  | [`moist_conv_mod`](https://eddy-stanford.github.io/MiMA/api/moist_conv_mod/) |
| `q_conv_col` | lon, lat | kg/m2/s | Water vapor path tendency from moist conv adj |  | [`moist_conv_mod`](https://eddy-stanford.github.io/MiMA/api/moist_conv_mod/) |
| `qdt_conv` | lon, lat, pfull | kg/kg/s | Spec humidity tendency from moist conv adj |  | [`moist_conv_mod`](https://eddy-stanford.github.io/MiMA/api/moist_conv_mod/) |
| `snow_conv` | lon, lat | kg/m2/s | Frozen precip rate from moist conv adj |  | [`moist_conv_mod`](https://eddy-stanford.github.io/MiMA/api/moist_conv_mod/) |
| `t_conv_col` | lon, lat | W/m2 | Column static energy tendency from moist conv adj |  | [`moist_conv_mod`](https://eddy-stanford.github.io/MiMA/api/moist_conv_mod/) |
| `tdt_conv` | lon, lat, pfull | K/s | Temperature tendency from moist conv adj |  | [`moist_conv_mod`](https://eddy-stanford.github.io/MiMA/api/moist_conv_mod/) |

<a id="module-moist"></a>

### `moist`

| field | dims | units | long_name | available with | source |
|---|---|---|---|---|---|
| `<tracer>` | lon, lat, pfull | &lt;tracer_units&gt; | &lt;tracer&gt; |  | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |
| `<tracer>_col` | lon, lat | &lt;tracer_units&gt; kg/m2 | column integrated &lt;tracer&gt; |  | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |
| `<tracer>dt_conv` | lon, lat, pfull | &lt;tracer_units&gt;/s | &lt;tracer&gt; total tendency from moist convection | `tracers_in_mca(n)` | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |
| `<tracer>dt_conv_col` | lon, lat | &lt;tracer_units&gt; kg/m2/s | &lt;tracer&gt; total path tendency from moist convection | `tracers_in_mca(n)` | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |
| `bmflag` | lon, lat | 1 | Betts-Miller flag | `do_bm` (moist_processes) | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |
| `cape` | lon, lat | J/kg | Convectively available potential energy |  | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |
| `cin` | lon, lat | J/kg | Convective inhibition |  | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |
| `entrop_ls` | lon, lat, pfull | 1/s | Entropy tendency from large-scale cond | `do_lsc` (moist_processes) | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |
| `gust_conv` | lon, lat | m/s | Gustiness from deep convection |  | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |
| `invtaubmq` | lon, lat | 1/s | Inverse humidity relaxation time | `do_bm` (moist_processes) | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |
| `invtaubmt` | lon, lat | 1/s | Inverse temperature relaxation time | `do_bm` (moist_processes) | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |
| `klzbs` | lon, lat | 1 | Betts-Miller level of zero buoyancy (model level index) | `do_bm` (moist_processes) | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |
| `prec_conv` | lon, lat | kg/m2/s | Precipitation rate |  | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |
| `prec_ls` | lon, lat | kg/m2/s | Precipitation rate from large-scale cond | `do_lsc` (moist_processes) | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |
| `precip` | lon, lat | kg/m2/s | Total precipitation rate |  | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |
| `q_conv_col` | lon, lat | kg/m2/s | Water vapor path tendency |  | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |
| `q_ls_col` | lon, lat | kg/m2/s | Water vapor path tendency from large-scale cond | `do_lsc` (moist_processes) | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |
| `qdt_conv` | lon, lat, pfull | kg/kg/s | Spec humidity tendency |  | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |
| `qdt_ls` | lon, lat, pfull | kg/kg/s | Spec humidity tendency from large-scale cond | `do_lsc` (moist_processes) | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |
| `qref` | lon, lat, pfull | kg/kg | Adjustment reference specific humidity profile | `do_bm` (moist_processes) | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |
| `rh` | lon, lat, pfull | percent | relative humidity |  | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |
| `rhsurf` | lon, lat | percent | Relative humidity at the lowest model level |  | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |
| `snow_conv` | lon, lat | kg/m2/s | Frozen precip rate |  | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |
| `snow_ls` | lon, lat | kg/m2/s | Frozen precip rate from large-scale cond | `do_lsc` (moist_processes) | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |
| `t_conv_col` | lon, lat | W/m2 | Column static energy tendency |  | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |
| `t_ls_col` | lon, lat | W/m2 | Column static energy tendency from large-scale cond | `do_lsc` (moist_processes) | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |
| `tdt_conv` | lon, lat, pfull | K/s | Temperature tendency |  | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |
| `tdt_ls` | lon, lat, pfull | K/s | Temperature tendency from large-scale cond | `do_lsc` (moist_processes) | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |
| `tref` | lon, lat, pfull | K | Adjustment reference temperature profile | `do_bm` (moist_processes) | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |
| `WVP` | lon, lat | kg/m2 | Column integrated water vapor |  | [`moist_processes_mod`](https://eddy-stanford.github.io/MiMA/api/moist_processes_mod/) |

<a id="module-radiation"></a>

### `radiation`

Which fields exist depends on `radiation_scheme` in `&radiation_nml`: **available with** lists the values for which the field is registered.

| field | dims | units | long_name | available with | source |
|---|---|---|---|---|---|
| `albedo_rad` | lon, lat | 1 | Surface albedo seen by the radiation | gray, rrtm | [`gray_radiation_mod`](https://eddy-stanford.github.io/MiMA/api/gray_radiation_mod/)<br>[`rrtm_radiation`](https://eddy-stanford.github.io/MiMA/api/rrtm_radiation/) |
| `coszen` | lon, lat | 1 | cosine of zenith angle | rrtm | [`rrtm_radiation`](https://eddy-stanford.github.io/MiMA/api/rrtm_radiation/) |
| `entrop_rad` | lon, lat, pfull | 1/s | Entropy production by radiation | gray | [`gray_radiation_mod`](https://eddy-stanford.github.io/MiMA/api/gray_radiation_mod/) |
| `lwdn_sfc` | lon, lat | W/m2 | LW flux down at surface | gray, rrtm | [`gray_radiation_mod`](https://eddy-stanford.github.io/MiMA/api/gray_radiation_mod/)<br>[`rrtm_radiation`](https://eddy-stanford.github.io/MiMA/api/rrtm_radiation/) |
| `lwnet_half` | lon, lat, phalf | W/m2 | Net LW flux on half levels (positive up) | gray | [`gray_radiation_mod`](https://eddy-stanford.github.io/MiMA/api/gray_radiation_mod/) |
| `lwup_sfc` | lon, lat | W/m2 | LW flux up at surface | gray, rrtm | [`gray_radiation_mod`](https://eddy-stanford.github.io/MiMA/api/gray_radiation_mod/)<br>[`rrtm_radiation`](https://eddy-stanford.github.io/MiMA/api/rrtm_radiation/) |
| `netrad_half` | lon, lat, phalf | W/m2 | Net radiative flux on half levels (positive up) | gray | [`gray_radiation_mod`](https://eddy-stanford.github.io/MiMA/api/gray_radiation_mod/) |
| `olr` | lon, lat | W/m2 | Outgoing longwave radiation at TOA | gray, rrtm | [`gray_radiation_mod`](https://eddy-stanford.github.io/MiMA/api/gray_radiation_mod/)<br>[`rrtm_radiation`](https://eddy-stanford.github.io/MiMA/api/rrtm_radiation/) |
| `ozone` | lon, lat, pfull | kg/kg | Ozone mass mixing ratio | rrtm | [`rrtm_radiation`](https://eddy-stanford.github.io/MiMA/api/rrtm_radiation/) |
| `swdn_toa` | lon, lat | W/m2 | SW flux down at TOA | gray, rrtm | [`gray_radiation_mod`](https://eddy-stanford.github.io/MiMA/api/gray_radiation_mod/)<br>[`rrtm_radiation`](https://eddy-stanford.github.io/MiMA/api/rrtm_radiation/) |
| `swnet_half` | lon, lat, phalf | W/m2 | Net SW flux on half levels (positive up) | gray | [`gray_radiation_mod`](https://eddy-stanford.github.io/MiMA/api/gray_radiation_mod/) |
| `swnet_sfc` | lon, lat | W/m2 | Net SW flux at surface (positive down) | gray, rrtm | [`gray_radiation_mod`](https://eddy-stanford.github.io/MiMA/api/gray_radiation_mod/)<br>[`rrtm_radiation`](https://eddy-stanford.github.io/MiMA/api/rrtm_radiation/) |
| `swnet_toa` | lon, lat | W/m2 | Net SW flux at TOA (positive down) | gray, rrtm | [`gray_radiation_mod`](https://eddy-stanford.github.io/MiMA/api/gray_radiation_mod/)<br>[`rrtm_radiation`](https://eddy-stanford.github.io/MiMA/api/rrtm_radiation/) |
| `tau_lw` | lon, lat, phalf | 1 | LW optical depth on half levels | gray | [`gray_radiation_mod`](https://eddy-stanford.github.io/MiMA/api/gray_radiation_mod/) |
| `tau_sw` | lon, lat, phalf | 1 | SW optical depth on half levels | gray | [`gray_radiation_mod`](https://eddy-stanford.github.io/MiMA/api/gray_radiation_mod/) |
| `tdt_lw` | lon, lat, pfull | K/s | Temperature tendency due to LW radiation | rrtm | [`rrtm_radiation`](https://eddy-stanford.github.io/MiMA/api/rrtm_radiation/) |
| `tdt_rad` | lon, lat, pfull | K/s | Temperature tendency due to radiation | gray, rrtm | [`gray_radiation_mod`](https://eddy-stanford.github.io/MiMA/api/gray_radiation_mod/)<br>[`rrtm_radiation`](https://eddy-stanford.github.io/MiMA/api/rrtm_radiation/) |
| `tdt_sw` | lon, lat, pfull | K/s | Temperature tendency due to SW radiation | rrtm | [`rrtm_radiation`](https://eddy-stanford.github.io/MiMA/api/rrtm_radiation/) |
| `thalf` | lon, lat, phalf | K | Temperature on half levels | rrtm | [`rrtm_radiation`](https://eddy-stanford.github.io/MiMA/api/rrtm_radiation/) |

<a id="module-simple_surface"></a>

### `simple_surface`

| field | dims | units | long_name | available with | source |
|---|---|---|---|---|---|
| `albedo` | lon, lat | 1 | surface albedo |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| `b_star` | lon, lat | m/s2 | buoyancy scale |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| `del_h` | lon, lat | 1 | Monin-Obukhov profile factor (T({z_ref_heat})-T_surf)/(T_atm-T_surf) |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| `del_m` | lon, lat | 1 | Monin-Obukhov profile factor u({z_ref_mom})/u_atm |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| `del_q` | lon, lat | 1 | Monin-Obukhov profile factor (q({z_ref_heat})-q_surf)/(q_atm-q_surf) |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| `drag_heat` | lon, lat | 1 | drag coeff for heat |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| `drag_moist` | lon, lat | 1 | drag coeff for moisture |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| `drag_mom` | lon, lat | 1 | drag coeff for momentum |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| `entrop_evap` | lon, lat | kg/m2/s/K | entropy source from evap |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| `entrop_lwflx` | lon, lat | W/m2/K | entropy source from LW flux |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| `entrop_shflx` | lon, lat | W/m2/K | entropy source from SH flux |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| `evap` | lon, lat | kg/m2/s | evaporation rate |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| `heat_capacity` | lon, lat | J/m2/K | mixed layer heat capacity |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| `lwflx` | lon, lat | W/m2 | net (down-up) longwave flux |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| `oflx` | lon, lat | W/m2 | prescribed ocean heat divergence |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| `rh_ref` | lon, lat | percent | relative humidity at {z_ref_heat} (100 q/q_sat) |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| `rough_heat` | lon, lat | m | surface roughness for heat |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| `rough_moist` | lon, lat | m | surface roughness for moisture |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| `rough_mom` | lon, lat | m | surface roughness for momentum |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| `shflx` | lon, lat | W/m2 | sensible heat flux |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| `t_atm` | lon, lat | K | temperature at btm level |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| `t_ref` | lon, lat | K | air temperature at {z_ref_heat} |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| `t_surf` | lon, lat | K | surface temperature |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| `tau_x` | lon, lat | N/m2 | zonal surface stress on the atmosphere (positive eastward) |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| `tau_y` | lon, lat | N/m2 | meridional surface stress on the atmosphere (positive northward) |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| `u_atm` | lon, lat | m/s | u wind component at btm level |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| `u_ref` | lon, lat | m/s | zonal wind at {z_ref_mom} |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| `u_star` | lon, lat | m/s | friction velocity |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| `v_atm` | lon, lat | m/s | v wind component at btm level |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| `v_ref` | lon, lat | m/s | meridional wind at {z_ref_mom} |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |
| `wind` | lon, lat | m/s | wind speed for flux calculations |  | [`simple_surface_mod`](https://eddy-stanford.github.io/MiMA/api/simple_surface_mod/) |

<a id="module-tracers"></a>

### `tracers`

| field | dims | units | long_name | available with | source |
|---|---|---|---|---|---|
| `bcemiss` | lon, lat; static | g/m2/s | black carbon emission | `nbcphobic > 0` | [`atmos_carbon_aerosol_mod`](https://eddy-stanford.github.io/MiMA/api/atmos_carbon_aerosol_mod/) |
| `ocemiss` | lon, lat; static | g/m2/s | organic carbon emission | `nbcphobic > 0` | [`atmos_carbon_aerosol_mod`](https://eddy-stanford.github.io/MiMA/api/atmos_carbon_aerosol_mod/) |
| `sf6emiss` | lon, lat; static | g/m2/s | SF6 emission | `nsf6 > 0` | [`atmos_sulfur_hex_mod`](https://eddy-stanford.github.io/MiMA/api/atmos_sulfur_hex_mod/) |
| `{tracer_ddep_names(n)}` | lon, lat | {tracer_units(n)} kg/m2/s | {tracer_ddep_longnames(n)} |  | [`atmos_tracer_utilities_mod`](https://eddy-stanford.github.io/MiMA/api/atmos_tracer_utilities_mod/) |
| `{tracer_wdep_names(n)}_cv` | lon, lat | {tracer_units(n)} kg/m2/s | {tracer_wdep_longnames(n)} in convective scheme |  | [`atmos_tracer_utilities_mod`](https://eddy-stanford.github.io/MiMA/api/atmos_tracer_utilities_mod/) |
| `{tracer_wdep_names(n)}_ls` | lon, lat | {tracer_units(n)} kg/m2/s | {tracer_wdep_longnames(n)} in large scale |  | [`atmos_tracer_utilities_mod`](https://eddy-stanford.github.io/MiMA/api/atmos_tracer_utilities_mod/) |

<a id="module-vert_diff"></a>

### `vert_diff`

| field | dims | units | long_name | available with | source |
|---|---|---|---|---|---|
| `diss_heat_vdif` | lon, lat | W/m2 | Integrated dissipative heating from vert diff |  | [`vert_diff_driver_mod`](https://eddy-stanford.github.io/MiMA/api/vert_diff_driver_mod/) |
| `entrop_vdif_kediss` | lon, lat, pfull | 1/s | Entropy tendency from vert diff kinetic energy dissipation |  | [`vert_diff_driver_mod`](https://eddy-stanford.github.io/MiMA/api/vert_diff_driver_mod/) |
| `entrop_vdif_sens` | lon, lat, pfull | 1/s | Entropy tendency from vert diff sensible heating |  | [`vert_diff_driver_mod`](https://eddy-stanford.github.io/MiMA/api/vert_diff_driver_mod/) |
| `evap_vdif` | lon, lat | kg/m2/s | Integrated moisture flux from vert diff |  | [`vert_diff_driver_mod`](https://eddy-stanford.github.io/MiMA/api/vert_diff_driver_mod/) |
| `qdt_vdif` | lon, lat, pfull | kg/kg/s | Spec humidity tendency from vert diff |  | [`vert_diff_driver_mod`](https://eddy-stanford.github.io/MiMA/api/vert_diff_driver_mod/) |
| `sens_vdif` | lon, lat | W/m2 | Integrated heat flux from vert diff |  | [`vert_diff_driver_mod`](https://eddy-stanford.github.io/MiMA/api/vert_diff_driver_mod/) |
| `tdt_diss_vdif` | lon, lat, pfull | K/s | Dissipative heating from vert diff |  | [`vert_diff_driver_mod`](https://eddy-stanford.github.io/MiMA/api/vert_diff_driver_mod/) |
| `tdt_vdif` | lon, lat, pfull | K/s | Temperature tendency from vert diff |  | [`vert_diff_driver_mod`](https://eddy-stanford.github.io/MiMA/api/vert_diff_driver_mod/) |
| `udt_vdif` | lon, lat, pfull | m/s2 | Zonal wind tendency from vert diff |  | [`vert_diff_driver_mod`](https://eddy-stanford.github.io/MiMA/api/vert_diff_driver_mod/) |
| `vdt_vdif` | lon, lat, pfull | m/s2 | Meridional wind tendency from vert diff |  | [`vert_diff_driver_mod`](https://eddy-stanford.github.io/MiMA/api/vert_diff_driver_mod/) |

<a id="module-vert_turb"></a>

### `vert_turb`

| field | dims | units | long_name | available with | source |
|---|---|---|---|---|---|
| `diff_m` | lon, lat, phalf | m2/s | vert diff coeff for momentum |  | [`vert_turb_driver_mod`](https://eddy-stanford.github.io/MiMA/api/vert_turb_driver_mod/) |
| `diff_t` | lon, lat, phalf | m2/s | vert diff coeff for temp |  | [`vert_turb_driver_mod`](https://eddy-stanford.github.io/MiMA/api/vert_turb_driver_mod/) |
| `gust` | lon, lat | m/s | wind gustiness in surface layer |  | [`vert_turb_driver_mod`](https://eddy-stanford.github.io/MiMA/api/vert_turb_driver_mod/) |
| `uwnd` | lon, lat, pfull | m/s | zonal wind on mass grid |  | [`vert_turb_driver_mod`](https://eddy-stanford.github.io/MiMA/api/vert_turb_driver_mod/) |
| `vwnd` | lon, lat, pfull | m/s | meridional wind on mass grid |  | [`vert_turb_driver_mod`](https://eddy-stanford.github.io/MiMA/api/vert_turb_driver_mod/) |
| `z_full` | lon, lat, pfull | m | geopotential height relative to surface at full levels |  | [`vert_turb_driver_mod`](https://eddy-stanford.github.io/MiMA/api/vert_turb_driver_mod/) |
| `z_half` | lon, lat, phalf | m | geopotential height relative to surface at half levels |  | [`vert_turb_driver_mod`](https://eddy-stanford.github.io/MiMA/api/vert_turb_driver_mod/) |
| `z_pbl` | lon, lat | m | depth of planetary boundary layer |  | [`vert_turb_driver_mod`](https://eddy-stanford.github.io/MiMA/api/vert_turb_driver_mod/) |

<!-- mimadoc:end -->
