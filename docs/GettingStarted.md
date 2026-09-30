[back to contents](README.md)

# Getting started with MiMA

This page explains how to compile MiMA and run the test case that ships with the repository.

* [Downloading the source](#downloading-the-source)
* [Dependencies](#dependencies)
  * [FMS](#fms)
  * [Installing FRE-NCtools](#installing-fre-nctools)
* [Compiling](#compiling)
* [Running the test case](#running-the-test-case)
* [Output](#output)
* [Restarting a run](#restarting-a-run)
* [Adding files to the build](#adding-files-to-the-build)

## Downloading the source

Clone the repository from GitHub:

```bash
git clone https://github.com/Eddy-Stanford/MiMA.git
cd MiMA
```

Tagged releases are listed on the [releases page](https://github.com/Eddy-Stanford/MiMA/releases) and in the [version history](Versions.md). MiMA is free to use, but please cite the relevant [references](README.md#references) in any publication that uses it.

## Dependencies

MiMA needs:

* a Fortran and a C compiler: GNU (`gfortran`/`gcc`/`clang`) or Intel oneAPI (`ifx`/`icx`, or the classic `ifort`)
* an MPI library (e.g. Open MPI, MPICH, Intel MPI) with the Fortran `mpi_f08` module
* netCDF, **both** the C library and the Fortran library (`netcdf-c` and `netcdf-fortran`)
* CMake ≥ 3.22
* the [FMS](https://github.com/NOAA-GFDL/FMS) library, release 2026.02 or later. You don't need to install it: CMake downloads and builds it if it can't find it (see [FMS](#fms)).

To combine the per-processor output files you will also need `mppnccombine` from FRE-NCtools, which is installed separately (see [Installing FRE-NCtools](#installing-fre-nctools)).

Typical ways to install them:

| Platform | Command |
|---|---|
| macOS (Homebrew) | `brew install gcc open-mpi netcdf netcdf-fortran cmake` |
| Ubuntu / Debian | `sudo apt install gfortran libopenmpi-dev openmpi-bin libnetcdf-dev libnetcdff-dev cmake` |
| HPC cluster | load the equivalent modules, e.g. `module load gcc openmpi netcdf-c netcdf-fortran cmake` (names vary between systems) |

CMake finds netCDF using `nc-config`/`nf-config` on your `PATH`. If netCDF is installed somewhere non-standard, point CMake at it by setting `NetCDF_ROOT` (or `NetCDF_C_ROOT` and `NetCDF_Fortran_ROOT` if they are installed separately), either as an environment variable or with `-DNetCDF_ROOT=/path/to/netcdf`.

### FMS

MiMA is built on NOAA-GFDL's Flexible Modeling System (FMS) library, which provides the parallel infrastructure, I/O, diagnostics and time management. CMake looks for an installed FMS 2026.02 or later that was built with 8-byte reals (FMS's `-D64BIT=ON`, which provides the `FMS::fms_r8` target). To use one, point CMake at its install prefix:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/fms
```

(or set `FMS_ROOT=/path/to/fms`). If none is found, CMake downloads the FMS 2026.02 source from GitHub during the configure step and builds it with MiMA. This needs network access the first time. Where there is none, e.g. on some HPC compute nodes, download [FMS 2026.02](https://github.com/NOAA-GFDL/FMS/archive/refs/tags/2026.02.tar.gz) elsewhere, unpack it, and pass its location:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DFETCHCONTENT_SOURCE_DIR_FMS=/path/to/FMS-2026.02
```

A downloaded FMS is built with OpenMP only if `MIMA_OPENMP` is on (see [Compiling](#compiling)), and it is not installed by `cmake --install`.

An installed FMS must also use FMS's default GFDL physical constants (`-DCONSTANTS=GFDL`, the default; Spack: `constants=GFDL`). MiMA checks this when it starts and stops if, for example, the GFS constants were chosen. It also needs `do_simple = .true.` in `&sat_vapor_pres_nml` (set in the shipped `input.nml`), and stops if it is missing.

With the Intel compilers, a downloaded FMS is compiled with FMS's own Intel flags, not MiMA's, so answers will not match builds that used MiMA's old bundled FMS.

### Installing FRE-NCtools

MiMA writes one output file per MPI process (see [Output](#output)). You join them with `mppnccombine`, which is part of NOAA-GFDL's [FRE-NCtools](https://github.com/NOAA-GFDL/FRE-NCtools) and is not included with MiMA. FRE-NCtools also provides `plevel.sh` for interpolating output to pressure levels.

FRE-NCtools isn't available from Homebrew, apt or conda-forge, so we recommend building it from source. It needs the same compilers and netCDF libraries as MiMA, plus `autoconf` and `automake` (`brew install autoconf automake` on macOS, `sudo apt install autoconf automake` on Ubuntu/Debian). To build it and install it into, for example, `~/fre-nctools`:

```bash
git clone --branch 2026.01.01 https://github.com/NOAA-GFDL/FRE-NCtools.git
cd FRE-NCtools
autoreconf -i
mkdir build && cd build
../configure --prefix=$HOME/fre-nctools
make -j 8
make install
```

`2026.01.01` is the latest release at the time of writing; see the [FRE-NCtools releases](https://github.com/NOAA-GFDL/FRE-NCtools/releases) for newer ones. Then add the tools to your `PATH` (e.g. in `~/.bashrc` or `~/.zshrc`):

```bash
export PATH=$HOME/fre-nctools/bin:$PATH
```

Check that it worked with `which mppnccombine`. The same steps work on HPC systems once the compiler and netCDF modules are loaded. See the [FRE-NCtools README](https://github.com/NOAA-GFDL/FRE-NCtools#readme) for more build options.

### Using the container

If you would rather not install the dependencies yourself, a development container with GNU compilers, Open MPI and netCDF is available as the Docker image `robcking/eddy_builder_dev:gnu_openmpi`. The repository's [`.devcontainer`](https://github.com/Eddy-Stanford/MiMA/blob/master/.devcontainer/devcontainer.json) configuration uses this image, so opening the repository in VS Code (with the Dev Containers extension) or GitHub Codespaces gives you a ready-to-build environment. Then follow the compile steps below as normal.

## Compiling

MiMA is built with CMake. From the top of the repository:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DINSTALL_EXEC=ON
cmake --build build -j 8
cmake --install build
```

The first command configures the build in the `build/` directory, the second compiles (using 8 parallel jobs), and the third installs the executable. The build produces the model executable `build/mima`.

The CMake options are:

| Option | Default | Meaning |
|---|---|---|
| `CMAKE_BUILD_TYPE` | `Debug` | Use `Release` for production runs (enables compiler optimisation). `Debug` builds are much slower. |
| `INSTALL_EXEC` | `OFF` | If `ON`, `cmake --install` creates a ready-to-run test case in `exec/` in the repository (see [below](#running-the-test-case)). If `OFF`, the executable is installed to `<prefix>/bin` in the usual CMake way (set the prefix with `-DCMAKE_INSTALL_PREFIX=...`). |
| `MIMA_OPENMP` | `OFF` on macOS, `ON` elsewhere | Compile with OpenMP. MiMA itself has no OpenMP code, so this matters only for a downloaded FMS, which then needs OpenMP for both C and Fortran (Apple's `clang` has none). |

### Choosing the compiler

CMake picks up the compilers from the `FC` and `CC` environment variables. The GNU and Intel compilers are both supported, and the correct compiler flags are chosen automatically. To use the Intel compilers, for example:

```bash
FC=ifx CC=icx cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DINSTALL_EXEC=ON
```

If you change compilers, delete the `build/` directory first so CMake starts afresh.

## Running the test case

With `-DINSTALL_EXEC=ON`, `cmake --install build` creates a run directory `exec/` containing everything needed for the test case:

```
exec/
├── mima            # model executable
├── input.nml       # namelists: all model parameters
├── diag_table      # which diagnostics to write, and how often
├── field_table     # tracers to advect
├── INPUT/          # input data (ozone, topography, land-sea mask)
└── RESTART/        # restart files are written here at the end of the run
```

These files are copied from the [`input/`](https://github.com/Eddy-Stanford/MiMA/tree/master/input) directory of the repository. **Note:** re-running `cmake --install build` overwrites `input.nml`, `diag_table` and `field_table` in `exec/`. For your own experiments, copy `exec/` (or `input/` plus the executable) to a separate run directory.

To run the model:

```bash
cd exec
ulimit -s unlimited
mpirun -n 4 ./mima
```

* `ulimit -s unlimited` removes the limit on the stack size, which MiMA needs. macOS doesn't allow an unlimited stack, so use `ulimit -s hard` there instead.
* `-n 4` sets the number of MPI processes. The number of processes must divide the number of latitudes (64 at T42) evenly. Use `mpiexec`, `srun`, etc., as appropriate on your system.
* MiMA reads `input.nml` automatically, so don't pass it on the command line (i.e. don't do `./mima < input.nml`).

As a rough guide, the test case runs at about 10 s per model day on 4 cores of a laptop, so the full year takes around an hour. To try things out more quickly, reduce `days` in `coupler_nml` (e.g. `days = 5`). The run length should be a whole multiple of the time step `dt_atmos` (500 s in the test case).

### The test case

The test case is defined entirely by the files in [`input/`](https://github.com/Eddy-Stanford/MiMA/tree/master/input):

* `input.nml`: This is the most important file. It sets all the input parameters within the various namelists of MiMA. Any variable not present in `input.nml` takes its (hard-coded) default value. This file completely defines the simulation you are running. See [Parameter settings](Parameters.md) for what the parameters mean.
* `diag_table`: A list of the diagnostics you would like in your output files. It doesn't change the simulation you are running. It only decides which variables are written, how frequently, and whether the output is averaged or instantaneous.
* `field_table`: A list of passive tracers you'd like to advect during the simulation. There are two types: grid or spectral tracers. To get the temporal evolution of a tracer (or its time average), add its name as a diagnostic output in `diag_table`.

The test run is one 360-day year (12 months of 30 days) with the following setup:
* T42 horizontal resolution (128 × 64 grid) with 40 vertical levels
* realistic topography and land-sea mask, interpolated from `INPUT/navy_topography.data.nc` and `INPUT/navy_pctwater.data.nc`
* RRTM radiation scheme, with 390 ppm CO<sub>2</sub>, ozone from `INPUT/ozone_1990.nc`, and a solar constant of 1370 W/m<sup>2</sup>
* seasonal cycle with a circular Earth-Sun orbit
* mixed-layer ocean with a meridional Q flux, plus zonally asymmetric Q fluxes (tropical warm pool, Gulf Stream, Kuroshio, …) to generate realistic stationary waves as in [Garfinkel et al. (2020)](https://doi.org/10.1175/JCLI-D-19-0181.1)
* surface albedo of 0.23, increasing to 0.8 in polar regions, with brighter Sahara, Gobi and Australian deserts
* Betts-Miller convection and large-scale condensation
* parameterized non-orographic gravity-wave drag (`cg_drag`)

## Output

MiMA writes one output file per MPI process for each file listed in `diag_table`, e.g. `atmos_daily.nc.0000`, `atmos_daily.nc.0001`, …. Each file holds a band of latitudes. Combine them into a single netCDF file with `mppnccombine` (see [Installing FRE-NCtools](#installing-fre-nctools)):

```bash
for f in atmos_daily atmos_avg atmos_davg atmos_dext; do
    mppnccombine -r $f.nc $f.nc.????
done
```

The `-r` flag removes the per-process files once they've been combined successfully. Run `mppnccombine` without arguments for its other options.

The test case produces:

| File | Contents |
|---|---|
| `atmos_daily.nc` | daily instantaneous wind, temperature, humidity and surface pressure |
| `atmos_avg.nc` | wind, temperature and humidity averaged over the whole run |
| `atmos_davg.nc` | daily-mean surface temperature and precipitation |
| `atmos_dext.nc` | daily maximum/minimum surface temperature and maximum precipitation |

The output is on the model's hybrid sigma levels. To interpolate it to pressure levels, use `plevel.sh` from FRE-NCtools on a combined file:

```bash
plevel.sh -a -i atmos_daily.nc -o atmos_daily_plev.nc
```

`-a` interpolates all fields. By default the output is on 17 standard levels from 1000 to 10 hPa; use `-p "100000 85000 ..."` to choose your own levels (in Pa). Run `plevel.sh` without arguments for all options. The interpolation needs `pk`, `bk` and `ps` in the file, which is why the test case writes them to `atmos_daily` (the pressure at the level interfaces is `pk + bk * ps`).

## Restarting a run

At the end of a run, MiMA writes its final state to `RESTART/`. To continue the simulation, move the restart files to `INPUT/` and run again:

```bash
mv RESTART/* INPUT/
mpirun -n 4 ./mima
```

The model detects the restart files in `INPUT/` and continues from the date stored in `INPUT/coupler.res`. Move or combine the output files from the previous segment first, because the new run overwrites them. Long simulations are usually run as a sequence of such segments, e.g. one year at a time.

## Adding files to the build

* If you work on your own version of MiMA, put each extension in a new file where possible, so as not to disturb the main branch and any other fork that might exist.
* When adding a source file, add it to the `CMakeLists.txt` in the same directory, so that it is compiled the next time you build.
