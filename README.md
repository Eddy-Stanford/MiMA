# Model of an idealized Moist Atmosphere (MiMA) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.597136.svg)](https://doi.org/10.5281/zenodo.597136)

MiMA is an intermediate-complexity General Circulation Model with interactive water vapor and full radiation (RRTM). It is built on the GFDL spectral dynamical core and the gray-radiation moist model of Frierson, Held and Zurita-Gotor (2006).

Full documentation is in the [`docs/`](docs/) folder and online at <https://eddy-stanford.github.io/MiMA/>.

## Quick start

You need a Fortran and C compiler (GNU or Intel), MPI, netCDF (C and Fortran libraries) and CMake ≥ 3.22. MiMA uses the [FMS](https://github.com/NOAA-GFDL/FMS) library, which CMake downloads and builds automatically if it can't find an installed copy. To combine the output you also need [FRE-NCtools](docs/GettingStarted.md#installing-fre-nctools). See [Getting started](docs/GettingStarted.md#dependencies) for how to install them, or use the provided [container](docs/GettingStarted.md#using-the-container).

Compile, and create a ready-to-run test case in `exec/`:

```bash
git clone https://github.com/Eddy-Stanford/MiMA.git
cd MiMA
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DINSTALL_EXEC=ON
cmake --build build -j 8
cmake --install build
```

Run the test case (one 360-day year at T42 resolution with 40 levels):

```bash
cd exec
ulimit -s unlimited   # on macOS use: ulimit -s hard
mpirun -n 4 ./mima
```

MiMA writes each output file listed in `diag_table` (e.g. `atmos_daily.nc`) as a single netCDF file.

See [Getting started](docs/GettingStarted.md) for details on the build options, the test configuration, output and restarts.

## Documentation

* [Getting started](docs/GettingStarted.md): dependencies, compiling, running the test case
* [Model configurations](docs/Configurations.md): radiation schemes, specified initial conditions, initial-condition noise
* [Parameter settings](docs/Parameters.md): default and recommended namelist values
* [Diagnostics](docs/Diagnostics.md): how to write a `diag_table`, and every diagnostic field the model can output
* [Fortran API reference](docs/FortranAPI.md): the modules and procedures of the source code, and how to document them
* [Version history](docs/Versions.md)
* [Migrating from v1 to v2.0](docs/Migration_v2.md)
* [References](docs/README.md#references)

See the 30 second trailer on [YouTube](https://www.youtube.com/watch?v=8UfaFnGtCrk "Model of an idealized Moist Atmosphere (MiMA)"):

[![MiMA thumbnail](https://img.youtube.com/vi/8UfaFnGtCrk/0.jpg)](https://www.youtube.com/watch?v=8UfaFnGtCrk "Model of an idealized Moist Atmosphere (MiMA)")

## Citing MiMA

MiMA is free to use, but please cite the relevant scientific work in any publication that uses it:

* [M Jucker and EP Gerber, 2017: *Untangling the annual cycle of the tropical tropopause layer with an idealized moist model*, Journal of Climate 30, 7339-7358](https://doi.org/10.1175/JCLI-D-17-0127.1) (MiMA v1.x)
* [Garfinkel et al., 2020: *The building blocks of Northern Hemisphere wintertime stationary waves*, Journal of Climate](https://doi.org/10.1175/JCLI-D-19-0181.1) (stationary-wave configuration, formerly referred to as v2.0)

Citation metadata for the code itself is in [`CITATION.cff`](CITATION.cff). Further references for the radiation schemes are listed in the [documentation](docs/README.md#references).

## License

MiMA is distributed under a GNU GPLv3 license. That means you have permission to use, modify, and distribute the code, even for commercial use. However, you must make your code publicly available under the same license. See [LICENSE.txt](LICENSE.txt) for more details.

AM2 is distributed under a GNU GPLv2 license. That means you have permission to use, modify, and distribute the code, even for commercial use. However, you must make your code publicly available under the same license.

RRTM/RRTMG: Copyright © 2002-2010, Atmospheric and Environmental Research, Inc. (AER, Inc.). This software may be used, copied, or redistributed as long as it is not sold and this copyright notice is reproduced on each copy made. This model is provided as is without any express or implied warranties.
