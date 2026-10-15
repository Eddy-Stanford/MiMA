# MiMA [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.597136.svg)](https://doi.org/10.5281/zenodo.597136)
Model of an idealized Moist Atmosphere

MiMA is an intermediate-complexity General Circulation Model with interactive water vapor and full radiation. It is based on the gray radiation model of [Frierson, Held, and Zurita-Gotor, JAS (2006)](https://doi.org/10.1175/JAS3753.1), which it still contains as a namelist option. The major development in MiMA is replacing the gray radiation scheme with a full radiative transfer code, the Rapid Radiative Transfer Model ([RRTM](http://rtweb.aer.com/rrtm_frame.html)) developed by AER.

MiMA is publicly available, but users are asked to cite the appropriate [references](#references) in any publication resulting from its use.

The source code is on [GitHub](https://github.com/Eddy-Stanford/MiMA).

## Contents

* [Getting started](GettingStarted.md): dependencies, compiling, running the test case
* [Model configurations](Configurations.md): radiation schemes, specified initial conditions, initial-condition noise
* [Parameter settings](Parameters.md): default and recommended namelist values
* [Diagnostics](Diagnostics.md): how to write a `diag_table`, and every diagnostic field the model can output
* [Fortran API reference](FortranAPI.md): the modules and procedures of the source code, and how to document them
* [Version history](Versions.md): main additions and changes
* [Migrating from v1 to v2.0](Migration_v2.md): what changed in v2.0 and how to convert a v1 setup
* [References](#references): required and relevant references
* [License](https://github.com/Eddy-Stanford/MiMA#license)

## References

MiMA

* [Jucker and Gerber, J Clim (2017)](https://doi.org/10.1175/JCLI-D-17-0127.1)
* [Garfinkel et al., J Clim (2020)](https://doi.org/10.1175/JCLI-D-19-0181.1)
* code:
    * latest version: [![DOI](https://zenodo.org/badge/36012278.svg)](https://zenodo.org/badge/latestdoi/36012278)
    * v1.1: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3637607.svg)](https://doi.org/10.5281/zenodo.3637607)
    * v1.0: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.321708.svg)](https://doi.org/10.5281/zenodo.321708)

Gray radiation model

* [Frierson, Held, Zurita-Gotor, JAS (2006)](https://doi.org/10.1175/JAS3753.1)
* [Frierson, JAS (2007)](https://doi.org/10.1175/JAS3935.1)
* [Frierson, Held, Zurita-Gotor, JAS (2007)](https://doi.org/10.1175/JAS3913.1)

RRTM

* [Mlawer et al., JGR (1997)](https://doi.org/10.1029/97JD00237)
* [Iacono et al., JGR (2000)](https://doi.org/10.1029/2000JD900091)
* [Iacono et al., JGR (2008)](https://doi.org/10.1029/2008JD009944)
* [Clough et al., JQSRT (2005)](https://doi.org/10.1016/j.jqsrt.2004.05.058)

Isca framework (which contains MiMA)

* [Vallis et al., GMD (2018)](https://doi.org/10.5194/gmd-11-843-2018)
