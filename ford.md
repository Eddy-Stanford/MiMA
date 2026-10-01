---
project: MiMA
version: 2.0.0
src_dir: ./src
exclude_dir: ./src/atmos_param/radiation/rrtm/rrtmg_lw
             ./src/atmos_param/radiation/rrtm/rrtmg_sw
fpp_extensions: F90
predocmark: >
docmark: !
predocmark_alt: @
docmark_alt: ^
sort: src
graph: false
search: false
extra_mods: fms_mod: https://noaa-gfdl.github.io/FMS/
            fms2_io_mod: https://noaa-gfdl.github.io/FMS/
            mpp_mod: https://noaa-gfdl.github.io/FMS/
            mpp_domains_mod: https://noaa-gfdl.github.io/FMS/
            constants_mod: https://noaa-gfdl.github.io/FMS/
            time_manager_mod: https://noaa-gfdl.github.io/FMS/
            diag_manager_mod: https://noaa-gfdl.github.io/FMS/
            field_manager_mod: https://noaa-gfdl.github.io/FMS/
            tracer_manager_mod: https://noaa-gfdl.github.io/FMS/
            sat_vapor_pres_mod: https://noaa-gfdl.github.io/FMS/
            topography_mod: https://noaa-gfdl.github.io/FMS/
            horiz_interp_mod: https://noaa-gfdl.github.io/FMS/
            time_interp_mod: https://noaa-gfdl.github.io/FMS/
            memutils_mod: https://noaa-gfdl.github.io/FMS/
            platform_mod: https://noaa-gfdl.github.io/FMS/
            netcdf: https://docs.unidata.ucar.edu/netcdf-fortran/current/
---

FORD settings for parsing MiMA's Fortran sources and their doc comments. FORD is used
only as a parser: the MkDocs hook `tools/fortran_api.py` turns the parsed modules into the
"Fortran API reference" pages of the documentation site (see `mkdocs.yml`). `extra_mods`
gives the links for the external modules a MiMA module uses.

The vendored AER RRTMG radiation code (`src/atmos_param/radiation/rrtm/rrtmg_*`) is
not included.
