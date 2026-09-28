[back to contents](README.md)

# Postprocessing

The [`postprocessing/`](https://github.com/Eddy-Stanford/MiMA/tree/master/postprocessing) directory contains tools for working with MiMA output.

## Combining output files: `mppnccombine`

MiMA writes one output file per MPI process (`atmos_daily.nc.0000`, `atmos_daily.nc.0001`, …). `mppnccombine` joins them into a single netCDF file. It is built by CMake when you configure with `-DBUILD_COMBINE=ON`, and installed next to `mima` (in `exec/` with `-DINSTALL_EXEC=ON`, otherwise in `<prefix>/bin`).

```bash
./mppnccombine -r atmos_daily.nc atmos_daily.nc.????
```

Options:

* `-r`: remove the per-process files after a successful combine
* `-a`: append to an existing output file
* `-v`: print progress information
* `-n #`: input file extensions start at `#` instead of `0000`

Run `./mppnccombine` without arguments for the full usage message.

## Interpolating to pressure levels: `plevel_interpolation`

MiMA's output is on the model's (hybrid) sigma levels. The GFDL `plevel` tool in [`postprocessing/plevel_interpolation/`](https://github.com/Eddy-Stanford/MiMA/tree/master/postprocessing/plevel_interpolation) interpolates it to pressure levels. It needs `bk`, `pk` and `ps` in the output file, which is why the test case's `diag_table` writes these to `atmos_daily`.

This tool is **not** built by CMake. It still uses the older `mkmf` build system: see [`postprocessing/plevel_interpolation/README`](https://github.com/Eddy-Stanford/MiMA/blob/master/postprocessing/plevel_interpolation/README) for how to compile `plev.x` and run it through `scripts/plevel.sh`. Alternatively, many users interpolate to pressure levels in Python (e.g. with `xarray`), using `p = pk + bk * ps` on the half levels.

## Turning output into input: `output_to_input.py`

[`postprocessing/output_to_input.py`](https://github.com/Eddy-Stanford/MiMA/blob/master/postprocessing/output_to_input.py) converts a MiMA output file into a file MiMA can read back as input (for example as an input field in `INPUT/`). It can select variables, add the cell-boundary coordinates `lonb`/`latb` from another file, duplicate the first or last time record, and compress the result. Run `python output_to_input.py -h` for the options.

Note that this script was written for Python 2 and uses `xarray`. It needs minor changes (the `print` statements) to run under Python 3.
