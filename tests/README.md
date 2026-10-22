# MiMA tests

A few short runs of the model that check what can be checked on any machine.
Every test compares the model with itself (or with an expected error message),
so nothing machine-specific is stored in git. The full configurations are in
`examples/` and `input/`; the cases here are copies of them with a 10-hour run.

| Test | Label | Proves | Passes when |
|---|---|---|---|
| `smoke_rrtm`, `smoke_gray`, `smoke_hs` | smoke | the shipped configurations run (`input/input.nml`, `examples/gray`, `examples/held_suarez`) | exit 0, a digest is written, no NaN or Inf in any output |
| `restart_rrtm` | self | a restart continues a run exactly: 2 x 10000 s against one 20000 s run (20 steps is not a multiple of the 9-step radiation interval) | the `RESTART/` digests are identical |
| `defaults_rrtm` | self | the code defaults are the shipped configuration (`cases/rrtm_defaults` sets only what has no usable default) | its digest is identical to `smoke_rrtm`'s |
| `guard_<name>` | guard | bad input stops with a clear message (one line of `guards.txt` each, 1 rank) | the message appears in the output |
| `compare_<case>` | ref | bit for bit against another build (only with `-DMIMA_TEST_REFERENCE`) | the digests are identical |

A **digest** (`<build>/tests/digests/<case>.digest`) has one line
`<file> <variable> <md5>` for every variable of every netCDF file the run wrote,
restarts included. The md5 is taken over the data printed by
`ncdump -v <variable> -p 9,17`, which shows every float and double exactly, so equal
digests mean bit-for-bit equal output. A digest is valid only for one compiler,
flags, FMS, MPI library and rank count, which is why none is committed.

## Running

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DMIMA_TESTING=ON
cmake --build build -j 4
ctest --test-dir build --output-on-failure        # all tests, under 2 minutes on 4 ranks
ctest --test-dir build -L smoke                   # only one label
ctest --test-dir build -R guard_ -V               # some tests, with their output
```

`MIMA_TEST_NPES` (default 4) sets the ranks per model run; it is capped at the
number of physical cores (`-DMIMA_TEST_CORES=logical` caps at the logical cores instead,
for small CI runners whose two cores show four threads). The tests need `mpirun`, `ncdump`, `md5sum` or `md5`,
bash and awk. `ctest` only calls `tests/run_case.sh`, which you can also run by hand:

```sh
tests/run_case.sh run build/mima rrtm build/tests 4       # -> build/tests/digests/rrtm.digest
tests/run_case.sh compare a.digest b.digest               # lists the differing variables
```

## Disk and CPU rules

`run_case.sh` enforces these; they are not just advice.

* Output goes to `<build>/tests/` (or anywhere under `$TMPDIR`); any other
  directory is refused. `input/INPUT/*` is linked into the run directory, not copied.
* A passing run leaves only its digest and the last 50 lines of its output
  (`<case>.log`). While it runs, a run directory holds about 120 MB (200 MB in
  the restart test).
* A failed run keeps its directory in `<build>/tests/run/` for you to inspect, and
  every later run refuses to start until you remove it (the message gives the
  `rm -rf` command). So at most one failed run is ever on disk.
* A run directory larger than 300 MB (`MIMA_TEST_MAX_MB`) is killed and deleted
  with "output too large"; a run longer than 900 s (`MIMA_TEST_TIMEOUT`) is killed.
* The ranks may not exceed the physical cores. Only one `run_case.sh` per user runs
  at a time: a second one, from another shell or another `ctest`, refuses at once
  (lock file `$TMPDIR/mima-tests-$USER.lock`) instead of waiting or overloading the
  machine. A lock left by a killed run is removed automatically. All tests share a
  CTest `RESOURCE_LOCK`, so `ctest -j` runs them one after another.

## Comparing with another commit

To check that a refactoring changes no answer, build the reference commit in
`$TMPDIR`, digest its runs with its own runner and cases, then compare:

```sh
ref=$TMPDIR/mima-ref && mkdir -p "$ref/src" && git archive <commit> | tar -x -C "$ref/src"
cmake -S "$ref/src" -B "$ref/build" -DCMAKE_BUILD_TYPE=Release -DMIMA_TESTING=ON
cmake --build "$ref/build" -j 4 && ctest --test-dir "$ref/build" -L smoke
cmake -B build -DMIMA_TEST_REFERENCE="$ref/build/tests/digests" && ctest --test-dir build -L ref
rm -rf "$ref"
```

Use the same compiler, flags and `MIMA_TEST_NPES` for both builds. For a commit
older than `tests/`, digest its binary with this tree's runner instead
(`tests/run_case.sh run "$ref/build/mima" rrtm "$ref/build/tests" 4`, and the same
for `gray` and `hs`), provided its namelists accept the same keys. After an
intended change of answers there is nothing to update: the reference is whatever
commit you name.

## Adding to the tests

* **A guard**: one line in `guards.txt` (the format is in its header); re-run cmake.
* **A case**: a directory in `cases/` with `input.nml` (a short run: whole steps,
  ideally 10 hours or less) and, if it differs from rrtm's, `diag_table` and
  `field_table` (missing files are taken from `cases/rrtm/`). Use one `atmos_avg`
  file with output frequency `-1` and packing 1, and keep the outputs small. Add it
  to the `foreach` lists in `CMakeLists.txt`.

## Continuous integration

`.github/workflows/tests.yml` builds MiMA and runs this suite on every push to `v2` and
`main` and on pull requests (FMS is built once and cached). The `ctest` step has a hard
budget of a few minutes (`TEST_BUDGET_MINUTES` in the workflow): if the suite outgrows it,
the job fails and the suite has to be trimmed, not the budget raised. The same recipe in a
4-CPU `ubuntu:24.04` container took 71 s for the tests and about 2 minutes for a cold build.
