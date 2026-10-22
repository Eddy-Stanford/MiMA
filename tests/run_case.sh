#!/usr/bin/env bash
#
# Run one MiMA test case and reduce its output to a digest. See tests/README.md.
#
#   run_case.sh run     <mima> <case> <out_dir> <npes> [<reference.digest>]
#   run_case.sh restart <mima> <case> <out_dir> <npes> <s1>,<s2>[,...]
#   run_case.sh guard   <mima> <case> <out_dir> <name> <edit>
#   run_case.sh compare <a.digest> <b.digest> [<file prefix>]
#
# run: run tests/cases/<case> and write <out_dir>/digests/<case>.digest (and .log, the
#   last 50 lines of the output); with a reference digest, the two must be identical.
# restart: run in segments of s1, s2, ... seconds, each from the last one's restart,
#   and once for their sum; the RESTART/ digests must be identical.
# guard: one edit, <group>:<key>=<value> (namelist), touch:<path> or rm:<path> (run
#   directory), then one step on 1 rank; the edit should make the model stop.
# compare: list the variables whose hashes differ; exit 1 if any do.
#
# Disk and CPU rules (tests/README.md), enforced here: output only in the build tree
# (the directory of <mima>) or $TMPDIR; INPUT/ files are links; a passing run directory
# is deleted, a failed one blocks the next run until it is removed; at most
# MIMA_TEST_MAX_MB (300) MB and MIMA_TEST_TIMEOUT (900) s per run; npes <= physical
# cores; one run_case.sh per user at a time (a PID lock file, never waited for).

set -euo pipefail
shopt -s nullglob

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
MAX_MB=${MIMA_TEST_MAX_MB:-300}
TIME_LIMIT=${MIMA_TEST_TIMEOUT:-900}
LOCK="${TMPDIR:-/tmp}"; LOCK="${LOCK%/}/mima-tests-${USER:-$(id -un)}.lock"
MPI_PID=""

die() { echo "run_case: $*" >&2; exit 1; }

cleanup() {
  [[ -z "$MPI_PID" ]] || kill "$MPI_PID" 2> /dev/null || true
  [[ "$(cat "$LOCK" 2> /dev/null)" != "$$" ]] || rm -f "$LOCK"
}

# One run per user at a time: refuse (do not wait) while the lock's PID is alive.
take_lock() {
  local pid
  if ! (set -o noclobber; echo $$ > "$LOCK") 2> /dev/null; then
    pid=$(cat "$LOCK" 2> /dev/null || true)
    if [[ -n "$pid" ]] && kill -0 "$pid" 2> /dev/null; then
      die "another MiMA test run (PID $pid) is active; not starting a second one (lock $LOCK)"
    fi
    echo "run_case: removing the stale lock $LOCK (PID $pid has ended)" >&2
    rm -f "$LOCK"
    (set -o noclobber; echo $$ > "$LOCK") 2> /dev/null || die "cannot take the lock $LOCK"
  fi
  trap cleanup EXIT
  trap 'exit 130' INT TERM
}

physical_cores() { sysctl -n hw.physicalcpu 2> /dev/null || nproc; }

abspath() {  # physical absolute path of $1, which need not exist yet
  local p="$1" rest=""
  while [[ ! -d "$p" ]]; do rest="/$(basename "$p")$rest"; p="$(dirname "$p")"; done
  echo "$(cd "$p" && pwd -P)$rest"
}

# Common start of every run: check the arguments, take the lock, check the output dir.
setup() {  # setup <mima> <out_dir> <npes>
  [[ -x "$1" ]] || die "not an executable: $1"
  EXE="$(abspath "$1")"
  [[ "$3" =~ ^[1-9][0-9]*$ ]] || die "npes must be a positive integer, not '$3'"
  (( $3 <= $(physical_cores) )) || die "$3 ranks exceed the $(physical_cores) physical cores"
  local build tmp; build="$(abspath "$(dirname "$EXE")")"; tmp="$(abspath "${TMPDIR:-/tmp}")"
  OUT="$(abspath "$2")"
  case "$OUT/" in
    "$build"/*|"$tmp"/*) ;;
    *) die "output dir $OUT is neither in the build tree ($build) nor in \$TMPDIR ($tmp)" ;;
  esac
  take_lock
  mkdir -p "$OUT/run" "$OUT/digests"
  local old; for old in "$OUT"/run/*; do
    die "a failed run is still in $old ($(du -sh "$old" | cut -f1)); look at it, then: rm -rf '$old'"
  done
}

# nml_set <file> <group> <key> <value>: set a scalar key in a namelist group,
# removing any earlier assignment of it; the group is appended if it is absent.
nml_set() {
  awk -v g="$2" -v k="$3" -v val="$4" '
    function scan(s,   i, c, q) {  # ncode: length before a ! comment; closed: a / ends the group
      q = ""; closed = 0
      for (i = 1; i <= length(s); i++) {
        c = substr(s, i, 1)
        if (q != "") { if (c == q) q = "" }
        else if (c == "\047" || c == "\"") q = c
        else if (c == "!") break
        else if (c == "/") closed = 1
      }
      ncode = i - 1
    }
    function drop(s,   re, t) {  # remove every "key = value[,]" from s
      re = "[^a-z0-9_%]" tolower(k) "[ \t]*=[ \t]*(\047[^\047]*\047|\"[^\"]*\"|[^ \t,/]+)[ \t]*,?"
      t = " " s
      while (match(tolower(t), re)) t = substr(t, 1, RSTART) substr(t, RSTART + RLENGTH)
      return substr(t, 2)
    }
    BEGIN { G = "&" tolower(g) }
    !ing && !done {
      s = $0; sub(/^[ \t]+/, "", s); c = substr(s, length(G) + 1, 1)
      if (tolower(substr(s, 1, length(G))) != G || (c != "" && c != " " && c != "\t")) { print; next }
      ing = done = 1; p = index(tolower($0), G) + length(G) - 1
      head = substr($0, 1, p) " " k " = " val ","; $0 = substr($0, p + 1)
    }
    !ing { print; next }
    { scan($0); print head drop(substr($0, 1, ncode)) substr($0, ncode + 1); head = ""
      if (closed) ing = 0 }
    END { if (!done) printf "\n&%s %s = %s /\n", g, k, val }' "$1" > "$1.new"
  mv "$1.new" "$1"
}

set_run_length() {  # set_run_length <input.nml> <seconds>
  local k; for k in months days hours minutes; do nml_set "$1" coupler_nml "$k" 0; done
  nml_set "$1" coupler_nml seconds "$2"
}

prepare() {  # prepare <run_dir> <case>: a run directory; files missing from a case come from rrtm
  local d="$1" c="$HERE/cases/$2" f
  [[ -f "$c/input.nml" ]] || die "no such case: $c"
  mkdir -p "$d/INPUT" "$d/RESTART"
  for f in "$REPO"/input/INPUT/*; do ln -s "$f" "$d/INPUT/"; done
  for f in input.nml diag_table field_table; do
    if [[ -f "$c/$f" ]]; then cp "$c/$f" "$d/"; else cp "$HERE/cases/rrtm/$f" "$d/"; fi
  done
}

# run_model <run_dir> <npes> <log>: mpirun in the background, killed if the run
# directory grows past MAX_MB or the run takes longer than TIME_LIMIT.
run_model() {
  local d="$1" t0=$SECONDS kb rc=0
  ( cd "$d" || exit 1; ulimit -s hard 2> /dev/null || true; export OMP_NUM_THREADS=1
    exec "${MIMA_MPIRUN:-mpirun}" -n "$2" "$EXE" ) > "$3" 2>&1 &
  MPI_PID=$!
  while kill -0 "$MPI_PID" 2> /dev/null; do
    sleep 1
    kb=$(du -sk "$d" | cut -f1)
    if (( kb > MAX_MB * 1024 )); then
      kill "$MPI_PID" 2> /dev/null; wait "$MPI_PID" || true
      du -sk "$d"/* "$d"/RESTART/* | sort -n | tail -5
      rm -rf "$d"; die "output too large: $d passed $MAX_MB MB (largest files above); deleted it"
    fi
    if (( SECONDS - t0 > TIME_LIMIT )); then
      kill "$MPI_PID" 2> /dev/null; wait "$MPI_PID" || true; die "killed after $TIME_LIMIT s; output kept in $d"
    fi
  done
  wait "$MPI_PID" || rc=$?
  MPI_PID=""
  return $rc
}

if command -v md5sum > /dev/null; then HASH=md5sum; else HASH=md5; fi

# digest <run_dir>: "<file> <variable> <md5>" for every variable of every netCDF file;
# the md5 is over the variable's data as printed by ncdump. Fails on NaN or Inf.
digest() {
  local d="$1" f v t bad=0
  t="$(mktemp "${TMPDIR:-/tmp}/mima-digest.XXXXXX")"
  for f in "$d"/*.nc "$d"/RESTART/*.nc; do
    for v in $(ncdump -h "$f" | awk '/^variables:/ { on = 1; next } /^(data:|\/\/ global)/ { on = 0 }
                                    on && NF >= 2 && $1 !~ /:/ { sub(/\(.*|;/, "", $2); print $2 }'); do
      ncdump -v "$v" -p 9,17 "$f" | sed -n '/^data:/,$p' > "$t" || { echo "ncdump failed: $f $v" >&2; bad=1; }
      if grep -qF -e NaN -e Infinity "$t"; then echo "run_case: NaN or Inf in ${f#"$d"/} $v" >&2; bad=1; fi
      echo "${f#"$d"/} $v $($HASH < "$t" | cut -c1-32)"
    done
  done > "$t.all"
  LC_ALL=C sort "$t.all"
  rm -f "$t" "$t.all"
  return $bad
}

cmd_compare() {  # compare <a.digest> <b.digest> [<file prefix>]
  local a="$1" b="$2" diffs
  [[ -f "$a" ]] || die "no digest $a"
  [[ -f "$b" ]] || die "no digest $b"
  diffs=$(awk -v p="${3:-}" '
    p != "" && index($1, p) != 1 { next }
    NR == FNR { ha[$1 " " $2] = $3; next }
    { k = $1 " " $2; seen[k] = 1
      if (!(k in ha)) print "  only in B: " k; else if (ha[k] != $3) print "  differs:   " k }
    END { for (k in ha) if (!(k in seen)) print "  only in A: " k }' "$a" "$b" | LC_ALL=C sort)
  if [[ -n "$diffs" ]]; then
    echo "DIFFER  A = $a"; echo "        B = $b"; echo "$diffs"
    echo "$(echo "$diffs" | wc -l | tr -d ' ') variable(s) differ"
    return 1
  fi
  echo "IDENTICAL  $(grep -c "^${3:-}" "$b") variables: $a = $b"
}

# run_segments <case> <name> <npes> [<s1>,<s2>,...]: run, digest, clean up.
run_segments() {
  local d="$OUT/run/$2" log="$OUT/digests/$2.log" segs="${4:--}" k=0 s len t0=$SECONDS
  rm -f "$OUT/digests/$2.digest"   # no stale digest if this run fails
  prepare "$d" "$1"
  : > "$log"
  for s in ${segs//,/ }; do
    k=$((k + 1))
    len="length from input.nml"
    if [[ $s != - ]]; then set_run_length "$d/input.nml" "$s"; len="$s s"; fi
    if (( k > 1 )); then mv "$d"/RESTART/* "$d/INPUT/"; rm -f "$d"/*.nc; fi
    if ! run_model "$d" "$3" "$d/stdout.$k"; then
      tail -30 "$d/stdout.$k"; die "$2 failed in segment $k; output kept in $d"
    fi
    { echo "== $2, segment $k ($len)"; tail -50 "$d/stdout.$k"; } >> "$log"
  done
  digest "$d" > "$OUT/digests/$2.digest" || die "$2 wrote NaN or Inf; output kept in $d"
  echo "$2: $(wc -l < "$OUT/digests/$2.digest" | tr -d ' ') variables in $((SECONDS - t0)) s" \
       "($(du -sh "$d" | cut -f1) of output) -> $OUT/digests/$2.digest"
  rm -rf "$d"
}

cmd_guard() {  # guard <mima> <case> <out_dir> <name> <edit>
  setup "$1" "$3" 1
  local d="$OUT/run/guard_$4" e="$5" rc=0
  prepare "$d" "$2"
  set_run_length "$d/input.nml" 500
  case "$e" in
    touch:*) : > "$d/${e#touch:}" ;;
    rm:*)    [[ $e != *..* ]] || die "bad edit '$e'"; rm -rf "${d:?}/${e#rm:}" ;;
    *:*=*)   local g="${e%%:*}" kv="${e#*:}"; nml_set "$d/input.nml" "$g" "${kv%%=*}" "${kv#*=}" ;;
    *)       die "unknown edit '$e'" ;;
  esac
  TIME_LIMIT=120
  run_model "$d" 1 "$d/stdout" || rc=$?
  grep -a -E 'FATAL|ERROR' "$d/stdout" | uniq | head -5 || true
  echo "--- last lines of the output:"; tail -15 "$d/stdout" | tr -d '\000'
  rm -rf "$d"
  (( rc != 0 )) || die "guard_$4: the model ran to the end; the guard did not stop it"
  echo "guard_$4: the model stopped (exit $rc)"
}

sub="${1:-}"; shift || true
case "$sub" in
  run)     (( $# == 4 || $# == 5 )) || die "usage: run <mima> <case> <out_dir> <npes> [<reference.digest>]"
           setup "$1" "$3" "$4"; run_segments "$2" "$2" "$4"
           if (( $# == 5 )); then cmd_compare "$5" "$OUT/digests/$2.digest"; fi ;;
  restart) (( $# == 5 )) || die "usage: restart <mima> <case> <out_dir> <npes> <s1>,<s2>[,...]"
           setup "$1" "$3" "$4"
           total=0; for s in ${5//,/ }; do total=$((total + s)); done
           run_segments "$2" "${2}_restart" "$4" "$5"
           run_segments "$2" "${2}_continuous" "$4" "$total"
           cmd_compare "$OUT/digests/${2}_continuous.digest" "$OUT/digests/${2}_restart.digest" RESTART/ ;;
  guard)   (( $# == 5 )) || die "usage: guard <mima> <case> <out_dir> <name> <edit>"; cmd_guard "$@" ;;
  compare) (( $# == 2 || $# == 3 )) || die "usage: compare <a.digest> <b.digest> [<file prefix>]"
           cmd_compare "$@" ;;
  *)       sed -n '3,8p' "$0"; exit 2 ;;
esac
