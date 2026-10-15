"""Check a diag_table against the fields registered in the source, and
optionally against the namelist settings of an input.nml."""

import io
import re
import sys

from . import conditions as cnd
from .fortran import split_top

PLACEHOLDER_RE = re.compile(r"\{[^}]*\}|<[^>]*>")
REDUCTIONS = {".true.", "mean", "average", "avg", ".false.", "none", "point", "rms",
              "max", "maximum", "min", "minimum", "sum", "cumsum"}
TIME_UNITS = {"seconds", "minutes", "hours", "days", "months", "years"}
# diag_manager appends these to the output name for these reductions
REDUCTION_SUFFIX = {"max": "max", "maximum": "max", "min": "min", "minimum": "min",
                    "sum": "sum", "cumsum": "sum"}


def name_matcher(pattern, ignore_case):
    """Regex for a registered name; placeholders ({expr}, <tracer>) match anything."""
    parts = PLACEHOLDER_RE.split(pattern)
    rx = "^" + ".+".join(re.escape(p) for p in parts) + "$"
    return re.compile(rx, re.I if ignore_case else 0)


def diag_table_lines(path):
    """(line number, [tokens]) for the non-blank lines of a diag_table, with
    comments after '#' removed and quotes stripped from the tokens."""
    with io.open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()
    for ln, line in enumerate(lines, 1):
        s = line.split("#", 1)[0].strip()
        if not s:
            continue
        toks = [t.strip().strip('"').strip("'").strip() for t in split_top(s, ",")]
        while toks and toks[-1] == "":
            toks.pop()
        yield ln, toks


def validate(path, rows, nml_values=None, defaults=None, out=None):
    """Report problems in a diag_table; returns 1 if there are any, else 0.
    Module names are case sensitive in diag_manager, field names are not.
    With nml_values ({group: {var: value}} from an input.nml) the conditions
    of each field are evaluated; variables not set there take the defaults
    ({(group, var): value})."""
    out = out or sys.stdout
    pats = [(name_matcher(r.module, False), name_matcher(r.field, True), r) for r in rows]
    defaults = defaults or {}
    bad = 0

    def value_of(group, name):
        if nml_values is None or group is None:
            return None
        g = group.lower()
        if name in nml_values.get(g, {}):
            return nml_values[g][name]
        return defaults.get((g, name))

    def report(ln, msg):
        out.write("%s:%d: %s\n" % (path, ln, msg))

    with io.open(path, encoding="utf-8", errors="replace") as fh:
        head = [fh.readline().strip() for _ in range(2)]
    # diag_manager reads the title and base date from lines 1 and 2 exactly
    if not head[0] or head[0].startswith("#"):
        bad += 1
        report(1, "line 1 must be the title (a quoted string)")
    if not re.match(r"^\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\b", head[1]):
        bad += 1
        report(2, "line 2 must be the base date: year month day hour minute second")

    files = {}
    outputs = {}
    for ln, toks in diag_table_lines(path):
        if ln <= 2:
            continue
        if len(toks) >= 6 and re.match(r"^-?\d+$", toks[1]):
            fname = re.sub(r"\.nc$", "", toks[0])
            if toks[2].lower() not in TIME_UNITS or toks[4].lower() not in TIME_UNITS:
                bad += 1
                report(ln, "file %s: units must be one of %s" % (fname, ", ".join(sorted(TIME_UNITS))))
            if "time" not in toks[5].lower():
                bad += 1
                report(ln, "file %s: the time axis name must contain 'time'" % fname)
            files[fname] = ln
            continue
        if len(toks) != 8:
            bad += 1
            report(ln, "not a file line (6 or more columns) or a field line (8 columns)")
            continue
        mod, fld, oname, fname, _, red, region, pack = toks
        fname = re.sub(r"\.nc$", "", fname)
        redl = red.lower()
        if fname.lower() != "null" and fname not in files:
            bad += 1
            report(ln, "%s/%s: file %s is not defined in this diag_table" % (mod, fld, fname))
        if redl not in REDUCTIONS and not re.match(r"^(pow|diurnal)\d+$", redl):
            bad += 1
            report(ln, "%s/%s: unknown reduction %s" % (mod, fld, red))
        if not re.match(r"^[1-8]$", pack):
            bad += 1
            report(ln, "%s/%s: packing must be 1 (double), 2 (float), 4 or 8" % (mod, fld))
        if region.lower() != "none" and len(region.split()) != 6:
            bad += 1
            report(ln, "%s/%s: region must be none or 6 numbers" % (mod, fld))
        eff = oname
        suffix = REDUCTION_SUFFIX.get(redl)
        if suffix and len(oname) >= 3 and oname[-3:].lower() != suffix:
            eff = oname + "_" + suffix
        if (fname, eff) in outputs:
            bad += 1
            report(ln, "%s/%s: output name %s already used in %s at line %d" % (
                mod, fld, eff, fname, outputs[(fname, eff)]))
        else:
            outputs[(fname, eff)] = ln

        hit = [r for pm, pfld, r in pats if pm.match(mod) and pfld.match(fld)]
        if not hit:
            bad += 1
            near = sorted({r.module for r in rows if r.field.lower() == fld.lower()})
            hint = (" (the field exists in module %s)" % ", ".join(near)) if near else ""
            report(ln, "%s/%s is not registered by any source file%s" % (mod, fld, hint))
            continue
        notes = []
        if not any(h.built for h in hit):
            notes.append("only in files that are not compiled")
        alts = []
        for h in hit:
            for a in h.alts:
                if a not in alts:
                    alts.append(a)
        label = cnd.alts_label(alts, md=False)
        if nml_values is not None:
            ok = cnd.alts_eval(alts, value_of)
            if ok is False:
                bad += 1
                report(ln, "%s/%s is not registered with these namelist settings; it needs %s"
                       % (mod, fld, label))
                continue
            if ok is None:
                notes.append("registered only when %s" % label)
        elif label:
            notes.append("registered only when %s" % label)
        if all(h.kind == "static" for h in hit) and redl not in (".false.", "none", "point"):
            notes.append("static field: use .false.")
        if any(h.send_status == "never-sent" for h in hit):
            notes.append("registered but never sent")
        if notes:
            report(ln, "%s/%s ok (%s)" % (mod, fld, "; ".join(notes)))
    return 1 if bad else 0
