"""Markdown for the generated parts of the documentation: the diagnostic
field tables (docs/Diagnostics.md) and the namelist tables
(docs/Parameters.md)."""

import os
import re
from collections import OrderedDict

from . import conditions as cnd
from .conditions import alt_label, md_code
from .diagnostics import merge_fields

SOURCE_URL = "https://github.com/Eddy-Stanford/MiMA/blob/main/"


def md_escape(s):
    s = str(s).replace("|", "\\|").replace("\n", " ")
    return s.replace("<", "&lt;").replace(">", "&gt;")


def source_link(path):
    return "[%s](%s%s)" % (os.path.basename(path), SOURCE_URL, path)


def anchor(module):
    return "module-" + re.sub(r"[^a-z0-9_-]", "", module.lower())


# --------------------------------------------------------------------------
# Diagnostics
# --------------------------------------------------------------------------

def module_availability(rows):
    """Factor the conditions of a module's rows.  Returns (common clauses,
    {id(row): alternatives without the common clauses}, selector), where
    selector is (var, groups) when every alternative of every row selects on
    the value of the same character variable (e.g. radiation_scheme)."""
    all_alts = [a for r in rows for a in r.alts]
    common = []
    if all_alts and all(all_alts):
        for c in all_alts[0]:
            if all(c in a for a in all_alts) and c not in common:
                common.append(c)
    cells = {id(r): [tuple(c for c in a if c not in common) for a in r.alts] for r in rows}
    keys, values = set(), set()
    for r in rows:
        for a in cells[id(r)]:
            sel, _ = cnd.split_selector(a)
            if sel is None:
                return common, cells, None
            keys.add((sel[0], sel[2].groups))
            values.add(sel[1])
    if len(keys) == 1 and len(values) > 1:
        return common, cells, keys.pop()
    return common, cells, None


def selector_values(alts):
    return sorted({cnd.split_selector(a)[0][1] for a in alts})


def selector_cell(alts):
    """'gray, rrtm' or, with further conditions, 'rrtm if `x` (group)'."""
    groups = OrderedDict()
    for a in alts:
        sel, extra = cnd.split_selector(a)
        groups.setdefault(extra, [])
        if sel[1] not in groups[extra]:
            groups[extra].append(sel[1])
    if () in groups:
        # values without further conditions make the conditional ones redundant
        base = set(groups[()])
        for k in list(groups):
            if k:
                groups[k] = [v for v in groups[k] if v not in base]
                if not groups[k]:
                    del groups[k]
    parts = []
    for extra in sorted(groups, key=lambda e: (len(e), alt_label(e))):
        parts.append(", ".join(sorted(groups[extra])) + (" if " + alt_label(extra) if extra else ""))
    return "<br>".join(parts)


def render_diagnostics(rows):
    """The field tables: a summary of the modules, one table per module and
    the fields whose metadata differ between registrations."""
    rows = [r for r in rows if r.built]
    mods = OrderedDict()
    for r in rows:
        mods.setdefault(r.module, []).append(r)
    info = {m: module_availability(rs) for m, rs in mods.items()}
    out = []
    w = out.append

    w("| module | fields | registered when | source |")
    w("|---|---:|---|---|")
    for m, rs in mods.items():
        common, cells, sel = info[m]
        when = alt_label(common) if common else ""
        if sel is not None:
            vals = sorted({v for r in rs for v in selector_values(cells[id(r)])})
            when = (when + " and " if when else "") + "%s (%s): %s" % (
                md_code(sel[0]), ", ".join(re.sub(r"_nml$", "", g) for v, g in sel[1]),
                ", ".join(vals))
        srcs = sorted({r.file for r in rs}, key=os.path.basename)
        w("| [`%s`](#%s) | %d | %s | %s |" % (
            m, anchor(m), len(merge_fields(rs)), when, ", ".join(source_link(s) for s in srcs)))
    w("")

    differs = []
    for m, rs in mods.items():
        common, cells, sel = info[m]
        w('<a id="%s"></a>' % anchor(m))
        w("")
        w("### `%s`" % m)
        w("")
        if common:
            w("Registered only when %s." % alt_label(common))
            w("")
        if sel is not None:
            w("Which fields exist depends on %s in `&%s`: **available with** lists the "
              "values for which the field is registered." % (
                  md_code(sel[0]), ", ".join(g for v, g in sel[1])))
            w("")
        w("| field | dims | units | long_name | available with | source |")
        w("|---|---|---|---|---|---|")
        for frs in merge_fields(rs).values():
            alts = []
            for r in frs:
                for a in cells[id(r)]:
                    if a not in alts:
                        alts.append(a)
            if any(len(a) == 0 for a in alts):
                avail = ""
            elif sel is not None:
                avail = selector_cell(alts)
            else:
                avail = "<br>or ".join(alt_label(a) for a in alts)

            def label(r):
                if sel is not None:
                    return ", ".join(selector_values(cells[id(r)]))
                return os.path.basename(r.file)

            cols, diff_cols = [], []
            for name, get in (("dims", lambda r: r.dims + ("; static" if r.kind == "static" else "")),
                              ("units", lambda r: r.units),
                              ("long_name", lambda r: r.long_name)):
                vals = []
                for r in frs:
                    if get(r) not in vals:
                        vals.append(get(r))
                if len(vals) == 1:
                    cols.append(md_escape(vals[0]))
                else:
                    diff_cols.append(name)
                    cols.append("**differs**<br>" + "<br>".join(
                        "%s: %s" % (label(r), md_escape(get(r))) for r in frs))
            if diff_cols:
                differs.append((m, frs[0].field, diff_cols, frs))
            locs = []
            for r in frs:
                loc = source_link(r.file)
                if r.send_status == "never-sent":
                    loc += " *(never sent)*"
                if loc not in locs:
                    locs.append(loc)
            fname = frs[0].field
            if len({r.field for r in frs}) > 1:
                fname = " / ".join(sorted({r.field for r in frs}))
            w("| %s | %s | %s | %s | %s | %s |" % (
                md_code(fname), cols[0], cols[1], cols[2], avail, "<br>".join(locs)))
        w("")
    if differs:
        w("## Fields with inconsistent metadata")
        w("")
        w("These fields are registered in more than one place with different metadata. "
          "They should be made consistent in the source.")
        w("")
        w("| module | field | differs in | registered in |")
        w("|---|---|---|---|")
        for m, f, dcols, frs in differs:
            w("| `%s` | %s | %s | %s |" % (m, md_code(f), ", ".join(dcols), ", ".join(
                sorted({source_link(r.file) for r in frs}))))
        w("")
    return "\n".join(out).rstrip("\n")


# --------------------------------------------------------------------------
# Namelists
# --------------------------------------------------------------------------

def render_namelist(variables):
    """The table of a namelist group.  Variables declared on the same line
    with the same documentation share a row."""
    rows = OrderedDict()
    for v in variables:
        key = (v.file, v.line, v.doc) if v.doc else id(v)
        rows.setdefault(key, []).append(v)
    out = ["| Variable | Type | Default | Description |", "|---|---|---|---|"]
    for vs in rows.values():
        types = []
        for v in vs:
            if v.type not in types:
                types.append(v.type)
        defaults = ", ".join(md_code(v.default) if v.default else "unset" for v in vs)
        out.append("| %s | %s | %s | %s |" % (
            ", ".join("`%s`" % v.spelling for v in vs), ", ".join(types), defaults,
            vs[0].doc.replace("|", "\\|")))
    return "\n".join(out)
