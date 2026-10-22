"""Namelist groups and their variables (type, default, documentation) from
the source, and a reader for input.nml files."""

import io
import re
from dataclasses import dataclass

from .fortran import doc_comment, mask_strings, norm_ws


@dataclass
class NamelistVar:
    group: str
    name: str          # lower case
    spelling: str      # as written in the namelist statement
    type: str          # e.g. 'real', 'character(len=64), dimension(10)'
    default: str       # Fortran text of the default ('' if unset)
    value: object      # Python value of a scalar default, or None
    doc: str           # the '!!' documentation
    file: str
    line: int
    module: str = ""   # module (or program) whose namelist statement it is


def fortran_value(text):
    """Python value of a scalar namelist/initialiser constant, or None."""
    t = norm_ws(text or "").rstrip(",").strip()
    if not t:
        return None
    tl = t.lower()
    if tl in (".true.", ".t.", "t", "true"):
        return True
    if tl in (".false.", ".f.", "f", "false"):
        return False
    m = re.match(r"^'((?:[^']|'')*)'$|^\"((?:[^\"]|\"\")*)\"$", t)
    if m:
        return m.group(1).replace("''", "'") if m.group(1) is not None else m.group(2).replace('""', '"')
    try:
        return float(re.sub(r"_\w+$", "", tl).replace("d", "e"))
    except ValueError:
        return None


def _default_text(sym, scope, evaluator, depth=0):
    """The default of a variable as Fortran text, with named parameters
    replaced by their values."""
    if sym.expr is None:
        return ""
    expr = norm_ws(sym.expr)
    if depth < 5 and re.match(r"^[A-Za-z_]\w*$", expr) and expr.lower() not in (".true.", ".false."):
        psym, pscope = evaluator.lookup(expr.lower(), scope)
        if psym is not None and psym.is_param and psym.expr is not None:
            return _default_text(psym, pscope, evaluator, depth + 1)
    return expr


def _declaration(name, scope, modtab):
    """(Symbol, scope) of a variable: in the scope chain, else in a module
    that the scope uses."""
    for s in scope.chain():
        if name in s.symbols:
            return s.symbols[name], s
    for s in scope.chain():
        for u in s.uses:
            only = s.use_only.get(u)
            ms = modtab.get(u)
            if ms is not None and name in ms.symbols and (only is None or name in only):
                return ms.symbols[name], ms
    return None, None


def _type_label(sym):
    t = re.sub(r"\s*([(),=*])\s*", r"\1", sym.spec.lower()).replace(",", ", ")
    if sym.dims and "dimension" not in t:
        t += ", dimension(%s)" % norm_ws(sym.dims)
    return t


def namelist_reference(parsed, evaluator):
    """{group: [NamelistVar in namelist order]} for every namelist group."""
    out = {}
    for pf in parsed:
        for nl in pf.namelist_stmts:
            rows = out.setdefault(nl.group, [])
            for v, spelling in zip(nl.vars, nl.spellings):
                sym, sscope = _declaration(v, nl.scope, evaluator.modtab)
                owner = nl.scope.module_scope()
                owner = owner.name if owner is not None and owner.kind != "file" else ""
                if sym is None:
                    rows.append(NamelistVar(nl.group, v, spelling, "", "", None, "", pf.relpath,
                                            nl.line, owner))
                    continue
                default = _default_text(sym, sscope, evaluator)
                rows.append(NamelistVar(
                    nl.group, v, spelling, _type_label(sym), default,
                    fortran_value(default) if not sym.dims else None,
                    doc_comment(pf.lines, sym.line), pf.relpath.replace("\\", "/"), sym.line, owner))
    return out


def namelist_defaults(reference):
    """{(group, var): Python value of the default} for evaluating conditions."""
    return {(g, v.name): v.value for g, vs in reference.items() for v in vs}


_SETTING = re.compile(r"([A-Za-z_]\w*)\s*(\([^)]*\))?\s*=")


def _namelist_groups(path):
    """(group, line, body) for each &group ... / of an input.nml, in order, with
    the comments removed. Line i of body is line `line + i` of the file."""
    with io.open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()
    out = []
    cur = None
    first = 0
    buf = []
    for n, ln in enumerate(lines, 1):
        m = mask_strings(ln.rstrip("\n"))
        k = m.find("!")
        if k >= 0:
            ln, m = ln[:k], m[:k]
        ln = ln.rstrip("\n")
        if cur is None:
            mg = re.match(r"^\s*&(\w+)", ln)
            if not mg:
                continue
            cur = mg.group(1).lower()
            if cur == "end":
                cur = None
                continue
            ln, m = ln[mg.end():], m[mg.end():]
            first, buf = n, []
        e = m.find("/")
        me = re.search(r"&end\b", m, re.I)
        if me is not None and (e < 0 or me.start() < e):
            e = me.start()
        if e >= 0:
            buf.append(ln[:e])
            out.append((cur, first, "\n".join(buf)))
            cur = None
        else:
            buf.append(ln)
    if cur is not None:
        out.append((cur, first, "\n".join(buf)))
    return out


def _settings(body):
    """(name, subscript, value text, offset) of each setting in a group body."""
    keys = list(_SETTING.finditer(mask_strings(body)))
    for i, mk in enumerate(keys):
        end = keys[i + 1].start() if i + 1 < len(keys) else len(body)
        raw = body[mk.end():end].strip().rstrip(",").strip()
        yield mk.group(1).lower(), mk.group(2), raw, mk.start()


def read_namelist_file(path):
    """Minimal input.nml reader: {group: {var: value}} for scalar settings.
    Array settings (x = 1, 2, 3 or x(2) = ...) are skipped."""
    groups = {}
    for group, _, body in _namelist_groups(path):
        vals = groups.setdefault(group, {})
        for name, sub, raw, _ in _settings(body):
            if sub or "," in mask_strings(raw):
                continue
            vals[name] = fortran_value(raw)
    return groups


def read_namelist_keys(path):
    """[(group, line, [(var, line), ...])] for each group of an input.nml, with
    every variable it sets, array settings included."""
    out = []
    for group, first, body in _namelist_groups(path):
        out.append((group, first, [(name, first + body[:off].count("\n"))
                                   for name, _, _, off in _settings(body)]))
    return out


# Groups that FMS 2026.02 reads (except its tests), which are not in MiMA's
# source. A group in an input.nml that neither reads is ignored at run time.
FMS_GROUPS = frozenset("""
    amip_interp_nml answers_nml astronomy_nml column_diagnostics_nml data_override_nml
    diag_integral_nml diag_manager_nml drifters_comm_nml ensemble_nml field_manager_nml
    fms2_io_nml fms_affinity_nml fms_nml gaussian_topog_nml get_cal_time_nml
    horiz_interp_spherical_nml interpolator_nml monin_obukhov_nml mpp_domains_nml mpp_nml
    offloading_io_nml sat_vapor_pres_nml time_interp_external_nml time_interp_nml
    topography_nml xgrid_nml""".split())

# Groups that make the model stop, with the reason.
FATAL_GROUPS = {
    "simple_surface_nml": "the v1 group simple_surface_nml stops the model; convert it to the "
                          "v2 surface namelists (see docs/Migration_v2.md)",
}


def check_namelist_file(path, reference, out):
    """Report the settings of an input.nml that the model would stop on: a
    variable that is not in its MiMA group, or a group in FATAL_GROUPS. A group
    that neither MiMA nor FMS reads gets a warning. Returns 1 on an error."""
    bad = 0
    for group, gline, keys in read_namelist_keys(path):
        if group in FATAL_GROUPS:
            bad += 1
            out.write("%s:%d: %s\n" % (path, gline, FATAL_GROUPS[group]))
        elif group in reference:
            known = {v.name for v in reference[group]}
            for name, line in keys:
                if name not in known:
                    bad += 1
                    out.write("%s:%d: %s is not a variable of %s; the model would stop with "
                              "\"Unknown namelist, or mistyped namelist variable\"\n"
                              % (path, line, name, group))
        elif group not in FMS_GROUPS:
            out.write("%s:%d: warning: no MiMA or FMS module reads %s, so its settings are "
                      "ignored\n" % (path, gline, group))
    return 1 if bad else 0
