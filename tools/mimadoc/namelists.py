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


def read_namelist_file(path):
    """Minimal input.nml reader: {group: {var: value}} for scalar settings.
    Array settings (x = 1, 2, 3 or x(2) = ...) are skipped."""
    with io.open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()
    groups = {}
    cur = None
    buf = []

    def flush():
        body = "\n".join(buf)
        masked = mask_strings(body)
        keys = list(re.finditer(r"([A-Za-z_]\w*)\s*(\([^)]*\))?\s*=", masked))
        vals = groups.setdefault(cur, {})
        for i, mk in enumerate(keys):
            end = keys[i + 1].start() if i + 1 < len(keys) else len(body)
            raw = body[mk.end():end].strip().rstrip(",").strip()
            if mk.group(2) or "," in mask_strings(raw):
                continue
            vals[mk.group(1).lower()] = fortran_value(raw)

    for ln in lines:
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
            buf = []
        e = m.find("/")
        me = re.search(r"&end\b", m, re.I)
        if me is not None and (e < 0 or me.start() < e):
            e = me.start()
        if e >= 0:
            buf.append(ln[:e])
            flush()
            cur = None
        else:
            buf.append(ln)
    if cur is not None:
        flush()
    return groups
