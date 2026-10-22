"""Evaluate the string and integer expressions passed to register_diag_field
(module and field names, long names, units, axes) from the declarations and
assignments in the source."""

import re

from .fortran import WriteValue, array_items, match_paren, mask_strings, norm_ws, \
    paren_group, split_top, unquote

# MiMA passes the atmosphere axes as an array: axes(1:3) = lon, lat, pfull.
AXIS_INDEX_NAMES = {1: "lon", 2: "lat", 3: "pfull", 4: "phalf"}
# Intrinsics that return (a transformation of) their string argument.
STRING_FUNCS = ("trim", "adjustl", "adjustr", "lowercase", "uppercase", "lcase", "ucase")


class Evaluation(object):
    """What an evaluation found out besides the value."""

    def __init__(self):
        self.expand = {}      # loop variable -> (array name, number of elements)
        self.notes = []


class Evaluator(object):
    def __init__(self, modtab, max_depth=8):
        self.modtab = modtab       # module name -> module Scope
        self.max_depth = max_depth

    # -- symbol lookup through the scope chain and USE'd modules
    def lookup(self, name, scope):
        for s in scope.chain():
            if name in s.symbols:
                return s.symbols[name], s
        seen = set()
        for s in scope.chain():
            for u in s.uses:
                if u in seen:
                    continue
                seen.add(u)
                ms = self.modtab.get(u)
                if ms is not None and name in ms.symbols and ms.symbols[name].is_param:
                    return ms.symbols[name], ms
        return None, None

    def last_assignment(self, name, scope, stmt_idx):
        """The last whole-variable assignment to name before stmt_idx in the
        procedure (or its host): ((index, rhs, subscript), scope)."""
        for s in scope.chain():
            prev = [a for a in s.assigns.get(name, ()) if a[0] < stmt_idx and a[2] is None]
            if prev:
                return prev[-1], s
            if s.kind in ("subroutine", "function"):
                break
        return None, None

    def tracer_role(self, name, scope):
        for s in scope.chain():
            if name in s.tracer_vars:
                return s.tracer_vars[name]
        return None

    # -- strings
    def eval_str(self, expr, scope, stmt_idx, env, info, depth=0):
        """(string, fully resolved).  Unresolved parts are rendered as {expr}
        (or <tracer> for tracer names).  env maps loop variables to values."""
        expr = expr.strip()
        if depth > self.max_depth:
            return "{%s}" % norm_ws(expr), False
        out = []
        ok = True
        for t in split_top(expr, "//"):
            s, r = self._term(t.strip(), scope, stmt_idx, env, info, depth)
            out.append(s)
            ok = ok and r
        return "".join(out), ok

    def _term(self, t, scope, stmt_idx, env, info, depth):
        if not t:
            return "", True
        if t.startswith("(") and not t.startswith("(/"):
            if match_paren(mask_strings(t), 0) == len(t) - 1:
                return self.eval_str(t[1:-1], scope, stmt_idx, env, info, depth + 1)
        mq = re.match(r"^(?:\w+_)?('([^']|'')*'|\"([^\"]|\"\")*\")$", t, re.S)
        if mq:
            return unquote(mq.group(1)), True
        mf = re.match(r"^(%s)\s*\(" % "|".join(STRING_FUNCS), t, re.I)
        if mf:
            inner, e = paren_group(t, mf.end() - 1)
            if e == len(t) - 1:
                s, r = self.eval_str(inner, scope, stmt_idx, env, info, depth + 1)
                if r:
                    s = _apply_string_func(mf.group(1).lower(), s)
                return s, r
        mi = re.match(r"^(\w+)$", t)
        if mi:
            name = mi.group(1).lower()
            if isinstance(env.get(name), str):
                return env[name], True
            return self._name(name, scope, stmt_idx, env, info, depth)
        ma = re.match(r"^(\w+)\s*\(\s*([^()]+?)\s*\)$", t)
        if ma:
            name = ma.group(1).lower()
            idx = ma.group(2).strip().lower()
            elems = self.literal_array(name, scope, stmt_idx, depth)
            if elems is not None:
                iv = self.eval_int(idx, scope, env)
                if iv is not None and 1 <= iv <= len(elems):
                    return elems[iv - 1], True
                if re.match(r"^\w+$", idx):
                    info.expand[idx] = (name, len(elems))
        return "{%s}" % norm_ws(t), False

    def _name(self, name, scope, stmt_idx, env, info, depth):
        sym, sscope = self.lookup(name, scope)
        # a local assignment wins for non-parameters
        if sym is None or not sym.is_param:
            asg, ascope = self.last_assignment(name, scope, stmt_idx)
            if asg is not None and isinstance(asg[1], WriteValue):
                return "{%s}" % asg[1], False
            if asg is not None:
                s, r = self.eval_str(asg[1], ascope, asg[0], env, info, depth + 1)
                return (self._fit(s, sym, info, name), True) if r else (s, False)
        role = self.tracer_role(name, scope)
        if role:
            return "<%s>" % role, False
        if sym is not None and sym.expr is not None and array_items(sym.expr) is None:
            s, r = self.eval_str(sym.expr, sscope, 0, env, info, depth + 1)
            if r:
                return self._fit(s, sym, info, name), True
        return "{%s}" % name, False

    def _fit(self, s, sym, info, name):
        """Truncate s to the declared length of the character variable."""
        if sym is None or sym.charlen is None:
            return s
        try:
            n = int(sym.charlen)
        except ValueError:
            return s
        if len(s.rstrip()) > n:
            info.notes.append("TRUNCATED: %s is character(len=%d) but value '%s' is longer"
                              % (name, n, s))
            return s[:n]
        return s

    def literal_array(self, name, scope, stmt_idx, depth=0):
        """The evaluated elements of a character array set from an array
        constructor, or None."""
        sym, sscope = self.lookup(name, scope)
        items = array_items(sym.expr) if sym is not None and sym.expr else None
        if items is None:
            asg, ascope = self.last_assignment(name, scope, stmt_idx)
            if asg is not None:
                items = array_items(asg[1])
                sscope = ascope
        if items is None:
            return None
        vals = []
        for item in items:
            s, r = self.eval_str(item, sscope, 0, {}, Evaluation(), depth + 1)
            if not r:
                return None
            vals.append(self._fit(s, sym, Evaluation(), name))
        return vals

    # -- integers
    def eval_int(self, expr, scope, env):
        expr = expr.strip().lower()
        if re.match(r"^[+-]?\d+$", expr):
            return int(expr)
        if isinstance(env.get(expr), int):
            return env[expr]
        ms = re.match(r"^size\s*\(\s*(\w+)\s*(\(\s*:\s*\))?\s*\)$", expr)
        if ms:
            arr = self.literal_array(ms.group(1), scope, 10 ** 9)
            return len(arr) if arr is not None else None
        if re.match(r"^\w+$", expr):
            sym, sscope = self.lookup(expr, scope)
            if sym is not None and sym.is_param and sym.expr:
                return self.eval_int(sym.expr, sscope, env)
        mo = re.match(r"^(\w+)\s*([+-])\s*(\w+)$", expr)
        if mo:
            a = self.eval_int(mo.group(1), scope, env)
            b = self.eval_int(mo.group(3), scope, env)
            if a is not None and b is not None:
                return a + b if mo.group(2) == "+" else a - b
        return None

    # -- axes
    def axes(self, axes, scope, stmt_idx):
        """(axis names such as 'lon, lat, pfull' or '', number of dimensions
        or None).  id_<name> variables are named axes."""
        a = norm_ws(axes or "")
        if not a:
            return "", None
        if a == "(scalar)":
            return "scalar", 0

        def indexed(items):
            try:
                return ", ".join(AXIS_INDEX_NAMES[int(i)] for i in items)
            except (KeyError, ValueError):
                return ""

        def named(items):
            out = []
            for it in items:
                it = it.strip()
                m = re.match(r"^\w+\s*\(\s*(\d+)\s*\)$", it)
                if m and int(m.group(1)) in AXIS_INDEX_NAMES:
                    out.append(AXIS_INDEX_NAMES[int(m.group(1))])
                    continue
                m = re.match(r"^id_(\w+)$", it, re.I)
                if not m:
                    return ""
                out.append(m.group(1).lower())
            return ", ".join(out)

        def value_of(name):
            sym, _ = self.lookup(name, scope)
            items = array_items(sym.expr) if sym is not None and sym.expr else None
            if items is None:
                asg, _ = self.last_assignment(name, scope, stmt_idx)
                if asg is not None:
                    items = array_items(asg[1])
            return items, sym

        def declared_size(sym):
            d = (sym.dims or "").strip() if sym is not None else ""
            if re.match(r"^\d+$", d):
                return int(d)
            m = re.match(r"^(\d+)\s*:\s*(\d+)$", d)
            return int(m.group(2)) - int(m.group(1)) + 1 if m else None

        m = re.match(r"^\w+\s*\(\s*(\d+)\s*:\s*(\d+)\s*\)$", a)
        if m:
            lo, hi = int(m.group(1)), int(m.group(2))
            return indexed(range(lo, hi + 1)), hi - lo + 1
        items = array_items(a)
        if items is not None:
            return named(items), len(items)
        m = re.match(r"^\w+\s*\(\s*(\w+)\s*\)$", a)          # axes(half)
        if m:
            items, sym = value_of(m.group(1).lower())
            if sym is not None and sym.expr and array_items(sym.expr) is not None:
                n = len(array_items(sym.expr))
            else:
                n = declared_size(sym) if sym is not None and re.match(
                    r"^\d+$", (sym.dims or "").strip()) else None
            return (indexed(items) if items else ""), n
        m = re.match(r"^(\w+)$", a)
        if m:
            items, sym = value_of(m.group(1).lower())
            n = declared_size(sym)
            if items:
                return named(items), n
            if sym is not None and re.match(r"^\d+$", (sym.dims or "").strip()):
                return indexed(range(1, n + 1)), n
            return "", n
        return "", None


def _apply_string_func(fn, s):
    if fn == "trim":
        return s.rstrip()
    if fn == "adjustl":
        return s.lstrip() + " " * (len(s) - len(s.lstrip()))
    if fn in ("lowercase", "lcase"):
        return s.lower()
    if fn in ("uppercase", "ucase"):
        return s.upper()
    return s
