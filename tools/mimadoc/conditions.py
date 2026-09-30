"""Conditions under which a diagnostic is registered.

The source parser records conditions as Fortran text.  Here they are parsed
once into a small tree, which is printed in a readable form
    trim(radiation_scheme) == ('rrtm')   ->  radiation_scheme = 'rrtm'
    (.not. (do_bm))                      ->  not do_bm
    x == ('a', 'b')    (select case)     ->  x in ('a', 'b')
and evaluated (three-valued) against the settings of an input.nml.

A Clause is one top-level '.and.' term, with the namelist group of each
namelist variable in it.  An alternative is a tuple of clauses that must all
hold (one call path to the registration); a field is registered when any of
its alternatives holds."""

import re
from collections import namedtuple

from .fortran import norm_ws


class Clause(namedtuple("Clause", "text groups node")):
    """text: readable form; groups: ((var, namelist group), ...); node: the
    parsed tree, or None if the text could not be parsed.  Clauses compare
    by text and groups only."""
    __slots__ = ()

    def _key(self):
        return (self.text, self.groups)

    def __eq__(self, other):
        return isinstance(other, Clause) and self._key() == other._key()

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(self._key())


# --------------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------------

_DOP_RE = r"\.(?:and|or|not|eqv|neqv|eq|ne|lt|le|gt|ge|true|false)\."
TOKEN_RE = re.compile(
    r"\s*(?:(?P<str>'(?:[^']|'')*'|\"(?:[^\"]|\"\")*\")"
    r"|(?P<dop>" + _DOP_RE + r")"
    r"|(?P<num>\d+(?:\.(?!(?:and|or|not|eqv|neqv|eq|ne|lt|le|gt|ge)\.)\d*)?(?:[eEdD][+-]?\d+)?)"
    r"|(?P<name>[A-Za-z_][\w%]*)"
    r"|(?P<op>==|/=|>=|<=|=|<|>|\(|\)|,))", re.I)

CMP_OPS = {"==": "=", "=": "=", ".eq.": "=", "/=": "/=", ".ne.": "/=",
           ".eqv.": "=", ".neqv.": "/=", "<": "<", ".lt.": "<", "<=": "<=",
           ".le.": "<=", ">": ">", ".gt.": ">", ">=": ">=", ".ge.": ">="}
CMP_NEG = {"=": "/=", "/=": "=", "<": ">=", ">=": "<", ">": "<=", "<=": ">"}
TRANSPARENT_FUNCS = {"trim", "adjustl"}
# a clause that selects on the value of a character variable
SELECTOR_RE = re.compile(r"^(\w+) = '([^']*)'$")


class ParseError(Exception):
    pass


def _tokens(text):
    toks = []
    pos = 0
    text = text.rstrip()
    while pos < len(text):
        m = TOKEN_RE.match(text, pos)
        if not m or m.end() == pos:
            raise ParseError(text)
        pos = m.end()
        kind = m.lastgroup
        val = m.group(kind)
        if kind == "dop":
            val = val.lower()
            if val in (".and.", ".or.", ".not."):
                kind, val = "kw", val.strip(".")
            elif val in (".true.", ".false."):
                kind = "bool"
            else:
                kind = "op"
        elif kind == "name" and val.lower() in ("and", "or", "not", "in"):
            kind, val = "kw", val.lower()
        toks.append((kind, val))
    return toks


class _Parser(object):
    """Recursive descent parser.  Nodes are tuples:
    ('or'|'and', [nodes]), ('not', node), ('cmp', op, left, right),
    ('in', left, [values]), ('lit', value, text), ('var', name),
    ('call', name, argument text), ('tuple', [nodes])."""

    def __init__(self, text):
        self.text = norm_ws(text)
        self.toks = _tokens(self.text)
        self.i = 0

    def peek(self):
        return self.toks[self.i] if self.i < len(self.toks) else (None, None)

    def take(self):
        t = self.peek()
        self.i += 1
        return t

    def expect(self, val):
        if self.take()[1] != val:
            raise ParseError(self.text)

    def parse(self):
        n = self.p_or()
        if self.i != len(self.toks):
            raise ParseError(self.text)
        return n

    def p_or(self):
        items = [self.p_and()]
        while self.peek() == ("kw", "or"):
            self.take()
            items.append(self.p_and())
        return items[0] if len(items) == 1 else ("or", items)

    def p_and(self):
        items = [self.p_not()]
        while self.peek() == ("kw", "and"):
            self.take()
            items.append(self.p_not())
        flat = []
        for it in items:
            flat.extend(it[1] if it[0] == "and" else [it])
        return flat[0] if len(flat) == 1 else ("and", flat)

    def p_not(self):
        if self.peek() == ("kw", "not"):
            self.take()
            return ("not", self.p_not())
        return self.p_cmp()

    def p_cmp(self):
        left = self.p_primary()
        k, v = self.peek()
        if k == "kw" and v == "in":
            self.take()
            right = self.p_primary()
            return ("in", left, right[1] if right[0] == "tuple" else [right])
        if k == "op" and v in CMP_OPS:
            self.take()
            right = self.p_primary()
            if right[0] == "tuple":
                return ("in", left, right[1])
            return ("cmp", CMP_OPS[v], left, right)
        if left[0] == "tuple":
            raise ParseError(self.text)
        return left

    def p_primary(self):
        k, v = self.take()
        if k is None:
            raise ParseError(self.text)
        if k == "op" and v == "(":
            first = self.p_or()
            if self.peek() == ("op", ","):
                items = [first]
                while self.peek() == ("op", ","):
                    self.take()
                    items.append(self.p_or())
                self.expect(")")
                return ("tuple", items)
            self.expect(")")
            return first
        if k == "str":
            s = v[1:-1].replace(v[0] * 2, v[0])
            return ("lit", s, "'%s'" % s.replace("'", "''"))
        if k == "num":
            try:
                val = float(v.lower().replace("d", "e"))
            except ValueError:
                val = None
            return ("lit", val, v)
        if k == "bool":
            return ("lit", v == ".true.", v)
        if k == "name":
            name = v.lower()
            if self.peek() != ("op", "("):
                return ("var", name)
            # function call or array element: keep the argument text
            start = self.i
            depth = 0
            while True:
                tk, tv = self.take()
                if tk is None:
                    raise ParseError(self.text)
                if tk == "op" and tv == "(":
                    depth += 1
                elif tk == "op" and tv == ")":
                    depth -= 1
                    if depth == 0:
                        break
            inner = self.toks[start + 1:self.i - 1]
            if name in TRANSPARENT_FUNCS and len(inner) == 1 and inner[0][0] == "name":
                return ("var", inner[0][1].lower())
            return ("call", name, _join_tokens(inner))
        raise ParseError(self.text)


def _join_tokens(toks):
    out = ""
    for k, v in toks:
        if v == ",":
            out += ", "
        elif k == "kw" or (k == "op" and v not in ("(", ")")):
            out += " %s " % v
        else:
            out += v
    return norm_ws(out)


def parse(text):
    """The tree of a condition, or None if it cannot be parsed."""
    try:
        return _Parser(text).parse()
    except ParseError:
        return None


# --------------------------------------------------------------------------
# Printing, variables, evaluation
# --------------------------------------------------------------------------

def show(n, prec=0):
    k = n[0]
    if k == "or":
        s = " or ".join(show(c, 1) for c in n[1])
        return "(%s)" % s if prec > 1 else s
    if k == "and":
        s = " and ".join(show(c, 2) for c in n[1])
        return "(%s)" % s if prec > 2 else s
    if k == "not":
        inner = n[1]
        if inner[0] == "not":
            return show(inner[1], prec)
        if inner[0] == "cmp":
            return show(("cmp", CMP_NEG[inner[1]], inner[2], inner[3]), prec)
        if inner[0] == "in":
            return "%s not in (%s)" % (show(inner[1], 5), ", ".join(show(v, 5) for v in inner[2]))
        return "not " + show(inner, 3)
    if k == "cmp":
        return "%s %s %s" % (show(n[2], 5), n[1], show(n[3], 5))
    if k == "in":
        if len(n[2]) == 1:
            return "%s = %s" % (show(n[1], 5), show(n[2][0], 5))
        return "%s in (%s)" % (show(n[1], 5), ", ".join(show(v, 5) for v in n[2]))
    if k == "lit":
        return n[2]
    if k == "var":
        return n[1]
    if k == "call":
        return "%s(%s)" % (n[1], n[2])
    return "?"


def variables(n, out=None):
    """Names of the variables in a tree, in order of appearance."""
    if out is None:
        out = []
    k = n[0]
    if k == "var":
        if n[1] not in out:
            out.append(n[1])
    elif k in ("or", "and", "tuple"):
        for c in n[1]:
            variables(c, out)
    elif k == "not":
        variables(n[1], out)
    elif k == "cmp":
        variables(n[2], out)
        variables(n[3], out)
    elif k == "in":
        variables(n[1], out)
        for c in n[2]:
            variables(c, out)
    elif k == "call":
        for t in re.findall(r"[a-z_]\w*", n[2].lower()):
            if t not in out:
                out.append(t)
    return out


def evaluate(n, value_of):
    """Three-valued evaluation: True, False or None (unknown).  value_of(name)
    returns a Python bool/str/float, or None when the value is unknown."""
    k = n[0]
    if k == "lit":
        return n[1]
    if k == "var":
        return value_of(n[1])
    if k in ("call", "tuple"):
        return None
    if k == "not":
        v = evaluate(n[1], value_of)
        return (not v) if isinstance(v, bool) else None
    if k in ("and", "or"):
        vals = [evaluate(c, value_of) for c in n[1]]
        vals = [v if isinstance(v, bool) else None for v in vals]
        if k == "and":
            if False in vals:
                return False
            return None if None in vals else True
        if True in vals:
            return True
        return None if None in vals else False
    if k in ("cmp", "in"):
        a = evaluate(n[2] if k == "cmp" else n[1], value_of)
        bs = [evaluate(n[3], value_of)] if k == "cmp" else [evaluate(v, value_of) for v in n[2]]
        if a is None or any(b is None for b in bs):
            return None
        if isinstance(a, str):
            a = a.rstrip()
            bs = [b.rstrip() if isinstance(b, str) else b for b in bs]
        op = n[1] if k == "cmp" else "="
        try:
            if op == "=":
                return any(a == b for b in bs)
            b = bs[0]
            return {"/=": a != b, "<": a < b, "<=": a <= b, ">": a > b, ">=": a >= b}[op]
        except TypeError:
            return None
    return None


# --------------------------------------------------------------------------
# Clauses and alternatives
# --------------------------------------------------------------------------

def clauses(text, group_of):
    """Split a Fortran condition into readable top-level '.and.' clauses.
    group_of(var) is the namelist group of a variable, or None."""
    text = norm_ws(text or "")
    if not text:
        return []
    n = parse(text)
    if n is None:
        names = re.findall(r"[a-z_]\w*", text.lower())
        return [Clause(text, _groups(names, group_of), None)]
    out = []
    for it in (n[1] if n[0] == "and" else [n]):
        c = Clause(show(it), _groups(variables(it), group_of), it)
        if c not in out:
            out.append(c)
    return out


def _groups(names, group_of):
    out = []
    for v in names:
        g = group_of(v)
        if g is not None and (v, g) not in out:
            out.append((v, g))
    return tuple(out)


def md_code(s):
    s = str(s)
    return "`" + s.replace("`", "'").replace("|", "\\|") + "`" if s else ""


def clause_label(clause, md=True):
    """A clause followed by the namelist groups of its variables."""
    groups = []
    for v, g in clause.groups:
        g = re.sub(r"_nml$", "", g)
        if g not in groups:
            groups.append(g)
    s = md_code(clause.text) if md else clause.text
    return s + (" (%s)" % ", ".join(groups) if groups else "")


def alt_label(alt, md=True):
    return " and ".join(clause_label(c, md) for c in alt)


def alts_label(alts, md=True, sep=" or "):
    """Readable form of a list of alternatives ('' if one of them is empty,
    i.e. the field is registered unconditionally)."""
    if not alts or any(len(a) == 0 for a in alts):
        return ""
    # x = 'a' or x = 'b'  ->  x in ('a', 'b')
    ms = [SELECTOR_RE.match(a[0].text) if len(a) == 1 else None for a in alts]
    if len(alts) > 1 and all(ms) and len({(m.group(1), a[0].groups) for m, a in zip(ms, alts)}) == 1:
        vals = sorted({m.group(2) for m in ms})
        text = "%s in (%s)" % (ms[0].group(1), ", ".join("'%s'" % v for v in vals))
        return clause_label(Clause(text, alts[0][0].groups, None), md)
    return sep.join(alt_label(a, md) for a in alts)


def alts_eval(alts, value_of):
    """True if some alternative holds, False if none can, None if unknown.
    value_of(group, name) gives the value of a namelist variable."""
    if not alts or any(len(a) == 0 for a in alts):
        return True
    res = []
    for a in alts:
        vals = []
        for c in a:
            groups = dict(c.groups)
            vals.append(evaluate(c.node, lambda name, g=groups: value_of(g.get(name), name))
                        if c.node is not None else None)
        if False in vals:
            res.append(False)
        elif None in vals:
            res.append(None)
        else:
            res.append(True)
    if True in res:
        return True
    return None if None in res else False


def split_selector(alt):
    """((var, value, clause), other clauses) if exactly one clause of alt has
    the form var = 'value', else (None, alt)."""
    sel = [(SELECTOR_RE.match(c.text), c) for c in alt]
    sel = [(m.group(1), m.group(2), c) for m, c in sel if m]
    if len(sel) != 1:
        return None, alt
    return sel[0], tuple(c for c in alt if c is not sel[0][2])
