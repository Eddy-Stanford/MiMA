"""Free-form Fortran source model: statements, scopes, declarations, block
conditions and doc comments.  Knows nothing about diag_manager; see
diagnostics.py and namelists.py for the users of this model."""

import io
import os
import re
from collections import defaultdict


# --------------------------------------------------------------------------
# Low-level text handling
# --------------------------------------------------------------------------

def mask_strings(s):
    """Copy of s with the contents of string literals replaced by '_' (quotes
    kept), so that indices line up with s."""
    out = []
    q = None
    i = 0
    n = len(s)
    while i < n:
        c = s[i]
        if q is None:
            if c in ("'", '"'):
                q = c
            out.append(c)
        else:
            if c == q:
                if i + 1 < n and s[i + 1] == q:      # doubled quote escape
                    out.append("__")
                    i += 2
                    continue
                q = None
                out.append(c)
            else:
                out.append("_")
        i += 1
    return "".join(out)


def comment_start(line, in_string=None):
    """Index of the '!' that starts a comment on a physical line (or -1), and
    the quote state at the end of the code part.  in_string is the quote
    character if the line starts inside a continued string."""
    q = in_string
    i = 0
    n = len(line)
    while i < n:
        c = line[i]
        if q is None:
            if c == "!":
                return i, None
            if c in ("'", '"'):
                q = c
        else:
            if c == q:
                if i + 1 < n and line[i + 1] == q:
                    i += 2
                    continue
                q = None
        i += 1
    return -1, q


def match_paren(masked, open_idx):
    """Index of the ')' matching '(' at open_idx in a masked string."""
    depth = 0
    for i in range(open_idx, len(masked)):
        c = masked[i]
        if c in "([":
            depth += 1
        elif c in ")]":
            depth -= 1
            if depth == 0:
                return i
    return -1


def paren_group(text, start):
    """(inner text, index of the closing paren) for the '(' at start."""
    e = match_paren(mask_strings(text), start)
    if e < 0:
        return text[start + 1:], len(text)
    return text[start + 1:e], e


def split_top(text, sep=","):
    """Split text on top-level separators (outside parens and strings)."""
    m = mask_strings(text)
    parts = []
    depth = 0
    last = 0
    i = 0
    while i < len(m):
        c = m[i]
        if c in "([":
            depth += 1
        elif c in ")]":
            depth -= 1
        elif depth == 0 and m.startswith(sep, i):
            parts.append(text[last:i])
            last = i + len(sep)
            i += len(sep)
            continue
        i += 1
    parts.append(text[last:])
    return parts


def norm_ws(s):
    return re.sub(r"\s+", " ", s).strip()


def unquote(lit):
    q = lit[0]
    return lit[1:-1].replace(q + q, q)


def array_items(expr):
    """Items of an array constructor '(/ a, b /)' or '[a, b]', else None."""
    e = norm_ws(expr or "")
    if e.startswith("(/") and e.endswith("/)"):
        return split_top(e[2:-2])
    if e.startswith("[") and e.endswith("]"):
        return split_top(e[1:-1])
    return None


# --------------------------------------------------------------------------
# Statements
# --------------------------------------------------------------------------

class Statement(object):
    __slots__ = ("text", "line", "segs", "is_cpp", "label")

    def __init__(self, text, line, segs, is_cpp=False, label=None):
        self.text = text          # joined code (original case)
        self.line = line          # first physical line (1-based)
        self.segs = segs          # [(offset in text, physical line)]
        self.is_cpp = is_cpp
        self.label = label        # numeric statement label, if any

    def line_at(self, offset):
        """Physical line of the character at offset in text."""
        ln = self.line
        for off, pl in self.segs:
            if off <= offset:
                ln = pl
            else:
                break
        return ln


def logical_statements(lines):
    """Join continuation lines, drop comments and split ';' statements."""
    stmts = []
    buf = ""
    segs = []
    start = None
    in_str = None
    cont = False
    for idx, raw in enumerate(lines, 1):
        line = raw.rstrip("\n\r")
        if not cont and line.lstrip().startswith("#"):
            stmts.append(Statement(line.lstrip(), idx, [(0, idx)], is_cpp=True))
            continue
        k, q = comment_start(line, in_str)
        code = line[:k] if k >= 0 else line
        if cont:
            s = code.lstrip()
            if s == "" and in_str is None:
                continue                      # blank / comment line inside continuation
            if s.startswith("&"):
                code = s[1:]
            elif in_str is None:
                code = " " + s
        else:
            if code.strip() == "":
                continue
            start = idx
        rs = code.rstrip()
        if rs.endswith("&"):
            segs.append((len(buf), idx))
            buf += rs[:-1]
            in_str = q
            cont = True
            continue
        segs.append((len(buf), idx))
        buf += code
        in_str = None
        cont = False
        _split_semicolons(buf, start, segs, stmts)
        buf = ""
        segs = []
    if buf.strip():
        _split_semicolons(buf, start, segs, stmts)
    return stmts


def _split_semicolons(buf, start, segs, stmts):
    m = mask_strings(buf)
    parts = []
    last = 0
    for i, c in enumerate(m):
        if c == ";":
            parts.append((last, i))
            last = i + 1
    parts.append((last, len(buf)))
    for a, b in parts:
        text = buf[a:b]
        if not text.strip():
            continue
        lead = len(text) - len(text.lstrip())
        text = text.strip()
        a2 = a + lead
        rel = []
        for off, pl in segs:
            if off >= b:
                break
            rel.append((max(0, off - a2), pl))
        label = None
        mlab = re.match(r"(\d+)\s+", text)
        if mlab:
            label = mlab.group(1)
            cut = mlab.end()
            text = text[cut:]
            rel = [(max(0, off - cut), pl) for off, pl in rel]
        line = rel[0][1] if rel else start
        for off, pl in rel:
            if off <= 0:
                line = pl
        stmts.append(Statement(text, line, rel or [(0, start)], label=label))


# --------------------------------------------------------------------------
# Doc comments (FORD style)
# --------------------------------------------------------------------------

def doc_comment(lines, lineno):
    """The '!!' documentation attached to physical line lineno (1-based): the
    trailing '!!' comment on that line and the comment-only '!!' lines that
    follow it.  Returns '' when there is none."""
    parts = []
    k, _ = comment_start(lines[lineno - 1].rstrip("\n"))
    trailing = lines[lineno - 1][k:].rstrip() if k >= 0 else ""
    if trailing.startswith("!!"):
        parts.append(trailing[2:].strip())
    i = lineno
    while i < len(lines):
        s = lines[i].strip()
        if not s.startswith("!!"):
            break
        parts.append(s[2:].strip())
        i += 1
    return norm_ws(" ".join(parts))


# --------------------------------------------------------------------------
# Declarations and scopes
# --------------------------------------------------------------------------

TYPE_DECL_RE = re.compile(
    r"^(integer|real|logical|character|complex|double\s+precision|type\s*\()", re.I)


class Symbol(object):
    __slots__ = ("name", "expr", "charlen", "is_param", "dims", "kind", "spec", "line")

    def __init__(self, name, expr, charlen, is_param, dims, kind, spec, line):
        self.name = name
        self.expr = expr          # initialiser text, or None
        self.charlen = charlen
        self.is_param = is_param
        self.dims = dims          # dimension text, or None
        self.kind = kind          # integer, real, logical, character, type, ...
        self.spec = spec          # the type spec, e.g. 'real(kind=8)', 'character(len=64)'
        self.line = line          # physical line of the entity


def parse_declaration(st, text=None):
    """Parse a type declaration statement.  Returns a list of Symbols, or None
    if the statement is not a declaration."""
    text = st.text if text is None else text
    base = len(st.text) - len(text)
    if "::" not in text or not TYPE_DECL_RE.match(text):
        return None
    m = mask_strings(text)
    dc = m.index("::")
    spec = text[:dc]
    spec_l = spec.lower()
    kind = spec_l.split("(")[0].split(",")[0].split("*")[0].strip()
    is_param = bool(re.search(r"\bparameter\b", spec_l))
    charlen = None
    if kind == "character":
        mm = re.match(r"character\s*\(\s*(?:len\s*=\s*)?([^,)]+)", spec, re.I)
        if mm:
            charlen = mm.group(1).strip()
        else:
            mm = re.match(r"character\s*\*\s*(\d+)", spec, re.I)
            charlen = mm.group(1) if mm else "1"
    gdims = None
    md = re.search(r"\bdimension\s*\(", spec, re.I)
    if md:
        o = md.end() - 1
        e = match_paren(mask_strings(spec), o)
        gdims = spec[o + 1:e] if e > 0 else None
    type_spec = norm_ws(split_top(spec)[0])
    syms = []
    off = dc + 2
    for ent in split_top(text[dc + 2:]):
        ent_off = off + len(ent) - len(ent.lstrip())
        off += len(ent) + 1
        ent = ent.strip()
        mm = re.match(r"(\w+)\s*", ent)
        if not mm:
            continue
        name = mm.group(1).lower()
        rest = ent[mm.end():]
        dims = gdims
        if rest.startswith("("):
            e = match_paren(mask_strings(rest), 0)
            dims = rest[1:e]
            rest = rest[e + 1:].lstrip()
        clen = charlen
        mm2 = re.match(r"\*\s*(\d+|\(\s*\*\s*\))\s*", rest)
        if mm2:
            clen = mm2.group(1)
            rest = rest[mm2.end():]
        expr = None
        if rest.startswith("=") and not rest.startswith("=>"):
            expr = rest[1:].strip()
        syms.append(Symbol(name, expr, clen, is_param, dims, kind, type_spec,
                           st.line_at(base + ent_off)))
    return syms


class Scope(object):
    def __init__(self, kind, name, parent, line):
        self.kind = kind        # 'file', 'module', 'program', 'subroutine', 'function'
        self.name = name
        self.parent = parent
        self.line = line
        self.symbols = {}
        self.assigns = defaultdict(list)   # name -> [(stmt index, rhs, lhs subscript)]
        self.uses = []
        self.use_only = {}                 # module -> set of names, or None (no ONLY)
        self.tracer_vars = {}              # var -> role ('tracer', 'tracer_longname', ...)
        self.guards = []                   # (stmt index, cond) for "if (c) return"
        self.calls = []                    # (callee, stmt index, conds, line)

    def module_scope(self):
        s = self
        while s is not None and s.kind not in ("module", "program", "file"):
            s = s.parent
        return s

    def proc_scope(self):
        s = self
        while s is not None and s.kind not in ("subroutine", "function"):
            s = s.parent
        return s

    def chain(self):
        s = self
        while s is not None:
            yield s
            s = s.parent


class Namelist(object):
    """A namelist group statement: namelist /group/ var, ..."""
    __slots__ = ("group", "vars", "spellings", "scope", "line")

    def __init__(self, group, vars, spellings, scope, line):
        self.group = group
        self.vars = vars              # lower case
        self.spellings = spellings    # as written in the namelist statement
        self.scope = scope
        self.line = line


class ExecStmt(object):
    """An executable statement with the conditions and loops around it."""
    __slots__ = ("index", "stmt", "scope", "conds", "loops", "body", "offset")

    def __init__(self, index, stmt, scope, conds, loops, body, offset):
        self.index = index        # index into ParsedFile.stmts
        self.stmt = stmt
        self.scope = scope
        self.conds = conds        # Fortran condition texts, all of which hold
        self.loops = loops        # enclosing DO headers, e.g. 'n = 1, 3'
        self.body = body          # the statement without a one-line IF prefix
        self.offset = offset      # offset of body within stmt.text

    def line_at(self, pos):
        """Physical line of body[pos]."""
        return self.stmt.line_at(self.offset + pos)


class ParsedFile(object):
    def __init__(self, path, relpath):
        self.path = path
        self.relpath = relpath
        self.lines = []
        self.stmts = []
        self.scopes = []
        self.execs = []             # ExecStmt
        self.namelist_stmts = []    # Namelist
        self.namelists = {}         # var -> group
        self.modules = []
        self.proc_defs = {}         # procedure name -> Scope


# --------------------------------------------------------------------------
# File parsing
# --------------------------------------------------------------------------

UNIT_START_RE = re.compile(
    r"^(?:(?:recursive|pure|elemental|impure)\s+|"
    r"(?:integer|real|logical|character|complex|double\s+precision|type\s*\([^)]*\))"
    r"(?:\s*\([^)]*\)|\s*\*\s*\d+)?\s+)*(subroutine|function)\s+(\w+)", re.I)
MODULE_RE = re.compile(r"^module\s+(\w+)\s*$", re.I)
PROGRAM_RE = re.compile(r"^program\s+(\w+)", re.I)
END_UNIT_RE = re.compile(r"^end\s*(subroutine|function|module|program)\b", re.I)
END_BARE_RE = re.compile(r"^end\s*$", re.I)
INTERFACE_RE = re.compile(r"^(abstract\s+)?interface\b", re.I)
END_INTERFACE_RE = re.compile(r"^end\s*interface\b", re.I)
CONSTRUCT_LABEL_RE = re.compile(r"^(\w+)\s*:(?!:)\s*")
IF_THEN_RE = re.compile(r"^if\s*\(", re.I)
ELSEIF_RE = re.compile(r"^else\s*if\s*\(", re.I)
ELSE_RE = re.compile(r"^else\s*(\w+)?\s*$", re.I)
ENDIF_RE = re.compile(r"^end\s*if\b", re.I)
DO_RE = re.compile(r"^do\b(?!\w)", re.I)
ENDDO_RE = re.compile(r"^end\s*do\b", re.I)
SELECT_RE = re.compile(r"^select\s*case\s*\(", re.I)
CASE_RE = re.compile(r"^case\b", re.I)
ENDSELECT_RE = re.compile(r"^end\s*select\b", re.I)
WHERE_RE = re.compile(r"^(where|forall)\s*\(", re.I)
ENDWHERE_RE = re.compile(r"^end\s*(where|forall)\b", re.I)
ELSEWHERE_RE = re.compile(r"^else\s*where\b", re.I)
NAMELIST_RE = re.compile(r"^namelist\s*/", re.I)
NAMELIST_GROUP_RE = re.compile(r"/\s*(\w+)\s*/")
USE_RE = re.compile(r"^use\b\s*(?:,\s*\w+\s*::)?\s*(\w+)", re.I)
CALL_RE = re.compile(r"^call\s+(\w+)", re.I)
ASSIGN_RE = re.compile(r"^(\w+)\s*(\([^=]*\))?\s*=(?!=)", re.I)
KEYWORD_STMT = re.compile(
    r"^(if|do|else|end|select|case|where|forall|call|return|go\s*to|print|write|"
    r"read|open|close|allocate|deallocate|nullify|use|implicit|public|private|"
    r"save|data|namelist|contains|cycle|exit|stop|format|include|interface|"
    r"module|subroutine|function|program|type|integer|real|logical|character|"
    r"complex|double|equivalence|common|external|intrinsic|optional|pointer|"
    r"target|parameter|entry|continue|inquire|rewind|backspace|endfile)\b", re.I)
# Units of internal-file WRITEs that are real output units, not character variables.
OUTPUT_UNITS = ("unit", "stdout", "logunit", "unit_log")
# get_tracer_names(model, n, name, longname, units): the variables it fills.
TRACER_NAME_ROLES = {"name": "tracer", "longname": "tracer_longname", "units": "tracer_units"}


def neg(c):
    """Fortran text of .not. c, simplified for a plain (non-compound) c."""
    c = norm_ws(c)
    m = re.match(r"^\.not\.\s*(.+)$", c, re.I)
    if m:
        inner = m.group(1).strip()
        if inner.startswith("(") and match_paren(mask_strings(inner), 0) == len(inner) - 1:
            inner = inner[1:-1].strip()
        if not re.search(r"\.(and|or|eqv|neqv)\.", inner, re.I):
            return "(%s)" % inner
    return "(.not. (%s))" % c


# Guards that only protect against double initialisation.
TRIVIAL_GUARD_RE = re.compile(r"^\s*\(?\s*module_is_initialized\s*\)?\s*$", re.I)


# Call-site conditions that are always true in the model: atmosphere PEs
# (MiMA runs only the atmosphere) and first initialisation.
TRIVIAL_COND_RE = re.compile(
    r"^\s*\(\s*(atm%pe|\.not\.\s*\(?\s*module_is_initialized\s*\)?)\s*\)\s*$", re.I)


def is_trivial_cond(c):
    return bool(TRIVIAL_COND_RE.match(c))


def guard_conds(guards, before=None):
    """Conditions implied by the 'if (c) return' guards before a statement."""
    return [neg(g) for (gi, g) in guards
            if (before is None or gi < before) and not TRIVIAL_GUARD_RE.match(g)]


class Block(object):
    __slots__ = ("kind", "cond", "prior", "do_label", "header", "cases")

    def __init__(self, kind, cond=None, header=None, do_label=None):
        self.kind = kind
        self.cond = cond        # condition of the current branch
        self.prior = []         # conditions of earlier IF branches
        self.do_label = do_label
        self.header = header    # DO header or SELECT CASE expression
        self.cases = []         # SELECT CASE: selectors of earlier cases

    def describe(self):
        parts = [neg(p) for p in self.prior]
        if self.cond is not None:
            parts.append("(%s)" % self.cond if self.kind != "select" else self.cond)
        return " .and. ".join(parts)


def _cpp_condition(stmt, cpp):
    """Update the stack of #if conditions for a preprocessor line."""
    text = stmt.text
    d = text[1:].strip().lower()
    arg = text.split(None, 1)[1].strip() if len(text.split(None, 1)) > 1 else "?"
    if d.startswith("ifdef"):
        cpp.append(["defined(%s)" % arg])
    elif d.startswith("ifndef"):
        cpp.append(["!defined(%s)" % arg])
    elif d.startswith("if"):
        cpp.append([text[1:].strip()[2:].strip()])
    elif d.startswith("elif") and cpp:
        cpp[-1] = ["!(%s)" % " && ".join(cpp[-1]), text[1:].strip()[4:].strip()]
    elif d.startswith("else") and cpp:
        cpp[-1] = ["!(%s)" % " && ".join(cpp[-1])]
    elif d.startswith("endif") and cpp:
        cpp.pop()


def _tracer_name_vars(body, scope):
    """Record the variables that get_tracer_names fills (tracer name, long
    name, units), so that strings built from them are shown as <tracer>."""
    o = body.find("(")
    if o < 0:
        return
    inner, _ = paren_group(body, o)
    args = [a.strip() for a in split_top(inner)]
    roles = ["tracer", "tracer_longname", "tracer_units"]
    pos = [a for a in args if not re.match(r"^\w+\s*=", a)]
    for k, a in enumerate(pos[2:5]):
        if re.match(r"^\w+$", a):
            scope.tracer_vars[a.lower()] = roles[k]
    for a in args:
        mk = re.match(r"^(name|longname|units)\s*=\s*(\w+)$", a, re.I)
        if mk:
            scope.tracer_vars[mk.group(2).lower()] = TRACER_NAME_ROLES[mk.group(1).lower()]


def parse_file(path, relpath):
    with io.open(path, "r", encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()
    pf = ParsedFile(path, relpath)
    pf.lines = lines
    pf.stmts = stmts = logical_statements(lines)
    root = Scope("file", os.path.basename(path), None, 0)
    scope = root
    pf.scopes.append(root)
    blocks = []
    cpp = []           # stack of [condition texts] of #if blocks
    in_interface = 0

    for si, st in enumerate(stmts):
        if st.is_cpp:
            _cpp_condition(st, cpp)
            continue
        text = st.text

        # labelled DO termination (e.g. "10 continue")
        if st.label:
            while blocks and blocks[-1].kind == "do" and blocks[-1].do_label == st.label:
                blocks.pop()

        body = text
        mlab = CONSTRUCT_LABEL_RE.match(body)
        if mlab and not re.match(r"^(\w+)\s*:\s*:", body) and \
                mlab.group(1).lower() not in ("case", "default"):
            body = body[mlab.end():]          # construct name, e.g. "outer: do i=1,n"
        low = body.lower()

        # ---- interface blocks: skip the procedure headers inside them
        if INTERFACE_RE.match(body):
            in_interface += 1
            continue
        if END_INTERFACE_RE.match(body):
            in_interface = max(0, in_interface - 1)
            continue
        if in_interface:
            continue

        # ---- program units
        mm = MODULE_RE.match(body)
        if mm and not low.startswith("module procedure"):
            scope = Scope("module", mm.group(1).lower(), scope, st.line)
            pf.scopes.append(scope)
            pf.modules.append(scope.name)
            blocks = []
            continue
        mm = PROGRAM_RE.match(body)
        if mm:
            scope = Scope("program", mm.group(1).lower(), scope, st.line)
            pf.scopes.append(scope)
            pf.modules.append(scope.name)
            blocks = []
            continue
        if END_UNIT_RE.match(body) or END_BARE_RE.match(body):
            if scope.parent is not None:
                scope = scope.parent
            blocks = []
            continue
        mm = UNIT_START_RE.match(body)
        if mm and not low.startswith("end"):
            scope = Scope(mm.group(1).lower(), mm.group(2).lower(), scope, st.line)
            pf.scopes.append(scope)
            pf.proc_defs.setdefault(scope.name, scope)
            blocks = []
            continue

        # ---- namelist /a/ x, y [/b/ z]
        if NAMELIST_RE.match(body):
            heads = list(NAMELIST_GROUP_RE.finditer(body))
            for k, mg in enumerate(heads):
                end = heads[k + 1].start() if k + 1 < len(heads) else len(body)
                grp = mg.group(1).lower()
                spellings = [v.strip() for v in body[mg.end():end].split(",")]
                spellings = [v for v in spellings if re.match(r"^\w+$", v)]
                names = [v.lower() for v in spellings]
                for v in names:
                    pf.namelists[v] = grp
                pf.namelist_stmts.append(Namelist(grp, names, spellings, scope, st.line))
            continue

        # ---- use statements
        mm = USE_RE.match(body)
        if mm and not re.match(r"^use\w", low):
            umod = mm.group(1).lower()
            scope.uses.append(umod)
            mo = re.search(r"\bonly\s*:(.*)$", body, re.I)
            if mo:
                names = set()
                for item in mo.group(1).split(","):
                    for nm in item.split("=>"):
                        nm = nm.strip().lower()
                        if re.match(r"^\w+$", nm):
                            names.add(nm)
                prev = scope.use_only.get(umod, set())
                if prev is not None:
                    scope.use_only[umod] = prev | names
            else:
                scope.use_only[umod] = None
            continue

        # ---- declarations
        decl = parse_declaration(st, body)
        if decl is not None:
            for sym in decl:
                scope.symbols[sym.name] = sym
            continue

        # ---- data statements (simple form)
        if re.match(r"^data\b", low):
            for mmd in re.finditer(r"(\w+)\s*/([^/]*)/", body):
                nm = mmd.group(1).lower()
                sym = scope.symbols.get(nm)
                if sym is None:
                    sym = Symbol(nm, None, None, False, None, "?", "", st.line)
                    scope.symbols[nm] = sym
                sym.expr = "(/" + mmd.group(2) + "/)"
            continue

        cond_now = [b.describe() for b in blocks if b.kind in ("if", "select", "where")]
        cond_now += ["#if " + " && ".join(c) for c in cpp]
        loops_now = [b.header for b in blocks if b.kind == "do"]

        # ---- block constructs
        if ELSEIF_RE.match(body):
            c, _ = paren_group(body, body.index("("))
            if blocks and blocks[-1].kind == "if":
                b = blocks[-1]
                if b.cond is not None:
                    b.prior.append(norm_ws(b.cond))
                b.cond = norm_ws(c)
            continue
        if ELSEWHERE_RE.match(body):
            continue
        if ELSE_RE.match(body):
            if blocks and blocks[-1].kind == "if":
                b = blocks[-1]
                if b.cond is not None:
                    b.prior.append(norm_ws(b.cond))
                b.cond = None
            continue
        ends = ((ENDIF_RE, "if"), (ENDDO_RE, "do"), (ENDSELECT_RE, "select"),
                (ENDWHERE_RE, "where"))
        kind_ended = next((k for rx, k in ends if rx.match(body)), None)
        if kind_ended:
            if blocks and blocks[-1].kind == kind_ended:
                blocks.pop()
            continue
        if SELECT_RE.match(body):
            c, _ = paren_group(body, body.index("("))
            blocks.append(Block("select", header=norm_ws(c)))
            continue
        if CASE_RE.match(body) and blocks and blocks[-1].kind == "select":
            b = blocks[-1]
            sel = norm_ws(body[4:])
            if sel.lower().startswith("default"):
                # the default branch: none of the listed cases
                listed = [s[1:-1].strip() if s.startswith("(") and s.endswith(")") else s
                          for s in b.cases]
                b.cond = "(.not. (%s == (%s)))" % (b.header, ", ".join(listed)) if listed else None
            else:
                b.cases.append(sel)
                b.cond = "%s == %s" % (b.header, sel)
            continue
        if DO_RE.match(body):
            mdl = re.match(r"^do\s+(\d+)\b\s*,?", body, re.I)
            hdr = norm_ws(body[mdl.end():] if mdl else body[2:])
            blocks.append(Block("do", header=hdr, do_label=mdl.group(1) if mdl else None))
            continue
        if WHERE_RE.match(body):
            c, e = paren_group(body, body.index("("))
            if body[e + 1:].strip() == "":
                blocks.append(Block("where", norm_ws(c)))
                continue

        # ---- if ... then, or a one-line if
        extra = []
        exec_body = body
        if IF_THEN_RE.match(body):
            c, e = paren_group(body, body.index("("))
            rest = body[e + 1:].strip()
            if rest.lower() == "then":
                blocks.append(Block("if", norm_ws(c)))
                continue
            extra = ["(%s)" % norm_ws(c)]
            exec_body = rest
            if re.match(r"^return\b", rest, re.I) and scope.kind in ("subroutine", "function") \
                    and not blocks:
                scope.guards.append((si, norm_ws(c)))
        conds = cond_now + extra
        pf.execs.append(ExecStmt(si, st, scope, conds, loops_now, exec_body,
                                 len(text) - len(exec_body)))

        # ---- call sites (for the call-chain conditions)
        mc = CALL_RE.match(exec_body)
        if mc:
            callee = mc.group(1).lower()
            scope.calls.append((callee, si, list(conds), st.line))
            if callee == "get_tracer_names":
                _tracer_name_vars(exec_body, scope)

        # ---- internal-file writes:  write (chvers, '(i2)') n  ->  chvers = {n}
        mw = re.match(r"^write\s*\(\s*(\w+)\s*,", exec_body, re.I)
        if mw:
            _, e = paren_group(exec_body, exec_body.index("("))
            var = mw.group(1).lower()
            if not re.match(r"^\d+$", var) and var not in OUTPUT_UNITS:
                scope.assigns[var].append((si, WriteValue(norm_ws(exec_body[e + 1:])), None))

        # ---- simple assignments
        ma = ASSIGN_RE.match(exec_body)
        if ma and not KEYWORD_STMT.match(exec_body):
            sub = ma.group(2)
            scope.assigns[ma.group(1).lower()].append(
                (si, exec_body[ma.end():].strip(), sub.strip() if sub else None))
    return pf


class WriteValue(str):
    """The value list of an internal-file WRITE: its text is only known at
    run time."""


# --------------------------------------------------------------------------
# Source files
# --------------------------------------------------------------------------

def find_sources(src):
    out = []
    for dirpath, dirnames, filenames in os.walk(src):
        dirnames.sort()
        for f in sorted(filenames):
            if f.endswith((".f90", ".F90")):
                out.append(os.path.join(dirpath, f))
    return out


def built_sources(root, src):
    """Paths (relative to root) of the Fortran files listed in the
    CMakeLists.txt files under src, or None if there are none."""
    built = set()
    found_any = False
    for dirpath, dirnames, filenames in os.walk(src):
        if "CMakeLists.txt" not in filenames:
            continue
        found_any = True
        with io.open(os.path.join(dirpath, "CMakeLists.txt"), encoding="utf-8",
                     errors="replace") as fh:
            for line in fh:
                s = line.split("#", 1)[0].strip()
                for tok in re.findall(r"[\w./${}-]+\.(?:f90|F90|f|F)\b", s):
                    if "${" in tok:
                        continue
                    for base in (dirpath, root):
                        cand = os.path.normpath(os.path.join(base, tok))
                        if os.path.exists(cand):
                            built.add(os.path.relpath(cand, root))
                            break
    return built if found_any else None


def parse_tree(root, src):
    """Parse every Fortran file under src.  Returns (parsed files, module name
    -> module Scope)."""
    parsed = [parse_file(f, os.path.relpath(f, root)) for f in find_sources(src)]
    modtab = {}
    for pf in parsed:
        for s in pf.scopes:
            if s.kind == "module":
                modtab.setdefault(s.name, s)
    return parsed, modtab
