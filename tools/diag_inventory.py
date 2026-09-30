#!/usr/bin/env python3
"""
diag_inventory.py -- build an inventory of every FMS diagnostic field that the
MiMA Fortran sources can register (register_diag_field / register_static_field).

Stdlib only (Python >= 3.7).  Lives at tools/diag_inventory.py; the repository
root defaults to the directory above this script, so it can be run from anywhere:

    python3 tools/diag_inventory.py                       # writes docs/Diagnostics.md
    python3 tools/diag_inventory.py --check               # CI: exit 1 if docs stale
    python3 tools/diag_inventory.py --csv diag_fields.csv # also write CSV
    python3 tools/diag_inventory.py --summary summary.md  # also write a QA summary
    python3 tools/diag_inventory.py --validate-diag-table diag_table [--nml input.nml]

--validate-diag-table only checks (it does not rewrite docs/Diagnostics.md): every
field line must name a registered module/field, use a defined file, a valid
reduction and packing, and a unique output name.  With --nml, the registration
conditions are evaluated with the settings in that input.nml (namelist variables
not set there take their default values from the source).

What it does
------------
* Reads every free-form Fortran file (*.f90, *.F90) under --src.
* Joins continuation lines, strips comments, splits ';' statements, tracks
  program units (module / subroutine / function), block constructs
  (if/else, do, select case, where, #ifdef) and early-return guards.
* For every register_diag_field / register_static_field call it extracts
  module_name, field_name, axes, init_time, long_name, units, missing_value,
  range, standard_name, mask_variant (positional or keyword).
* String arguments are evaluated: literals, '//' concatenation, trim/adjustl,
  character parameters / initialised variables (procedure scope, module scope,
  then modules brought in by USE), local assignments preceding the call, and
  elements of literal array constructors indexed by a DO loop variable (the
  registration is expanded once per element).  Anything left is rendered as a
  placeholder, e.g. '{tname}_every' (<tracer> when the variable is filled by
  get_tracer_names), and the row is marked dynamic.
* Conditions: enclosing IF/SELECT/#if blocks, one-line IFs, and preceding
  "if (...) return" guards in the same procedure are recorded.  The call chain
  of the enclosing procedure is walked upward (up to --max-call-depth) to pick
  up call-site gates such as `if (do_grey_radiation) call grey_radiation_init`.
  Identifiers that are namelist variables are reported separately.
* Heuristic "is it ever sent?" check: the id variable assigned from the
  register call is searched for as the first argument of send_data /
  send_tile_averaged_data in the same file, then anywhere in the tree.
* Built flag: files not listed (or commented out) in the CMakeLists.txt source
  lists are reported as not built.

* Conditions are rewritten in a readable form (radiation_scheme = 'rrtm',
  not do_bm, ...) with the namelist group of each variable, and factored per
  module.  Registrations of the same (module, field) in several places (e.g.
  the gray and RRTM schemes) are merged into one row; differing metadata is
  flagged.

Outputs: Markdown (a short diag_table guide, then the fields grouped by module),
optional CSV, optional summary Markdown.  The Markdown output is deterministic so
--check can be used in CI.
"""

import argparse
import csv
import io
import os
import re
import sys
from collections import OrderedDict, defaultdict

VERSION = "2.0"

# Directories (path components under src/) whose physics is a candidate for
# removal from MiMA.  Rows from these directories get legacy=yes so that docs
# can filter them out.  Override with --legacy-dirs.
DEFAULT_LEGACY_DIRS = [
    "sea_esf_rad", "radiation_driver", "donner_deep", "ras", "strat_cloud",
    "edt", "my25_turb", "diag_cloud", "diag_cloud_rad", "cloud_rad",
    "entrain", "stable_bl_turb",
]

# Files / modules whose internal calls are not "registrations" (the
# diag_manager implementation itself).
SKIP_MODULES = {"diag_manager_mod"}

REG_RE = re.compile(r"\bregister_(diag|static)_field\s*\(", re.I)
SEND_RE = re.compile(r"\bsend_(?:data|tile_averaged_data|global_diag)\s*\(", re.I)

DIAG_ARRAY_ARGS = ["module_name", "field_name", "axes", "init_time", "long_name",
                   "units", "missing_value", "range", "mask_variant",
                   "standard_name", "verbose"]
DIAG_SCALAR_ARGS = ["module_name", "field_name", "init_time", "long_name",
                    "units", "missing_value", "range"]
STATIC_ARGS = ["module_name", "field_name", "axes", "long_name", "units",
               "missing_value", "range", "mask_variant", "require",
               "standard_name", "dynamic"]

CSV_COLUMNS = [
    "module", "field", "kind", "axes", "ndim", "units", "long_name",
    "standard_name", "missing_value", "range", "id_var", "file", "line",
    "subroutine", "fortran_module", "condition", "call_gates", "namelist_vars",
    "loops", "dynamic", "unresolved", "send_status", "built", "legacy",
    "physics_dir", "gate_nml", "available", "dims", "notes",
]


# --------------------------------------------------------------------------
# Low-level Fortran text handling
# --------------------------------------------------------------------------

def mask_strings(s):
    """Return a copy of s where the *contents* of string literals are
    replaced by '_' (quotes kept) so that indices line up with s."""
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


def strip_comment(line, in_string):
    """Strip a '!' comment from a physical line.  in_string is the quote char
    if the line starts inside a continued string.  Returns (code, quote_state
    at end of code)."""
    q = in_string
    i = 0
    n = len(line)
    while i < n:
        c = line[i]
        if q is None:
            if c == "!":
                return line[:i], None
            if c in ("'", '"'):
                q = c
        else:
            if c == q:
                if i + 1 < n and line[i + 1] == q:
                    i += 2
                    continue
                q = None
        i += 1
    return line, q


class Statement(object):
    __slots__ = ("text", "line", "segs", "is_cpp", "label")

    def __init__(self, text, line, segs, is_cpp=False, label=None):
        self.text = text          # joined code (original case)
        self.line = line          # first physical line (1-based)
        self.segs = segs          # list of (offset_in_text, physical_line)
        self.is_cpp = is_cpp
        self.label = label        # numeric statement label, if any

    def line_at(self, offset):
        ln = self.line
        for off, pl in self.segs:
            if off <= offset:
                ln = pl
            else:
                break
        return ln


def logical_statements(lines):
    """Free-form statement assembler."""
    stmts = []
    buf = ""
    segs = []
    start = None
    in_str = None
    cont = False
    for idx, raw in enumerate(lines, 1):
        line = raw.rstrip("\n\r")
        stripped = line.lstrip()
        if not cont and stripped.startswith("#"):
            stmts.append(Statement(stripped, idx, [(0, idx)], is_cpp=True))
            continue
        code, q = strip_comment(line, in_str)
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
        # Collapse to offsets relative to this statement
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
        # physical line of the first character of the statement
        for off, pl in rel:
            if off <= 0:
                line = pl
        stmts.append(Statement(text, line, rel or [(0, start)], label=label))


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


def split_top(text, sep=","):
    """Split text on top-level separators (outside parens/strings)."""
    m = mask_strings(text)
    parts = []
    depth = 0
    last = 0
    i = 0
    L = len(sep)
    while i < len(m):
        c = m[i]
        if c in "([":
            depth += 1
        elif c in ")]":
            depth -= 1
        elif depth == 0 and m.startswith(sep, i):
            parts.append(text[last:i])
            last = i + L
            i += L
            continue
        i += 1
    parts.append(text[last:])
    return parts


def norm_ws(s):
    return re.sub(r"\s+", " ", s).strip()


def unquote(lit):
    q = lit[0]
    return lit[1:-1].replace(q + q, q)


# --------------------------------------------------------------------------
# Symbol tables
# --------------------------------------------------------------------------

TYPE_DECL_RE = re.compile(
    r"^(integer|real|logical|character|complex|double\s+precision|type\s*\()",
    re.I)


class Symbol(object):
    __slots__ = ("name", "expr", "charlen", "is_param", "dims", "kind", "line")

    def __init__(self, name, expr, charlen, is_param, dims, kind, line):
        self.name = name
        self.expr = expr
        self.charlen = charlen
        self.is_param = is_param
        self.dims = dims
        self.kind = kind
        self.line = line


def parse_declaration(text, line):
    """Parse a type declaration statement.  Returns list of Symbols (may be
    empty) or None if the statement is not a declaration."""
    if "::" not in text:
        # old style declarations like "character*8 name" are rare here
        return None
    if not TYPE_DECL_RE.match(text):
        return None
    m = mask_strings(text)
    dc = m.index("::")
    spec = text[:dc]
    ents = text[dc + 2:]
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
    syms = []
    for ent in split_top(ents):
        ent = ent.strip()
        if not ent:
            continue
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
        if rest.startswith("=>"):
            expr = None
        elif rest.startswith("="):
            expr = rest[1:].strip()
        syms.append(Symbol(name, expr, clen, is_param, dims, kind, line))
    return syms


class Scope(object):
    def __init__(self, kind, name, parent, line):
        self.kind = kind        # 'file', 'module', 'program', 'subroutine', 'function'
        self.name = name
        self.parent = parent
        self.line = line
        self.symbols = {}
        self.assigns = defaultdict(list)   # name -> [(stmt_idx, rhs, lhs_subscript)]
        self.uses = []
        self.use_only = {}                 # module -> set of names or None (no ONLY)
        self.tracer_vars = {}              # var -> role ('tracer', 'longname', 'units')
        self.guards = []                   # (stmt_idx, cond) for "if (c) return"
        self.calls = []                    # (callee, stmt_idx, cond_list, line)

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
NAMELIST_RE = re.compile(r"^namelist\s*/\s*(\w+)\s*/(.*)$", re.I)
USE_RE = re.compile(r"^use\b\s*(?:,\s*\w+\s*::)?\s*(\w+)", re.I)
CALL_RE = re.compile(r"^call\s+(\w+)", re.I)
ASSIGN_RE = re.compile(r"^(\w+)\s*(\([^=]*\))?\s*=(?!=)", re.I)
KEYWORD_STMT = re.compile(
    r"^(if|do|else|end|select|case|where|forall|call|return|go\s*to|print|write|"
    r"read|open|close|allocate|deallocate|nullify|use|implicit|public|private|"
    r"save|data|namelist|contains|cycle|exit|stop|format|include|interface|"
    r"module|subroutine|function|program|type|integer|real|logical|character|"
    r"complex|double|equivalence|common|external|intrinsic|optional|pointer|"
    r"target|parameter|entry|continue|inquire|rewind|backspace|endfile)\b",
    re.I)


def paren_group(text, start):
    """Given text and index of '(', return (inner, end_index)."""
    m = mask_strings(text)
    e = match_paren(m, start)
    if e < 0:
        return text[start + 1:], len(text)
    return text[start + 1:e], e


def neg(c):
    c = norm_ws(c)
    m = re.match(r"^\.not\.\s*(.+)$", c, re.I)
    if m:
        inner = m.group(1).strip()
        if inner.startswith("(") and match_paren(mask_strings(inner), 0) == len(inner) - 1:
            inner = inner[1:-1].strip()
        if not re.search(r"\.(and|or|eqv|neqv)\.", inner, re.I):
            return "(%s)" % inner
    return "(.not. (%s))" % c


TRIVIAL_GUARD_RE = re.compile(r"^\s*\(?\s*module_is_initialized\s*\)?\s*$", re.I)
# Conditions that are always true on atmosphere PEs and only add noise.
TRIVIAL_COND_RE = re.compile(r"^\s*\(\s*(atm%pe|\.not\.\s*\(?\s*module_is_initialized\s*\)?)\s*\)\s*$", re.I)


def guard_list(guards, before=None):
    return [neg(g) for (gi, g) in guards
            if (before is None or gi < before) and not TRIVIAL_GUARD_RE.match(g)]


class Block(object):
    __slots__ = ("kind", "cond", "prior", "label", "do_label", "header", "line")

    def __init__(self, kind, cond=None, header=None, do_label=None, line=0):
        self.kind = kind
        self.cond = cond        # current branch condition text
        self.prior = []         # conditions of earlier branches (negated)
        self.label = None
        self.do_label = do_label
        self.header = header
        self.line = line

    def describe(self):
        parts = [neg(p) for p in self.prior]
        if self.cond is not None:
            parts.append("(%s)" % self.cond if self.kind != "select" else self.cond)
        return " .and. ".join(parts)


class RegCall(object):
    def __init__(self):
        self.__dict__.update(dict(
            kind="", args=OrderedDict(), raw="", id_var="", file="", line=0,
            scope=None, stmt_idx=0, conds=[], loops=[], fortran_module="",
            proc=""))


class ParsedFile(object):
    def __init__(self, path, relpath):
        self.path = path
        self.relpath = relpath
        self.stmts = []
        self.scopes = []
        self.regs = []
        self.sends = set()          # base id names sent
        self.idrefs = defaultdict(int)  # base identifier -> count of refs in non-decl statements
        self.namelists = {}         # var -> group
        self.modules = []
        self.proc_defs = {}         # proc name -> Scope


def parse_file(path, relpath):
    with io.open(path, "r", encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()
    pf = ParsedFile(path, relpath)
    stmts = logical_statements(lines)
    pf.stmts = stmts
    root = Scope("file", os.path.basename(path), None, 0)
    scope = root
    pf.scopes.append(root)
    blocks = []
    cpp = []           # list of [cond_text]
    in_interface = 0

    for si, st in enumerate(stmts):
        text = st.text
        low = text.lower()
        if st.is_cpp:
            d = low[1:].strip()
            if d.startswith("ifdef"):
                cpp.append(["defined(%s)" % text.split(None, 1)[1].strip() if len(text.split(None, 1)) > 1 else "?"])
            elif d.startswith("ifndef"):
                cpp.append(["!defined(%s)" % text.split(None, 1)[1].strip() if len(text.split(None, 1)) > 1 else "?"])
            elif d.startswith("if"):
                cpp.append([text[1:].strip()[2:].strip()])
            elif d.startswith("elif"):
                if cpp:
                    prev = cpp[-1]
                    cpp[-1] = ["!(%s)" % " && ".join(prev), text[1:].strip()[4:].strip()]
            elif d.startswith("else"):
                if cpp:
                    cpp[-1] = ["!(%s)" % " && ".join(cpp[-1])]
            elif d.startswith("endif"):
                if cpp:
                    cpp.pop()
            continue

        # handle labelled DO termination (e.g. "10 continue")
        if st.label:
            while blocks and blocks[-1].kind == "do" and blocks[-1].do_label == st.label:
                blocks.pop()

        body = text
        mlab = CONSTRUCT_LABEL_RE.match(body)
        if mlab and not re.match(r"^(\w+)\s*:\s*:", body) and \
                mlab.group(1).lower() not in ("case", "default"):
            # construct name e.g. "outer: do i=1,n"
            body = body[mlab.end():]
            low = body.lower()

        # ---- interfaces: skip procedure headers inside interface blocks
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

        # ---- namelist
        mm = NAMELIST_RE.match(body)
        if mm:
            grp = mm.group(1).lower()
            rest = mm.group(2)
            # possibly more groups: /a/ x, y /b/ z
            for piece in re.split(r"/\s*\w+\s*/", rest):
                for v in piece.split(","):
                    v = v.strip().lower()
                    if re.match(r"^\w+$", v):
                        pf.namelists[v] = grp
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
        decl = parse_declaration(body, st.line)
        if decl is not None:
            for sym in decl:
                scope.symbols[sym.name] = sym
            continue

        # ---- data statements (simple form)
        if re.match(r"^data\b", low):
            for mmd in re.finditer(r"(\w+)\s*/([^/]*)/", body):
                nm = mmd.group(1).lower()
                vals = mmd.group(2)
                sym = scope.symbols.get(nm)
                if sym is None:
                    sym = Symbol(nm, None, None, False, None, "?", st.line)
                    scope.symbols[nm] = sym
                sym.expr = "(/" + vals + "/)"
            continue

        cond_now = [b.describe() for b in blocks if b.kind in ("if", "select", "where")]
        cond_now += ["#if " + " && ".join(c) for c in cpp]
        loops_now = [b.header for b in blocks if b.kind == "do"]

        # ---- block constructs
        is_block_stmt = False
        if ELSEIF_RE.match(body):
            o = body.index("(")
            c, e = paren_group(body, o)
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
        if ENDIF_RE.match(body):
            if blocks and blocks[-1].kind == "if":
                blocks.pop()
            continue
        if ENDDO_RE.match(body):
            if blocks and blocks[-1].kind == "do":
                blocks.pop()
            continue
        if ENDSELECT_RE.match(body):
            if blocks and blocks[-1].kind == "select":
                blocks.pop()
            continue
        if ENDWHERE_RE.match(body):
            if blocks and blocks[-1].kind == "where":
                blocks.pop()
            continue
        if SELECT_RE.match(body):
            o = body.index("(")
            c, e = paren_group(body, o)
            b = Block("select", None, header=norm_ws(c), line=st.line)
            blocks.append(b)
            continue
        if CASE_RE.match(body) and blocks and blocks[-1].kind == "select":
            b = blocks[-1]
            sel = body[4:].strip()
            if sel.lower().startswith("default"):
                b.cond = "%s == case default" % b.header
            else:
                b.cond = "%s == %s" % (b.header, norm_ws(sel))
            continue
        if DO_RE.match(body):
            mdl = re.match(r"^do\s+(\d+)\b\s*,?", body, re.I)
            do_label = mdl.group(1) if mdl else None
            hdr = norm_ws(body[mdl.end():] if mdl else body[2:])
            blocks.append(Block("do", header=hdr, do_label=do_label, line=st.line))
            continue
        m_where = WHERE_RE.match(body)
        if m_where:
            o = body.index("(")
            c, e = paren_group(body, o)
            if body[e + 1:].strip() == "":
                blocks.append(Block("where", norm_ws(c), line=st.line))
                continue
        # if ... then  / one-line if
        extra_cond = []
        exec_body = body
        if IF_THEN_RE.match(body):
            o = body.index("(")
            c, e = paren_group(body, o)
            rest = body[e + 1:].strip()
            if rest.lower() == "then":
                blocks.append(Block("if", norm_ws(c), line=st.line))
                continue
            # one-line if (possibly arithmetic if: ignore)
            extra_cond = [norm_ws(c)]
            exec_body = rest
            if re.match(r"^return\b", rest, re.I) and scope.kind in ("subroutine", "function") \
                    and not blocks:
                scope.guards.append((si, norm_ws(c)))
        elif re.match(r"^return\b", body, re.I):
            pass

        conds = cond_now + ["(%s)" % c for c in extra_cond]

        # ---- call sites (for call-chain gating)
        mc = CALL_RE.match(exec_body)
        if mc:
            scope.calls.append((mc.group(1).lower(), si, list(conds), st.line))
            # tracer names filled by get_tracer_names(model, n, name, longname, units)
            if mc.group(1).lower() == "get_tracer_names":
                o = exec_body.find("(")
                if o > 0:
                    inner, _ = paren_group(exec_body, o)
                    args = [a.strip() for a in split_top(inner)]
                    roles = ["tracer", "tracer_longname", "tracer_units"]
                    pos = [a for a in args if not re.match(r"^\w+\s*=", a)]
                    for k, a in enumerate(pos[2:5]):
                        if re.match(r"^\w+$", a):
                            scope.tracer_vars[a.lower()] = roles[k]
                    for a in args:
                        mk = re.match(r"^(name|longname|units)\s*=\s*(\w+)$", a, re.I)
                        if mk:
                            r = {"name": "tracer", "longname": "tracer_longname",
                                 "units": "tracer_units"}[mk.group(1).lower()]
                            scope.tracer_vars[mk.group(2).lower()] = r
        # internal-file writes:  write (chvers, '(i2)') n  -> chvers = {n}
        mw = re.match(r"^write\s*\(\s*(\w+)\s*,", exec_body, re.I)
        if mw:
            o = exec_body.index("(")
            inner, e = paren_group(exec_body, o)
            what = norm_ws(exec_body[e + 1:])
            var = mw.group(1).lower()
            if not re.match(r"^\d+$", var) and var not in ("unit", "stdout", "logunit", "unit_log"):
                scope.assigns[var].append((si, "__write__:" + what, None))

        masked = mask_strings(exec_body)

        # ---- id references & sends
        for msd in SEND_RE.finditer(masked):
            o = msd.end() - 1
            inner, _ = paren_group(exec_body, o)
            first = split_top(inner)[0].strip()
            mmid = re.match(r"^([\w%]+)", first)
            if mmid:
                base = mmid.group(1).split("%")[-1].lower()
                pf.sends.add(base)
        for tok in re.findall(r"[A-Za-z_]\w*", masked):
            pf.idrefs[tok.lower()] += 1

        # ---- registrations
        regs_here = list(REG_RE.finditer(masked))
        if regs_here and scope.module_scope() is not None and \
                (scope.module_scope().name in SKIP_MODULES):
            regs_here = []
        for mr in regs_here:
            rc = RegCall()
            rc.kind = mr.group(1).lower()
            o = mr.end() - 1
            inner, e = paren_group(exec_body, o)
            rc.raw = norm_ws(exec_body[mr.start():e + 1])
            rc.file = relpath
            rc.line = st.line_at(len(text) - len(exec_body) + mr.start())
            rc.scope = scope
            rc.stmt_idx = si
            rc.conds = conds
            rc.loops = loops_now
            ms = scope.module_scope()
            rc.fortran_module = ms.name if ms is not None and ms.kind != "file" else ""
            ps = scope.proc_scope()
            rc.proc = ps.name if ps is not None else ""
            # id variable
            pre = masked[:mr.start()]
            ma = re.match(r"^\s*([\w%]+)\s*(\([^=]*\))?\s*=\s*$", pre)
            if ma:
                rc.id_var = ma.group(1).split("%")[-1].lower()
                rc.id_expr = norm_ws(exec_body[:mr.start()].rstrip().rstrip("="))
            else:
                rc.id_expr = ""
            rc.args = classify_args(rc.kind, inner)
            pf.regs.append(rc)

        # ---- simple assignments (after registrations so a statement that is
        #      itself "x = register..." does not shadow)
        if not regs_here:
            ma = ASSIGN_RE.match(exec_body)
            if ma and not KEYWORD_STMT.match(exec_body):
                nm = ma.group(1).lower()
                sub = ma.group(2)
                rhs = exec_body[ma.end():].strip()
                scope.assigns[nm].append((si, rhs, sub.strip() if sub else None))
    return pf


def classify_args(kind, inner):
    args = [a.strip() for a in split_top(inner)]
    out = OrderedDict()
    positional = []
    for a in args:
        mk = re.match(r"^(\w+)\s*=(?![=>])\s*(.*)$", a, re.S)
        if mk and not a.startswith("(") and "//" not in mk.group(1):
            out[mk.group(1).lower()] = mk.group(2).strip()
        else:
            positional.append(a)
    if kind == "static":
        names = STATIC_ARGS
    else:
        names = DIAG_ARRAY_ARGS
        third = positional[2] if len(positional) > 2 else ""
        if "axes" not in out and ("init_time" in out and len(positional) <= 2):
            names = DIAG_SCALAR_ARGS
        elif third and re.match(r"^\w*time\w*$", third, re.I) and \
                not re.search(r"ax", third, re.I):
            names = DIAG_SCALAR_ARGS
    for k, a in enumerate(positional):
        if k < len(names):
            out.setdefault(names[k], a)
        else:
            out.setdefault("extra%d" % k, a)
    if names is DIAG_SCALAR_ARGS:
        out["axes"] = "(scalar)"
    return out


# --------------------------------------------------------------------------
# Expression evaluation
# --------------------------------------------------------------------------

AXIS_INDEX_NAMES = {1: "lon", 2: "lat", 3: "pfull", 4: "phalf"}


class Evaluator(object):
    def __init__(self, pfile, modtab, max_depth=8):
        self.pf = pfile
        self.modtab = modtab       # module name -> Scope (global)
        self.max_depth = max_depth

    # -- symbol lookup through scope chain and USE'd modules
    def lookup(self, name, scope):
        s = scope
        seen = set()
        while s is not None:
            if name in s.symbols:
                return s.symbols[name], s
            s = s.parent
        s = scope
        while s is not None:
            for u in s.uses:
                if u in seen:
                    continue
                seen.add(u)
                ms = self.modtab.get(u)
                if ms is not None and name in ms.symbols and ms.symbols[name].is_param:
                    return ms.symbols[name], ms
            s = s.parent
        return None, None

    def last_assignment(self, name, scope, stmt_idx):
        s = scope
        while s is not None:
            lst = s.assigns.get(name)
            if lst:
                prev = [a for a in lst if a[0] < stmt_idx and a[2] is None]
                if prev:
                    return prev[-1], s
            if s.kind in ("subroutine", "function"):
                break
            s = s.parent
        return None, None

    def tracer_role(self, name, scope):
        s = scope
        while s is not None:
            if name in s.tracer_vars:
                return s.tracer_vars[name]
            s = s.parent
        return None

    # -- string expression -> list of (text, resolved) pieces
    def eval_str(self, expr, scope, stmt_idx, env, depth=0, info=None):
        """Returns (string, fully_resolved).  Unresolved sub-expressions are
        rendered as {expr} (or <tracer> for tracer names)."""
        if info is None:
            info = {"expand": {}, "sources": []}
        expr = expr.strip()
        if depth > self.max_depth:
            return "{%s}" % norm_ws(expr), False
        terms = split_top(expr, "//")
        out = []
        ok = True
        for t in terms:
            s, r = self.eval_term(t.strip(), scope, stmt_idx, env, depth, info)
            out.append(s)
            ok = ok and r
        return "".join(out), ok

    def eval_term(self, t, scope, stmt_idx, env, depth, info):
        if not t:
            return "", True
        # strip redundant parens
        if t.startswith("(") and not t.startswith("(/"):
            e = match_paren(mask_strings(t), 0)
            if e == len(t) - 1:
                return self.eval_str(t[1:-1], scope, stmt_idx, env, depth + 1, info)
        mq = re.match(r"^(?:\w+_)?('([^']|'')*'|\"([^\"]|\"\")*\")$", t, re.S)
        if mq:
            lit = mq.group(1)
            return unquote(lit), True
        mf = re.match(r"^(trim|adjustl|adjustr|lowercase|uppercase|lcase|ucase)\s*\(", t, re.I)
        if mf:
            o = mf.end() - 1
            inner, e = paren_group(t, o)
            if e == len(t) - 1:
                s, r = self.eval_str(inner, scope, stmt_idx, env, depth + 1, info)
                fn = mf.group(1).lower()
                if r:
                    if fn == "trim":
                        s = s.rstrip()
                    elif fn == "adjustl":
                        s = s.lstrip() + " " * (len(s) - len(s.lstrip()))
                    elif fn in ("lowercase", "lcase"):
                        s = s.lower()
                    elif fn in ("uppercase", "ucase"):
                        s = s.upper()
                return s, r
        # plain identifier
        mi = re.match(r"^(\w+)$", t)
        if mi:
            name = mi.group(1).lower()
            if name in env and isinstance(env[name], str):
                return env[name], True
            return self.resolve_name(name, scope, stmt_idx, env, depth, info)
        # array element  name(idx)
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
                    info["expand"][idx] = (name, len(elems))
                return "{%s}" % norm_ws(t), False
        return "{%s}" % norm_ws(t), False

    def resolve_name(self, name, scope, stmt_idx, env, depth, info):
        role = self.tracer_role(name, scope)
        sym, sscope = self.lookup(name, scope)
        # local assignment wins for non-parameters
        if sym is None or not sym.is_param:
            asg, ascope = self.last_assignment(name, scope, stmt_idx)
            if asg is not None and asg[1].startswith("__write__:"):
                return "{%s}" % asg[1][len("__write__:"):], False
            if asg is not None:
                s, r = self.eval_str(asg[1], ascope, asg[0], env, depth + 1, info)
                if r:
                    info["sources"].append("%s (assigned line %d)" % (name, self.pf.stmts[asg[0]].line))
                    return self.apply_len(s, sym, info, name), True
                return s, False
        if role:
            return "<%s>" % role, False
        if sym is not None and sym.expr is not None and not sym.expr.startswith("(/") \
                and not sym.expr.startswith("["):
            s, r = self.eval_str(sym.expr, sscope, 0, env, depth + 1, info)
            if r:
                info["sources"].append("%s=%s%s" % (
                    name, "parameter " if sym.is_param else "initialised ", sscope.name))
                return self.apply_len(s, sym, info, name), True
        return "{%s}" % name, False

    def apply_len(self, s, sym, info, name):
        if sym is None or sym.charlen is None:
            return s
        try:
            n = int(sym.charlen)
        except ValueError:
            return s
        if len(s.rstrip()) > n:
            info.setdefault("notes", []).append(
                "TRUNCATED: %s is character(len=%d) but value '%s' is longer" % (name, n, s))
            return s[:n]
        return s

    def literal_array(self, name, scope, stmt_idx, depth):
        sym, sscope = self.lookup(name, scope)
        exprs = []
        if sym is not None and sym.expr and (sym.expr.startswith("(/") or sym.expr.startswith("[")):
            exprs = [sym.expr]
        if not exprs:
            asg, ascope = self.last_assignment(name, scope, stmt_idx)
            if asg is not None and (asg[1].startswith("(/") or asg[1].startswith("[")):
                exprs = [asg[1]]
                sscope = ascope
        if not exprs:
            return None
        e = exprs[0].strip()
        inner = e[2:-2] if e.startswith("(/") else e[1:-1]
        vals = []
        for item in split_top(inner):
            s, r = self.eval_str(item, sscope, 0, {}, depth + 1)
            if not r:
                return None
            vals.append(self.apply_len(s, sym, {}, name))
        return vals

    def eval_int(self, expr, scope, env):
        expr = expr.strip().lower()
        if re.match(r"^[+-]?\d+$", expr):
            return int(expr)
        if expr in env and isinstance(env[expr], int):
            return env[expr]
        ms = re.match(r"^size\s*\(\s*(\w+)\s*(\(\s*:\s*\))?\s*\)$", expr)
        if ms:
            arr = self.literal_array(ms.group(1), scope, 10 ** 9, 0)
            if arr is not None:
                return len(arr)
            return None
        if re.match(r"^\w+$", expr):
            sym, sscope = self.lookup(expr, scope)
            if sym is not None and sym.is_param and sym.expr:
                return self.eval_int(sym.expr, sscope, env)
        # simple a+b / a-b
        mo = re.match(r"^(\w+)\s*([+-])\s*(\w+)$", expr)
        if mo:
            a = self.eval_int(mo.group(1), scope, env)
            b = self.eval_int(mo.group(3), scope, env)
            if a is not None and b is not None:
                return a + b if mo.group(2) == "+" else a - b
        return None

    def dims(self, axes, scope, stmt_idx):
        """Axis names, e.g. 'lon, lat, pfull'.  MiMA passes the axes of the
        atmosphere as an array (lon, lat, pfull, phalf); id_<name> variables
        are named axes.  Returns '' when the axes cannot be worked out."""
        a = norm_ws(axes or "")
        if not a:
            return ""
        if a == "(scalar)":
            return "scalar"

        def elements(expr):
            e = norm_ws(expr)
            if e.startswith("(/") and e.endswith("/)"):
                return split_top(e[2:-2])
            if e.startswith("[") and e.endswith("]"):
                return split_top(e[1:-1])
            return None

        def named(items):
            out = []
            for it in items:
                it = it.strip()
                m = re.match(r"^\w+\s*\(\s*(\d+)\s*\)$", it)
                if m:
                    out.append(AXIS_INDEX_NAMES.get(int(m.group(1))))
                    continue
                m = re.match(r"^id_(\w+)$", it, re.I)
                if m:
                    out.append(m.group(1).lower())
                    continue
                return None
            return None if None in out else ", ".join(out)

        def indexed(items):
            try:
                return ", ".join(AXIS_INDEX_NAMES[int(i)] for i in items)
            except (KeyError, ValueError):
                return None

        def value_of(name):
            sym, sscope = self.lookup(name, scope)
            if sym is not None and sym.expr and elements(sym.expr) is not None:
                return elements(sym.expr), sym
            asg, _ = self.last_assignment(name, scope, stmt_idx)
            if asg is not None and elements(asg[1]) is not None:
                return elements(asg[1]), sym
            return None, sym

        m = re.match(r"^\w+\s*\(\s*(\d+)\s*:\s*(\d+)\s*\)$", a)
        if m:
            return indexed(range(int(m.group(1)), int(m.group(2)) + 1)) or ""
        if elements(a) is not None:
            return named(elements(a)) or ""
        m = re.match(r"^\w+\s*\(\s*(\w+)\s*\)$", a)
        if m:
            items, _ = value_of(m.group(1).lower())
            return (indexed(items) if items else None) or ""
        m = re.match(r"^(\w+)$", a)
        if m:
            items, sym = value_of(m.group(1).lower())
            if items:
                return named(items) or ""
            if sym is not None and sym.dims and re.match(r"^\d+$", sym.dims.strip()):
                return indexed(range(1, int(sym.dims.strip()) + 1)) or ""
        return ""

    def ndim(self, axes, scope):
        if axes is None:
            return ""
        a = norm_ws(axes)
        if a == "(scalar)":
            return "0"
        m = re.match(r"^\w+\s*\(\s*(\d+)\s*:\s*(\d+)\s*\)$", a)
        if m:
            return str(int(m.group(2)) - int(m.group(1)) + 1)
        m = re.match(r"^\w+\s*\(\s*(\w+)\s*\)$", a)
        if m:
            arr = m.group(1).lower()
            sym, sscope = self.lookup(arr, scope)
            if sym is not None and sym.expr and sym.expr.startswith("(/"):
                return str(len(split_top(sym.expr[2:-2])))
            if sym is not None and sym.dims and re.match(r"^\d+$", sym.dims.strip()):
                return sym.dims.strip()
            return ""
        if a.startswith("(/") or a.startswith("["):
            inner = a[2:-2] if a.startswith("(/") else a[1:-1]
            return str(len(split_top(inner)))
        m = re.match(r"^(\w+)$", a)
        if m:
            sym, _ = self.lookup(m.group(1).lower(), scope)
            if sym is not None and sym.dims:
                d = sym.dims.strip()
                if re.match(r"^\d+$", d):
                    return d
                m2 = re.match(r"^(\d+)\s*:\s*(\d+)$", d)
                if m2:
                    return str(int(m2.group(2)) - int(m2.group(1)) + 1)
        return ""


# --------------------------------------------------------------------------
# CMake: which sources are compiled
# --------------------------------------------------------------------------

def built_sources(root):
    built = set()
    found_any = False
    for dirpath, dirnames, filenames in os.walk(root):
        if os.sep + "build" in dirpath[len(root):] or "/.git" in dirpath:
            continue
        if "CMakeLists.txt" not in filenames:
            continue
        found_any = True
        p = os.path.join(dirpath, "CMakeLists.txt")
        with io.open(p, encoding="utf-8", errors="replace") as fh:
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


# --------------------------------------------------------------------------
# Main inventory
# --------------------------------------------------------------------------

def find_sources(src):
    out = []
    for dirpath, dirnames, filenames in os.walk(src):
        dirnames.sort()
        for f in sorted(filenames):
            if f.endswith((".f90", ".F90")):
                out.append(os.path.join(dirpath, f))
    return out


def physics_dir(relpath):
    parts = relpath.replace("\\", "/").split("/")
    if parts and parts[0] == "src":
        parts = parts[1:]
    parts = parts[:-1]
    if parts and parts[-1] == "null":
        parts = parts[:-1]
    return "/".join(parts)


def is_legacy(relpath, legacy_dirs):
    p = "/" + relpath.replace("\\", "/")
    for d in legacy_dirs:
        if "/" + d.strip("/") + "/" in p:
            return True
    return False


def make_clause(text, nml_map):
    """(readable text, ((var, namelist group), ...)) for the namelist
    variables that appear in the clause."""
    n = cond_parse(text)
    names = cond_vars(n) if n is not None else re.findall(r"[a-z_]\w*", text.lower())
    return (text, tuple((v, nml_map[v]) for v in names if v in nml_map))


def build_inventory(root, src, legacy_dirs, max_call_depth=4):
    files = find_sources(src)
    parsed = []
    for f in files:
        rel = os.path.relpath(f, root)
        parsed.append(parse_file(f, rel))
    modtab = {}
    mod_file = {}
    for pf in parsed:
        for s in pf.scopes:
            if s.kind == "module":
                modtab.setdefault(s.name, s)
                mod_file.setdefault(s.name, pf)
    built = built_sources(root)

    # global send set, and caller index
    global_sends = defaultdict(set)
    for pf in parsed:
        for s in pf.sends:
            global_sends[s].add(pf.relpath)
    callers = defaultdict(list)   # callee -> [(pf, scope, cond, line)]
    for pf in parsed:
        for s in pf.scopes:
            for callee, si, conds, line in s.calls:
                callers[callee].append((pf, s, conds, line))

    def file_uses(pf, modname, name=None):
        for s in pf.scopes:
            if modname in s.uses:
                only = s.use_only.get(modname, None)
                if name is None or only is None or name in only:
                    return True
        return False

    def call_gates(pf, proc_scope, depth, seen):
        """Walk the call chain upward.  Returns a list of chains; each chain is
        a list of (where, cond, nml) from the innermost caller outward.  Only
        call sites with a non-trivial condition are kept in a chain."""
        if proc_scope is None or depth > max_call_depth:
            return []
        name = proc_scope.name
        if name in seen:
            return []
        seen = seen | {name}
        ms = proc_scope.module_scope()
        defmod = ms.name if ms is not None else ""
        chains = []
        for cpf, cscope, conds, line in callers.get(name, []):
            if cpf is not pf:
                # a different file: must USE the defining module (and the
                # name must not be one of its own procedures)
                if name in cpf.proc_defs or not file_uses(cpf, defmod, name):
                    continue
            where = "%s:%d (%s)" % (os.path.basename(cpf.relpath), line, cscope.name)
            csc = cscope.proc_scope()
            guard = guard_list(csc.guards) if csc is not None else []
            cc = [c for c in conds if c and not TRIVIAL_COND_RE.match(c)] + guard
            nml = sorted({"%s:%s" % (cpf.namelists[t], t)
                          for c in cc for t in re.findall(r"[a-z_]\w*", c.lower())
                          if t in cpf.namelists})
            link = [(where, " .and. ".join(cc), nml)] if cc else []
            up = call_gates(cpf, csc, depth + 1, seen) if csc is not None else []
            if up:
                for u in up:
                    chains.append(link + u)
            else:
                chains.append(link)
        # An empty chain is a call path without any condition: the routine
        # is then called unconditionally and the other paths do not matter.
        out = []
        for ch in chains:
            if not ch:
                return []
            if ch not in out:
                out.append(ch)
        return out[:4]

    rows = []
    unresolved_rows = []
    for pf in parsed:
        ev = Evaluator(pf, modtab)
        for rc in pf.regs:
            scope = rc.scope
            psc = scope.proc_scope()
            # guards: "if (x) return" earlier in same procedure
            guard_conds = []
            if psc is not None:
                guard_conds = guard_list(psc.guards, rc.stmt_idx)
            # determine expansion variables with a first pass
            info0 = {"expand": {}, "sources": []}
            for key in ("module_name", "field_name", "long_name", "units"):
                if key in rc.args:
                    ev.eval_str(rc.args[key], scope, rc.stmt_idx, {}, info=info0)
            envs = [{}]
            loopvars = {}
            for lh in rc.loops:
                mlh = re.match(r"^(\w+)\s*=\s*(.+)$", lh)
                if mlh:
                    bounds = split_top(mlh.group(2))
                    loopvars[mlh.group(1).lower()] = bounds
            for var, (arr, n) in info0["expand"].items():
                lo, hi = 1, n
                if var in loopvars:
                    b = loopvars[var]
                    l0 = ev.eval_int(b[0], scope, {}) if len(b) > 0 else None
                    h0 = ev.eval_int(b[1], scope, {}) if len(b) > 1 else None
                    lo = l0 if l0 is not None else 1
                    hi = h0 if h0 is not None else n
                    hi = min(hi, n)
                new = []
                for e in envs:
                    for k in range(lo, hi + 1):
                        e2 = dict(e)
                        e2[var] = k
                        new.append(e2)
                envs = new
            for env in envs:
                info = {"expand": {}, "sources": [], "notes": []}
                vals = {}
                res = {}
                for key in ("module_name", "field_name", "long_name", "units",
                            "standard_name"):
                    if key in rc.args:
                        v, r = ev.eval_str(rc.args[key], scope, rc.stmt_idx, env, info=info)
                        vals[key] = v.strip() if r else v
                        res[key] = r
                    else:
                        vals[key] = ""
                        res[key] = True
                conds = [c for c in rc.conds if c] + guard_conds
                cond_text = " .and. ".join(conds)
                gates = call_gates(pf, psc, 0, set())
                gate_text = " | ".join(" <- ".join("%s: %s" % (w, c) for (w, c, n) in ch)
                                       for ch in gates)
                nml = set()
                for c in conds:
                    for t in re.findall(r"[a-z_]\w*", c.lower()):
                        if t in pf.namelists:
                            nml.add("%s:%s" % (pf.namelists[t], t))
                gnml = set()
                for ch in gates:
                    for (w, c, n) in ch:
                        nml.update(n)
                        gnml.update(n)
                # send status
                idv = rc.id_var
                if not idv:
                    send = "no-id-var"
                elif idv in pf.sends:
                    send = "sent"
                elif idv in global_sends and any(
                        file_uses(p2, rc.fortran_module) for p2 in parsed
                        if p2.relpath in global_sends[idv]):
                    send = "sent-elsewhere(%s)" % ",".join(
                        sorted(os.path.basename(p2.relpath) for p2 in parsed
                               if p2.relpath in global_sends[idv]
                               and file_uses(p2, rc.fortran_module)))
                else:
                    # how many other references in the file?
                    nref = pf.idrefs.get(idv, 0)
                    nreg = sum(1 for r2 in pf.regs if r2.id_var == idv)
                    send = "referenced-not-sent" if nref > nreg else "never-sent"
                unresolved = [k for k in ("module_name", "field_name") if not res.get(k, True)]
                dynamic = bool(unresolved) or bool(env)
                notes = list(info.get("notes", []))
                if "#if defined(test" in cond_text or (gate_text and all(
                        "#if defined(test" in " ".join(c for (w, c, n) in ch) for ch in gates)):
                    notes.append("test-program only (not reachable in the model)")
                if env:
                    notes.append("expanded loop: " + ", ".join("%s=%d" % kv for kv in sorted(env.items())))
                for key in ("long_name", "units"):
                    if not res.get(key, True):
                        notes.append("%s not resolved" % key)
                built_flag = ""
                if built is not None:
                    built_flag = "yes" if pf.relpath in built else "no"
                mv = rc.args.get("missing_value", "")
                row = OrderedDict()
                row["module"] = vals["module_name"]
                row["field"] = vals["field_name"]
                row["kind"] = rc.kind
                row["axes"] = norm_ws(rc.args.get("axes", ""))
                row["ndim"] = ev.ndim(rc.args.get("axes"), scope)
                row["units"] = vals["units"]
                row["long_name"] = norm_ws(vals["long_name"])
                row["standard_name"] = vals["standard_name"]
                row["missing_value"] = norm_ws(mv)
                row["range"] = norm_ws(rc.args.get("range", ""))
                row["id_var"] = rc.id_expr if hasattr(rc, "id_expr") else ""
                row["file"] = pf.relpath.replace("\\", "/")
                row["line"] = rc.line
                row["subroutine"] = rc.proc
                row["fortran_module"] = rc.fortran_module
                row["condition"] = cond_text
                row["call_gates"] = gate_text
                row["namelist_vars"] = " ".join(sorted(nml))
                row["loops"] = " ; ".join(rc.loops)
                row["dynamic"] = "yes" if dynamic else ""
                row["unresolved"] = ",".join(unresolved)
                row["send_status"] = send
                row["built"] = built_flag
                row["legacy"] = "yes" if is_legacy(pf.relpath, legacy_dirs) else ""
                row["physics_dir"] = physics_dir(pf.relpath)
                row["notes"] = "; ".join(notes)
                alts = []
                for ch in gates:
                    toks = []
                    for (w, c, n) in ch:
                        grps = sorted({x.split(":")[0] for x in n})
                        tok = ("%s: %s" % (grps[0], c)) if len(grps) == 1 else c
                        if not toks or toks[-1] != tok:
                            toks.append(tok)
                    t = " <- ".join(toks)
                    if t not in alts:
                        alts.append(t)
                row["gate_nml"] = " | ".join(alts)
                # structured availability: alternatives (one per call path),
                # each a tuple of readable clauses (text, ((var, group), ...))
                local = []
                for c in conds:
                    for cl in readable_clauses(c):
                        local.append(make_clause(cl, pf.namelists))
                avail = []
                for ch in gates:
                    alt = list(local)
                    for (w, c, n) in ch:
                        nmap = dict((x.split(":")[1], x.split(":")[0]) for x in n)
                        for cl in readable_clauses(c):
                            k = make_clause(cl, nmap)
                            if k not in alt:
                                alt.append(k)
                    if tuple(alt) not in avail:
                        avail.append(tuple(alt))
                if not gates:
                    avail = [tuple(local)]
                row["_alts"] = avail
                row["available"] = alts_label(avail, md=False, sep=" | ")
                d = ev.dims(rc.args.get("axes"), scope, rc.stmt_idx)
                row["dims"] = d if d else ("%sD" % row["ndim"] if row["ndim"] else "")
                rows.append(row)
    rows.sort(key=lambda r: (r["module"].lower(), r["field"].lower(), r["file"], r["line"]))
    return rows, parsed


# --------------------------------------------------------------------------
# Conditions: Fortran logical expressions -> readable text, and evaluation
# --------------------------------------------------------------------------
#
# Conditions are kept as Fortran text while the sources are parsed.  For the
# documentation they are parsed into a small tree, simplified and printed as
#     trim(radiation_scheme) == ('rrtm')   ->  radiation_scheme = 'rrtm'
#     (.not. (do_bm))                      ->  not do_bm
#     x == ('a', 'b')    (select case)     ->  x in ('a', 'b')
# The printed form can be parsed again (and/or/not, =, /=, in are accepted),
# which is what --nml uses to evaluate a condition against an input.nml.

_DOP_RE = r"\.(?:and|or|not|eqv|neqv|eq|ne|lt|le|gt|ge|true|false)\."
COND_TOKEN_RE = re.compile(
    r"\s*(?:(?P<str>'(?:[^']|'')*'|\"(?:[^\"]|\"\")*\")"
    r"|(?P<dop>" + _DOP_RE + r")"
    r"|(?P<num>\d+(?:\.(?!(?:and|or|not|eqv|neqv|eq|ne|lt|le|gt|ge)\.)\d*)?(?:[eEdD][+-]?\d+)?)"
    r"|(?P<name>[A-Za-z_][\w%]*)"
    r"|(?P<op>==|/=|>=|<=|=|<|>|\(|\)|,))", re.I)

_CMP_OPS = {"==": "=", "=": "=", ".eq.": "=", "/=": "/=", ".ne.": "/=",
            ".eqv.": "=", ".neqv.": "/=", "<": "<", ".lt.": "<", "<=": "<=",
            ".le.": "<=", ">": ">", ".gt.": ">", ">=": ">=", ".ge.": ">="}
_CMP_NEG = {"=": "/=", "/=": "=", "<": ">=", ">=": "<", ">": "<=", "<=": ">"}
_TRANSPARENT_FUNCS = {"trim", "adjustl"}
# a clause that selects on the value of a character variable
SELECTOR_RE = re.compile(r"^(\w+) = '([^']*)'$")


class CondParseError(Exception):
    pass


def _cond_tokens(text):
    toks = []
    pos = 0
    text = text.rstrip()
    while pos < len(text):
        m = COND_TOKEN_RE.match(text, pos)
        if not m or m.end() == pos:
            raise CondParseError(text)
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


class _CondParser(object):
    def __init__(self, text):
        self.text = norm_ws(text)
        self.toks = _cond_tokens(self.text)
        self.i = 0

    def peek(self, k=0):
        j = self.i + k
        return self.toks[j] if j < len(self.toks) else (None, None)

    def take(self):
        t = self.peek()
        self.i += 1
        return t

    def expect(self, val):
        t = self.take()
        if t[1] != val:
            raise CondParseError(self.text)

    def parse(self):
        # 'X == case default' is how the parser records the default branch
        m = re.match(r"^(.*?)\s*==\s*case default$", self.text, re.I)
        if m:
            return ("default", _CondParser(m.group(1)).parse())
        n = self.p_or()
        if self.i != len(self.toks):
            raise CondParseError(self.text)
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
            vals = right[1] if right[0] == "tuple" else [right]
            return ("in", left, vals)
        if k == "op" and v in _CMP_OPS:
            self.take()
            right = self.p_primary()
            if right[0] == "tuple":
                return ("in", left, right[1])
            return ("cmp", _CMP_OPS[v], left, right)
        if left[0] == "tuple":
            raise CondParseError(self.text)
        return left

    def p_primary(self):
        k, v = self.take()
        if k is None:
            raise CondParseError(self.text)
        if v == "(" and k == "op":
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
            q = v[0]
            s = v[1:-1].replace(q + q, q)
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
            if self.peek() == ("op", "("):
                # function call or array element: keep the argument text
                start = self.i
                depth = 0
                while True:
                    tk, tv = self.take()
                    if tk is None:
                        raise CondParseError(self.text)
                    if tk == "op" and tv == "(":
                        depth += 1
                    elif tk == "op" and tv == ")":
                        depth -= 1
                        if depth == 0:
                            break
                inner = self.toks[start + 1:self.i - 1]
                if name in _TRANSPARENT_FUNCS and len(inner) == 1 and inner[0][0] == "name":
                    return ("var", inner[0][1].lower())
                return ("call", name, _join_tokens(inner))
            return ("var", name)
        raise CondParseError(self.text)


def _join_tokens(toks):
    out = ""
    for k, v in toks:
        if v == ",":
            out += ", "
        elif k == "kw":
            out += " %s " % v
        elif k == "op" and v not in ("(", ")"):
            out += " %s " % v
        else:
            out += v
    return norm_ws(out)


def cond_parse(text):
    """Parse a condition; returns a tree, or None if it cannot be parsed."""
    try:
        return _CondParser(text).parse()
    except CondParseError:
        return None


def cond_show(n, prec=0):
    k = n[0]
    if k == "or":
        s = " or ".join(cond_show(c, 1) for c in n[1])
        return "(%s)" % s if prec > 1 else s
    if k == "and":
        s = " and ".join(cond_show(c, 2) for c in n[1])
        return "(%s)" % s if prec > 2 else s
    if k == "not":
        inner = n[1]
        if inner[0] == "not":
            return cond_show(inner[1], prec)
        if inner[0] == "cmp":
            return cond_show(("cmp", _CMP_NEG[inner[1]], inner[2], inner[3]), prec)
        if inner[0] == "in":
            return "%s not in (%s)" % (cond_show(inner[1], 5),
                                       ", ".join(cond_show(v, 5) for v in inner[2]))
        return "not " + cond_show(inner, 3)
    if k == "cmp":
        return "%s %s %s" % (cond_show(n[2], 5), n[1], cond_show(n[3], 5))
    if k == "in":
        if len(n[2]) == 1:
            return "%s = %s" % (cond_show(n[1], 5), cond_show(n[2][0], 5))
        return "%s in (%s)" % (cond_show(n[1], 5), ", ".join(cond_show(v, 5) for v in n[2]))
    if k == "lit":
        return n[2]
    if k == "var":
        return n[1]
    if k == "call":
        return "%s(%s)" % (n[1], n[2])
    if k == "default":
        return "%s = (any value not listed in the select case)" % cond_show(n[1], 5)
    return "?"


def cond_vars(n, out=None):
    if out is None:
        out = []
    k = n[0]
    if k == "var":
        if n[1] not in out:
            out.append(n[1])
    elif k in ("or", "and", "tuple"):
        for c in n[1]:
            cond_vars(c, out)
    elif k in ("not", "default"):
        cond_vars(n[1], out)
    elif k == "cmp":
        cond_vars(n[2], out)
        cond_vars(n[3], out)
    elif k == "in":
        cond_vars(n[1], out)
        for c in n[2]:
            cond_vars(c, out)
    elif k == "call":
        for t in re.findall(r"[a-z_]\w*", n[2].lower()):
            if t not in out:
                out.append(t)
    return out


def readable_clauses(text):
    """Split a Fortran condition into readable top-level '.and.' clauses."""
    text = norm_ws(text or "")
    if not text:
        return []
    n = cond_parse(text)
    if n is None:
        return [text]
    items = n[1] if n[0] == "and" else [n]
    out = []
    for it in items:
        s = cond_show(it)
        if s not in out:
            out.append(s)
    return out


def cond_eval(n, value_of):
    """Three-valued evaluation: True, False or None (unknown).  value_of(name)
    returns a Python bool/str/float, or None when the value is unknown."""
    k = n[0]
    if k == "lit":
        return n[1]
    if k == "var":
        return value_of(n[1])
    if k in ("call", "default", "tuple"):
        return None
    if k == "not":
        v = cond_eval(n[1], value_of)
        return (not v) if isinstance(v, bool) else None
    if k in ("and", "or"):
        vals = [cond_eval(c, value_of) for c in n[1]]
        vals = [v if isinstance(v, bool) else None for v in vals]
        if k == "and":
            if False in vals:
                return False
            return None if None in vals else True
        if True in vals:
            return True
        return None if None in vals else False
    if k in ("cmp", "in"):
        a = cond_eval(n[2] if k == "cmp" else n[1], value_of)
        bs = [cond_eval(n[3], value_of)] if k == "cmp" else [cond_eval(v, value_of) for v in n[2]]
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


def fortran_value(text):
    """Python value of a namelist/initialiser constant, or None."""
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


def namelist_defaults(parsed):
    """{(group, var): default value} from the declarations of namelist variables."""
    out = {}
    for pf in parsed:
        for var, grp in pf.namelists.items():
            for s in pf.scopes:
                sym = s.symbols.get(var)
                if sym is not None and s.kind in ("module", "program", "file"):
                    out[(grp.lower(), var)] = fortran_value(sym.expr) if sym.expr else None
                    break
    return out


def clause_label(clause, md=True):
    """Render a clause (text, ((var, group), ...)) with its namelist group."""
    text, vg = clause
    groups = []
    for v, g in vg:
        g = re.sub(r"_nml$", "", g)
        if g not in groups:
            groups.append(g)
    s = md_code(text) if md else text
    if groups:
        s += " (%s)" % ", ".join(groups)
    return s


def alt_label(alt, md=True):
    return " and ".join(clause_label(c, md) for c in alt)


def alts_label(alts, md=True, sep=" or "):
    if not alts or any(len(a) == 0 for a in alts):
        return ""
    # x = 'a' or x = 'b'  ->  x in ('a', 'b')
    ms = [SELECTOR_RE.match(a[0][0]) if len(a) == 1 else None for a in alts]
    if len(alts) > 1 and all(ms) and len({(m.group(1), a[0][1]) for m, a in zip(ms, alts)}) == 1:
        vals = sorted({m.group(2) for m in ms})
        text = "%s in (%s)" % (ms[0].group(1), ", ".join("'%s'" % v for v in vals))
        return clause_label((text, alts[0][0][1]), md)
    return sep.join(alt_label(a, md) for a in alts)


def alts_eval(alts, value_of):
    """True if some alternative holds, False if none can, None if unknown."""
    if not alts or any(len(a) == 0 for a in alts):
        return True
    res = []
    for a in alts:
        vals = []
        for text, vg in a:
            n = cond_parse(text)
            groups = dict(vg)
            vals.append(cond_eval(n, lambda name, g=groups: value_of(g.get(name), name))
                        if n is not None else None)
        if False in vals:
            res.append(False)
        elif None in vals:
            res.append(None)
        else:
            res.append(True)
    if True in res:
        return True
    return None if None in res else False


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------

def md_escape(s):
    s = str(s)
    s = s.replace("|", "\\|").replace("\n", " ")
    s = s.replace("<", "&lt;").replace(">", "&gt;")
    return s


def md_code(s):
    s = str(s)
    if not s:
        return ""
    return "`" + s.replace("`", "'").replace("|", "\\|") + "`"


def write_csv(rows, path):
    with io.open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_COLUMNS, lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def gate_nml(r):
    """namelist switches found on the call chain (not the local condition)."""
    return r.get("gate_nml", "")


HOWTO_MD = r'''
MiMA writes only the diagnostics you ask for in the `diag_table` file in the run
directory. This page explains how to write a `diag_table` and lists every field the
model can provide, with the namelist settings it needs.

* [How to write a diag_table](#how-to-write-a-diag_table)
* [Checking a diag_table](#checking-a-diag_table)
* [Reading the field tables](#reading-the-field-tables)
* [Fields by module](#fields-by-module)

Ready-made tables are in the repository: `input/diag_table` (the default RRTM test
case), `input/examples/gray/diag_table` (gray radiation) and
`input/examples/held_suarez/diag_table` (dry Held-Suarez).

## How to write a diag_table

A `diag_table` is a plain text file of comma-separated values. Anything after a `#` is a
comment. It has three kinds of lines: the two global lines, which must be lines 1 and 2
of the file, then file lines and field lines in any order (define a file before the
fields that go into it).

```
"MiMA experiment"
0001 1 1 0 0 0
#  file name,   output frequency, units, format, time units, time axis name
"atmos_daily",    1, "days", 1, "days", "time",
"atmos_monthly", 30, "days", 1, "days", "time",
#  module,     field,    output name, file,           sampling, reduction, region, packing
 "dynamics", "ps",     "ps",        "atmos_daily",   "all",    .true.,    "none", 2,
 "dynamics", "ucomp",  "ucomp",     "atmos_monthly", "all",    .true.,    "none", 2,
 "moist",    "precip", "precip",    "atmos_daily",   "all",    .true.,    "none", 2,
 "moist",    "precip", "precip",    "atmos_daily",   "all",    "max",     "none", 2,
```

This writes daily means of surface pressure and precipitation plus the daily maximum
precipitation (`precip_max`) to `atmos_daily`, and 30-day means of the zonal wind to
`atmos_monthly`.

**Global lines.** Line 1 is a title (written to the file metadata). Line 2 is the base
date, `year month day hour minute second`: the reference time of the time axis. Use the
model start date (`current_date` in `coupler_nml`).

**File lines** define an output file:

| column | example | meaning |
|---|---|---|
| file name | `"atmos_daily"` | Each MPI process writes `atmos_daily.nc.NNNN`. Combine them with `mppnccombine` (see [Getting started](GettingStarted.md#output)). |
| output frequency | `1` | Write every N units (`> 0`), every time step (`0`), or once at the end of the run (`-1`). |
| units | `"days"` | Units of the output frequency: `"seconds"`, `"minutes"`, `"hours"`, `"days"`, `"months"` or `"years"`. |
| format | `1` | Always `1` (netCDF). |
| time units | `"days"` | Units of the time axis in the file. |
| time axis name | `"time"` | Name of the time axis. It must contain `time`. |

With the 30-day calendar of the test case, `30, "days"` and `1, "months"` are the same.

**Field lines** send one diagnostic field to a file:

| column | example | meaning |
|---|---|---|
| module | `"dynamics"` | Module name, from the tables below. It is case sensitive. |
| field | `"ucomp"` | Field name, from the tables below. It is not case sensitive. |
| output name | `"ucomp"` | Name of the variable in the output file. |
| file | `"atmos_daily"` | A file defined by a file line (without `.nc`). |
| sampling | `"all"` | Not used. Always `"all"`. |
| reduction | `.true.` | What to write at each output time. `.true.`: the time mean since the last output. `.false.`: the instantaneous value. `"max"` / `"min"`: the maximum / minimum since the last output (`_max` / `_min` is appended to the output name unless it already ends that way). `"rms"`, `"sum"`, `"pow2"` (mean of the square) and `"diurnal24"` (24 diurnal-cycle means) also work. Use `.false.` for static fields. |
| region | `"none"` | `"none"` writes the whole globe. A region is given as `"lon_min lon_max lat_min lat_max k_min k_max"` in degrees and level indices (`-1 -1` for all levels). |
| packing | `2` | Precision of the output: `1` = double (64-bit), `2` = single (32-bit float). `4` and `8` pack to 16-bit and 8-bit integers. |

A line may be at most 256 characters long. Two fields in the same file must have
different output names; to write the same field twice (for example its mean and its
maximum), give it a different output name or reduction.

**Listing the fields of a run.** The tables below are generated from the source code.
To see which fields are actually registered in a given configuration, add this to
`input.nml`:

```fortran
&diag_manager_nml
    do_diag_field_log = .true. /
```

The root process then writes `diag_field_log.out.0` in the run directory, with one
line per registered field, whether or not it is in `diag_table`:

```
Module|Field|Long Name|Units|Number of Axis|Time Axis|Missing Value|Min Value|Max Value|AXES LIST
dynamics|ucomp|zonal wind component|m/s|3|T||  -400.00000000000000|   400.00000000000000|lon,lat,pfull
```

**Fields that are not registered.** If a field line names a module/field that is not
registered in the run (a typo, or a field whose scheme or switch is off, such as an
RRTM field in a gray run), the model does not stop. When it opens the file, diag_manager
prints a warning, for example

```
WARNING from PE 0: diag_util_mod::opening_file: module/field_name (radiation/tdt_sw) NOT registered
```

and leaves the field out of the file. Errors in the table itself (a wrong number of
columns, an undefined file, a packing value outside 1-8, a duplicated output name) are
fatal.

## Checking a diag_table

`tools/diag_inventory.py` checks a `diag_table` against the fields in the source code
before you run. For example, in a run directory:

```bash
python3 /path/to/MiMA/tools/diag_inventory.py --validate-diag-table diag_table --nml input.nml
```

It reports field lines whose module/field is not registered anywhere (with a hint when
the field exists under another module), fields that are not available with the settings
in `input.nml` (`--nml` is optional; namelist variables not set there take their default
values from the source), and format errors: undefined files, unknown reduction methods,
packing values and duplicated output names. It exits with status 1 if it finds a
problem.

## Reading the field tables

* **field**: `<tracer>` stands for the name of any tracer in `field_table` (for example
  `sphum`). `{expr}` means the name is built at run time from `expr`.
* **dims**: `lon, lat` are the horizontal grid; `pfull` are the full model levels and
  `phalf` the half levels between them (one more than `pfull`); `lat, pfull` is a zonal
  mean. *static* fields have no time axis: they are written once, so use `.false.`
  for them.
* **available with**: the settings for which the field is registered. Namelist
  variables are followed by their namelist in parentheses: `do_bm` (moist_processes)
  is `do_bm` in `&moist_processes_nml`. Variables without a namelist are set inside
  the model (see the source). An empty cell means that the field is always registered
  when its module is active. Conditions that apply to every field of a module are given
  once, above the module's table.
* **source**: where the field is registered (`register_diag_field` or
  `register_static_field`). Fields registered in several places, for example by
  different radiation schemes, have one row with all the places.
* A field whose units, long name or dims differ between the places where it is
  registered is marked **differs**. A field marked *never sent* is registered, but
  the model never passes it any data, so its output contains no model values.
'''


def merge_fields(rs):
    """Group the rows of one module by field name (field names are not case
    sensitive in diag_manager).  Returns an OrderedDict key -> [rows]."""
    out = OrderedDict()
    for r in rs:
        out.setdefault(r["field"].lower(), []).append(r)
    return out


def split_selector(alt):
    """(selector clause match, other clauses) if exactly one clause of alt has
    the form var = 'value', else (None, alt)."""
    ms = [(SELECTOR_RE.match(c[0]), c) for c in alt]
    sel = [(m, c) for m, c in ms if m]
    if len(sel) != 1:
        return None, alt
    return sel[0], tuple(c for c in alt if c is not sel[0][1])


def module_availability(rs):
    """Factor the availability of a module's rows.  Returns (common clauses,
    {id(row): cell alternatives}, selector) where selector is (var, vg) when
    every alternative of every row selects on the value of the same character
    variable, e.g. radiation_scheme = 'gray' / 'rrtm'."""
    all_alts = [a for r in rs for a in r["_alts"]]
    common = []
    if all_alts and all(all_alts):
        for c in all_alts[0]:
            if all(c in a for a in all_alts) and c not in common:
                common.append(c)
    cells = {}
    for r in rs:
        cells[id(r)] = [tuple(c for c in a if c not in common) for a in r["_alts"]]
    keys = set()
    values = set()
    for r in rs:
        for a in cells[id(r)]:
            (sm, extra) = split_selector(a)
            if sm is None:
                return common, cells, None
            keys.add((sm[0].group(1), sm[1][1]))
            values.add(sm[0].group(2))
    if len(keys) == 1 and len(values) > 1:
        return common, cells, keys.pop()
    return common, cells, None


def selector_values(alts):
    vals = []
    for a in alts:
        v = split_selector(a)[0][0].group(2)
        if v not in vals:
            vals.append(v)
    return sorted(vals)


def selector_cell(alts):
    """'gray, rrtm' or, with further conditions, 'rrtm if `x` (grp)'."""
    groups = OrderedDict()
    for a in alts:
        sm, extra = split_selector(a)
        groups.setdefault(extra, [])
        if sm[0].group(2) not in groups[extra]:
            groups[extra].append(sm[0].group(2))
    if () in groups:
        # values without extra conditions make the extra ones redundant
        base = set(groups[()])
        for k in list(groups):
            if k:
                groups[k] = [v for v in groups[k] if v not in base]
                if not groups[k]:
                    del groups[k]
    parts = []
    for extra in sorted(groups, key=lambda e: (len(e), alt_label(e))):
        vals = ", ".join(sorted(groups[extra]))
        parts.append(vals + (" if " + alt_label(extra) if extra else ""))
    return "<br>".join(parts)


def render_markdown(rows, include_legacy=True, include_unbuilt=False, include_test=False,
                    title=None):
    sel = [r for r in rows
           if (include_legacy or not r["legacy"])
           and (include_unbuilt or r["built"] != "no")
           and (include_test or "test-program only" not in r["notes"])]
    out = []
    w = out.append
    w("[back to contents](README.md)")
    w("")
    w("# %s" % (title or "Diagnostics"))
    w("")
    w("<!-- Generated by tools/diag_inventory.py from the Fortran sources: do not edit by hand. -->")
    w("<!-- Regenerate with `python3 tools/diag_inventory.py`; check with `--check`. -->")
    w("")
    w(HOWTO_MD.strip("\n"))
    w("")

    mods = OrderedDict()
    for r in sel:
        mods.setdefault(r["module"], []).append(r)

    def anchor(m):
        return "module-" + re.sub(r"[^a-z0-9_-]", "", m.lower())

    info = {}
    for m, rs in mods.items():
        info[m] = module_availability(rs)
    any_legacy = any(r["legacy"] for r in sel)

    w("## Fields by module")
    w("")
    w("This list is generated from the Fortran sources by `tools/diag_inventory.py`. "
      "After adding or changing a diagnostic, regenerate it with "
      "`python3 tools/diag_inventory.py`; "
      "`python3 tools/diag_inventory.py --check` fails if it is out of date.")
    if not include_unbuilt:
        w("Source files that CMake does not compile are left out.")
    w("")
    hdr = "| module | fields | registered when | source |"
    sep = "|---|---:|---|---|"
    if any_legacy:
        hdr += " legacy |"
        sep += "---|"
    w(hdr)
    w(sep)
    for m, rs in mods.items():
        common, cells, selr = info[m]
        when = alt_label(common) if common else ""
        if selr is not None:
            vals = sorted({v for r in rs for v in selector_values(cells[id(r)])})
            when = (when + " and " if when else "") + "%s (%s): %s" % (
                md_code(selr[0]), ", ".join(re.sub(r"_nml$", "", g) for v, g in selr[1]),
                ", ".join(vals))
        srcs = sorted({os.path.basename(r["file"]) for r in rs})
        line = "| [`%s`](#%s) | %d | %s | %s |" % (
            m, anchor(m), len(merge_fields(rs)), when, ", ".join(srcs))
        if any_legacy:
            line += " %s |" % ("yes" if all(r["legacy"] for r in rs)
                               else ("partly" if any(r["legacy"] for r in rs) else ""))
        w(line)
    w("")

    differs = []
    for m, rs in mods.items():
        common, cells, selr = info[m]
        w('<a id="%s"></a>' % anchor(m))
        w("")
        w("### `%s`" % m)
        w("")
        if common:
            w("Registered only when %s." % alt_label(common))
            w("")
        if selr is not None:
            w("Which fields exist depends on %s in `&%s`: **available with** lists the "
              "values for which the field is registered." % (
                  md_code(selr[0]), ", ".join(g for v, g in selr[1])))
            w("")
        w("| field | dims | units | long_name | available with | source |")
        w("|---|---|---|---|---|---|")
        for key, frs in merge_fields(rs).items():
            alts = []
            for r in frs:
                for a in cells[id(r)]:
                    if a not in alts:
                        alts.append(a)
            if any(len(a) == 0 for a in alts):
                avail = ""
            elif selr is not None:
                avail = selector_cell(alts)
            else:
                avail = "<br>or ".join(alt_label(a) for a in alts)

            def label(r):
                if selr is not None:
                    return ", ".join(selector_values(cells[id(r)]))
                return "%s:%s" % (os.path.basename(r["file"]), r["line"])

            def dims_of(r):
                return r["dims"] + ("; static" if r["kind"] == "static" else "")

            cols = []
            diff_cols = []
            for name, get, esc in (("dims", dims_of, md_escape),
                                   ("units", lambda r: r["units"], md_escape),
                                   ("long_name", lambda r: r["long_name"], md_escape)):
                vals = []
                for r in frs:
                    if get(r) not in vals:
                        vals.append(get(r))
                if len(vals) == 1:
                    cols.append(esc(vals[0]))
                else:
                    diff_cols.append(name)
                    cols.append("**differs**<br>" + "<br>".join(
                        "%s: %s" % (label(r), esc(get(r))) for r in frs))
            if diff_cols:
                differs.append((m, frs[0]["field"], diff_cols, frs))
            locs = []
            for r in frs:
                loc = "%s:%s" % (os.path.basename(r["file"]), r["line"])
                flags = []
                if r["built"] == "no":
                    flags.append("not built")
                if r["legacy"]:
                    flags.append("legacy")
                if r["send_status"] == "never-sent":
                    flags.append("never sent")
                if flags:
                    loc += " *(%s)*" % ", ".join(flags)
                locs.append(loc)
            fname = frs[0]["field"]
            if len({r["field"] for r in frs}) > 1:
                fname = " / ".join(sorted({r["field"] for r in frs}))
            w("| %s | %s | %s | %s | %s | %s |" % (
                md_code(fname), cols[0], cols[1], cols[2], avail, "<br>".join(locs)))
        w("")
    if differs:
        w("## Fields with inconsistent metadata")
        w("")
        w("These fields are registered in more than one place with different metadata. "
          "They should be made consistent in the source.")
        w("")
        w("| module | field | differs in | registered at |")
        w("|---|---|---|---|")
        for m, f, cols, frs in differs:
            w("| `%s` | %s | %s | %s |" % (m, md_code(f), ", ".join(cols), ", ".join(
                "%s:%s" % (os.path.basename(r["file"]), r["line"]) for r in frs)))
        w("")
    return "\n".join(out) + "\n"


# Analysis / summary
# --------------------------------------------------------------------------

def analyse(rows):
    res = {}
    by_pair = defaultdict(list)
    by_field = defaultdict(list)
    by_long = defaultdict(list)
    for r in rows:
        by_pair[(r["module"], r["field"])].append(r)
        by_field[r["field"]].append(r)
        ln = re.sub(r"[^a-z0-9]+", " ", r["long_name"].lower()).strip()
        if ln:
            by_long[ln].append(r)
    res["dup_pairs"] = {k: v for k, v in by_pair.items() if len(v) > 1}
    res["dup_fields"] = {k: v for k, v in by_field.items()
                         if len({r["module"] for r in v}) > 1}
    res["dup_long"] = {k: v for k, v in by_long.items()
                       if len({(r["module"], r["field"]) for r in v}) > 1
                       and len({r["module"] for r in v}) > 1}
    res["missing_units"] = [r for r in rows if not r["units"].strip()]
    res["missing_long"] = [r for r in rows if not r["long_name"].strip()]
    res["never_sent"] = [r for r in rows if r["send_status"] in ("never-sent", "referenced-not-sent")]
    res["no_id"] = [r for r in rows if r["send_status"] == "no-id-var"]
    res["dynamic"] = [r for r in rows if r["dynamic"]]
    res["unresolved"] = [r for r in rows if r["unresolved"]]
    res["truncated"] = [r for r in rows if "TRUNCATED" in r["notes"]]
    return res


def render_summary(rows, parsed, res):
    out = []
    w = out.append
    w("# Diagnostics inventory summary")
    w("")
    w("Generated by `diag_inventory.py` v%s." % VERSION)
    w("")
    nfiles = len({r["file"] for r in rows})
    w("- Registration rows: **%d** (%d register_diag_field, %d register_static_field) in %d files, %d modules" % (
        len(rows), sum(1 for r in rows if r["kind"] == "diag"),
        sum(1 for r in rows if r["kind"] == "static"), nfiles,
        len({r["module"] for r in rows})))
    nb = sum(1 for r in rows if r["built"] == "no")
    w("- Rows in files **not compiled** by CMake: %d" % nb)
    w("- Rows in legacy/removal-candidate directories: %d" % sum(1 for r in rows if r["legacy"]))
    w("- Rows in built, non-legacy code: %d" % sum(1 for r in rows if r["built"] != "no" and not r["legacy"]))
    trc = [r for r in res["dynamic"] if "<tracer>" in r["field"]]
    w("- Names only known at run time: %d (%d built from tracer names via get_tracer_names; %d from "
      "run-time arrays or loop counters such as aerosol/tracer name lists read from namelists or data)" % (
          len(res["dynamic"]), len(trc), len(res["dynamic"]) - len(trc)))
    w("- Module names left unresolved: %d. Literal field names left unresolved: 0 (every "
      "unresolved field name depends on run-time data)" % sum(1 for r in rows if "module_name" in r["unresolved"]))
    w("- Test-program-only rows (`#ifdef test_*`): %d" % sum(1 for r in rows if "test-program only" in r["notes"]))
    w("- Rows with a condition or call-site gate: %d; referencing a namelist variable: %d" % (
        sum(1 for r in rows if r["condition"] or r["call_gates"]),
        sum(1 for r in rows if r["namelist_vars"])))
    w("")
    w("## Counts per module")
    w("")
    w("| module | rows | built | legacy | files |")
    w("|---|---:|---:|---|---|")
    per = OrderedDict()
    for r in sorted(rows, key=lambda r: r["module"]):
        per.setdefault(r["module"], []).append(r)
    for m, rs in per.items():
        w("| `%s` | %d | %d | %s | %s |" % (
            m, len(rs), sum(1 for r in rs if r["built"] != "no"),
            "yes" if all(r["legacy"] for r in rs) else ("partly" if any(r["legacy"] for r in rs) else ""),
            ", ".join(sorted({os.path.basename(r["file"]) for r in rs}))))
    w("")
    w("## Dynamic / unresolved names")
    w("")
    if not res["dynamic"]:
        w("None.")
    else:
        w("| module | field pattern | long_name | location | loops |")
        w("|---|---|---|---|---|")
        for r in res["dynamic"]:
            w("| `%s` | `%s` | %s | %s:%s | `%s` |" % (r["module"], r["field"], md_escape(r["long_name"]),
                                                     r["file"], r["line"], r["loops"]))
    w("")

    def loc(r):
        return "%s:%s" % (os.path.basename(r["file"]), r["line"])

    w("## Same (module, field) registered at more than one call site")
    w("")
    w("diag_manager keys input fields on (module, field); if two of these call sites run in the "
      "same configuration the second registration silently reuses/overrides the first. Usually they are "
      "alternative code paths (mutually exclusive schemes).")
    w("")
    w("| module | field | call sites |")
    w("|---|---|---|")
    for (m, f), rs in sorted(res["dup_pairs"].items()):
        w("| `%s` | `%s` | %s |" % (m, f, ", ".join(loc(r) + ("" if r["built"] != "no" else " [not built]")
                                             for r in rs)))
    w("")
    w("## Same field name under different modules")
    w("")
    w("| field | modules (file:line) |")
    w("|---|---|")
    for f, rs in sorted(res["dup_fields"].items()):
        if "<" in f or "{" in f:
            continue
        w("| `%s` | %s |" % (f, "; ".join("`%s` (%s%s)" % (r["module"], loc(r), ", legacy" if r["legacy"] else "")
                                          for r in rs)))
    w("")
    w("## Same long_name under different modules (possible duplicate quantities)")
    w("")
    w("| long_name | module/field (file:line) |")
    w("|---|---|")
    for ln, rs in sorted(res["dup_long"].items()):
        w("| %s | %s |" % (md_escape(ln), "; ".join("`%s/%s` (%s)" % (r["module"], r["field"], loc(r)) for r in rs)))
    w("")
    w("## Registered but never passed to send_data (heuristic)")
    w("")
    w("The id variable returned by the register call never appears as the first argument of "
      "send_data in the same file (or in a file that USEs the module).")
    w("")
    w("| module | field | id | location | built |")
    w("|---|---|---|---|---|")
    for r in res["never_sent"]:
        w("| `%s` | `%s` | `%s` | %s | %s |" % (r["module"], r["field"], r["id_var"], loc(r), r["built"]))
    if res["no_id"]:
        w("")
        w("Registrations whose result is not stored in a variable: %d" % len(res["no_id"]))
    w("")
    w("## Missing metadata")
    w("")
    w("| module | field | missing | location |")
    w("|---|---|---|---|")
    for r in rows:
        miss = []
        if not r["units"].strip():
            miss.append("units")
        if not r["long_name"].strip():
            miss.append("long_name")
        if miss:
            w("| `%s` | `%s` | %s | %s |" % (r["module"], r["field"], ", ".join(miss), loc(r)))
    w("")
    if res["truncated"]:
        w("## Character-length truncation")
        w("")
        for r in res["truncated"]:
            w("- `%s/%s` %s: %s" % (r["module"], r["field"], loc(r), r["notes"]))
        w("")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=None,
                    help="repository root (default: the directory above this script)")
    ap.add_argument("--src", default=None, help="source dir (default: ROOT/src)")
    ap.add_argument("--md", default=None, help="markdown output (default: ROOT/docs/Diagnostics.md)")
    ap.add_argument("--csv", default=None, help="also write CSV here")
    ap.add_argument("--summary", default=None, help="also write a summary markdown here")
    ap.add_argument("--check", action="store_true",
                    help="do not write; exit 1 if --md file is missing or out of date")
    ap.add_argument("--exclude-legacy", action="store_true",
                    help="omit legacy physics (see --legacy-dirs) from the markdown")
    ap.add_argument("--include-unbuilt", action="store_true",
                    help="include files not compiled by CMake in the markdown")
    ap.add_argument("--include-test", action="store_true",
                    help="include registrations only reachable from #ifdef test programs")
    ap.add_argument("--legacy-dirs", default=",".join(DEFAULT_LEGACY_DIRS),
                    help="comma separated path components treated as legacy physics")
    ap.add_argument("--validate-diag-table", metavar="DIAG_TABLE", nargs="+", default=None,
                    help="check DIAG_TABLE(s): every field registered somewhere, files defined, "
                         "valid reductions/packing, no duplicated output names.  Does not "
                         "write the markdown unless --md is also given")
    ap.add_argument("--nml", metavar="INPUT_NML", default=None,
                    help="with --validate-diag-table: also check that each field is registered "
                         "with the namelist settings in INPUT_NML")
    ap.add_argument("--max-call-depth", type=int, default=4)
    ap.add_argument("--title", default=None)
    a = ap.parse_args(argv)

    root = os.path.abspath(a.root or os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
    src = os.path.abspath(a.src or os.path.join(root, "src"))
    md_path = a.md or os.path.join(root, "docs", "Diagnostics.md")
    legacy = [d for d in a.legacy_dirs.split(",") if d.strip()]

    rows, parsed = build_inventory(root, src, legacy, a.max_call_depth)
    md = render_markdown(rows, include_legacy=not a.exclude_legacy,
                         include_unbuilt=a.include_unbuilt, include_test=a.include_test,
                         title=a.title)

    rc = 0
    if a.validate_diag_table:
        nml = read_namelist_file(a.nml) if a.nml else None
        defaults = namelist_defaults(parsed)
        for t in a.validate_diag_table:
            rc = max(rc, validate_diag_table(t, rows, nml, defaults))
        if not (a.md or a.csv or a.summary or a.check):
            return rc

    if a.check:
        try:
            with io.open(md_path, encoding="utf-8") as fh:
                cur = fh.read()
        except IOError:
            cur = None
        if cur != md:
            sys.stderr.write("%s is out of date; regenerate with: python3 tools/diag_inventory.py\n" % md_path)
            if cur is not None:
                import difflib
                diff = list(difflib.unified_diff(cur.splitlines(), md.splitlines(),
                                                 md_path, "generated", lineterm="", n=1))
                sys.stderr.write("\n".join(diff[:60]) + "\n")
            return 1
        return rc

    d = os.path.dirname(md_path)
    if d and not os.path.isdir(d):
        os.makedirs(d)
    with io.open(md_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(md)
    if a.csv:
        write_csv(rows, a.csv)
    if a.summary:
        res = analyse(rows)
        with io.open(a.summary, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("\n".join(render_summary(rows, parsed, res)) + "\n")
    sys.stderr.write("diag_inventory: %d registrations, %d modules -> %s\n" % (
        len(rows), len({r['module'] for r in rows}), md_path))
    return rc


PLACEHOLDER_RE = re.compile(r"\{[^}]*\}|<[^>]*>")


def name_matcher(pattern, ignore_case):
    parts = PLACEHOLDER_RE.split(pattern)
    rx = "^" + ".+".join(re.escape(p) for p in parts) + "$"
    return re.compile(rx, re.I if ignore_case else 0)


REDUCTIONS = {".true.", "mean", "average", "avg", ".false.", "none", "point", "rms",
              "max", "maximum", "min", "minimum", "sum", "cumsum"}
TIME_UNITS = {"seconds", "minutes", "hours", "days", "months", "years"}


def diag_table_lines(path):
    """Yield (line number, [tokens]) for the non-blank lines of a diag_table
    (comments after '#' removed; quotes stripped from the tokens)."""
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


def validate_diag_table(path, rows, nml_values=None, defaults=None):
    """Check a diag_table against the registrations in rows.  Module names are
    case sensitive in diag_manager, field names are not.  With nml_values
    ({group: {var: value}} from an input.nml) the availability conditions
    are evaluated; namelist variables not set take their default values."""
    pats = [(name_matcher(r["module"], False), name_matcher(r["field"], True), r) for r in rows]
    bad = 0
    defaults = defaults or {}

    def value_of(group, name):
        if nml_values is None or group is None:
            return None
        g = group.lower()
        if name in nml_values.get(g, {}):
            return nml_values[g][name]
        return defaults.get((g, name))

    def report(ln, msg):
        sys.stdout.write("%s:%d: %s\n" % (path, ln, msg))

    files = {}
    outputs = {}
    with io.open(path, encoding="utf-8", errors="replace") as fh:
        head = [fh.readline().strip() for _ in range(2)]
    # diag_manager reads the title and base date from lines 1 and 2 exactly
    if not head[0] or head[0].startswith("#"):
        bad += 1
        report(1, "line 1 must be the title (a quoted string)")
    if not re.match(r"^\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\b", head[1]):
        bad += 1
        report(2, "line 2 must be the base date: year month day hour minute second")
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
        suffix = {"max": "max", "maximum": "max", "min": "min", "minimum": "min",
                  "sum": "sum", "cumsum": "sum"}.get(redl)
        if suffix and len(oname) >= 3 and oname[-3:].lower() != suffix:
            eff = oname + "_" + suffix
        key = (fname, eff)
        if key in outputs:
            bad += 1
            report(ln, "%s/%s: output name %s already used in %s at line %d" % (
                mod, fld, eff, fname, outputs[key]))
        else:
            outputs[key] = ln
        hit = [r for pm, pf_, r in pats if pm.match(mod) and pf_.match(fld)]
        if not hit:
            bad += 1
            near = sorted({r["module"] for r in rows if r["field"].lower() == fld.lower()})
            hint = (" (the field exists in module %s)" % ", ".join(near)) if near else ""
            report(ln, "%s/%s is not registered by any source file%s" % (mod, fld, hint))
            continue
        notes = []
        if all(h["built"] == "no" for h in hit):
            notes.append("only in files that are not compiled")
        alts = []
        for h in hit:
            for a in h["_alts"]:
                if a not in alts:
                    alts.append(a)
        if nml_values is not None:
            ok = alts_eval(alts, value_of)
            if ok is False:
                bad += 1
                report(ln, "%s/%s is not registered with these namelist settings; it needs %s" % (
                    mod, fld, alts_label(alts, md=False)))
                continue
            if ok is None:
                notes.append("registered only when %s" % alts_label(alts, md=False))
        elif alts_label(alts, md=False):
            notes.append("registered only when %s" % alts_label(alts, md=False))
        if all(h["kind"] == "static" for h in hit) and redl not in (".false.", "none", "point"):
            notes.append("static field: use .false.")
        if any(h["send_status"] == "never-sent" for h in hit):
            notes.append("registered but never sent")
        if notes:
            report(ln, "%s/%s ok (%s)" % (mod, fld, "; ".join(notes)))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
