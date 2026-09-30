"""The diagnostic fields that the model can register: every call of
register_diag_field / register_static_field in the source, with its evaluated
names and metadata, the conditions under which it runs, and whether its id is
ever passed to send_data."""

import os
import re
from collections import OrderedDict, defaultdict
from dataclasses import dataclass, field
from typing import List, Tuple

from . import conditions as cnd
from .evaluate import Evaluation, Evaluator
from .fortran import guard_conds, is_trivial_cond, mask_strings, norm_ws, paren_group, split_top

REG_RE = re.compile(r"\bregister_(diag|static)_field\s*\(", re.I)
SEND_RE = re.compile(r"\bsend_(?:data|tile_averaged_data|global_diag)\s*\(", re.I)

# Positional argument order of the register functions.
DIAG_ARRAY_ARGS = ["module_name", "field_name", "axes", "init_time", "long_name",
                   "units", "missing_value", "range", "mask_variant",
                   "standard_name", "verbose"]
DIAG_SCALAR_ARGS = ["module_name", "field_name", "init_time", "long_name",
                    "units", "missing_value", "range"]
STATIC_ARGS = ["module_name", "field_name", "axes", "long_name", "units",
               "missing_value", "range", "mask_variant", "require",
               "standard_name", "dynamic"]

# How far up the call chain to look for the conditions of a registration,
# and how many call paths to keep.
MAX_CALL_DEPTH = 4
MAX_CALL_PATHS = 4


@dataclass
class Registration:
    """One registration (one row per element of an expanded loop)."""
    module: str
    field: str
    kind: str                     # 'diag' or 'static'
    axes: str                     # the axes argument as written
    dims: str                     # e.g. 'lon, lat, pfull'
    ndim: str
    units: str
    long_name: str
    standard_name: str
    missing_value: str
    range: str
    id_var: str
    file: str                     # path relative to the repository root
    line: int
    subroutine: str
    fortran_module: str
    condition: str                # Fortran text of the local condition
    call_gates: str               # conditions found up the call chain
    namelist_vars: List[str]      # 'group:var'
    loops: str
    dynamic: bool                 # name only known at run time
    unresolved: List[str]
    send_status: str              # sent, sent-elsewhere(...), referenced-not-sent, never-sent, no-id-var
    built: bool
    notes: List[str]
    alts: List[Tuple[cnd.Clause, ...]] = field(default_factory=list)

    @property
    def available(self):
        return cnd.alts_label(self.alts, md=False, sep=" | ")


def classify_args(kind, inner):
    """Map the arguments of a register call to their names."""
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
        # the scalar form has the init time where the array form has the axes
        if "axes" not in out and ("init_time" in out and len(positional) <= 2):
            names = DIAG_SCALAR_ARGS
        elif third and re.match(r"^\w*time\w*$", third, re.I) and not re.search(r"ax", third, re.I):
            names = DIAG_SCALAR_ARGS
    for k, a in enumerate(positional):
        out.setdefault(names[k] if k < len(names) else "extra%d" % k, a)
    if names is DIAG_SCALAR_ARGS:
        out["axes"] = "(scalar)"
    return out


class _Call(object):
    """A register call found in an executable statement."""

    def __init__(self, pf, ex, match):
        self.pf = pf
        self.ex = ex
        self.kind = match.group(1).lower()
        inner, _ = paren_group(ex.body, match.end() - 1)
        self.args = classify_args(self.kind, inner)
        self.line = ex.line_at(match.start())
        pre = mask_strings(ex.body)[:match.start()]
        ma = re.match(r"^\s*([\w%]+)\s*(\([^=]*\))?\s*=\s*$", pre)
        self.id_var = ma.group(1).split("%")[-1].lower() if ma else ""
        self.id_expr = norm_ws(ex.body[:match.start()].rstrip().rstrip("=")) if ma else ""
        ms = ex.scope.module_scope()
        self.fortran_module = ms.name if ms is not None and ms.kind != "file" else ""
        ps = ex.scope.proc_scope()
        self.proc = ps.name if ps is not None else ""


class Inventory(object):
    """All registrations in a parsed source tree."""

    def __init__(self, parsed, modtab, built):
        self.parsed = parsed
        self.modtab = modtab
        self.built = built
        self.ev = Evaluator(modtab)
        # namelist groups of each module's variables, for USE'd variables
        self.module_namelists = defaultdict(dict)
        for pf in parsed:
            for nl in pf.namelist_stmts:
                ms = nl.scope.module_scope()
                for v in nl.vars:
                    self.module_namelists[ms.name][v] = nl.group
        self.calls = []
        self.sends = defaultdict(set)        # file -> id names sent
        self.idrefs = defaultdict(lambda: defaultdict(int))
        for pf in parsed:
            for ex in pf.execs:
                masked = mask_strings(ex.body)
                for msd in SEND_RE.finditer(masked):
                    inner, _ = paren_group(ex.body, msd.end() - 1)
                    mmid = re.match(r"^([\w%]+)", split_top(inner)[0].strip())
                    if mmid:
                        self.sends[pf.relpath].add(mmid.group(1).split("%")[-1].lower())
                for tok in re.findall(r"[A-Za-z_]\w*", masked):
                    self.idrefs[pf.relpath][tok.lower()] += 1
                for mr in REG_RE.finditer(masked):
                    self.calls.append(_Call(pf, ex, mr))
        self.callers = defaultdict(list)     # callee -> [(pf, scope, conds, line)]
        for pf in parsed:
            for s in pf.scopes:
                for callee, si, conds, line in s.calls:
                    self.callers[callee].append((pf, s, conds, line))

    # -- helpers
    def group_of(self, pf, scope):
        """Function giving the namelist group of a variable seen from scope:
        the file's own namelists, then those of the modules it uses."""
        def get(var):
            if var in pf.namelists:
                return pf.namelists[var]
            for s in scope.chain():
                for u in s.uses:
                    only = s.use_only.get(u)
                    if (only is None or var in only) and var in self.module_namelists.get(u, {}):
                        return self.module_namelists[u][var]
            return None
        return get

    @staticmethod
    def file_uses(pf, modname, name=None):
        for s in pf.scopes:
            if modname in s.uses:
                only = s.use_only.get(modname, None)
                if name is None or only is None or name in only:
                    return True
        return False

    def call_paths(self, pf, proc_scope, depth=0, seen=frozenset()):
        """The conditions on the call paths into a procedure.  Returns a list
        of paths; each path is a list of (where, [conds], group_of) from the
        innermost caller outward, keeping only call sites with conditions.
        An empty list means the procedure is called unconditionally (or is
        never called)."""
        if proc_scope is None or depth > MAX_CALL_DEPTH or proc_scope.name in seen:
            return []
        name = proc_scope.name
        seen = seen | {name}
        ms = proc_scope.module_scope()
        defmod = ms.name if ms is not None else ""
        paths = []
        for cpf, cscope, conds, line in self.callers.get(name, []):
            if cpf is not pf:
                # another file must USE the defining module (and not define
                # a procedure of the same name itself)
                if name in cpf.proc_defs or not self.file_uses(cpf, defmod, name):
                    continue
            where = "%s:%d (%s)" % (os.path.basename(cpf.relpath), line, cscope.name)
            csc = cscope.proc_scope()
            cc = [c for c in conds if c and not is_trivial_cond(c)] + (guard_conds(csc.guards) if csc is not None else [])
            link = [(where, cc, self.group_of(cpf, cscope))] if cc else []
            up = self.call_paths(cpf, csc, depth + 1, seen) if csc is not None else []
            if up:
                paths.extend(link + u for u in up)
            else:
                paths.append(link)
        out = []
        for p in paths:
            if not p:
                return []           # an unconditional path: the others do not matter
            key = [(w, c) for w, c, g in p]
            if key not in [[(w, c) for w, c, g in q] for q in out]:
                out.append(p)
        return out[:MAX_CALL_PATHS]

    def send_status(self, call):
        pf = call.pf
        idv = call.id_var
        if not idv:
            return "no-id-var"
        if idv in self.sends[pf.relpath]:
            return "sent"
        elsewhere = sorted(os.path.basename(p.relpath) for p in self.parsed
                           if idv in self.sends[p.relpath]
                           and self.file_uses(p, call.fortran_module))
        if elsewhere:
            return "sent-elsewhere(%s)" % ",".join(elsewhere)
        nreg = sum(1 for c in self.calls if c.pf is pf and c.id_var == idv)
        return "referenced-not-sent" if self.idrefs[pf.relpath].get(idv, 0) > nreg else "never-sent"

    # -- the inventory
    def registrations(self):
        rows = []
        for call in self.calls:
            rows.extend(self._rows(call))
        rows.sort(key=lambda r: (r.module.lower(), r.field.lower(), r.file, r.line))
        return rows

    def _loop_envs(self, call):
        """One environment per element of the literal arrays indexed by a DO
        variable in the names (the registration is expanded once per element)."""
        ex = call.ex
        info = Evaluation()
        for key in ("module_name", "field_name", "long_name", "units"):
            if key in call.args:
                self.ev.eval_str(call.args[key], ex.scope, ex.index, {}, info)
        loopvars = {}
        for lh in ex.loops:
            mlh = re.match(r"^(\w+)\s*=\s*(.+)$", lh)
            if mlh:
                loopvars[mlh.group(1).lower()] = split_top(mlh.group(2))
        envs = [{}]
        for var, (arr, n) in info.expand.items():
            lo, hi = 1, n
            if var in loopvars:
                b = loopvars[var]
                l0 = self.ev.eval_int(b[0], ex.scope, {}) if len(b) > 0 else None
                h0 = self.ev.eval_int(b[1], ex.scope, {}) if len(b) > 1 else None
                lo = l0 if l0 is not None else 1
                hi = min(h0 if h0 is not None else n, n)
            envs = [dict(e, **{var: k}) for e in envs for k in range(lo, hi + 1)]
        return envs

    def _rows(self, call):
        pf, ex = call.pf, call.ex
        scope = ex.scope
        psc = scope.proc_scope()
        conds = [c for c in ex.conds if c] + (guard_conds(psc.guards, ex.index) if psc else [])
        group_of = self.group_of(pf, scope)
        paths = self.call_paths(pf, psc)
        send = self.send_status(call)

        local = [c for t in conds for c in cnd.clauses(t, group_of)]
        alts = []
        for p in paths:
            alt = list(local)
            for w, cc, g in p:
                for t in cc:
                    for c in cnd.clauses(t, g):
                        if c not in alt:
                            alt.append(c)
            if tuple(alt) not in alts:
                alts.append(tuple(alt))
        if not paths:
            alts = [tuple(local)]
        nml = sorted({"%s:%s" % (g, v) for a in alts for c in a for v, g in c.groups})
        gates = " | ".join(" <- ".join("%s: %s" % (w, " .and. ".join(cc)) for w, cc, g in p)
                           for p in paths)
        dims, ndim = self.ev.axes(call.args.get("axes"), scope, ex.index)
        ndim = "" if ndim is None else str(ndim)

        rows = []
        for env in self._loop_envs(call):
            info = Evaluation()
            vals, res = {}, {}
            for key in ("module_name", "field_name", "long_name", "units", "standard_name"):
                if key in call.args:
                    v, r = self.ev.eval_str(call.args[key], scope, ex.index, env, info)
                    vals[key] = v.strip() if r else v
                    res[key] = r
                else:
                    vals[key], res[key] = "", True
            unresolved = [k for k in ("module_name", "field_name") if not res[k]]
            notes = list(info.notes)
            if env:
                notes.append("expanded loop: " + ", ".join("%s=%d" % kv for kv in sorted(env.items())))
            notes += ["%s not resolved" % k for k in ("long_name", "units") if not res[k]]
            rows.append(Registration(
                module=vals["module_name"], field=vals["field_name"], kind=call.kind,
                axes=norm_ws(call.args.get("axes", "")),
                dims=dims or ("%sD" % ndim if ndim else ""), ndim=ndim,
                units=vals["units"], long_name=norm_ws(vals["long_name"]),
                standard_name=vals["standard_name"],
                missing_value=norm_ws(call.args.get("missing_value", "")),
                range=norm_ws(call.args.get("range", "")), id_var=call.id_expr,
                file=pf.relpath.replace("\\", "/"), line=call.line, subroutine=call.proc,
                fortran_module=call.fortran_module, condition=" .and. ".join(conds),
                call_gates=gates, namelist_vars=nml, loops=" ; ".join(ex.loops),
                dynamic=bool(unresolved) or bool(env), unresolved=unresolved,
                send_status=send,
                built=self.built is None or pf.relpath in self.built,
                notes=notes, alts=alts))
        return rows


def merge_fields(rows):
    """Group the rows of one module by field name (field names are not case
    sensitive in diag_manager): OrderedDict lowercase name -> [rows]."""
    out = OrderedDict()
    for r in rows:
        out.setdefault(r.field.lower(), []).append(r)
    return out
