"""Command line interface: python3 tools/mimadoc {generate,check,validate,list}."""

import argparse
import csv
import difflib
import io
import json
import os
import sys

from . import diag_table, markers, render
from .diagnostics import Inventory
from .evaluate import Evaluator
from .fortran import built_sources, parse_tree
from .namelists import namelist_defaults, namelist_reference, read_namelist_file

# Pages with generated regions, relative to the repository root.
PAGES = ["docs/Diagnostics.md", "docs/Parameters.md"]

CSV_COLUMNS = [
    "module", "field", "kind", "dims", "ndim", "units", "long_name", "standard_name",
    "missing_value", "range", "id_var", "file", "line", "subroutine", "fortran_module",
    "condition", "call_gates", "namelist_vars", "loops", "dynamic", "unresolved",
    "send_status", "built", "available", "notes",
]


class Model(object):
    """Everything the documentation is generated from."""

    def __init__(self, root):
        self.root = root
        src = os.path.join(root, "src")
        parsed, modtab = parse_tree(root, src)
        self.inventory = Inventory(parsed, modtab, built_sources(root, src))
        self.rows = self.inventory.registrations()
        self.namelists = namelist_reference(parsed, Evaluator(modtab))

    def region(self, kind, arg):
        if kind == "diagnostics":
            return render.render_diagnostics(self.rows)
        if kind == "namelist" and arg in self.namelists:
            return render.render_namelist(self.namelists[arg])
        raise KeyError("unknown region %s %s" % (kind, arg or ""))


def default_root():
    return os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        os.pardir, os.pardir))


def _read(path):
    with io.open(path, encoding="utf-8") as fh:
        return fh.read()


def problems(model, pages):
    """Problems with the regions themselves (not their contents)."""
    out = []
    seen = {}
    for rel, text in pages.items():
        for kind, arg in markers.regions(text):
            if kind == "namelist":
                if arg not in model.namelists:
                    out.append("%s: no namelist group %s in the source" % (rel, arg))
                seen[arg] = rel
            elif kind != "diagnostics":
                out.append("%s: unknown region kind %s" % (rel, kind))
    # once the namelist reference is generated it must be complete
    if seen:
        for g in sorted(set(model.namelists) - set(seen)):
            out.append("namelist group %s has no <!-- mimadoc:namelist %s --> region" % (g, g))
        for g in sorted(seen):
            for v in model.namelists.get(g, []):
                if not v.doc:
                    out.append("%s:%d: %s in %s has no '!!' documentation" % (
                        v.file, v.line, v.name, g))
    for r in model.rows:
        if not r.built:
            out.append("%s:%d: %s/%s is in a file that CMake does not compile (left out)" % (
                r.file, r.line, r.module, r.field))
    return out


def cmd_generate(model, args, check):
    pages = {rel: _read(os.path.join(model.root, rel)) for rel in PAGES}
    issues = problems(model, pages)
    stale = []
    for rel, text in pages.items():
        new = markers.fill(text, model.region)
        if new == text:
            continue
        if check:
            stale.append(rel)
            diff = difflib.unified_diff(text.splitlines(), new.splitlines(), rel, "generated",
                                        lineterm="", n=1)
            sys.stderr.write("\n".join(list(diff)[:60]) + "\n")
        else:
            with io.open(os.path.join(model.root, rel), "w", encoding="utf-8", newline="\n") as fh:
                fh.write(new)
            sys.stderr.write("mimadoc: updated %s\n" % rel)
    for p in issues:
        sys.stderr.write("mimadoc: %s\n" % p)
    if check:
        for rel in stale:
            sys.stderr.write("mimadoc: %s is out of date; run: python3 tools/mimadoc generate\n" % rel)
        return 1 if stale or issues else 0
    return 0


def cmd_validate(model, args):
    nml = read_namelist_file(args.nml) if args.nml else None
    defaults = namelist_defaults(model.namelists)
    rc = 0
    for t in args.diag_table:
        rc = max(rc, diag_table.validate(t, model.rows, nml, defaults))
    return rc


def _row_dict(r):
    d = {k: getattr(r, k) for k in CSV_COLUMNS}
    for k in ("namelist_vars", "unresolved"):
        d[k] = " ".join(d[k]) if k == "namelist_vars" else ",".join(d[k])
    d["notes"] = "; ".join(r.notes)
    return d


def cmd_list(model, args):
    rows = [_row_dict(r) for r in model.rows]
    if args.csv:
        with io.open(args.csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=CSV_COLUMNS, lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
    if args.json:
        data = {"fields": rows,
                "namelists": {g: [v.__dict__ for v in vs] for g, vs in sorted(model.namelists.items())}}
        with io.open(args.json, "w", encoding="utf-8") as fh:
            json.dump(data, fh, indent=1, sort_keys=True)
            fh.write("\n")
    if not (args.csv or args.json):
        for r in rows:
            sys.stdout.write("%s/%s\t%s\t%s\t%s\n" % (r["module"], r["field"], r["dims"],
                                                      r["units"], r["available"]))
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="mimadoc",
        description="Generate and check MiMA's diagnostics and namelist reference from the "
                    "Fortran sources, and check diag_table files.")
    ap.add_argument("--root", default=None,
                    help="repository root (default: two directories above this package)")
    sub = ap.add_subparsers(dest="cmd")
    sub.add_parser("generate", help="rewrite the generated regions of %s" % ", ".join(PAGES))
    sub.add_parser("check", help="exit 1 if a generated region is out of date or incomplete")
    p = sub.add_parser("validate", help="check diag_table files against the source")
    p.add_argument("diag_table", nargs="+")
    p.add_argument("--nml", metavar="INPUT_NML",
                   help="also check that each field is registered with these namelist settings")
    p = sub.add_parser("list", help="list the registered fields (and namelists)")
    p.add_argument("--csv", help="write the fields as CSV")
    p.add_argument("--json", help="write the fields and namelists as JSON")
    args = ap.parse_args(argv)
    if not args.cmd:
        ap.print_help()
        return 2
    model = Model(os.path.abspath(args.root or default_root()))
    if args.cmd in ("generate", "check"):
        return cmd_generate(model, args, check=args.cmd == "check")
    if args.cmd == "validate":
        return cmd_validate(model, args)
    return cmd_list(model, args)
