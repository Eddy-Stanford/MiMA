"""MkDocs hook: the Fortran API reference as pages of the documentation site.

FORD (configured by ford.md) parses the Fortran sources and their doc
comments; this hook turns the result into one Markdown page per module (and
for the coupler_main program) under api/, listed in the navigation under
"Fortran API reference" and indexed by the site search. Enabled in
mkdocs.yml with `hooks: [tools/fortran_api.py]`."""

import contextlib
import io
import os
import pathlib
import re
import textwrap

import posixpath
import sys

import ford
from mkdocs.structure.files import File

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from mimadoc.cli import Model  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
SOURCE_URL = "https://github.com/Eddy-Stanford/MiMA/blob/main/"
SITE_URL = "https://eddy-stanford.github.io/MiMA/"
# Links to the published site, as written in docs/*.md and in doc comments
# (so that they also work on GitHub and in the source):
#   https://eddy-stanford.github.io/MiMA/              -> README.md
#   https://eddy-stanford.github.io/MiMA/Page/#anchor  -> Page.md#anchor
#   https://eddy-stanford.github.io/MiMA/api/mod/#x    -> api/mod.md#x
SITE_LINK_RE = re.compile(re.escape(SITE_URL) + r"((?:api/)?\w+/)?(#[\w-]+)?(?=[)\s>\]]|$)")
NAV_TITLE = "Fortran API reference"
OVERVIEW = "FortranAPI.md"

# Navigation groups, by the first matching prefix of the source path.
GROUPS = [
    ("Dynamical core", ("src/atmos_spectral/",)),
    ("Radiation", ("src/atmos_param/radiation/",)),
    ("Physics", ("src/atmos_param/",)),
    ("Coupler and surface", ("src/coupler/", "src/atmos_coupled/")),
    ("Shared", ("src/atmos_shared/", "src/mima_shared/")),
]

_units = []      # [(name, kind, entity, relpath)] of the last parse
_pages = {}      # api/<name>.md -> Markdown
_type_pages = {} # derived type name -> module page that defines it
_diag_rows = []  # mimadoc registrations (for the links to docs/Diagnostics.md)


def _diag_anchor(module):
    from mimadoc.render import anchor
    return anchor(module)


# --------------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------------

def _parse():
    """Parse the sources with FORD; returns [(name, kind, entity, relpath)]."""
    os.environ["PATH"] = os.path.dirname(os.sys.executable) + os.pathsep + os.environ["PATH"]
    text = (ROOT / "ford.md").read_text()
    with contextlib.redirect_stdout(io.StringIO()):
        docs, settings = ford.load_settings(text, ROOT, "ford.md")
        settings, docs = ford.parse_arguments({"quiet": True}, docs, settings, ROOT)
        project = ford.fortran_project.Project(settings)
        project.correlate()
    units = [(m.name, "module", m) for m in project.modules] + \
            [(p.name, "program", p) for p in project.programs]
    out = []
    for name, kind, ent in units:
        rel = os.path.relpath(str(ent.parent.path), ROOT).replace(os.sep, "/")
        out.append((name, kind, ent, rel))
    return out, settings.extra_mods


def _group(relpath):
    for title, prefixes in GROUPS:
        if relpath.startswith(prefixes):
            return title
    return "Other"


# --------------------------------------------------------------------------
# Markdown
# --------------------------------------------------------------------------

def _doc(ent):
    """The doc comment of an entity as Markdown."""
    lines = list(getattr(ent, "doc_list", None) or [])
    return textwrap.dedent("\n".join(lines)).strip()


def _summary(ent):
    first = _doc(ent).split("\n\n")[0]
    return re.sub(r"\s+", " ", first)


def _type(v, link=False):
    """Fortran type of a variable or argument, e.g. real, dimension(:, :).
    With link, MiMA's derived types link to their definition."""
    t = v.vartype
    if v.vartype in ("type", "class") and v.proto:
        name = re.sub(r"<[^>]+>", "", str(getattr(v.proto[0], "name", v.proto[0]))).strip()
        page = _type_pages.get(name.lower())
        if link and page:
            name = "[%s](%s.md#%s)" % (name, page, name.lower())
        t = "%s(%s)" % (v.vartype, name)
    elif v.vartype == "character" and v.strlen:
        t = "character(len=%s)" % v.strlen
    elif v.kind:
        t = "%s(kind=%s)" % (v.vartype, v.kind)
    attrs = [a for a in (v.attribs or []) if a.lower() not in ("optional",)]
    if v.dimension:
        attrs.append("dimension%s" % v.dimension)
    return ", ".join([t] + attrs)


def _cell(text):
    return re.sub(r"\s+", " ", text).replace("|", "\\|").strip()


def _code(s):
    return "`%s`" % s


def _grouped(variables):
    """Variables with the same type, intent and doc (declared in one
    statement) share a row, at the position of the first of them."""
    rows = []
    index = {}
    for v in variables:
        key = (_type(v, link=True), getattr(v, "intent", ""),
               bool(getattr(v, "optional", False)), _doc(v))
        if key[3] and key in index:
            rows[index[key]][1].append(v)
        else:
            index[key] = len(rows)
            rows.append((key, [v]))
    return rows


def _arguments(proc):
    rows = _grouped(proc.args)
    if not rows:
        return []
    out = ["| Argument | Type | Intent | Description |", "|---|---|---|---|"]
    for (typ, intent, optional, doc), vs in rows:
        out.append("| %s | %s | %s | %s |" % (
            ", ".join(_code(v.name) for v in vs), _cell(typ),
            (intent or "") + (", optional" if optional else ""), _cell(doc)))
    return out


def _signature(proc, kind):
    args = ", ".join(a.name for a in proc.args)
    if kind == "function":
        rv = proc.retvar
        rname = getattr(rv, "name", rv)
        res = "" if rname == proc.name else " result(%s)" % rname
        rtype = _type(rv) if hasattr(rv, "vartype") else ""
        return "%s function %s(%s)%s" % (rtype, proc.name, args, res)
    return "subroutine %s(%s)" % (proc.name, args)


def _procedure(proc, kind, level):
    out = ["%s `%s`" % ("#" * level, proc.name), ""]
    out += ["```fortran", _signature(proc, kind).strip(), "```", ""]
    if _doc(proc):
        out += [_doc(proc), ""]
    args = _arguments(proc)
    if args:
        out += args + [""]
    return out


def _interface(iface, level):
    out = ["%s `%s`" % ("#" * level, iface.name), "", "Generic interface.", ""]
    if _doc(iface):
        out += [_doc(iface), ""]
    for mp in iface.modprocs or []:
        proc = getattr(mp, "procedure", None)
        if proc is None or isinstance(proc, str):
            continue
        kind = "function" if proc.obj == "function" else "subroutine"
        out += _procedure(proc, kind, level + 1)
    return out


def _variables(variables, header):
    rows = _grouped(variables)
    out = ["| %s | Type | Value | Description |" % header, "|---|---|---|---|"]
    for (typ, intent, optional, doc), vs in rows:
        vals = [v.initial for v in vs if v.initial not in (None, "")]
        if vals and len(set(vals)) == 1:
            values = _code(vals[0]) + (" (each)" if len(vs) > 1 else "")
        else:
            values = ", ".join(_code(x) for x in vals)
        out.append("| %s | %s | %s | %s |" % (
            ", ".join(_code(v.name) for v in vs), _cell(typ), _cell(values), _cell(doc)))
    return out


def _module_link(u, names, extra_mods):
    name = getattr(u, "name", u)
    if name in names:
        return "[`%s`](%s.md)" % (name, name)
    if name in extra_mods:
        return "[`%s`](%s)" % (name, extra_mods[name])
    return _code(name)


def _render(name, kind, ent, rel, names, extra_mods):
    public = (lambda xs: [x for x in xs if getattr(x, "permission", "public") == "public"]) \
        if kind == "module" else (lambda xs: list(xs))
    out = ["# `%s`" % name, ""]
    if _doc(ent):
        out += [_doc(ent), ""]
    facts = ["**%s** in [%s](%s%s)" % (kind.capitalize(), rel, SOURCE_URL, rel)]
    groups = [n.name.lower() for n in ent.namelists or []]
    if groups:
        facts.append("**Namelist:** " + ", ".join(
            "[`%s`](../Parameters.md#%s)" % (g, g) for g in groups))
    diag = sorted({r.module for r in _diag_rows if r.fortran_module == name.lower()})
    if diag:
        facts.append("**Diagnostics:** " + ", ".join(
            "[`%s`](../Diagnostics.md#%s)" % (d, _diag_anchor(d)) for d in diag))
    uses = sorted({getattr(u, "name", u) for u in ent.uses or []})
    if uses:
        facts.append("**Uses:** " + ", ".join(_module_link(u, names, extra_mods) for u in uses))
    out += ["  \n".join(facts), ""]

    types = public(ent.types or [])
    ifaces = [i for i in public(ent.interfaces or []) if i.generic or i.modprocs]
    subs = public(ent.subroutines or [])
    funcs = public(ent.functions or [])
    variables = public(ent.variables or []) if kind == "module" else []
    if kind == "program":
        subs, funcs = list(ent.subroutines or []), list(ent.functions or [])

    if variables:
        out += ["## Variables and parameters", ""] + _variables(variables, "Name") + [""]
    if types:
        out += ["## Derived types", ""]
        for t in types:
            out += ["### `%s`" % t.name, ""]
            if _doc(t):
                out += [_doc(t), ""]
            if t.variables:
                out += _variables(t.variables, "Component") + [""]
    if ifaces:
        out += ["## Interfaces", ""]
        for i in ifaces:
            out += _interface(i, 3)
    if funcs:
        out += ["## Functions", ""]
        for f in funcs:
            out += _procedure(f, "function", 3)
    if subs:
        out += ["## Subroutines", ""] if kind == "module" else ["## Internal procedures", ""]
        for s in subs:
            out += _procedure(s, "subroutine", 3)
    return "\n".join(out).rstrip() + "\n"


def _index():
    """The module list appended to the overview page."""
    out = ["", "## Modules", ""]
    for title, _ in GROUPS + [("Other", ())]:
        units = [u for u in _units if _group(u[3]) == title]
        if not units:
            continue
        out += ["### %s" % title, "", "| Module | Summary |", "|---|---|"]
        for name, kind, ent, rel in units:
            out.append("| [`%s`](api/%s.md) | %s |" % (name, name, _cell(_summary(ent))))
        out.append("")
    return "\n".join(out)


# --------------------------------------------------------------------------
# MkDocs events
# --------------------------------------------------------------------------

def on_config(config):
    global _units, _pages
    units, extra_mods = _parse()
    model = Model(str(ROOT))
    _diag_rows[:] = [r for r in model.rows if r.built]
    order = {title: k for k, (title, _) in enumerate(GROUPS)}
    _units = sorted(units, key=lambda u: (order.get(_group(u[3]), 99), u[0]))
    names = {u[0] for u in _units}
    _type_pages.clear()
    for name, kind, ent, rel in _units:
        for t in ent.types or []:
            _type_pages[t.name.lower()] = name
    _pages = {"api/%s.md" % name: _render(name, kind, ent, rel, names, extra_mods)
              for name, kind, ent, rel in _units}
    # the overview page becomes a section with the module pages, by group
    section = [{"Overview": OVERVIEW}]
    for title, _ in GROUPS + [("Other", ())]:
        units = [u for u in _units if _group(u[3]) == title]
        if units:
            section.append({title: [{u[0]: "api/%s.md" % u[0]} for u in units]})
    nav = []
    for item in config["nav"]:
        if isinstance(item, dict) and item.get(NAV_TITLE) == OVERVIEW:
            nav.append({NAV_TITLE: section})
        else:
            nav.append(item)
    config["nav"] = nav
    return config


def on_files(files, config):
    for path, content in _pages.items():
        files.append(File.generated(config, path, content=content))
    return files


def _relative_links(markdown, page):
    """Links to the published site become relative links to the pages."""
    here = posixpath.dirname(page.file.src_uri) or "."

    def sub(m):
        path = (m.group(1) or "").rstrip("/")
        target = (path + ".md") if path else "README.md"
        return posixpath.relpath(target, here) + (m.group(2) or "")
    return SITE_LINK_RE.sub(sub, markdown)


def on_page_markdown(markdown, page, config, files):
    if page.file.src_uri == OVERVIEW:
        markdown = markdown.rstrip("\n") + "\n" + _index()
    return _relative_links(markdown, page)
