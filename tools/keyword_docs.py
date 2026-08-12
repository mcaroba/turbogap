#!/usr/bin/env python3
"""Generate the TurboGAP input-keyword reference from the source that parses it.

The keyword table has exactly one authority, src/read_files.f90, because that
is the only place where a keyword exists at all: input_parameters holds fields,
and the two sets do not line up.  Roughly a third of the fields are internal
state with no keyword, several keywords write more than one field, a few write
a field named differently from themselves (nd_n_samples writes xrd_n_samples),
and the aliases and the deprecations live only here.  So the prose lives here
too, in a !> block above the branch that handles the keyword.

Everything a machine can already read is NOT repeated in that block:

  group        the enclosing read_options_* subroutine
  fields       the params%<field> the branch assigns
  type         the field's declaration in src/types.f90
  default      the field's initialiser in src/types.f90
  side effects the other params%<field> = the branch assigns

A default written into a comment is a default that will eventually disagree
with the code.  Harvesting it from types.f90 means the table cannot drift.

Outputs
  docs/keywords.md     the reference, for the repo and the wiki
  docs/keywords.html   the same as a standalone searchable page
  src/keyword_help.f90 a generated module backing `turbogap --help [mode]`

Usage
  tools/keyword_docs.py               write all three
  tools/keyword_docs.py --check       exit 1 if any output is stale or any
                                      keyword is undocumented (for CI/hooks)
  tools/keyword_docs.py --list-undocumented
"""

import argparse
import html
import os
import re
import sys
from collections import OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
READ_FILES = os.path.join(ROOT, "src", "read_files.f90")
TYPES = os.path.join(ROOT, "src", "types.f90")
DOC_MD = os.path.join(ROOT, "docs", "keywords.md")
DOC_HTML = os.path.join(ROOT, "docs", "keywords.html")
GEN_F90 = os.path.join(ROOT, "src", "keyword_help.f90")

MODES = ["predict", "md", "mc", "soap"]
#  What `turbogap --help <topic>` accepts. The modes filter the input-file
#  listing; "gap" selects the potential file instead, which has no modes.
HELP_TOPICS = MODES + ["gap"]

INPUT_FILE = "the input file"
GAP_FILE = "the potential (.gap) file"

# read_options_<suffix> -> (title, blurb).  Anything not listed still appears,
# under its bare subroutine suffix.
GROUPS = OrderedDict([
    ("general", ("General", "The atoms, the species and how chatty the run is.")),
    ("control", ("Run control", "What is computed, and how a relaxation is driven.")),
    ("md", ("Molecular dynamics", "Time stepping, thermostat and barostat.")),
    ("nested", ("Nested sampling", "Nested-sampling walks.")),
    ("mc", ("Monte Carlo", "Move types, acceptance and the grand-canonical ensemble.")),
    ("vdw", ("Van der Waals", "Tkatchenko-Scheffler and many-body dispersion.")),
    ("estat", ("Electrostatics", "Charge-charge interactions and their damping.")),
    ("exp", ("Experimental data (MAD)", "XPS, pair distributions, structure factors and diffraction.")),
    ("stopping", ("Electronic stopping and EPH", "Radiation cascades and electron-phonon coupling.")),
    ("local_properties", ("Local properties", "Per-atom quantities predicted alongside the energy.")),
    ("output", ("Output", "What is written, and how often.")),
    ("gpu", ("GPU", "Accepted by the CPU build too, so one deck runs on both.")),

    # .gap blocks. The key is the descriptor name that follows gap_beg.
    ("soap_turbo", ("soap_turbo", "The SOAP descriptor and the GAP fitted on it, "
                    "one block per central species. Also carries the local-property "
                    "models -- Hirshfeld volumes, charges, core-electron binding "
                    "energies -- that ride on the same descriptor.")),
    ("distance_2b", ("distance_2b", "A two-body GAP on the interatomic distance, "
                     "one block per species pair.")),
    ("angle_3b", ("angle_3b", "A three-body GAP on the two distances and the angle "
                  "between them, one block per centre-and-pair combination.")),
    ("core_pot", ("core_pot", "A tabulated pair potential added on top of the GAPs, "
                  "one block per species pair. Splined on read.")),
])

#  Which read_gap_hypers block writes into which derived type in types.f90.
GAP_BLOCKS = OrderedDict([
    ("soap_turbo", "soap_turbo"),
    ("distance_2b", "distance_2b"),
    ("angle_3b", "angle_3b"),
    ("core_pot", "core_pot"),
])


class Keyword(object):
    def __init__(self, names, group, line, section=INPUT_FILE, type_name=""):
        self.names = names          # every spelling that reaches this branch
        self.group = group
        self.section = section      # which file the keyword belongs in
        self.type_name = type_name  # derived type in types.f90 holding its fields
        self.line = line            # 1-based line of the branch header
        self.desc = []              # prose lines
        self.units = ""
        self.modes = []             # empty means every mode
        self.needs = []
        self.see = []
        self.status = ""            # "", "ignored" or "deprecated"
        self.type_override = ""     # only where the automatic join cannot work
        self.default_override = ""  # ditto
        self.fields = []            # params% fields the branch assigns
        self.documented = False

    @property
    def name(self):
        return self.names[0]

    @property
    def aliases(self):
        return self.names[1:]


# ---------------------------------------------------------------- types.f90

def parse_types(path):
    """derived type name (lowercased) -> {component: (type, default)}.

    Every type in the file, not only input_parameters: the .gap keywords write
    into soap_turbo, distance_2b, angle_3b and core_pot, and their defaults
    have to come from the same place for the same reason.
    """
    lines = open(path).read().split("\n")
    spans, open_at, open_name = [], None, None
    for i, l in enumerate(lines):
        m = re.match(r"\s*type\s+(\w+)\s*$", l)
        if m and open_at is None:
            open_at, open_name = i, m.group(1)
            continue
        if re.match(r"\s*end type\b", l) and open_at is not None:
            spans.append((open_name.lower(), open_at, i))
            open_at = open_name = None
    if not any(n == "input_parameters" for n, _, _ in spans):
        sys.exit("could not find `type input_parameters` in " + path)
    return dict((name, _parse_type_body(lines[a + 1:b])) for name, a, b in spans)


def _parse_type_body(raw_lines):
    decl = re.compile(
        r"^\s*(real\(dp\)|integer|logical|character\*\d+|type\([\w ]+\))"
        r"\s*(,[^:]*)?::\s*(.*)$")
    out = {}
    for raw in raw_lines:
        if not raw.strip() or raw.strip().startswith("!"):
            continue
        m = decl.match(raw)
        if not m:
            continue
        ftype, attrs, rest = m.group(1), (m.group(2) or ""), m.group(3)
        allocatable = "allocatable" in attrs
        dim = ""
        dm = re.search(r"dimension\(([^)]*)\)", attrs)
        if dm:
            dim = dm.group(1)
        for part in split_top_level(rest):
            part = part.strip()
            if not part:
                continue
            nm = re.match(r"(\w+)", part)
            if not nm:
                continue
            name = nm.group(1)
            shape = ""
            sm = re.match(r"\w+\(([^)]*)\)", part)
            if sm:
                shape = sm.group(1)
            default = ""
            if "=" in part and "=>" not in part:
                default = part.split("=", 1)[1].strip()
            out[name.lower()] = (pretty_type(ftype, allocatable, shape or dim),
                                 pretty_default(default))
    return out


def split_top_level(text):
    parts, depth, cur = [], 0, ""
    for ch in text:
        if ch in "([":
            depth += 1
        elif ch in ")]":
            depth -= 1
        if ch == "," and depth == 0:
            parts.append(cur)
            cur = ""
        else:
            cur += ch
    parts.append(cur)
    return parts


def pretty_type(ftype, allocatable, shape):
    base = {"real(dp)": "real", "integer": "integer", "logical": "logical"}.get(ftype)
    if base is None:
        if ftype.startswith("character"):
            base = "string"
        else:
            base = ftype
    if allocatable and not shape:
        return base + " list"
    if shape:
        if allocatable:
            return base + " list"
        return "%s(%s)" % (base, shape)
    return base


def pretty_default(default):
    if not default:
        return ""
    d = default.strip()
    d = re.sub(r"^\.(true|false)\.$", lambda m: m.group(1), d)
    d = re.sub(r'^"(.*)"$', r"\1", d)
    d = re.sub(r"^'(.*)'$", r"\1", d)
    if d.startswith("reshape("):
        return "identity"
    # 1.d-6 -> 1e-6, 300.d0 -> 300.0
    if re.match(r"^-?[\d.]+[dD][-+]?\d+$", d):
        mant, exp = re.split(r"[dD]", d)
        if int(exp) == 0:
            d = mant + "0" if mant.endswith(".") else mant
            if "." not in d:
                d += ".0"
        else:
            d = "%se%d" % (mant, int(exp))
    elif re.match(r"^-?\d+\.$", d):
        d += "0"
    if d.startswith("(/") and d.endswith("/)"):
        d = "[" + d[2:-2].strip() + "]"
    if len(d) > 34:
        d = d[:31] + "..."
    return d


# ----------------------------------------------------------- read_files.f90

DOC_RE = re.compile(r"^\s*!>\s?(.*)$")
BRANCH_RE = re.compile(r"^\s*(?:else\s+)?if\s*\(\s*(?:trim\()?keyword")
KW_RE = re.compile(r"keyword\s*==\s*['\"]([^'\"]+)['\"]")
FIELD_RE = re.compile(r"params\s*%\s*(\w+)")


def parse_read_files(path):
    lines = open(path).read().split("\n")

    # subroutine spans
    spans, cur = [], None
    for i, l in enumerate(lines):
        m = re.match(r"\s*subroutine\s+(\w+)", l)
        if m:
            cur = (m.group(1), i)
        elif re.match(r"\s*end subroutine", l) and cur:
            spans.append((cur[0], cur[1], i))
            cur = None

    keywords = []
    for sub, a, b in spans:
        if not sub.startswith("read_options_"):
            continue
        group = sub[len("read_options_"):]
        body = lines[a:b]
        # branch starts, relative to a
        starts = [j for j, l in enumerate(body) if BRANCH_RE.match(l)]
        for n, j in enumerate(starts):
            end = starts[n + 1] if n + 1 < len(starts) else len(body)
            header = body[j]
            k = j
            while header.rstrip().endswith("&") and k + 1 < len(body):
                k += 1
                header += body[k]
            names = KW_RE.findall(header)
            if not names:
                continue
            kw = Keyword(names, group, a + j + 1)
            text = "\n".join(body[j:end])
            seen = []
            for f in FIELD_RE.findall(text):
                if f.lower() not in [s.lower() for s in seen]:
                    seen.append(f)
            kw.fields = seen
            attach_doc(kw, body, j)
            keywords.append(kw)
    return keywords


HYPERS_RE = re.compile(r"(\w+)_hypers\s*\([^)]*\)((?:\s*%\s*\w+(?:\([^)]*\))?)+)")


def parse_gap_hypers(path):
    """The .gap keywords, out of read_gap_hypers in the same file.

    The file's own structure is the grouping: `gap_beg <descriptor>` opens a
    block and `gap_end` closes it, and read_gap_hypers mirrors that with one
    `do while (keyword /= "gap_end")` loop per descriptor. A branch counts as
    a keyword of a block when it lies inside that block's loop, which also
    drops the two gap_beg dispatchers -- those are file syntax, not settings,
    and are described in the section preamble instead.
    """
    lines = open(path).read().split("\n")
    try:
        a = next(i for i, l in enumerate(lines)
                 if l.strip().startswith("subroutine read_gap_hypers"))
    except StopIteration:
        return []
    b = next(i for i, l in enumerate(lines[a:], a) if l.strip() == "end subroutine")
    body = lines[a:b]

    block_re = re.compile(r'^(\s*)(?:else\s+)?if\s*\(\s*keyword\s*==\s*"(%s)"\s*\)\s*then'
                          % "|".join(GAP_BLOCKS))
    loop_re = re.compile(r"^(\s*)do while \(")

    # Every do-while the dispatcher opens belongs to its block. soap_turbo has
    # two: a pre-scan that finds n_species, because the array-valued keywords
    # cannot be allocated without it, and then the real loop over the rest.
    loops, block, block_indent, loop_indent = [], None, -1, -1
    for j, l in enumerate(body):
        m = block_re.match(l)
        if m:
            block, block_indent, loop_indent = m.group(2), len(m.group(1)), -1
            continue
        m = loop_re.match(l)
        if not m or block is None:
            continue
        indent = len(m.group(1))
        if indent <= block_indent:
            continue
        if loop_indent < 0:
            loop_indent = indent
        elif indent != loop_indent:
            continue        # a loop nested inside one of this block's branches
        end = len(body)
        for k in range(j + 1, len(body)):
            s = body[k]
            if s.strip().startswith("end do") and len(s) - len(s.lstrip()) == indent:
                end = k
                break
        loops.append((block, j, end))

    keywords, seen = [], set()
    for block, lo, hi in loops:
        starts = [j for j in range(lo, hi) if BRANCH_RE.match(body[j])]
        for n, j in enumerate(starts):
            end = starts[n + 1] if n + 1 < len(starts) else hi
            header = body[j]
            k = j
            while header.rstrip().endswith("&") and k + 1 < len(body):
                k += 1
                header += body[k]
            names = KW_RE.findall(header)
            if not names:
                continue
            kw = Keyword(names, block, a + j + 1,
                         section=GAP_FILE, type_name=GAP_BLOCKS[block])
            fields = []
            for m in HYPERS_RE.finditer("\n".join(body[j:end])):
                path_ = re.sub(r"\([^)]*\)", "", m.group(2)).strip("%")
                # Components of a nested type (local_property_models%zeta) are
                # an implementation detail of the block, not a second keyword.
                if "%" in path_ or path_ in fields:
                    continue
                fields.append(path_)
            kw.fields = fields
            attach_doc(kw, body, j)
            if (block, kw.name) in seen:
                continue        # already claimed by an earlier pass of the block
            seen.add((block, kw.name))
            keywords.append(kw)
    return keywords


def attach_doc(kw, body, j):
    """Read the !> block sitting immediately above branch line j."""
    block = []
    i = j - 1
    while i >= 0:
        m = DOC_RE.match(body[i])
        if not m:
            break
        block.append(m.group(1).rstrip())
        i -= 1
    if not block:
        return
    block.reverse()
    kw.documented = True
    for line in block:
        m = re.match(r"@(\w+)\s*(.*)$", line.strip())
        if not m:
            if line.strip() or kw.desc:
                kw.desc.append(line.strip())
            continue
        tag, val = m.group(1).lower(), m.group(2).strip()
        if tag == "kw":
            declared = [x.strip() for x in val.split(",") if x.strip()]
            actual = [x.strip() for x in kw.names]
            if declared != actual:
                sys.stderr.write(
                    "%s:%d: @kw says %s but the branch tests %s\n"
                    % (READ_FILES, kw.line, ", ".join(declared), ", ".join(actual)))
                sys.exit(2)
        elif tag == "type":
            kw.type_override = val
        elif tag == "default":
            kw.default_override = val
        elif tag == "units":
            kw.units = val
        elif tag == "modes":
            if kw.section == GAP_FILE:
                sys.exit("%s:%d: @modes on a .gap keyword; the potential file is "
                         "read the same way in every mode" % (READ_FILES, kw.line))
            kw.modes = [x.strip() for x in re.split(r"[,\s]+", val) if x.strip()]
            for m2 in kw.modes:
                if m2 not in MODES:
                    sys.exit("%s:%d: unknown mode %r (expected %s)"
                             % (READ_FILES, kw.line, m2, "/".join(MODES)))
        elif tag == "needs":
            kw.needs = [x.strip() for x in re.split(r"[,\s]+", val) if x.strip()]
        elif tag == "see":
            kw.see = [x.strip() for x in re.split(r"[,\s]+", val) if x.strip()]
        elif tag == "status":
            if val not in ("ignored", "deprecated"):
                sys.exit("%s:%d: @status must be ignored or deprecated"
                         % (READ_FILES, kw.line))
            kw.status = val
        else:
            sys.exit("%s:%d: unknown tag @%s" % (READ_FILES, kw.line, tag))
    while kw.desc and not kw.desc[-1]:
        kw.desc.pop()


# ------------------------------------------------------------------ joining

def enrich(keywords, types):
    """Attach type/default from types.f90 and work out the side effects.

    `types` maps a derived type name to its components. A keyword's fields
    live in input_parameters unless it named another type -- the .gap ones do.
    """
    by_field = {}
    for kw in keywords:
        if kw.fields:
            by_field.setdefault((kw.section, kw.group if kw.section == GAP_FILE else "",
                                 kw.fields[0].lower()), kw.name)
    for kw in keywords:
        members = types.get(kw.type_name or "input_parameters", {})
        kw.type = ""
        kw.default = ""
        kw.primary = ""
        kw.also = []
        if kw.status == "ignored":
            kw.type = kw.type_override or "any"
            continue
        # The primary field is the one that shares the keyword's name when
        # there is one, otherwise the first assigned.
        low = [f.lower() for f in kw.fields]
        primary = None
        for n in kw.names:
            if n.lower() in low:
                primary = kw.fields[low.index(n.lower())]
                break
        if primary is None and kw.fields:
            primary = kw.fields[0]
        kw.primary = primary or ""
        if primary and primary.lower() in members:
            kw.type, kw.default = members[primary.lower()]
        if kw.type_override:
            kw.type = kw.type_override
        if kw.default_override:
            kw.default = "" if kw.default_override == "none" else kw.default_override
        for f in kw.fields:
            if primary and f.lower() == primary.lower():
                continue
            other = by_field.get((kw.section,
                                  kw.group if kw.section == GAP_FILE else "",
                                  f.lower()))
            kw.also.append((f, other if other and other != kw.name else ""))


def applies(kw, mode):
    return not kw.modes or mode in kw.modes


def grouped(keywords):
    out = OrderedDict()
    for key in list(GROUPS) + sorted({k.group for k in keywords} - set(GROUPS)):
        rows = [k for k in keywords if k.group == key]
        if rows:
            out[key] = rows
    return out


def sectioned(keywords):
    """[(section, OrderedDict(group -> rows))], input file first."""
    out = []
    for section in (INPUT_FILE, GAP_FILE):
        rows = [k for k in keywords if k.section == section]
        if rows:
            out.append((section, grouped(rows)))
    return out


def group_title(key):
    return GROUPS.get(key, (key.replace("_", " ").title(), ""))[0]


def group_blurb(key):
    return GROUPS.get(key, ("", ""))[1]


# ------------------------------------------------------------------ outputs

SECTION_INTRO = {
    INPUT_FILE: [
        "One `keyword = value` per line. A keyword this list does not contain",
        "aborts the run.",
    ],
    GAP_FILE: [
        "The fitted potential, named by `pot_file` in the input file. It is a",
        "sequence of blocks, each opened by `gap_beg <descriptor>` and closed by",
        "`gap_end`, holding one `keyword = value` per line. `gap_beg` and `gap_end`",
        "are file syntax rather than settings and are not listed below. Unlike the",
        "input file, a keyword the block does not recognise is skipped in silence,",
        "so a misspelling here shows up as a wrong answer rather than an error.",
        "These files are written by the fitting code, not by hand.",
    ],
}


def render_md(keywords):
    o = []
    o.append("# TurboGAP keywords")
    o.append("")
    o.append("Generated by `tools/keyword_docs.py` from `src/read_files.f90` and")
    o.append("`src/types.f90`. Do not edit this file: edit the `!>` block above the")
    o.append("keyword in `src/read_files.f90` and regenerate with `make docs`.")
    o.append("")
    o.append("## Contents")
    o.append("")
    for section, groups in sectioned(keywords):
        o.append("**%s** &mdash; %d keywords"
                 % (section[0].upper() + section[1:],
                    sum(len(r) for r in groups.values())))
        o.append("")
        for key, rows in groups.items():
            o.append("- [%s](#%s) (%d)"
                     % (group_title(key), anchor(group_title(key)), len(rows)))
        o.append("")

    for section, groups in sectioned(keywords):
        o.append("# Keywords for %s" % section)
        o.append("")
        o.extend(SECTION_INTRO.get(section, []))
        o.append("")
        for key, rows in groups.items():
            o.append("## " + group_title(key))
            o.append("")
            if group_blurb(key):
                o.append(group_blurb(key))
                o.append("")
            o.append("| Keyword | Type | Default | Units | Modes | Description | Depends on |"
                     if section == INPUT_FILE else
                     "| Keyword | Type | Default | Units | Description | Depends on |")
            o.append("| --- | --- | --- | --- | --- | --- | --- |"
                     if section == INPUT_FILE else
                     "| --- | --- | --- | --- | --- | --- |")
            for kw in sorted(rows, key=lambda k: k.name):
                cells = [md_names(kw), kw.type or "", md_code(kw.default), kw.units or ""]
                if section == INPUT_FILE:
                    cells.append(modes_text(kw))
                cells += [md_desc(kw), md_deps(kw)]
                o.append("| %s |" % " | ".join(cells))
            o.append("")
    return "\n".join(o) + "\n"


def anchor(title):
    return re.sub(r"[^a-z0-9]+", "-", title.lower()).strip("-")


def md_code(text):
    return "`%s`" % text if text else ""


def md_names(kw):
    s = "`%s`" % kw.name
    if kw.aliases:
        s += "<br>" + ", ".join("`%s`" % a for a in kw.aliases)
    return s


def md_desc(kw):
    body = " ".join(kw.desc) if kw.desc else "_undocumented_"
    if kw.status == "ignored":
        body = "**Accepted and ignored.** " + body
    elif kw.status == "deprecated":
        body = "**Deprecated.** " + body
    return body.replace("|", "\\|")


def md_deps(kw):
    bits = []
    for n in kw.needs:
        bits.append("needs `%s`" % n)
    for field, other in kw.also:
        bits.append("sets `%s`" % (other or field))
    for n in kw.see:
        bits.append("see `%s`" % n)
    return "; ".join(bits)


def modes_text(kw):
    return "all" if not kw.modes else ", ".join(kw.modes)


HTML_HEAD = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>TurboGAP input keywords</title>
<style>
/* Headings are set in the monospace face and the prose in a sans one: the
   subject here is a list of literal strings typed into a file, so the face
   the keyword is written in is the one that leads. Neutrals are pulled a few
   degrees towards the verdigris accent so the greys read as chosen. */
:root {
  --bg:      #fafbfa;
  --panel:   #f0f4f3;
  --raised:  #ffffff;
  --fg:      #16211f;
  --muted:   #5b6c69;
  --faint:   #8a9996;
  --rule:    #dbe4e1;
  --rule-hi: #c3d2ce;
  --accent:  #0d6f6a;
  --accent-soft: #e2efed;
  --flag:    #7d5312;
  --flag-bg: #f8efdc;
  --shadow:  0 1px 2px rgba(16, 40, 37, .05);
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --bg:      #0f1514;
    --panel:   #161f1d;
    --raised:  #131b1a;
    --fg:      #e2eae8;
    --muted:   #93a5a1;
    --faint:   #6d807c;
    --rule:    #232f2d;
    --rule-hi: #32403d;
    --accent:  #5cc4bb;
    --accent-soft: #16302e;
    --flag:    #e0b978;
    --flag-bg: #2a2214;
    --shadow:  none;
  }
}
:root[data-theme="dark"] {
  --bg:      #0f1514;
  --panel:   #161f1d;
  --raised:  #131b1a;
  --fg:      #e2eae8;
  --muted:   #93a5a1;
  --faint:   #6d807c;
  --rule:    #232f2d;
  --rule-hi: #32403d;
  --accent:  #5cc4bb;
  --accent-soft: #16302e;
  --flag:    #e0b978;
  --flag-bg: #2a2214;
  --shadow:  none;
}

* { box-sizing: border-box; }

html { scroll-behavior: smooth; scroll-padding-top: 5rem; }
@media (prefers-reduced-motion: reduce) { html { scroll-behavior: auto; } }

body {
  margin: 0;
  background: var(--bg);
  color: var(--fg);
  font: 15px/1.6 ui-sans-serif, system-ui, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
  -webkit-font-smoothing: antialiased;
}

code, .mono, h1, h2, th, .eyebrow {
  font-family: ui-monospace, SFMono-Regular, "SF Mono", Menlo, Consolas, "Liberation Mono", monospace;
}

.page {
  display: grid;
  grid-template-columns: 15rem minmax(0, 1fr);
  gap: 0 3rem;
  max-width: 84rem;
  margin: 0 auto;
  padding: 0 1.5rem 6rem;
  align-items: start;
}
@media (max-width: 60rem) {
  .page { grid-template-columns: minmax(0, 1fr); gap: 0; }
}

/* --- masthead -------------------------------------------------------- */
.masthead { grid-column: 1 / -1; padding: 3rem 0 1.75rem; }
.eyebrow {
  font-size: 11px; letter-spacing: .16em; text-transform: uppercase;
  color: var(--accent); margin: 0 0 .7rem;
}
h1 {
  font-size: clamp(1.6rem, 1.1rem + 1.6vw, 2.3rem);
  font-weight: 600; letter-spacing: -.02em; margin: 0;
  text-wrap: balance;
}
.lede {
  color: var(--muted); margin: .8rem 0 0; max-width: 62ch;
}
.lede code { font-size: .92em; }

/* --- section rail ---------------------------------------------------- */
.rail {
  position: sticky; top: 1.5rem;
  display: flex; flex-direction: column; gap: .1rem;
  padding-bottom: 2rem;
}
.rail a {
  display: flex; justify-content: space-between; gap: .75rem; align-items: baseline;
  text-decoration: none; color: var(--muted);
  font-size: 13px; padding: .3rem .6rem;
  border-left: 2px solid var(--rule); border-radius: 0 4px 4px 0;
}
.rail a:hover { color: var(--fg); background: var(--panel); border-left-color: var(--accent); }
.rail a:focus-visible { outline: 2px solid var(--accent); outline-offset: 2px; }
.rail .n { font-family: ui-monospace, monospace; font-variant-numeric: tabular-nums;
           font-size: 11.5px; color: var(--faint); }
@media (max-width: 60rem) {
  .rail {
    position: static; flex-direction: row; flex-wrap: wrap; gap: .4rem;
    margin-bottom: 1.5rem;
  }
  .rail a { border-left: 0; border: 1px solid var(--rule); border-radius: 999px; padding: .2rem .7rem; }
  .rail .n { display: none; }
}

/* --- filter bar ------------------------------------------------------ */
.filters {
  position: sticky; top: 0; z-index: 5;
  display: flex; gap: .6rem; flex-wrap: wrap; align-items: center;
  padding: .85rem 0;
  background: var(--bg);
  border-bottom: 1px solid var(--rule);
  margin-bottom: .25rem;
}
input[type=search], select {
  font: inherit; font-size: 14px;
  padding: .4rem .65rem;
  color: var(--fg); background: var(--raised);
  border: 1px solid var(--rule-hi); border-radius: 6px;
  box-shadow: var(--shadow);
}
input[type=search] { flex: 1 1 18rem; min-width: 10rem; }
input[type=search]:focus-visible, select:focus-visible {
  outline: 2px solid var(--accent); outline-offset: 1px; border-color: var(--accent);
}
.tally {
  font-family: ui-monospace, monospace; font-variant-numeric: tabular-nums;
  font-size: 12.5px; color: var(--faint); margin-left: auto;
}

/* --- groups ---------------------------------------------------------- */
section { padding-top: 2.4rem; }
h2 {
  font-size: 1rem; font-weight: 600; letter-spacing: -.005em;
  margin: 0; display: flex; align-items: baseline; gap: .6rem;
}
h2::after {
  content: ""; flex: 1; height: 1px; background: var(--rule);
}
.blurb { color: var(--muted); font-size: 13.5px; margin: .45rem 0 1rem; max-width: 60ch; }

.scroll { overflow-x: auto; }
table { border-collapse: collapse; width: 100%; font-size: 13.5px; }
th, td { text-align: left; vertical-align: top; padding: .55rem .7rem; }
th {
  font-size: 10.5px; font-weight: 600; letter-spacing: .09em; text-transform: uppercase;
  color: var(--faint); border-bottom: 1px solid var(--rule-hi); white-space: nowrap;
  padding-bottom: .4rem;
}
tbody tr { border-bottom: 1px solid var(--rule); }
tbody tr:hover { background: var(--panel); }
tbody tr:last-child { border-bottom: 0; }

td.kw { white-space: nowrap; padding-left: .7rem; border-left: 2px solid transparent; }
tr.flagged td.kw { border-left-color: var(--flag); }
td.kw code { font-weight: 600; color: var(--fg); }
td.kw .alias { display: block; margin-top: .15rem; }
td.kw .alias code { font-weight: 400; color: var(--faint); font-size: 12px; }

td.ty, td.df, td.un, td.mo { white-space: nowrap; color: var(--muted); font-size: 12.5px; }
td.df code, td.ty { font-variant-numeric: tabular-nums; }
td.mo { color: var(--faint); }
td.de { min-width: 24rem; }
td.dp { color: var(--muted); font-size: 12.5px; min-width: 10rem; }
td.dp .dep { display: block; }

code {
  font-size: 12.5px;
  background: var(--panel); border-radius: 3px; padding: .08em .32em;
}
td.dp code, td.de code { background: none; padding: 0; color: var(--accent); }
a { color: var(--accent); }

.flag {
  display: inline-block; font-family: ui-monospace, monospace;
  font-size: 10px; letter-spacing: .08em; text-transform: uppercase;
  background: var(--flag-bg); color: var(--flag);
  border-radius: 3px; padding: .1rem .35rem; margin-right: .4rem;
  vertical-align: 1px;
}

.railhead {
  font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
  font-size: 10px; letter-spacing: .14em; text-transform: uppercase;
  color: var(--faint); margin: 1.1rem 0 .35rem .6rem;
}
.rail .railhead:first-child { margin-top: 0; }
@media (max-width: 60rem) { .railhead { display: none; } }

.part { padding-top: 3rem; }
.part > h2 {
  font-size: 1.35rem; letter-spacing: -.015em; margin-bottom: .5rem;
}
.part > h2::after { display: none; }
.part > .intro {
  color: var(--muted); font-size: 13.5px; max-width: 64ch; margin: 0 0 .5rem;
  border-left: 2px solid var(--accent); padding-left: .85rem;
}
.part.hide { display: none; }

tr.hide, section.hide { display: none; }
.empty {
  grid-column: 1 / -1; color: var(--faint); font-size: 14px;
  padding: 3rem 0; text-align: center; display: none;
}
.empty.show { display: block; }

footer {
  grid-column: 1 / -1; margin-top: 3.5rem; padding-top: 1.25rem;
  border-top: 1px solid var(--rule); color: var(--faint); font-size: 12.5px;
}
</style>
"""


def render_html(keywords):
    parts = sectioned(keywords)
    o = [HTML_HEAD, '<div class="page">']

    o.append('<header class="masthead">')
    o.append('<p class="eyebrow">TurboGAP &middot; keyword reference</p>')
    o.append("<h1>Keywords</h1>")
    o.append('<p class="lede">%d keywords across %s. Types and defaults are read out '
             'of <code>src/types.f90</code>, and the rest out of the branch in '
             '<code>src/read_files.f90</code> that parses the keyword, so nothing here '
             'is transcribed by hand.</p>'
             % (len(keywords), " and ".join(esc(s) for s, _ in parts)))
    o.append("</header>")

    o.append('<nav class="rail" aria-label="Groups">')
    for section, groups in parts:
        o.append('<p class="railhead">%s</p>' % esc(section))
        for key, rows in groups.items():
            o.append('<a href="#%s"><span>%s</span><span class="n">%d</span></a>'
                     % (anchor(group_title(key)), esc(group_title(key)), len(rows)))
    o.append("</nav>")

    o.append("<main>")
    o.append('<div class="filters">')
    o.append('<input type="search" id="q" aria-label="Filter keywords" '
             'placeholder="Filter by name or description…" '
             'autocomplete="off" spellcheck="false">')
    o.append('<select id="mode" aria-label="Mode">'
             '<option value="">every mode</option>')
    for m in MODES:
        o.append('<option value="%s">turbogap %s</option>' % (m, m))
    o.append("</select>")
    o.append('<span class="tally" id="tally"></span>')
    o.append("</div>")

    for section, groups in parts:
        gap = section == GAP_FILE
        o.append('<div class="part">')
        o.append("<h2>Keywords for %s</h2>" % esc(section))
        intro = " ".join(SECTION_INTRO.get(section, []))
        intro = re.sub(r"`([^`]+)`", r"<code>\1</code>", esc(intro).replace("&#x27;", "'"))
        if intro:
            o.append('<p class="intro">%s</p>' % intro)
        for key, rows in groups.items():
            o.append('<section id="%s">' % anchor(group_title(key)))
            o.append("<h2>%s</h2>" % esc(group_title(key)))
            if group_blurb(key):
                o.append('<p class="blurb">%s</p>' % esc(group_blurb(key)))
            head = ("<th>Keyword</th><th>Type</th><th>Default</th><th>Units</th>"
                    + ("" if gap else "<th>Modes</th>")
                    + "<th>Description</th><th>Depends on</th>")
            o.append('<div class="scroll"><table><thead><tr>%s</tr></thead><tbody>' % head)
            for kw in sorted(rows, key=lambda k: k.name):
                hay = " ".join([kw.name] + kw.aliases + kw.desc).lower()
                o.append('<tr%s data-h="%s" data-m="%s">' % (
                    ' class="flagged"' if kw.status else "",
                    esc(hay),
                    "" if gap else esc(" ".join(kw.modes) or " ".join(MODES))))
                names = "<code>%s</code>" % esc(kw.name)
                if kw.aliases:
                    names += '<span class="alias">%s</span>' % ", ".join(
                        "<code>%s</code>" % esc(a) for a in kw.aliases)
                o.append('<td class="kw">%s</td>' % names)
                o.append('<td class="ty">%s</td>' % esc(kw.type))
                o.append('<td class="df">%s</td>'
                         % (("<code>%s</code>" % esc(kw.default)) if kw.default else ""))
                o.append('<td class="un">%s</td>' % esc(kw.units))
                if not gap:
                    o.append('<td class="mo">%s</td>' % esc(modes_text(kw)))
                o.append('<td class="de">%s</td>' % html_desc(kw))
                o.append('<td class="dp">%s</td>' % html_deps(kw))
                o.append("</tr>")
            o.append("</tbody></table></div></section>")
        o.append("</div>")
    o.append('<p class="empty" id="empty">Nothing matches that filter.</p>')
    o.append("</main>")

    o.append("<footer>Generated by <code>tools/keyword_docs.py</code>. "
             "To change an entry, edit the <code>!&gt;</code> block above the keyword in "
             "<code>src/read_files.f90</code> and run <code>make docs</code>.</footer>")
    o.append("</div>")
    o.append(HTML_SCRIPT)
    o.append("</body>\n</html>")
    return "\n".join(o) + "\n"


def esc(text):
    return html.escape(text or "", quote=True)


def html_desc(kw):
    body = esc(" ".join(kw.desc)) if kw.desc else '<em>undocumented</em>'
    if kw.status:
        body = '<span class="tag">%s</span>%s' % (kw.status, body)
    return body


def html_deps(kw):
    bits = []
    for n in kw.needs:
        bits.append("needs <code>%s</code>" % esc(n))
    for field, other in kw.also:
        bits.append("sets <code>%s</code>" % esc(other or field))
    for n in kw.see:
        bits.append("see <code>%s</code>" % esc(n))
    return "".join('<span class="dep">%s</span>' % b for b in bits)


HTML_SCRIPT = """<script>
(function () {
  var q = document.getElementById('q'),
      mode = document.getElementById('mode'),
      tally = document.getElementById('tally'),
      empty = document.getElementById('empty'),
      rows = Array.prototype.slice.call(document.querySelectorAll('tbody tr')),
      secs = Array.prototype.slice.call(document.querySelectorAll('section')),
      parts = Array.prototype.slice.call(document.querySelectorAll('.part'));

  function apply() {
    var needle = q.value.trim().toLowerCase(), m = mode.value, shown = 0;
    rows.forEach(function (r) {
      // A row with no modes belongs to the potential file, which the mode
      // filter says nothing about, so leave it visible.
      var modes = r.dataset.m;
      var ok = (!needle || r.dataset.h.indexOf(needle) !== -1) &&
               (!m || !modes || modes.split(' ').indexOf(m) !== -1);
      r.classList.toggle('hide', !ok);
      if (ok) { shown++; }
    });
    secs.forEach(function (s) {
      s.classList.toggle('hide', !s.querySelector('tbody tr:not(.hide)'));
    });
    parts.forEach(function (p) {
      p.classList.toggle('hide', !p.querySelector('tbody tr:not(.hide)'));
    });
    empty.classList.toggle('show', shown === 0);
    tally.textContent = shown === rows.length
      ? rows.length + ' keywords'
      : shown + ' of ' + rows.length;
  }

  q.addEventListener('input', apply);
  mode.addEventListener('change', apply);
  apply();
})();
</script>
"""


# --------------------------------------------------------- generated Fortran

def f_str(text):
    """A Fortran character literal, doubling embedded apostrophes."""
    return "'" + text.replace("'", "''") + "'"


def wrap(text, width):
    words, lines, cur = text.split(), [], ""
    for w in words:
        if cur and len(cur) + 1 + len(w) > width:
            lines.append(cur)
            cur = w
        else:
            cur = (cur + " " + w) if cur else w
    if cur:
        lines.append(cur)
    return lines or [""]


def render_f90(keywords):
    parts = sectioned(keywords)
    o = []
    o.append("! Generated by tools/keyword_docs.py from src/read_files.f90 and")
    o.append("! src/types.f90. Do not edit: edit the !> block above the keyword in")
    o.append("! src/read_files.f90 and run `make docs`.")
    o.append("")
    o.append("module keyword_help")
    o.append("")
    o.append("   implicit none")
    o.append("")
    o.append("   private")
    o.append("   public :: print_keyword_help, keyword_help_topics")
    o.append("")
    o.append("contains")
    o.append("")
    o.append("!  What `turbogap --help <topic>` accepts, for the error message.")
    o.append("   function keyword_help_topics() result(text)")
    o.append("      character(len=64) :: text")
    o.append("      text = " + f_str("/".join(HELP_TOPICS)))
    o.append("   end function keyword_help_topics")
    o.append("")
    o.append("!  Print the keyword reference. An empty topic prints everything; a mode")
    o.append("!  prints the input-file keywords that do something in it; \"gap\" prints")
    o.append("!  the potential-file keywords instead.")
    o.append("   subroutine print_keyword_help(topic)")
    o.append("")
    o.append("      character(len=*), intent(in) :: topic")
    o.append("")
    o.append("      logical :: every")
    o.append("      logical :: gap_only")
    o.append("      character(len=len(topic)) :: mode")
    o.append("")
    o.append("      every = (len_trim(topic) == 0)")
    o.append("      gap_only = (topic == " + f_str("gap") + ")")
    o.append("!     A mode never selects the potential file, so blank it out there and")
    o.append("!     the per-keyword guards below need no second variable.")
    o.append("      mode = topic")
    o.append("      if (gap_only) mode = " + f_str("") )
    o.append("")
    o.append("      write (*, '(A)') " + f_str("TurboGAP keywords"))
    o.append("      write (*, '(A)') " + f_str(""))
    o.append("      if (every) then")
    o.append("         write (*, '(A)') " + f_str(
        "  Everything. `turbogap --help <%s>` narrows this down." % "|".join(HELP_TOPICS)))
    o.append("      else if (gap_only) then")
    o.append("         write (*, '(A)') " + f_str(
        "  Keywords of the potential (.gap) file."))
    o.append("      else")
    o.append("         write (*, '(A)') " + f_str(
        "  Input-file keywords that do something in mode: ") + "//trim(topic)")
    o.append("      end if")
    o.append("      write (*, '(A)') " + f_str(""))
    o.append("      write (*, '(A)') " + f_str(
        "  Full reference, with the cross-references: docs/keywords.html"))
    o.append("      write (*, '(A)') " + f_str(""))

    for section, groups in parts:
        gap = section == GAP_FILE
        o.append("")
        o.append("!  ======== %s" % section)
        o.append("      if (%s) then" % ("every .or. gap_only" if gap
                                         else "every .or. .not. gap_only"))
        o.append("         write (*, '(A)') " + f_str(""))
        o.append("         write (*, '(A)') " + f_str(
            "##### KEYWORDS FOR %s #####" % section.upper()))
        o.append("         write (*, '(A)') " + f_str(""))
        for line in wrap(" ".join(SECTION_INTRO.get(section, []))
                         .replace("`", ""), 72):
            o.append("         write (*, '(A)') " + f_str("  " + line))
        o.append("         write (*, '(A)') " + f_str(""))
        o.append("      end if")

        for key, rows in groups.items():
            rows = sorted(rows, key=lambda k: k.name)
            o.append("")
            o.append("!     ---- %s" % group_title(key))
            o.append("      if (%s) then" % section_guard(section, rows))
            o.append("         write (*, '(A)') "
                     + f_str("=== %s ===" % group_title(key).upper()))
            o.append("         write (*, '(A)') " + f_str(""))
            o.append("      end if")
            for kw in rows:
                emit_keyword(o, kw)

    o.append("")
    o.append("   end subroutine print_keyword_help")
    o.append("")
    o.append("end module keyword_help")
    text = "\n".join(o) + "\n"
    # Nothing downstream would tell us which line was too long, only that the
    # module does not compile, so check here where the offender is nameable.
    for n, line in enumerate(text.split("\n"), 1):
        if len(line) > 130:
            sys.exit("%s:%d: generated line is %d columns, over Fortran's 132:\n  %s"
                     % (os.path.relpath(GEN_F90, ROOT), n, len(line), line[:100]))
    return text


def section_guard(section, rows):
    """Show a .gap group whenever the potential file is shown; show an
    input-file group when its modes overlap the one asked for."""
    if section == GAP_FILE:
        return "every .or. gap_only"
    if any(not k.modes for k in rows):
        return "every .or. .not. gap_only"
    modes = sorted({m for k in rows for m in k.modes})
    return "every .or. (" + " .or. ".join("mode == " + f_str(m) for m in modes) + ")"


def emit_keyword(o, kw):
    if kw.section == GAP_FILE:
        guard = "every .or. gap_only"
    elif not kw.modes:
        guard = "every .or. .not. gap_only"
    else:
        guard = "every .or. (" + " .or. ".join(
            "mode == " + f_str(m) for m in kw.modes) + ")"
    o.append("      if (%s) then" % guard)
    head = "  " + kw.name
    if kw.aliases:
        head += " (or " + ", ".join(kw.aliases) + ")"
    meta = []
    if kw.type:
        meta.append(kw.type)
    if kw.default:
        meta.append("default " + kw.default)
    if kw.units:
        meta.append(kw.units)
    if meta:
        meta = "[" + ", ".join(meta) + "]"
        # Free-form Fortran stops at column 132 and the write statement costs 26
        # of them, so a long name and a long default go on separate lines rather
        # than making the module fail to compile.
        if len(head) + 2 + len(meta) <= 74:
            head = head.ljust(38) + " " + meta
        else:
            o.append("         write (*, '(A)') " + f_str(head))
            head = "      " + meta
    o.append("         write (*, '(A)') " + f_str(head))
    body = " ".join(kw.desc) if kw.desc else "(undocumented)"
    if kw.status == "ignored":
        body = "Accepted and ignored. " + body
    elif kw.status == "deprecated":
        body = "Deprecated. " + body
    for line in wrap(body, 68):
        o.append("         write (*, '(A)') " + f_str("      " + line))
    deps = []
    if kw.needs:
        deps.append("needs " + ", ".join(kw.needs))
    for field, other in kw.also:
        deps.append("sets " + (other or field))
    if kw.see:
        deps.append("see " + ", ".join(kw.see))
    if deps:
        for n, line in enumerate(wrap("; ".join(deps), 62)):
            lead = "      -> " if n == 0 else "         "
            o.append("         write (*, '(A)') " + f_str(lead + line))
    o.append("         write (*, '(A)') " + f_str(""))
    o.append("      end if")


# --------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    ap.add_argument("--check", action="store_true",
                    help="do not write; exit 1 if an output is stale or a keyword is undocumented")
    ap.add_argument("--list-undocumented", action="store_true")
    args = ap.parse_args()

    types = parse_types(TYPES)
    keywords = parse_read_files(READ_FILES) + parse_gap_hypers(READ_FILES)
    enrich(keywords, types)

    missing = [k for k in keywords if not k.documented]
    if args.list_undocumented:
        for k in missing:
            print("%s:%d: %s (%s, %s)"
                  % (READ_FILES, k.line, k.name, k.group, k.section))
        print("%d of %d undocumented" % (len(missing), len(keywords)))
        return 0

    outputs = [(DOC_MD, render_md(keywords)),
               (DOC_HTML, render_html(keywords)),
               (GEN_F90, render_f90(keywords))]

    if args.check:
        bad = 0
        for k in missing:
            sys.stderr.write("%s:%d: keyword '%s' has no !> block\n"
                             % (READ_FILES, k.line, k.name))
            bad += 1
        for path, text in outputs:
            old = open(path).read() if os.path.exists(path) else None
            if old != text:
                sys.stderr.write("%s is out of date; run `make docs`\n"
                                 % os.path.relpath(path, ROOT))
                bad += 1
        if bad:
            return 1
        print("keyword docs up to date (%d keywords, all documented)" % len(keywords))
        return 0

    for path, text in outputs:
        d = os.path.dirname(path)
        if d and not os.path.isdir(d):
            os.makedirs(d)
        with open(path, "w") as fh:
            fh.write(text)
        print("wrote %s" % os.path.relpath(path, ROOT))
    if missing:
        sys.stderr.write("warning: %d of %d keywords have no !> block\n"
                         % (len(missing), len(keywords)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
