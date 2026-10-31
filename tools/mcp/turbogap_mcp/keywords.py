"""The keyword table, as tools/keyword_docs.py parses it out of src/read_files.f90.

Nothing here restates a keyword: names, types, defaults, modes and prose all come
from the parser, so this cannot drift from the Fortran.
"""

import difflib
import importlib.util
import re

INPUT = "input"
GAP = "gap"

# Read by turbogap_setup.f90 before the main pass, which skips them; that pass
# lists them on its comment-skipping line, and nothing else documents them.
PREPASS_TYPES = {"pot_file": "string", "n_species": "integer"}
PREPASS_DOC = {
    "pot_file": "Path to the potential (.gap) file, or none. Read before the main pass; "
                "the first occurrence wins and the keyword must be written in lower case.",
    "n_species": "Number of species; every per-species list is read this many values wide. "
                 "Read before the main pass, defaults to 1, and must be written in lower case.",
}

SKIP_LINE_RE = re.compile(r"keyword_first\s*==\s*'#'.*?\bcycle\b", re.S)
SKIPPED_KW_RE = re.compile(r"keyword\s*==\s*'(\w+)'")
COUNT_RES = [
    re.compile(r"read_parameters\s*\(\s*unit\s*,\s*iostatus\s*,\s*([^,]+?)\s*,"),
    re.compile(r"cjunk\s*,\s*cjunk\s*,\s*\(\s*params[^=]*?,\s*\w+\s*=\s*1\s*,\s*([^)]+?)\s*\)"),
    re.compile(r"cjunk\s*,\s*cjunk\s*,\s*params\s*%\s*\w+\s*\(\s*1\s*:\s*([^)]+?)\s*\)"),
]


def load_keyword_docs(root):
    spec = importlib.util.spec_from_file_location(
        "keyword_docs", root / "tools" / "keyword_docs.py")
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {root}/tools/keyword_docs.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def source_with_line_starts(path):
    """The text of a source file and the offset at which each of its lines starts."""
    lines = path.read_text().split("\n")
    starts = []
    pos = 0
    for line in lines:
        starts.append(pos)
        pos += len(line) + 1
    text = "\n".join(lines)
    return text, starts


def prepass_names(text):
    m = SKIP_LINE_RE.search(text)
    if not m:
        return []
    return SKIPPED_KW_RE.findall(m.group(0))


def count_expression(branch):
    flat = re.sub(r"&\s*\n\s*&?", "", branch)
    for pattern in COUNT_RES:
        m = pattern.search(flat)
        if m:
            return m.group(1).strip()
    return ""


def fixed_count(type_name):
    m = re.match(r"\w+\((\d+)\)$", type_name)
    if m:
        return m.group(1)
    return ""


class KeywordTable:
    def __init__(self, root):
        self.root = root
        kd = load_keyword_docs(root)
        self.modes = list(kd.MODES)
        types = kd.parse_types(kd.TYPES)
        keywords = kd.parse_read_files(kd.READ_FILES) + kd.parse_gap_hypers(kd.READ_FILES)
        kd.enrich(keywords, types)
        self.input_keywords = [k for k in keywords if k.section == kd.INPUT_FILE]
        self.gap_keywords = [k for k in keywords if k.section == kd.GAP_FILE]
        text, starts = source_with_line_starts(root / "src" / "read_files.f90")
        self.counts = self._count_expressions(text, starts)
        self.prepass = prepass_names(text)
        self.by_name = {}
        for kw in self.input_keywords:
            for name in kw.names:
                self.by_name[name.lower()] = kw

    def _count_expressions(self, text, starts):
        lines = sorted(k.line for k in self.input_keywords)
        out = {}
        for kw in self.input_keywords:
            later = [n for n in lines if n > kw.line]
            end = starts[later[0] - 1] if later else len(text)
            expr = count_expression(text[starts[kw.line - 1]:end]) or fixed_count(kw.type)
            if expr:
                out[kw.name] = expr
        return out

    def is_known(self, name):
        return name in self.by_name or name in self.prepass

    def canonical(self, name):
        if name in self.by_name:
            return self.by_name[name].name
        return name

    def type_of(self, name):
        if name in self.by_name:
            return self.by_name[name].type
        return PREPASS_TYPES.get(name, "any")

    def field_owner(self, field):
        """The input keyword whose primary field is params%<field>."""
        for kw in self.input_keywords:
            if kw.primary and kw.primary.lower() == field.lower():
                return kw
        return None

    def describe(self, kw, file=INPUT):
        out = {
            "name": kw.name,
            "aliases": list(kw.aliases),
            "file": file,
            "group": kw.group,
            "type": kw.type,
            "default": kw.default,
            "units": kw.units,
            "modes": list(kw.modes) or ["all"],
            "needs": list(kw.needs),
            "see": list(kw.see),
            "status": kw.status,
            "description": " ".join(kw.desc),
        }
        if file == INPUT and kw.name in self.counts:
            out["count"] = self.counts[kw.name]
        return out

    def describe_prepass(self, name):
        return {
            "name": name, "aliases": [], "file": INPUT, "group": "general",
            "type": PREPASS_TYPES.get(name, "any"), "default": "1" if name == "n_species" else "",
            "units": "", "modes": ["all"], "needs": [], "see": [], "status": "",
            "description": PREPASS_DOC.get(name, ""),
        }

    def lookup(self, name, file=INPUT):
        name = name.strip().lower()
        if file == GAP:
            return [self.describe(k, GAP) for k in self.gap_keywords if name in k.names]
        if name in self.by_name:
            return [self.describe(self.by_name[name])]
        if name in self.prepass:
            return [self.describe_prepass(name)]
        return []

    def search(self, query="", mode="", file=INPUT, limit=25):
        """Name matches first, then description matches; each as a short row."""
        words = query.lower().split()
        pool = self.gap_keywords if file == GAP else self.input_keywords
        rows = []
        for kw in pool:
            if mode and file == INPUT and kw.modes and mode not in kw.modes:
                continue
            names = " ".join(kw.names).lower()
            desc = " ".join(kw.desc).lower()
            if not all(w in names or w in desc for w in words):
                continue
            score = sum(w in names for w in words)
            rows.append((-score, kw.name, kw))
        rows.sort(key=lambda r: (r[0], r[1]))
        out = []
        for _, _, kw in rows[:limit]:
            summary = " ".join(kw.desc)
            out.append({
                "name": kw.name,
                "group": kw.group,
                "type": kw.type,
                "default": kw.default,
                "summary": summary[:160] + ("..." if len(summary) > 160 else ""),
            })
        return {"total": len(rows), "shown": len(out), "keywords": out}

    def close_names(self, name, n=3):
        pool = list(self.by_name) + list(self.prepass)
        return difflib.get_close_matches(name, pool, n=n, cutoff=0.7)
