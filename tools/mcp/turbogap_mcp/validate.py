"""Check a TurboGAP input file the way src/read_files.f90 will read it.

Each line is read list-directed: `read (unit, *) keyword` for the first token,
then `cjunk, cjunk, value(s)`. So the second token is discarded whatever it is,
values end at a blank, a comma or an unquoted slash, a list keeps reading into
the next line until it has its count, and a keyword glued to its `=` is one
unknown token.
"""

import re
from dataclasses import asdict, dataclass

from .keywords import KeywordTable

REQUIRED = ["atoms_file", "pot_file"]

INT_RE = re.compile(r"^[+-]?\d+$")
REAL_RE = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eEdDqQ][+-]?\d+)?$")
SPECIAL_REAL_RE = re.compile(r"^[+-]?(inf|infinity|nan)$", re.I)
LOGICAL_RE = re.compile(r"^\.?[tTfF]")
REPEAT_RE = re.compile(r"^(\d+)\*(.*)$")
COUNT_TERM_RE = re.compile(r"params\s*%\s*(\w+)|\b(n_species)\b")


@dataclass
class Issue:
    line: int
    keyword: str
    severity: str
    message: str


@dataclass
class Entry:
    line: int
    raw_keyword: str
    keyword: str
    separator: str
    values: list
    slash_cut: str


def split_list_directed(text):
    """Tokens of a list-directed record, and the token an unquoted slash cut
    short ("" when no slash cut one). Quotes are kept off the token."""
    tokens = []
    i = 0
    n = len(text)
    while i < n:
        ch = text[i]
        if ch in " \t,":
            i += 1
            continue
        if ch == "/":
            return tokens, ""
        if ch in "'\"":
            quote = ch
            i += 1
            token = ""
            while i < n:
                if text[i] == quote:
                    if i + 1 < n and text[i + 1] == quote:
                        token += quote
                        i += 2
                        continue
                    i += 1
                    break
                token += text[i]
                i += 1
            tokens.append(token)
            continue
        start = i
        while i < n and text[i] not in " \t,/":
            i += 1
        token = text[start:i]
        if i < n and text[i] == "/":
            rest = text[i:].split()[0] if text[i:].split() else ""
            tokens.append(token)
            return tokens, token + rest
        tokens.append(token)
    return tokens, ""


def expand_repeats(tokens):
    out = []
    for token in tokens:
        m = REPEAT_RE.match(token)
        if m and m.group(2):
            out.extend([m.group(2)] * int(m.group(1)))
        else:
            out.append(token)
    return out


def parse_entries(text):
    """(entries, issues) for every non-comment line."""
    entries = []
    issues = []
    for number, line in enumerate(text.splitlines(), start=1):
        tokens, slash_cut = split_list_directed(line)
        if not tokens:
            continue
        first = tokens[0]
        if first[:1] in ("#", "!"):
            continue
        if "=" in first and first != "=":
            issues.append(Issue(number, first, "error",
                                f"'{first}' is read as one token, so TurboGAP aborts on an unknown "
                                f"keyword; write it as '{first.split('=')[0]} = {first.split('=', 1)[1]}'"))
            continue
        separator = tokens[1] if len(tokens) > 1 else ""
        entries.append(Entry(number, first, first.lower(), separator,
                             expand_repeats(tokens[2:]), slash_cut))
    return entries, issues


def significant(values):
    """Values up to an inline comment. Fortran never sees the comment unless it
    is still short of values, which check_values reports on its own."""
    for i, token in enumerate(values):
        if token[:1] in ("!", "#"):
            return values[:i]
    return values


def base_type(type_name):
    for base in ("integer", "real", "logical", "string"):
        if type_name.startswith(base):
            return base
    return ""


def value_ok(base, token):
    if base == "integer":
        return bool(INT_RE.match(token))
    if base == "real":
        return bool(REAL_RE.match(token) or SPECIAL_REAL_RE.match(token))
    if base == "logical":
        return bool(LOGICAL_RE.match(token))
    return True


def integer_setting(table, entries, field):
    """Value of the integer keyword writing params%<field>, from the input or
    its default; None when neither gives one."""
    owner = table.field_owner(field)
    if owner is None:
        return None
    value = None
    for entry in entries:
        if table.canonical(entry.keyword) == owner.name and entry.values:
            value = entry.values[0]
    if value is None:
        value = owner.default
    if value and INT_RE.match(value):
        return int(value)
    return None


def n_species_setting(entries):
    value = 1
    for entry in entries:
        if entry.raw_keyword == "n_species" and entry.values and INT_RE.match(entry.values[0]):
            return int(entry.values[0])
    return value


def evaluate_count(table, entries, expr):
    def term(m):
        if m.group(2):
            return str(n_species_setting(entries))
        value = integer_setting(table, entries, m.group(1))
        if value is None:
            raise LookupError(m.group(1))
        return str(value)

    try:
        arithmetic = COUNT_TERM_RE.sub(term, expr)
    except LookupError:
        return None
    if not re.match(r"^[\d\s*+-]+$", arithmetic):
        return None
    return eval(arithmetic)


def check_known(table, entry):
    if table.is_known(entry.keyword):
        return []
    hint = ""
    close = table.close_names(entry.keyword)
    if close:
        hint = "; did you mean " + ", ".join(close) + "?"
    return [Issue(entry.line, entry.raw_keyword, "error",
                  f"unknown keyword: TurboGAP aborts on it{hint}")]


def check_shape(entry):
    if entry.separator == "":
        return [Issue(entry.line, entry.raw_keyword, "error",
                      "no '= value': the read carries on into the next line")]
    if entry.separator != "=":
        return [Issue(entry.line, entry.raw_keyword, "error",
                      f"expected '=' after the keyword, found '{entry.separator}'; "
                      "the second token is discarded whatever it is, so the value is misread")]
    if not entry.values:
        return [Issue(entry.line, entry.raw_keyword, "error",
                      "no value: the read carries on into the next line")]
    return []


def check_prepass_case(table, entry):
    if entry.keyword in table.prepass and entry.raw_keyword != entry.keyword:
        return [Issue(entry.line, entry.raw_keyword, "error",
                      f"{entry.keyword} is read in a first pass that does not lower-case the "
                      "keyword, so this spelling is skipped there; write it in lower case")]
    return []


def check_slash(entry):
    if entry.slash_cut:
        kept = entry.values[-1] if entry.values else ""
        return [Issue(entry.line, entry.raw_keyword, "error",
                      f"an unquoted '/' ends a list-directed record, so '{entry.slash_cut}' "
                      f"is read as '{kept}'; quote the value")]
    return []


def check_values(table, entries, entry):
    type_name = table.type_of(entry.keyword)
    base = base_type(type_name)
    name = table.canonical(entry.keyword)
    issues = []
    expr = table.counts.get(name)
    want = 1
    if expr:
        want = evaluate_count(table, entries, expr)
    elif "list" in type_name or "(" in type_name:
        want = None
    if want is not None and len(entry.values) < want:
        issues.append(Issue(entry.line, entry.raw_keyword, "error",
                            f"needs {want} value(s) ({expr or 'scalar'}), has {len(entry.values)}; "
                            "Fortran keeps reading into the next line for the rest"))
    checked = significant(entry.values) if want is None else entry.values[:want]
    for token in checked:
        if not value_ok(base, token):
            issues.append(Issue(entry.line, entry.raw_keyword, "error",
                                f"'{token}' is not a valid {base}"))
    if base == "string" and want == 1 and len(significant(entry.values)) > 1:
        issues.append(Issue(entry.line, entry.raw_keyword, "warning",
                            f"only '{entry.values[0]}' is read; quote the value if it has blanks or commas"))
    return issues


def check_mode(table, entry, mode):
    kw = table.by_name.get(entry.keyword)
    if kw is None or not kw.modes or mode in kw.modes:
        return []
    return [Issue(entry.line, entry.raw_keyword, "warning",
                  f"documented for mode(s) {', '.join(kw.modes)}, not {mode}")]


def check_status(table, entry):
    kw = table.by_name.get(entry.keyword)
    if kw is None or not kw.status:
        return []
    return [Issue(entry.line, entry.raw_keyword, "warning", f"this keyword is {kw.status}")]


def check_needs(table, entry, present):
    kw = table.by_name.get(entry.keyword)
    if kw is None:
        return []
    missing = [n for n in kw.needs if table.canonical(n) not in present]
    if not missing:
        return []
    return [Issue(entry.line, entry.raw_keyword, "warning",
                  f"needs {', '.join(missing)}, which the input does not set")]


def check_duplicates(table, entries):
    issues = []
    first_line = {}
    for entry in entries:
        name = table.canonical(entry.keyword)
        if not table.is_known(entry.keyword):
            continue
        if name in first_line:
            winner = "the first" if name in table.prepass else "this later"
            issues.append(Issue(entry.line, entry.raw_keyword, "warning",
                                f"also set on line {first_line[name]}; {winner} value is used"))
        else:
            first_line[name] = entry.line
    return issues


def check_required(present):
    return [Issue(0, name, "error", f"{name} is required") for name in REQUIRED if name not in present]


def validate(table: KeywordTable, text, mode):
    if mode not in table.modes:
        return [Issue(0, "", "error", f"unknown mode '{mode}'; expected one of {', '.join(table.modes)}")]
    entries, issues = parse_entries(text)
    present = {table.canonical(e.keyword) for e in entries if table.is_known(e.keyword)}
    for entry in entries:
        found = check_known(table, entry)
        if found:
            issues += found
            continue
        issues += check_prepass_case(table, entry)
        shape = check_shape(entry)
        if shape:
            issues += shape
            continue
        issues += check_slash(entry)
        issues += check_values(table, entries, entry)
        issues += check_mode(table, entry, mode)
        issues += check_status(table, entry)
        issues += check_needs(table, entry, present)
    issues += check_duplicates(table, entries)
    issues += check_required(present)
    issues.sort(key=lambda i: (i.line, i.severity != "error"))
    return issues


def as_dicts(issues):
    return [asdict(i) for i in issues]
