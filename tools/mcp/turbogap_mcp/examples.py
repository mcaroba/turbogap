"""The regression cases as worked examples: each one's input and the recipe to run it."""

import re

CONF_RE = re.compile(r"""^([A-Z_][A-Z0-9_]*)=(?:"([^"]*)"|'([^']*)'|(\S*))""")


def cases_dir(root):
    return root / "tests" / "regression" / "cases"


def parse_case_conf(text):
    """KEY=value lines and the leading comment. case.conf is meant to be sourced,
    but every key in it is a plain assignment, so it is read, never run."""
    values = {}
    about = []
    for line in text.splitlines():
        if line.startswith("#") and not values:
            about.append(line.lstrip("#").strip())
            continue
        m = CONF_RE.match(line.strip())
        if m:
            values[m.group(1)] = next(g for g in m.groups()[1:] if g is not None)
    return values, " ".join(a for a in about if a)


def list_cases(root):
    out = []
    for conf in sorted(cases_dir(root).glob("*/case.conf")):
        values, about = parse_case_conf(conf.read_text())
        out.append({"name": conf.parent.name, "mode": values.get("MODE", ""),
                    "ranks": int(values.get("RANKS", "1") or 1), "about": about[:160]})
    return out


def read_case(root, name):
    directory = cases_dir(root) / name
    if not re.match(r"^[\w.-]+$", name) or not (directory / "case.conf").is_file():
        raise ValueError(f"no regression case named '{name}'")
    values, about = parse_case_conf((directory / "case.conf").read_text())
    links = []
    for item in values.get("LINKS", "").split():
        source, _, target = item.partition(":")
        links.append({"source": source, "target": target or source})
    return {
        "name": name,
        "about": about,
        "mode": values.get("MODE", ""),
        "ranks": int(values.get("RANKS", "1") or 1),
        "data": values.get("DATA", ""),
        "links": links,
        "needs": values.get("NEEDS", "").split(),
        "outputs": values.get("OUTPUTS", "").split(),
        "reference": values.get("REFERENCE", "baseline"),
        "input": (directory / "input").read_text(),
        "note": "links are relative to the case's DATA directory in the turbogap_tests repository",
    }
