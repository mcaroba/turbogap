"""Compact summaries of TurboGAP output. The files stay on the host; only one
frame of a trajectory, or the decimated thermo.log columns, comes back."""

import math
import re
import shlex
from typing import Any

MAX_ARRAY_ATOMS = 1000
MAX_SERIES_POINTS = 200

# Frames can change size (MC insertions), so frames are walked, not counted.
FRAME_AWK = r"""
left == 0 { n = $1 + 0; left = n + 2; frame++; if (want < 0 || frame == want) buf = "" }
{ if (want < 0 || frame == want) buf = buf $0 "\n"; left-- }
want > 0 && frame > want { exit }
END { print frame; printf "%s", buf }
"""


def frame_script(path, frame):
    """Bash that prints the number of frames walked, then the frame asked for:
    frame is 0-based, negative counts from the end (only -1 is supported)."""
    want = -1 if frame < 0 else frame + 1
    return f"awk -v want={want} {shlex.quote(FRAME_AWK)} {shlex.quote(path)}\n"


def parse_header(line):
    """key=value pairs of an extxyz comment line; quoted values kept whole."""
    out = {}
    for key, quoted, bare in re.findall(r'(\w+)=(?:"([^"]*)"|(\S+))', line):
        out[key] = quoted if quoted else bare
    return out


def columns(properties):
    """{name: (first column, width)} from an extxyz Properties string."""
    parts = properties.split(":")
    out = {}
    start = 0
    for i in range(0, len(parts) - 2, 3):
        width = int(parts[i + 2])
        out[parts[i]] = (start, width)
        start += width
    return out


def numbers(text):
    return [float(x) for x in text.split()]


def summarise_frame(text, arrays=False):
    lines = text.splitlines()
    if len(lines) < 2:
        return {"error": "no frame"}
    n = int(lines[0].split()[0])
    header = parse_header(lines[1])
    cols = columns(header.get("Properties", "species:S:1:pos:R:3"))
    atoms = [line.split() for line in lines[2:2 + n]]
    out: dict[str, Any] = {"n_atoms": n}
    if "species" in cols:
        first = cols["species"][0]
        counts = {}
        for a in atoms:
            counts[a[first]] = counts.get(a[first], 0) + 1
        out["species"] = counts
    for key in ("energy", "time", "step"):
        if key in header:
            out[key] = float(header[key])
    for key in ("virial", "stress", "Lattice"):
        if key in header:
            out[key.lower()] = numbers(header[key])
    if "forces" in cols:
        first, width = cols["forces"]
        forces = [[float(x) for x in a[first:first + width]] for a in atoms]
        norms = [math.sqrt(sum(c * c for c in f)) for f in forces]
        out["force_max"] = max(norms) if norms else 0.0
        out["force_rms"] = math.sqrt(sum(x * x for x in norms) / len(norms)) if norms else 0.0
        out["force_sum"] = [sum(f[k] for f in forces) for k in range(3)]
        if arrays and n <= MAX_ARRAY_ATOMS:
            out["forces"] = forces
    if arrays and n <= MAX_ARRAY_ATOMS and "pos" in cols:
        first, width = cols["pos"]
        out["positions"] = [[float(x) for x in a[first:first + width]] for a in atoms]
    if arrays and n > MAX_ARRAY_ATOMS:
        out["arrays_omitted"] = f"{n} atoms is over {MAX_ARRAY_ATOMS}; read the file on the host"
    return out


def thermo_names(header, n_columns):
    """Column names. The lattice block is unnamed when nothing follows it, so
    beyond the named columns the rest are numbered."""
    names = header.lstrip("#").split()
    if len(names) == n_columns:
        return names
    named = names[:names.index("Pressure") + 1] if "Pressure" in names else names[:n_columns]
    return named + [f"col_{i + 1}" for i in range(len(named), n_columns)]


def summarise_thermo(text):
    lines = [line for line in text.splitlines() if line.strip()]
    header = next((line for line in lines if line.startswith("#")), "")
    rows = [numbers(line) for line in lines if not line.startswith("#")]
    if not rows:
        return {"rows": 0}
    names = thermo_names(header, len(rows[-1]))
    stride = max(1, math.ceil(len(rows) / MAX_SERIES_POINTS))
    kept = rows[::stride]
    if kept[-1] is not rows[-1]:
        kept.append(rows[-1])
    series = {name: [r[i] for r in kept if i < len(r)] for i, name in enumerate(names)}
    return {"rows": len(rows), "stride": stride, "last": dict(zip(names, rows[-1])), "series": series}
