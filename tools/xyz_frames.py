#!/usr/bin/env python3
"""Read an extended-XYZ trajectory, honouring its Properties line.

The column layout is not fixed: whether velocities, local energies, charges or
fix_atoms appear depends on the deck, so a reader that counts columns from the
left gets a different quantity for a different input and says nothing about it.
This resolves each field by name.
"""

import re


def parse_properties(comment):
    """[(name, type, count)] from the Properties=... entry of a comment line."""
    match = re.search(r'Properties=(\S+)', comment)
    if not match:
        raise ValueError("no Properties entry in the comment line")
    parts = match.group(1).split(":")
    out = []
    for i in range(0, len(parts) - 2, 3):
        out.append((parts[i], parts[i + 1], int(parts[i + 2])))
    return out


def column_of(properties, name):
    """First column index of a named field, 0-based."""
    index = 0
    for field, _, count in properties:
        if field == name:
            return index
        index += count
    raise KeyError(f"{name} is not in this trajectory: {[p[0] for p in properties]}")


def scalar(comment, key):
    match = re.search(rf'\b{key}=(-?[\d.EeDd+-]+)', comment)
    if not match:
        raise KeyError(f"{key} is not in the comment line")
    return float(match.group(1).replace("D", "E").replace("d", "e"))


def frames(path):
    """Yield (comment, properties, rows) for each frame."""
    with open(path, encoding="utf-8") as handle:
        lines = handle.read().splitlines()
    index = 0
    while index < len(lines):
        if not lines[index].strip():
            index += 1
            continue
        count = int(lines[index].split()[0])
        comment = lines[index + 1]
        rows = [line.split() for line in lines[index + 2:index + 2 + count]]
        yield comment, parse_properties(comment), rows
        index += 2 + count


def max_force(path):
    """Largest force component in the final frame, and that frame's energy."""
    last = None
    for comment, properties, rows in frames(path):
        last = (comment, properties, rows)
    if last is None:
        raise ValueError(f"{path} holds no frames")
    comment, properties, rows = last
    start = column_of(properties, "forces")
    biggest = 0.0
    for row in rows:
        for value in row[start:start + 3]:
            biggest = max(biggest, abs(float(value)))
    return biggest, scalar(comment, "energy")


def energies(path):
    return [scalar(comment, "energy") for comment, _, _ in frames(path)]
