#!/usr/bin/env python3
"""Check that every column of a thermo.log sits under its own heading.

thermo.log is a fixed-width table whose header is assembled from string
literals while the values are written with Fortran format specifiers. The two
are in different places and nothing makes them agree, so a column added on one
side and not the other shifts every heading after it -- which is worse than no
heading at all, because it still looks like a label.

The values are right-justified in their fields, so a heading belongs to a
column when the two end in the same character position.

    check_thermo_header.py <thermo.log> [more.log ...]
"""

import sys


def token_ends(line):
    """End position of each whitespace-separated token, 0-based exclusive."""
    ends = []
    in_token = False
    for index, character in enumerate(line):
        if character.isspace():
            if in_token:
                ends.append(index)
                in_token = False
        else:
            in_token = True
    if in_token:
        ends.append(len(line))
    return ends


def check(path):
    lines = [line.rstrip("\n") for line in open(path, encoding="utf-8")]
    header = next((line for line in lines if line.startswith("#")), None)
    if header is None:
        print(f"{path}: no header line")
        return 1
    data = [line for line in lines if line and not line.startswith("#")]
    if not data:
        print(f"{path}: no data rows")
        return 1

    # The leading '#' marks the line as a comment; it is not a heading.
    head_ends = token_ends(" " + header[1:])
    row_ends = token_ends(data[0])

    # Trailing columns may be deliberately unnamed -- the nine lattice columns
    # are, when nothing follows them -- so the headings that exist are checked
    # against the columns they sit over, and a shortfall is reported rather
    # than treated as a misalignment.
    if len(head_ends) > len(row_ends):
        print(f"{path}: {len(head_ends)} headings over only {len(row_ends)} columns")
        return 1
    unnamed = len(row_ends) - len(head_ends)

    bad = [(i, h, r) for i, (h, r) in enumerate(zip(head_ends, row_ends), start=1) if h != r]
    if bad:
        for index, head_end, row_end in bad[:6]:
            print(f"{path}: column {index} heading ends at {head_end}, values end at {row_end}")
        print(f"{path}: {len(bad)} of {len(row_ends)} columns misaligned")
        return 1

    note = f", {unnamed} unnamed" if unnamed else ""
    print(f"{path}: {len(row_ends)} columns, all aligned{note}")
    return 0


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    return max(check(path) for path in sys.argv[1:])


sys.exit(main())
