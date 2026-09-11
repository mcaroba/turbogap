#!/usr/bin/env python3
"""Compare two output files numerically, within a tolerance.

The bit-exact diff the regression suite uses for the CPU binary cannot be used
for a device build: the SOAP batch decomposition and the cuBLAS reductions sum
the same terms in a different order, so the last digits move on a run that is
entirely correct. This compares number against number instead, and still
requires every non-numeric token -- species names, keywords, the key of a
key=value pair -- to match exactly, so a structural change is still a failure.

An extended-xyz comment line is a run of key=value pairs, some of them quoted
and split across whitespace, so a token is reduced to its numeric part before
being compared.

  compare_tol.py <reference> <test> [--rtol R] [--atol A] [--max-report N]

Exit 0 when every difference is within tolerance.
"""

import argparse
import re
import sys

# Fortran writes 1.0D+00 as readily as 1.0E+00.
EXPONENT = re.compile(r"([0-9])[dD]([-+]?[0-9])")


def split_token(token):
    """Return (label, number) for a token, with number None if there is none.

    'energy=-47309.16' is a labelled number; 'virial="18177.89' is one too,
    with a quote left over from the opening of a quoted list.
    """
    label, _, value = token.rpartition("=")
    value = value.strip('"')
    try:
        return label, float(EXPONENT.sub(r"\1e\2", value))
    except ValueError:
        return token, None


def deviations(reference, test):
    """Yield (column, ref, test, absolute, relative) for each differing token.

    absolute and relative are None where the tokens are not numbers.
    """
    ref_tokens, test_tokens = reference.split(), test.split()
    if len(ref_tokens) != len(test_tokens):
        yield (0, f"{len(ref_tokens)} tokens", f"{len(test_tokens)} tokens", None, None)
        return
    for column, (one, two) in enumerate(zip(ref_tokens, test_tokens), start=1):
        if one == two:
            continue
        label_a, a = split_token(one)
        label_b, b = split_token(two)
        if a is None or b is None or label_a != label_b:
            yield (column, one, two, None, None)
            continue
        absolute = abs(b - a)
        yield (column, one, two, absolute, absolute / abs(a) if a else float("inf"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("reference")
    parser.add_argument("test")
    parser.add_argument("--rtol", type=float, default=1e-6)
    parser.add_argument("--atol", type=float, default=1e-6)
    parser.add_argument("--max-report", type=int, default=8)
    args = parser.parse_args()

    reference = open(args.reference, encoding="utf-8", errors="replace").read().splitlines()
    test = open(args.test, encoding="utf-8", errors="replace").read().splitlines()

    if len(reference) != len(test):
        print(f"      line count differs: reference {len(reference)}, test {len(test)}")
        return 1

    reported, worst_abs, worst_rel, total = 0, 0.0, 0.0, 0
    for number, (one, two) in enumerate(zip(reference, test), start=1):
        for column, a, b, absolute, relative in deviations(one, two):
            if absolute is not None and absolute <= args.atol + args.rtol * abs(float(split_token(a)[1])):
                continue
            total += 1
            if absolute is not None:
                worst_abs = max(worst_abs, absolute)
                worst_rel = max(worst_rel, relative)
            if reported < args.max_report:
                reported += 1
                shown = "not comparable" if absolute is None else f"abs {absolute:.3e}, rel {relative:.3e}"
                print(f"      line {number} col {column}: {a} vs {b}  ({shown})")
    if total:
        print(f"      {total} value(s) outside rtol={args.rtol:g} atol={args.atol:g}; "
              f"worst absolute {worst_abs:.3e}, worst relative {worst_rel:.3e}")
        return 1
    return 0


sys.exit(main())
