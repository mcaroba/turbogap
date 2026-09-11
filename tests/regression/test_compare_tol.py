#!/usr/bin/env python3
"""Self-test for compare_tol.py, the comparator the --gpu suite decides with.

A comparator that silently passes everything turns the device suite into a
formality, and nothing else would notice. These cases are the ones that
distinguish "the last digit moved" from "the answer changed".

    python3 tests/regression/test_compare_tol.py
"""

import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.join(HERE, "compare_tol.py")

XYZ_HEAD = "2\nProperties=species:S:1:pos:R:3 energy=-100.00000000\n"


def run(reference, test, *flags):
    with tempfile.TemporaryDirectory() as directory:
        paths = []
        for name, text in (("ref", reference), ("test", test)):
            path = os.path.join(directory, name)
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(text)
            paths.append(path)
        done = subprocess.run([sys.executable, TOOL, *paths, *flags],
                              capture_output=True, text=True)
        return done.returncode, done.stdout


CASES = [
    ("identical files pass",
     "C 1.0 2.0 3.0\n", "C 1.0 2.0 3.0\n", (), 0),

    ("a difference below the tolerance passes",
     "C 1.00000000 2.0 3.0\n", "C 1.00000005 2.0 3.0\n", (), 0),

    ("a difference above the tolerance fails",
     "C 1.00000000 2.0 3.0\n", "C 1.00100000 2.0 3.0\n", (), 1),

    # The absolute term is what carries a force component near zero, where any
    # relative tolerance is meaningless.
    ("a tiny absolute difference on a near-zero value passes",
     "C 0.00000001 2.0 3.0\n", "C 0.00000002 2.0 3.0\n", (), 0),
    ("...and fails once the absolute tolerance is tightened",
     "C 0.00000001 2.0 3.0\n", "C 0.00000002 2.0 3.0\n",
     ("--rtol", "0", "--atol", "0"), 1),

    # A structural change has to fail whatever the tolerance is.
    ("a changed species fails at any tolerance",
     "C 1.0 2.0 3.0\n", "O 1.0 2.0 3.0\n", ("--rtol", "1e9", "--atol", "1e9"), 1),
    ("a changed Properties string fails",
     XYZ_HEAD, XYZ_HEAD.replace("pos:R:3", "pos:R:4"), ("--rtol", "1e9"), 1),
    ("a different number of columns fails",
     "C 1.0 2.0 3.0\n", "C 1.0 2.0\n", ("--rtol", "1e9", "--atol", "1e9"), 1),
    ("a different number of lines fails",
     "C 1.0 2.0 3.0\n", "C 1.0 2.0 3.0\nC 1.0 2.0 3.0\n", ("--rtol", "1e9"), 1),

    # key=value is how an extended-xyz comment line carries its numbers.
    ("a key=value number is compared as a number",
     "energy=-100.00000000\n", "energy=-100.00000005\n", (), 0),
    ("...and its key still has to match",
     "energy=-100.0\n", "enthalpy=-100.0\n", ("--rtol", "1e9"), 1),
    ("a quoted list entry is compared as a number",
     'virial="1.00000000 2.0\n', 'virial="1.00000005 2.0\n', (), 0),

    # Fortran writes both spellings of an exponent.
    ("D and E exponents compare equal",
     "1.5D+02\n", "1.5E+02\n", ("--rtol", "0", "--atol", "0"), 0),
]


def main():
    failures = 0
    for name, reference, test, flags, expected in CASES:
        code, output = run(reference, test, *flags)
        if code == expected:
            print(f"  ok    {name}")
        else:
            failures += 1
            print(f"  FAIL  {name}: expected exit {expected}, got {code}")
            print("        " + output.strip().replace("\n", "\n        "))
    print(f"\n{len(CASES) - failures}/{len(CASES)} passed")
    return 1 if failures else 0


sys.exit(main())
