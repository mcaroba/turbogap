#!/usr/bin/env python3
"""Prove a reformat changed only layout, by comparing token streams.

`diff -w` answers this for Fortran, where the formatter only ever moves
whitespace. It does not answer it for C: clang-format also breaks long lines,
and breaking a #define means ADDING a `\\` continuation -- a real new token in
a file the commit message calls a pure reformat.

So compare tokens, not text. Anything this reports is a genuine change to the
token stream and has to be looked at; anything it does not report is layout.

The clang-format adoption on this branch reported exactly one insertion per
affected .cu file, both of them the `\\` from splitting

    #define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

across two lines. Benign, and worth having confirmed rather than assumed.

Note this is why .clang-format sets ReflowComments: false. With reflow on,
re-wrapping a long // comment inserts and removes comment markers, every file
reports differences, and the check stops being able to tell you anything.

Usage:
    tools/check_reformat_only.py <git-ref> <file> [<file> ...]
    tools/check_reformat_only.py HEAD src/*.cu src/*.cc
"""

import difflib
import re
import subprocess
import sys

TOKEN = re.compile(r"[A-Za-z_][A-Za-z0-9_]*|[0-9]+\.?[0-9]*|\S")


def tokens(text):
    return TOKEN.findall(text)


def main(argv):
    if len(argv) < 3:
        print(__doc__)
        return 2
    ref, paths = argv[1], argv[2:]

    rc = 0
    for path in paths:
        try:
            old = subprocess.run(
                ["git", "show", "%s:%s" % (ref, path)],
                capture_output=True, text=True, check=True,
            ).stdout
        except subprocess.CalledProcessError:
            print("  SKIP    %s (not in %s)" % (path, ref))
            continue

        with open(path) as fh:
            new = fh.read()

        a, b = tokens(old), tokens(new)
        if a == b:
            print("  ok      %-36s %d tokens, layout only" % (path, len(a)))
            continue

        rc = 1
        print("  CHANGED %-36s %d -> %d tokens" % (path, len(a), len(b)))
        sm = difflib.SequenceMatcher(None, a, b, autojunk=False)
        shown = 0
        for tag, i1, i2, j1, j2 in sm.get_opcodes():
            if tag == "equal":
                continue
            shown += 1
            if shown > 10:
                print("      ... more")
                break
            print("      %-8s old %r" % (tag, " ".join(a[i1:i2])[:80]))
            print("               new %r" % (" ".join(b[j1:j2])[:80]))

    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv))
