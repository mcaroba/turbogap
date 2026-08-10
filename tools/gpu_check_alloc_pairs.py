#!/usr/bin/env python3
"""Check that every device buffer is released by the matching deallocator.

The wrappers in src/gpu/gpu_memory.cu are two distinct allocation families, and
mixing them is undefined behaviour rather than a leak:

    gpu_malloc_async(a_d, n, stream)    -> hipMallocAsync   (stream-ordered pool)
    gpu_malloc_all_blocking(a_d, n)     -> hipMalloc        (ordinary device heap)

    gpu_free_async(a_d, stream)         -> hipFreeAsync
    gpu_free(a_d)                       -> hipFree

A pointer from the pool allocator handed to hipFree is not merely untidy: the
allocation belongs to a memory pool the driver still owns. This is a defect
that has already been found twice in this tree by hand -- virial_2b_d and
virial_core_pot_d were released with plain gpu_free despite being
hipMallocAsync allocations, with the intended fix half-written in a trailing
comment beside them (gpu_fixes_handoff.md 6e, and two further sites found
during the gap_backend extraction).

Neither a compiler nor compute-sanitizer reliably reports this: the free
usually "works", and the corruption surfaces later somewhere unrelated. It is
a static property of the source, so check it statically.

What this reports, in order of severity:

  MISMATCH   allocated by one family, freed by the other.
  UNFREED    allocated and never passed to any deallocator in the tree.
  FREE-ONLY  freed but never allocated here (often fine -- the buffer is
             allocated in .cu code or passed in -- but worth an eye).

Exit status is 1 if any MISMATCH is found, so it can gate a build.

Usage:
    tools/gpu_check_alloc_pairs.py [src ...]        # default: src/*.f90
    tools/gpu_check_alloc_pairs.py --quiet          # mismatches only
"""

import glob
import os
import re
import sys
from collections import defaultdict

POOL_ALLOC = ("gpu_malloc_async", "gpu_malloc_neighbors")
PLAIN_ALLOC = ("gpu_malloc_all_blocking",)
POOL_FREE = ("gpu_free_async", "gpu_free_neighbors")
PLAIN_FREE = ("gpu_free",)

FAMILY = {}
for n in POOL_ALLOC:
    FAMILY[n] = ("alloc", "pool")
for n in PLAIN_ALLOC:
    FAMILY[n] = ("alloc", "plain")
for n in POOL_FREE:
    FAMILY[n] = ("free", "pool")
for n in PLAIN_FREE:
    FAMILY[n] = ("free", "plain")

# `call gpu_free(x)` and `call gpu_free_async(x, s)` -- longest name first, or
# gpu_free matches the prefix of gpu_free_async and every async free is
# misreported as a plain one. This is the same shape as the time_mpi /
# time_mpi_positions ordering trap in tools/bundle_timers.py.
NAMES = sorted(FAMILY, key=len, reverse=True)
CALL = re.compile(r"\bcall\s+(%s)\s*\(" % "|".join(NAMES), re.I)


def first_arg(text, open_paren):
    """The first actual argument, from the '(' at open_paren.

    A regex cannot do this: the argument may itself contain parentheses
    (j2_index_d(n_dim_idx)) and a lazy character class either stops inside it
    or swallows the closing paren of the call. Swallowing it is the worse
    failure and it is silent -- gpu_free(Qs_d) yields the name "Qs_d)" while
    gpu_malloc_async(Qs_d, n, s) yields "Qs_d", so one buffer becomes two
    entries and a genuine MISMATCH is reported as an UNFREED plus a FREE-ONLY.
    Track the depth instead.
    """
    depth = 0
    out = []
    for ch in text[open_paren:]:
        if ch == "(":
            depth += 1
            if depth == 1:
                continue
        elif ch == ")":
            depth -= 1
            if depth == 0:
                break
        elif ch == "," and depth == 1:
            break
        out.append(ch)
    return "".join(out).strip()


def strip_comment(line):
    """Drop a trailing Fortran comment, honouring quotes."""
    out, q = [], None
    for ch in line:
        if q:
            if ch == q:
                q = None
        elif ch in "\"'":
            q = ch
        elif ch == "!":
            break
        out.append(ch)
    return "".join(out)


def scan(paths):
    events = defaultdict(list)  # buffer name -> [(kind, family, file, line, fn)]
    for path in paths:
        with open(path, errors="replace") as fh:
            for n, raw in enumerate(fh, 1):
                code = strip_comment(raw)
                m = CALL.search(code)
                if not m:
                    continue
                fn = m.group(1).lower()
                buf = first_arg(code, m.end() - 1)
                if not buf:
                    continue
                # normalise soap_turbo_hypers(i)%W_d -> %W_d so the same
                # component across loop indices is one buffer
                buf = re.sub(r"^.*%", "%", buf)
                kind, fam = FAMILY[fn]
                events[buf].append((kind, fam, path, n, fn))
    return events


def main(argv):
    quiet = "--quiet" in argv
    paths = [a for a in argv[1:] if not a.startswith("--")]
    if not paths:
        root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
        paths = sorted(glob.glob(os.path.join(root, "src", "*.f90")))
    if not paths:
        print("no sources found", file=sys.stderr)
        return 2

    events = scan(paths)
    mismatch, unfreed, freeonly = [], [], []

    for buf, evs in sorted(events.items()):
        allocs = [e for e in evs if e[0] == "alloc"]
        frees = [e for e in evs if e[0] == "free"]
        afam = {e[1] for e in allocs}
        ffam = {e[1] for e in frees}

        # Any free family the allocations do not cover is a mismatch. Asking
        # only whether the two SETS are disjoint misses the partial case: a
        # buffer freed with gpu_free_async at one site and gpu_free at another
        # intersects, and the plain free is still wrong.
        if allocs and frees and (ffam - afam):
            bad = [e for e in frees if e[1] not in afam]
            mismatch.append((buf, allocs, bad))
        elif allocs and not frees:
            unfreed.append((buf, allocs))
        elif frees and not allocs:
            freeonly.append((buf, frees))

    def show(label, rows, detail):
        if not rows:
            return
        print("\n%s (%d)" % (label, len(rows)))
        for row in rows:
            print("  %s" % row[0])
            for e in detail(row):
                print("      %-22s %s:%d" % (e[4], os.path.basename(e[2]), e[3]))

    show(
        "MISMATCH -- pool allocation freed with the plain deallocator, or vice versa",
        mismatch,
        lambda r: r[1] + r[2],
    )
    if not quiet:
        show("UNFREED -- allocated, never passed to a deallocator here", unfreed, lambda r: r[1])
        show("FREE-ONLY -- freed, never allocated here", freeonly, lambda r: r[1])

    print(
        "\n%d device buffers seen across %d files: %d mismatched, %d unfreed, %d free-only"
        % (len(events), len(paths), len(mismatch), len(unfreed), len(freeonly))
    )
    return 1 if mismatch else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
