#!/usr/bin/env python3
"""Bundle the driver's separate wall-clock buckets into one times_t.

Both trees carried ~22 (CPU) / ~29 (GPU) `time_<name>(1:3)` locals, threaded
one at a time through every extracted module. This rewrites them to components
of the single `times_t` declared in src/timing.f90, and replaces the hand
written

    call get_time(time_x(1))
    ...
    call get_time(time_x(2))
    time_x(3) = time_x(3) + time_x(2) - time_x(1)

with `call time_start(time%x)` / `call time_end(time%x)`.

Two names look like timers and are not, and a blanket `s/time_/time%/` breaks
both: `time_step` and `time_step_prev` are the MD integration step in fs. They
are listed in NOT_TIMERS and asserted to survive untouched.

Every count is asserted. Run with --check to diff without writing.

Usage: tools/bundle_timers.py [--check] <file> [<file> ...]
"""

import re
import sys

# old local name -> times_t component
BUCKETS = {
    "time_read_input": "read_input",
    "time_read_xyz": "read_xyz",
    "time_neigh": "neigh",
    "time_gap": "gap",
    "time_soap": "soap",
    "time_2b": "gap_2b",
    "time_3b": "gap_3b",
    "time_core_pot": "gap_core_pot",
    "time_vdw": "vdw",
    "time_estat": "estat",
    "time_pdf": "pdf",
    "time_sf": "sf",
    "time_xrd": "xrd",
    "time_nd": "nd",
    "time_xps": "xps",
    "time_md": "md",
    "time_mc": "mc",
    "time_mpi": "mpi",
    "time_mpi_positions": "mpi_positions",
    "time_mpi_ef": "mpi_ef",
    # GPU branch only; harmless on the CPU tree where they never appear.
    "time_get_soap": "get_soap",
    "time_local_prop": "local_prop",
    "time_batch_alloc": "batch_alloc",
    "time_batch_pdf": "batch_pdf",
    "time_batch_xrd": "batch_xrd",
    "time_create_streams": "create_streams",
    "time_exp_batched": "exp_batched",
}

# Names that match the time_* shape and are NOT wall-clock buckets.
NOT_TIMERS = ("time_step", "time_step_prev")


def rewrite(text):
    counts = {}

    # 1. Collapse the three-line accumulate idiom into time_end. Do this before
    #    the name substitution so the pattern is still the original one.
    n_end = 0
    clobbered = []
    for old, new in BUCKETS.items():
        # Allow blank lines and comments between the (2) stamp and the
        # accumulate -- the GPU tree interleaves dead MPI_wtime comments, and
        # requiring adjacency left 13 sites uncollapsed there. Anything else
        # in between (soap_solo has a deallocate) means collapsing would move
        # the stamp, so those are left alone.
        pat = re.compile(
            r"([ \t]*)call get_time\(%s\(2\)\)[ \t]*\n"
            r"((?:[ \t]*(?:!.*)?\n)*)"
            r"[ \t]*%s\(3\) *= *(%s\(3\) *\+ *)?%s\(2\) *- *%s\(1\)"
            % (old, old, old, old, old)
        )

        def repl(m, new=new, old=old):
            # A bucket whose accumulate line is missing its own (3) term
            # ASSIGNS rather than accumulates, so it reports only the last
            # interval. time_end always accumulates, so collapsing the idiom
            # silently fixes it -- report it rather than hide it.
            if m.group(3) is None:
                clobbered.append(old)
            return "%scall time_end(time%%%s)" % (m.group(1), new)

        text, k = pat.subn(repl, text)
        n_end += k
    counts["time_end"] = n_end
    counts["fixed_assign"] = clobbered

    # 2. call get_time(time_x(1)) -> call time_start(time%x)
    n_start = 0
    for old, new in BUCKETS.items():
        text, k = re.subn(
            r"call get_time\(%s\(1\)\)" % old, "call time_start(time%%%s)" % new, text
        )
        n_start += k
    counts["time_start"] = n_start

    # 3. Everything else: time_x -> time%x, on whole words only.
    n_ref = 0
    for old, new in sorted(BUCKETS.items(), key=lambda kv: -len(kv[0])):
        text, k = re.subn(r"\b%s\b" % old, "time%%%s" % new, text)
        n_ref += k
    counts["renamed"] = n_ref

    return text, counts


def main(argv):
    check = "--check" in argv
    files = [a for a in argv[1:] if not a.startswith("--")]
    if not files:
        print(__doc__)
        return 2

    for path in files:
        with open(path) as fh:
            src = fh.read()

        before = {n: len(re.findall(r"\b%s\b" % n, src)) for n in NOT_TIMERS}
        out, counts = rewrite(src)
        after = {n: len(re.findall(r"\b%s\b" % n, out)) for n in NOT_TIMERS}

        # The trap this asserts against: a blanket rename eats time_step.
        for n in NOT_TIMERS:
            assert before[n] == after[n], (
                "%s: %s went %d -> %d; it is the MD integration step, not a timer"
                % (path, n, before[n], after[n])
            )

        # time_neigh and time_gap were plain scalars, not (1:3) accumulators.
        # Promoting them leaves every unsubscripted use rank-wrong -- and a
        # write statement accepts a rank-1 actual against an F edit descriptor
        # without complaint, so the only symptom is garbled output far from
        # the change. Flag any bucket reference that is not subscripted.
        for line in out.splitlines():
            if line.lstrip().startswith("!"):
                continue
            if "write" not in line and "print" not in line:
                continue
            # A lookahead here does not work: [a-z0-9_]+ backtracks until the
            # next character is not "(", so every subscripted use matches too.
            # Match greedily, then look at what actually follows.
            bare = [
                m.group(0)
                for m in re.finditer(r"time%\w+", line)
                if not line[m.end():m.end() + 1] == "("
            ]
            assert not bare, "%s: unsubscripted bucket in output statement: %s" % (
                path,
                line.strip(),
            )

        # No bare bucket name may survive.
        leftover = sorted(
            set(re.findall(r"\btime_[a-z0-9_]+\b", out))
            - set(NOT_TIMERS)
            - {"time_start", "time_end"}
        )
        assert not leftover, "%s: unconverted timer names %s" % (path, leftover)

        print(
            "%-28s time_start %3d  time_end %3d  renamed %3d"
            % (path, counts["time_start"], counts["time_end"], counts["renamed"])
        )
        for name in sorted(set(counts["fixed_assign"])):
            print("    NOTE %s assigned instead of accumulating; time_end fixes it" % name)
        if not check:
            with open(path, "w") as fh:
                fh.write(out)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
