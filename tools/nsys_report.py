#!/usr/bin/env python3
"""Turn an nsys SQLite export into the numbers that decide what to optimise.

    tools/nsys_report.py profiling/CO_predict/CO_predict.sqlite
    tools/nsys_report.py --top 30 --plot phases.png <sqlite>

`nsys stats` already prints per-kernel and per-API summaries, and
tools/profile_gpu.sh saves those as CSV. This exists for the four questions
those summaries cannot answer, all of which are about *time nothing is
happening* rather than time something is:

  1. What fraction of the run is the GPU actually busy? A kernel summary tells
     you the ranking of the kernels; it cannot tell you they collectively
     occupy 30% of the wall clock, which is the number that says whether to
     tune a kernel or to stop leaving the device idle.

  2. Within one phase -- SOAP, 2b, 3b, the batched pdf loop -- how much of that
     phase is device work? This is what the NVTX ranges are for. A phase that
     is 90% idle is a phase whose cost is host-side or launch-bound, and
     tuning its kernels is wasted effort.

  3. Do any two kernels ever run at the same time? TurboGAP creates one stream
     per OpenMP thread and hands it to the batched loops, so the question
     "are the streams doing anything" has a yes/no answer in the data: if peak
     concurrency is 1, every batch ran one after another regardless of how many
     streams were created.

  4. Where does the wall clock go between GPU operations? Summed gaps, split by
     the phase they fall in, distinguish "the host is computing" from "the host
     is launching thousands of tiny kernels" -- opposite fixes.

Reads the export with the standard library's sqlite3 and needs no packages, so
it works on a machine that only has the .sqlite file. --plot additionally needs
matplotlib (tools/setup_profiling_env.sh installs it).
"""

import argparse
import sqlite3
import sys
from collections import defaultdict

NS = 1e-9


# --------------------------------------------------------------------- schema
#
# nsys renames and adds tables between versions, and an export contains only the
# tables the trace actually collected. Every query here goes through these two
# helpers so a missing table degrades that one section instead of raising.
def tables(db):
    return {r[0] for r in db.execute(
        "SELECT name FROM sqlite_master WHERE type IN ('table','view')")}


def columns(db, table):
    return {r[1] for r in db.execute(f'PRAGMA table_info("{table}")')}


def string_map(db):
    if 'StringIds' not in tables(db):
        return {}
    return {i: s for i, s in db.execute('SELECT id, value FROM StringIds')}


# ------------------------------------------------------------------- intervals
def merge(intervals):
    """Union of [start, end) intervals, as a sorted, disjoint list.

    Overlap is the point: kernels on different streams overlap in time, so
    summing durations double-counts and can exceed the wall clock. Busy time is
    the measure of the union, never the sum.
    """
    if not intervals:
        return []
    intervals = sorted(intervals)
    out = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= out[-1][1]:
            if e > out[-1][1]:
                out[-1][1] = e
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def measure(intervals):
    return sum(e - s for s, e in intervals)


def overlap(intervals, lo, hi):
    """Measure of `intervals` (already merged) inside [lo, hi)."""
    total = 0
    for s, e in intervals:
        if e <= lo:
            continue
        if s >= hi:
            break
        total += min(e, hi) - max(s, lo)
    return total


def peak_concurrency(intervals):
    """Largest number of intervals overlapping at any instant, and how long the
    run spends at each level."""
    events = []
    for s, e in intervals:
        events.append((s, 1))
        events.append((e, -1))
    events.sort()
    depth = 0
    last = None
    at_level = defaultdict(float)
    peak = 0
    for t, d in events:
        if last is not None and t > last:
            at_level[depth] += t - last
        depth += d
        peak = max(peak, depth)
        last = t
    return peak, at_level


# ---------------------------------------------------------------- data loading
def gpu_ops(db):
    """Every device-side operation: (start, end, stream, kind, name, bytes)."""
    names = string_map(db)
    tabs = tables(db)
    ops = []

    if 'CUPTI_ACTIVITY_KIND_KERNEL' in tabs:
        cols = columns(db, 'CUPTI_ACTIVITY_KIND_KERNEL')
        # demangledName is the readable one; shortName exists on older exports.
        namecol = 'demangledName' if 'demangledName' in cols else 'shortName'
        stream = 'streamId' if 'streamId' in cols else 'NULL'
        for s, e, st, n in db.execute(
                f'SELECT start, end, {stream}, {namecol} '
                f'FROM CUPTI_ACTIVITY_KIND_KERNEL'):
            ops.append((s, e, st, 'kernel', names.get(n, f'<{n}>'), 0))

    if 'CUPTI_ACTIVITY_KIND_MEMCPY' in tabs:
        cols = columns(db, 'CUPTI_ACTIVITY_KIND_MEMCPY')
        stream = 'streamId' if 'streamId' in cols else 'NULL'
        size = 'bytes' if 'bytes' in cols else '0'
        kind = 'copyKind' if 'copyKind' in cols else '-1'
        for s, e, st, b, k in db.execute(
                f'SELECT start, end, {stream}, {size}, {kind} '
                f'FROM CUPTI_ACTIVITY_KIND_MEMCPY'):
            ops.append((s, e, st, f'memcpy{MEMCPY_KIND.get(k, k)}', '', b or 0))

    if 'CUPTI_ACTIVITY_KIND_MEMSET' in tabs:
        cols = columns(db, 'CUPTI_ACTIVITY_KIND_MEMSET')
        stream = 'streamId' if 'streamId' in cols else 'NULL'
        size = 'bytes' if 'bytes' in cols else '0'
        for s, e, st, b in db.execute(
                f'SELECT start, end, {stream}, {size} '
                f'FROM CUPTI_ACTIVITY_KIND_MEMSET'):
            ops.append((s, e, st, 'memset', '', b or 0))

    ops.sort()
    return ops


# CUPTI's copyKind enum. Named here rather than joined from ENUM_CUDA_MEMCPY_OPER
# because that table is not present in every export.
MEMCPY_KIND = {1: ' HtoD', 2: ' DtoH', 3: ' HtoH', 4: ' DtoD',
               8: ' HtoD(p)', 9: ' DtoH(p)', 10: ' DtoD(p)'}


def api_calls(db):
    """Host-side CUDA runtime calls: (start, end, name)."""
    tabs = tables(db)
    if 'CUPTI_ACTIVITY_KIND_RUNTIME' not in tabs:
        return []
    names = string_map(db)
    cols = columns(db, 'CUPTI_ACTIVITY_KIND_RUNTIME')
    namecol = 'nameId' if 'nameId' in cols else 'name'
    return [(s, e, names.get(n, f'<{n}>'))
            for s, e, n in db.execute(
                f'SELECT start, end, {namecol} FROM CUPTI_ACTIVITY_KIND_RUNTIME')]


def nvtx_ranges(db):
    """Our phase markers: (start, end, text, thread). Push/pop ranges only."""
    tabs = tables(db)
    if 'NVTX_EVENTS' not in tabs:
        return []
    names = string_map(db)
    cols = columns(db, 'NVTX_EVENTS')
    tid = 'globalTid' if 'globalTid' in cols else 'NULL'
    out = []
    for s, e, text, tid_v, textId in db.execute(
            f'SELECT start, end, text, {tid}, textId FROM NVTX_EVENTS'):
        if e is None or s is None:
            continue                      # a mark, not a range
        label = text if text else names.get(textId)
        if not label:
            continue
        out.append((s, e, label, tid_v))
    out.sort()
    return out


# -------------------------------------------------------------------- printing
def hdr(title):
    print()
    print(title)
    print('-' * len(title))


def secs(ns):
    return f'{ns * NS:10.4f}'


def pct(part, whole):
    return f'{100.0 * part / whole:5.1f}%' if whole else '    -'


def human_bytes(b):
    for unit in ('B', 'KB', 'MB', 'GB', 'TB'):
        if abs(b) < 1024 or unit == 'TB':
            return f'{b:7.1f} {unit}'
        b /= 1024.0


def report(path, top, plot):
    db = sqlite3.connect(f'file:{path}?mode=ro', uri=True)

    ops = gpu_ops(db)
    if not ops:
        print(f'{path}: no CUDA activity in this export.', file=sys.stderr)
        print('Was the run traced with -t cuda?', file=sys.stderr)
        return 1

    t0 = min(o[0] for o in ops)
    t1 = max(o[1] for o in ops)
    span = t1 - t0

    busy_iv = merge([(o[0], o[1]) for o in ops])
    busy = measure(busy_iv)
    kernel_iv = merge([(o[0], o[1]) for o in ops if o[3] == 'kernel'])

    print(f'nsys report: {path}')
    print(f'{len(ops)} GPU operations over {span * NS:.3f} s '
          f'(first launch to last completion)')

    # ---------------------------------------------------------------- headline
    hdr('Device occupancy')
    print(f'  wall clock (GPU span)   {secs(span)} s')
    print(f'  GPU busy (any op)       {secs(busy)} s   {pct(busy, span)}')
    print(f'  GPU busy (kernels only) {secs(measure(kernel_iv))} s   '
          f'{pct(measure(kernel_iv), span)}')
    print(f'  GPU IDLE                {secs(span - busy)} s   {pct(span - busy, span)}')
    print()
    print('  Idle is the budget for overlap. A run that is mostly idle is not')
    print('  short of FLOPs -- it is waiting on the host, on launches, or on')
    print('  synchronisation, and none of those is fixed by a faster kernel.')

    # ------------------------------------------------------------ concurrency
    streams = sorted({o[2] for o in ops if o[2] is not None})
    peak, at_level = peak_concurrency([(o[0], o[1]) for o in ops])
    hdr('Stream concurrency')
    print(f'  distinct streams used   {len(streams)}   {streams[:12]}'
          f'{" ..." if len(streams) > 12 else ""}')
    print(f'  peak concurrent ops     {peak}')
    for lvl in sorted(at_level):
        if lvl == 0 or at_level[lvl] <= 0:
            continue
        print(f'    {lvl} op(s) at once      {secs(at_level[lvl])} s   '
              f'{pct(at_level[lvl], span)}')
    if peak <= 1:
        print()
        print('  Peak concurrency 1: nothing on the device ever overlapped.')
        print('  Extra streams, however many were created, bought nothing on')
        print('  this run.')

    # ---------------------------------------------------------------- kernels
    agg = defaultdict(lambda: [0, 0])         # name -> [count, total ns]
    for s, e, _st, kind, name, _b in ops:
        if kind == 'kernel':
            a = agg[name]
            a[0] += 1
            a[1] += e - s
    hdr(f'Top {top} kernels (by summed duration; overlapping streams double-count)')
    print(f'  {"total s":>10} {"%span":>6} {"count":>8} {"mean us":>10}  kernel')
    for name, (n, tot) in sorted(agg.items(), key=lambda kv: -kv[1][1])[:top]:
        print(f'  {secs(tot)} {pct(tot, span)} {n:8d} {tot / n / 1000.0:10.1f}  '
              f'{name[:78]}')

    # -------------------------------------------------------------- transfers
    hdr('Host <-> device transfers')
    xfer = defaultdict(lambda: [0, 0, 0])     # kind -> [count, ns, bytes]
    for s, e, _st, kind, _n, b in ops:
        if kind.startswith('memcpy') or kind == 'memset':
            a = xfer[kind]
            a[0] += 1
            a[1] += e - s
            a[2] += b
    if xfer:
        print(f'  {"total s":>10} {"%span":>6} {"count":>8} {"volume":>12} '
              f'{"GB/s":>8}  kind')
        for kind, (n, ns, b) in sorted(xfer.items(), key=lambda kv: -kv[1][1]):
            bw = (b / (ns * NS) / 1e9) if ns else 0.0
            print(f'  {secs(ns)} {pct(ns, span)} {n:8d} {human_bytes(b)} '
                  f'{bw:8.1f}  {kind}')
    else:
        print('  none recorded')

    # ------------------------------------------------------------------ phases
    ranges = nvtx_ranges(db)
    outer, inner = [], []
    if ranges:
        # Only outermost ranges per thread, so a phase and the sub-phases inside
        # it are not both charged the same device time.
        by_thread = defaultdict(list)
        for s, e, label, tid in ranges:
            by_thread[tid].append((s, e, label))
        for tid, rs in by_thread.items():
            rs.sort()
            end_of_open = 0
            for s, e, label in rs:
                (outer if s >= end_of_open else inner).append((s, e, label))
                end_of_open = max(end_of_open, e) if s >= end_of_open else end_of_open

        def phase_table(title, rows, note=None):
            hdr(title)
            if note:
                print(f'  {note}')
            agg = defaultdict(lambda: [0, 0, 0])   # label -> [count, wall, gpu]
            for s, e, label in rows:
                a = agg[label]
                a[0] += 1
                a[1] += e - s
                a[2] += overlap(busy_iv, s, e)
            print(f'  {"wall s":>10} {"%span":>6} {"gpu s":>10} {"gpu%":>6} '
                  f'{"count":>7}  phase')
            for label, (n, wall, g) in sorted(agg.items(), key=lambda kv: -kv[1][1]):
                print(f'  {secs(wall)} {pct(wall, span)} {secs(g)} '
                      f'{pct(g, wall)} {n:7d}  {label}')
            return agg

        outer_agg = phase_table(
            'Phases (NVTX, outermost only)',
            outer,
            'gpu% is the fraction of that phase during which the device was busy.')
        if inner:
            phase_table('Phases (nested)', inner)
    else:
        outer_agg = {}
        hdr('Phases (NVTX)')
        print('  No NVTX ranges in this export.')
        print('  Either the binary was not built with PROFILE=1, or nsys was not')
        print('  run with -t nvtx. tools/profile_gpu.sh checks the first.')

    # ------------------------------------------------------------- host stalls
    calls = api_calls(db)
    if calls:
        hdr(f'Top {top} CUDA API calls (host side, by summed duration)')
        agg = defaultdict(lambda: [0, 0])
        for s, e, name in calls:
            a = agg[name]
            a[0] += 1
            a[1] += e - s
        print(f'  {"total s":>10} {"%span":>6} {"count":>8} {"mean us":>10}  call')
        for name, (n, tot) in sorted(agg.items(), key=lambda kv: -kv[1][1])[:top]:
            print(f'  {secs(tot)} {pct(tot, span)} {n:8d} {tot / n / 1000.0:10.1f}  '
                  f'{name}')
        print()
        # Allocator churn is the specific thing worth naming. cudaMalloc and
        # cudaFree serialise against the whole context: they are not stream
        # operations, so a per-batch malloc/free pair inside a loop the streams
        # were meant to overlap is a barrier no amount of streaming removes.
        alloc = sum(t for nm, (_c, t) in agg.items()
                    if 'Malloc' in nm or 'Free' in nm)
        if alloc:
            print(f'  of which allocation ({"cudaMalloc/cudaFree".ljust(20)}) '
                  f'{secs(alloc)} s  {pct(alloc, span)}')
            print('  cudaMalloc and cudaFree are context-wide barriers, not stream')
            print('  operations. Allocation inside a loop serialises that loop no')
            print('  matter how many streams it uses.')

    # -------------------------------------------------------------------- gaps
    hdr('Gaps between GPU operations')
    gaps = []
    for (_s0, e0), (s1, _e1) in zip(busy_iv, busy_iv[1:]):
        gaps.append((e0, s1 - e0))
    gaps.sort(key=lambda g: -g[1])
    total_gap = sum(g[1] for g in gaps)
    print(f'  {len(gaps)} gaps totalling {total_gap * NS:.3f} s '
          f'({pct(total_gap, span).strip()} of the span)')
    if gaps:
        print(f'  largest {min(10, len(gaps))}:')
        labels = [(s, e, l) for s, e, l in (outer if ranges else [])]
        for at, dur in gaps[:10]:
            where = next((l for s, e, l in labels if s <= at < e), '-')
            print(f'    {dur * NS:9.4f} s  at t+{(at - t0) * NS:8.4f} s  in {where}')

    if plot:
        make_plot(plot, span, busy, outer_agg)

    db.close()
    return 0


def make_plot(path, span, busy, phases):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('\n--plot needs matplotlib: run tools/setup_profiling_env.sh',
              file=sys.stderr)
        return
    if not phases:
        print('\n--plot needs NVTX phases; nothing to draw.', file=sys.stderr)
        return
    rows = sorted(phases.items(), key=lambda kv: kv[1][1])[-15:]
    labels = [r[0] for r in rows]
    wall = [r[1][1] * NS for r in rows]
    gpu = [r[1][2] * NS for r in rows]
    fig, ax = plt.subplots(figsize=(9, 0.4 * len(rows) + 2))
    ax.barh(labels, wall, label='wall', color='#c9d5e3')
    ax.barh(labels, gpu, label='GPU busy', color='#2b6cb0')
    ax.set_xlabel('seconds')
    ax.set_title(f'TurboGAP phases  (span {span * NS:.2f} s, '
                 f'GPU busy {100 * busy / span:.0f}%)')
    ax.legend(loc='lower right')
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    print(f'\nwrote {path}')


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('sqlite', help='the .sqlite produced by `nsys export`')
    p.add_argument('--top', type=int, default=15, help='rows per table (default 15)')
    p.add_argument('--plot', metavar='PNG', help='write a phase bar chart')
    a = p.parse_args()
    return report(a.sqlite, a.top, a.plot)


if __name__ == '__main__':
    sys.exit(main())
