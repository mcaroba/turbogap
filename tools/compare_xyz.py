#!/usr/bin/env python3
"""Compare energies/forces/virial between two TurboGAP extended-XYZ trajectories.

Pure standard library on purpose: this runs on cluster nodes and in CI where
numpy is often not installed, and a comparison tool that cannot run is worse
than no comparison tool.

Usage: compare_xyz.py <reference.xyz> <test.xyz> [--ftol F] [--etol E]
Exit status is 0 if every checked quantity is within tolerance, 1 otherwise.
"""
import re
import sys

DEF_FTOL = 1e-4   # absolute, on forces
DEF_ETOL = 1e-6   # relative, on energies and virial


def parse_frames(path):
    frames = []
    with open(path) as f:
        while True:
            line = f.readline()
            if not line:
                break
            if not line.strip():
                continue
            nat = int(line.split()[0])
            comment = f.readline()
            props = dict(re.findall(r'(\w+)=("[^"]*"|\S+)', comment))
            rows = [f.readline().split() for _ in range(nat)]
            frames.append((props, rows))
    return frames


def getf(props, key):
    v = props.get(key)
    return None if v is None else float(v.strip('"'))


def getvec(props, key):
    v = props.get(key)
    return None if v is None else [float(x) for x in v.strip('"').split()]


def rms(xs):
    return (sum(x * x for x in xs) / len(xs)) ** 0.5 if xs else 0.0


def main(argv):
    if len(argv) < 3:
        print(__doc__)
        return 2
    a_path, b_path = argv[1], argv[2]
    ftol, etol = DEF_FTOL, DEF_ETOL
    for i, arg in enumerate(argv):
        if arg == "--ftol":
            ftol = float(argv[i + 1])
        elif arg == "--etol":
            etol = float(argv[i + 1])

    A, B = parse_frames(a_path), parse_frames(b_path)
    print("frames: %s=%d  %s=%d" % (a_path, len(A), b_path, len(B)))
    if not A or not B:
        print("RESULT: FAIL (no frames parsed)")
        return 1

    failed = False
    for i in range(min(len(A), len(B))):
        pa, ca = A[i]
        pb, cb = B[i]
        print("\n--- frame %d ---" % i)

        for k in ("energy", "energy_soap", "energy_2b", "energy_3b",
                  "energy_core_pot", "energy_vdw"):
            va, vb = getf(pa, k), getf(pb, k)
            if va is None or vb is None:
                continue
            d = abs(va - vb)
            rel = d / max(abs(va), 1e-30)
            ok = rel < etol or d < 1e-6
            print("  %s %-16s %20.8f %20.8f  absdiff=%.3e rel=%.2e"
                  % ("OK " if ok else "DIFF", k, va, vb, d, rel))
            failed |= not ok

        va, vb = getvec(pa, "virial"), getvec(pb, "virial")
        if va and vb and len(va) == len(vb):
            d = max(abs(x - y) for x, y in zip(va, vb))
            rel = d / max(max(abs(x) for x in va), 1e-30)
            ok = rel < etol
            print("  %s %-16s maxabsdiff=%.6e rel=%.2e"
                  % ("OK " if ok else "DIFF", "virial", d, rel))
            if not ok:
                print("       ref=%s" % va)
                print("       new=%s" % vb)
            failed |= not ok

        # columns: species, x, y, z, fx, fy, fz, ...
        try:
            fa = [[float(v) for v in r[4:7]] for r in ca]
            fb = [[float(v) for v in r[4:7]] for r in cb]
        except (IndexError, ValueError):
            fa = fb = []
        if fa and len(fa) == len(fb):
            diffs = [abs(x - y) for ra, rb in zip(fa, fb)
                     for x, y in zip(ra, rb)]
            md = max(diffs)
            fmax = max(abs(x) for r in fa for x in r)
            ok = md < ftol
            print("  %s %-16s maxabsdiff=%.6e (|F|max=%.4f) rms=%.3e"
                  % ("OK " if ok else "DIFF", "forces", md, fmax, rms(diffs)))
            failed |= not ok

        try:
            la = [float(r[7]) for r in ca]
            lb = [float(r[7]) for r in cb]
            md = max(abs(x - y) for x, y in zip(la, lb))
            print("  %s %-16s maxabsdiff=%.6e"
                  % ("OK " if md < 1e-6 else "DIFF", "local_energy", md))
            failed |= md >= 1e-6
        except (IndexError, ValueError):
            pass

    print("\nRESULT:", "FAIL" if failed else "PASS")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
