#!/usr/bin/env python3
"""Self-test for tools/compare_xyz.py.

Checks that the comparison actually discriminates: identical trajectories
must pass, and a trajectory perturbed in energy, force or virial must fail.
Without this a silently broken comparator would make the whole GPU
regression suite pass unconditionally.

Runs anywhere; needs no GPU and no TurboGAP build.
"""
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
COMPARE = os.path.join(HERE, "..", "..", "tools", "compare_xyz.py")

HEADER = (
    'Properties=species:S:1:pos:R:3:forces:R:3:local_energy:R:1:hirshfeld_v:R:1 '
    'Lattice="10.0 0.0 0.0 0.0 10.0 0.0 0.0 0.0 10.0" '
    'energy={energy} energy_soap=-2.0 energy_2b=0.5 energy_3b=0.25 '
    'virial="{v} 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0"'
)


def write_xyz(path, energy=-1.5, fx=0.125, virial=1.0, local=-0.75):
    with open(path, "w") as f:
        f.write("2\n")
        f.write(HEADER.format(energy="%.8f" % energy, v="%.8f" % virial) + "\n")
        f.write("C 0.0 0.0 0.0 %.8f 0.0 0.0 %.8f 1.0\n" % (fx, local))
        f.write("O 1.1 0.0 0.0 %.8f 0.0 0.0 %.8f 1.0\n" % (-fx, local))


def compare(a, b):
    r = subprocess.run([sys.executable, COMPARE, a, b],
                       capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def main():
    failures = []
    with tempfile.TemporaryDirectory() as d:
        ref = os.path.join(d, "ref.xyz")
        write_xyz(ref)

        cases = [
            ("identical",        dict(),                 0),
            ("energy perturbed", dict(energy=-1.4),      1),
            ("forces perturbed", dict(fx=0.5),           1),
            ("virial perturbed", dict(virial=2.0),       1),
        ]
        for name, kwargs, want in cases:
            other = os.path.join(d, "other.xyz")
            write_xyz(other, **kwargs)
            rc, out = compare(ref, other)
            got = 0 if rc == 0 else 1
            ok = got == want
            print("%-18s expected %s, got %s  %s"
                  % (name,
                     "PASS" if want == 0 else "FAIL",
                     "PASS" if got == 0 else "FAIL",
                     "ok" if ok else "<-- WRONG"))
            if not ok:
                failures.append((name, out))

    if failures:
        print("\ncompare_xyz.py self-test FAILED")
        for name, out in failures:
            print("--- %s ---" % name)
            print(out)
        return 1
    print("\ncompare_xyz.py self-test passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
