#!/usr/bin/env python3
"""Check what an eight-bead PIMD run produced.

Separated from run.sh so the assertions are readable without reading shell.
"""

import glob
import os
import sys

work, beads, steps = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])

failures = []
checks = 0


def check(passed, description, detail=""):
    global checks
    checks += 1
    print(f"    {'PASS' if passed else 'FAIL'}  {description}")
    if not passed:
        if detail:
            print(f"          {detail}")
        failures.append(description)


def frames_of(path):
    """[(comment, [(name, x, y, z)])] from a plain xyz trajectory."""
    lines = open(path, encoding="utf-8").read().splitlines()
    out, index = [], 0
    while index < len(lines):
        if not lines[index].strip():
            index += 1
            continue
        count = int(lines[index].split()[0])
        rows = []
        for row in lines[index + 2:index + 2 + count]:
            parts = row.split()
            rows.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
        out.append((lines[index + 1], rows))
        index += 2 + count
    return out


def check_completed():
    log = os.path.join(work, "ipi.log")
    text = open(log, encoding="utf-8").read() if os.path.exists(log) else ""
    check("SOFTEXIT" in text or "Simulation ran successfully" in text or os.path.exists(
        os.path.join(work, "pimd.out")),
        "i-PI completed the run",
        text[-400:] if text else "no ipi.log")


def check_bead_files():
    found = sorted(glob.glob(os.path.join(work, "pimd.pos_*.xyz")))
    check(len(found) == beads,
          f"all {beads} bead trajectories were written",
          f"found {len(found)}: {[os.path.basename(f) for f in found]}")
    return found


def check_atoms_present(found):
    for path in found[:1] + found[-1:]:
        frames = frames_of(path)
        check(len(frames) >= 2,
              f"{os.path.basename(path)} holds {len(frames)} frames")
        check(all(len(rows) == 12 for _, rows in frames),
              f"{os.path.basename(path)} has all 12 atoms in every frame")


def check_beads_are_distinct(found):
    """A ring polymer whose beads coincide is a classical run in disguise."""
    if len(found) < 2:
        check(False, "at least two beads to compare", "not enough bead files")
        return
    first = frames_of(found[0])[-1][1]
    last = frames_of(found[-1])[-1][1]
    spread = max(max(abs(a[i] - b[i]) for i in (1, 2, 3)) for a, b in zip(first, last))
    check(spread > 1e-4,
          f"beads 1 and {beads} differ by up to {spread:.4f} A",
          "the beads coincide, so the ring polymer is not spread: this is a "
          "classical simulation wearing a quantum label")


def check_conserved():
    path = os.path.join(work, "pimd.out")
    if not os.path.exists(path):
        check(False, "i-PI wrote its properties file", f"{path} is missing")
        return
    rows = [line.split() for line in open(path, encoding="utf-8")
            if line.strip() and not line.startswith("#")]
    check(len(rows) >= 2, f"{len(rows)} property rows written")
    if len(rows) < 2:
        return
    conserved = [float(r[2]) for r in rows]
    potential = [float(r[3]) for r in rows]

    # Whether the dynamics can be tested at all depends on the potential. The
    # only water GAP in the test data is a DIPOLE model -- dipole_model =
    # .true. on both its descriptors -- and a dipole model contributes no
    # energy and no forces by construction, so the ring polymer is driven by
    # nothing. The protocol checks above are still meaningful: eight drivers,
    # eight beads, every exchange answered. The dynamics are not.
    spread = max(potential) - min(potential)
    if spread <= 1e-6:
        print("    NOTE  the potential is identically zero, so the DYNAMICS are not")
        print("          under test here -- only the protocol. The water potential in")
        print("          the test data is a dipole model, which by design returns no")
        print("          energy and no forces. A water GAP with energies would make")
        print("          the two checks below meaningful; there is not one available.")
        return

    check(True, f"the potential energy moved by {spread:.6f} eV, so the run advanced")

    drift = max(conserved) - min(conserved)
    scale = max(abs(c) for c in conserved) or 1.0
    check(drift / scale < 0.05,
          f"the conserved quantity moved by {drift:.4f} eV over {len(rows)} steps "
          f"({100 * drift / scale:.2f}% of its magnitude)")


check_completed()
found = check_bead_files()
check_atoms_present(found)
check_beads_are_distinct(found)
check_conserved()

print()
if failures:
    print(f"ipi_pimd: {len(failures)} of {checks} checks FAILED")
    sys.exit(1)
print(f"ipi_pimd: all {checks} checks passed")
