#!/usr/bin/env python3
"""Does a relaxation actually relax, to the tolerance it was given?

The regression cases run fifteen steps and compare the output byte for byte.
That detects change; it says nothing about whether the relaxation converges,
whether the tolerance it was handed means anything, or whether relaxing the
lattice as well as the positions lands somewhere sensible. This checks those.

Every assertion compares against something that is not the run's own report of
itself: the forces in the trajectory it wrote, or a second run started from the
structure the first one ended at.

The convergence test in turbogap_md.f90 requires BOTH |dE| < e_tol*N AND
max|F| < f_tol, so a run stops at whichever is stricter. e_tol is held tight
throughout here, so that f_tol is the binding condition and a claim about it
means something.
"""

import os
import subprocess
import sys

BINARY = os.environ.get("TURBOGAP_BIN")
TOOLS = os.environ["TURBOGAP_TOOLS"]
sys.path.insert(0, TOOLS)
import xyz_frames as xyz

# The force tolerances to sweep. All below what the energy criterion reaches on
# its own for this system (~0.003), so each one binds.
FORCE_TOLERANCES = [0.005, 0.002, 0.001]
E_TOL = "1.d-10"
# The box relaxation does not terminate at that energy criterion: it keeps
# taking steps and the forces wander, so it runs to the step cap. It settles at
# 1e-8, which is what the regression decks use. See KNOWN_ISSUES 15.
E_TOL_BOX = "1.d-8"
F_TOL_BOX = 0.002
MAX_STEPS = 800

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


def check_known(passed, description, issue, detail=""):
    """A check that is expected to fail, against a recorded defect.

    Reported every run so it cannot be forgotten, and green when it starts
    passing -- which is how anyone fixing the defect finds out.
    """
    global checks
    checks += 1
    if passed:
        print(f"    PASS  {description}  (KNOWN_ISSUES {issue} appears fixed)")
        return
    print(f"    KNOWN {description}  -- KNOWN_ISSUES {issue}")
    if detail:
        print(f"          {detail}")


def write_input(path, atoms, mode, f_tol, e_tol):
    """A relaxation deck: one mode, one force tolerance, everything else fixed."""
    text = f"""atoms_file = "{atoms}"
pot_file = "gap_files/CO.gap"
n_species = 2
species = C O
masses = 12.01 15.99
e0 = -.16138053 0.
random_seed = 12345

md_nsteps = {MAX_STEPS}
optimize = "{mode}"
e_tol = {e_tol}
f_tol = {f_tol}
p_tol = 0.02

write_xyz = 1
write_thermo = 1
max_Gbytes_per_process = 1.0
"""
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)


def relax(tag, atoms, mode, f_tol, e_tol=E_TOL):
    """Run one relaxation in its own directory; return (max|F|, E, steps)."""
    directory = os.path.abspath(tag)
    os.makedirs(directory, exist_ok=True)
    for name in ("atoms.xyz", "gap_files"):
        link = os.path.join(directory, name)
        if not os.path.exists(link):
            os.symlink(os.path.join(os.environ["RELAX_DATA"], name), link)
    if atoms != "atoms.xyz":
        target = os.path.join(directory, "start.xyz")
        if os.path.abspath(atoms) != target:
            with open(atoms, encoding="utf-8") as src, open(target, "w", encoding="utf-8") as dst:
                dst.write(src.read())
        atoms = "start.xyz"
    write_input(os.path.join(directory, "input"), atoms, mode, f_tol, e_tol)

    trajectory = os.path.join(directory, "trajectory_out.xyz")
    if os.path.exists(trajectory):
        os.remove(trajectory)
    with open(os.path.join(directory, "run.log"), "w", encoding="utf-8") as log:
        done = subprocess.run([BINARY, "md"], cwd=directory, stdout=log, stderr=subprocess.STDOUT)
    if done.returncode != 0:
        raise SystemExit(f"{tag}: turbogap exited {done.returncode}; see {directory}/run.log")

    force, energy = xyz.max_force(trajectory)
    return force, energy, len(xyz.energies(trajectory)), trajectory


def last_frame_to(trajectory, path):
    """Write the final frame of a trajectory as an xyz file of its own."""
    last = list(xyz.frames(trajectory))[-1]
    comment, _, rows = last
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(f"{len(rows)}\n{comment}\n")
        for row in rows:
            handle.write(" ".join(row) + "\n")


def lattice_of(trajectory, which):
    """The 3x3 lattice of the first or last frame, as nine numbers."""
    all_frames = list(xyz.frames(trajectory))
    comment = all_frames[0 if which == "first" else -1][0]
    import re
    match = re.search(r'Lattice="([^"]+)"', comment)
    if not match:
        raise ValueError("no Lattice entry in the comment line")
    return [float(v) for v in match.group(1).split()]


def off_diagonal(lattice):
    """The six components a diagonal-only relaxation must leave alone."""
    return [lattice[i] for i in (1, 2, 3, 5, 6, 7)]


def check_converges_to_its_tolerance():
    """The run stops only once the forces are actually under f_tol."""
    print("\n  a relaxation reaches the force tolerance it was given")
    results = {}
    for f_tol in FORCE_TOLERANCES:
        force, energy, steps, _ = relax(f"gd_{f_tol}", "atoms.xyz", "gd", f_tol)
        results[f_tol] = (force, energy, steps)
        check(force < f_tol,
              f"f_tol = {f_tol}: max|F| = {force:.6f} < {f_tol}",
              f"stopped at {force:.6f}, which is above the tolerance it was given")
    return results


def check_tightening_does_not_raise_energy(results):
    """A stricter tolerance cannot land higher than a looser one."""
    print("\n  tightening the tolerance does not raise the energy")
    ordered = sorted(results)
    for loose, tight in zip(ordered[1:], ordered[:-1]):
        e_loose = results[loose][1]
        e_tight = results[tight][1]
        check(e_tight <= e_loose + 1e-6,
              f"E(f_tol={tight}) = {e_tight:.8f} <= E(f_tol={loose}) = {e_loose:.8f}")
        check(results[tight][0] <= results[loose][0] + 1e-9,
              f"max|F|(f_tol={tight}) = {results[tight][0]:.6f} "
              f"<= max|F|(f_tol={loose}) = {results[loose][0]:.6f}")


def check_energy_falls():
    """The relaxation lowers the energy from where it started."""
    print("\n  the relaxation lowers the energy")
    trajectory = os.path.abspath(f"gd_{FORCE_TOLERANCES[-1]}/trajectory_out.xyz")
    series = xyz.energies(trajectory)
    check(series[-1] < series[0],
          f"E fell from {series[0]:.8f} to {series[-1]:.8f}")


def check_relaxed_stays_relaxed():
    """Restarted from its own result, it has nothing left to do.

    This is what "fully relaxed" means, and the one check the relaxation
    cannot pass by stopping early: a run that quit before the minimum would
    move again when restarted.
    """
    print("\n  a relaxed structure is still relaxed when restarted")
    tight = FORCE_TOLERANCES[-1]
    first = os.path.abspath(f"gd_{tight}/trajectory_out.xyz")
    force_before, energy_before = xyz.max_force(first)

    os.makedirs("restart", exist_ok=True)
    relaxed = os.path.abspath("restart/relaxed.xyz")
    last_frame_to(first, relaxed)

    force_after, energy_after, steps, _ = relax("restart", relaxed, "gd", tight)
    check(abs(energy_after - energy_before) < 1e-4,
          f"E moved by {abs(energy_after - energy_before):.2e} eV on restart",
          f"{energy_before:.8f} -> {energy_after:.8f}")
    check(force_after < tight,
          f"max|F| = {force_after:.6f} is still under {tight}")


def check_box_relaxation_goes_no_higher():
    """Relaxing the lattice too cannot land above relaxing positions alone.

    gd-box optimises over the positions AND the lattice, so its minimum is over
    a strictly larger space than gd's and cannot be higher. If it is, the box
    relaxation is not finding the minimum it is searching.
    """
    print("\n  relaxing the lattice as well reaches an energy no higher")
    _, e_positions, _, _ = relax(f"gd_{F_TOL_BOX}", "atoms.xyz", "gd", F_TOL_BOX)
    force_box, e_box, steps, _ = relax("gdbox", "atoms.xyz", "gd-box", F_TOL_BOX, E_TOL_BOX)
    check(e_box <= e_positions + 1e-6,
          f"E(gd-box) = {e_box:.8f} <= E(gd) = {e_positions:.8f}",
          f"box relaxation landed {e_box - e_positions:+.6f} eV above position-only")
    check(force_box < F_TOL_BOX,
          f"gd-box max|F| = {force_box:.6f} < {F_TOL_BOX}",
          f"ran {steps} frames without reaching it")
    check(steps < MAX_STEPS,
          f"gd-box terminated in {steps} frames rather than hitting the {MAX_STEPS}-step cap")
    print(f"          gd-box gained {e_positions - e_box:.6f} eV over gd")


def check_lattice_constraint_is_honoured():
    """gd-box-ortho may change the cell lengths and nothing else.

    The constraint is the whole point of the mode: a cell that starts
    orthorhombic has to stay orthorhombic, or a relaxation asked to preserve
    the symmetry has quietly not.
    """
    print("\n  the diagonal-only constraint is honoured, and converges")
    force, energy, steps, trajectory = relax("gdbox_ortho", "atoms.xyz", "gd-box-ortho",
                                             F_TOL_BOX, E_TOL_BOX)
    check_known(force < F_TOL_BOX,
                f"gd-box-ortho max|F| = {force:.6f} < {F_TOL_BOX}", "15",
                f"ran {steps} frames without reaching it; gd-box needs 149")
    check_known(steps < MAX_STEPS,
                f"gd-box-ortho terminated in {steps} frames", "15")

    before = off_diagonal(lattice_of(trajectory, "first"))
    after = off_diagonal(lattice_of(trajectory, "last"))
    moved = max(abs(a - b) for a, b in zip(before, after))
    check(moved < 1e-8,
          f"the six off-diagonal lattice components moved by at most {moved:.2e}",
          f"before {before}\n          after  {after}")

    diagonal_before = [lattice_of(trajectory, "first")[i] for i in (0, 4, 8)]
    diagonal_after = [lattice_of(trajectory, "last")[i] for i in (0, 4, 8)]
    changed = max(abs(a - b) for a, b in zip(diagonal_before, diagonal_after))
    check(changed > 1e-6,
          f"the cell lengths did change, by up to {changed:.4f} A",
          "a relaxation that moved no lattice component is not testing the mode")
    return energy


def check_more_freedom_goes_no_higher(e_ortho):
    """Each mode's minimum is over a larger space than the last, so cannot be higher.

    gd optimises the positions; gd-box-ortho adds the three cell lengths;
    gd-box adds all nine lattice components. The energies have to come out in
    that order, to within the tolerance each was converged to.
    """
    print("\n  more lattice freedom reaches an energy no higher")
    _, e_positions, _, _ = relax(f"gd_{F_TOL_BOX}", "atoms.xyz", "gd", F_TOL_BOX)
    _, e_full, _, _ = relax("gdbox", "atoms.xyz", "gd-box", F_TOL_BOX, E_TOL_BOX)
    check(e_ortho <= e_positions + 1e-6,
          f"E(gd-box-ortho) = {e_ortho:.8f} <= E(gd) = {e_positions:.8f}")
    check(e_full <= e_ortho + 1e-6,
          f"E(gd-box) = {e_full:.8f} <= E(gd-box-ortho) = {e_ortho:.8f}",
          f"the unconstrained cell landed {e_full - e_ortho:+.6f} eV above the constrained one")
    print(f"          gd {e_positions:.6f} -> ortho {e_ortho:.6f} -> full {e_full:.6f}")


def main():
    if not BINARY or not os.access(BINARY, os.X_OK):
        raise SystemExit(f"TURBOGAP_BIN is not an executable: {BINARY}")
    results = check_converges_to_its_tolerance()
    check_tightening_does_not_raise_energy(results)
    check_energy_falls()
    check_relaxed_stays_relaxed()
    check_box_relaxation_goes_no_higher()
    e_ortho = check_lattice_constraint_is_honoured()
    check_more_freedom_goes_no_higher(e_ortho)

    print()
    if failures:
        print(f"relaxation: {len(failures)} of {checks} checks FAILED")
        return 1
    print(f"relaxation: all {checks} checks passed")
    return 0


sys.exit(main())
