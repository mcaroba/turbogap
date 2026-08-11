#!/usr/bin/env python3
"""Finite-difference check of the Debye XRD/ND forces and virial.

The experimental term is isolated by running each geometry twice, once at the
deck's energy scale and once at zero. Everything but the XRD contribution is
identical between the two and cancels, so what is left is exactly the family
under test -- no dependence on the GAP forces being right, and no need to
print an extra diagnostic from the Fortran.

Forces:  F_a = -dE/dr_a, checked by displacing one atom one Cartesian
         component at a time.
Virial:  W_kl = -dE/deps_kl, checked by applying a symmetric homogeneous
         strain to the cell and every position. This is the convention
         TurboGAP writes out (its stress is -virial/volume).

Usage:  fd_gradient.py <workdir> --bin <turbogap>
"""
import argparse
import re
import shutil
import subprocess
import sys

import numpy as np


def _need(match, what):
    if match is None:
        raise SystemExit(what)
    return match


def read_frame(path):
    lines = open(path).read().splitlines()
    n = int(lines[0].split()[0])
    return n, lines[1], lines[2:2 + n]


def write_frame(path, n, comment, body):
    with open(path, "w") as fh:
        fh.write(f"{n:>8d}\n{comment}\n" + "\n".join(body) + "\n")


class Runner:
    """Runs turbogap in a staging directory and reads back what it produced."""

    def __init__(self, binary, workdir, atoms, scale, mode="predict", family="xrd"):
        self.binary, self.workdir, self.atoms, self.mode = binary, workdir, atoms, mode
        self.scale = scale
        # Which per-family energy in the output frame this deck is exercising.
        # The forces and the virial are read from the same frame whatever it
        # is, because contribution() differences two energy scales; only the
        # energy has to be named, and naming the family rather than the total
        # keeps the finite difference away from the GAP energy's round-off.
        self.family = family
        self.deck = open(f"{workdir}/input").read()
        self.n, self.comment, self.body = read_frame(f"{workdir}/{atoms}")
        self.lattice = np.array(
            [float(v) for v in _need(re.search(r'Lattice="([^"]+)"', self.comment),
                                     "no Lattice in the reference frame").group(1).split()]
        ).reshape(3, 3)
        self.calls = 0

    def deck_at_scale(self, scale):
        return re.sub(r"(?m)^\s*exp_energy_scales\s*=.*$",
                      f"exp_energy_scales = {scale}", self.deck)

    def geometry(self, strain=None, atom=None, dim=0, h=0.0):
        """Write the working atoms file: reference geometry, optionally with a
        homogeneous strain applied to cell and positions, or one atom nudged."""
        comment, body = self.comment, list(self.body)
        if strain is not None:
            d = np.eye(3) + strain
            lat = self.lattice @ d.T
            comment = re.sub(r'Lattice="[^"]+"',
                             'Lattice="' + " ".join(f"{float(v):.12f}" for v in lat.flatten()) + '"',
                             comment)
            out = []
            for line in body:
                w = line.split()
                p = np.array([float(v) for v in w[1:4]]) @ d.T
                w[1:4] = [f"{float(v):.14f}" for v in p]
                out.append(" ".join(w))
            body = out
        if atom is not None:
            w = body[atom].split()
            w[1 + dim] = f"{float(w[1 + dim]) + h:.14f}"
            body[atom] = " ".join(w)
        write_frame(f"{self.workdir}/{self.atoms}", self.n, comment, body)

    def run(self, scale, label):
        open(f"{self.workdir}/input", "w").write(self.deck_at_scale(scale))
        r = subprocess.run([self.binary, self.mode], cwd=self.workdir,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        self.calls += 1
        if r.returncode != 0:
            print(r.stdout.decode()[-2000:])
            raise SystemExit(f"turbogap failed at {label}")
        lines = open(f"{self.workdir}/trajectory_out.xyz").read().splitlines()
        c = lines[1]
        m = _need(re.search(rf"energy_{self.family}=(\S+)", c),
                  f"no energy_{self.family} in the output frame; "
                  "does the deck set exp_energies?")
        energy = float(m.group(1))
        virial = np.array([float(v) for v in _need(
            re.search(r'virial="([^"]+)"', c), "no virial in the output frame").group(1).split()]).reshape(3, 3)
        forces = np.array([[float(v) for v in ln.split()[4:7]] for ln in lines[2:2 + self.n]])
        return energy, forces, virial

    def contribution(self, label, **geom):
        """Energy, forces and virial of the experimental term alone."""
        self.geometry(**geom)
        e0, f0, v0 = self.run("0.0", f"{label} (scale 0)")
        self.geometry(**geom)
        e1, f1, v1 = self.run(self.scale, label)
        return e1 - e0, f1 - f0, v1 - v0


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("workdir")
    p.add_argument("--bin", required=True)
    p.add_argument("--atoms", default="atoms.xyz")
    p.add_argument("--scale", default="100.0",
                   help="energy scale of the deck under test")
    p.add_argument("--family", default="xrd", choices=("xrd", "nd", "pdf", "sf"),
                   help="which energy_* in the output frame the deck drives")
    p.add_argument("--h", type=float, default=1e-3, help="displacement, Angstrom")
    p.add_argument("--strain", type=float, default=1e-4)
    p.add_argument("--atoms-to-check", type=int, default=3)
    p.add_argument("--seed", type=int, default=0)
    # An escape hatch for a deck whose energy is known not to be a function of
    # the cell alone -- nothing in the tree needs it today. The pdf route used
    # to: its pair distribution is normalised by the number density N/V, and
    # the virial did not carry the resulting volume term. tests/pdf_virial now
    # checks that it does.
    p.add_argument("--skip-virial", action="store_true",
                   help="check the forces only, not the virial")
    # The reported energy carries 8 decimals, so at E ~ 1e3 and h = 1e-3 the
    # finite difference is only good to a few parts in 1e6.
    p.add_argument("--tol", type=float, default=1e-4)
    a = p.parse_args()

    shutil.copy(f"{a.workdir}/{a.atoms}", f"{a.workdir}/atoms_reference.xyz")
    r = Runner(a.bin, a.workdir, a.atoms, a.scale, family=a.family)
    deck = r.deck
    status = 0
    try:
        e_ref, f_ref, v_ref = r.contribution("reference")
        f_scale = max(float(np.abs(f_ref).max()), 1e-12)
        v_scale = max(float(np.abs(v_ref).max()), 1e-12)
        print(f"    energy_{a.family} = {e_ref:.8f} eV   "
              f"max|F| = {f_scale:.4e} eV/A   max|W| = {v_scale:.4e} eV")

        drift = float(np.abs(f_ref.sum(axis=0)).max())
        print(f"    sum of forces = {drift:.2e} eV/A (translation invariance)")
        if drift > 1e-6 * f_scale * r.n:
            print("    FAIL: forces do not sum to zero")
            status = 1

        rng = np.random.default_rng(a.seed)
        picks = rng.choice(r.n, size=min(a.atoms_to_check, r.n), replace=False)

        print(f"\n    forces, h = {a.h} A")
        print(f"    {'atom':>5} {'dim':>3} {'analytic':>15} {'finite diff':>15} {'rel':>10}")
        worst = 0.0
        for atom in sorted(int(i) for i in picks):
            for dim in range(3):
                e = [r.contribution(f"atom {atom} dim {dim} {s:+}",
                                    atom=atom, dim=dim, h=s * a.h)[0]
                     for s in (+1.0, -1.0)]
                fd = -(e[0] - e[1]) / (2.0 * a.h)
                an = f_ref[atom, dim]
                rel = abs(fd - an) / f_scale
                worst = max(worst, rel)
                print(f"    {atom:>5} {dim:>3} {an:>15.8f} {fd:>15.8f} {rel:>10.2e}")
        print(f"    worst relative deviation: {worst:.2e}  (tolerance {a.tol:.1e})")
        if worst > a.tol:
            print("    FAIL: forces disagree with finite differences")
            status = 1

        if a.skip_virial:
            print("\n    virial: SKIPPED (--skip-virial)")
        else:
            print(f"\n    virial, strain = {a.strain}")
            print(f"    {'kl':>5} {'analytic':>15} {'finite diff':>15} {'rel':>10}")
            worst = 0.0
            for k, l in [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]:
                e = []
                for s in (+1.0, -1.0):
                    st = np.zeros((3, 3))
                    st[k, l] += s * a.strain / 2.0
                    st[l, k] += s * a.strain / 2.0
                    e.append(r.contribution(f"strain {k}{l} {s:+}", strain=st)[0])
                fd = -(e[0] - e[1]) / (2.0 * a.strain)
                an = v_ref[k, l]
                rel = abs(fd - an) / v_scale
                worst = max(worst, rel)
                print(f"    {k}{l:>4} {an:>15.6f} {fd:>15.6f} {rel:>10.2e}")
            print(f"    worst relative deviation: {worst:.2e}  (tolerance {a.tol:.1e})")
            if worst > a.tol:
                print("    FAIL: virial disagrees with finite differences")
                status = 1
    finally:
        shutil.move(f"{a.workdir}/atoms_reference.xyz", f"{a.workdir}/{a.atoms}")
        open(f"{a.workdir}/input", "w").write(deck)

    print(f"\n    {r.calls} turbogap invocations")
    return status


if __name__ == "__main__":
    sys.exit(main())
