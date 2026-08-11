#!/usr/bin/env python3
"""Build the synthetic fit targets the virial legs are driven by.

The virial test needs a target the prediction does not already sit on: if
y_exp were the prediction, every residual would be zero and the whole gradient
would vanish, which no amount of wrong algebra could fail. An affine distortion
of the prediction gives a residual with structure across the whole grid while
staying on the same abscissa, so no interpolation enters and the test says
nothing about the experimental reader.

That needs the prediction first, and on this branch the prediction cannot be
had without a target: `perform%pdf` is `do_pair_distribution .and. valid_pdf`,
and `valid_pdf` is set only by a dataset labelled `pair_distribution` being
declared. So the run that produces the prediction is itself driven by a flat
seed, whose values never reach anything -- the pattern turbogap writes out does
not depend on what it is being fitted to.

Usage:
  make_targets.py <prediction.dat> <target.dat>        distort a prediction
  make_targets.py --flat <x0> <x1> <n> <target.dat>    seed a grid to bootstrap
"""
import sys

import numpy as np

SLOPE, OFFSET = 0.85, 0.05


def main():
    if sys.argv[1] == "--flat":
        x0, x1, n, dst = float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
        d = np.column_stack([np.linspace(x0, x1, n), np.ones(n)])
    else:
        src, dst = sys.argv[1], sys.argv[2]
        d = np.loadtxt(src)
        d[:, 1] = SLOPE * d[:, 1] + OFFSET
    np.savetxt(dst, d, fmt="%.10f")


if __name__ == "__main__":
    main()
