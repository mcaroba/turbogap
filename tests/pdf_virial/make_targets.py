#!/usr/bin/env python3
"""Turn a predicted pattern into a synthetic fit target.

The virial test needs a target the prediction does not already sit on: if
y_exp were the prediction, every residual would be zero and the whole gradient
would vanish, which no amount of wrong algebra could fail. An affine distortion
of the prediction gives a residual with structure across the whole grid while
staying on the same abscissa, so no interpolation enters and the test says
nothing about the experimental reader.

Usage: make_targets.py <prediction.dat> <target.dat>
"""
import sys

import numpy as np

SLOPE, OFFSET = 0.85, 0.05


def main():
    src, dst = sys.argv[1], sys.argv[2]
    d = np.loadtxt(src)
    d[:, 1] = SLOPE * d[:, 1] + OFFSET
    np.savetxt(dst, d, fmt="%.10f")


if __name__ == "__main__":
    main()
