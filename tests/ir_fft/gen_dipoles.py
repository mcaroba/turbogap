#!/usr/bin/env python3
"""A synthetic dipole trajectory with water-like structure, for the ir_fft checks.

Five damped modes at wavenumbers a water spectrum has features at, each with an
independent random phase per Cartesian component, plus white noise and a large
STATIC OFFSET. The offset is the point of the last term: it exercises
subtract_mean, and without it the mean-subtraction path and the gradient's
prefix-sum correction would both be tested on a series whose mean is already
zero -- which is to say, not tested.

Seeded, so the comparison with the Python is on identical numbers.
"""
import sys
import numpy as np

out = sys.argv[1] if len(sys.argv) > 1 else "dipoles.dat"
rng = np.random.default_rng(20260824)

T, dt = 6000, 1.0                      # frames, fs
CM = 33356.40952                       # cm^-1 per (1/fs)
t = np.arange(T) * dt
mu = np.zeros((T, 3))
for nu, amp, tau in [(600., 1.0, 300.), (1650., 0.6, 250.), (3400., 0.9, 180.),
                     (3600., 0.4, 150.), (150., 1.4, 500.)]:
    w = 2 * np.pi * nu / CM            # rad/fs
    for d in range(3):
        ph = rng.uniform(0, 2 * np.pi)
        # Re-excited every 2 ps, so the correlation decays without the whole
        # series decaying to nothing over 6 ps.
        mu[:, d] += amp * np.cos(w * t + ph) * np.exp(-((t % 2000) / tau))
mu += 0.15 * rng.standard_normal((T, 3))
mu += np.array([2.5, -1.0, 0.7])

np.savetxt(out, mu, fmt="%.17e")
print(f"wrote {out}: {T} frames at {dt} fs")
