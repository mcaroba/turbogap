# The Lorentz-polarization factor

`xrd_lorentz_polarization = .true.` multiplies the predicted diffraction
pattern by

```
LP(theta) = ( 1 + P cos^2(2 theta) ) / ( sin^2(theta) cos(theta) )
```

with `P = xrd_lp_polarization` the degree of polarization of the incident
beam: 1 for an unpolarized source, `cos^2(2 theta_M)` behind a monochromator
set at `2 theta_M`. The scattering angle comes from the q grid and
`xrd_wavelength`, since the grid holds `2 sin(theta)/lambda` whatever
`q_units` the deck was written in.

Two parts of the grid have no usable angle and are set to zero rather than to
an infinity: below `xrd_lp_sin_theta_min`, where the Lorentz factor diverges
as the forward direction is approached, and at or beyond the Ewald limit
`sin(theta) >= 1`, which cannot be measured at that wavelength at all. On a
grid given in `q` that second cut bites earlier than it looks: at Cu K-alpha
it removes everything above `Q = 4 pi / lambda ~ 8.16 1/A`.

The factor applies to whatever `xrd_output` asked for. It corrects a raw
powder intensity, so `xrd_output = 'xrd'` is what it is meant for; applied to
`F(q)` it is just a per-q rescaling of an already normalised quantity.

It applies to neutron diffraction too, on the same keyword. The polarization
half of the expression is an X-ray notion — for neutrons `P = 0` leaves the
Lorentz factor alone — so `xrd_lp_polarization` is a modelling choice there
rather than a property of the beam, and one run predicting both XRD and ND
weights them with the same `P`.

## Both routes apply it

The two routes to a pattern (see [xrd_debye.md](xrd_debye.md)) apply the
factor in different places, and both keep the gradient consistent with the
energy. `LP` is a weight per q sample that does not depend on the positions,
so the chain rule carries it through untouched: `dy/dr` picks up exactly the
same number `y` does, and the MAD forces stay the forces of the MAD energy.

**Debye route** — `calculate_xrd_debye` scales `y_xrd` and the gradient
channel by `lp` together, right after the affine map that implements
`xrd_output`.

**pdf route** — `calculate_xrd` scales the reduced pattern, and hands `lp` to
`get_structure_factor_forces` and `get_structure_factor_forces_matrix`, which
fold it into `all_scattering_factors`: the per-q weight that already turns
`dS_ab/dr` into `dy/dr`. Three details are worth knowing:

- It is applied **after** the MPI reduce, not inside
  `get_xrd_from_partial_structure_factors`. Each rank builds only its own
  slice `[q_beg, q_end]` of the grid, so weighting inside that routine would
  work but would put the factor in the middle of the q decomposition; doing it
  once on the whole pattern is what makes the MPI leg of the test meaningful.
- Being downstream of the branch, it covers the `get_xrd_explicit` fallback
  (`structure_factor_from_pdf = .false.`) as well.
- The self term and the `+ 1` that `xrd_output` conventions add are constants
  of the geometry, so they are inside the weight, not outside it, and they
  contribute nothing to the gradient either way.

`lp` stays allocated and equal to one when the keyword is off, so the force
routines need no second call form. Multiplying by `1.0d0` is exact, and the
regression suite confirms bit-identical output with the factor off.

## Keywords

| keyword                    | default     | meaning                                        |
| -------------------------- | ----------- | ---------------------------------------------- |
| `xrd_lorentz_polarization` | `.false.`   | apply the powder LP factor                     |
| `xrd_lp_polarization`      | `1.0`       | `P` in `1 + P cos^2 2theta`                    |
| `xrd_lp_sin_theta_min`     | `1e-3`      | below this `sin(theta)`, LP is set to zero     |
| `xrd_wavelength`           | `1.5405981` | sets the angle each q sample corresponds to    |

## Tests

- `tests/xrd_debye/run.sh` — the factor on the Debye route, against a pattern
  rebuilt from scratch in numpy at `P = 0.8`, plus finite-difference forces
  and virial.
- `tests/xrd_lp_pdf/run.sh` — the factor on the pdf route: the weight itself
  over four decks including MPI, and finite-difference forces and virial. That
  the factor reaches the virial at all is worth a word: writing this test is
  what turned up the pdf route's virial disagreeing with a strain derivative of
  its own energy, for reasons older than the factor. `tests/pdf_virial` covers
  that.
