# XRD and ND from the Debye scattering equation

TurboGAP has two ways to predict an X-ray or neutron diffraction pattern. This
document describes the second one, added by `xrd_debye = .true.`, and how it
differs from the one that was there before.

## The two routes

**The pdf route** (the default, `xrd_debye = .false.`) is the chain described
at the top of `compute_exp_spectra`: bin the pair distances into partial pair
distribution functions `g_ab(r)`, Fourier transform those into partial
structure factors `S_ab(q)`, then weight by the form factors to get `I(q)`.
Its cost is set by the number of pairs inside `pair_distribution_rcut` and by
`pair_distribution_n_samples`, so it scales linearly in the number of atoms
and is the right choice for a large cell.

**The Debye route** goes straight from the positions:

```
I(q) = 1/N sum_i sum_j f_i(q) f_j(q) sinc(q r_ij) w(r_ij)
```

summed over every pair in the cell, `i = j` included — those terms contribute
`f_i(q)^2`, the self term the pdf route has to add back by hand. There is no
binning, so nothing to converge in `r`: no `pair_distribution_n_samples`, no
kernel width, no interpolation onto a grid. The price is that it is
`O(N^2 n_samples)`, which is fine at a few thousand atoms and hopeless at a
hundred thousand.

Because the Debye route needs neither the pair distribution nor the structure
factor, `read_input_file` switches both back off when `xrd_debye` is on —
`do_xrd` and the `exp_labels` block turn them on as a side effect, and leaving
them on would make the keyword a pure slowdown. A deck that asks for either in
as many words keeps it:

```
do_xrd               = .true.
xrd_debye            = .true.
do_pair_distribution = .true.    ! computed anyway, because you asked
```

`xrd_rcut` means something different on the two routes. On the pdf route it
widens the neighbour list. On the Debye route no neighbour list is involved:
it is the cutoff of the Lorch window below, and `turbogap_setup` deliberately
does not widen `rcut_max` for it.

## Periodicity

The Debye sum runs over the atoms as given, with no minimum-image convention
and no periodic images. That is exact for a cluster or a nanoparticle. For a
periodic cell it is an approximation, and truncating a sum of sinusoids puts
ringing into `I(q)`; `structure_factor_window = .true.` applies the Lorch
window `sinc(pi r / r_cut)` out to `r_cut = xrd_rcut` and drops pairs beyond
it, which is what damps that ringing. The window's own derivative is carried
through into the forces, so turning it on does not cost gradient consistency.

## The Lorentz-polarization factor

`xrd_lorentz_polarization = .true.` multiplies the predicted pattern by

```
LP(theta) = ( 1 + P cos^2(2 theta) ) / ( sin^2(theta) cos(theta) )
```

with `P = xrd_lp_polarization` the degree of polarization of the incident
beam: 1 for an unpolarized source, `cos^2(2 theta_M)` behind a monochromator
set at `2 theta_M`. The scattering angle comes from the q grid and
`xrd_wavelength`, since the grid holds `2 sin(theta)/lambda`.

Two parts of the grid have no usable angle and are set to zero rather than to
an infinity: below `xrd_lp_sin_theta_min`, where the Lorentz factor diverges
as the forward direction is approached, and above the Ewald limit
`sin(theta) >= 1`, which cannot be measured at that wavelength at all.

The factor is applied on the Debye route only, and to whatever `xrd_output`
asked for. It corrects a raw powder intensity, so `xrd_output = 'xrd'` is what
it is meant for; applied to `F(q)` it is just a per-q rescaling of an already
normalised quantity. It is a multiplicative weight per q sample either way, so
the energy and the gradient stay consistent with each other.

The pdf route does not apply it. Making it consistent there means touching
`get_structure_factor_forces` and the self-term bookkeeping in
`get_xrd_from_partial_structure_factors`, and no test covers that today.

## Output conventions

`xrd_output` (and `nd_output`) mean exactly what they mean on the pdf route,
whose `y` is the interference part `I - sum_i c_i f_i^2`:

| `xrd_output`        | predicted quantity                            |
| ------------------- | --------------------------------------------- |
| `xrd`               | `I(q)`, the raw intensity per atom            |
| `F(q)`, `i(q)`      | `(I - sum_i c_i f_i^2) / (sum_i c_i f_i)^2 + 1` |
| `q*F(q)`, `q*i(q)`  | `2 pi q (I - sum_i c_i f_i^2) / (sum_i c_i f_i)^2` |

Each is an affine map of `I(q)` with coefficients that do not depend on the
positions, which is why `calculate_xrd_debye` can carry the gradient through
all of them by scaling one channel per q sample.

## Energies and forces

Identical in form to the pdf route:

```
E = 1/2 s sum_l ( y_l - y_exp_l )^2,   F_a = -s sum_l ( y_l - y_exp_l ) dy_l/dr_a
```

with `s = exp_energy_scales` for this observable. The virial is
`sum_i F_i (x) r_i`, symmetrised — valid here because the energy depends on
the positions only through distances between atoms of the cell, with no
periodic images in the sum, so that expression telescopes into the pair form.

Under MPI the double sum is split over the outer atom index. The pattern is
reduced; the per-atom derivatives are not, because each rank's slice is
already complete for the atoms it owns — that rank looped over every `j`.

## Keywords

| keyword                    | default     | meaning                                        |
| -------------------------- | ----------- | ---------------------------------------------- |
| `xrd_debye`                | `.false.`   | use the Debye sum for XRD and ND               |
| `xrd_lorentz_polarization` | `.false.`   | apply the powder LP factor (Debye route only)  |
| `xrd_lp_polarization`      | `1.0`       | `P` in `1 + P cos^2 2theta`                    |
| `xrd_lp_sin_theta_min`     | `1e-3`      | below this `sin(theta)`, LP is set to zero     |
| `xrd_rcut`                 | `4.0`       | Lorch window cutoff when `xrd_debye` is on     |
| `structure_factor_window`  | `.true.`    | apply the Lorch window                         |

## Tests

- `tests/xrd_debye/run.sh` — the acceptance test. Checks the predicted pattern
  against an independent Debye sum written from scratch in numpy, across the
  output conventions, the window, the LP factor, neutrons and MPI; then checks
  the forces and the virial against central finite differences of the energy.
  This is what says the physics is right.
- `tests/regression/cases/xrd_debye_mad{,_mpi2}` — golden characterization
  cases. They pin the output against drift and assert nothing about
  correctness; the baseline binary predates the keyword and exits on it.
