# Cartesian gradient of the dipole

`d(mu_a)/d(r_jb)` by central differences, for a GAP dipole model.

```sh
tools/dipole_gradient/dipole_gradient.py \
    --gap gap_files/water_dipole.gap \
    --atoms frame.xyz \
    --turbogap bin/turbogap
```

Writes `dipole_gradient.dat` (one line per atom, per `a`, per `b`) and prints
the translational sum rule as a correctness check.

## Why this is finite difference and not analytic

A GAP dipole model predicts `mu_i = dE_i/dr_i`, so its Cartesian gradient is a
**second** derivative of the fitted scalar. Expanding through the descriptor
`q`:

```
d(mu_ia)/d(r_jb) = sum_de (d2E_i/dq_d dq_e) (dq_d/dr_ia) (dq_e/dr_jb)   [A]
                 + sum_d  (dE_i/dq_d) (d2q_d / dr_ia dr_jb)             [B]
```

**Term A** is computable from quantities TurboGAP already forms. For the
`dot_product` kernel with exponent `zeta`,

```
d2E_i/dq dq = delta^2 zeta (zeta-1) sum_s alpha_s (q_i.q_s)^(zeta-2) q_s q_s^T
```

so term A is `delta^2 zeta (zeta-1) sum_s alpha_s k_s^(zeta-2) (q_s . dq/dr_ia)
(q_s . dq/dr_jb)`, i.e. a contraction of the sparse set against
`soap_cart_der`. It costs roughly `n_sparse` times the existing force
contraction, because the sparse index no longer collapses.

**Term B** requires `d2q/dr dr`, the second derivative of the SOAP descriptor.
`soap_turbo` does not provide it — `get_soap` returns `soap` and
`soap_cart_der` and nothing beyond. Producing it means differentiating
`soap_turbo_radial.f90`, `soap_turbo_angular.f90` and `soap_turbo_compress.f90`
a second time and threading a rank-higher array through `get_soap`. That is a
real project, not an afternoon.

So there is **no exact analytic gradient available today**, and term A alone is
not the gradient — nothing bounds term B in general. Central differences need
none of that machinery and are exact to `O(h^2)`.

## Cost

`6N+1` configurations. TurboGAP evaluates every frame of an xyz in a single
invocation, so this is **one run**, not `6N` runs. Fine for validation and for
small systems; not a production route for MD on large cells.

## Choosing the step

`dipole=` is printed to `F16.8`, which puts a floor under the differencing:
rounding contributes `~5e-9 / 2h` per element while truncation contributes
`O(h^2)`. Measured on a 2-water frame (largest `|J|` element 28.7), using the
translational sum rule as the error estimate:

| h      | max sum-rule residual | relative |
| ------ | --------------------- | -------- |
| 1e-4   | 1.0e-4                | 3.5e-6   |
| 5e-4   | **2.0e-5**            | **7.0e-7** |
| 1e-3   | 5.5e-5                | 1.9e-6   |
| 5e-3   | 1.3e-3                | 4.6e-5   |
| 1e-2   | 5.3e-3                | 1.8e-4   |

`5e-4` is the optimum and is the default. The `1/h` branch below it and the
`h^2` branch above it are both clearly visible, which is itself evidence the
differencing is behaving.

## The sum rule

`mu` is invariant under a rigid translation of the whole system, so

```
sum_j d(mu_a)/d(r_jb) = 0    for every a, b
```

This is a property of the model, not of the differencing, so it is a genuine
check: a residual far above the rounding floor means something is wrong. The
script reports it on every run.

## If you need this in production (MAD / IR forces)

Experimental forces from an IR spectrum — matching a DACF-derived spectrum to
measurement — need `d(mu)/d(r)` at every MD step, where `6N+1` evaluations per
step is not viable. Two routes:

1. **Implement `d2q/dr dr` in soap_turbo** and then term A + term B
   analytically. Exact and fast to evaluate, but the differentiation work is
   substantial and the memory for a `(3, 3, n_soap, n_pairs)` object is large.

2. **Reformulate the dipole so its gradient is first-order.** If the model
   predicts per-atom charges `q_i` instead, with `mu = sum_i q_i r_i`, then

   ```
   d(mu_a)/d(r_jb) = q_j delta_ab + sum_i (dq_i/dr_jb) r_ia
   ```

   and `dq_i/dr_jb` is exactly `local_properties_cart_der`, which TurboGAP
   **already computes** — it is what the `atomic_charge` local property feeds to
   electrostatics. No second derivatives anywhere.

   The catch is that `sum_i q_i r_i` is origin-dependent unless the system is
   neutral, and it is a less expressive model than the `dE_i/dr_i` form. A
   hybrid (charges plus a per-atom vector correction) is the usual compromise.

Route 2 reuses machinery that exists and is already exercised by the
electrostatics path; route 1 is the one to take only if the accuracy of the
`dE_i/dr_i` form turns out to be needed. Worth deciding before building either.
