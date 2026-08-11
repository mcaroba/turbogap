# Dipole acceptance test

Checks that TurboGAP reproduces QUIP's dipole prediction for the same model.

```sh
make DEBUG=0 all
tests/dipole/run.sh
```

Non-zero exit means TurboGAP and QUIP disagree by more than the tolerance.
Missing data is a skip, not a failure.

| variable              | meaning                                          |
| --------------------- | ------------------------------------------------ |
| `TURBOGAP_BIN`        | binary under test (default `bin/turbogap`)       |
| `TURBOGAP_DATA_ROOT`  | take the systems from here and fetch nothing |
| `TURBOGAP_TESTS_DIR`  | where the data repository is cloned          |
| `TURBOGAP_DIPOLE_TOL` | absolute tolerance (default `1e-6`)              |
| `TURBOGAP_RANKS`      | MPI ranks (default 1)                            |
| `TURBOGAP_KEEP`       | keep the staging directory                       |

## Why this is not in tests/regression

That suite compares byte-for-byte against `baseline/turbogap.e6eb1aa`, a binary
frozen before the refactor started. It has no dipole support, so it cannot
produce a reference here at all, and `diff` admits no tolerance. This test asks
a different question — do two independent codes agree on the same model — and
so needs both a real reference and a tolerance.

## What is being compared

The model is a GAP whose local "energy" is a fictitious scalar fitted so that

    mu_i = dE_i/dr_i    (derivative w.r.t. the central atom's *own* position)

is the local dipole, and `mu = sum_i mu_i`. The test compares, over 200 frames:

- per-atom `local_dipole` against QUIP
- per-frame `dipole` against QUIP
- TurboGAP's `dipole` against the sum of its own `local_dipole` columns

It does **not** compare against `mu=` in the QUIP output. That is the DFT target
the model was fitted to; comparing against it would measure the quality of the
fit rather than whether the two codes agree.

Observed agreement is ~1e-8, which is the `F16.8` print precision of
`trajectory_out.xyz`, not a numerical difference: both codes call the same
`soap_turbo` library. The 1e-6 default leaves room without being vacuous.

## Test data

Under `$TURBOGAP_DATA_ROOT/water_dipole/` (~8 MB):

| file                                       | what it is                          |
| ------------------------------------------ | ----------------------------------- |
| `water_dipole.xml`                         | the QUIP dipole GAP                 |
| `water_dipole.xml.sparseX.*`               | its two sparse sets                 |
| `water_dipole_tnep_converted_test.xyz`     | 200 frames, 2 waters, 20 A cubic    |
| `out_quip_predict`                         | QUIP's predictions (the reference)  |
| `train.sh`                                 | the `gap_fit` call that made it     |

`out_quip_predict` was produced with

```sh
quip atoms_filename=water_dipole_tnep_converted_test.xyz \
     param_filename=water_dipole.xml e \
     calc_args="dipole=dipole local_dipole=local_dipole" > out_quip_predict
```

The model is two `soap_turbo` descriptors (`central_index` 1 and 2, i.e. H and
O), `zeta=4`, `dot_product`, `delta=0.1`, `n_sparse=500`, `compress_mode=trivial`.
`run.sh` converts the xml to TurboGAP's `.gap` format on each run with
`tools/quip_xml_to_gap/make_gap_files.py --dipole-model`, rather than checking in
a derived file that could drift from the xml.

## A useful pre-check

QUIP's output also carries the model's fictitious `energy=` per frame. Running
`turbogap predict` on the converted `.gap` *without* the `--dipole-model` flag
should reproduce that as `energy_soap=`. That isolates "does TurboGAP's
soap_turbo agree with QUIP's on this model" from "is the dipole code right", and
was how this implementation was validated before a line of it was written.
