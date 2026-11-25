# XRD and ND do not depend on atom order

The partial structure factors are indexed by species id, so slot `k` must be
scattered with species `k`'s form factor or scattering length. The slot names
used to come from the first atom to raise the slot count. If the species did
not first appear in the order of `species =`, Te was scattered with Sb's form
factor (GST listed Ge, Sb, Te with `species = Ge Te Sb`), and the partial
structure factor files overwrote one another.

This test runs the xrd_mad cell (`species = C O`) twice: atoms in file order,
C first, and with the O atoms moved to the front. The O-first run is done at
one rank and at several. Every pdf, sf, xrd and nd output must agree to 1e-6,
and the partials must be written as C_C, C_O and O_O. On 9ffd623, before the
fix, every O-first leg fails.

A species can be listed and absent from the structure, as when hybrid MC
removes it. Its partials have no atoms to be normalised by, and 0/0 there made
every total NaN. The last leg runs the cell with its O removed as
`species = O C` (O absent and first) against `species = C`: the totals and the
C_C partials must agree, with no NaN in any output, and so must every energy
term: with `species = C` the potential's O descriptors match no listed species,
and the device build counted them as another species pair's (GST without Sb,
`species = Ge Te`: 2b -10840 eV against -4450). On e28f1d1, before these
fixes, the leg fails on the host build (NaN) and on the device build (3b).

```sh
make DEBUG=0 all
tests/xrd_species_order/run.sh
```

Non-zero exit means a leg failed. Missing test data is a skip, not a failure.

| variable             | meaning                                      |
| -------------------- | -------------------------------------------- |
| `TURBOGAP_BIN`       | binary under test (default `bin/turbogap`)   |
| `TURBOGAP_DATA_ROOT` | take the systems from here and fetch nothing |
| `TURBOGAP_TESTS_DIR` | where the data repository is cloned          |
| `TURBOGAP_RANKS`     | ranks for the MPI leg (default 3)            |
| `TURBOGAP_KEEP`      | keep the staging directory                   |
