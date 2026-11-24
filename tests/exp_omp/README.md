# The threaded pdf and diffraction forces do not depend on the thread count

On the host, the pdf kernel density estimate, the structure-factor and
diffraction forces (`get_structure_factor_forces_matrix`) and the pdf forces
(`get_pair_distribution_forces`) form each pair's kernel or force in parallel
and sum them in pair order afterwards, so they are bit-identical to the serial
loop at any thread count.

This test runs the 897-atom carbon cell of the mad_sf_matrix case, which has
enough pairs to cross several of the pdf's scratch chunks, at 1 thread, 4
threads and 2 ranks x 2 threads:

- predict: every pdf, sf, xrd and nd file must be identical across thread
  counts, and agree to 1e-7 across ranks, which reduce in another order;
- xrd and pdf MAD: the step-0 forces and virial must agree to 1e-6.

The rest of TurboGAP is not bit-reproducible across thread counts, hence the
tolerance on the MAD legs; a race in the threaded loops is an error of order
one. With the order-preserving sums made parallel on purpose, every leg fails.

```sh
make DEBUG=0 all
tests/exp_omp/run.sh
```

Non-zero exit means a leg failed. Missing test data is a skip, not a failure.

| variable             | meaning                                      |
| -------------------- | -------------------------------------------- |
| `TURBOGAP_BIN`       | binary under test (default `bin/turbogap`)   |
| `TURBOGAP_DATA_ROOT` | take the systems from here and fetch nothing |
| `TURBOGAP_TESTS_DIR` | where the data repository is cloned          |
| `TURBOGAP_KEEP`      | keep the staging directory                   |
