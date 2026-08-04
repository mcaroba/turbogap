# Known issues found while building the regression suite

These are pre-existing defects in master, found while establishing the
refactor baseline. They are **not** caused by the refactor and are **not**
fixed on this branch — fixing them would change results and so cannot ride
along with commits whose contract is bit-identical output.

---

## 1. `vdw_type = ts+mbd` segfaults in `predict` and `mc` modes

**Found:** 2026-08-04, on `e6eb1aa` (master HEAD) with `-fcheck=bounds`.

**Reproducer** (any system with a `hirshfeld_v` local property):

```sh
cd tests/regression && sed 's/^vdw_type = mbd/vdw_type = ts+mbd/' \
  cases/vdw_mbd_cell_mpi2/input > /tmp/in
# stage p4_cell.xyz + gap_files alongside, then:
turbogap predict          # -> SIGSEGV
```

Bounds-checked build reports:

```
At line 2364 of file src/turbogap.f90
Fortran runtime error: Array bound mismatch for dimension 1 of
array 'this_energies_vdw_corr' (1/108)
```

**Cause.** Two guards on the same expression must agree, and in `predict`
they both take the branch that assumes the other one ran:

* `turbogap.f90:909` sets `md_istep = -1`. It is only incremented inside
  `if( params%do_md )`, so in `predict` and `mc` it stays `-1`.
* `types.f90:139` defaults `mbd_correction_freq = 100`, so
  `modulo(-1, 100) = 99`.
* `turbogap.f90:2023` allocates `this_energies_vdw_corr(1:n_sites)` only when
  `vdw_type == "ts+mbd" .and. modulo(md_istep, mbd_correction_freq) == 0`.
  With `md_istep = -1` this is false, so the array is never allocated.
* `turbogap.f90:2279` branches on the *same* expression and therefore takes
  its `else` at line 2360, which at 2364 evaluates
  `this_energies_vdw = this_energies_vdw + this_energies_vdw_corr` —
  reading the array that was never allocated.

`md` mode is unaffected: step 0 gives `modulo(0, 100) == 0`, so the array is
allocated on the first step and the `else` branch at later steps reuses it,
which is the intended design.

**Why no fix here.** The obvious repair — allocate whenever
`vdw_type == "ts+mbd"` — is not behaviour-preserving. The current code
deallocates, reallocates and re-zeros every `mbd_correction_freq` steps, and
the `else` branch depends on the value *persisting* between those points.
Making allocation unconditional would re-zero the correction every step and
silently change MD results. Deciding what the correction should do when there
is no MD step counter is a question for whoever owns the vdW code, not
something a refactor commit should answer.

**Consequence for the suite.** `cases/vdw_mbd_cell_mpi2` uses `vdw_type = mbd`
instead, and `ts+mbd` is covered from `md` mode instead of `predict`, by
`cases/vdw_tsmbd_md{,_mpi2}` — `md` is unaffected by this bug, so the
correction machinery *is* exercised, including the reuse branch that is
unreachable in `predict`. What remains uncovered is specifically `ts+mbd`
under `predict` and `mc`, which is the crash itself. Once the bug is fixed,
switch `vdw_mbd_cell_mpi2` back to `ts+mbd`.

---

## 2. Negative "Miscellaneous" time in the timing report

Cosmetic. The end-of-run summary prints e.g. `Miscellaneous: -3.252 seconds`
because the miscellaneous bucket is computed as total minus the sum of the
measured buckets, and some buckets are double-counted. Harmless, but it makes
the printed totals useless for judging refactor performance — use the
wall-clock times reported by `run.sh` instead.
