# Known issues

Defects in the vdW path found while building the refactor baseline. Issues 1
and 2 are fixed on this branch, each in its own commit and each separate from
the pure-refactor commits, whose contract is bit-identical output. Issue 3 is
open and cosmetic.

---

## 1. `vdw_type = ts+mbd` segfaulted outside MD — FIXED

**Found** 2026-08-04 on `e6eb1aa` with `-fcheck=bounds`; **fixed** the same day.

`turbogap predict` and `turbogap mc` with `vdw_type = ts+mbd` died with:

```
At line 2364 of file src/turbogap.f90
Fortran runtime error: Array bound mismatch for dimension 1 of
array 'this_energies_vdw_corr' (1/108)
```

**Cause.** The ts+mbd scheme runs full MBD every `mbd_correction_freq` steps
and, in between, reuses a stored correction on top of a cheap TS evaluation.
Which half of that cycle a call takes was decided by
`modulo(md_istep, mbd_correction_freq) == 0`, written out at **six** separate
sites that all had to agree.

`md_istep` starts at `-1` and is only incremented under `if( params%do_md )`,
so in `predict` and `mc` it stays `-1`. With the default
`mbd_correction_freq = 100`, `modulo(-1, 100)` is `99`. Every one of the six
sites therefore chose the *reuse* half — including the site that allocates the
correction buffers, which is only reached on the *recompute* half. The consumer
then read buffers that had never been allocated. A second, quieter instance of
the same mistake: `mbd_ts_scaling` is initialised under `md_istep == 0`, so in
`predict` and `mc` it was read uninitialised.

**Fix.** The reuse half only means anything during MD, where there was a
previous step to have stored something. Outside MD there is no cycle, so every
call recomputes. The predicate is now computed once, named, and shared:

```fortran
is_correction_step = ( md_istep < 0 ) .or. &
                     ( modulo(md_istep, params%mbd_correction_freq) == 0 )
```

and the initialisation guard became `md_istep <= 0`. For `md_istep >= 0` both
are exactly the old expressions, so **the MD path was unchanged by this fix** —
verified at the time by `vdw_tsmbd_md` and `vdw_tsmbd_md_mpi2` remaining
bit-identical to the baseline. (Issue 2, committed after this one, then changed
ts+mbd forces deliberately, so those two cases are golden-compared today.)

Naming the predicate once is part of the fix, not tidying: the bug was six
copies of one condition disagreeing about whether the buffers existed.

**Coverage.** `vdw_tsmbd_predict`, `vdw_tsmbd_predict_mpi2` and `vdw_tsmbd_mc`.
These use `REFERENCE=golden` because the baseline cannot produce a reference —
it crashes. 1 and 2 ranks agree exactly on energy and to ~1e-8 on forces, the
latter being ordinary reduction-order non-associativity, which is why the
expected outputs are stored per rank count.

---

## 2. `ts+mbd` corrected energies but not forces — FIXED

The correction was applied to energies and the virial but not to forces, so
`ts+mbd` returned MBD energies with TS-scaled forces. Uncommenting the dormant
line would not have fixed it: `this_forces_vdw_corr` was snapshotted as the TS
result and never updated to the MBD one, and the step that forms
`MBD - TS_scaled` and lands on MBD did not exist for forces at all. The fix
mirrors the energy path exactly, at all three points.

**Measured on the 108-atom P4 cell**, against `mbd_correction_freq = 1`
(true MBD at every step) as reference:

| | max abs force error at first reuse step | after 10 reuse steps |
|---|---|---|
| without the correction | 1.9e-05 (0.33% of \|F\|max) | — |
| with the correction | 9.0e-08 (0.00%) | 8.8e-06 (0.15%) |

A factor of ~200 at the first reuse step. The corrected forces after *ten*
reuse steps are still better than the uncorrected forces after *one*. The
residual grows smoothly as the frozen correction ages, which is the expected
behaviour of the scheme and what `mbd_correction_freq` trades against cost.

On a recompute step `ts+mbd` now reproduces `mbd` exactly, forces included —
previously only the energies matched.

The cost argument the scheme exists for, same system, 10 MD steps:
`mbd_correction_freq = 1` takes 102.0 s, `= 1000` takes 8.6 s — **11.9x**.

**Consequence for the suite.** This deliberately changes ts+mbd results, so
`vdw_tsmbd_md` and `vdw_tsmbd_md_mpi2` could no longer be compared against the
baseline binary and became `REFERENCE=golden` like the other ts+mbd cases.
That is a genuine loss of safety for those two, and the price of the change.
It is contained: the five ts+mbd cases were the only ones that moved, and
`ts`, `mbd` and every non-vdW case remain bit-identical to the baseline.

Still dormant, and untouched: the `S_xyz_inv` strain term
(`forces_vdw = forces_vdw - local_virial_vdw_diag_corr * S_xyz_inv`) and the
`this_local_virial_vdw_diag` correction. Those are a separate correction from
the one fixed here, and `local_virial_vdw_diag_corr` already has a live role
feeding the scaling update, so enabling them is its own question.

---

## 3. Negative "Miscellaneous" time in the timing report — OPEN, cosmetic

The end-of-run summary prints e.g. `Miscellaneous: -3.252 seconds`, because
that bucket is total minus the sum of the measured buckets and some buckets are
double-counted. Harmless, but it makes the printed totals useless for judging
refactor performance — use the wall-clock `run.sh` reports instead.
