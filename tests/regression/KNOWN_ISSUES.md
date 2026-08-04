# Known issues

Defects in the vdW path found while building the refactor baseline. Issue 1 is
fixed on this branch; 2 and 3 are open and deliberately untouched, because
fixing them changes results and cannot ride on a commit whose contract is
bit-identical output.

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
are exactly the old expressions, so **the MD path is unchanged** — verified by
`vdw_tsmbd_md` and `vdw_tsmbd_md_mpi2` remaining bit-identical to the baseline.

Naming the predicate once is part of the fix, not tidying: the bug was six
copies of one condition disagreeing about whether the buffers existed.

**Coverage.** `vdw_tsmbd_predict`, `vdw_tsmbd_predict_mpi2` and `vdw_tsmbd_mc`.
These use `REFERENCE=golden` because the baseline cannot produce a reference —
it crashes. 1 and 2 ranks agree exactly on energy and to ~1e-8 on forces, the
latter being ordinary reduction-order non-associativity, which is why the
expected outputs are stored per rank count.

---

## 2. `ts+mbd` corrects energies but not forces — OPEN

On a recompute step the correction is applied to energies and to the virial,
but the corresponding force lines are commented out in the source:

```fortran
!  this_forces_vdw = this_forces_vdw + state%this_forces_vdw_corr
```

The result is that `ts+mbd` returns **MBD energies with TS-scaled forces**.
Confirmed by running the P4 dimer under `mbd` and under `ts+mbd`: energies match
to every printed digit, forces do not (e.g. 0.0867 vs 0.0800 on atom 1).

This is pre-existing and not a consequence of the issue-1 fix — MD has always
behaved this way, and the MD cases are bit-identical across that fix. It is
flagged because inconsistent energies and forces will not conserve energy in
MD, so anyone using `ts+mbd` for dynamics should know. The surrounding
`! REMOVE EVERYTHING UNDER THIS WHEN YOU ARE DONE MERGING MBD` comment suggests
the MBD merge was simply never finished.

Deciding whether to enable the force correction is a question for whoever owns
the vdW code. It would change every `ts+mbd` result, so it cannot be done on
this branch.

---

## 3. Negative "Miscellaneous" time in the timing report — OPEN, cosmetic

The end-of-run summary prints e.g. `Miscellaneous: -3.252 seconds`, because
that bucket is total minus the sum of the measured buckets and some buckets are
double-counted. Harmless, but it makes the printed totals useless for judging
refactor performance — use the wall-clock `run.sh` reports instead.
