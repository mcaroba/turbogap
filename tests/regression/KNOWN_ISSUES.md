# Known issues

Findings from building the refactor baseline and from the driver work that
followed. Issues 1 and 4 are genuine defects and are fixed on this branch.
Issue 2 was reported as a defect, fixed, and then the fix was reverted: the
original behaviour is correct. Issue 3 is open and cosmetic. Issue 5 is a
latent trap that the current call structure makes harmless.

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

## 2. `ts+mbd` corrects energies but not forces — NOT A DEFECT

**Resolved 2026-08-04 by the owner of the vdW code: the uncorrected forces are
correct, and the force correction must stay disabled.**

This was reported as a defect and fixed in `113d0b6`; that commit has been
reverted. Do not re-enable the commented-out force lines. The description below
is kept because the observation itself is real and will be noticed again — what
was wrong was the conclusion drawn from it, not the measurement.

Note for anyone revisiting this: `113d0b6` justified enabling the correction by
measuring against `mbd_correction_freq = 1` and finding the corrected forces
~200x closer to it. That agreement is not evidence of correctness — whatever
makes the uncorrected forces right also makes `freq = 1` the wrong reference to
score them against. Re-deriving that measurement is not a reason to reopen this.

With the revert in place `vdw_tsmbd_md` and `vdw_tsmbd_md_mpi2` are back to
`REFERENCE=baseline` and are bit-identical to the pre-refactor binary, so this
path has stronger coverage than it did while the fix was applied.

---

### Original report (retained for context)

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

*(That decision has since been made — see the resolution at the top of this
section. The correction stays off.)*

---

## 3. Negative "Miscellaneous" time in the timing report — OPEN, cosmetic

The end-of-run summary prints e.g. `Miscellaneous: -3.252 seconds`, because
that bucket is total minus the sum of the measured buckets and some buckets are
double-counted. Harmless, but it makes the printed totals useless for judging
refactor performance — use the wall-clock `run.sh` reports instead.

---

## 4. `mpi_reduce` read uninitialised memory for exp-spectra forces — FIXED

**Found** 2026-08-04 while reading the reduce block for the `contrib_on` change.
Not introduced by it, and not fixed by it.

A contribution family owns a slot in `all_forces` if its *slot* predicate holds:

```fortran
contrib_on(C_PDF) = allocated(this_energies_pdf) .and. params%valid_pdf
```

but its forces are written into that slot under a **stricter** predicate:

```fortran
if( params%do_forces .and. params%exp_forces )then
   all_forces(1:3, 1:n_sites, counter2) = this_forces_pdf(1:3, 1:n_sites)
```

So with `do_forces` set and `exp_forces` unset, a valid `pdf`, `sf`, `xrd` or
`nd` family claims a slice of `all_forces` that is **never written** before

```fortran
call mpi_reduce(all_forces, all_this_forces, 3*n_sites*counter2, ...)
```

`all_forces` is `allocate`d and never zeroed, so that slice is whatever was on
the heap. It is reduced across ranks and unpacked, though under the same
stricter predicate, so the garbage is not read back into `forces_pdf`. The
observable risk is therefore not a wrong force but a trap: a signalling NaN or
a denormal in that slice would fault or stall inside `mpi_reduce`, in a spot
with no obvious connection to the exp-spectra code.

**Why it has not bitten.** Until 2026-08-04 no regression case exercised
`pdf`, `sf`, `xrd` or `nd` at all. The `xrd_predict` case now pins exactly this
configuration — `do_forces` on, `exp_forces` off, all three families valid — so
the slot bookkeeping is protected even though the uninitialised values
themselves are not observable in the compared outputs.

**Fix.** The pack loop now clears the slot of any family that owns one but
contributes no forces:

```fortran
if( contrib(i_contrib)%forces )then
   all_forces(1:3, 1:n_sites, i_contrib) = contrib(i_contrib)%f_src(1:3, 1:n_sites)
   all_virial(1:3, 1:3, i_contrib) = contrib(i_contrib)%v_src(1:3, 1:3)
else if( params%do_forces )then
   all_forces(1:3, 1:n_sites, i_contrib) = 0.d0
   all_virial(1:3, 1:3, i_contrib) = 0.d0
end if
```

Clearing only the affected slots rather than zeroing the whole array keeps the
cost off the slots that are about to be overwritten anyway.

This section previously claimed the fix could not ride on a bit-exact commit
because it changes what `mpi_reduce` sees. That was wrong, and the suite says
so: all 18 cases pass **unchanged**. Each slot reduces independently, and the
garbage was never unpacked, so replacing it with zeros moves nothing that is
read. The defect was a trap, not a wrong number — which is exactly why it
needed fixing and exactly why no test could have caught it.


---

## 5. `mpi_reduce` leaves the receive buffer undefined on non-root ranks — OPEN, by design

Noticed while fixing #4; recorded so it is not mistaken for a new defect.

`mpi_reduce` writes its receive buffer **only on the root rank**. The unpack
loop is not rank-guarded, so on every non-root rank

```fortran
contrib(i_contrib)%e_dst(1:n_sites) = all_this_energies(1:n_sites, i_contrib)
```

copies undefined heap into `energies_soap`, `energies_2b` and the rest.

This is harmless as the code stands: rank 0 owns the MD integration and
broadcasts positions and velocities back out (`turbogap.f90`, the
`mpi_bcast(positions, ...)` calls), so nothing downstream reads a non-root
rank's copy. It is recorded because the same trap as #4 applies -- a
signalling NaN or denormal landing in those arrays could fault in the "Add up
all the energy terms" arithmetic that follows, on a rank whose answer nobody
wants. `mpi_allreduce`, or guarding the unpack with `rank == 0`, would both
close it; which one is right depends on whether any future consumer needs
these per-rank, so it is left to whoever changes that.