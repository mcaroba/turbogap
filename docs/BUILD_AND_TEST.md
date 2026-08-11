# Build and test: the one-liners

Everything here is copy-pasteable. `$M` and `$G` are the two trees:

```sh
export M=/u/74/zarrout1/unix/work/cpu_vs_gpu_tests/turbogap_master_2026
export G=/u/74/zarrout1/unix/work/cpu_vs_gpu_tests/turbogap_gpu_commit_mahti
```

---

## 0. Once per clone

```sh
tools/setup_dev_env.sh                 # venv + pinned fprettify + git hook
export PATH="$HOME/.venvs/turbogap-tools/bin:$PATH"
tools/setup_dev_env.sh --check         # exits non-zero if anything is missing
```

`tests/fetch_test_data.sh` runs automatically on the first suite run and clones
`https://github.com/TiganyZ/turbogap_tests` to `<repo>/../turbogap_tests`.

---

## 1. Build

```sh
# CPU
cd $M && make -j12                                    # ~26 s

# GPU -- DEBUG defaults to 1, which is a 2.1x slowdown. Always pass DEBUG=0.
cd $G && export HOP_ROOT=/u/74/zarrout1/unix/work/hop \
      && export TURBOGAP_ARCH=Aalto_gfortran_openblas_hip_cuda \
      && make -j12 DEBUG=0                            # ~15 s
```

### Dependencies -- REGENERATE AFTER ADDING OR REMOVING A SOURCE FILE

```sh
python3 tools/gen_fortran_deps.py . > makefiles/Makefile.deps
```

The Makefile expresses no Fortran module dependencies of its own, so without
this `make -j` races and fails with `Cannot open module file`. Regenerate
whenever the set of modules changes -- adding `turbogap_vdw.f90` and removing
the three `gpu_var_*` modules both required it.

Two traps in that file, both of which report success while being wrong:

* **The `include` must come LAST in the Makefile.** Placed before the first
  rule, a generated dependency line becomes make's default goal: `make` builds
  a single object and **exits 0**. Only counting objects catches it.
* **Duplicate module definitions must be gone first.** The generator fails loudly
  if two files define the same module name.

### Full rebuild

```sh
make clean       # objects, mods, binary; keeps the build dirs
make deepclean   # removes build/ bin/ include/ lib/ entirely
```

After `make clean` on the GPU tree, rebuild **serially once** if `-j` fails --
that is the symptom of a stale `Makefile.deps`, not a flaky build.

---

## 2. Test

```sh
# CPU: 20 cases, bit-exact against a frozen baseline binary, ~8 min
cd $M && tests/regression/run.sh
cd $M && tests/regression/run.sh vdw_ts co_predict     # named cases
cd $M && tests/regression/run.sh --list

# GPU: 4 cases, compared against the CPU build. REF_BIN must be ABSOLUTE --
# it is resolved from the staging directory, not from $PWD.
cd $G && export HOP_ROOT=/u/74/zarrout1/unix/work/hop \
      && export TURBOGAP_ARCH=Aalto_gfortran_openblas_hip_cuda \
      && TURBOGAP_REF_BIN=$M/bin/turbogap tests/gpu/run_regression.sh
```

Expected: CPU **20 passed, 0 failed**; GPU **3 passed, 0 failed, 1 xfail**
(`XRD_mad`, the known ~1e-5 `local_energy` drift).

Three suites sit outside `tests/regression`, because they ask whether an
answer is *right* rather than whether it is *unchanged*, and the frozen
baseline cannot run their inputs at all:

```sh
cd $M && tests/dipole/run.sh        # dipole prediction against QUIP
cd $M && tests/xrd_debye/run.sh     # Debye XRD against an independent sum,
                                    # and its forces against finite differences
cd $M && tests/xrd_lp_pdf/run.sh    # the Lorentz-polarization factor on the
                                    # pdf route, pattern and forces
```

### Do not let a rebuild race a running suite

`run.sh` reads `bin/turbogap` as it goes, so rebuilding mid-run gives a
meaningless mixed result. Stage the binary first:

```sh
cd $M && make -j12 && cp bin/turbogap /tmp/tg_under_test \
      && TURBOGAP_BIN=/tmp/tg_under_test tests/regression/run.sh
```

### Useful environment

| | |
|---|---|
| `TURBOGAP_BIN` | binary under test |
| `TURBOGAP_REF_BIN` | reference binary (**absolute path**) |
| `TURBOGAP_KEEP=1` | keep the staging directory for inspection |
| `TURBOGAP_BLESS=1` | regenerate a golden case's `expected/` -- deliberately, never to green a red suite |
| `TURBOGAP_TIME_TOL` | fail if test/ref wall-clock exceeds this ratio |

---

## 3. Checks that are not the suite

```sh
# device allocator/deallocator pairing (static; exits 1 on a mismatch)
cd $G && python3 tools/gpu_check_alloc_pairs.py
cd $G && python3 tools/gpu_check_alloc_pairs.py --quiet     # mismatches only

# compute-sanitizer, on the SMALL cell by default
cd $G && tools/gpu_memcheck.sh                       # memcheck
cd $G && tools/gpu_memcheck.sh --tool initcheck      # uninitialised device reads
cd $G && tools/gpu_memcheck.sh --leak --keep
cd $G && tools/gpu_memcheck.sh --case CO_predict     # slow: 10-100x

# formatting, over everything (rare; never on one branch alone)
pre-commit run --all-files
```

### Cross-branch divergence, the number the merge is measured by

```sh
cd /u/74/zarrout1/unix/work/cpu_vs_gpu_tests
for f in $(ls turbogap_master_2026/src/*.f90 | xargs -n1 basename); do
  m=turbogap_master_2026/src/$f; g=turbogap_gpu_commit_mahti/src/$f
  [ -f "$g" ] || continue
  d=$(diff -w --ignore-blank-lines "$m" "$g" | grep -c '^[<>]')
  [ "$d" -gt 0 ] && printf '%-22s %5s %5s  %6s\n' "$f" "$(wc -l < $m)" "$(wc -l < $g)" "$d"
done | sort -k4 -rn
```

Use `diff -w --ignore-blank-lines`. **`git diff -B` is `--break-rewrites`, not
"ignore blank lines"** -- it reported ~19,000 non-whitespace changes in a pure
reformat.

Files that are byte-identical across the trees (the merge scoreboard):

```sh
cd /u/74/zarrout1/unix/work/cpu_vs_gpu_tests
for f in $(ls turbogap_master_2026/src/*.f90 | xargs -n1 basename); do
  cmp -s turbogap_master_2026/src/$f turbogap_gpu_commit_mahti/src/$f && echo "  $f"
done
```

---

## 4. Committing

fprettify runs as a pre-commit hook and **rewrites files during `git commit`,
which then ABORTS**. That is the normal flow, not a failure:

```sh
git add -A && git commit -m "..."          # may abort, having reformatted
git diff -w --ignore-blank-lines           # MUST be empty: whitespace only
make -j12 DEBUG=0 && tests/.../run.sh      # rebuild and retest the new bytes
git add -A && git commit -m "..."          # now it lands
```

Commit message bodies with parentheses in shell-quoted heredocs break `ssh`
one-liners -- write the message to a file and use `git commit -F`.

---

## 5. Debugging

```sh
# gdb cannot attach on alt (ptrace_scope = 2). -fcheck=bounds is the tool.
# Patch the DEBUG=1 branch of makefiles/Makefile.Aalto_* to add it, then restore.
make DEBUG=1

# run one case by hand, keeping stderr
cd $(mktemp -d) && cp $G/../turbogap_tests/vdw_P/p4_dimer.xyz atoms.xyz \
  && ln -s $G/../turbogap_tests/vdw_P/gap_files . && printf '...' > input \
  && $G/bin/turbogap predict > out.log 2> err.log; echo "exit=$?"; cat err.log
```

**A bare `stop` in Fortran exits 0.** Three separate places reported failure and
exited successfully. Never trust an exit code alone; check for the output the
run was supposed to produce.

**The GPU binary is not bit-reproducible run to run** (three builds of an
unchanged tree gave three md5s, and one unchanged binary gives trajectories
differing at ~1e-12). So the bit-exact contract cannot be used there: do each
change on master first, where it holds, and port the verified result.

---

## `mem_fraction`: sizing `max_Gbytes_per_process` from the node

`max_Gbytes_per_process` is the divisor `get_number_of_atom_pairs` uses to split
the SOAP descriptor loop: it keeps splitting until an estimated batch fits
inside it. Its default is `1.0`, a number chosen with no machine in mind. On a
256 GB node that asks for roughly sixty times more batches than necessary, and
every extra batch repeats the whole per-batch setup — `assign_species_multiplicity`,
the neighbour-subset gather, the allocation, the scatter back — for the same
total work.

So it is now sized from the node instead:

```
mem_fraction = 0.25        ! the default; set it to 0 to turn sizing off
```

```
 Host memory budget:                    |
 node has       31.0 GB (MemTotal)      |
 budget =  0.250 x total / 1 rank(s)    |
 max_Gbytes_per_process =      7.761 GB |
```

`gpu_memory_budget_init` in `src/gpu_context.f90` does it, called from
`turbogap.f90` immediately after `read_input_and_gap_files`. The GPU branch
calls the same name in the same place, where it budgets from the device rather
than from the node — the driver line is identical on both branches and neither
driver learns which implementation it got.

Four decisions, each of which is somewhere this could be wrong in a way that
looks like a code bug:

* **An explicit `max_Gbytes_per_process` in the input wins.** Someone who set
  that keyword did so because the estimate is wrong for their system.
  `read_files` records that the input named it (`params%max_Gbytes_set`) and the
  sizing then leaves it alone and says so.

* **The cgroup limit is read before `/proc/meminfo`.** Under slurm, `--mem` puts
  the job in a cgroup whose limit is a fraction of the node, and `MemTotal`
  reports the whole node regardless. Budgeting from `MemTotal` on a shared node
  is a reliable way to be OOM-killed while the code believes it left three
  quarters of memory spare. The limit that applies is not the one at the root of
  the hierarchy — slurm nests each job step — so it reads `/proc/self/cgroup` and
  walks *up*, taking the first numeric limit, which is the tightest that applies.

* **It uses the total, not what is free right now.** Free memory moves between
  runs and between MD steps. See the next section for why that matters more than
  it sounds like it should.

* **The fraction is well under 1.** The SOAP batch is not the only thing
  resident. The global neighbour arrays (`rjs`, `thetas`, `phis`, `xyz` and the
  list itself, ~52 bytes per pair over *all* pairs, not just this batch), the
  positions, forces and velocities, and the descriptor outputs are all live at
  the same time, and none of them is counted by the estimate this keyword
  bounds.

### Why the regression suite pins it

The batch count changes the order the per-batch force contributions are summed,
so it moves the last bits of any force on an atom that neighbours a batch
boundary. Turning this on unpinned took **8 of the 21 cases red** — forces
differing at ~1e-8, virial in the eighth digit. Nothing was wrong with the
numbers; what was wrong is that the answers had become a function of how much
RAM the machine has, so a reference recorded on one node would not reproduce on
another.

`tests/regression/run.sh` therefore appends `max_Gbytes_per_process = 1.0` to
every staged input, for the same reason and in the same way it pins
`random_seed`: a comparison against stored text is only meaningful if everything
that moves the text is held still. Production runs get the automatic sizing; the
suite gets a number. With that in place the suite is back to **21 passed, 0
failed**.

This is worth remembering as a general property rather than a quirk of this
keyword: **anything that changes the batch decomposition changes the last digits.**
That includes `mem_fraction`, an explicit `max_Gbytes_per_process`, and the rank
count.

---

## The SOAP batch memory model on this branch, and what it costs the suite

`get_number_of_atom_pairs` sized a batch as `n_atom_pairs_in * n_max * k_max *
150`, one constant standing in for every array the descriptor path allocates.
Only the three `cnk` derivative arrays actually go as `n_max*k_max`; the rest go
as `n_soap`, as `n_max`, as `n_species`, or are flat — so no constant tracks them
across potentials, and it is blind to SOAP compression, which takes `n_soap` from
324 to 72 for `carbon.gap`. Nothing was counted **per site** either, though `cnk`
is `16*k_max*n_max` per site.

`soap_batch_memory_model` in `src/neighbors.f90` enumerates them instead, one
line per allocation, with `SOAP_BATCH_SAFETY = 1.10` for what it cannot see.

**These coefficients are this branch's and must not be copied from the GPU
branch.** This branch calls `get_derivatives`, so the coefficient and derivative
arrays are live host arrays here; on the GPU branch they are dummies and the real
allocations are on the device. The numbers differ (23,804 vs 25,760 bytes per
pair) and copying either way under-counts the other side's memory.

For `carbon.gap` this gives **23,804 bytes per pair against the 54,000 assumed** —
the loop was split about 2.3x more finely than the memory required.

Measured, 32,785 atoms, 2 MD steps, 6 ranks: **49.70 -> 47.43 s (-4.6%)**.
Smaller than the same change on the GPU branch (-14.9% at 1M), and the reason is
worth knowing: the CPU spends its time on the descriptor FLOPs themselves, so
per-batch overhead is a smaller share of the total. The gain grows with the batch
count, so it is larger at a million atoms and at low rank counts.

### It turns the suite red, and that is a policy question

`tests/regression/run.sh` compares output **text, bit-exactly**, against a frozen
pre-refactor binary (`baseline/turbogap.e6eb1aa`). Changing the estimate changes
the batch count, which regroups the per-batch partial sums, which moves the last
digits. Result: **13 passed, 8 failed** — all of them the SOAP-carrying cases,
all last-digit (`energy_soap` differs by 1.1e-5 on 65828, i.e. 1.7e-10 relative).

Pinning `max_Gbytes_per_process` in the staging step does NOT fix this, the way
it fixed the `mem_fraction` case: there the keyword moved, here the estimate that
consumes it moved, and the baseline binary carries the old one.

So there are three options and they are not equivalent:

1. **Regenerate the baseline.** Cheapest, but it discards the property the suite
   was built for — bit-exactness against the pre-refactor code — and that
   property has already caught real defects.
2. **Give the estimate a legacy switch** and have the suite set it. Keeps both
   properties; costs a keyword that exists only to keep a frozen comparison
   valid, and carries the old model forever.
3. **Compare numerically rather than textually** for these cases. The real fix:
   the suite is asserting text equality where it means numerical agreement. It is
   also the largest change, and it weakens the check for everything else unless
   the tolerance is chosen per quantity.

Not decided here. The change is **left uncommitted on this branch** pending that
decision, because committing something that turns the suite red is not a call to
make on someone else's test policy.

---

## GPU vs this branch, measured

Same input (`large_systems/diamond_1M/input.md`, `gpu_mem_fraction` stripped for
the CPU runs), 2 MD steps, everything else equal, each run alone on the box.

| atoms | GPU, 1 rank | CPU, 1 rank | CPU, 6 ranks | GPU vs 1 core | GPU vs node |
|---|---|---|---|---|---|
| 13,824 | 14.70 s | 68.65 s | 21.78 s | 4.67x | 1.48x |
| 32,785 | 32.31 s | 160.66 s | 48.28 s | 4.97x | 1.49x |
| 124,959 | 121.50 s | 608.35 s | 178.41 s | 5.01x | 1.47x |

The ratio is flat across a 9x range in system size: the GPU is **~5x one core**
and **~1.5x the whole six-core node**. The CPU's own parallel efficiency is 57%
at 6 ranks (608.35 / 178.41 = 3.41x on 6 cores).

**What this does and does not say.** `alt` pairs an RTX A2000 -- 26 SMs, 5.7 GB,
an entry-level workstation card -- with a 6-core desktop i5-11600. Both ends are
low-end, so 1.5x characterises *this box*. On Mahti or LUMI the CPU side is a
64-128 core node and the GPU an A100 or MI250X; that ratio is a different
measurement and cannot be inferred from this one. What does carry over is the
per-phase split and the scaling, which is what the optimisation work was aimed
at.

For scale, the GPU side moved during this work: 125k went 207.3 -> 121.5 s, and
the million-atom case 3583 -> 1770 s.
