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
# CPU: 18 cases, bit-exact against a frozen baseline binary, ~7 min
cd $M && tests/regression/run.sh
cd $M && tests/regression/run.sh vdw_ts co_predict     # named cases
cd $M && tests/regression/run.sh --list

# GPU: 4 cases, compared against the CPU build. REF_BIN must be ABSOLUTE --
# it is resolved from the staging directory, not from $PWD.
cd $G && export HOP_ROOT=/u/74/zarrout1/unix/work/hop \
      && export TURBOGAP_ARCH=Aalto_gfortran_openblas_hip_cuda \
      && TURBOGAP_REF_BIN=$M/bin/turbogap tests/gpu/run_regression.sh
```

Expected: CPU **18 passed, 0 failed**; GPU **3 passed, 0 failed, 1 xfail**
(`XRD_mad`, the known ~1e-5 `local_energy` drift).

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

## 2b. Profile

`docs/PROFILING.md` is the whole workflow. The two commands:

```sh
tools/setup_profiling_env.sh --check    # what would stop a profile, and the fix
tools/profile_gpu.sh CO_predict         # build PROFILE=1, run under nsys, analyse
```

The cases come from `tests/gpu/cases.sh`, which `tests/gpu/run_regression.sh`
reads too, so a profile is always of an input the suite checks. **Never profile
a `DEBUG=1` build** -- the Makefile refuses it, because `-G` is a 2.1x slowdown
that also reorders the kernel ranking.

`make PROFILE=1` and `make OPENMP=1` each get their own object tree
(`build-profile/`, `build-omp/`), so neither can silently reuse the other's
objects or invalidate `bin/turbogap` while the suite is reading it.

## 2c. Device memory, and systems too big for the suite

Every run now ends with a `GPUmem` block: what the device has, what this rank
held, and its peak. An allocation that does not fit retries once (after draining
the stream, because `hipFreeAsync` is stream-ordered) and then reports the size,
the budget and the keyword that fixes it, instead of a bare "out of memory".

```
gpu_mem_fraction = 0.8      # size max_Gbytes_per_process from the card
```

Default 0 = off, so existing inputs are unchanged. The 1.0 GB default for
`max_Gbytes_per_process` was chosen with no device in mind. An explicit
`gpu_n_batches` is a **floor**: the automatic sizing raises it, never lowers it.

The diamond ladder (13,824 → 1,000,000 atoms) lives outside both repos in
`../large_systems/diamond_1M/`; see its `README.md` and `docs/PROFILING.md` §6.

```sh
tools/run_scaling_ladder.sh              # wall time, host RSS, device peak per size
tools/profile_gpu.sh diamond_125k        # one rung, under nsys
```

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
