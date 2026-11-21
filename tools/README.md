# TurboGAP tools

Author(s): Miguel A. Caro, Tigany Zarrouk, Mikhail Kuklin, Richard Jana, Jan Kloppenburg

Only tools the build, the tests, the hooks or a developer check use are kept
here. `docs/DEVELOPMENT.md` says when each is run.

| tool | what it does |
|---|---|
| `gen_fortran_deps.py` | regenerates `makefiles/Makefile.deps` and `Makefile.deps.gpu`; run it whenever the set of modules changes |
| `keyword_docs.py` | builds `docs/keywords.{md,html}` and `src/keyword_help.f90` from the `!>` blocks in `read_files.f90` (`make docs`) |
| `setup_dev_env.sh` | installs the pinned fprettify/clang-format venv and the pre-commit hook |
| `check_reformat_only.py` | confirms a formatter changed only layout in CUDA/C++ files, by comparing tokens |
| `split_fortran_vars.sh` | reshapes a comma-batched declaration block into one variable per line |
| `install_kokkos.sh` | builds Kokkos for the `KOKKOS=1` device build |
| `verify_kokkos.sh` | runs the regression suite on the Kokkos build against the CUDA build, less the CUDA build's own noise |
| `gpu_check_alloc_pairs.py` | checks every device buffer is freed by the call matching its allocation |
| `gpu_memcheck.sh` | runs a case under compute-sanitizer |
| `xyz_frames.py` | reads extended-XYZ trajectories by field name; used by `tests/relaxation` |
| `quip_xml_to_gap/` | converts a QUIP XML potential into the `.gap` files TurboGAP reads |

See the `README.md` inside a subdirectory for that tool's own documentation.
