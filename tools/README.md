# TurboGAP tools

Author(s): Miguel A. Caro, Tigany Zarrouk, Mikhail Kuklin, Richard Jana, Jan Kloppenburg

Pre- and post-processing aids for TurboGAP. Only tools the build, the tests or
the hooks depend on are kept here; the one-off analysis scripts live in
`archive/cpu_repo/tools/` beside the repo.

| tool | what it does |
|---|---|
| `gen_fortran_deps.py` | regenerates `makefiles/Makefile.deps`; run it whenever the set of modules changes |
| `keyword_docs.py` | builds `docs/keywords.{md,html}` and `src/keyword_help.f90` from the `!>` blocks in `read_files.f90` (`make docs`) |
| `setup_dev_env.sh` | installs the pinned fprettify/clang-format venv and the pre-commit hook |
| `quip_xml_to_gap/` | converts a QUIP XML potential into the `.gap` files TurboGAP reads |
| `shear_xyz.py` | applies a shear to an extended-XYZ cell, for the regression cases |
| `split_fortran_vars.sh` | reshapes a comma-batched declaration block into one variable per line |

See the `README.md` inside a subdirectory for that tool's own documentation.
