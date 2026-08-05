#!/usr/bin/env python3
"""Generate Fortran module dependencies for the TurboGAP Makefile.

Scans the source files the Makefile actually compiles, maps `module <name>` to the
object that defines it, and emits `build/x.o: build/y.o` ordering rules so that a
parallel build cannot compile a consumer before its .mod file exists.

Prints the dependency block on stdout. Diagnostics go to stderr.
"""
import re
import sys
from pathlib import Path

REPO = Path(sys.argv[1]).resolve()

# Mirrors the SRC_* lists in the Makefile. Each entry: (filename, directory).
# Search order matters only for locating the file on disk.
SEARCH_DIRS = [
    'src',
    'src/third_party/bussi_thermostat',
    'src/soap_turbo/src',
    'src/stopping',
]

SRC_BASE = ['kinds.f90']
SRC_TP_BT = ['resamplekin.f90', 'fortran_cuda_interfaces.f90']
SRC_ST = ['soap_turbo_functions.f90', 'soap_turbo_radial.f90', 'soap_turbo_angular.f90',
          'soap_turbo.f90', 'soap_turbo_compress.f90']
SRC = ['timing.f90', 'splines.f90', 'types.f90', 'neighbors.f90', 'gap.f90', 'vdw.f90',
       'local_properties.f90', 'exp_utils.f90', 'xyz.f90', 'md.f90', 'mc.f90',
       'read_files.f90', 'gap_backend_gpu.f90', 'gap_interface.f90', 'mpi.f90',
       'exp_interface.f90', 'turbogap_exp.f90', 'turbogap_setup.f90', 'electrostatics.f90']
# Compiled straight into the binary by the $(BIN_DIR)/% rule, not into an object.
MAIN = ['turbogap.f90']

ALL_OBJ_SRC = SRC_BASE + SRC_TP_BT + SRC_ST + SRC

# Modules provided by the compiler / MPI / OpenMP, not by this tree.
EXTERNAL = {
    'iso_c_binding', 'iso_fortran_env', 'ieee_arithmetic', 'ieee_exceptions',
    'ieee_features', 'mpi', 'mpi_f08', 'omp_lib', 'omp_lib_kinds', 'cudafor',
    'openacc',
}

RE_MODULE = re.compile(r'^\s*module\s+([a-z_][a-z0-9_]*)\s*$', re.I)
RE_MODULE_PROC = re.compile(r'^\s*module\s+(procedure|subroutine|function)\b', re.I)
RE_USE = re.compile(r'^\s*use\s*(?:,\s*intrinsic\s*)?(?:::)?\s*([a-z_][a-z0-9_]*)', re.I)


def locate(fname):
    for d in SEARCH_DIRS:
        p = REPO / d / fname
        if p.is_file():
            return p
    return None


def scan(path):
    """Return (modules_defined, modules_used) for a Fortran source file."""
    defines, uses = set(), set()
    text = path.read_text(encoding='utf-8', errors='replace')
    for raw in text.splitlines():
        line = raw.split('!', 1)[0]  # strip comments; no '!' inside strings in practice here
        if not line.strip():
            continue
        if RE_MODULE_PROC.match(line):
            continue
        m = RE_MODULE.match(line)
        if m:
            defines.add(m.group(1).lower())
            continue
        u = RE_USE.match(line)
        if u:
            name = u.group(1).lower()
            if name not in EXTERNAL:
                uses.add(name)
    return defines, uses


def main():
    obj_of_module = {}      # module name -> object stem
    uses_of_obj = {}        # object stem -> set of module names
    dupes = []

    for fname in ALL_OBJ_SRC:
        path = locate(fname)
        if path is None:
            print(f'ERROR: cannot locate {fname}', file=sys.stderr)
            return 1
        stem = Path(fname).stem
        defines, uses = scan(path)
        uses_of_obj[stem] = uses
        for mod in defines:
            if mod in obj_of_module and obj_of_module[mod] != stem:
                dupes.append((mod, obj_of_module[mod], stem))
            obj_of_module[mod] = stem

    if dupes:
        for mod, a, b in dupes:
            print(f'ERROR: module {mod} defined by both {a}.f90 and {b}.f90', file=sys.stderr)
        return 1

    unresolved = set()
    lines = []
    for fname in ALL_OBJ_SRC:
        stem = Path(fname).stem
        deps = set()
        for mod in sorted(uses_of_obj[stem]):
            owner = obj_of_module.get(mod)
            if owner is None:
                unresolved.add(mod)
            elif owner != stem:
                deps.add(owner)
        if deps:
            rhs = ' '.join(f'$(BUILD_DIR)/{d}.o' for d in sorted(deps))
            lines.append(f'$(BUILD_DIR)/{stem}.o: {rhs}')

    # turbogap.f90 is compiled by the $(BIN_DIR)/% rule, which already lists every
    # object group as a prerequisite, so it needs no extra ordering. Report its
    # uses only as a cross-check.
    main_path = locate(MAIN[0])
    _, main_uses = scan(main_path)
    main_unresolved = {m for m in main_uses if m not in obj_of_module}

    print('# ---- Fortran module dependencies (generated) ----')
    print('# Regenerate with: tools/gen_fortran_deps.py . > makefiles/Makefile.deps')
    print('# Without these, a parallel build (make -j) can compile a consumer before')
    print('# the .mod file it USEs has been written, which fails intermittently.')
    print()
    for line in lines:
        print(line)

    print(f'\ngenerated {len(lines)} dependency rules '
          f'over {len(ALL_OBJ_SRC)} objects', file=sys.stderr)
    if unresolved:
        print(f'NOTE: modules USEd but not defined in the build set (assumed external): '
              f'{sorted(unresolved)}', file=sys.stderr)
    if main_unresolved:
        print(f'NOTE: turbogap.f90 USEs undefined-in-build-set: {sorted(main_unresolved)}',
              file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
