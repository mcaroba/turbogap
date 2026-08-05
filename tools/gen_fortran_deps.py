#!/usr/bin/env python3
"""Generate Fortran module dependencies for the TurboGAP Makefile.

The Makefile lists which sources are compiled but expresses no dependency between
an object and the modules it USEs, so `make -j` can compile a consumer before the
producer has written its .mod file. This reads the Makefile's own SRC_* lists and
pattern rules, maps each `module <name>` to the object defining it, and emits the
ordering rules.

    python3 tools/gen_fortran_deps.py . > makefiles/Makefile.deps

Works on both the CPU and GPU trees: the source lists, the directories searched,
and which object groups are actually linked are all read from the Makefile rather
than hardcoded, so this stays correct as those lists change.

Exits non-zero if two files define the same module, which would make the mapping
ambiguous.
"""
import re
import sys
from pathlib import Path

# Modules provided by the compiler / MPI / OpenMP rather than by this tree.
EXTERNAL = {
    'iso_c_binding', 'iso_fortran_env', 'ieee_arithmetic', 'ieee_exceptions',
    'ieee_features', 'mpi', 'mpi_f08', 'omp_lib', 'omp_lib_kinds', 'cudafor',
    'openacc',
}

# Groups that hold Fortran sources. SRC_CUDA / SRC_CC are .cu/.cc and irrelevant.
GROUPS = ['SRC_BASE', 'SRC_TP_BT', 'SRC_ST', 'SRC_STOP', 'SRC']

RE_MODULE = re.compile(r'^\s*module\s+([a-z_][a-z0-9_]*)\s*$', re.I)
RE_MODULE_PROC = re.compile(r'^\s*module\s+(procedure|subroutine|function)\b', re.I)
RE_USE = re.compile(r'^\s*use\s*(?:,\s*intrinsic\s*)?(?:::)?\s*([a-z_][a-z0-9_]*)', re.I)
RE_OBJ_RULE = re.compile(r'^\$\(BUILD_DIR\)/%\.o:\s*(\S+)/%\.f90')


def read_makefile(path):
    """Join continuation lines and strip comments."""
    lines = []
    buf = ''
    for raw in path.read_text(encoding='utf-8').splitlines():
        buf += raw
        if buf.rstrip().endswith('\\'):
            buf = buf.rstrip()[:-1] + ' '
            continue
        lines.append(buf)
        buf = ''
    if buf:
        lines.append(buf)
    return lines


def parse_var(lines, name):
    """Return the .f90 entries assigned to a make variable, comments stripped."""
    pat = re.compile(r'^\s*' + re.escape(name) + r'\s*:?=\s*(.*)$')
    for line in lines:
        m = pat.match(line)
        if m:
            body = m.group(1).split('#', 1)[0]
            return [t for t in body.split() if t.endswith('.f90')]
    return []


def main():
    repo = Path(sys.argv[1] if len(sys.argv) > 1 else '.').resolve()
    mk = repo / 'Makefile'
    if not mk.is_file():
        print(f'ERROR: no Makefile at {mk}', file=sys.stderr)
        return 1
    lines = read_makefile(mk)

    search_dirs = [m.group(1) for line in lines for m in [RE_OBJ_RULE.match(line)] if m]
    if not search_dirs:
        print('ERROR: found no $(BUILD_DIR)/%.o: <dir>/%.f90 pattern rules', file=sys.stderr)
        return 1

    # Only generate for groups whose objects are actually a prerequisite of
    # something. The GPU tree defines SRC_STOP but never links OBJ_STOP.
    rule_text = '\n'.join(l for l in lines if l.startswith(('libturbogap:', '$(BIN_DIR)/%:')))
    groups, skipped = [], []
    for g in GROUPS:
        files = parse_var(lines, g)
        if not files:
            continue
        objvar = '$(' + g.replace('SRC', 'OBJ', 1) + ')'
        (groups if objvar in rule_text else skipped).append((g, files, objvar))

    sources = [f for _, files, _ in groups for f in files]
    if not sources:
        print('ERROR: no source files found in the Makefile SRC_* lists', file=sys.stderr)
        return 1

    def locate(fname):
        for d in search_dirs:
            p = repo / d / fname
            if p.is_file():
                return p
        return None

    obj_of_module, uses_of_obj, dupes = {}, {}, []
    for fname in sources:
        path = locate(fname)
        if path is None:
            print(f'ERROR: cannot locate {fname} under {search_dirs}', file=sys.stderr)
            return 1
        stem = Path(fname).stem
        defines, uses = set(), set()
        for raw in path.read_text(encoding='utf-8', errors='replace').splitlines():
            line = raw.split('!', 1)[0]
            if not line.strip() or RE_MODULE_PROC.match(line):
                continue
            m = RE_MODULE.match(line)
            if m:
                defines.add(m.group(1).lower())
                continue
            u = RE_USE.match(line)
            if u and u.group(1).lower() not in EXTERNAL:
                uses.add(u.group(1).lower())
        uses_of_obj[stem] = uses
        for mod in defines:
            if mod in obj_of_module and obj_of_module[mod] != stem:
                dupes.append((mod, obj_of_module[mod], stem))
            obj_of_module[mod] = stem

    if dupes:
        for mod, a, b in dupes:
            print(f'ERROR: module {mod} defined by both {a}.f90 and {b}.f90 -- '
                  f'the mapping is ambiguous; remove or rename one', file=sys.stderr)
        return 1

    shadow = sorted(set(obj_of_module) & EXTERNAL)
    if shadow:
        print(f'ERROR: in-tree modules shadow compiler-provided names: {shadow}',
              file=sys.stderr)
        return 1

    unresolved, out = set(), []
    for fname in sources:
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
            out.append(f'$(BUILD_DIR)/{stem}.o: {rhs}')

    print('# ---- Fortran module dependencies (generated, do not edit) ----')
    print('# python3 tools/gen_fortran_deps.py . > makefiles/Makefile.deps')
    print('#')
    print('# Compiling a file that USEs a module needs that module\'s .mod to exist')
    print('# already. Without these rules `make -j` races and fails with')
    print('# "Cannot open module file".')
    print(f'# Covers: {", ".join(g for g, _, _ in groups)}')
    if skipped:
        print(f'# Not covered (objects are in no rule\'s prerequisites, so never '
              f'built): {", ".join(g for g, _, _ in skipped)}')
    print()
    for line in out:
        print(line)

    print(f'generated {len(out)} rules over {len(sources)} objects '
          f'from {mk}', file=sys.stderr)
    print(f'  groups covered: {[g for g, _, _ in groups]}', file=sys.stderr)
    if skipped:
        print(f'  groups skipped (never linked): {[g for g, _, _ in skipped]}',
              file=sys.stderr)
    if unresolved:
        print(f'  NOTE: USEd but defined nowhere in the build set, assumed '
              f'external: {sorted(unresolved)}', file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
