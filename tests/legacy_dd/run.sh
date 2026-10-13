#!/usr/bin/env bash
#
# Does dd_type = legacy give the same energies, forces, virial and local
# properties as the replicated run of the same input?
#
# Each leg runs one input twice, replicated and decomposed, and compares
# trajectory_out.xyz number by number. The domains sum each atom's terms in a
# different order, so the comparison is a tolerance, not a diff.
#
# The legs are chosen for the ways a decomposition goes wrong:
#   - 2 domains along an axis: both neighbours are the same domain, so a halo
#     atom arrives twice and must be kept once;
#   - 3 along an axis: the two neighbours differ;
#   - 2x2x2: edge and corner atoms reach a domain through two forwards, and
#     their forces have to come back the same way;
#   - several ranks per domain, block and cyclic, so the halo forces of
#     non-root ranks have to reach the root;
#   - a sheared cell, where a face's normal is not a cell vector;
#   - five MD steps, so rank 0 integrates gathered forces and redistributes.
#
# Environment:
#   TURBOGAP_BIN        binary under test (default: <repo>/bin/turbogap)
#   TURBOGAP_DATA_ROOT  directory holding the test systems (see ../data_root.sh)
#   TURBOGAP_PYTHON     interpreter for the comparison
#   TURBOGAP_KEEP       keep the staging directory for inspection
#   TURBOGAP_RTOL, TURBOGAP_ATOL   comparison tolerances (default 1e-7 each)
#
# Usage: run.sh [leg ...]      (no argument runs all of them)

set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/../.." && pwd)

BIN=${TURBOGAP_BIN:-$repo/bin/turbogap}
PYTHON=${TURBOGAP_PYTHON:-python3}
# shellcheck source=../data_root.sh
. "$here/../data_root.sh"
# shellcheck source=../mpi_oversubscribe.sh
. "$here/../mpi_oversubscribe.sh"
DATA=$DATA_ROOT/CO
WORK=${TMPDIR:-/tmp}/turbogap_legacy_dd.$$
COMPARE=$repo/tests/regression/compare_tol.py

# Per-atom values are written to 8-10 significant digits; a reordered sum moves
# the last of them. A device binary does not reproduce itself even replicated,
# so it needs its own, looser, tolerance.
RTOL=${TURBOGAP_RTOL:-1e-7}
ATOL=${TURBOGAP_ATOL:-1e-7}

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 2
}

check_setup() {
  [ -x "$BIN" ] || die "binary under test not found or not executable: $BIN"
  command -v mpirun >/dev/null || die "mpirun not found; dd_type legacy needs an MPI build"
  command -v "$PYTHON" >/dev/null || die "python interpreter not found: $PYTHON"
  if [ ! -e "$DATA/gap_files/CO.gap" ]; then
    printf 'SKIP: missing test data %s\n' "$DATA/gap_files/CO.gap"
    exit 0
  fi
}

# The 897-atom cell with c tilted by 0.3 a, atoms carried along with the cell.
make_sheared() {
  "$PYTHON" - "$DATA/atoms_897.xyz" "$1" <<'EOF'
import re
import sys

src, dst = sys.argv[1], sys.argv[2]
lines = open(src).read().splitlines()
n = int(lines[0])
lat = [float(x) for x in re.search(r'Lattice="([^"]*)"', lines[1]).group(1).split()]
a, b, c = lat[0:3], lat[3:6], lat[6:9]
c2 = [c[k] + 0.3 * a[k] for k in range(3)]
out = [lines[0], lines[1].replace(re.search(r'Lattice="([^"]*)"', lines[1]).group(1),
                                   " ".join(f"{x:.10f}" for x in a + b + c2))]
for line in lines[2:2 + n]:
    f = line.split()
    r = [float(x) for x in f[1:4]]
    frac_c = r[2] / c[2]
    r2 = [r[k] + 0.3 * a[k] * frac_c for k in range(3)]
    out.append(" ".join([f[0]] + [f"{x:.10f}" for x in r2] + f[4:]))
open(dst, "w").write("\n".join(out) + "\n")
EOF
}

# write_input <dir> <md steps, 0 for a single point> [extra lines]
write_input() {
  local dir=$1 steps=$2 extra=${3:-}
  cat >"$dir/input" <<EOF
atoms_file = "atoms.xyz"
pot_file = "gap_files/CO.gap"
n_species = 2
species = C O
masses = 12.01 15.99
e0 = -.16138053 0.
EOF
  if [ "$steps" -gt 0 ]; then
    printf 'md_nsteps = %s\nthermostat = berendsen\nwrite_xyz = 1\nrandom_seed = 12345\n' "$steps" >>"$dir/input"
  fi
  if [ -n "$extra" ]; then
    printf '%s\n' "$extra" >>"$dir/input"
  fi
}

# run_one <dir> <mode> <ranks>
run_one() {
  local dir=$1 mode=$2 ranks=$3
  (cd "$dir" && OMP_NUM_THREADS=1 mpirun $MPI_OVERSUBSCRIBE -np "$ranks" "$BIN" "$mode" >out.log 2>&1)
}

# leg <name> <atoms file> <md steps> <ranks> <dd_grid> [affinity]
leg() {
  local name=$1 atoms=$2 steps=$3 ranks=$4 grid=$5 affinity=${6:-block}
  local mode=predict
  [ "$steps" -gt 0 ] && mode=md
  printf '== %s: %s ranks, grid %s, %s ==\n' "$name" "$ranks" "$grid" "$affinity"

  local ok=1 kind
  for kind in ref dd; do
    mkdir -p "$WORK/$name/$kind"
    ln -sf "$DATA/gap_files" "$WORK/$name/$kind/gap_files"
    cp "$atoms" "$WORK/$name/$kind/atoms.xyz"
  done
  write_input "$WORK/$name/ref" "$steps"
  write_input "$WORK/$name/dd" "$steps" "dd_type = legacy
dd_grid = $grid
dd_grid_affinity = $affinity"

  run_one "$WORK/$name/ref" "$mode" "$ranks" || ok=0
  run_one "$WORK/$name/dd" "$mode" "$ranks" || ok=0
  if [ "$ok" -eq 1 ] && grep -q 'dd_type' "$WORK/$name/dd/out.log" &&
     "$PYTHON" "$COMPARE" "$WORK/$name/ref/trajectory_out.xyz" "$WORK/$name/dd/trajectory_out.xyz" \
       --rtol "$RTOL" --atol "$ATOL"; then
    printf '    PASS\n'
    pass=$((pass + 1))
  else
    printf '    FAIL (logs in %s)\n' "$WORK/$name"
    tail -5 "$WORK/$name/dd/out.log" | sed 's/^/      /'
    fail=$((fail + 1))
  fi
}

# refused <name> <ranks> <extra line>: the run must stop with a dd_ error.
refused() {
  local name=$1 ranks=$2 extra=$3
  printf '== %s: refused ==\n' "$name"
  mkdir -p "$WORK/$name"
  ln -sf "$DATA/gap_files" "$WORK/$name/gap_files"
  cp "$DATA/atoms_897.xyz" "$WORK/$name/atoms.xyz"
  write_input "$WORK/$name" 0 "dd_type = legacy
$extra"
  if ! run_one "$WORK/$name" predict "$ranks" && grep -q 'ERROR: dd' "$WORK/$name/out.log"; then
    printf '    PASS\n'
    pass=$((pass + 1))
  else
    printf '    FAIL: ran, or stopped for another reason (log in %s)\n' "$WORK/$name"
    fail=$((fail + 1))
  fi
}

run_leg() {
  local small=$DATA/atoms_897.xyz
  case $1 in
  x2) leg x2 "$small" 0 2 "2 1 1" ;;
  x3) leg x3 "$small" 0 3 "3 1 1" ;;
  xy4) leg xy4 "$small" 0 4 "2 2 1" ;;
  xyz8) leg xyz8 "$small" 0 8 "2 2 2" ;;
  one3) leg one3 "$small" 0 3 "1 1 1" ;;
  block4) leg block4 "$small" 0 4 "2 1 1" block ;;
  cyclic6) leg cyclic6 "$small" 0 6 "1 3 1" cyclic ;;
  sheared) make_sheared "$WORK/sheared.xyz" && leg sheared "$WORK/sheared.xyz" 0 4 "2 1 2" ;;
  md8) leg md8 "$DATA/atoms_7176_vel.xyz" 5 8 "2 2 2" ;;
  narrow) refused narrow 4 "dd_grid = 4 1 1" ;;
  vdw) refused vdw 2 "dd_grid = 2 1 1
vdw_type = ts" ;;
  *) die "no such leg: $1" ;;
  esac
}

check_setup
mkdir -p "$WORK"
if [ -z "${TURBOGAP_KEEP:-}" ]; then
  trap 'rm -rf "$WORK"' EXIT
else
  trap 'printf "staging kept in %s\n" "$WORK"' EXIT
fi

pass=0
fail=0
legs=${*:-x2 x3 xy4 xyz8 one3 block4 cyclic6 sheared md8 narrow vdw}
for l in $legs; do
  run_leg "$l"
done

printf '\npassed: %d   failed: %d\n' "$pass" "$fail"
[ "$fail" -eq 0 ]
