#!/usr/bin/env bash
#
# Do the XPS forces and virial agree with a finite difference of the XPS
# energy?
#
# The term is isolated by running every geometry at the deck's energy scale
# and at zero and differencing (fd_gradient.py --family xps): the GAP, and the
# core_electron_be model the spectrum is built from, are the same in both and
# cancel. The system is the 512-atom diamond cell of xps_mad_diamond, fitted
# against its sp3 spectrum with the default squared_diff similarity.
#
# Environment:
#   TURBOGAP_BIN        binary under test (default: <repo>/bin/turbogap)
#   TURBOGAP_DATA_ROOT  directory holding the test systems (see ../data_root.sh)
#   TURBOGAP_PYTHON     interpreter for the reference; needs numpy
#   TURBOGAP_KEEP       keep the staging directory for inspection
#
# Usage: run.sh [leg ...]      (no argument runs all of them)

set -u

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/../.." && pwd)

BIN=${TURBOGAP_BIN:-$repo/bin/turbogap}
PYTHON=${TURBOGAP_PYTHON:-python3}
# shellcheck source=../data_root.sh
. "$here/../data_root.sh"
DATA=$DATA_ROOT/xps_opt
WORK=${TMPDIR:-/tmp}/turbogap_xps_fd.$$

FD=$here/../xrd_debye/fd_gradient.py

# The spectrum is broadened by sigma = 0.2 eV, so the energy is strongly
# curved: the central-difference error falls 4.8e-4, 1.2e-4, 3.0e-5, 7.4e-6 for
# h = 2e-3, 1e-3, 5e-4, 2.5e-4 -- second order, so the gradient is right and h
# is what has to be small. The scale is raised from the deck's 2.48 so that the
# energy, printed to 8 decimals, still resolves that h.
H=2.5e-4
STRAIN=1e-4
TOL=1e-4

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 2
}

[ -x "$BIN" ] || die "binary under test not found or not executable: $BIN"
command -v "$PYTHON" >/dev/null || die "python interpreter not found: $PYTHON"
"$PYTHON" -c 'import numpy' 2>/dev/null ||
  die "$PYTHON has no numpy, which the reference needs (set TURBOGAP_PYTHON)"

for f in atoms.xyz gap_files/CO.gap xps_diamond_sp3_285.35.dat; do
  if [ ! -e "$DATA/$f" ]; then
    printf 'SKIP: missing test data %s\n' "$DATA/$f"
    exit 0
  fi
done

mkdir -p "$WORK"
if [ -z "${TURBOGAP_KEEP:-}" ]; then
  trap 'rm -rf "$WORK"' EXIT
else
  trap 'printf "staging kept in %s\n" "$WORK"' EXIT
fi

pass=0
fail=0

# leg <name> <ranks>
leg() {
  local name=$1 ranks=$2
  local dir=$WORK/$name
  printf '== %s ==\n' "$name"

  mkdir -p "$dir"
  ln -sf "$DATA/gap_files" "$dir/gap_files"
  ln -sf "$DATA/xps_diamond_sp3_285.35.dat" "$dir/exp.dat"
  cp "$DATA/atoms.xyz" "$dir/atoms.xyz"
  cat >"$dir/input" <<DECK
atoms_file = 'atoms.xyz'
pot_file = 'gap_files/CO.gap'
n_species = 2
species = C O
masses = 12.01 15.99
e0 = -.16138053 0.
n_exp = 1
exp_labels = 'xps'
exp_data_files = 'exp.dat'
exp_n_samples = 601
exp_energy_scales = 248.0
exp_energies = .true.
exp_forces = .true.
xps_sigma = 0.2
DECK

  if "$PYTHON" "$FD" "$dir" --bin "$BIN" --family xps --scale 248.0 --ranks "$ranks" \
       --h "$H" --strain "$STRAIN" --tol "$TOL" --atoms-to-check 3 --pick largest; then
    printf '    PASS\n'
    pass=$((pass + 1))
  else
    printf '    FAIL\n'
    fail=$((fail + 1))
  fi
}

run_leg() {
  case $1 in
  serial) leg serial 1 ;;
  mpi2) leg mpi2 2 ;;
  *) die "no such leg: $1 (have serial mpi2)" ;;
  esac
}

legs=${*:-serial mpi2}
for l in $legs; do
  run_leg "$l"
done

printf '\npassed: %d   failed: %d\n' "$pass" "$fail"
[ "$fail" -eq 0 ]
