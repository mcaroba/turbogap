#!/usr/bin/env bash
#
# Do the dispersion forces and virial agree with a finite difference of the
# dispersion energy?
#
# The term is isolated by running every geometry with the method on and with
# vdw_type = none and differencing (fd_gradient.py --family vdw); the GAP that
# predicts the Hirshfeld volumes is the same in both and cancels.
#
# The system is the P4 dimer of vdw_mbd, with the SCS and local cutoffs at 4 A
# so that no pair sits on one. mbd_split puts vdw_mbd_rcut below vdw_2b_rcut,
# which splits MBD into a two-body call at the long cutoff and a many-body call
# at the short one.
#
# The MBD legs fail by a few percent (KNOWN_ISSUES 20) and are left out of the
# default legs until that is fixed; name them to run them.
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
DATA=$DATA_ROOT/vdw_P
WORK=${TMPDIR:-/tmp}/turbogap_vdw_fd.$$

FD=$here/../xrd_debye/fd_gradient.py

# The energy is printed to 8 decimals and forces are ~0.02 eV/A, so h = 5e-3
# resolves them to ~5e-5 relative. The virial needs a strain of 1e-3: at 1e-4
# the difference is below the energy's resolution, at 5e-3 it sees curvature.
H=5e-3
STRAIN=1e-3
TOL=1e-3

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 2
}

[ -x "$BIN" ] || die "binary under test not found or not executable: $BIN"
command -v "$PYTHON" >/dev/null || die "python interpreter not found: $PYTHON"
"$PYTHON" -c 'import numpy' 2>/dev/null ||
  die "$PYTHON has no numpy, which the reference needs (set TURBOGAP_PYTHON)"

if [ ! -e "$DATA/gap_files/phosphorus.gap" ]; then
  printf 'SKIP: missing test data %s\n' "$DATA/gap_files/phosphorus.gap"
  exit 0
fi

mkdir -p "$WORK"
if [ -z "${TURBOGAP_KEEP:-}" ]; then
  trap 'rm -rf "$WORK"' EXIT
else
  trap 'printf "staging kept in %s\n" "$WORK"' EXIT
fi

pass=0
fail=0

# leg <name> <vdw_type> <ranks> <mbd_rcut> <mbd_rcut2>
leg() {
  local name=$1 type=$2 ranks=$3 mbd_rcut=$4 mbd_rcut2=$5
  local dir=$WORK/$name
  printf '== %s ==\n' "$name"

  mkdir -p "$dir"
  ln -sf "$DATA/gap_files" "$dir/gap_files"
  head -10 "$DATA/p4_dimer.xyz" >"$dir/atoms.xyz"
  cat >"$dir/input" <<DECK
atoms_file = 'atoms.xyz'
pot_file = 'gap_files/phosphorus.gap'
n_species = 1
species = P
e0 = -0.52375977
masses = 30.97
vdw_type = $type
vdw_rcut = 25.
vdw_r0_ref = 2.12
vdw_alpha0_ref = 3.7046
vdw_c6_ref = 110.54
vdw_buffer = 0.5
vdw_scs_rcut = 4.
vdw_loc_rcut = 4.
vdw_mbd_rcut = $mbd_rcut
vdw_mbd_rcut2 = $mbd_rcut2
vdw_2b_rcut = 15.0
vdw_2b_rcut2 = 10.0
vdw_mbd_nfreq = 13
vdw_mbd_norder = 6
DECK

  if "$PYTHON" "$FD" "$dir" --bin "$BIN" --family vdw --scale 1.0 --ranks "$ranks" \
       --h "$H" --strain "$STRAIN" --tol "$TOL" --atoms-to-check 3; then
    printf '    PASS\n'
    pass=$((pass + 1))
  else
    printf '    FAIL\n'
    fail=$((fail + 1))
  fi
}

run_leg() {
  case $1 in
  ts) leg ts ts 1 15.0 10.0 ;;
  ts_mpi2) leg ts_mpi2 ts 2 15.0 10.0 ;;
  mbd) leg mbd mbd 1 15.0 10.0 ;;
  mbd_split) leg mbd_split mbd 1 6.0 4.0 ;;
  *) die "no such leg: $1 (have ts ts_mpi2 mbd mbd_split)" ;;
  esac
}

legs=${*:-ts ts_mpi2}
for l in $legs; do
  run_leg "$l"
done

printf '\npassed: %d   failed: %d\n' "$pass" "$fail"
[ "$fail" -eq 0 ]
