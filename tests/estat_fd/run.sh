#!/usr/bin/env bash
#
# Do the electrostatic forces and virial agree with a finite difference of the
# electrostatic energy, for each method and on more than one rank?
#
# The term is isolated by running every geometry with the method on and with
# estat_method = "none" and differencing (fd_gradient.py --family estat); the
# GAP that predicts the charges is the same in both and cancels.
#
# The system is a nine-atom C/Li cluster in a 30 A box with estat_rcut = 8 A,
# larger than any pair in it. "direct" has a sharp cutoff, so a displacement
# that carries a pair across it has no derivative; here none can.
#
# The mpi3 legs split nine atoms three ways, so every rank but the first holds
# a slice that starts past atom 1. That is where a local index read as a global
# one goes wrong, and one rank never shows it.
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
DATA=$DATA_ROOT/CLi
WORK=${TMPDIR:-/tmp}/turbogap_estat_fd.$$

FD=$here/../xrd_debye/fd_gradient.py

# The energy is read as the sum of per-atom energies written to 8 decimals, so
# over nine atoms a central difference at h = 1e-3 is good to ~1e-7 eV/A,
# three orders under the forces here.
H=1e-3
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

if [ ! -e "$DATA/gap_files/CCLi_estat_ljrep.gap" ]; then
  printf 'SKIP: missing test data %s\n' "$DATA/gap_files/CCLi_estat_ljrep.gap"
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

# leg <name> <estat_method> <ranks> [extra deck lines]
leg() {
  local name=$1 method=$2 ranks=$3 extra=${4:-}
  local dir=$WORK/$name
  printf '== %s ==\n' "$name"

  mkdir -p "$dir"
  ln -sf "$DATA/gap_files" "$dir/gap_files"
  cp "$here/cluster.xyz" "$dir/atoms.xyz"
  cat >"$dir/input" <<EOF
atoms_file = 'atoms.xyz'
pot_file = 'gap_files/CCLi_estat_ljrep.gap'
n_species = 2
species = C Li
masses = 12.01 6.94
e0 = 0. 0.
estat_method = "$method"
estat_rcut = 8.0
estat_dsf_alpha = 0.12
estat_damped = .true.
EOF
  [ -n "$extra" ] && printf '%s\n' "$extra" >>"$dir/input"

  if "$PYTHON" "$FD" "$dir" --bin "$BIN" --family estat --scale 1.0 --ranks "$ranks" \
       --h "$H" --strain "$STRAIN" --tol "$TOL" --atoms-to-check 3; then
    printf '    PASS\n'
    pass=$((pass + 1))
  else
    printf '    FAIL\n'
    fail=$((fail + 1))
  fi
}

GSF='estat_tsf = .true.
estat_sp = .true.
estat_gsf = .true.'
SELF='estat_self_energy_correction = .true.'

run_leg() {
  case $1 in
  gsf) leg gsf gsf 1 "$GSF" ;;
  gsf_mpi3) leg gsf_mpi3 gsf 3 "$GSF" ;;
  gsf_self) leg gsf_self gsf 1 "$GSF
$SELF" ;;
  dsf) leg dsf dsf 1 ;;
  dsf_mpi3) leg dsf_mpi3 dsf 3 ;;
  direct) leg direct direct 1 ;;
  direct_mpi3) leg direct_mpi3 direct 3 ;;
  *) die "no such leg: $1 (have gsf gsf_mpi3 gsf_self dsf dsf_mpi3 direct direct_mpi3)" ;;
  esac
}

legs=${*:-gsf gsf_mpi3 gsf_self dsf dsf_mpi3 direct direct_mpi3}
for l in $legs; do
  run_leg "$l"
done

printf '\npassed: %d   failed: %d\n' "$pass" "$fail"
[ "$fail" -eq 0 ]
