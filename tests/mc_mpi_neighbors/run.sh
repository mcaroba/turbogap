#!/usr/bin/env bash
#
# Are the neighbour lists a Monte Carlo run carries still right after the walk
# has changed the number of atoms -- and are they right on every rank?
#
# An MC step can insert an atom, remove one, swap two species or displace one,
# and each of those invalidates the neighbour lists. `turbogap mc` rebuilds them
# every step and, under MPI, each rank rebuilds its own slice from a broadcast
# copy of the new positions. There is a lot to get wrong: an atom count that is
# stale on one rank, a rebuild skipped, a slice boundary that moved when
# n_sites did.
#
# None of it shows up in a diff-based regression case, because a wrong
# neighbour list is perfectly reproducible. What it does show up in is the
# energy: an atom whose neighbours are missing has the wrong force. So the test
# is a round trip.
#
#   1. Run `turbogap mc` on N ranks with mc_write_xyz, so mc_all.xyz records
#      every accepted configuration together with the forces and local energies
#      the MC path computed for it.
#   2. Feed mc_all.xyz back to `turbogap predict` on one rank. Prediction reads
#      concatenated frames, and it builds every neighbour list from scratch for
#      each of them, which is the reference.
#   3. Compare per-atom forces and local energies, frame by frame.
#
# Agreement to round-off means the lists the walk carried were the lists a cold
# start would have built. Disagreement localises to the atom.
#
# The legs are rank counts. Leg r1 is the control: it shares the whole
# neighbour path with the MPI legs but none of the communication, so an r1
# failure is a rebuild bug and an r1 pass with an r2 failure is a
# communication bug.
#
# Environment:
#   TURBOGAP_BIN        binary under test (default: <repo>/bin/turbogap)
#   TURBOGAP_DATA_ROOT  directory holding the test systems. Default: the
#                       turbogap_tests clone beside this repository, which
#                       tests/fetch_test_data.sh creates on the first run
#                       (override where it goes with TURBOGAP_TESTS_DIR).
#   TURBOGAP_PYTHON     interpreter for the comparison (default: python3).
#                       Standard library only.
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
WORK=${TMPDIR:-/tmp}/turbogap_mc_mpi_neighbors.$$

# One neighbour missing from one atom moves that atom's force by O(0.1) eV/A,
# five orders of magnitude above this. The tolerance is here only to absorb the
# reduction order changing with the rank count, and the F16.8 the coordinates
# make the round trip in.
FORCE_TOL=1e-5
ENERGY_TOL=1e-6

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 2
}

[ -x "$BIN" ] || die "binary under test not found or not executable: $BIN"
command -v "$PYTHON" >/dev/null || die "python interpreter not found: $PYTHON"

for f in atoms.xyz gap_files/CO.gap co_molecule.xyz; do
  if [ ! -e "$DATA/$f" ]; then
    printf 'SKIP: missing test data %s\n' "$DATA/$f"
    printf '      (run tests/fetch_test_data.sh, or see tests/mc_mpi_neighbors/README.md)\n'
    exit 0
  fi
done

mkdir -p "$WORK"
if [ -z "${TURBOGAP_KEEP:-}" ]; then
  trap 'rm -rf "$WORK"' EXIT
else
  trap 'printf "staging kept in %s\n" "$WORK"' EXIT
fi

# The walk. Every move type that changes the atom list is in it, and
# mc_write_xyz puts every step in mc_all.xyz rather than every write_xyz-th.
# mc_min_dist is loose enough that insertions land close to something, which is
# the case where a stale neighbour list would matter.
deck() {
  cat <<EOF
atoms_file = 'atoms.xyz'
pot_file = 'gap_files/CO.gap'
n_species = 2
species = C O
masses = 12.01 15.99
e0 = -.16138053 0.
random_seed = 12345
max_Gbytes_per_process = 1.0

mc_nsteps = 14
n_mc_types = 4
mc_types = 'insertion' 'removal' 'move' 'swap'
mc_acceptance = 1 1 1 1
mc_move_max = 0.4

n_mc_mu = 1
mc_mu = 0.0
mc_species = 'O'
mc_min_dist = 1.0
mc_max_dist = 2.5

n_mc_swaps = 1
mc_swaps = 'C' 'O'

mc_write_xyz = .true.
EOF
}

# The same walk exchanging whole CO molecules instead of single atoms. Two
# atoms arrive or leave at once, so every array the walk carries changes length
# by two, and the removal takes atoms that need not be adjacent once earlier
# removals have compacted the list -- both places a neighbour list could be
# left describing the wrong number of atoms.
#
# No swap: swapping the species of an atom that belongs to a molecule would
# turn that molecule into something that is no longer the thing mc_mu prices.
mol_deck() {
  cat <<EOF
atoms_file = 'atoms.xyz'
pot_file = 'gap_files/CO.gap'
n_species = 2
species = C O
masses = 12.01 15.99
e0 = -.16138053 0.
random_seed = 12345
max_Gbytes_per_process = 1.0

mc_nsteps = 14
n_mc_types = 3
mc_types = 'insertion' 'removal' 'move'
mc_acceptance = 2 1 1
mc_move_max = 0.4

n_mc_mu = 1
mc_mu = -1.0
mc_species = 'CO'
mc_molecule_files = 'co_molecule.xyz'
mc_mu_reference = 'e0'
mc_min_dist = 1.2
mc_max_dist = 3.0

mc_write_xyz = .true.
EOF
}

# The prediction. Same potential, and nothing that would make it write anything
# but forces and local energies.
predict_deck() {
  cat <<EOF
atoms_file = 'mc_all.xyz'
pot_file = 'gap_files/CO.gap'
n_species = 2
species = C O
masses = 12.01 15.99
e0 = -.16138053 0.
max_Gbytes_per_process = 1.0
EOF
}

pass=0
fail=0

leg() {
  local ranks=$1
  local dir=$WORK/$3
  printf '== %s ==\n' "$3"

  if [ "$ranks" -gt 1 ] && ! command -v mpirun >/dev/null; then
    printf '    SKIP: mpirun not available\n'
    return 0
  fi

  mkdir -p "$dir"
  ln -sf "$DATA/gap_files" "$dir/gap_files"
  ln -sf "$DATA/atoms.xyz" "$dir/atoms.xyz"
  ln -sf "$DATA/co_molecule.xyz" "$dir/co_molecule.xyz"
  "$2" >"$dir/input"

  if [ "$ranks" -gt 1 ]; then
    (cd "$dir" && mpirun -np "$ranks" "$BIN" mc) >"$dir/mc.out" 2>&1
  else
    (cd "$dir" && "$BIN" mc) >"$dir/mc.out" 2>&1
  fi
  if [ $? -ne 0 ]; then
    printf '    FAIL: the mc run failed\n'
    tail -12 "$dir/mc.out" | sed 's/^/      /'
    fail=$((fail + 1))
    return 0
  fi

  # Predict always on one rank: it is the reference, so it must not inherit
  # whatever the MPI path might be doing wrong.
  predict_deck >"$dir/input"
  (cd "$dir" && "$BIN" predict) >"$dir/predict.out" 2>&1
  if [ $? -ne 0 ]; then
    printf '    FAIL: the prediction over mc_all.xyz failed\n'
    tail -12 "$dir/predict.out" | sed 's/^/      /'
    fail=$((fail + 1))
    return 0
  fi

  if "$PYTHON" "$here/compare_frames.py" "$dir/mc_all.xyz" "$dir/trajectory_out.xyz" \
       --force-tol "$FORCE_TOL" --energy-tol "$ENERGY_TOL"; then
    printf '    PASS\n'
    pass=$((pass + 1))
  else
    printf '    FAIL\n'
    fail=$((fail + 1))
  fi
}

legs=${*:-r1 r2 r4 molr1 molr2 molr4}
for l in $legs; do
  case $l in
  r1) leg 1 deck r1 ;;
  r2) leg 2 deck r2 ;;
  r4) leg 4 deck r4 ;;
  molr1) leg 1 mol_deck molr1 ;;
  molr2) leg 2 mol_deck molr2 ;;
  molr4) leg 4 mol_deck molr4 ;;
  *) die "no such leg: $l (have r1 r2 r4 molr1 molr2 molr4)" ;;
  esac
done

printf '\npassed: %d   failed: %d\n' "$pass" "$fail"
[ "$fail" -eq 0 ]
