#!/usr/bin/env bash
#
# Does randomize_velocities draw from the distribution its keyword names?
#
# velocity_distribution = "maxwell" has to give Maxwell-Boltzmann velocities,
# and "uniform" has to keep giving the historical draw, which is uniform on a
# box and then scaled to an exact kinetic energy.
#
# A temperature check cannot tell the two apart: both are constructed to have
# the right mean kinetic energy, which is exactly why "uniform" sat under
# mc_hamiltonian unnoticed. So each leg checks the *shape* of the distribution,
# through the excess kurtosis of the mass-scaled components -- 0 for a
# Gaussian, -1.2 for a box. See check_draw.py.
#
# The distinction is not cosmetic. With mc_hamiltonian the momenta are part of
# the state the Metropolis test accepts or rejects, and detailed balance holds
# only if they come from the distribution the test assumes.
#
# Environment:
#   TURBOGAP_BIN        binary under test (default: <repo>/bin/turbogap)
#   TURBOGAP_DATA_ROOT  directory holding the test systems. Default: the
#                       turbogap_tests clone beside this repository, which
#                       tests/fetch_test_data.sh creates on the first run
#                       (override where it goes with TURBOGAP_TESTS_DIR).
#   TURBOGAP_PYTHON     interpreter for the check (default: python3). Standard
#                       library only.
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
WORK=${TMPDIR:-/tmp}/turbogap_velocity_draw.$$

TEMP=600

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 2
}

[ -x "$BIN" ] || die "binary under test not found or not executable: $BIN"
command -v "$PYTHON" >/dev/null || die "python interpreter not found: $PYTHON"

for f in atoms.xyz gap_files/CO.gap; do
  if [ ! -e "$DATA/$f" ]; then
    printf 'SKIP: missing test data %s\n' "$DATA/$f"
    printf '      (run tests/fetch_test_data.sh, or see tests/velocity_draw/README.md)\n'
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

# md_nsteps = 1 with optimize = vv: the frame written at md_istep = 0 carries
# the velocities as drawn, because velocity Verlet's first step leaves them
# alone and the trajectory writer emits positions_prev and velocities, which are
# synchronous. md_nsteps = 0 would be the obvious thing to ask for and divides
# by zero in the progress bar.
# masses have to be in the output for the mass-scaling, hence write_masses.
leg() {
  local name=$1 dist=$2 shape=$3
  local dir=$WORK/$name
  printf '== %s ==\n' "$name"

  mkdir -p "$dir"
  ln -sf "$DATA/gap_files" "$dir/gap_files"
  ln -sf "$DATA/atoms.xyz" "$dir/atoms.xyz"
  cat >"$dir/input" <<EOF
atoms_file = 'atoms.xyz'
pot_file = 'gap_files/CO.gap'
n_species = 2
species = C O
masses = 12.01 15.99
e0 = -.16138053 0.
random_seed = 12345
max_Gbytes_per_process = 1.0

md_nsteps = 1
md_step = 0.5
optimize = 'vv'
t_beg = $TEMP
t_end = $TEMP
randomize_velocities = .true.
velocity_distribution = '$dist'
write_xyz = 1
write_masses = .true.
EOF

  (cd "$dir" && "$BIN" md) >"$dir/run.log" 2>&1
  if [ $? -ne 0 ]; then
    printf '    FAIL: the run failed\n'
    tail -12 "$dir/run.log" | sed 's/^/      /'
    fail=$((fail + 1))
    return 0
  fi

  if "$PYTHON" "$here/check_draw.py" "$dir/trajectory_out.xyz" \
       --shape "$shape" --temperature "$TEMP"; then
    printf '    PASS\n'
    pass=$((pass + 1))
  else
    printf '    FAIL\n'
    fail=$((fail + 1))
  fi
}

# A deck naming a distribution that does not exist has to stop, not fall back
# to a default and sample something the user did not ask for.
leg_reject() {
  local dir=$WORK/reject
  printf '== reject ==\n'
  mkdir -p "$dir"
  ln -sf "$DATA/gap_files" "$dir/gap_files"
  ln -sf "$DATA/atoms.xyz" "$dir/atoms.xyz"
  cat >"$dir/input" <<EOF
atoms_file = 'atoms.xyz'
pot_file = 'gap_files/CO.gap'
n_species = 2
species = C O
masses = 12.01 15.99
e0 = -.16138053 0.
md_nsteps = 1
randomize_velocities = .true.
velocity_distribution = 'gaussian'
EOF
  (cd "$dir" && "$BIN" md) >"$dir/run.log" 2>&1
  local rc=$?
  local msg
  msg=$(grep -m1 'Invalid velocity_distribution' "$dir/run.log" | sed 's/^ *//')
  # A non-zero exit on its own would not do: the deck has to be rejected *for
  # this reason*, not because something else fell over first.
  if [ $rc -eq 0 ] || [ -z "$msg" ]; then
    printf '    FAIL: expected a rejection naming velocity_distribution (exit %d)\n' "$rc"
    tail -6 "$dir/run.log" | sed 's/^/      /'
    fail=$((fail + 1))
  else
    printf '    PASS   %s\n' "$msg"
    pass=$((pass + 1))
  fi
}

legs=${*:-maxwell uniform reject}
for l in $legs; do
  case $l in
  maxwell) leg maxwell maxwell gaussian ;;
  uniform) leg uniform uniform box ;;
  reject) leg_reject ;;
  *) die "no such leg: $l (have maxwell uniform reject)" ;;
  esac
done

printf '\npassed: %d   failed: %d\n' "$pass" "$fail"
[ "$fail" -eq 0 ]
