#!/usr/bin/env bash
#
# Acceptance test for the Debye scattering route (xrd_debye = .true.).
#
# This is deliberately NOT part of tests/regression, for the same reason
# tests/dipole is not: that suite compares byte-for-byte against a frozen
# baseline binary which does not know the keyword and exits on it, and its
# comparator is plain diff, which has no tolerance. Those cases
# (xrd_debye_mad, xrd_debye_mad_mpi2) pin the Debye output against drift and
# say nothing about whether it is right. This one asks whether it is right:
#
#   forward   the predicted pattern against an independent Debye sum written
#             from scratch in numpy, over the four xrd_output conventions,
#             with and without the Lorch window, with and without the
#             Lorentz-polarization factor, for X-rays and for neutrons, and
#             under MPI
#   gradient  the forces and the virial against central finite differences of
#             the energy, with the experimental term isolated by differencing
#             two runs at different energy scales
#
# Environment:
#   TURBOGAP_BIN        binary under test (default: <repo>/bin/turbogap)
#   TURBOGAP_DATA_ROOT  directory holding the test systems
#                       (default: $HOME/work/cpu_vs_gpu_tests/input)
#   TURBOGAP_RANKS      ranks for the MPI leg (default: 2)
#   TURBOGAP_PYTHON     interpreter for the reference (default: python3). The
#                       reference needs numpy; on a machine where the system
#                       python3 has none, point this at a venv that does.
#   TURBOGAP_KEEP       keep the staging directory for inspection
#
# Usage: run.sh [forward|gradient]      (no argument runs both)

set -u

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/../.." && pwd)

BIN=${TURBOGAP_BIN:-$repo/bin/turbogap}
DATA_ROOT=${TURBOGAP_DATA_ROOT:-$HOME/work/cpu_vs_gpu_tests/input}
DATA=$DATA_ROOT/xrd_mad
RANKS=${TURBOGAP_RANKS:-2}
PYTHON=${TURBOGAP_PYTHON:-python3}
WORK=${TMPDIR:-/tmp}/turbogap_xrd_debye.$$

which=${1:-both}

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 2
}

[ -x "$BIN" ] || die "binary under test not found or not executable: $BIN"
command -v "$PYTHON" >/dev/null || die "python interpreter not found: $PYTHON"
"$PYTHON" -c 'import numpy' 2>/dev/null ||
  die "$PYTHON has no numpy, which the reference needs (set TURBOGAP_PYTHON)"

# A case whose data is absent is skipped, not failed -- same contract as
# tests/regression/run.sh.
for f in atoms.xyz gap_files/CO.gap xrd_glassy_carbon_zeng_2017.fq; do
  if [ ! -e "$DATA/$f" ]; then
    printf 'SKIP: missing test data %s\n' "$DATA/$f"
    printf '      (set TURBOGAP_DATA_ROOT, or see tests/xrd_debye/README.md)\n'
    exit 0
  fi
done

mkdir -p "$WORK"
if [ -z "${TURBOGAP_KEEP:-}" ]; then
  trap 'rm -rf "$WORK"' EXIT
else
  trap 'printf "staging kept in %s\n" "$WORK"' EXIT
fi

ln -sf "$DATA/gap_files" "$WORK/gap_files"
ln -sf "$DATA/xrd_glassy_carbon_zeng_2017.fq" "$WORK/xrd_glassy_carbon_zeng_2017.fq"
cp "$DATA/atoms.xyz" "$WORK/atoms.xyz"

NSAMPLES=51
QMIN=1.0
QMAX=10.0

pass=0
fail=0

# The species block every deck needs.
preamble() {
  cat <<'EOF'
atoms_file = 'atoms.xyz'
pot_file = 'gap_files/CO.gap'
n_species = 2
species = C O
masses = 12.01 15.99
e0 = -.16138053 0.
EOF
}

# forward <name> <keywords> <prediction file> <reference args...>
forward() {
  local name=$1 keywords=$2 pred=$3
  shift 3
  printf '== forward: %s ==\n' "$name"
  { preamble
    printf 'xrd_debye = .true.\n'
    printf 'q_range_min = %s\nq_range_max = %s\nq_units = %s\n' "$QMIN" "$QMAX" "'q'"
    printf 'xrd_n_samples = %s\n' "$NSAMPLES"
    printf '%s\n' "$keywords"
  } >"$WORK/input"

  rm -f "$WORK/xrd_prediction.dat" "$WORK/nd_prediction.dat"
  if ! (cd "$WORK" && "$BIN" predict) >"$WORK/run.log" 2>&1; then
    printf '    FAIL: turbogap exited non-zero\n'
    tail -12 "$WORK/run.log" | sed 's/^/      /'
    fail=$((fail + 1))
    return
  fi
  if [ ! -f "$WORK/$pred" ]; then
    printf '    FAIL: %s was not written\n' "$pred"
    fail=$((fail + 1))
    return
  fi
  if "$PYTHON" "$here/debye_reference.py" "$WORK/atoms.xyz" "$WORK/$pred" \
       --n "$NSAMPLES" --qmin "$QMIN" --qmax "$QMAX" "$@"; then
    printf '    PASS\n'
    pass=$((pass + 1))
  else
    printf '    FAIL\n'
    fail=$((fail + 1))
  fi
}

if [ "$which" = forward ] || [ "$which" = both ]; then
  # The raw intensity, no window: the Debye equation with nothing on top.
  forward "raw intensity" \
    "do_xrd = .true.
write_xrd = .true.
xrd_output = 'xrd'
structure_factor_window = .false." \
    xrd_prediction.dat --output xrd

  # The two normalised conventions, which have to mean the same thing here as
  # they do on the pdf/sf route.
  forward "q*F(q)" \
    "do_xrd = .true.
write_xrd = .true.
xrd_output = 'q*F(q)'
structure_factor_window = .false." \
    xrd_prediction.dat --output 'q*F(q)'

  forward "F(q) with Lorch window" \
    "do_xrd = .true.
write_xrd = .true.
xrd_output = 'F(q)'
structure_factor_window = .true.
xrd_rcut = 6.0" \
    xrd_prediction.dat --output 'F(q)' --window --rcut 6.0

  # The Lorentz-polarization factor, at a polarization that is neither 0 nor 1
  # so a dropped or defaulted parameter shows up.
  forward "Lorentz-polarization, P = 0.8" \
    "do_xrd = .true.
write_xrd = .true.
xrd_output = 'xrd'
structure_factor_window = .false.
xrd_wavelength = 1.5405981
xrd_lorentz_polarization = .true.
xrd_lp_polarization = 0.8" \
    xrd_prediction.dat --output xrd --lp --polarization 0.8

  # Neutrons: the same sum with q-independent scattering lengths.
  forward "neutron diffraction" \
    "do_nd = .true.
write_nd = .true.
nd_output = 'xrd'
structure_factor_window = .false." \
    nd_prediction.dat --output xrd --neutron

  # Under MPI the double sum is split over the outer atom index and reduced.
  if command -v mpirun >/dev/null && [ "$RANKS" -gt 1 ]; then
    printf '== forward: q*F(q) with window, %s ranks ==\n' "$RANKS"
    { preamble
      printf 'xrd_debye = .true.\n'
      printf 'q_range_min = %s\nq_range_max = %s\nq_units = %s\n' "$QMIN" "$QMAX" "'q'"
      printf 'xrd_n_samples = %s\n' "$NSAMPLES"
      cat <<'EOF'
do_xrd = .true.
write_xrd = .true.
xrd_output = 'q*F(q)'
structure_factor_window = .true.
xrd_rcut = 6.0
EOF
    } >"$WORK/input"
    rm -f "$WORK/xrd_prediction.dat"
    if (cd "$WORK" && mpirun -np "$RANKS" "$BIN" predict) >"$WORK/run.log" 2>&1 &&
       "$PYTHON" "$here/debye_reference.py" "$WORK/atoms.xyz" "$WORK/xrd_prediction.dat" \
         --n "$NSAMPLES" --qmin "$QMIN" --qmax "$QMAX" \
         --output 'q*F(q)' --window --rcut 6.0; then
      printf '    PASS\n'
      pass=$((pass + 1))
    else
      printf '    FAIL\n'
      tail -12 "$WORK/run.log" | sed 's/^/      /'
      fail=$((fail + 1))
    fi
  else
    printf '== forward: MPI leg SKIPPED (no mpirun, or TURBOGAP_RANKS=1) ==\n'
  fi
fi

if [ "$which" = gradient ] || [ "$which" = both ]; then
  printf '== gradient: forces and virial against finite differences ==\n'
  { preamble
    cat <<'EOF'
do_xrd = .true.
xrd_debye = .true.
structure_factor_window = .true.
xrd_rcut = 6.0
xrd_output = 'q*F(q)'

q_range_min = 1.0
q_range_max = 10.0
q_units = 'q'

n_exp = 1
exp_labels = 'xrd'
exp_data_files = 'xrd_glassy_carbon_zeng_2017.fq'
exp_n_samples = 51
exp_energy_scales = 100.0
exp_energies = .true.
exp_forces = .true.
EOF
  } >"$WORK/input"

  if "$PYTHON" "$here/fd_gradient.py" "$WORK" --bin "$BIN" --scale 100.0; then
    printf '    PASS\n'
    pass=$((pass + 1))
  else
    printf '    FAIL\n'
    fail=$((fail + 1))
  fi
fi

printf '\npassed: %d   failed: %d\n' "$pass" "$fail"
[ "$fail" -eq 0 ]
