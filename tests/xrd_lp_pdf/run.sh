#!/usr/bin/env bash
#
# Acceptance test for the Lorentz-polarization factor on the pdf/sf route,
# i.e. xrd_lorentz_polarization = .true. with xrd_debye off, where the pattern
# comes from the pair distribution through the partial structure factors.
#
# The Debye route has its own test (tests/xrd_debye), which checks the factor
# against a pattern rebuilt from scratch in numpy. Here the pattern itself is
# not the question -- the pdf chain is what the xrd_mad regression case pins.
# The question is the weight applied on top of it and, above all, whether the
# gradient carries the same weight as the energy:
#
#   forward   the same prediction with the factor off and on, checked against
#             LP(theta) recomputed from the q grid, for the raw intensity and
#             for a normalised output, for X-rays and for neutrons, and under
#             MPI where the q grid is split across ranks and reduced before
#             the factor is applied
#   gradient  the forces and virial of an LP-weighted MAD run against central
#             finite differences of its energy
#
# Environment:
#   TURBOGAP_BIN        binary under test (default: <repo>/bin/turbogap)
#   TURBOGAP_DATA_ROOT  directory holding the test systems. Default: the
#                       turbogap_tests clone beside this repository, which
#                       tests/fetch_test_data.sh creates on the first run
#                       (override where it goes with TURBOGAP_TESTS_DIR).
#   TURBOGAP_RANKS      ranks for the MPI leg (default: 2)
#   TURBOGAP_PYTHON     interpreter for the reference (default: python3). It
#                       needs numpy; on a machine whose system python3 has
#                       none, point this at a venv that does.
#   TURBOGAP_KEEP       keep the staging directory for inspection
#
# Usage: run.sh [forward|gradient]      (no argument runs both)

set -u

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/../.." && pwd)

BIN=${TURBOGAP_BIN:-$repo/bin/turbogap}
# shellcheck source=../data_root.sh
. "$here/../data_root.sh"
DATA=$DATA_ROOT/xrd_mad
RANKS=${TURBOGAP_RANKS:-2}
PYTHON=${TURBOGAP_PYTHON:-python3}
WORK=${TMPDIR:-/tmp}/turbogap_xrd_lp_pdf.$$

which=${1:-both}

WAVELENGTH=1.5405981
# Neither 0 nor 1, so a dropped or defaulted polarization shows up.
POLARIZATION=0.8
# The xrd_mad case uses 100. LP multiplies the pattern by up to ~117 over this
# q range, and the residual enters the energy squared, so the same scale here
# puts the virial past the field width of the xyz writer and the frame comes
# back with two numbers run together. Nothing about the gradient depends on
# the scale -- it divides out of the comparison.
GRADIENT_SCALE=1.0

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
    printf '      (run tests/fetch_test_data.sh, or see tests/xrd_lp_pdf/README.md)\n'
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

# The species block and the pdf -> structure factor chain the pattern comes
# from. Cut down from the tutorial values, as the xrd_mad case is.
preamble() {
  cat <<EOF
atoms_file = 'atoms.xyz'
pot_file = 'gap_files/CO.gap'
n_species = 2
species = C O
masses = 12.01 15.99
e0 = -.16138053 0.

do_pair_distribution        = .true.
pair_distribution_kde_sigma =   0.1
pair_distribution_n_samples =  201
pair_distribution_partial   = .true.
pair_distribution_rcut      =   6.0
r_range_min                 =   0.1
r_range_max                 =   6.0

do_structure_factor         = .true.
structure_factor_from_pdf   = .true.
structure_factor_window     = .true.

q_range_min                 = $QMIN
q_range_max                 = $QMAX
q_units                     = 'q'
structure_factor_n_samples  = $NSAMPLES
xrd_wavelength              = $WAVELENGTH
EOF
}

# run <ranks> <deck body> <prediction file> <destination>
run_pattern() {
  local ranks=$1 body=$2 pred=$3 dest=$4
  { preamble; printf '%s\n' "$body"; } >"$WORK/input"
  rm -f "$WORK/$pred"
  if [ "$ranks" -gt 1 ]; then
    (cd "$WORK" && mpirun -np "$ranks" "$BIN" predict) >"$WORK/run.log" 2>&1
  else
    (cd "$WORK" && "$BIN" predict) >"$WORK/run.log" 2>&1
  fi
  if [ $? -ne 0 ]; then
    printf '    FAIL: turbogap exited non-zero\n'
    tail -12 "$WORK/run.log" | sed 's/^/      /'
    return 1
  fi
  if [ ! -f "$WORK/$pred" ]; then
    printf '    FAIL: %s was not written\n' "$pred"
    return 1
  fi
  mv "$WORK/$pred" "$WORK/$dest"
}

# forward <name> <ranks> <deck body> <prediction file>
#
# Runs the same deck twice, with the factor off and on, and hands both
# patterns to the reference.
forward() {
  local name=$1 ranks=$2 body=$3 pred=$4
  printf '== forward: %s ==\n' "$name"

  if ! run_pattern "$ranks" "$body" "$pred" raw.dat; then
    fail=$((fail + 1))
    return
  fi
  if ! run_pattern "$ranks" "$body
xrd_lorentz_polarization = .true.
xrd_lp_polarization = $POLARIZATION" "$pred" lp.dat; then
    fail=$((fail + 1))
    return
  fi

  if "$PYTHON" "$here/lp_reference.py" "$WORK/raw.dat" "$WORK/lp.dat" \
       --wavelength "$WAVELENGTH" --polarization "$POLARIZATION" --q-units q; then
    printf '    PASS\n'
    pass=$((pass + 1))
  else
    printf '    FAIL\n'
    fail=$((fail + 1))
  fi
}

if [ "$which" = forward ] || [ "$which" = both ]; then
  # The raw powder intensity, which is what the factor is meant to correct.
  forward "raw intensity" 1 \
    "do_xrd = .true.
write_xrd = .true.
xrd_output = 'xrd'" \
    xrd_prediction.dat

  # A normalised output: the factor multiplies whatever xrd_output asked for,
  # after the affine map, not before it.
  forward "q*F(q)" 1 \
    "do_xrd = .true.
write_xrd = .true.
xrd_output = 'q*F(q)'" \
    xrd_prediction.dat

  # Neutrons take the same weight. The polarization part of it is an X-ray
  # notion, so P is a modelling choice there and not a property of the beam;
  # what is checked here is only that nd is weighted like xrd.
  forward "neutron diffraction" 1 \
    "do_nd = .true.
write_nd = .true.
nd_output = 'xrd'" \
    nd_prediction.dat

  # Under MPI each rank builds its own slice of the q grid and the pattern is
  # reduced before the factor is applied, so a factor applied per rank rather
  # than once would show up here as a multiple count.
  if command -v mpirun >/dev/null && [ "$RANKS" -gt 1 ]; then
    forward "raw intensity, $RANKS ranks" "$RANKS" \
      "do_xrd = .true.
write_xrd = .true.
xrd_output = 'xrd'" \
      xrd_prediction.dat
  else
    printf '== forward: MPI leg SKIPPED (no mpirun, or TURBOGAP_RANKS=1) ==\n'
  fi
fi

if [ "$which" = gradient ] || [ "$which" = both ]; then
  printf '== gradient: LP-weighted MAD forces against finite differences ==\n'
  { preamble
    cat <<EOF
do_xrd                      = .true.
xrd_output                  = 'q*F(q)'

xrd_lorentz_polarization = .true.
xrd_lp_polarization = $POLARIZATION

n_exp = 1
exp_labels = 'xrd'
exp_data_files = 'xrd_glassy_carbon_zeng_2017.fq'
exp_n_samples = $NSAMPLES
exp_energy_scales = $GRADIENT_SCALE
exp_energies = .true.
exp_forces = .true.
EOF
  } >"$WORK/input"

  # fd_gradient.py is deck-agnostic: it differences two runs at different
  # energy scales, so it isolates whatever experimental family the deck turns
  # on. Shared with tests/xrd_debye rather than copied.
  #
  # Forces and virial both, since the factor multiplies the whole pattern and
  # so has to reach both. What the strain leg says about the pdf route in
  # general -- and why it needs a smaller strain than the Debye route does --
  # belongs to tests/pdf_virial, which covers it without the factor.
  if "$PYTHON" "$here/../xrd_debye/fd_gradient.py" "$WORK" --bin "$BIN" \
       --scale "$GRADIENT_SCALE" --strain 1e-5; then
    printf '    PASS\n'
    pass=$((pass + 1))
  else
    printf '    FAIL\n'
    fail=$((fail + 1))
  fi
fi

printf '\npassed: %d   failed: %d\n' "$pass" "$fail"
[ "$fail" -eq 0 ]
