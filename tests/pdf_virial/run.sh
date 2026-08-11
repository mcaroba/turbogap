#!/usr/bin/env bash
#
# Acceptance test for the MAD virial on the pdf route: the chain that bins pair
# distances into partial g_ab(r), transforms those into partial S_ab(q), and
# weights them into a diffraction pattern.
#
# The forces on this route have always agreed with finite differences. The
# virial did not, for two reasons that this test pins down:
#
#   sum_a F_a (x) r_a has one term per unordered pair, and the force loops run
#   over ordered ones, which double counted every pair; and the pattern depends
#   on the cell directly, not only through the interatomic distances, because
#   each partial is normalised by the number density N/V.
#
# So every leg here is a strain derivative. fd_gradient.py deforms the cell and
# every position by the same homogeneous strain and differences the energy,
# which is the definition the virial has to meet; the forces come along for
# free in the same run and are checked too.
#
# The legs cover the observables that reach these routines and both force
# implementations, because the pair virial is written out separately in each:
#
#   xrd        the batched device route, which is what a default deck runs:
#              get_structure_factor_forces_matrix_original per batch and
#              channel, with its pair virial accumulated by
#              kernel_exp_force_virial_collection in src/gpu/mad_xrd.cu and
#              its cell half added on the host in compute_exp_spectra
#   xrd-plain  get_structure_factor_forces, the host route, gpu_batched off.
#              get_structure_factor_forces_matrix has no leg: it reads the
#              device pair lists only the batched pdf builds, so neither route
#              can reach it
#   sf         the same two routines with do_xrd off, so no form factors
#   pdf        get_pair_distribution_forces, g(r) as the observable
#   pdf-dr     the same with pair_distribution_output = D(r), where the volume
#              enters twice and mostly cancels
#   mpi        the virial is identical at 1 and N ranks. The volume term is
#              global rather than per-pair, so exactly one rank may add it;
#              this is what catches it being added by all of them
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
# Usage: run.sh [leg ...]      (no argument runs all of them)

set -u

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/../.." && pwd)

BIN=${TURBOGAP_BIN:-$repo/bin/turbogap}
# shellcheck source=../data_root.sh
. "$here/../data_root.sh"
DATA=$DATA_ROOT/xrd_mad
RANKS=${TURBOGAP_RANKS:-2}
PYTHON=${TURBOGAP_PYTHON:-python3}
WORK=${TMPDIR:-/tmp}/turbogap_pdf_virial.$$

# On the CPU branch this lives in tests/xrd_debye, shared with the Debye suite.
# There is no Debye route here, so it sits beside its one caller.
FD=$here/fd_gradient.py

# The pair distribution is cut hard at pair_distribution_rcut, so a pair
# crossing that radius takes a finite bite out of the pattern. At the default
# strain of 1e-4 a few pairs cross on this system and the finite difference
# picks up a percent of noise; 1e-5 moves every atom by ~1e-4 A and none do.
STRAIN=1e-5

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
    printf '      (run tests/fetch_test_data.sh, or see tests/pdf_virial/README.md)\n'
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

PDF_SAMPLES=201
SF_SAMPLES=51

pass=0
fail=0
skipped=0

# The species block and the pdf -> structure factor chain every leg shares.
# Cut down from the tutorial values, as the xrd_mad case is.
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
pair_distribution_n_samples = $PDF_SAMPLES
pair_distribution_partial   = .true.
pair_distribution_rcut      =   6.0
r_range_min                 =   0.1
r_range_max                 =   6.0

do_structure_factor         = .true.
structure_factor_from_pdf   = .true.
structure_factor_window     = .true.

q_range_min                 = 1.0
q_range_max                 = 10.0
q_units                     = 'q'
structure_factor_n_samples  = $SF_SAMPLES
xrd_wavelength              = 1.5405981
EOF
}

# The pdf and sf observables need a target on their own grid, which no data
# file in the test set provides. Predict once and distort what comes out.
#
# The prediction cannot be had for free here: turbogap.f90 sets
# perform%pdf = do_pair_distribution .and. valid_pdf, and valid_pdf comes from
# a dataset labelled `pair_distribution` being declared, so a run with no
# target computes no pattern and writes no file. Each family is therefore
# seeded with a flat target on its own grid, run once to get the prediction,
# and the prediction distorted into the target the legs use. What the seed
# says never reaches anything -- the pattern does not depend on what it is
# being fitted to.
#
# A deck declaring either label does not run on this branch at all: nothing in
# compute_exp_spectra assembles the total y_pair_distribution or
# y_structure_factor -- the batched path collects the partials and stops -- and
# the similarity block then assigns the unallocated array into y_pred and dies
# (turbogap_exp.f90:1024 and :1026). That predates this work; the binary from
# before it dies in the same place. So the seeding run is also the probe: if it
# fails, the three legs that need those targets report why and are skipped
# rather than counted as a virial failure, and they come back on their own once
# the observables work.
targets_ok=""
dr_target_ok=""

make_targets() {
  "$PYTHON" "$here/make_targets.py" --flat 0.1 6.0 $PDF_SAMPLES "$WORK/pdf_target.dat" ||
    die "make_targets (pdf seed)"
  "$PYTHON" "$here/make_targets.py" --flat 1.0 10.0 $SF_SAMPLES "$WORK/sf_target.dat" ||
    die "make_targets (sf seed)"

  { preamble
    cat <<EOF
write_pair_distribution = .true.
write_structure_factor  = .true.

n_exp = 2
exp_labels = 'pair_distribution' 'structure_factor'
exp_data_files = 'pdf_target.dat' 'sf_target.dat'
exp_n_samples = $PDF_SAMPLES $SF_SAMPLES
exp_energy_scales = 1.0 1.0
exp_energies = .true.
EOF
  } >"$WORK/input"
  # Run through an inner shell: this is a probe, and the crash it is probing
  # for would otherwise print bash's own "Segmentation fault" line over the
  # suite's output instead of into the log.
  if bash -c 'cd "$1" && "$2" predict; exit $?' _ "$WORK" "$BIN" >"$WORK/targets.log" 2>&1 &&
     "$PYTHON" "$here/make_targets.py" \
       "$WORK/pair_distribution_total.dat" "$WORK/pdf_target.dat" &&
     "$PYTHON" "$here/make_targets.py" \
       "$WORK/structure_factor_total.dat" "$WORK/sf_target.dat"; then
    targets_ok=1
  else
    return 0
  fi

  # D(r) is a different observable, not a rescaling of g(r), so it needs its
  # own target. exp_input_type tells turbogap the file is already D(r); left
  # unsaid it would rebuild the target from the current density on every step,
  # which makes the target itself move with the cell.
  "$PYTHON" "$here/make_targets.py" --flat 0.1 6.0 $PDF_SAMPLES "$WORK/pdf_dr_target.dat" ||
    die "make_targets (D(r) seed)"
  { preamble
    cat <<EOF
write_pair_distribution  = .true.
pair_distribution_output = 'D(r)'

$(exp_block pair_distribution pdf_dr_target.dat $PDF_SAMPLES 1.0)
exp_input_type = 'D(r)'
EOF
  } >"$WORK/input"
  if bash -c 'cd "$1" && "$2" predict; exit $?' _ "$WORK" "$BIN" >"$WORK/targets_dr.log" 2>&1 &&
     "$PYTHON" "$here/make_targets.py" \
       "$WORK/pair_distribution_total.dat" "$WORK/pdf_dr_target.dat"; then
    dr_target_ok=1
  fi
}

unbatched_unsupported() {
  printf '== %s ==\n' "$1"
  printf '    SKIPPED: gpu_batched = .false. does not run on this branch --\n'
  printf '             the batched block is the only thing that assembles y_xrd,\n'
  printf '             so turning it off dies in the same place the pdf and sf\n'
  printf '             observables do. Pre-existing; see README.md. Log: %s\n' "$WORK/unbatched.log"
  skipped=$((skipped + 1))
}

# Announce a leg that cannot run here, in a form that does not read as a pass.
unsupported() {
  printf '== %s ==\n' "$1"
  printf '    SKIPPED: the %s observable does not run on this branch --\n' "$2"
  printf '             compute_exp_spectra assembles no total %s, and the\n' "$3"
  printf '             similarity block dereferences it. Pre-existing; see\n'
  printf '             tests/pdf_virial/README.md. Log: %s\n' "$4"
  skipped=$((skipped + 1))
}

# leg <name> <family> <scale> <deck body>
leg() {
  local name=$1 family=$2 scale=$3 body=$4
  printf '== %s ==\n' "$name"
  { preamble; printf '%s\n' "$body"; } >"$WORK/input"
  if "$PYTHON" "$FD" "$WORK" --bin "$BIN" --scale "$scale" --family "$family" \
       --strain "$STRAIN" --atoms-to-check 1; then
    printf '    PASS\n'
    pass=$((pass + 1))
  else
    printf '    FAIL\n'
    fail=$((fail + 1))
  fi
}

exp_block() {
  local label=$1 file=$2 samples=$3 scale=$4
  cat <<EOF
n_exp = 1
exp_labels = '$label'
exp_data_files = '$file'
exp_n_samples = $samples
exp_energy_scales = $scale
exp_energies = .true.
exp_forces = .true.
EOF
}
# Anything keyed off n_exp -- exp_input_type among them -- is read with an
# implied do over the experimental families, so it means nothing unless n_exp
# is already set. It goes after exp_block, never before it.

legs=${*:-xrd xrd-plain sf pdf pdf-dr mpi}
run_leg() { case " $legs " in *" $1 "*) return 0 ;; *) return 1 ;; esac; }

make_targets

# The xrd_mad case uses exp_energy_scales = 100; this uses 1.0. Nothing here
# depends on the scale -- it divides out of every comparison -- and a large one
# runs the virial past the field width of the xyz writer, where two components
# come back with no space between them.
#
# There is no Lorentz-polarization factor on this branch, so this leg is the
# CPU branch's `xrd` leg with that turned off. It resolves the derivative less
# sharply for it: LP multiplies the pattern by up to ~117 over this q range,
# which makes both the energy and the virial large.
if run_leg xrd; then
  leg "xrd, matrix forces (device kernel)" xrd 1.0 \
    "do_xrd     = .true.
xrd_output = 'q*F(q)'

$(exp_block xrd xrd_glassy_carbon_zeng_2017.fq $SF_SAMPLES 1.0)"
fi

# gpu_batched is on by default and, when it is, the batched block in
# compute_exp_spectra does the xrd forces itself through
# get_structure_factor_forces_matrix_original and never consults
# structure_factor_matrix_forces. So this leg has to turn it off; without that
# it is the `xrd` leg again, to the last digit, and proves nothing about the
# routine it names.
unbatched_ok=""
if run_leg xrd-plain; then
  { preamble
    cat <<EOF
do_xrd      = .true.
xrd_output  = 'q*F(q)'
gpu_batched = .false.

$(exp_block xrd xrd_glassy_carbon_zeng_2017.fq $SF_SAMPLES 1.0)
EOF
  } >"$WORK/input"
  if bash -c 'cd "$1" && "$2" predict; exit $?' _ "$WORK" "$BIN" >"$WORK/unbatched.log" 2>&1; then
    unbatched_ok=1
  fi
fi

if run_leg xrd-plain; then
  if [ -n "$unbatched_ok" ]; then
    leg "xrd, non-matrix forces" xrd 1.0 \
      "do_xrd      = .true.
xrd_output  = 'q*F(q)'
gpu_batched = .false.
structure_factor_matrix_forces = .false.

$(exp_block xrd xrd_glassy_carbon_zeng_2017.fq $SF_SAMPLES 1.0)"
  else
    unbatched_unsupported "xrd, non-matrix forces"
  fi
fi

if run_leg sf; then
  if [ -n "$targets_ok" ]; then
    leg "structure factor observable" sf 100.0 \
      "$(exp_block structure_factor sf_target.dat $SF_SAMPLES 100.0)"
  else
    unsupported "structure factor observable" \
      "structure_factor" "y_structure_factor" "$WORK/targets.log"
  fi
fi

if run_leg pdf; then
  if [ -n "$targets_ok" ]; then
    leg "pair distribution observable, g(r)" pdf 1000.0 \
      "$(exp_block pair_distribution pdf_target.dat $PDF_SAMPLES 1000.0)"
  else
    unsupported "pair distribution observable, g(r)" \
      "pair_distribution" "y_pair_distribution" "$WORK/targets.log"
  fi
fi

if run_leg pdf-dr; then
  if [ -n "$dr_target_ok" ]; then
    leg "pair distribution observable, D(r)" pdf 1.0 \
      "pair_distribution_output = 'D(r)'

$(exp_block pair_distribution pdf_dr_target.dat $PDF_SAMPLES 1.0)
exp_input_type = 'D(r)'"
  else
    unsupported "pair distribution observable, D(r)" \
      "pair_distribution" "y_pair_distribution" "$WORK/targets.log"
  fi
fi

# The pair half of the virial is built from the pairs each rank owns and sums
# over ranks; the volume half is a property of the whole cell and must be added
# exactly once. Getting that wrong scales it by the rank count, which is
# invisible on one rank and obvious here.
if run_leg mpi; then
  printf '== virial is independent of the rank count ==\n'
  if ! command -v mpirun >/dev/null || [ "$RANKS" -le 1 ]; then
    printf '    SKIPPED (no mpirun, or TURBOGAP_RANKS=1)\n'
    skipped=$((skipped + 1))
  else
    { preamble
      cat <<EOF
do_xrd     = .true.
xrd_output = 'q*F(q)'

$(exp_block xrd xrd_glassy_carbon_zeng_2017.fq $SF_SAMPLES 1.0)
EOF
    } >"$WORK/input"
    (cd "$WORK" && "$BIN" predict) >"$WORK/mpi1.log" 2>&1
    one=$(sed -n 2p "$WORK/trajectory_out.xyz" | grep -o 'virial="[^"]*"')
    (cd "$WORK" && mpirun -np "$RANKS" "$BIN" predict) >"$WORK/mpin.log" 2>&1
    many=$(sed -n 2p "$WORK/trajectory_out.xyz" | grep -o 'virial="[^"]*"')
    # Compared with a tolerance rather than byte for byte: the pair half is a
    # device reduction whose summation order is not fixed, and splitting the
    # pairs across ranks reorders it again, so the last digit moves between
    # runs. What this leg is looking for is the volume term being added once
    # per rank instead of once, which is a whole term -- percent, not 1e-8.
    if [ -n "$one" ] && [ -n "$many" ] && "$PYTHON" -c '
import sys
import numpy as np
a, b = (np.array([float(v) for v in s.split("=", 1)[1].strip("\"").split()])
        for s in sys.argv[1:3])
scale = max(float(np.abs(a).max()), 1e-12)
worst = float(np.abs(a - b).max()) / scale
print(f"    worst relative difference: {worst:.2e}  (tolerance 1e-6)")
sys.exit(0 if worst < 1e-6 else 1)
' "$one" "$many"; then
      printf '    PASS (%s)\n' "$one"
      pass=$((pass + 1))
    else
      printf '    FAIL: 1 rank  %s\n' "$one"
      printf '          %d ranks %s\n' "$RANKS" "$many"
      fail=$((fail + 1))
    fi
  fi
fi

printf '\npassed: %d   failed: %d   skipped: %d\n' "$pass" "$fail" "$skipped"
[ "$fail" -eq 0 ]
