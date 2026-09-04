#!/usr/bin/env bash
#
# poly3operator radial expansion: the coefficients a turbogap build actually
# produces, against Gauss-Legendre quadrature of the same integral in real128.
#
# The quadrature is not an approximation. The integrand is a polynomial of
# degree at most alpha+8 on each panel between its breakpoints -- rj-width, rj,
# rcut_soft, rj+width -- and 16-node Gauss-Legendre is exact through degree 31.
# What the comparison measures is therefore the closed form's own conditioning,
# dominated by the change of basis: W has entries of 400 at alpha_max 4 and
# 5e5 at 8, for a result of order 1.
#
# Nothing here compares the code to itself: the reference integrates the
# density numerically instead of by parts, the tabulated S and W are checked
# against the analytic overlap matrix rather than trusted, and the derivative
# is checked by an h-scan that asserts second-order convergence rather than a
# tolerance.
#
# Two configurations, because one covers less than half the arithmetic:
#
#   GST   alpha_max 4, atom_sigma_r_scaling 0, amplitude_scaling 1
#   CO    alpha_max 8, atom_sigma_r_scaling 0.1, amplitude_scaling 2
#
# The CO potential is fitted with poly3gauss; this switches its basis in a
# scratch copy, which makes its energies meaningless and its descriptor exactly
# the thing under test. Without it the sigma-scaling terms in the density width
# and in the amplitude derivative never run -- and an early version of this
# test passed a reference that had dropped them both.
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/../.." && pwd)
work=$(mktemp -d)
trap 'if [ "${TURBOGAP_KEEP:-0}" = "1" ]; then echo "kept $work"; else rm -rf "$work"; fi' EXIT

BIN=${TURBOGAP_BIN:-$repo/bin/turbogap}
DATA=${TURBOGAP_POLY3OP_DATA:-$repo/../turbogap_tests/GST}
CO_DATA=${TURBOGAP_POLY3OP_CO_DATA:-$repo/../turbogap_tests/CO}
F90=${TURBOGAP_F90:-mpif90}
# soap_turbo_radial.f90 also holds the LAPACK-based orthonormalisation, so the
# module needs a LAPACK at link time even though this test uses the tabulated
# one. Which library provides it is a per-machine answer: -llapack -lblas on a
# Debian box, -lopenblas on Roihu.
LIBS=${TURBOGAP_TEST_LIBS:--llapack -lblas}
TOL=${TURBOGAP_POLY3OP_TOL:-1e-8}

if [ ! -x "$BIN" ]; then
  echo "SKIP: no binary at $BIN (set TURBOGAP_BIN)"
  exit 0
fi

echo "--- building poly3operatorverify ---"
$F90 -O2 -J "$work" -c "$repo/src/soap_turbo/src/soap_turbo_functions.f90" -o "$work/f.o"
$F90 -O2 -J "$work" -c "$repo/src/soap_turbo/src/soap_turbo_radial.f90" -o "$work/r.o"
$F90 -O2 -I "$work" "$here/poly3operatorverify.f90" "$work/f.o" "$work/r.o" \
     -o "$work/poly3operatorverify" $LIBS

# Read a scalar hyper out of a .gap file: first value of the first match.
get() { grep -m1 "^$2 *=" "$1" | sed 's/.*= *//' | tr -d '"' | awk '{print $1}'; }

failed=0

check() {   # check <label> <datadir> <atoms> <gapfile basename> <input body> <flip basis>
  local label=$1 data=$2 atoms=$3 gapname=$4 body=$5 flip=$6
  local run="$work/$label"

  if [ ! -d "$data/gap_files" ] || [ ! -f "$data/$atoms" ]; then
    echo "SKIP $label: no data at $data"
    return
  fi

  mkdir -p "$run/gap_files"
  cp "$data/$atoms" "$run/atoms.xyz"
  find "$data/gap_files/" -maxdepth 1 -type f -exec cp {} "$run/gap_files/" \;
  if [ "$flip" = yes ]; then
    sed -i 's/basis = "poly3gauss"/basis = "poly3operator"/' "$run/gap_files/$gapname"
  fi
  printf '%s\n' "$body" > "$run/input"

  local gapfile="$run/gap_files/$gapname"
  local basis alpha_max rcut_hard buffer atom_sigma sigma_scaling amp_scaling enh cw rcut_soft
  basis=$(get "$gapfile" basis)
  if [ "$basis" != poly3operator ]; then
    echo "SKIP $label: $gapname is basis=$basis"
    return
  fi
  alpha_max=$(get "$gapfile" n_max)
  rcut_hard=$(get "$gapfile" rcut)
  buffer=$(get "$gapfile" buffer)
  atom_sigma=$(get "$gapfile" atom_sigma_r)
  sigma_scaling=$(get "$gapfile" atom_sigma_r_scaling)
  amp_scaling=$(get "$gapfile" amplitude_scaling)
  enh=$(get "$gapfile" radial_enhancement)
  cw=$(get "$gapfile" central_weight)
  rcut_soft=$(awk -v a="$rcut_hard" -v b="$buffer" 'BEGIN{printf "%.12f", a-b}')

  echo
  echo "############ $label ############"
  echo "  alpha_max=$alpha_max rcut_hard=$rcut_hard rcut_soft=$rcut_soft atom_sigma=$atom_sigma"
  echo "  atom_sigma_scaling=$sigma_scaling amplitude_scaling=$amp_scaling radial_enhancement=$enh central_weight=$cw"

  ( cd "$run" && TURBOGAP_DUMP_RADIAL=1 "$BIN" predict > out.log 2> err.log ) || {
    echo "FAIL $label: turbogap predict exited non-zero"; tail -20 "$run/err.log"; failed=1; return; }
  if [ ! -s "$run/radial_exp_coeff_dump.dat" ]; then
    echo "FAIL $label: no radial_exp_coeff_dump.dat was written"; failed=1; return
  fi

  "$work/poly3operatorverify" "$run/radial_exp_coeff_dump.dat" \
      "$alpha_max" "$rcut_hard" "$rcut_soft" "$atom_sigma" "$sigma_scaling" \
      "$amp_scaling" "$enh" "$cw" "$TOL" || failed=1
}

check GST_alpha4 "$DATA" atoms_897.xyz soap_turbo_pot_v8_3.gap \
'atoms_file = "atoms.xyz"
pot_file = "gap_files/soap_turbo_pot_v8_3.gap"
n_species = 3
species = Ge Sb Te
masses = 72.64 121.76 127.6
e0 = 0.0 0.0 0.0
random_seed = 12345
write_xyz = 1' no

check CO_alpha8_sigma_scaled "$CO_DATA" atoms_7176.xyz CO.gap \
'atoms_file = "atoms.xyz"
pot_file = "gap_files/CO.gap"
n_species = 2
species = C O
masses = 12.01 15.99
e0 = -.16138053 0.
random_seed = 12345
write_xyz = 1' yes

echo
if [ "$failed" = 0 ]; then
  echo "poly3operator: PASS"
else
  echo "poly3operator: FAIL"
  exit 1
fi
