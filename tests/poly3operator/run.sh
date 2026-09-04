#!/usr/bin/env bash
#
# poly3operator radial expansion: the coefficients a turbogap build actually
# produces, against Gauss-Legendre quadrature of the same integral in real128.
#
# The build under test writes radial_exp_coeff_dump.dat (TURBOGAP_DUMP_RADIAL)
# for a GST cell, whose 5.5 A cutoff and 0.5 A buffer put roughly a third of
# the neighbours in the buffer region, where the filtered density raises the
# integrand from cubic to degree six. poly3operatorverify then recomputes every
# nonzero coefficient by quadrature and compares.
#
# Nothing here compares the code to itself: the reference integrates the density
# numerically instead of by parts, and the tabulated basis matrices are checked
# against the analytic overlap matrix rather than taken on trust.
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/../.." && pwd)
work=$(mktemp -d)
trap 'if [ "${TURBOGAP_KEEP:-0}" = "1" ]; then echo "kept $work"; else rm -rf "$work"; fi' EXIT

BIN=${TURBOGAP_BIN:-$repo/bin/turbogap}
DATA=${TURBOGAP_POLY3OP_DATA:-$repo/../turbogap_tests/GST}
F90=${TURBOGAP_F90:-mpif90}

if [ ! -x "$BIN" ]; then
  echo "SKIP: no binary at $BIN (set TURBOGAP_BIN)"
  exit 0
fi
if [ ! -d "$DATA/gap_files" ]; then
  echo "SKIP: no GST test data at $DATA (set TURBOGAP_POLY3OP_DATA)"
  exit 0
fi
if [ ! -f "$DATA/atoms_897.xyz" ]; then
  echo "SKIP: $DATA/atoms_897.xyz is missing"
  exit 0
fi

# The hypers the verifier needs, read out of the potential rather than hardcoded.
gapfile=$(ls "$DATA"/gap_files/*.gap | head -1)
get() { grep -m1 "^$1 *=" "$gapfile" | sed 's/.*= *//' | tr -d '"' | awk '{print $1}'; }
basis=$(get basis)
if [ "$basis" != "poly3operator" ]; then
  echo "SKIP: $gapfile is basis=$basis, not poly3operator"
  exit 0
fi
alpha_max=$(get n_max)
rcut_hard=$(get rcut)
buffer=$(get buffer)
atom_sigma=$(get atom_sigma_r)
amplitude_scaling=$(get amplitude_scaling)
radial_enhancement=$(get radial_enhancement)
central_weight=$(get central_weight)
rcut_soft=$(awk -v a="$rcut_hard" -v b="$buffer" 'BEGIN{printf "%.12f", a-b}')

echo "potential: $gapfile"
echo "  alpha_max=$alpha_max rcut_hard=$rcut_hard rcut_soft=$rcut_soft atom_sigma=$atom_sigma"
echo "  amplitude_scaling=$amplitude_scaling radial_enhancement=$radial_enhancement central_weight=$central_weight"
echo

echo "--- building poly3operatorverify ---"
$F90 -O2 -J "$work" -c "$repo/src/soap_turbo/src/soap_turbo_functions.f90" -o "$work/f.o"
$F90 -O2 -J "$work" -c "$repo/src/soap_turbo/src/soap_turbo_radial.f90" -o "$work/r.o"
$F90 -O2 -I "$work" "$here/poly3operatorverify.f90" "$work/f.o" "$work/r.o" \
     -o "$work/poly3operatorverify" -llapack -lblas

echo "--- running $BIN on the 897-atom GST cell ---"
run=$work/run
mkdir -p "$run"
cp "$DATA/atoms_897.xyz" "$run/atoms.xyz"
ln -s "$DATA/gap_files" "$run/gap_files"
cat > "$run/input" <<EOF
atoms_file = 'atoms.xyz'
pot_file = '$(basename "$gapfile" | sed 's|^|gap_files/|')'
n_species = 3
species = Ge Sb Te
masses = 72.64 121.76 127.6
e0 = 0.0 0.0 0.0
write_forces = .true.
write_local_energies = .true.
EOF
( cd "$run" && TURBOGAP_DUMP_RADIAL=1 "$BIN" predict > out.log 2> err.log ) || {
  echo "FAIL: turbogap predict exited non-zero"; tail -20 "$run/err.log"; exit 1; }
if [ ! -s "$run/radial_exp_coeff_dump.dat" ]; then
  echo "FAIL: no radial_exp_coeff_dump.dat was written"
  exit 1
fi

echo
"$work/poly3operatorverify" "$run/radial_exp_coeff_dump.dat" \
    "$alpha_max" "$rcut_hard" "$rcut_soft" "$atom_sigma" \
    "$amplitude_scaling" "$radial_enhancement" "$central_weight" \
    "${TURBOGAP_POLY3OP_TOL:-1e-9}"
