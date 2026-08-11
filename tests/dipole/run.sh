#!/usr/bin/env bash
#
# Dipole acceptance test: does TurboGAP reproduce QUIP's dipole prediction?
#
# This is deliberately NOT part of tests/regression. That suite compares
# byte-for-byte against baseline/turbogap.e6eb1aa, a binary frozen before the
# refactor which knows nothing about dipole models and cannot produce a
# reference here at all; and its comparator is plain `diff`, which has no
# tolerance. This test asks a different question -- do two independent codes
# agree on the same model -- and so needs a tolerance and a real comparator.
#
# The reference is QUIP's own output, checked in as produced by
#   quip atoms_filename=water_dipole_tnep_converted_test.xyz \
#        param_filename=water_dipole.xml e \
#        calc_args="dipole=dipole local_dipole=local_dipole" > out_quip_predict
# on 200 frames of 2-water clusters.
#
# Environment:
#   TURBOGAP_BIN        binary under test (default: <repo>/bin/turbogap)
#   TURBOGAP_DATA_ROOT  directory holding the test systems
#                       (default: $HOME/work/cpu_vs_gpu_tests/input)
#   TURBOGAP_DIPOLE_TOL absolute tolerance (default: 1e-6)
#   TURBOGAP_RANKS      MPI ranks to run with (default: 1)
#   TURBOGAP_KEEP       keep the staging directory for inspection
#
# Usage: run.sh

set -u

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/../.." && pwd)

BIN=${TURBOGAP_BIN:-$repo/bin/turbogap}
DATA_ROOT=${TURBOGAP_DATA_ROOT:-$HOME/work/cpu_vs_gpu_tests/input}
DATA=$DATA_ROOT/water_dipole
TOL=${TURBOGAP_DIPOLE_TOL:-1e-6}
RANKS=${TURBOGAP_RANKS:-1}
WORK=${TMPDIR:-/tmp}/turbogap_dipole.$$

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 2
}

[ -x "$BIN" ] || die "binary under test not found or not executable: $BIN"

# A case whose data is absent is skipped, not failed -- same contract as
# tests/regression/run.sh.
for f in water_dipole.xml out_quip_predict water_dipole_tnep_converted_test.xyz; do
  if [ ! -e "$DATA/$f" ]; then
    printf 'SKIP: missing test data %s\n' "$DATA/$f"
    printf '      (set TURBOGAP_DATA_ROOT, or see tests/dipole/README.md)\n'
    exit 0
  fi
done

if [ "$RANKS" -gt 1 ] && ! command -v mpirun >/dev/null; then
  printf 'SKIP: mpirun not available (asked for %d ranks)\n' "$RANKS"
  exit 0
fi

mkdir -p "$WORK" || die "cannot create $WORK"
if [ -z "${TURBOGAP_KEEP:-}" ]; then
  trap 'rm -rf "$WORK"' EXIT
else
  trap 'printf "staging kept in %s\n" "$WORK"' EXIT
fi

printf 'binary   %s\n' "$BIN"
printf 'data     %s\n' "$DATA"
printf 'ranks    %s\n\n' "$RANKS"

# The model ships as a QUIP xml. Convert it on each run rather than relying on a
# derived .gap that could drift from the xml it came from -- but the converter
# needs BeautifulSoup, which is not everywhere, so fall back to a pre-converted
# gap_files/ beside the data when it is missing.
( cd "$WORK" && cp "$DATA/water_dipole.xml" . && cp "$DATA"/water_dipole.xml.sparseX.* . ) \
  || die "could not stage the model"

if python3 -c 'import bs4' >/dev/null 2>&1; then
  ( cd "$WORK" && python3 "$repo/tools/quip_xml_to_gap/make_gap_files.py" \
      water_dipole.xml water_dipole.gap --dipole-model ) >"$WORK/convert.log" 2>&1 \
    || { sed 's/^/    /' "$WORK/convert.log"; die "xml -> gap conversion failed"; }
elif [ -d "$DATA/gap_files" ]; then
  printf 'NOTE: python bs4 not available; using pre-converted %s/gap_files\n\n' "$DATA"
  cp -r "$DATA/gap_files" "$WORK/gap_files" || die "could not copy pre-converted gap_files"
else
  printf 'SKIP: python bs4 is not installed and %s/gap_files does not exist,\n' "$DATA"
  printf '      so the model cannot be converted. pip install beautifulsoup4 lxml\n'
  exit 0
fi

grep -q 'dipole_model = .true.' "$WORK/gap_files/water_dipole.gap" \
  || die "the .gap in use has no dipole_model flag"

cp "$DATA/water_dipole_tnep_converted_test.xyz" "$WORK/atoms.xyz"
cp "$here/input" "$WORK/input"

python3 "$here/make_reference.py" "$DATA/out_quip_predict" "$WORK/reference.dat" \
  || die "could not extract the QUIP reference"

if [ "$RANKS" -gt 1 ]; then
  ( cd "$WORK" && mpirun -np "$RANKS" "$BIN" predict ) >"$WORK/run.log" 2>&1
else
  ( cd "$WORK" && "$BIN" predict ) >"$WORK/run.log" 2>&1
fi
rc=$?
if [ $rc -ne 0 ]; then
  printf 'turbogap failed (exit %d):\n' "$rc"
  tail -20 "$WORK/run.log" | sed 's/^/    /'
  exit 1
fi

python3 "$here/compare_dipole.py" "$WORK/reference.dat" "$WORK/trajectory_out.xyz" "$TOL"
