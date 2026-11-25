#!/usr/bin/env bash
#
# The pdf -> partial structure factor -> XRD/ND chain must not depend on the
# order the atoms are listed in.
#
# The partials are indexed by species id, so the form factor or scattering
# length of slot k has to be that of species k. It used to be taken from
# whichever atom first raised the slot count, which named the slots wrongly
# whenever the species did not first appear in `species` order: for
# `species = C O` with an O first in the file, C was scattered as O, and the
# partial structure factor files overwrote one another under repeated names.
#
# The same melted graphite-O cell as the xrd_mad case, run with its atoms in
# file order (C first) and with the O atoms moved to the front, at one rank and
# at several. Every pdf, sf, xrd and nd output must agree, and the partial
# structure factors must come out as C_C, C_O and O_O.
#
# A species can also be listed in `species` and absent from the structure, as
# when hybrid MC removes it. Its partials have no atoms to normalise by; 0/0
# there used to make every total NaN. The cell with its O removed, run as
# `species = O C` (O absent and first) and as `species = C`, must give the
# same totals and C_C partials, with no NaN anywhere -- and the same energy,
# term by term: with `species = C` the potential's O descriptors match no
# listed species, and the device build used to count them as another pair's.
#
# Environment:
#   TURBOGAP_BIN        binary under test (default: <repo>/bin/turbogap)
#   TURBOGAP_DATA_ROOT  directory holding the test systems (see tests/data_root.sh)
#   TURBOGAP_RANKS      ranks for the MPI leg (default: 3)
#   TURBOGAP_KEEP       keep the staging directory for inspection

set -u

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/../.." && pwd)

BIN=${TURBOGAP_BIN:-$repo/bin/turbogap}
# shellcheck source=../data_root.sh
. "$here/../data_root.sh"
# shellcheck source=../mpi_oversubscribe.sh
. "$here/../mpi_oversubscribe.sh"
DATA=$DATA_ROOT/xrd_mad
RANKS=${TURBOGAP_RANKS:-3}
WORK=${TMPDIR:-/tmp}/turbogap_xrd_species_order.$$

# Reordering changes only the order pairs are summed in; the outputs carry 8
# decimals.
TOL=1e-6
OUTPUTS="xrd_prediction.dat nd_prediction.dat structure_factor_total.dat pair_distribution_total.dat
structure_factor_C_C_prediction.dat structure_factor_C_O_prediction.dat structure_factor_O_O_prediction.dat
pair_distribution_C_C_prediction.dat pair_distribution_C_O_prediction.dat pair_distribution_O_O_prediction.dat"
# What a C-only run writes either way.
OUTPUTS_C_ONLY="xrd_prediction.dat nd_prediction.dat structure_factor_total.dat pair_distribution_total.dat
structure_factor_C_C_prediction.dat pair_distribution_C_C_prediction.dat"

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 2
}

species_c_o() {
  printf "n_species = 2\nspecies = C O\nmasses = 12.01 15.99\ne0 = -.16138053 0.\n"
}

species_o_c() {
  printf "n_species = 2\nspecies = O C\nmasses = 15.99 12.01\ne0 = 0. -.16138053\n"
}

species_c() {
  printf "n_species = 1\nspecies = C\nmasses = 12.01\ne0 = -.16138053\n"
}

# write_input <dir> [species function, default species_c_o]
write_input() {
  {
    echo "atoms_file = 'atoms.xyz'"
    echo "pot_file = 'gap_files/CO.gap'"
    "${2:-species_c_o}"
  } >"$1/input"
  cat >>"$1/input" <<EOF

do_pair_distribution        = .true.
pair_distribution_kde_sigma =   0.1
pair_distribution_n_samples =  201
pair_distribution_partial   = .true.
pair_distribution_rcut      =   6.0
r_range_min                 =   0.1
r_range_max                 =   6.0
write_pair_distribution     = .true.

do_structure_factor         = .true.
structure_factor_from_pdf   = .true.
structure_factor_window     = .true.
write_structure_factor      = .true.
q_range_min                 = 1.0
q_range_max                 = 10.0
q_units                     = 'q'
structure_factor_n_samples  = 51

do_xrd                      = .true.
write_xrd                   = .true.
xrd_output                  = 'q*F(q)'
do_nd                       = .true.
write_nd                    = .true.
EOF
}

# stage <dir> <xyz> [species function]
stage() {
  mkdir -p "$1"
  ln -sf "$DATA/gap_files" "$1/gap_files"
  cp "$2" "$1/atoms.xyz"
  write_input "$1" "${3:-species_c_o}"
}

# run <dir> <ranks>
run() {
  local status
  if [ "$2" -gt 1 ]; then
    (cd "$1" && mpirun $MPI_OVERSUBSCRIBE -np "$2" "$BIN" predict) >"$1/run.log" 2>&1
  else
    (cd "$1" && "$BIN" predict) >"$1/run.log" 2>&1
  fi
  status=$?
  if [ $status -ne 0 ]; then
    printf '    FAIL: turbogap exited %d in %s\n' "$status" "$1"
    tail -12 "$1/run.log" | sed 's/^/      /'
    return 1
  fi
}

# max |y_a - y_b| over two (x, y) files with the same x grid; empty on a mismatch.
max_diff() {
  awk '/^#/ || NF < 2 { next }
       FNR == NR { m++; x[m] = $1; y[m] = $2; next }
       { k++; if (k > m || x[k] != $1) bad = 1; d = $2 - y[k]; if (d < 0) d = -d; if (d > w) w = d }
       END { if (bad || k != m || m == 0) print ""; else printf "%.3e\n", w }' "$1" "$2"
}

# compare <reference dir> <test dir> [outputs, default OUTPUTS]
compare() {
  local ok=0 f d
  for f in ${3:-$OUTPUTS}; do
    if [ ! -f "$1/$f" ] || [ ! -f "$2/$f" ]; then
      printf '    FAIL: %s missing\n' "$f"
      ok=1
      continue
    fi
    d=$(max_diff "$1/$f" "$2/$f")
    if [ -z "$d" ]; then
      printf '    FAIL: %s has a different grid\n' "$f"
      ok=1
    elif awk -v d="$d" -v t="$TOL" 'BEGIN { exit !(d > t) }'; then
      printf '    FAIL: %s differs by %s\n' "$f" "$d"
      ok=1
    fi
  done
  return $ok
}

# partials_named <dir>
partials_named() {
  local got
  got=$(cd "$1" && ls structure_factor_*_prediction.dat | tr '\n' ' ')
  if [ "$got" != "structure_factor_C_C_prediction.dat structure_factor_C_O_prediction.dat structure_factor_O_O_prediction.dat " ]; then
    printf '    FAIL: partial structure factors written as: %s\n' "$got"
    return 1
  fi
}

# no_nan <dir>
no_nan() {
  local bad
  bad=$(cd "$1" && grep -il nan -- *.dat | tr '\n' ' ')
  if [ -n "$bad" ]; then
    printf '    FAIL: NaN in %s\n' "$bad"
    return 1
  fi
}

# energy_terms <dir>: "name value" for each energy term in the run's log
energy_terms() {
  awk '/ (SOAP|2b|3b|core_pot|Total) energy:/ { print $1, $(NF - 2) }' "$1/run.log"
}

# same_energies <reference dir> <dir>
same_energies() {
  energy_terms "$1" >"$1/energy_terms.txt"
  energy_terms "$2" >"$2/energy_terms.txt"
  if [ ! -s "$1/energy_terms.txt" ]; then
    printf '    FAIL: no energies in %s\n' "$1/run.log"
    return 1
  fi
  awk -v t="$TOL" 'FNR == NR { e[$1] = $2; next }
       { d = $2 - e[$1]; if (d < 0) d = -d
         if (d > t) { printf "    FAIL: %s energy %s, %s with species = C\n", $1, $2, e[$1]; bad = 1 } }
       END { exit bad }' "$1/energy_terms.txt" "$2/energy_terms.txt"
}

# absent_leg <dir> <reference dir>: O listed, none in the cell
absent_leg() {
  local ok=0
  run "$1" 1 || return 1
  no_nan "$1" || ok=1
  compare "$2" "$1" "$OUTPUTS_C_ONLY" || ok=1
  same_energies "$2" "$1" || ok=1
  return $ok
}

# reordered_leg <dir> <ranks>: both checks run, so a failure shows every symptom.
reordered_leg() {
  local ok=0
  run "$1" "$2" || return 1
  partials_named "$1" || ok=1
  compare "$WORK/c_first" "$1" || ok=1
  return $ok
}

[ -x "$BIN" ] || die "binary under test not found or not executable: $BIN"
for f in atoms.xyz gap_files/CO.gap; do
  if [ ! -e "$DATA/$f" ]; then
    printf 'SKIP: missing test data %s\n' "$DATA/$f"
    exit 0
  fi
done

mkdir -p "$WORK"
if [ -z "${TURBOGAP_KEEP:-}" ]; then
  trap 'rm -rf "$WORK"' EXIT
else
  trap 'printf "staging kept in %s\n" "$WORK"' EXIT
fi

{
  head -n 2 "$DATA/atoms.xyz"
  awk 'NR > 2 && $1 == "O"' "$DATA/atoms.xyz"
  awk 'NR > 2 && $1 != "O"' "$DATA/atoms.xyz"
} >"$WORK/o_first.xyz"
[ "$(awk 'NR == 3 { print $1 }' "$WORK/o_first.xyz")" = O ] || die "could not put an O atom first"

{
  awk 'NR > 2 && $1 != "O"' "$DATA/atoms.xyz" | wc -l
  sed -n 2p "$DATA/atoms.xyz"
  awk 'NR > 2 && $1 != "O"' "$DATA/atoms.xyz"
} >"$WORK/c_only.xyz"

pass=0
fail=0

stage "$WORK/c_first" "$DATA/atoms.xyz"
stage "$WORK/o_first_1" "$WORK/o_first.xyz"
stage "$WORK/o_first_n" "$WORK/o_first.xyz"
stage "$WORK/c_only_c" "$WORK/c_only.xyz" species_c
stage "$WORK/c_only_o_c" "$WORK/c_only.xyz" species_o_c

printf '== reference: atoms in file order (C first), 1 rank ==\n'
if run "$WORK/c_first" 1 && partials_named "$WORK/c_first"; then
  printf '    PASS\n'
  pass=$((pass + 1))
else
  fail=$((fail + 1))
  printf 'xrd_species_order: no reference, nothing compared\n'
  exit 1
fi

printf '== O atoms first, 1 rank: every output as in file order ==\n'
if reordered_leg "$WORK/o_first_1" 1; then
  printf '    PASS\n'
  pass=$((pass + 1))
else
  fail=$((fail + 1))
fi

printf '== O atoms first, %d ranks: every output as in file order ==\n' "$RANKS"
if reordered_leg "$WORK/o_first_n" "$RANKS"; then
  printf '    PASS\n'
  pass=$((pass + 1))
else
  fail=$((fail + 1))
fi

printf '== O listed but absent (species = O C), against species = C ==\n'
if run "$WORK/c_only_c" 1 && no_nan "$WORK/c_only_c" && absent_leg "$WORK/c_only_o_c" "$WORK/c_only_c"; then
  printf '    PASS\n'
  pass=$((pass + 1))
else
  fail=$((fail + 1))
fi

printf '\nxrd_species_order: %d passed, %d failed\n' "$pass" "$fail"
[ "$fail" -eq 0 ]
