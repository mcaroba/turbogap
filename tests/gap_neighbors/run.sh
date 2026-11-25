#!/usr/bin/env bash
#
# The GAP terms walk only the pairs within their own cutoffs.
#
# A pdf, XRD, ND or vdW term widens the neighbour list far past the
# potential's cutoffs. The GAP terms then get the pairs of it within theirs:
# SOAP within the largest GAP cutoff, and the 2b, 3b and core-potential terms
# each within their own. Each site keeps its neighbours in list order, so the
# sums run over the same pairs in the same order as on the whole list.
#
# The xrd_mad cell with the CO potential (SOAP and 2b at 4.5 A, 3b at 2 A,
# core potential to 1 A) is predicted
#   - with the GAP alone: the 3b and core-potential lists are cut;
#   - with a pdf to 9 A on top, writing no forces, at one rank and at several:
#     the SOAP and 2b lists are cut too;
#   - with a potential of its 3b descriptors alone.
# expected/ holds what master 8a7df2b, which walked the whole list in every
# term, wrote for the first and the last. Every number must agree with it to
# TOL; the pdf moves the last printed digit even on master, so not bit for bit.
# A cut that drops one pair moves the energies by orders of magnitude more.
#
# Environment:
#   TURBOGAP_BIN        binary under test (default: <repo>/bin/turbogap)
#   TURBOGAP_DATA_ROOT  directory holding the test systems (see tests/data_root.sh)
#   TURBOGAP_RANKS      ranks for the MPI leg (default: 3)
#   TURBOGAP_KEEP       keep the staging directory for inspection
#   TURBOGAP_BLESS      write expected/ from TURBOGAP_BIN instead of checking;
#                       only from a build that does not cut the lists

set -u

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/../.." && pwd)

BIN=${TURBOGAP_BIN:-$repo/bin/turbogap}
# shellcheck source=../data_root.sh
. "$here/../data_root.sh"
# shellcheck source=../mpi_oversubscribe.sh
. "$here/../mpi_oversubscribe.sh"
DATA=$DATA_ROOT/xrd_mad
EXPECTED=$here/expected
RANKS=${TURBOGAP_RANKS:-3}
WORK=${TMPDIR:-/tmp}/turbogap_gap_neighbors.$$
PDF_RCUT=9.0
TOL=1e-6

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 2
}

gap_deck() {
  cat <<EOF
atoms_file = 'atoms.xyz'
pot_file = 'gap_files/CO.gap'
n_species = 2
species = C O
masses = 12.01 15.99
e0 = -.16138053 0.
write_xyz = 1
EOF
}

pdf_deck() {
  gap_deck
  cat <<EOF
do_pair_distribution = .true.
pair_distribution_n_samples = 201
pair_distribution_kde_sigma = 0.1
pair_distribution_rcut = $PDF_RCUT
r_range_min = 0.1
r_range_max = $PDF_RCUT
EOF
}

three_body_deck() {
  gap_deck | sed "s#gap_files/CO.gap#three_body.gap#"
}

# The angle_3b blocks of CO.gap; their paths stay relative to gap_files.
three_body_gap() {
  awk '/^gap_beg/ { keep = ($2 == "angle_3b") } keep { print } /^gap_end/ { keep = 0 }' "$DATA/gap_files/CO.gap"
}

# run <dir> <deck function> <ranks>
run() {
  local status
  mkdir -p "$1"
  ln -sf "$DATA/gap_files" "$1/gap_files"
  cp "$DATA/atoms.xyz" "$1/atoms.xyz"
  three_body_gap >"$1/three_body.gap"
  "$2" >"$1/input"
  if [ "$3" -gt 1 ]; then
    (cd "$1" && mpirun $MPI_OVERSUBSCRIBE -np "$3" "$BIN" predict) >"$1/run.log" 2>&1
  else
    (cd "$1" && "$BIN" predict) >"$1/run.log" 2>&1
  fi
  status=$?
  if [ $status -ne 0 ]; then
    printf '    FAIL: turbogap exited %d in %s\n' "$status" "$1"
    tail -12 "$1/run.log" | sed 's/^/      /'
    return 1
  fi
  if [ ! -s "$1/trajectory_out.xyz" ]; then
    printf '    FAIL: no trajectory_out.xyz in %s\n' "$1"
    return 1
  fi
}

# max_diff <a> <b>: largest difference of the numbers two extended-xyz files
# write in the same places; "layout" if any other word differs.
max_diff() {
  awk '{ gsub(/[="]/, " ") }
       FNR == NR { line[FNR] = $0; next }
       { n = split(line[FNR], a, " "); m = split($0, b, " ")
         if (n != m) { bad = 1; exit }
         for (i = 1; i <= n; i++) {
           if (a[i] ~ /^-?[0-9.]+([eEdD][-+]?[0-9]+)?$/) { d = a[i] - b[i]; if (d < 0) d = -d; if (d > w) w = d }
           else if (a[i] != b[i]) bad = 1
         } }
       END { if (bad) print "layout"; else printf "%.3e\n", w }' "$1" "$2"
}

# leg <name> <deck function> <ranks> <expected file>
leg() {
  local d
  printf '== %s, %s rank(s), against %s ==\n' "$1" "$3" "$(basename "$4")"
  run "$WORK/$1_$3" "$2" "$3" || return 1
  d=$(max_diff "$4" "$WORK/$1_$3/trajectory_out.xyz")
  if [ "$d" = layout ] || awk -v d="$d" -v t="$TOL" 'BEGIN { exit !(d > t) }'; then
    printf '    FAIL: differs by %s\n' "$d"
    return 1
  fi
  printf '    largest difference %s\n' "$d"
}

# bless: expected/ from BIN, which must be a build that walks the whole list.
bless() {
  mkdir -p "$EXPECTED"
  run "$WORK/gap" gap_deck 1 || exit 1
  run "$WORK/three_body" three_body_deck 1 || exit 1
  cp "$WORK/gap/trajectory_out.xyz" "$EXPECTED/gap.xyz"
  cp "$WORK/three_body/trajectory_out.xyz" "$EXPECTED/three_body.xyz"
  printf 'wrote %s/gap.xyz and three_body.xyz from %s\n' "$EXPECTED" "$BIN"
  exit 0
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

[ -n "${TURBOGAP_BLESS:-}" ] && bless
[ -f "$EXPECTED/gap.xyz" ] && [ -f "$EXPECTED/three_body.xyz" ] || die "no expected/ (see TURBOGAP_BLESS)"

pass=0
fail=0
for args in "gap gap_deck 1 gap.xyz" "pdf pdf_deck 1 gap.xyz" "pdf pdf_deck $RANKS gap.xyz" \
  "three_body three_body_deck 1 three_body.xyz"; do
  read -r name deck ranks expected <<<"$args"
  if leg "$name" "$deck" "$ranks" "$EXPECTED/$expected"; then
    printf '    PASS\n'
    pass=$((pass + 1))
  else
    fail=$((fail + 1))
  fi
done

printf '\npassed: %d   failed: %d\n' "$pass" "$fail"
[ "$fail" -eq 0 ]
