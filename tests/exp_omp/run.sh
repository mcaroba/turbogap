#!/usr/bin/env bash
#
# The pdf and the pdf, structure-factor and diffraction forces are threaded with
# OpenMP on the host. Each pair's kernel and force are formed in parallel and
# summed in pair order afterwards, so the thread count must not change them.
#
# The 897-atom carbon cell of the mad_sf_matrix case, with enough pairs to cross
# several of the pdf's scratch chunks, at 1 thread, at 4 and at 2 ranks x 2:
#
#   predict   every pdf, sf, xrd and nd file: identical across thread counts,
#             and within 1e-7 across ranks, which reduce in a different order
#   xrd MAD   step-0 forces and virial with the xrd forces on
#   pdf MAD   the same with a pdf target
#
# The rest of the code is not bit-reproducible across thread counts, so the MAD
# legs compare to 1e-6: a race in the threaded loops is an error of order one.
#
# Environment:
#   TURBOGAP_BIN        binary under test (default: <repo>/bin/turbogap)
#   TURBOGAP_DATA_ROOT  directory holding the test systems (see tests/data_root.sh)
#   TURBOGAP_KEEP       keep the staging directory for inspection

set -u

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/../.." && pwd)

BIN=${TURBOGAP_BIN:-$repo/bin/turbogap}
# shellcheck source=../data_root.sh
. "$here/../data_root.sh"
# shellcheck source=../mpi_oversubscribe.sh
. "$here/../mpi_oversubscribe.sh"
DATA=$DATA_ROOT/C
XRD_TARGET=xrd_ta-C_Shang_2021_27_GPa_1000C_hand.fq
WORK=${TMPDIR:-/tmp}/turbogap_exp_omp.$$

RANK_TOL=1e-7
FORCE_TOL=1e-6
LAYOUTS="1x4 2x2"

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 2
}

species() {
  cat <<EOF
atoms_file = 'atoms.xyz'
pot_file = 'gap_files/carbon.gap'
n_species = 1
species = C
masses = 12.01
e0 = -.16138053
EOF
}

observables() {
  cat <<EOF
do_pair_distribution = .true.
write_pair_distribution = .true.
pair_distribution_partial = .true.
pair_distribution_n_samples = 301
pair_distribution_kde_sigma = 0.1
r_range_min = 0.1
r_range_max = 12.0
pair_distribution_rcut = 12.6
do_structure_factor = .true.
write_structure_factor = .true.
structure_factor_from_pdf = .true.
structure_factor_matrix = .true.
structure_factor_matrix_forces = .true.
structure_factor_window = .true.
structure_factor_n_samples = 301
q_units = "q"
q_range_min = 0.1
q_range_max = 10.0
EOF
}

deck_predict() {
  species
  observables
  cat <<EOF
do_xrd = .true.
write_xrd = .true.
xrd_output = 'q*i(q)'
do_nd = .true.
write_nd = .true.
nd_output = 'q*F(q)'
EOF
}

# One NVE step: the step-0 frame carries the forces and the virial.
deck_mad() {
  species
  observables
  cat <<EOF
md_nsteps = 1
md_step = 1.0
write_xyz = 1
n_exp = 1
exp_energies = .true.
exp_forces = .true.
exp_similarity_type = 'squared_diff'
exp_energy_scales_beg = 1000
exp_energy_scales_end = 1000
EOF
  if [ "$1" = xrd ]; then
    cat <<EOF
do_xrd = .true.
write_xrd = .true.
xrd_output = 'q*i(q)'
exp_labels = 'xrd'
exp_data_files = '$XRD_TARGET'
exp_input_type = 'q*i(q)'
EOF
  else
    cat <<EOF
exp_labels = 'pair_distribution'
exp_data_files = 'target_pdf.dat'
EOF
  fi
}

# run <dir> <layout RxT> <mode> <deck text>
run() {
  local dir=$1 ranks=${2%x*} threads=${2#*x} mode=$3
  mkdir -p "$dir"
  ln -sf "$DATA/gap_files" "$dir/gap_files"
  cp "$DATA/atoms_897.xyz" "$dir/atoms.xyz"
  cp "$DATA/$XRD_TARGET" "$dir/"
  [ -f "$WORK/target_pdf.dat" ] && cp "$WORK/target_pdf.dat" "$dir/"
  printf '%s\n' "$4" >"$dir/input"
  if [ "$ranks" -gt 1 ]; then
    (cd "$dir" && OMP_NUM_THREADS=$threads mpirun $MPI_OVERSUBSCRIBE -np "$ranks" "$BIN" "$mode") >"$dir/run.log" 2>&1
  else
    (cd "$dir" && OMP_NUM_THREADS=$threads "$BIN" "$mode") >"$dir/run.log" 2>&1
  fi
  local status=$?
  if [ $status -ne 0 ]; then
    printf '    FAIL: turbogap exited %d in %s\n' "$status" "$dir"
    tail -8 "$dir/run.log" | sed 's/^/      /'
    return 1
  fi
}

# max |y_a - y_b| over two (x, y) files on one grid; empty if the grids differ.
max_diff() {
  awk '/^#/ || NF < 2 { next }
       FNR == NR { m++; x[m] = $1; y[m] = $2; next }
       { k++; if (k > m || x[k] != $1) bad = 1; d = $2 - y[k]; if (d < 0) d = -d; if (d > w) w = d }
       END { if (bad || k != m || m == 0) print ""; else printf "%.3e\n", w }' "$1" "$2"
}

# max |difference| over the forces and the virial of the first frame of two xyz files.
frame0_diff() {
  awk 'function cols(line,   p, n, i, f, s) {
         match(line, /Properties=[^ ]*/); p = substr(line, RSTART + 11, RLENGTH - 11)
         n = split(p, f, ":"); s = 1
         for (i = 1; i <= n; i += 3) { if (f[i] == "forces") return s; s += f[i + 2] }
         return 0 }
       function virial(line, v,   t) {
         match(line, /virial="[^"]*"/); t = substr(line, RSTART + 8, RLENGTH - 9); return split(t, v, " ") }
       FNR == 1 { nat = $1; file++ }
       FNR == 2 { if (file == 1) { fc = cols($0); nv = virial($0, va) } else { virial($0, vb) } next }
       FNR > 2 && FNR <= nat + 2 {
         for (c = 0; c < 3; c++) { if (file == 1) a[FNR, c] = $(fc + c); else { d = $(fc + c) - a[FNR, c]; if (d < 0) d = -d; if (d > w) w = d } } }
       END { for (i = 1; i <= nv; i++) { d = va[i] - vb[i]; if (d < 0) d = -d; if (d > w) w = d }
             if (fc == 0) print ""; else printf "%.3e\n", w }' "$1" "$2"
}

over() {
  awk -v d="$1" -v t="$2" 'BEGIN { exit !(d > t) }'
}

[ -x "$BIN" ] || die "binary under test not found or not executable: $BIN"
for f in atoms_897.xyz gap_files/carbon.gap $XRD_TARGET; do
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

pass=0
fail=0

printf '== predict: pdf, sf, xrd and nd files ==\n'
if run "$WORK/predict_1x1" 1x1 predict "$(deck_predict)"; then
  outputs=$(cd "$WORK/predict_1x1" && ls pair_distribution_*.dat structure_factor_*.dat xrd_prediction.dat nd_prediction.dat)
  awk 'NF == 2 && $1 !~ /^#/ { print $1, 0.9 * $2 + 0.1 }' "$WORK/predict_1x1/pair_distribution_total.dat" \
    | head -n 301 >"$WORK/target_pdf.dat"
  for layout in $LAYOUTS; do
    ok=0
    run "$WORK/predict_$layout" "$layout" predict "$(deck_predict)" || ok=1
    for f in $outputs; do
      if [ "${layout%x*}" = 1 ]; then
        cmp -s "$WORK/predict_1x1/$f" "$WORK/predict_$layout/$f" || { printf '    FAIL: %s differs at %s\n' "$f" "$layout"; ok=1; }
      else
        d=$(max_diff "$WORK/predict_1x1/$f" "$WORK/predict_$layout/$f")
        if [ -z "$d" ] || over "$d" $RANK_TOL; then
          printf '    FAIL: %s differs by %s at %s\n' "$f" "${d:-grid}" "$layout"
          ok=1
        fi
      fi
    done
    if [ $ok = 0 ]; then
      printf '    PASS %s (%d files)\n' "$layout" "$(echo $outputs | wc -w)"
      pass=$((pass + 1))
    else
      fail=$((fail + 1))
    fi
  done
else
  fail=$((fail + 1))
fi

for kind in xrd pdf; do
  printf '== %s MAD: step-0 forces and virial ==\n' "$kind"
  if ! run "$WORK/${kind}_1x1" 1x1 md "$(deck_mad $kind)"; then
    fail=$((fail + 1))
    continue
  fi
  for layout in $LAYOUTS; do
    if run "$WORK/${kind}_$layout" "$layout" md "$(deck_mad $kind)"; then
      d=$(frame0_diff "$WORK/${kind}_1x1/trajectory_out.xyz" "$WORK/${kind}_$layout/trajectory_out.xyz")
      if [ -n "$d" ] && ! over "$d" $FORCE_TOL; then
        printf '    PASS %s (max |d| %s)\n' "$layout" "$d"
        pass=$((pass + 1))
      else
        printf '    FAIL %s: forces or virial differ by %s\n' "$layout" "${d:-unreadable}"
        fail=$((fail + 1))
      fi
    else
      fail=$((fail + 1))
    fi
  done
done

printf '\nexp_omp: %d passed, %d failed\n' "$pass" "$fail"
[ "$fail" -eq 0 ]
