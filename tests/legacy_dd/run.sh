#!/usr/bin/env bash
#
# Does dd_type = legacy give the same energies, forces, virial and local
# properties as the replicated run of the same input?
#
# Each leg runs one input twice, replicated and decomposed, and compares
# trajectory_out.xyz number by number. The domains sum each atom's terms in a
# different order, so the comparison is a tolerance, not a diff.
#
# The legs are chosen for the ways a decomposition goes wrong:
#   - 2 domains along an axis: both neighbours are the same domain, so a halo
#     atom arrives twice and must be kept once;
#   - 3 along an axis: the two neighbours differ;
#   - 2x2x2: edge and corner atoms reach a domain through two forwards, and
#     their forces have to come back the same way;
#   - several ranks per domain, block and cyclic, so the halo forces of
#     non-root ranks have to reach the root;
#   - a sheared cell, where a face's normal is not a cell vector;
#   - the grid chosen automatically, and equal-volume cuts;
#   - a slab with vacuum, where count balancing moves the cuts, and one cut
#     into equal slabs, two of which hold no atom at all, owned or halo;
#   - MD: five steps; twenty, which reuse the domains between rebuilds; and
#     twenty with a tiny skin, which rebuild nearly every step;
#   - dd_integrate = domain: Berendsen, Bussi and NVE, a tiny skin (atoms
#     migrate at every step), empty domains, and a thermostat it does not
#     cover, which must integrate on rank 0 instead and say so;
#   - TS van der Waals, which reads its neighbours' Hirshfeld volumes, so the
#     halo copies have to carry them;
#   - electrostatics, whose neutralising mean charge is a sum over all atoms
#     (on a 24k cell: this potential's rcut_max is 20 A);
#   - the pdf -> S(q) -> XRD chain and the Debye XRD, computed on the
#     whole-structure view while the potential stays decomposed;
#   - the dipole model, whose total dipole is a sum over all atoms, and the
#     MAD-IR bias it drives through 40 MD steps;
#   - a grid too fine for the cell, which must fall back to replicated.
#
# Under dd_type rcb (TURBOGAP_DD_TYPE=rcb, as tests/rcb_dd/run.sh sets it)
# the same legs run, a dd_grid giving only the number of domains, without the
# fallback leg, since no cell is too narrow for it, and with legs of their
# own: a cluster in vacuum, a pore, prime rank counts, and a droplet in thin
# vapour integrated on the domains, where count and cost disagree and the boxes
# must be redrawn by pair count.
#
# Environment:
#   TURBOGAP_BIN        binary under test (default: <repo>/bin/turbogap)
#   TURBOGAP_DATA_ROOT  directory holding the test systems (see ../data_root.sh)
#   TURBOGAP_PYTHON     interpreter for the comparison
#   TURBOGAP_KEEP       keep the staging directory for inspection
#   TURBOGAP_RTOL, TURBOGAP_ATOL   comparison tolerances (default 1e-7 each)
#   TURBOGAP_MD_ATOL    absolute tolerance of the MD legs (default 2e-6)
#   TURBOGAP_DD_TYPE    legacy (default) or rcb
#
# Usage: run.sh [leg ...]      (no argument runs all of them)

set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/../.." && pwd)

BIN=${TURBOGAP_BIN:-$repo/bin/turbogap}
PYTHON=${TURBOGAP_PYTHON:-python3}
# shellcheck source=../data_root.sh
. "$here/../data_root.sh"
# shellcheck source=../mpi_oversubscribe.sh
. "$here/../mpi_oversubscribe.sh"
DATA=$DATA_ROOT/CO
WORK=${TMPDIR:-/tmp}/turbogap_legacy_dd.$$
COMPARE=$repo/tests/regression/compare_tol.py

# Per-atom values are written to 8-10 significant digits; a reordered sum moves
# the last of them. A device binary does not reproduce itself even replicated,
# so it needs its own, looser, tolerance.
DD_TYPE=${TURBOGAP_DD_TYPE:-legacy}
RTOL=${TURBOGAP_RTOL:-1e-7}
ATOL=${TURBOGAP_ATOL:-1e-7}
# An MD frame's virial has components of 1e3-1e4 eV carrying ~1e-6 eV of
# round-off, so a small off-diagonal one fails a relative test. Control: the
# replicated run alone moves them by up to 1.2e-6 eV between 2 and 4 ranks
# (integ_bussi's input), 4e-7 between 4 and 8 (md20_skin's) and 2e-7 between 3
# and 6 (integ_empty's).
MD_ATOL=${TURBOGAP_MD_ATOL:-2e-6}

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 2
}

check_setup() {
  [ -x "$BIN" ] || die "binary under test not found or not executable: $BIN"
  command -v mpirun >/dev/null || die "mpirun not found; dd_type legacy needs an MPI build"
  command -v "$PYTHON" >/dev/null || die "python interpreter not found: $PYTHON"
  if [ ! -e "$DATA/gap_files/CO.gap" ]; then
    printf 'SKIP: missing test data %s\n' "$DATA/gap_files/CO.gap"
    exit 0
  fi
}

# make_slab <out> <factor>: the 897-atom cell with c stretched by factor, the
# rest of it vacuum.
make_slab() {
  "$PYTHON" - "$DATA/atoms_897.xyz" "$1" "$2" <<'EOF'
import re
import sys

src, dst, factor = sys.argv[1], sys.argv[2], float(sys.argv[3])
lines = open(src).read().splitlines()
old = re.search(r'Lattice="([^"]*)"', lines[1]).group(1)
lat = [float(x) for x in old.split()]
lat[6:9] = [factor * x for x in lat[6:9]]
lines[1] = lines[1].replace(old, " ".join(f"{x:.10f}" for x in lat))
open(dst, "w").write("\n".join(lines) + "\n")
EOF
}

# make_sphere <out> <keep> <radius>: the 7176-atom cell with velocities, keeping
# the atoms inside (keep = in) or outside (out) a sphere of radius times the
# shortest cell edge about the centre, or (vapour) all those inside a sphere a
# quarter of the way along a and every eighth of the rest: a droplet to one
# side, whose vapour atoms have few neighbours.
make_sphere() {
  "$PYTHON" - "$DATA/atoms_7176_vel.xyz" "$1" "$2" "$3" <<'EOF'
import re
import sys

src, dst, keep, radius = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
lines = open(src).read().splitlines()
n = int(lines[0])
lat = [float(x) for x in re.search(r'Lattice="([^"]*)"', lines[1]).group(1).split()]
box = [lat[0], lat[4], lat[8]]
r = radius * min(box)
centre = [b / 2 for b in box]
if keep == "vapour":
    centre[0] = box[0] / 4
out = []
n_outside = 0
for line in lines[2:2 + n]:
    x = [float(v) for v in line.split()[1:4]]
    inside = sum(((x[k] - centre[k] + box[k] / 2) % box[k] - box[k] / 2) ** 2 for k in range(3)) < r * r
    if not inside:
        n_outside += 1
    if keep == "vapour":
        if inside or n_outside % 8 == 0:
            out.append(line)
    elif inside == (keep == "in"):
        out.append(line)
open(dst, "w").write(f"{len(out)}\n{lines[1]}\n" + "\n".join(out) + "\n")
EOF
}

# The 897-atom cell with c tilted by 0.3 a, atoms carried along with the cell.
make_sheared() {
  "$PYTHON" - "$DATA/atoms_897.xyz" "$1" <<'EOF'
import re
import sys

src, dst = sys.argv[1], sys.argv[2]
lines = open(src).read().splitlines()
n = int(lines[0])
lat = [float(x) for x in re.search(r'Lattice="([^"]*)"', lines[1]).group(1).split()]
a, b, c = lat[0:3], lat[3:6], lat[6:9]
c2 = [c[k] + 0.3 * a[k] for k in range(3)]
out = [lines[0], lines[1].replace(re.search(r'Lattice="([^"]*)"', lines[1]).group(1),
                                   " ".join(f"{x:.10f}" for x in a + b + c2))]
for line in lines[2:2 + n]:
    f = line.split()
    r = [float(x) for x in f[1:4]]
    frac_c = r[2] / c[2]
    r2 = [r[k] + 0.3 * a[k] * frac_c for k in range(3)]
    out.append(" ".join([f[0]] + [f"{x:.10f}" for x in r2] + f[4:]))
open(dst, "w").write("\n".join(out) + "\n")
EOF
}

# write_input <dir> <md steps, 0 for a single point> [extra lines]
write_input() {
  local dir=$1 steps=$2 extra=${3:-}
  cat >"$dir/input" <<EOF
atoms_file = "atoms.xyz"
pot_file = "gap_files/CO.gap"
n_species = 2
species = C O
masses = 12.01 15.99
e0 = -.16138053 0.
EOF
  if [ "$steps" -gt 0 ]; then
    printf 'md_nsteps = %s\nthermostat = berendsen\nwrite_xyz = 1\nrandom_seed = 12345\n' "$steps" >>"$dir/input"
  fi
  if [ -n "$extra" ]; then
    printf '%s\n' "$extra" >>"$dir/input"
  fi
}

# run_one <dir> <mode> <ranks>
run_one() {
  local dir=$1 mode=$2 ranks=$3
  (cd "$dir" && OMP_NUM_THREADS=1 mpirun $MPI_OVERSUBSCRIBE -np "$ranks" "$BIN" "$mode" >out.log 2>&1)
}

# leg <name> <atoms file> <md steps> <ranks> <dd_grid, or auto> [affinity] [extra lines for both runs]
#     [extra lines for the decomposed run only]
leg() {
  local name=$1 atoms=$2 steps=$3 ranks=$4 grid=$5 affinity=${6:-block} extra=${7:-} dd_extra=${8:-}
  local mode=predict grid_line="" atol=$ATOL
  [ "$steps" -gt 0 ] && mode=md && atol=$MD_ATOL
  [ "$grid" != auto ] && grid_line="dd_grid = $grid"
  printf '== %s: %s ranks, grid %s, %s ==\n' "$name" "$ranks" "$grid" "$affinity"

  local ok=1 kind
  last_ok=0
  for kind in ref dd; do
    mkdir -p "$WORK/$name/$kind"
    ln -sf "$DATA/gap_files" "$WORK/$name/$kind/gap_files"
    cp "$atoms" "$WORK/$name/$kind/atoms.xyz"
  done
  write_input "$WORK/$name/ref" "$steps" "$extra"
  write_input "$WORK/$name/dd" "$steps" "$extra
dd_type = $DD_TYPE
$grid_line
dd_grid_affinity = $affinity
$dd_extra"

  run_one "$WORK/$name/ref" "$mode" "$ranks" || ok=0
  run_one "$WORK/$name/dd" "$mode" "$ranks" || ok=0
  # rebuilds is reported only by a run that stayed decomposed to the end.
  if [ "$ok" -eq 1 ] && grep -q 'rebuilds,' "$WORK/$name/dd/out.log" &&
     "$PYTHON" "$COMPARE" "$WORK/$name/ref/trajectory_out.xyz" "$WORK/$name/dd/trajectory_out.xyz" \
       --rtol "$RTOL" --atol "$atol"; then
    printf '    PASS  %s\n' "$(grep -a 'rebuilds,' "$WORK/$name/dd/out.log" | sed 's/^ *dd_type [a-z]*: //')"
    pass=$((pass + 1))
    last_ok=1
  else
    printf '    FAIL (logs in %s)\n' "$WORK/$name"
    tail -5 "$WORK/$name/dd/out.log" | sed 's/^/      /'
    fail=$((fail + 1))
  fi
}

# full_leg <name> <mode> <ranks> <atoms> <gap dir> <extra file or -> <input>:
# like leg, for an input that is not the CO one; the input names its own
# potential, which is linked in as gap_files.
full_leg() {
  local name=$1 mode=$2 ranks=$3 atoms=$4 gaps=$5 extra_file=$6 body=$7
  printf '== %s: %s ranks, grid auto ==\n' "$name" "$ranks"
  local ok=1 kind atol=$ATOL
  [ "$mode" = md ] && atol=$MD_ATOL
  last_ok=0
  for kind in ref dd; do
    mkdir -p "$WORK/$name/$kind"
    ln -sfn "$gaps" "$WORK/$name/$kind/gap_files"
    cp "$atoms" "$WORK/$name/$kind/atoms.xyz"
    [ "$extra_file" != - ] && ln -sf "$extra_file" "$WORK/$name/$kind/"
    printf '%s\n' "$body" >"$WORK/$name/$kind/input"
  done
  printf 'dd_type = %s\n' "$DD_TYPE" >>"$WORK/$name/dd/input"
  run_one "$WORK/$name/ref" "$mode" "$ranks" || ok=0
  run_one "$WORK/$name/dd" "$mode" "$ranks" || ok=0
  if [ "$ok" -eq 1 ] && grep -q 'rebuilds,' "$WORK/$name/dd/out.log" &&
     "$PYTHON" "$COMPARE" "$WORK/$name/ref/trajectory_out.xyz" "$WORK/$name/dd/trajectory_out.xyz" \
       --rtol "$RTOL" --atol "$atol"; then
    printf '    PASS  %s\n' "$(grep -a 'rebuilds,' "$WORK/$name/dd/out.log" | sed 's/^ *dd_type [a-z]*: //')"
    pass=$((pass + 1))
    last_ok=1
  else
    printf '    FAIL (logs in %s)\n' "$WORK/$name"
    tail -5 "$WORK/$name/dd/out.log" | sed 's/^/      /'
    fail=$((fail + 1))
  fi
}

# tile_water <out> <nx> <ny> <nz>: the first water frame of the dipole data,
# repeated periodically. Its molecules are far apart, which is fine: the point
# is a cell large enough to split, with the dipole model on.
tile_water() {
  "$PYTHON" - "$DATA_ROOT/water_dipole/water_dipole_tnep_converted_test.xyz" "$1" "$2" "$3" "$4" <<'EOF'
import re
import sys

src, dst = sys.argv[1], sys.argv[2]
reps = [int(x) for x in sys.argv[3:6]]
lines = open(src).read().splitlines()
n = int(lines[0])
lat = [float(x) for x in re.search(r'Lattice="([^"]*)"', lines[1]).group(1).split()]
vec = [lat[0:3], lat[3:6], lat[6:9]]
atoms = [l.split() for l in lines[2:2 + n]]
out = []
for i in range(reps[0]):
    for j in range(reps[1]):
        for k in range(reps[2]):
            shift = [i * vec[0][d] + j * vec[1][d] + k * vec[2][d] for d in range(3)]
            for a in atoms:
                out.append(f"{a[0]} " + " ".join(f"{float(a[1 + d]) + shift[d]:.10f}" for d in range(3)))
big = [v * r for row, r in zip(vec, reps) for v in row]
head = f'Lattice="{" ".join(f"{x:.10f}" for x in big)}" Properties=species:S:1:pos:R:3 pbc="T T T"'
open(dst, "w").write(f"{len(out)}\n{head}\n" + "\n".join(out) + "\n")
EOF
}

# A two-peak IR spectrum, as tests/mad_ir/md_run.sh makes one.
make_ir_exp() {
  "$PYTHON" - "$1" <<'EOF'
import math
import sys

with open(sys.argv[1], "w") as f:
    nu = 0.0
    while nu <= 4000.0:
        f.write("%10.2f %16.8f\n" % (nu, math.exp(-((nu - 1600.0) / 60.0) ** 2)
                                    + 0.8 * math.exp(-((nu - 3400.0) / 120.0) ** 2)))
        nu += 2.0
EOF
}

WATER_COMMON='atoms_file = "atoms.xyz"
pot_file = "gap_files/water_dipole.gap"
n_species = 2
species = H O
masses = 1.008 15.999
e0 = 0. 0.
random_seed = 12345
write_xyz = 1'

MADIR_INPUT="$WATER_COMMON
soap_radial_legacy_filter = .false.
md_nsteps = 40
md_step = 1.0
thermostat = none
n_exp = 1
exp_labels = ir
exp_data_files = \"ir_exp.dat\"
exp_energy_scales = 1.0e-6
exp_forces = .true.
exp_energies = .true.
ir_stride = 1
ir_nu_max = 4000.0
ir_nu_min = 400.0
ir_resolution = 3000.0
ir_lag_factor = 2"

ESTAT_INPUT='atoms_file = "atoms.xyz"
pot_file = "gap_files/CCLi_estat_ljrep.gap"
n_species = 2
species = C Li
masses = 12.01 6.94
e0 = 0. 0.
estat_method = "gsf"
estat_rcut = 10.0
estat_dsf_alpha = 0.12
estat_damped = .true.
estat_tsf = .true.
estat_sp = .true.
estat_gsf = .true.'

XRD_COMMON='atoms_file = "atoms.xyz"
pot_file = "gap_files/CO.gap"
n_species = 2
species = C O
masses = 12.01 15.99
e0 = -.16138053 0.
n_exp = 1
exp_labels = "xrd"
exp_data_files = "xrd_glassy_carbon_zeng_2017.fq"
exp_n_samples = 201
exp_energy_scales = 100.0
exp_energies = .true.
exp_forces = .true.
md_nsteps = 3
thermostat = berendsen
t_beg = 1000
t_end = 1000
random_seed = 12345
write_xyz = 1
q_range_min = 1.0
q_range_max = 10.0
xrd_output = "q*F(q)"
do_xrd = .true.
structure_factor_window = .true.'

XRD_CHAIN="$XRD_COMMON
do_pair_distribution = .true.
pair_distribution_kde_sigma = 0.1
pair_distribution_n_samples = 201
pair_distribution_partial = .true.
pair_distribution_rcut = 6.0
r_range_min = 0.1
r_range_max = 6.0
do_structure_factor = .true.
structure_factor_from_pdf = .true."

XRD_DEBYE="$XRD_COMMON
xrd_debye = .true.
xrd_rcut = 6.0"

# integ <name> <atoms> <md steps> <ranks> <extra lines> [rank0]: MD with
# dd_integrate = domain against the replicated run. The last argument says
# the run should have integrated on rank 0 instead, with its note.
integ() {
  local name=$1 atoms=$2 steps=$3 ranks=$4 extra=$5 expect=${6:-domain}
  leg "$name" "$atoms" "$steps" "$ranks" auto block "$extra" "dd_integrate = domain"
  local log=$WORK/$name/dd/out.log
  [ "$last_ok" -eq 1 ] || return 0
  if [ "$expect" = domain ] && ! grep -q 'integrated on the domains' "$log"; then
    printf '    FAIL: did not integrate on the domains\n'
    pass=$((pass - 1))
    fail=$((fail + 1))
  elif [ "$expect" = rank0 ] && ! grep -q 'is integrated on rank 0' "$log"; then
    printf '    FAIL: no note that it integrated on rank 0\n'
    pass=$((pass - 1))
    fail=$((fail + 1))
  fi
}

# redrawn <name>: the leg's rcb run must have redrawn its boxes by pair count.
redrawn() {
  [ "$last_ok" -eq 1 ] || return 0
  if ! grep -Eq 'boxes redrawn by pair count at rebuilds: [1-9]' "$WORK/$1/dd/out.log"; then
    printf '    FAIL: the boxes were never redrawn\n'
    pass=$((pass - 1))
    fail=$((fail + 1))
  fi
}

# fallback <name> <ranks> <dd lines>: a decomposition the cell cannot take
# must say so, run replicated, and match the replicated run.
fallback() {
  local name=$1 ranks=$2 dd_lines=$3
  printf '== %s: falls back ==\n' "$name"
  local kind
  for kind in ref dd; do
    mkdir -p "$WORK/$name/$kind"
    ln -sf "$DATA/gap_files" "$WORK/$name/$kind/gap_files"
    cp "$DATA/atoms_897.xyz" "$WORK/$name/$kind/atoms.xyz"
  done
  write_input "$WORK/$name/ref" 0
  write_input "$WORK/$name/dd" 0 "dd_type = legacy
$dd_lines"
  if run_one "$WORK/$name/ref" predict "$ranks" && run_one "$WORK/$name/dd" predict "$ranks" &&
     grep -q 'evaluat.* replicated' "$WORK/$name/dd/out.log" &&
     "$PYTHON" "$COMPARE" "$WORK/$name/ref/trajectory_out.xyz" "$WORK/$name/dd/trajectory_out.xyz" \
       --rtol "$RTOL" --atol "$ATOL"; then
    printf '    PASS  %s\n' "$(grep -a -m1 'dd_type legacy:' "$WORK/$name/dd/out.log" | sed 's/^ *dd_type [a-z]*: //')"
    pass=$((pass + 1))
  else
    printf '    FAIL (logs in %s)\n' "$WORK/$name"
    fail=$((fail + 1))
  fi
}

run_leg() {
  local small=$DATA/atoms_897.xyz
  case $1 in
  x2) leg x2 "$small" 0 2 "2 1 1" ;;
  x3) leg x3 "$small" 0 3 "3 1 1" ;;
  xy4) leg xy4 "$small" 0 4 "2 2 1" ;;
  xyz8) leg xyz8 "$small" 0 8 "2 2 2" ;;
  one3) leg one3 "$small" 0 3 "1 1 1" ;;
  block4) leg block4 "$small" 0 4 "2 1 1" block ;;
  cyclic6) leg cyclic6 "$small" 0 6 "1 3 1" cyclic ;;
  sheared) make_sheared "$WORK/sheared.xyz" && leg sheared "$WORK/sheared.xyz" 0 4 "2 1 2" ;;
  auto8) leg auto8 "$small" 0 8 auto ;;
  volume4) leg volume4 "$small" 0 4 "2 2 1" block "dd_balance = volume" ;;
  slab4) make_slab "$WORK/slab.xyz" 2 && leg slab4 "$WORK/slab.xyz" 0 4 "1 1 4" ;;
  empty6) make_slab "$WORK/slab3.xyz" 3 && leg empty6 "$WORK/slab3.xyz" 3 6 "1 1 6" block "dd_balance = volume" ;;
  md8) leg md8 "$DATA/atoms_7176_vel.xyz" 5 8 "2 2 2" ;;
  md20) leg md20 "$DATA/atoms_7176_vel.xyz" 20 8 auto ;;
  md20_skin) leg md20_skin "$DATA/atoms_7176_vel.xyz" 20 8 auto block "neighbors_buffer = 0.01" ;;
  integ8) integ integ8 "$DATA/atoms_7176_vel.xyz" 30 8 "" ;;
  integ_bussi) integ integ_bussi "$DATA/atoms_7176_vel.xyz" 20 4 "thermostat = bussi
t_beg = 600
t_end = 600" ;;
  integ_nve) integ integ_nve "$DATA/atoms_7176_vel.xyz" 20 6 "thermostat = none" ;;
  integ_skin) integ integ_skin "$DATA/atoms_7176_vel.xyz" 20 8 "neighbors_buffer = 0.01" ;;
  integ_empty) make_slab "$WORK/slab3i.xyz" 3 && integ integ_empty "$WORK/slab3i.xyz" 10 6 "dd_grid = 1 1 6
dd_balance = volume" ;;
  integ_gle) integ integ_gle "$DATA/atoms_7176_vel.xyz" 5 4 "thermostat = langevin
tau_t = 100." rank0 ;;
  narrow) fallback narrow 4 "dd_grid = 4 1 1" ;;
  estat) full_leg estat predict 2 "$DATA_ROOT/C/atoms_24219.xyz" "$DATA_ROOT/CLi/gap_files" - "$ESTAT_INPUT" ;;
  xrd_view) full_leg xrd_view md 4 "$DATA_ROOT/C/atoms_7176.xyz" "$DATA_ROOT/xrd_mad/gap_files" \
    "$DATA_ROOT/xrd_mad/xrd_glassy_carbon_zeng_2017.fq" "$XRD_CHAIN" ;;
  xrd_debye_view) full_leg xrd_debye_view md 4 "$DATA_ROOT/C/atoms_7176.xyz" "$DATA_ROOT/xrd_mad/gap_files" \
    "$DATA_ROOT/xrd_mad/xrd_glassy_carbon_zeng_2017.fq" "$XRD_DEBYE" ;;
  dipole) tile_water "$WORK/water222.xyz" 2 2 2 &&
    full_leg dipole predict 8 "$WORK/water222.xyz" "$DATA_ROOT/water_dipole/gap_files" - "$WATER_COMMON" ;;
  mad_ir) tile_water "$WORK/water221.xyz" 2 2 1 && make_ir_exp "$WORK/ir_exp.dat" &&
    full_leg mad_ir md 4 "$WORK/water221.xyz" "$DATA_ROOT/water_dipole/gap_files" "$WORK/ir_exp.dat" "$MADIR_INPUT" ;;
  vdw_ts) leg vdw_ts "$DATA/atoms_7176.xyz" 0 4 auto block "vdw_type = ts
vdw_rcut = 6.0" ;;
  prime5) leg prime5 "$DATA/atoms_7176.xyz" 0 5 auto ;;
  prime7) leg prime7 "$DATA/atoms_7176_vel.xyz" 10 7 auto ;;
  cluster) make_sphere "$WORK/cluster.xyz" in 0.42 && leg cluster "$WORK/cluster.xyz" 0 8 auto ;;
  pore) make_sphere "$WORK/pore.xyz" out 0.35 && leg pore "$WORK/pore.xyz" 0 6 auto ;;
  integ_droplet) make_sphere "$WORK/droplet.xyz" vapour 0.3 &&
    integ integ_droplet "$WORK/droplet.xyz" 20 8 "neighbors_buffer = 0.01" && redrawn integ_droplet ;;
  *) die "no such leg: $1" ;;
  esac
}

check_setup
mkdir -p "$WORK"
if [ -z "${TURBOGAP_KEEP:-}" ]; then
  trap 'rm -rf "$WORK"' EXIT
else
  trap 'printf "staging kept in %s\n" "$WORK"' EXIT
fi

pass=0
fail=0
last_ok=0
common="x2 x3 xy4 xyz8 one3 block4 cyclic6 sheared auto8 volume4 slab4 md8 md20 md20_skin integ8 integ_bussi
integ_nve integ_skin integ_empty integ_gle empty6 vdw_ts estat xrd_view xrd_debye_view dipole mad_ir"
if [ "$DD_TYPE" = rcb ]; then
  legs=${*:-$common prime5 prime7 cluster pore integ_droplet}
else
  legs=${*:-$common narrow}
fi
for l in $legs; do
  run_leg "$l"
done

printf '\npassed: %d   failed: %d\n' "$pass" "$fail"
[ "$fail" -eq 0 ]
