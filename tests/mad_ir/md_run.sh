#!/usr/bin/env bash
#
# The MAD IR bias, end to end through an actual MD run.
#
# irverify and madverify check the mathematics; this checks the wiring, which
# is a different question and fails in different ways. It asks:
#
#   1. is the ensemble sized and reported from the keywords as intended
#   2. is the bias withheld until the ensemble is full, and applied after
#   3. does the bias conserve momentum -- sum_j f_j must vanish, because a
#      rigid translation cannot change a dipole. This is the check that a
#      per-atom scatter or an MPI reduction gone wrong cannot survive.
#   4. does a restart resume the dipole history
#   5. is a restart written for different sizing refused rather than adopted
#   6. do the input guards fire: no dipole model, legacy filter seed on,
#      sampling too coarse for the requested wavenumber, missing spectrum
#
# The dipole model contributes no energy and no forces, so in this case the
# MAD bias is the ONLY force in the system. That makes both its onset and its
# net force directly readable from the trajectory rather than inferred.
#
# Usage:  ./md_run.sh [path/to/turbogap] [path/to/water_dipole/input]

set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="${1:-$(cd "$HERE/../.." && pwd)}"
CASE="${2:-$HOME/work/cpu_vs_gpu_tests/input/water_dipole}"
BIN="$ROOT/bin/turbogap"
D="${TMPDIR:-/tmp}/mad_ir_md_$$"

[ -x "$BIN" ] || {
  echo "no binary at $BIN; run make first" >&2
  exit 1
}
[ -d "$CASE/gap_files" ] || {
  echo "SKIP: no dipole case at $CASE" >&2
  exit 0
}

fail=0
check() { # name  expected-substring  file
  if grep -qF "$2" "$3"; then
    echo "  PASS  $1"
  else
    echo "  FAIL  $1  (expected to find: $2)"
    fail=$((fail + 1))
  fi
}

rm -rf "$D"
mkdir -p "$D"
cd "$D"
cp -rL "$CASE/gap_files" .
head -8 "$CASE"/*_test.xyz >atoms.xyz

# A synthetic experimental spectrum. Its shape is irrelevant to the wiring;
# what matters is that it is read, restricted to the fitted range, and used.
# Two plain columns and no header: read_exp_data counts lines and reads two
# reals from each, so a comment line is a runtime error rather than a skip.
python3 - <<'PY'
import math
with open("ir_exp.dat","w") as f:
    nu = 0.0
    while nu <= 4000.0:
        f.write("%10.2f %16.8f\n" % (nu,
            math.exp(-((nu-1600.0)/60.0)**2) + 0.8*math.exp(-((nu-3400.0)/120.0)**2)))
        nu += 2.0
PY

write_input() {
  cat >input <<EOF
atoms_file = "atoms.xyz"
pot_file = "${1:-gap_files/water_dipole.gap}"
n_species = 2
species = H O
masses = 1.008 15.999
e0 = 0. 0.
random_seed = 12345
soap_radial_legacy_filter = ${2:-.false.}

do_md = .true.
md_nsteps = ${3:-60}
md_step = 1.0
thermostat = none
write_xyz = 1

# The IR spectrum is an experimental observable like any other: named in
# exp_labels, weighted by exp_energy_scales, and its gradient reaches the
# forces only because exp_forces is set.
n_exp = 1
exp_labels = ir
exp_data_files = "${4:-ir_exp.dat}"
exp_energy_scales = 1.0e-6
exp_forces = .true.
exp_energies = .true.

ir_stride = ${5:-1}
ir_nu_max = 4000.0
ir_nu_min = 400.0
ir_resolution = ${6:-3000.0}
ir_lag_factor = 2
ir_write_spectrum = .true.
EOF
}

echo "==> 1/2/3. sizing, onset and momentum"
write_input
"$BIN" md >run1.log 2>&1
check "ensemble sized from the keywords" "ensemble size:               24" run1.log
check "fresh ensemble reported" "fresh ensemble" run1.log
[ -f ir_restart.dat ] || {
  echo "  FAIL  restart file not written"
  fail=$((fail + 1))
}
[ -f ir_spectrum.dat ] || {
  echo "  FAIL  spectrum file not written"
  fail=$((fail + 1))
}

python3 - <<'PY'
import re, sys
frames=[]; lines=open("trajectory_out.xyz").read().split("\n"); i=0
while i < len(lines) and lines[i].strip():
    n=int(lines[i]); frames.append((lines[i+1],[lines[i+2+k].split() for k in range(n)])); i+=n+2
props=re.search(r"Properties=(\S+)", frames[0][0]).group(1).split(":")
col=0; fcol=None
for k in range(0,len(props),3):
    if props[k].startswith("force"): fcol=col
    col+=int(props[k+2])
def stats(idx):
    F=[[float(r[fcol+j]) for j in range(3)] for r in frames[idx][1]]
    return (max(max(abs(v) for v in r) for r in F),
            max(abs(sum(r[j] for r in F)) for j in range(3)))
nat=len(frames[0][1])
# The trajectory writes forces to a fixed number of decimals, so the net force
# read back cannot be smaller than the rounding of the terms that made it. The
# floor is nat*ulp/2; comparing against anything smaller would be testing the
# output format rather than the physics.
dec=0
for tok in (r[fcol+j] for r in frames[-1][1] for j in range(3)):
    if "." in tok:
        dec=max(dec, len(tok.split(".")[1].rstrip("0123456789Ee+-")) or len(tok.split(".")[1]))
ulp=10.0**(-dec)
floor=nat*ulp
before=max(stats(i)[0] for i in range(0,20))
after =max(stats(i)[0] for i in range(30,len(frames)))
net   =max(stats(i)[1] for i in range(30,len(frames)))
ok=True
print("  %s  no bias before the ensemble is full (max|f| = %.3e)" %
      ("PASS" if before==0.0 else "FAIL", before)); ok &= before==0.0
print("  %s  bias applied once it is full (max|f| = %.3e)" %
      ("PASS" if after>0.0 else "FAIL", after)); ok &= after>0.0
good = net <= floor
print("  %s  momentum conserved (max|sum f| = %.3e; trajectory prints %d decimals,"
      " so the floor is %.3e)" % ("PASS" if good else "FAIL", net, dec, floor))
ok &= good
sys.exit(0 if ok else 1)
PY
[ $? -eq 0 ] || fail=$((fail + 1))

echo "==> 4. restart resumes"
write_input gap_files/water_dipole.gap .false. 3
"$BIN" md >run2.log 2>&1
check "history resumed" "resumed, frames:             24" run2.log

echo "==> 5. restart with different sizing refused"
write_input gap_files/water_dipole.gap .false. 3 ir_exp.dat 1 1500.0
"$BIN" md >run3.log 2>&1
check "sizing mismatch refused" "restart was written for" run3.log
check "fresh ensemble started instead" "fresh ensemble" run3.log

echo "==> 6. input guards"
write_input gap_files/water_dipole.gap .true. 3
"$BIN" md >run4.log 2>&1
check "legacy filter seed refused" "soap_radial_legacy_filter = .false." run4.log

write_input gap_files/water_dipole.gap .false. 3 nope.dat
"$BIN" md >run5.log 2>&1
check "missing spectrum refused" "nope.dat" run5.log

write_input gap_files/water_dipole.gap .false. 3 ir_exp.dat 8
"$BIN" md >run6.log 2>&1
check "coarse sampling refused" "sample more often" run6.log

cp -r gap_files gap_nodip
sed -i "s/dipole_model *= *\.true\./dipole_model = .false./" gap_nodip/water_dipole.gap
write_input gap_nodip/water_dipole.gap .false. 3
"$BIN" md >run7.log 2>&1
check "missing dipole model refused" "needs a dipole model" run7.log

echo
if [ "$fail" -eq 0 ]; then
  echo "==> all MAD IR MD checks passed"
  cd /
  rm -rf "$D"
  exit 0
else
  echo "==> $fail MAD IR MD check(s) FAILED; working directory kept at $D"
  exit 1
fi
