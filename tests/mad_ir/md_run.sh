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
# Absolute, because the run happens in a scratch directory and a relative root
# stops resolving the moment we cd into it. BIN in the environment overrides,
# so a variant tree (bin-dbg, bin-gle, ...) can be tested without rebuilding
# the one everything else names.
ROOT="$(cd "${1:-$HERE/../..}" && pwd)"
CASE="${2:-$HOME/work/cpu_vs_gpu_tests/input/water_dipole}"
BIN="${BIN:-$ROOT/bin/turbogap}"
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
# A bare PASS/FAIL, for checks that are not a grep over a log. Both go through
# the same counter as check(), so a failure here cannot be reported as a pass.
pass() { echo "  PASS  $1"; }
bad()  { echo "  FAIL  $1"; fail=$((fail + 1)); }

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

echo "==> 7. auxiliary variables (ir_acf_mode = exponential)"
#
# The same bias with the running correlation carried as auxiliary degrees of
# freedom instead of re-averaged from the buffer. auxverify checks the
# mathematics; what is checked here is that the keywords reach it, that the
# resulting force is still a force -- momentum conserved, because a rigid
# translation still cannot change a dipole -- and that the restart carries the
# filter rather than silently starting it again.
rm -f ir_restart.dat
cat >input <<EOF
atoms_file = "atoms.xyz"
pot_file = "gap_files/water_dipole.gap"
n_species = 2
species = H O
masses = 1.008 15.999
e0 = 0. 0.
random_seed = 12345
soap_radial_legacy_filter = .false.

do_md = .true.
md_nsteps = 60
md_step = 1.0
thermostat = none
write_xyz = 1

n_exp = 1
exp_labels = ir
exp_data_files = "ir_exp.dat"
exp_energy_scales = 1.0e-6
exp_forces = .true.
exp_energies = .true.

ir_stride = 1
ir_nu_max = 4000.0
ir_nu_min = 400.0
ir_resolution = 3000.0
ir_lag_factor = 2
ir_acf_mode = exponential
ir_tau_mem = 20.0
ir_write_spectrum = .true.
EOF
"$BIN" md >run8.log 2>&1
check "exponential mode accepted" "= exponential" run8.log
check "the memory constant is read" "ir_tau_mem" run8.log
[ -f ir_spectrum.dat ] || {
  echo "  FAIL  spectrum not written in exponential mode"
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
dec=0
for tok in (r[fcol+j] for r in frames[-1][1] for j in range(3)):
    if "." in tok:
        dec=max(dec, len(tok.split(".")[1]))
floor=nat*10.0**(-dec)
before=max(stats(i)[0] for i in range(0,20))
after =max(stats(i)[0] for i in range(30,len(frames)))
net   =max(stats(i)[1] for i in range(30,len(frames)))
ok=True
print("  %s  no bias while the filter is still charging (max|f| = %.3e)" %
      ("PASS" if before==0.0 else "FAIL", before)); ok &= before==0.0
print("  %s  bias applied once it has (max|f| = %.3e)" %
      ("PASS" if after>0.0 else "FAIL", after)); ok &= after>0.0
good = net <= floor
print("  %s  momentum conserved (max|sum f| = %.3e, floor %.3e)" %
      ("PASS" if good else "FAIL", net, floor)); ok &= good
sys.exit(0 if ok else 1)
PY
[ $? -eq 0 ] || fail=$((fail + 1))

# The filter is state: a second run must take it up, and must refuse one
# written under a different memory constant.
"$BIN" md >run9.log 2>&1
check "the filter is resumed" "resumed, frames:" run9.log
sed -i "s/ir_tau_mem = 20.0/ir_tau_mem = 40.0/" input
"$BIN" md >runA.log 2>&1
check "a different tau_mem is refused" "ir_tau_mem=" runA.log
check "and a fresh filter is started" "fresh ensemble" runA.log

echo "==> 8. auxiliary variables on the PREDICTION path"
#
# do_ir has its own setup routine, and an option wired into the bias path but
# not into that one would be accepted, printed back, and then quietly ignored.
# So this asks for more than "it ran": the spectrum a prediction run writes
# under the exponential estimator must actually DIFFER from the block one.
# Identical output here would mean ir_acf_mode never reached the estimator.
rm -f ir_restart.dat ir_spectrum.dat
cat >input <<EOF
atoms_file = "atoms.xyz"
pot_file = "gap_files/water_dipole.gap"
n_species = 2
species = H O
masses = 1.008 15.999
e0 = 0. 0.
random_seed = 12345
soap_radial_legacy_filter = .false.

do_md = .true.
md_nsteps = 200
md_step = 1.0
thermostat = none
write_xyz = 200

do_ir = .true.
ir_stride = 1
ir_nu_min = 400.0
ir_nu_max = 4000.0
ir_lag_factor = 2
ir_acf_mode = exponential
ir_tau_mem = 30.0
EOF
"$BIN" md >runB1.log 2>&1
check "prediction accepts exponential" "= exponential" runB1.log
if [ -f ir_spectrum.dat ]; then
  cp ir_spectrum.dat spec_exp.dat
else
  bad "no spectrum from the exponential prediction"
fi
sed -i "s/ir_acf_mode = exponential/ir_acf_mode = block/" input
"$BIN" md >runB2.log 2>&1
if [ -f spec_exp.dat ] && [ -f ir_spectrum.dat ]; then
  if cmp -s spec_exp.dat ir_spectrum.dat; then
    bad "block and exponential predictions are identical; the mode is ignored"
  else
    pass "the estimator reaches the prediction path"
  fi
fi

# And the same refusal the bias path gives, from the same shared check.
sed -i "s/ir_acf_mode = block/ir_acf_mode = exponential/; s/ir_tau_mem = 30.0/ir_tau_mem = 0.4/" input
"$BIN" md >runB3.log 2>&1
check "a tau_mem below the sampling interval is refused" "retain nothing" runB3.log

echo "==> 9. the dissimilarity columns in thermo.log"
#
# The MAD energy is 1/2 * exp_energy_scales * sum_k wgt_k (I_k - I_k^exp)^2, so
# the reported Dissimilarity must satisfy E_exp = 0.5 * scale * Dissimilarity
# exactly. That identity is the check: it ties the new column to the energy the
# force is actually the gradient of, and no independent reimplementation of the
# residual could satisfy it by accident.
#
# It matters because exp_energy_scales is RAMPED. E_exp alone therefore mixes
# how far the spectrum is from the experiment with how much that is currently
# being charged for, and can fall while the agreement gets worse.
rm -f ir_restart.dat
write_input
"$BIN" md >runC1.log 2>&1
head -1 thermo.log | grep -q "Dissimilarity" && pass "Dissimilarity column present" || bad "no Dissimilarity column"
head -1 thermo.log | grep -q "Rel_error"     && pass "Rel_error column present"     || bad "no Rel_error column"
# An IR-only deck switches params%do_exp back off, so before this it reported
# neither its MAD energy nor its mismatch.
head -1 thermo.log | grep -q "E_exp" && pass "E_exp reported for an IR-only run" || bad "no E_exp column"

python3 - <<'PY'
import sys
rows=[l.split() for l in open("thermo.log") if l.strip() and not l.startswith("#")]
scale=1.0e-6            # exp_energy_scales in the deck above, not ramped
eexp=[float(r[5]) for r in rows]
diss=[float(r[7]) for r in rows]
rel =[float(r[8]) for r in rows]
ok=True
# E_exp = 0.5 * scale * Dissimilarity, to the printed precision of E_exp (8 dp)
err=max(abs(e-0.5*scale*d) for e,d in zip(eexp,diss))
good = err <= 1e-8
print("  %s  E_exp == 0.5 * scale * Dissimilarity (max err %.1e)" %
      ("PASS" if good else "FAIL", err)); ok &= good
# Zero while the ensemble is still filling -- there is no spectrum yet, so
# there is no mismatch -- and positive once there is one. Same onset the bias
# itself has, and checked the same way.
n_zero = sum(1 for d in diss if d == 0.0)
tail = diss[n_zero:]
good = n_zero > 0 and len(tail) > 0 and all(d > 0.0 for d in tail)
print("  %s  zero while filling (%d rows), positive after (%d rows)" %
      ("PASS" if good else "FAIL", n_zero, len(tail))); ok &= good
# Rel_error is dimensionless and, for a model with no vibrational bands against
# a synthetic spectrum that has them, should be of order one rather than tiny.
rtail = rel[n_zero:]
good = all(0.0 < r < 10.0 for r in rtail) and max(rtail) > 0.1
print("  %s  Rel_error is dimensionless and O(1) (%.3f .. %.3f once active)" %
      ("PASS" if good else "FAIL", min(rtail), max(rtail))); ok &= good
sys.exit(0 if ok else 1)
PY
[ $? -eq 0 ] || fail=$((fail + 1))

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
