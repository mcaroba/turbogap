#!/usr/bin/env bash
#
# The generalized Langevin thermostat, end to end through an actual MD run.
#
# gleverify checks the mathematics; this checks the wiring, which is a
# different question and fails in different ways. It asks:
#
#   1. do the keywords reach the thermostat, and does the report say what was
#      actually built (kernel, ns, the timescales the matrix implies)
#   2. does a run sample the target temperature. The dipole potential used here
#      contributes NO energy and NO forces, so the atoms are free particles and
#      the O step is exact -- which makes <T> = t_beg an identity rather than
#      something that holds to within a splitting error. Any systematic offset
#      is a defect, not a discretisation.
#   3. does a matrix-file kernel give the ns it was given, and still thermalise
#   4. is the bath written, resumed, and REFUSED when it belongs to another run
#   5. do the input guards fire: gle without a drift matrix, a covariance
#      without a drift matrix, a matrix file that is not square, and an A that
#      violates the fluctuation-dissipation condition
#   6. is a seeded run reproducible
#
# Usage:  ./md_run.sh [path/to/turbogap] [path/to/water_dipole/input]

set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
# Absolute, because the checks below run from a scratch directory and a
# relative root stops resolving the moment we cd into it.
ROOT="$(cd "${1:-$HERE/../..}" && pwd)"
CASE="${2:-$HOME/work/cpu_vs_gpu_tests/input/water_dipole}"
# BIN in the environment overrides, so a variant tree (bin-dbg, bin-gle, ...)
# can be tested without rebuilding the one everything else names.
BIN="${BIN:-$ROOT/bin/turbogap}"
D="${TMPDIR:-/tmp}/gle_md_$$"

[ -x "$BIN" ] || {
  echo "no binary at $BIN; run make first, or set BIN=" >&2
  exit 1
}
[ -d "$CASE/gap_files" ] || {
  echo "SKIP: no dipole case at $CASE" >&2
  exit 0
}

fail=0
pass() { echo "  PASS  $1"; }
bad()  { echo "  FAIL  $1"; fail=$((fail + 1)); }

check() { # name  expected-substring  file
  if grep -qF "$2" "$3"; then pass "$1"; else
    bad "$1  (expected to find: $2)"
  fi
}
check_not() { # name  forbidden-substring  file
  if grep -qF "$2" "$3"; then bad "$1  (found: $2)"; else pass "$1"; fi
}

rm -rf "$D"; mkdir -p "$D"; cd "$D"
cp -rL "$CASE/gap_files" .
head -8 "$CASE"/*_test.xyz >atoms.xyz

# An ns = 2 kernel. The off-diagonal coupling is antisymmetric and the
# diagonal non-negative, which is what makes A + A^T positive semi-definite and
# so satisfies the fluctuation-dissipation condition against C = kB T I.
cat >gle_A.dat <<'EOF'
# a two-auxiliary-variable kernel, fs^-1
  0.002   0.05    0.02
 -0.05    0.02    0.0
 -0.02    0.0     0.004
EOF

# The same matrix with one number missing, so the count is not a square.
cat >gle_A_ragged.dat <<'EOF'
  0.002   0.05    0.02
 -0.05    0.02    0.0
 -0.02    0.0
EOF

# A negative friction on the physical momentum: A + A^T has a negative
# eigenvalue, so the noise variance the propagator needs does not exist.
cat >gle_A_bad.dat <<'EOF'
 -0.01    0.05
 -0.05    0.02
EOF

# C at the wrong order for a 3x3 A.
cat >gle_C_small.dat <<'EOF'
 0.02585
EOF

write_input() { # 1 thermostat  2 nsteps  3 extra lines  4 seed
  cat >input <<EOF
atoms_file = "atoms.xyz"
pot_file = "gap_files/water_dipole.gap"
n_species = 2
species = H O
masses = 1.008 15.999
e0 = 0. 0.
random_seed = ${4:-12345}

do_md = .true.
md_nsteps = ${2:-15000}
md_step = 1.0
write_xyz = ${2:-15000}

thermostat = ${1:-langevin}
tau_t = 20.0
t_beg = 300.0
t_end = 300.0
${3:-}
EOF
}

mean_temp() { awk 'NR>2{n++; s+=$3} END{if(n>0) printf "%.3f", s/n}' thermo.log; }

# Is the mean temperature within tol per cent of the target? The spread of a
# single sample is sqrt(2/Ndof)*T = 110 K here, and the samples decorrelate on
# tau_t, so 15000 steps at tau_t = 20 fs give a standard error near 4 K. 5 per
# cent is about four of those: loose enough not to flicker, tight enough that a
# thermostat off by a factor of two in the noise cannot pass.
temp_ok() { # name target tol_percent
  local t; t=$(mean_temp)
  awk -v t="$t" -v tgt="$2" -v tol="$3" -v n="$1" '
    BEGIN{ d = (t-tgt)/tgt*100; a = d<0?-d:d;
      if (a<=tol) printf "  PASS  %s  (mean T = %.1f K, %+.1f%%)\n", n, t, d;
      else        { printf "  FAIL  %s  (mean T = %.1f K, %+.1f%%, tol %g%%)\n", n, t, d, tol; exit 1 } }'
  [ $? -eq 0 ] || fail=$((fail + 1))
}

echo "==> 1/2. langevin: report, and the free-particle temperature"
write_input langevin 15000
"$BIN" md >run1.log 2>&1
check "the thermostat reports itself"  "Generalized Langevin thermostat" run1.log
check "kernel named"                   "langevin"                       run1.log
check "no auxiliary variables"         "auxiliary DOF:                0" run1.log
check "friction from tau_t"            "fastest mode:           20.0000" run1.log
check "fresh bath on a first run"      "fresh bath"                     run1.log
[ -f gle_restart.dat ] && pass "bath written" || bad "bath not written"
temp_ok "samples the target temperature" 300 5

echo
echo "==> 3. gle from a matrix file"
write_input gle 15000 'gle_a_file = "gle_A.dat"'
"$BIN" md >run2.log 2>&1
check "kernel read from the file"      "auxiliary DOF:                2" run2.log
check "named gle"                      "*) kernel:                     gle" run2.log
temp_ok "still samples the target"     300 5

echo
echo "==> 4. the bath is resumed, and refused when it is not ours"
# The bath from run 3 has ns = 2; rerun with it in place and it must be taken up.
write_input gle 200 'gle_a_file = "gle_A.dat"'
"$BIN" md >run3.log 2>&1
check "resumed from the restart file"  "bath resumed from"              run3.log

# The same file offered to a run with no auxiliary variables describes a
# different kernel, and adopting it would run one bath under another's
# propagator.
write_input langevin 200
"$BIN" md >run4.log 2>&1
check "refused for a different ns"     "describes a different run"      run4.log
check "and says it starts fresh"       "fresh bath"                     run4.log

echo
echo "==> 5. input guards"
write_input gle 50
"$BIN" md >run5.log 2>&1
check "gle without a drift matrix"     "thermostat = gle needs gle_a_file" run5.log

write_input none 50 'gle_c_file = "gle_C_small.dat"'
"$BIN" md >run6.log 2>&1
check "covariance without a drift"     "gle_c_file was given without gle_a_file" run6.log

write_input gle 50 'gle_a_file = "gle_A_ragged.dat"'
"$BIN" md >run7.log 2>&1
check "a matrix file that is not square" "is not square"                run7.log

write_input gle 50 'gle_a_file = "gle_A_bad.dat"'
"$BIN" md >run8.log 2>&1
check "an A that breaks fluctuation-dissipation" "not positive semi-definite" run8.log

write_input gle 50 'gle_a_file = "gle_A.dat"
gle_c_file = "gle_C_small.dat"'
"$BIN" md >run9.log 2>&1
check "a C of the wrong order"         "but A sets the order"           run9.log

echo
echo "==> 6. a seeded run is reproducible"
rm -f gle_restart.dat
write_input langevin 300 'gle_restart = .false.' 777
"$BIN" md >runA.log 2>&1
cp trajectory_out.xyz traj_a.xyz
write_input langevin 300 'gle_restart = .false.' 777
"$BIN" md >runB.log 2>&1
if cmp -s traj_a.xyz trajectory_out.xyz; then
  pass "same seed, identical trajectory"
else
  bad "same seed, trajectories differ"
fi

echo
echo "==> 7. the bath lives on rank 0 only"
# compute_md runs the whole MD step inside IF (rank == 0) and broadcasts
# positions and velocities afterwards, so the auxiliary variables and the
# random draws belong to one rank and there is nothing to reduce. If that is
# true, a seeded run on two ranks is BIT-IDENTICAL to the same run on one; if
# some rank ever drew its own noise, this is the check that would catch it.
if command -v mpirun >/dev/null 2>&1; then
  rm -f gle_restart.dat
  write_input gle 300 'gle_a_file = "gle_A.dat"
gle_restart = .false.' 999
  "$BIN" md >runC.log 2>&1
  cp trajectory_out.xyz traj_serial.xyz
  if mpirun -np 2 --oversubscribe "$BIN" md >runD.log 2>&1 ||
     mpirun -np 2 "$BIN" md >runD.log 2>&1; then
    if cmp -s traj_serial.xyz trajectory_out.xyz; then
      pass "1 rank and 2 ranks agree bit for bit"
    else
      bad "1 rank and 2 ranks disagree"
    fi
  else
    echo "  SKIP  mpirun present but the run failed; binary may be serial"
  fi
else
  echo "  SKIP  no mpirun"
fi

echo
echo "==> 8. the thermostat's energy ledger in thermo.log"
# A generalized Langevin run does not conserve E_pot + E_kin and is not meant
# to. What it conserves is E_pot + E_kin - E_thermo + E_cmrem, and with this
# potential -- which has no energy and no forces -- that is an IDENTITY, not an
# approximation: every joule that entered or left the atoms did so through the
# thermostat or through the centre-of-mass removal, and both are accounted.
# So the column is checked for EXACT constancy, not for small drift.
#
# E_cmrem is in there because remove_cm_vel discards kinetic energy every step
# while a stochastic thermostat keeps pumping energy back into the drift it
# removes. Left out, E_cons drifts by ~1 eV over 300 steps and looks like an
# integrator fault.
rm -f gle_restart.dat
write_input langevin 300 'gle_restart = .false.' 4242
"$BIN" md >runE.log 2>&1
head -1 thermo.log | grep -q "E_thermo" && pass "E_thermo column present" || bad "no E_thermo column"
head -1 thermo.log | grep -q "E_cons"   && pass "E_cons column present"   || bad "no E_cons column"

python3 - <<'PY'
import sys
rows=[l.split() for l in open("thermo.log") if l.strip() and not l.startswith("#")]
ek  =[float(r[3]) for r in rows]
eth =[float(r[6]) for r in rows]
ecm =[float(r[7]) for r in rows]
ec  =[float(r[8]) for r in rows]
ok=True
drift=max(abs(c-ec[0]) for c in ec)
span =max(ek)-min(ek)
moved=abs(eth[-1])
# The printed precision is 8 decimals, so anything at or below 1e-8 is the
# format, not the physics.
good = drift <= 1e-8
print("  %s  E_cons is constant (drift %.1e eV, while E_kin spans %.3f eV and the"
      " thermostat moved %.3f eV)" % ("PASS" if good else "FAIL", drift, span, moved))
ok &= good
# and the column really is the combination it claims to be
err=max(abs(k-t+c-e) for k,t,c,e in zip(ek,eth,ecm,ec))
good = err <= 1e-7
print("  %s  E_cons == E_kin - E_thermo + E_cmrem (max err %.1e)" %
      ("PASS" if good else "FAIL", err)); ok &= good
# a thermostat that did nothing would pass the above trivially
good = moved > 0.1
print("  %s  the thermostat actually moved energy" % ("PASS" if good else "FAIL")); ok &= good
sys.exit(0 if ok else 1)
PY
[ $? -eq 0 ] || fail=$((fail + 1))

echo
echo "==> 9. the header names exactly as many columns as the data has"
#
# This is the invariant the dynamic header exists to hold. It used to fail
# silently: the header was a fixed string while the columns were already
# conditional, so a do_exp run wrote a column the header did not name and
# everything after it was misread by one place -- in xrd_mad the column labelled
# "Pressure" actually held E_exp, and the pressure had no label at all.
#
# Checked over the combinations that switch columns on and off, because the bug
# only appears when two optional groups are on at once and no single-flag test
# would have caught it.
cat >hdr_A.dat <<'EOF'
  0.002   0.05    0.02
 -0.05    0.02    0.0
 -0.02    0.0     0.004
EOF
python3 - <<'PY'
import math
with open("hdr_ir.dat","w") as f:
    nu = 0.0
    while nu <= 4000.0:
        f.write("%10.2f %16.8f\n" % (nu, math.exp(-((nu-1600.0)/60.0)**2)))
        nu += 2.0
PY
# An ir observable needs the descriptor second derivatives, which the legacy
# radial filter cannot supply -- without this line the run aborts.
IRBLK='soap_radial_legacy_filter = .false.
n_exp = 1
exp_labels = ir
exp_data_files = "hdr_ir.dat"
exp_energy_scales = 1.0e-6
exp_forces = .true.
exp_energies = .true.
ir_stride = 1
ir_nu_max = 4000.0
ir_nu_min = 400.0
ir_resolution = 3000.0
ir_lag_factor = 2'
GLEBLK='gle_a_file = "hdr_A.dat"
gle_restart = .false.'

hdr_combo() { # 1 name  2 thermostat  3 extra  4 expected columns
  # DELETE thermo.log FIRST. A run that aborts leaves the previous combo's file
  # in place, and reading that compares a stale header against its own stale
  # data -- which agree, so a combination that never ran reports a pass. That
  # happened here: the IR combos aborted on the legacy radial filter and the
  # check went green on the file left by the run before.
  rm -f thermo.log
  write_input "$2" 40 "$3"
  "$BIN" md >"hdr_$1.log" 2>&1
  if [ ! -f thermo.log ]; then
    bad "$1: the run produced no thermo.log ($(grep -m1 -i error "hdr_$1.log" | cut -c1-60))"
    return
  fi
  local h d
  h=$(head -1 thermo.log | sed 's/^#//' | wc -w)
  d=$(sed -n 2p thermo.log | wc -w)
  if [ "$h" != "$d" ]; then
    bad "$1: header names $h columns but data has $d"
  elif [ -n "${4:-}" ] && [ "$h" != "$4" ]; then
    # The count is also asserted outright, so that a combination silently
    # losing a whole column group still fails even though it stays consistent.
    bad "$1: consistent at $h columns but $4 were expected"
  else
    pass "$1: header names $h columns, data has $d"
  fi
}
hdr_combo "plain"        none     ""                       6
hdr_combo "gle"          gle      "$GLEBLK"                9
hdr_combo "IR"           none     "$IRBLK"                 9
hdr_combo "gle+IR"       gle      "$GLEBLK
$IRBLK"                                                    12
hdr_combo "write_lv+gle" gle      "write_lv = .true.
$GLEBLK"                                                   18
hdr_combo "write_lv+gle+IR" gle   "write_lv = .true.
$GLEBLK
$IRBLK"                                                    21
# write_lv ALONE is deliberately left as it always was: 6 names for 15 columns,
# the nine lattice columns unnamed. Naming them would change thermo.log for
# every existing write_lv deck, and relax_gd and its three variants compare
# against the frozen baseline binary, which cannot produce a new header. They
# are named as soon as anything follows them, because only then does leaving
# them out misalign the names that do appear -- which the two cases above check.
rm -f thermo.log
write_input none 40 'write_lv = .true.'
"$BIN" md >hdr_lv.log 2>&1
h=$(head -1 thermo.log | sed 's/^#//' | wc -w)
d=$(sed -n 2p thermo.log | wc -w)
if [ "$h" = "6" ] && [ "$d" = "15" ]; then
  pass "write_lv alone keeps its historical header ($h names, $d columns)"
else
  bad "write_lv alone changed: $h names, $d columns (expected 6 and 15)"
fi

echo
if [ "$fail" -eq 0 ]; then
  echo "==> gle/md_run.sh: all checks passed"
  cd /; rm -rf "$D"
else
  echo "==> gle/md_run.sh: $fail CHECK(S) FAILED   (logs kept in $D)"
  exit 1
fi
