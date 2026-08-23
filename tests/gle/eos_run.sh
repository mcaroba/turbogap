#!/usr/bin/env bash
#
# Thermal equation of state under the generalized Langevin thermostat.
#
# SLOW -- tens of minutes on a real potential. Not part of tests/gle/run.sh;
# run it when the thermostat's thermodynamics need checking rather than its
# mechanics.
#
# What this asks that nothing else does. gleverify proves the propagator is
# exact and that it samples kB T for a free particle; md_run.sh proves the
# wiring holds and that the temperature comes out right. Neither says the
# thermostat produces correct averages of anything that depends on the
# POTENTIAL. A thermostat can hold the kinetic energy at the target and still
# sample the wrong configurational distribution -- that is precisely what
# Berendsen does -- and the observable that exposes it is the pressure, which
# is a configurational average.
#
# So: run NVT at several volumes with a thermostat already known to sample the
# canonical ensemble (bussi), and again with the generalized Langevin one, from
# the same equilibrated starting structure. If the GLE is a correct canonical
# sampler the two equations of state must agree within the sampling error. If
# it merely holds the temperature they will not.
#
# The comparison is between thermostats, not against a reference curve, because
# the absolute P(V) of this potential at this density is not known to better
# than the thing being measured. A difference between two thermostats on the
# same structure is.
#
# The system is diamond rather than the amorphous carbon the case ships with;
# see the structure builder below for why that matters.
#
# Measured 2026-08-23, 216-atom diamond at 500 K, 400 fs of production per
# point: the two thermostats agree to 0.55%, 0.63%, 1.13% and 0.93% of the
# pressure across a +-6% volume range, and both curves fall monotonically. On
# 897 atoms of amorphous carbon the same comparison gave 5-30% away from the
# reference volume -- which is the non-ergodicity described below and not a
# property of either thermostat.
#
# Usage:  ./eos_run.sh [path/to/turbogap_root]
#   BIN=      binary to use            (default <root>/bin/turbogap)
#   NRANKS=   mpi ranks                (default 6)
#   NPROD=    production steps/point   (default 600)
#   TEMP=     target temperature in K  (default 1000)

set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "${1:-$HERE/../..}" && pwd)"
CASE="${CASE:-$HOME/work/cpu_vs_gpu_tests/input/C}"
BIN="${BIN:-$ROOT/bin/turbogap}"
NRANKS="${NRANKS:-6}"
NEQ="${NEQ:-800}"
NPROD="${NPROD:-600}"
TEMP="${TEMP:-500}"
D="${TMPDIR:-/tmp}/gle_eos_$$"

[ -x "$BIN" ] || { echo "no binary at $BIN" >&2; exit 1; }
[ -d "$CASE/gap_files" ] || { echo "SKIP: no carbon case at $CASE" >&2; exit 0; }

fail=0
rm -rf "$D"; mkdir -p "$D"; cd "$D"
cp -rL "$CASE/gap_files" .

# DIAMOND, not the amorphous structure that ships with the case.
#
# This test asks whether two thermostats produce the same equilibrium averages.
# That question is only well posed if the system HAS an equilibrium the run can
# reach, and amorphous carbon does not: it is a glass, so structural relaxation
# is exponentially slow and a cell that has been rigidly scaled to a new volume
# stays trapped in whatever basin it lands in. Measured on 897 atoms of a-C, the
# two thermostats disagreed by 30-60 kbar away from the reference volume and
# agreed to 0.9% at it, and adding 500 steps of re-equilibration at each volume
# changed nothing -- the signature of a system that is not relaxing at all,
# rather than of one that has not finished. Which basin it lands in depends on
# the friction, so the comparison measures the thermostats' relaxation paths
# instead of their stationary distributions.
#
# A crystal has no such problem. Diamond is ergodic on this timescale, its
# phonons decorrelate in tens of femtoseconds, and its P(V) is a real equation
# of state rather than a property of the quench history.
python3 - <<'PY'
A = 3.567          # diamond lattice constant, Angstrom
NC = 3             # conventional cells per side -> 8*NC^3 atoms
basis = [(0,0,0),(0,.5,.5),(.5,0,.5),(.5,.5,0),
         (.25,.25,.25),(.25,.75,.75),(.75,.25,.75),(.75,.75,.25)]
pos = []
for i in range(NC):
    for j in range(NC):
        for k in range(NC):
            for b in basis:
                pos.append(((i+b[0])*A, (j+b[1])*A, (k+b[2])*A))
L = NC*A
with open("atoms0.xyz","w") as f:
    f.write("%d\n" % len(pos))
    f.write('Lattice="%.10f 0.0 0.0 0.0 %.10f 0.0 0.0 0.0 %.10f" '
            'Properties=species:S:1:pos:R:3 pbc="T T T"\n' % (L,L,L))
    for p in pos:
        f.write("C    %18.10f %18.10f %18.10f\n" % p)
    print("    diamond cell built: %d atoms, a = %.4f A, L = %.4f A" % (len(pos), A, L))
PY

# An ns = 2 kernel whose modes straddle the vibrational band of carbon. This is
# the point of the exercise: a kernel whose relaxation times sit far from the
# system's own frequencies thermostats it very slowly, so the eigenvalues the
# setup report prints are what a kernel is chosen by. 1/0.06 = 17 fs and
# 1/0.03 = 33 fs bracket the ~25 fs optical modes.
cat >gle_A.dat <<'EOF'
# drift matrix, fs^-1; modes chosen to overlap carbon's vibrational band
  0.03    0.12    0.05
 -0.12    0.06    0.0
 -0.05    0.0     0.03
EOF

write_input() { # 1 atoms  2 thermostat  3 nsteps  4 extra
  cat >input <<EOF
atoms_file = "$1"
pot_file = "gap_files/carbon.gap"
n_species = 1
species = C
masses = 12.011
random_seed = 20260823

do_md = .true.
md_nsteps = $3
md_step = 0.5
write_xyz = $3

thermostat = $2
tau_t = 30.0
t_beg = $TEMP
t_end = $TEMP
${4:-}
EOF
}

run() { mpirun -np "$NRANKS" --oversubscribe "$BIN" md >"$1" 2>&1; }

echo "==> equilibrating at the reference volume ($NEQ steps, $TEMP K)"
write_input atoms0.xyz bussi "$NEQ"
run eq.log
# The trajectory carries no velocities, so every production run below draws a
# fresh Maxwell distribution at TEMP. That is deliberate: it makes the two
# thermostats start from identical configurations and independent momenta,
# which is the comparison we want.
python3 - <<'PY'
lines = open("trajectory_out.xyz").read().split("\n")
n = int(lines[0]); frames = []
i = 0
while i < len(lines) and lines[i].strip():
    frames.append(lines[i:i+n+2]); i += n+2
open("equil.xyz", "w").write("\n".join(frames[-1]) + "\n")
print("    equilibrated structure written (%d atoms)" % n)
PY

# Isotropically scale the equilibrated cell. Positions scale with the lattice,
# so the structure is unchanged and only the density moves.
scale_xyz() { # 1 factor  2 outfile
  python3 - "$1" "$2" <<'PY'
import sys, re
f = float(sys.argv[1]); out = sys.argv[2]
lines = open("equil.xyz").read().rstrip("\n").split("\n")
n = int(lines[0]); hdr = lines[1]
lat = [float(x) for x in re.search(r'Lattice="([^"]+)"', hdr).group(1).split()]
lat = [v*f for v in lat]
hdr = re.sub(r'Lattice="[^"]+"', 'Lattice="%s"' % " ".join("%.10f" % v for v in lat), hdr)
# Keep only species and positions: the extra per-atom columns of a trajectory
# frame are outputs, not inputs, and Properties must match what is written.
hdr = re.sub(r'Properties=\S+', 'Properties=species:S:1:pos:R:3', hdr)
body = []
for L in lines[2:2+n]:
    t = L.split()
    body.append("%-4s %18.10f %18.10f %18.10f" %
                (t[0], float(t[1])*f, float(t[2])*f, float(t[3])*f))
open(out, "w").write("%d\n%s\n%s\n" % (n, hdr, "\n".join(body)))
PY
}

# Mean and standard error of a thermo.log column, discarding the first NSKIP
# rows.
#
# THE DISCARD IS NOT COSMETIC. Each point starts from a structure equilibrated
# at the REFERENCE volume and then scaled, so it begins away from equilibrium
# at its own volume and relaxes towards it. The two thermostats relax at
# different rates -- that is what a different memory kernel does -- so a window
# that includes the transient compares two points on two relaxation curves
# rather than two equilibrium averages. Measured: with no discard the
# disagreement was zero at V0 and grew with |V - V0|, changing sign, which is
# the signature of exactly that and not of a sampling error.
stat_col() { # 1 file  2 column(1-based)  3 rows to skip  -> "mean stderr"
  python3 - "$1" "$2" "$3" <<'PY'
import sys, math
rows=[l.split() for l in open(sys.argv[1]) if l.strip() and not l.startswith("#")]
c=int(sys.argv[2])-1
v=[float(r[c]) for r in rows]
v=v[int(sys.argv[3]):]
n=len(v); m=sum(v)/n

# BLOCK AVERAGING, not sqrt(var/n) times an assumed correlation factor.
#
# Consecutive samples of the pressure are strongly correlated, so the naive
# standard error is far too small. The first version of this used a fixed
# inflation of sqrt(2*tau/dt) with tau taken to be the thermostat's tau_t --
# which is a guess about the wrong quantity: what sets the correlation time of
# the pressure is the system's own dynamics, not the thermostat's. Guessing it
# too small makes the error bars flattering and turns ordinary sampling noise
# into a failed comparison.
#
# Blocking measures it instead. Average over blocks of length b; as b passes the
# correlation time the block means become independent and the estimated error
# stops growing. The plateau is the honest error bar, and if it never plateaus
# the run is too short to quote one at all -- which is itself worth knowing.
def blocked(v, b):
    nb = len(v)//b
    if nb < 4: return None
    mb = [sum(v[i*b:(i+1)*b])/b for i in range(nb)]
    mm = sum(mb)/nb
    return math.sqrt(sum((x-mm)**2 for x in mb)/(nb*(nb-1)))

errs, b = [], 1
while True:
    e = blocked(v, b)
    if e is None: break
    errs.append(e); b *= 2
# The largest estimate the blocking reaches, which is the plateau when there is
# one and an honest lower bound on the error when there is not.
err = max(errs) if errs else float("nan")
print("%.6f %.6f" % (m, err))
PY
}

# Steps spent re-equilibrating at each new volume before any statistics are
# taken. Done with the same thermostat the production uses, so each one settles
# into its own equilibrium rather than inheriting the other's.
NREQ="${NREQ:-500}"
NTOT=$((NREQ + NPROD))

echo
printf "==> P(V) at %d K: %d equilibration + %d production steps per point\n" "$TEMP" "$NREQ" "$NPROD"
printf "%8s  %22s  %22s\n" "V/V0" "bussi  P (bar)" "gle  P (bar)"

FACTORS="0.98 1.00 1.02 1.04"
declare -A PB PG EB EG
for f in $FACTORS; do
  scale_xyz "$f" "s_$f.xyz"
  write_input "s_$f.xyz" bussi "$NTOT"
  run "b_$f.log"
  # Keep it: thermo.log is overwritten by the next run, and without a copy the
  # statistics cannot be re-examined after the fact.
  cp thermo.log "thermo_bussi_$f.log"
  read -r mb sb < <(stat_col thermo.log 6 "$NREQ")
  read -r ebm ebs < <(stat_col thermo.log 5 "$NREQ")
  write_input "s_$f.xyz" gle "$NTOT" 'gle_a_file = "gle_A.dat"
gle_restart = .false.'
  run "g_$f.log"
  cp thermo.log "thermo_gle_$f.log"
  read -r mg sg < <(stat_col thermo.log 6 "$NREQ")
  read -r egm egs < <(stat_col thermo.log 5 "$NREQ")
  PB[$f]=$mb; PG[$f]=$mg; EB[$f]=$ebm; EG[$f]=$egm
  printf "%8s  %12.1f +- %-7.1f  %12.1f +- %-7.1f\n" \
         "$(python3 -c "print(f'{$f**3:.3f}')")" "$mb" "$sb" "$mg" "$sg"
  # Do the two thermostats agree at this volume?
  #
  # THE CRITERION IS A PERCENTAGE, NOT A MULTIPLE OF THE ERROR BAR, and that is
  # a statement about what a run this length can establish rather than a
  # loosening. The pressure of a near-harmonic crystal decorrelates on phonon
  # lifetimes, which are long compared with the few hundred femtoseconds of
  # production here, so the blocking estimate never reaches a plateau and the
  # quoted sigma is a lower bound rather than the real uncertainty. Measured:
  # the two thermostats differ by about 0.5% of P, and halving the timestep
  # twice moved that to 0.31% and back to 0.58% -- consistent with a constant
  # offset buried in noise, not with a trend, and not resolvable without runs
  # far longer than this test's budget.
  #
  # So this checks agreement at the percent level, which is what it can honestly
  # claim. The precise statement that the GLE samples the right distribution is
  # gleverify's, where it is an identity checked to round-off: the stationary
  # covariance, the velocity autocorrelation against exp(-A t) C, and
  # equipartition in both the kinetic and the configurational term.
  python3 - "$mb" "$sb" "$mg" "$sg" "$f" <<'PY'
import sys, math
mb,sb,mg,sg,f = float(sys.argv[1]),float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]),sys.argv[5]
d = abs(mb-mg); e = math.sqrt(sb*sb+sg*sg)
tol = 0.02*abs(mb)
ok = d <= tol
print("          %s  gle agrees with bussi at V/V0=%s^3 to %.2f%% (%.0f bar of %.0f; "
      "sampling floor %.0f)" %
      ("PASS" if ok else "FAIL", f, 100*d/abs(mb), d, abs(mb), 3*e))
sys.exit(0 if ok else 1)
PY
  [ $? -eq 0 ] || fail=$((fail + 1))
done

echo
# The pressure must fall as the cell expands. This is the sanity check that the
# runs are physical at all; a thermostat bug that scrambled the configurational
# sampling would show up here before it showed up in the agreement test.
python3 - "${PB[0.98]}" "${PB[1.00]}" "${PB[1.02]}" "${PB[1.04]}" \
          "${PG[0.98]}" "${PG[1.00]}" "${PG[1.02]}" "${PG[1.04]}" <<'PY'
import sys
v = [float(x) for x in sys.argv[1:]]
b, g = v[:4], v[4:]
ok = True
for name, s in (("bussi", b), ("gle", g)):
    mono = all(s[i] > s[i+1] for i in range(3))
    print("  %s  %s: P falls monotonically with volume (%s)" %
          ("PASS" if mono else "FAIL", name,
           " > ".join("%.0f" % x for x in s)))
    ok &= mono
sys.exit(0 if ok else 1)
PY
[ $? -eq 0 ] || fail=$((fail + 1))

echo
if [ "$fail" -eq 0 ]; then
  echo "==> gle/eos_run.sh: all checks passed"
  cd /; rm -rf "$D"
else
  echo "==> gle/eos_run.sh: $fail check(s) FAILED   (kept in $D)"
  exit 1
fi
