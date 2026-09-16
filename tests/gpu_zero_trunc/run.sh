#!/bin/bash
# Author: Tigany Zarrouk (tiganyzarrouk@gmail.com)
#
# The device local-property predictor clamps negative values at zero, and drops
# the gradient of the sites it clamped, exactly as the host one does.
#
# Both halves need checking and only one of them is obvious. The value clamp
# shows up in the output column; the gradient zeroing shows up nowhere unless
# something uses the gradient LINEARLY. vdW-TS does not -- it builds C6 from
# V^2, so dC6/dr carries a factor of V and vanishes on its own at a clamped
# site, and a device build with the gradient zeroing deleted passes a vdW
# comparison unchanged. The XPS bias does use it linearly:
#
#     f_j -= sum_i lambda_i dq_i/dr_j
#
# so this drives the bias instead. Measured discrimination: deleting the call
# to gpu_zero_trunc_der moves the forces by 1.9 eV/A on this case, against
# 4e-8 when it is there.
#
# The clamp does not fire on a sane potential -- a binding energy is ~285 eV --
# so the case shifts V0 to make it fire on 466 of 512 sites, and puts the
# target spectrum over a window straddling zero so clamped and unclamped sites
# are both inside it. The numbers are not physical; the agreement is.
#
# Host against device in the same session, never against a stored number: the
# device does not reproduce itself run to run.

set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)

PYTHON=${TURBOGAP_PYTHON:-python3}
BIN=${TURBOGAP_BIN:-$repo/bin/turbogap}
GPU_BIN=${TURBOGAP_GPU_BIN:-$repo/bin-gpu/turbogap}
WORK=${TURBOGAP_WORK:-/tmp/gpu_zero_trunc.$$}

# V0 for the carbon core-electron model. The fitted value is 290 eV; the GP
# part contributes about -3.7 to -6.8, so a V0 in that range straddles zero.
V0_SHIFTED=5.0
NSITES=512
TOL_Q=1.0e-6
TOL_F=1.0e-6

. "$repo/tests/data_root.sh"

say() { printf '\n== %s\n' "$1"; }
skip() { printf '    SKIP: %s\n' "$1"; exit 0; }
die() { printf '    FAIL: %s\n' "$1" >&2; exit 1; }

check_preconditions() {
    [ -x "$BIN" ] || skip "no host binary at $BIN"
    [ -x "$GPU_BIN" ] || skip "no device binary at $GPU_BIN (make TURBOGAP_ARCH=... BUILD_TAG_EXTRA=-gpu)"
    [ -d "$DATA_ROOT/xps_opt" ] || skip "no xps_opt data under $DATA_ROOT"
}

stage_potential() {
    rm -rf "$WORK"
    mkdir -p "$WORK"
    cp -rL "$DATA_ROOT/xps_opt/gap_files" "$WORK/gap_files"
    # First block only: carbon. Oxygen keeps its 536 eV and never clamps, so
    # the run also shows that clamping is per site rather than per call.
    sed -i "0,/^local_property_v0s = 290.0/s//local_property_v0s = $V0_SHIFTED/" \
        "$WORK/gap_files/CO.gap"
    grep -q "local_property_v0s = $V0_SHIFTED" "$WORK/gap_files/CO.gap" \
        || die "could not shift the carbon V0; has the gap file changed?"
}

stage_target() {
    "$PYTHON" - "$WORK/target.dat" <<'PY'
import math, sys
with open(sys.argv[1], "w") as f:
    for i in range(501):
        x = -5.0 + 10.0 * i / 500.0
        f.write("%.6f %.6f\n" % (x, math.exp(-((x - 0.5) / 0.6) ** 2)))
PY
}

stage_input() {
    cat > "$WORK/input" <<INPUT
atoms_file = "atoms.xyz"
pot_file = "gap_files/CO.gap"
n_species = 2
species = C O
masses = 12.01 15.99
e0 = -.16138053 0.
random_seed = 12345
xps_sigma = 0.4

n_exp = 1
exp_labels = "xps"
exp_data_files = "target.dat"
exp_n_samples = 501
exp_energy_scales = 100.0
exp_forces = .true.

md_nsteps = 1
write_xyz = 1
INPUT
}

run_one() {
    local dir=$1 bin=$2
    mkdir -p "$WORK/$dir"
    cp "$WORK/input" "$WORK/$dir/"
    ln -sf "$DATA_ROOT/xps_opt/atoms.xyz" "$WORK/$dir/"
    ln -sf "$WORK/gap_files" "$WORK/$dir/gap_files"
    ln -sf "$WORK/target.dat" "$WORK/$dir/"
    ( cd "$WORK/$dir" && "$bin" md > run.log 2>&1 ) \
        || { tail -20 "$WORK/$dir/run.log"; die "$dir run"; }
}

compare() {
    "$PYTHON" - "$WORK/cpu/trajectory_out.xyz" "$WORK/gpu/trajectory_out.xyz" \
               "$NSITES" "$TOL_Q" "$TOL_F" <<'PY'
import sys

path_h, path_d, n, tol_q, tol_f = sys.argv[1:6]
n, tol_q, tol_f = int(n), float(tol_q), float(tol_f)


def frame(path):
    body = open(path).read().split("\n")[2:2 + n]
    rows = [l.split() for l in body if len(l.split()) >= 12]
    if len(rows) != n:
        sys.exit("    FAIL: %s has %d sites, expected %d" % (path, len(rows), n))
    f = [[float(r[7]), float(r[8]), float(r[9])] for r in rows]
    q = [float(r[11]) for r in rows]
    return f, q


fh, qh = frame(path_h)
fd, qd = frame(path_d)

clamped_h = [i for i, v in enumerate(qh) if v == 0.0]
clamped_d = [i for i, v in enumerate(qd) if v == 0.0]
print("    clamped: %d of %d on the host, %d on the device" % (len(clamped_h), n, len(clamped_d)))

# Without this the rest is vacuous: two backends that both clamp nothing agree
# trivially, and so do two that clamp everything.
if not 0 < len(clamped_h) < n:
    sys.exit("    FAIL: the clamp fired on %d of %d sites; it must fire on some "
             "and not all, or nothing is tested" % (len(clamped_h), n))
if clamped_h != clamped_d:
    sys.exit("    FAIL: the two backends clamped different sites")

dq = max(abs(a - b) for a, b in zip(qh, qd))
df = max(abs(a - b) for x, y in zip(fh, fd) for a, b in zip(x, y))
mf = max(abs(a) for x in fh for a in x)
print("    core_electron_be: max |device - host| = %.3e" % dq)
print("    forces          : max |F| = %.4f eV/A, max |device - host| = %.3e eV/A (%.1e rel)"
      % (mf, df, df / mf))
if dq > tol_q:
    sys.exit("    FAIL: the clamped values disagree")
if df > tol_f:
    sys.exit("    FAIL: the forces disagree -- the gradient of a clamped site "
             "is being kept on one backend and dropped on the other")
print("    PASS: the device clamps the same sites and drops their gradients")
PY
}

main() {
    say "zero_trunc on the device: same values, same dropped gradients"
    check_preconditions
    stage_potential
    stage_target
    stage_input
    run_one cpu "$BIN"
    run_one gpu "$GPU_BIN"
    grep -q "floored at zero" "$WORK/cpu/run.log" \
        || die "the host did not report flooring; the case no longer clamps"
    grep -q "floored at zero" "$WORK/gpu/run.log" \
        || die "the device did not report flooring"
    compare
    [ -n "${TURBOGAP_KEEP:-}" ] || rm -rf "$WORK"
}

main "$@"
