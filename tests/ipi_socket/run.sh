#!/usr/bin/env bash
#
# Validation of the i-PI driver, src/ipi_socket.f90 and src/ipi_driver.f90.
#
# In `turbogap ipi` the coordinates arrive over a socket in bohr with a cell
# matrix whose index order is i-PI's, and the energy, forces and virial go back
# in Hartree. Every one of those is a convention that can be got wrong in a way
# that still runs: a transposed cell, a row/column swap in the positions, a
# missing unit factor, forces returned for the previous geometry. None of them
# raises an error. They produce a trajectory that is merely wrong.
#
# So the check is that the socket path and the file path give the SAME forces
# for the SAME configuration, where the file path is `turbogap predict`, which
# the regression suite already pins.
#
#   1. AGREEMENT. Two geometries, driven over the socket by a reference server
#      that is not this code -- tests/ipi_socket/refserver.py, written from the
#      wire protocol -- against `turbogap predict` on the same two. Energy,
#      forces and virial must agree to the precision trajectory_out.xyz is
#      written with, and nothing looser.
#
#      Two geometries and not one: a driver that answered the first exchange
#      correctly and then lost a step would pass a single-shot test. The second
#      geometry is displaced enough to move every force and to exercise the
#      Verlet-skin accounting between exchanges.
#
#   2. DISCRIMINATION. The test cell is TRICLINIC, and the run also asks what
#      `predict` gives for the same atoms in the TRANSPOSED cell. If those two
#      agreed, the cell convention would not be under test at all and check 1
#      would pass with the index order reversed. The run fails if they agree.
#      A cubic cell -- the obvious thing to test with -- has exactly that
#      defect, which is why it is not used here.
#
#   3. HANDSHAKE ORDER. refserver.py asserts NEEDINIT before INIT, READY before
#      each geometry, and HAVEDATA only after one has been sent. The last is
#      the one that matters: a client answering HAVEDATA on its first STATUS
#      hands back forces for the configuration in atoms_file, one exchange out
#      of step with i-PI forever after, and nothing downstream can tell.
#
#   4. THE SUPERCELL. Everything above runs a second time in a cell too small
#      to hold one cutoff sphere, so TurboGAP replicates it and `indices` is no
#      longer 1. That case has its own convention: a_box, b_box and c_box are
#      then the SUPERCELL vectors, while i-PI's cell is always the primitive
#      one. Assigning i-PI's cell straight into them cost hours here -- the
#      periodic images land at half spacing on top of the real atoms, and the
#      run does not fail, it grinds and then dies somewhere inside the GAP.
#      Phase 1 cannot see it, because there indices IS 1 and the two
#      conventions coincide.
#
#      The phase asserts that a supercell was actually built. Without that it
#      would quietly degenerate into a second copy of phase 1 the moment the
#      cutoff or the test cell changed.
#
# Usage:  ./run.sh [path/to/turbogap_root]

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "${1:-$HERE/../..}" && pwd)"
repo="$ROOT"
# shellcheck source=/dev/null
. "$ROOT/tests/data_root.sh"

TG="${TURBOGAP_BIN:-$ROOT/bin/turbogap}"
PY="${TURBOGAP_PYTHON:-python3}"
CO="$DATA_ROOT/CO"

[ -x "$TG" ] || { echo "ERROR: no turbogap binary at $TG"; exit 2; }
[ -d "$CO/gap_files" ] || { echo "ERROR: no CO potential under $CO"; exit 2; }
if ! "$PY" -c "import numpy" 2>/dev/null; then
    echo "SKIP: $PY has no numpy; set TURBOGAP_PYTHON to one that does"
    exit 0
fi

WORK="$(mktemp -d "${TMPDIR:-/tmp}/tg_ipi.XXXXXX")"
SOCK="tgtest$$"
cleanup() { rm -rf "$WORK"; rm -f "/tmp/ipi_$SOCK"* ; }
trap cleanup EXIT
cd "$WORK"
ln -s "$CO/gap_files" gap_files

cat > input_common <<IN
pot_file = "gap_files/CO.gap"
n_species = 2
species = C O
masses = 12.01 15.99
e0 = -.16138053 0.
random_seed = 12345
IN

echo "== building the test configurations =="
"$PY" "$HERE/make_configs.py" "$CO/atoms_897.xyz"

# ---------------------------------------------------------------- one phase
# $1 = tag, $2.. = the geometries to drive over the socket
run_phase () {
    local tag=$1; shift
    local addr="${SOCK}_$tag"
    echo "== $tag: driving turbogap ipi from the reference server =="
    { cat input_common; printf 'atoms_file = "%s_0.xyz"\nipi_address = %s\n' \
        "$tag" "'UNIX:$addr'"; } > input
    "$PY" "$HERE/refserver.py" "$addr" "$@" \
          --out "sock_$tag" --emit-effective "eff_$tag" > "server_$tag.log" 2>&1 &
    local srv=$!
    # The client does not retry, so the socket has to exist before it starts.
    for _ in $(seq 1 100); do [ -S "/tmp/ipi_$addr" ] && break; sleep 0.1; done
    [ -S "/tmp/ipi_$addr" ] || { echo "ERROR: the reference server never bound its socket"
                                 cat "server_$tag.log"; exit 1; }
    if ! "$TG" ipi > "ipi_$tag.log" 2>&1; then
        echo "ERROR: turbogap ipi failed in phase $tag"; tail -30 "ipi_$tag.log"
        cat "server_$tag.log"; exit 1
    fi
    if ! wait "$srv"; then
        echo "ERROR: the reference server reported a protocol violation in phase $tag"
        cat "server_$tag.log"; exit 1
    fi
    sed -n 's/^/  /p' "server_$tag.log"

    echo "== $tag: the same geometries through turbogap predict =="
    local n=0 f
    : > "eff_${tag}_all.xyz"
    for f in "$@"; do cat "eff_${tag}_${n}.xyz" >> "eff_${tag}_all.xyz"; n=$((n + 1)); done
    { cat input_common; printf 'atoms_file = "eff_%s_all.xyz"\nwrite_xyz = 1\n' "$tag"; } > input
    "$TG" predict > "predict_$tag.log" 2>&1 || {
        echo "ERROR: predict failed in phase $tag"; tail -20 "predict_$tag.log"; exit 1; }
    mv trajectory_out.xyz "ref_$tag.xyz"
}

run_phase big  big_0.xyz big_1.xyz

echo "== the discrimination control: the same atoms in the transposed cell =="
"$PY" - <<'PYEOF'
import numpy as np
src = open("eff_big_0.xyz").read().splitlines()
key = 'Lattice="'
i = src[1].index(key) + len(key)
j = src[1].index('"', i)
cell = np.array([float(x) for x in src[1][i:j].split()]).reshape(3, 3)
src[1] = src[1][:i] + " ".join(f"{v:.17g}" for v in cell.T.ravel()) + src[1][j:]
open("big_0T.xyz", "w").write("\n".join(src) + "\n")
PYEOF
{ cat input_common; printf 'atoms_file = "big_0T.xyz"\nwrite_xyz = 1\n'; } > input
"$TG" predict > predictT.log 2>&1 || { echo "ERROR: transposed predict failed"; tail -20 predictT.log; exit 1; }
mv trajectory_out.xyz refT.xyz

run_phase small small_0.xyz small_1.xyz

echo "== confirming the small cell really did build a supercell =="
if ! grep -q "smaller" "ipi_small.log"; then
    echo "  FAIL: no supercell was built, so the small-cell phase is only a"
    echo "        second copy of the large-cell one. Shrink the test cell."
    exit 1
fi
if grep -q "smaller" "ipi_big.log"; then
    echo "  FAIL: the large-cell phase built a supercell too, so the two"
    echo "        phases no longer test different code paths."
    exit 1
fi
echo "  ok: supercell in the small-cell phase, none in the large-cell phase"

echo "== comparing =="
"$PY" "$HERE/compare.py"

echo "PASS: the socket path and the file path agree, with and without a supercell"
