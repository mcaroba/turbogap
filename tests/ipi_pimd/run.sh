#!/usr/bin/env bash
#
# Path-integral MD: i-PI driving turbogap over a socket, eight beads.
#
# tests/ipi_socket checks the wire protocol against a reference server written
# for the purpose -- units, cell convention, one exchange after another. This
# checks the thing that protocol exists for: a real i-PI, a real ring polymer,
# and eight beads each driven by its own turbogap.
#
# What it asserts, and why each one can fail on its own:
#
#   1. The run completes. Eight drivers against one server is where a driver
#      that connects but never answers, or answers out of turn, shows up.
#   2. All eight bead trajectories are written, with every atom in each.
#   3. The beads are DISTINCT. A ring polymer whose beads coincide is a
#      classical simulation wearing a quantum label -- the commonest silent
#      failure in a PIMD setup, and invisible in any per-bead check.
#   4. i-PI's conserved quantity does not run away.
#
# Needs i-PI (tools/setup_dev_env.sh installs it) and the water potential from
# the test-data repository.
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/../.." && pwd)

PYTHON=${TURBOGAP_PYTHON:-python3}
BIN=${TURBOGAP_BIN:-$repo/bin/turbogap}
IPI=${IPI_BIN:-i-pi}
BEADS=8
STEPS=6

# shellcheck source=../data_root.sh
. "$repo/tests/data_root.sh"

# A water potential that returns ENERGIES. The water_dipole model in the test
# data is a dipole model -- dipole_model = .true. on both its descriptors --
# which by construction returns none, so a ring polymer driven by it does not
# move and the dynamics are untestable. Point TURBOGAP_WATER_POT at a directory
# holding gap_files/ for a water GAP with an energy model.
DATA=${TURBOGAP_WATER_POT:-$DATA_ROOT/water_dipole}
POT=${TURBOGAP_WATER_GAP:-water_dipole.gap}

WORK=""
cleanup() {
    [ -n "$WORK" ] || return 0
    pkill -f "turbogap ipi" 2>/dev/null || true
    if [ -n "${TURBOGAP_KEEP:-}" ]; then
        printf 'staging kept in %s\n' "$WORK"
    else
        rm -rf "$WORK"
    fi
}
trap cleanup EXIT

check_inputs() {
    [ -x "$BIN" ] || { printf 'ERROR: no binary at %s\n' "$BIN" >&2; exit 2; }
    command -v "$IPI" >/dev/null || {
        printf 'SKIP: i-pi is not on PATH; run tools/setup_dev_env.sh\n'; exit 0; }
    [ -f "$DATA/gap_files/$POT" ] || {
        printf 'SKIP: no water potential at %s/gap_files/%s\n' "$DATA" "$POT"; exit 0; }
}

write_structure() {
    "$PYTHON" "$here/make_water.py" "$WORK/water.xyz" > /dev/null
    # i-PI reads its own initial configuration; the comment line it wants is a
    # cell in its own notation, not an extended-xyz Lattice.
    "$PYTHON" - "$WORK/water.xyz" "$WORK/init.xyz" <<'PYEOF'
import sys
source, target = sys.argv[1], sys.argv[2]
lines = open(source).read().splitlines()
count = int(lines[0])
body = lines[2:2 + count]
with open(target, "w") as handle:
    handle.write(f"{count}\n")
    handle.write("# CELL(abcABC): 12.0 12.0 12.0 90.0 90.0 90.0 positions{angstrom}\n")
    for row in body:
        handle.write(row + "\n")
PYEOF
}

write_ipi_input() {
    cat > "$WORK/input.xml" <<XMLEOF
<simulation verbosity='low'>
  <output prefix='pimd'>
    <properties stride='1' filename='out'>
      [ step, time{picosecond}, conserved{electronvolt},
        potential{electronvolt}, kinetic_cv{electronvolt} ]
    </properties>
    <!-- No bead= attribute: i-PI then writes one file per bead, which is what
         lets the run check that the ring polymer is actually spread. -->
    <trajectory filename='pos' stride='1' format='xyz'> positions{angstrom} </trajectory>
  </output>
  <total_steps> $STEPS </total_steps>
  <prng> <seed> 12345 </seed> </prng>
  <ffsocket name='turbogap' mode='unix' pbc='false'>
    <address> ${SOCKET} </address>
    <latency> 0.01 </latency>
  </ffsocket>
  <system>
    <initialize nbeads='$BEADS'>
      <file mode='xyz'> init.xyz </file>
      <velocities mode='thermal' units='kelvin'> 300 </velocities>
    </initialize>
    <forces> <force forcefield='turbogap'/> </forces>
    <motion mode='dynamics'>
      <dynamics mode='nvt'>
        <timestep units='femtosecond'> 0.25 </timestep>
        <thermostat mode='pile_l'> <tau units='femtosecond'> 100 </tau> </thermostat>
      </dynamics>
    </motion>
    <ensemble> <temperature units='kelvin'> 300 </temperature> </ensemble>
  </system>
</simulation>
XMLEOF
}

write_turbogap_input() {
    cat > "$WORK/input" <<INPEOF
atoms_file = "water.xyz"
pot_file = "gap_files/${POT}"
n_species = 2
species = H O
masses = 1.008 15.999
e0 = 0. 0.
ipi_address = UNIX:${SOCKET}
max_Gbytes_per_process = 1.0
INPEOF
}

start_drivers() {
    local i
    for i in $(seq 1 "$BEADS"); do
        ( cd "$WORK" && "$BIN" ipi > "driver_$i.log" 2>&1 ) &
    done
}

run_simulation() {
    ( cd "$WORK" && "$IPI" input.xml > ipi.log 2>&1 ) &
    local server=$!
    sleep 3
    start_drivers
    wait "$server"
}

check_results() {
    "$PYTHON" "$here/pimdcheck.py" "$WORK" "$BEADS" "$STEPS"
}

check_inputs
WORK=$(mktemp -d "${TMPDIR:-/tmp}/turbogap_ipi_pimd.XXXXXX")
SOCKET="tgpimd_$$"
printf 'binary %s\ni-pi   %s\nbeads  %d\n\n' "$BIN" "$(command -v "$IPI")" "$BEADS"
ln -s "$DATA/gap_files" "$WORK/gap_files"
write_structure
write_turbogap_input
write_ipi_input
run_simulation
check_results
