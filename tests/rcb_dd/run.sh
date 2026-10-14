#!/usr/bin/env bash
#
# Does dd_type = rcb give the same results as the replicated run? The
# legacy_dd legs with the grid replaced by recursive coordinate bisection,
# and the legs only rcb can take; see ../legacy_dd/run.sh.
#
# Usage: run.sh [leg ...]

set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
TURBOGAP_DD_TYPE=rcb exec "$here/../legacy_dd/run.sh" "$@"
