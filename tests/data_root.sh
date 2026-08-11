# Where the test systems live. Sourced by every suite's run.sh.
#
# The systems are not in the source tree -- they are trajectories and fitted
# potentials, shared by the CPU and GPU branches and changing on a different
# cadence from the code, so they have a repository of their own:
#
#   https://github.com/TiganyZ/turbogap_tests
#
# tests/fetch_test_data.sh clones it, and this file calls that when the clone
# is not there yet, so a first `run.sh` on a fresh checkout works.
#
# The caller must have set `repo` (the top of the source tree) beforehand.
# Sets DATA_ROOT.
#
#   TURBOGAP_DATA_ROOT   take the systems from here and fetch nothing. For a
#                        checkout of the data at some other path, or a subset
#                        assembled by hand.
#   TURBOGAP_TESTS_DIR   where the data repository is, or should be, cloned.
#                        Default: <repo>/../turbogap_tests, the same default
#                        fetch_test_data.sh uses, so the two agree without
#                        being told anything.

DATA_ROOT=${TURBOGAP_DATA_ROOT:-${TURBOGAP_TESTS_DIR:-$repo/../turbogap_tests}}

if [ ! -d "$DATA_ROOT" ] && [ -z "${TURBOGAP_DATA_ROOT:-}" ]; then
  "$repo/tests/fetch_test_data.sh" "$DATA_ROOT" || {
    printf 'ERROR: could not fetch the test data into %s\n' "$DATA_ROOT" >&2
    exit 2
  }
fi
