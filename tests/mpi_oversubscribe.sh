# Does mpirun need telling that more ranks than slots is on purpose?
#
# OpenMPI counts a slot per physical core, so a two-core machine -- a hosted
# CI runner, a laptop -- refuses `-np 3` outright:
#
#   There are not enough slots available in the system to satisfy the 3
#   slots that were requested by the application
#
# The multi-rank cases here are about which rank owns which atom, not about
# speed, so oversubscribing is the right answer rather than skipping them.
# MPICH has no such flag and needs none, so the probe leaves this empty there.
#
# The probe runs the flag rather than reading `mpirun --help`, which does not
# list it (OpenMPI 4.1 names it only in the error message you get without it).
#
# Sets MPI_OVERSUBSCRIBE, to be used unquoted: mpirun $MPI_OVERSUBSCRIBE -np N

MPI_OVERSUBSCRIBE=""
if command -v mpirun > /dev/null 2>&1 &&
   mpirun --oversubscribe -np 1 true > /dev/null 2>&1; then
  MPI_OVERSUBSCRIBE="--oversubscribe"
fi
