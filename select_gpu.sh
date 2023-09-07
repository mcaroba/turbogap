#!/bin/bash
GPUSID="4 5 2 3 6 7 0 1"
GPUSID=(${GPUSID})
if [ ${#GPUSID[@]} -gt 0 -a -n "${SLURM_NTASKS_PER_NODE}" ]; then
if [ ${#GPUSID[@]} -gt $SLURM_NTASKS_PER_NODE ]; then
export ROCR_VISIBLE_DEVICES=${GPUSID[$(($SLURM_LOCALID))]}
else
export ROCR_VISIBLE_DEVICES=${GPUSID[$(($SLURM_LOCALID / ($SLURM_NTASKS_PER_NODE / ${#GPUSID[@]})))]}
fi fi
exec $*