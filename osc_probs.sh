#!/bin/bash
# Run osc_probs.sbatch locally (no SLURM) by iterating over array indices
# in series. Each iteration runs the sbatch body as a regular bash script,
# with SLURM_ARRAY_TASK_ID set so the names[]/params[] lookup works.
#
# Update NUM_TASKS to match the --array=0-N in osc_probs.sbatch.

set -eu

SBATCH_FILE="osc_probs.sbatch"
NUM_TASKS=26      # matches --array=0-25
THREADS=8         # OMP_NUM_THREADS for each task

for i in $(seq 0 $((NUM_TASKS - 1))); do
  echo
  echo "############### local task $i / $((NUM_TASKS - 1)) ###############"
  SLURM_ARRAY_TASK_ID=$i SLURM_CPUS_PER_TASK=$THREADS bash "$SBATCH_FILE"
done
