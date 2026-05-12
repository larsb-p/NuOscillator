#!/bin/bash
#SBATCH --job-name=osc_probs
#SBATCH --account=def-blairt2k
#SBATCH --time=0-23:59:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=16
#SBATCH --array=0-6
#SBATCH --output=logs/osc_probs_%A_%a.out
#SBATCH --error=logs/osc_probs_%A_%a.err

set -eu

echo "===== TASK $SLURM_ARRAY_TASK_ID START ====="
date

source /home/lbp/links/scratch/plotting/NuOscillator/build/Linux/bin/setup.NuOscillator.sh
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# --- 26 parameter sets ---
# Names line up with params; index by SLURM_ARRAY_TASK_ID.
names=(
  setA_n-3
  setA_n-2
  setA_n-1
  setA_n0
  setA_n1
  setA_n2
  setA_n3  
)

params=(
  # 0: setA_n-3
  "sin2theta12=0.307 sin2theta23=0.561 sin2theta13=0.022 dm21=7.53e-5 dm32=2.494e-3 dcp=-1.601 rho=2.6 Ye=0.5 gamma0=2.0e-15 n=-3.0 E0=1.0"
  # 1: setA_n-2
  "sin2theta12=0.307 sin2theta23=0.561 sin2theta13=0.022 dm21=7.53e-5 dm32=2.494e-3 dcp=-1.601 rho=2.6 Ye=0.5 gamma0=5.0e-15 n=-2.0 E0=1.0"
  # 2: setA_n-1
  "sin2theta12=0.307 sin2theta23=0.561 sin2theta13=0.022 dm21=7.53e-5 dm32=2.494e-3 dcp=-1.601 rho=2.6 Ye=0.5 gamma0=3.0e-14 n=-1.0 E0=1.0"
  # 3: setA_n0
  "sin2theta12=0.307 sin2theta23=0.561 sin2theta13=0.022 dm21=7.53e-5 dm32=2.494e-3 dcp=-1.601 rho=2.6 Ye=0.5 gamma0=2.0e-14 n=0.0 E0=1.0"
  # 4: setA_n1
  "sin2theta12=0.307 sin2theta23=0.561 sin2theta13=0.022 dm21=7.53e-5 dm32=2.494e-3 dcp=-1.601 rho=2.6 Ye=0.5 gamma0=5.0e-14 n=1.0 E0=1.0"
  # 5: setA_n2
  "sin2theta12=0.307 sin2theta23=0.561 sin2theta13=0.022 dm21=7.53e-5 dm32=2.494e-3 dcp=-1.601 rho=2.6 Ye=0.5 gamma0=1.0e-14 n=2.0 E0=1.0"
  # 6: setA_n3
  "sin2theta12=0.307 sin2theta23=0.561 sin2theta13=0.022 dm21=7.53e-5 dm32=2.494e-3 dcp=-1.601 rho=2.6 Ye=0.5 gamma0=1.0e-14 n=3.0 E0=1.0"
)

NAME=${names[$SLURM_ARRAY_TASK_ID]}
PARAMS=${params[$SLURM_ARRAY_TASK_ID]}

OUTBASE=/home/lbp/links/scratch/plotting/NuOscillator/output_t2k_energyscale_variations
OUTDIR=$OUTBASE/$NAME
mkdir -p "$OUTDIR" && cd "$OUTDIR"

CONFIG=/home/lbp/links/scratch/plotting/NuOscillator/NuOscillatorConfigs/Unbinned_NuSQUIDSLinear.yaml

echo "Task $SLURM_ARRAY_TASK_ID ($NAME) in $OUTDIR"
echo "Params: $PARAMS"

SingleOscProbCalcer "$CONFIG" $PARAMS

echo "===== TASK $SLURM_ARRAY_TASK_ID END ====="
date
