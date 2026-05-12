#!/bin/bash
#SBATCH --job-name=osc_probs
#SBATCH --account=def-blairt2k
#SBATCH --time=0-23:59:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=16
#SBATCH --array=0-33
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
  setA
  setB
  setA_g0
  setB_g0
  setA_g0_dcp0
  setA_g0_dcp_pos
  setA_g0_th13_low
  setA_g0_th23_max
  setA_g0_dm32_low
  setA_g0_dm32_high
  setA_g0_th12_low
  setA_g0_th12_high
  setA_g0_dm21_low
  setA_g0_dm21_high
  setA_g0_dcp_opp
  setB_g0_dcp_neg
  setB_g0_dcp_pos
  setB_g0_th13_low
  setB_g0_th23_max
  setB_g0_dm32_low
  setB_g0_dm32_high
  setB_g0_th12_low
  setB_g0_th12_high
  setB_g0_dm21_low
  setB_g0_dm21_high
  setB_g0_dcp_opp
  setA_g0_r0
  setA_r0
  setB_g0_r0
  setB_r0
  setA_all_osc_params_zero
  setA_all_osc_params_zero_vacuum
  set_all_osc_params_zero
  set_all_osc_params_zero_vacuum
)

params=(
  # 0: setA — T2K Set A
  "sin2theta12=0.307 sin2theta23=0.561 sin2theta13=0.022 dm21=7.53e-5 dm32=2.494e-3 dcp=-1.601 rho=2.6 Ye=0.5 gamma0=5.0e-15 n=-2.0 E0=1.0"
  # 1: setB — T2K Set B
  "sin2theta12=0.307 sin2theta23=0.450 sin2theta13=0.022 dm21=7.53e-5 dm32=2.494e-3 dcp=0.0    rho=2.6 Ye=0.5 gamma0=5.0e-15 n=-2.0 E0=1.0"
  # 2: setA_g0
  "sin2theta12=0.307 sin2theta23=0.561 sin2theta13=0.022 dm21=7.53e-5 dm32=2.494e-3 dcp=-1.601 rho=2.6 Ye=0.5 gamma0=0.0     n=-2.0 E0=1.0"
  # 3: setB_g0
  "sin2theta12=0.307 sin2theta23=0.450 sin2theta13=0.022 dm21=7.53e-5 dm32=2.494e-3 dcp=0.0    rho=2.6 Ye=0.5 gamma0=0.0     n=-2.0 E0=1.0"
  # 4: setA_g0_dcp0
  "sin2theta12=0.307 sin2theta23=0.561 sin2theta13=0.022 dm21=7.53e-5 dm32=2.494e-3 dcp=0.0    rho=2.6 Ye=0.5 gamma0=0.0     n=-2.0 E0=1.0"
  # 5: setA_g0_dcp_pos
  "sin2theta12=0.307 sin2theta23=0.561 sin2theta13=0.022 dm21=7.53e-5 dm32=2.494e-3 dcp=1.601  rho=2.6 Ye=0.5 gamma0=0.0     n=-2.0 E0=1.0"
  # 6: setA_g0_th13_low
  "sin2theta12=0.307 sin2theta23=0.561 sin2theta13=0.018  dm21=7.53e-5 dm32=2.494e-3 dcp=-1.601 rho=2.6 Ye=0.5 gamma0=0.0     n=-2.0 E0=1.0"
  # 7: setA_g0_th23_max
  "sin2theta12=0.307 sin2theta23=0.5 sin2theta13=0.022 dm21=7.53e-5 dm32=2.494e-3 dcp=-1.601 rho=2.6 Ye=0.5 gamma0=0.0     n=-2.0 E0=1.0"
  # 8: setA_g0_dm32_low
  "sin2theta12=0.307 sin2theta23=0.561 sin2theta13=0.022 dm21=7.53e-5 dm32=2.4e-3   dcp=-1.601 rho=2.6 Ye=0.5 gamma0=0.0     n=-2.0 E0=1.0"
  # 9: setA_g0_dm32_high
  "sin2theta12=0.307 sin2theta23=0.561 sin2theta13=0.022 dm21=7.53e-5 dm32=2.6e-3   dcp=-1.601 rho=2.6 Ye=0.5 gamma0=0.0     n=-2.0 E0=1.0"
  # 10: setA_g0_th12_low
  "sin2theta12=0.293 sin2theta23=0.561 sin2theta13=0.022 dm21=7.53e-5 dm32=2.494e-3 dcp=-1.601 rho=2.6 Ye=0.5 gamma0=0.0     n=-2.0 E0=1.0"
  # 11: setA_g0_th12_high
  "sin2theta12=0.321 sin2theta23=0.561 sin2theta13=0.022 dm21=7.53e-5 dm32=2.494e-3 dcp=-1.601 rho=2.6 Ye=0.5 gamma0=0.0     n=-2.0 E0=1.0"
  # 12: setA_g0_dm21_low
  "sin2theta12=0.307 sin2theta23=0.561 sin2theta13=0.022 dm21=7.0e-5  dm32=2.494e-3 dcp=-1.601 rho=2.6 Ye=0.5 gamma0=0.0     n=-2.0 E0=1.0"
  # 13: setA_g0_dm21_high
  "sin2theta12=0.307 sin2theta23=0.561 sin2theta13=0.022 dm21=8.0e-5  dm32=2.494e-3 dcp=-1.601 rho=2.6 Ye=0.5 gamma0=0.0     n=-2.0 E0=1.0"
  # 14: setA_g0_dcp_opp
  "sin2theta12=0.307 sin2theta23=0.561 sin2theta13=0.022 dm21=7.53e-5 dm32=2.494e-3 dcp=3.141592654 rho=2.6 Ye=0.5 gamma0=0.0 n=-2.0 E0=1.0"
  # 15: setB_g0_dcp_neg
  "sin2theta12=0.307 sin2theta23=0.450 sin2theta13=0.022 dm21=7.53e-5 dm32=2.494e-3 dcp=-1.601 rho=2.6 Ye=0.5 gamma0=0.0     n=-2.0 E0=1.0"
  # 16: setB_g0_dcp_pos
  "sin2theta12=0.307 sin2theta23=0.450 sin2theta13=0.022 dm21=7.53e-5 dm32=2.494e-3 dcp=1.601  rho=2.6 Ye=0.5 gamma0=0.0     n=-2.0 E0=1.0"
  # 17: setB_g0_th13_low
  "sin2theta12=0.307 sin2theta23=0.450 sin2theta13=0.018  dm21=7.53e-5 dm32=2.494e-3 dcp=0.0    rho=2.6 Ye=0.5 gamma0=0.0     n=-2.0 E0=1.0"
  # 18: setB_g0_th23_max
  "sin2theta12=0.307 sin2theta23=0.5 sin2theta13=0.022 dm21=7.53e-5 dm32=2.494e-3 dcp=0.0    rho=2.6 Ye=0.5 gamma0=0.0     n=-2.0 E0=1.0"
  # 19: setB_g0_dm32_low
  "sin2theta12=0.307 sin2theta23=0.450 sin2theta13=0.022 dm21=7.53e-5 dm32=2.4e-3   dcp=0.0    rho=2.6 Ye=0.5 gamma0=0.0     n=-2.0 E0=1.0"
  # 20: setB_g0_dm32_high
  "sin2theta12=0.307 sin2theta23=0.450 sin2theta13=0.022 dm21=7.53e-5 dm32=2.6e-3   dcp=0.0    rho=2.6 Ye=0.5 gamma0=0.0     n=-2.0 E0=1.0"
  # 21: setB_g0_th12_low
  "sin2theta12=0.293 sin2theta23=0.450 sin2theta13=0.022 dm21=7.53e-5 dm32=2.494e-3 dcp=0.0    rho=2.6 Ye=0.5 gamma0=0.0     n=-2.0 E0=1.0"
  # 22: setB_g0_th12_high
  "sin2theta12=0.321 sin2theta23=0.450 sin2theta13=0.022 dm21=7.53e-5 dm32=2.494e-3 dcp=0.0    rho=2.6 Ye=0.5 gamma0=0.0     n=-2.0 E0=1.0"
  # 23: setB_g0_dm21_low
  "sin2theta12=0.307 sin2theta23=0.450 sin2theta13=0.022 dm21=7.0e-5  dm32=2.494e-3 dcp=0.0    rho=2.6 Ye=0.5 gamma0=0.0     n=-2.0 E0=1.0"
  # 24: setB_g0_dm21_high
  "sin2theta12=0.307 sin2theta23=0.450 sin2theta13=0.022 dm21=8.0e-5  dm32=2.494e-3 dcp=0.0    rho=2.6 Ye=0.5 gamma0=0.0     n=-2.0 E0=1.0"
  # 25: setB_g0_dcp_opp
  "sin2theta12=0.307 sin2theta23=0.450 sin2theta13=0.022 dm21=7.53e-5 dm32=2.494e-3 dcp=3.141592654 rho=2.6 Ye=0.5 gamma0=0.0 n=-2.0 E0=1.0"
  # 26: setA_g0_r0
  "sin2theta12=0.307 sin2theta23=0.561 sin2theta13=0.022 dm21=7.53e-5 dm32=2.494e-3 dcp=-1.601 rho=0.0 Ye=0.5 gamma0=0.0 n=-2.0 E0=1.0"
  # 27: setA_r0
  "sin2theta12=0.307 sin2theta23=0.561 sin2theta13=0.022 dm21=7.53e-5 dm32=2.494e-3 dcp=-1.601 rho=0.0 Ye=0.5 gamma0=5.0e-15 n=-2.0 E0=1.0"
  # 28: setB_g0_r0
  "sin2theta12=0.307 sin2theta23=0.450 sin2theta13=0.0220 dm21=7.53e-5 dm32=2.494e-3 dcp=0.0 rho=0.0 Ye=0.5 gamma0=0.0 n=-2.0 E0=1.0"
  # 29: setB_r0
  "sin2theta12=0.307 sin2theta23=0.450 sin2theta13=0.0220 dm21=7.53e-5 dm32=2.494e-3 dcp=0.0 rho=0.0 Ye=0.5 gamma0=5.0e-15 n=-2.0 E0=1.0"
  # 30: setA_all_osc_params_zero
  "sin2theta12=0.0 sin2theta23=0.0 sin2theta13=0.0 dm21=0.0 dm32=0.0 dcp=0.0 rho=2.6 Ye=0.5 gamma0=5.0e-15 n=-2.0 E0=1.0"
  # 31: setA_all_osc_params_zero_vacuum
  "sin2theta12=0.0 sin2theta23=0.0 sin2theta13=0.0 dm21=0.0 dm32=0.0 dcp=0.0 rho=0.0 Ye=0.5 gamma0=5.0e-15 n=-2.0 E0=1.0"
  # 32: set_all_osc_params_zero
  "sin2theta12=0.0 sin2theta23=0.0 sin2theta13=0.0 dm21=0.0 dm32=0.0 dcp=0.0 rho=2.6 Ye=0.5 gamma0=0.0 n=-2.0 E0=1.0" 
  # 33: set_all_osc_params_zero_vacuum
  "sin2theta12=0.0 sin2theta23=0.0 sin2theta13=0.0 dm21=0.0 dm32=0.0 dcp=0.0 rho=0.0 Ye=0.5 gamma0=0.0 n=-2.0 E0=1.0"
)

NAME=${names[$SLURM_ARRAY_TASK_ID]}
PARAMS=${params[$SLURM_ARRAY_TASK_ID]}

OUTBASE=/home/lbp/links/scratch/plotting/NuOscillator/output
OUTDIR=$OUTBASE/$NAME
mkdir -p "$OUTDIR" && cd "$OUTDIR"

CONFIG=/home/lbp/links/scratch/plotting/NuOscillator/NuOscillatorConfigs/Unbinned_NuSQUIDSLinear.yaml

echo "Task $SLURM_ARRAY_TASK_ID ($NAME) in $OUTDIR"
echo "Params: $PARAMS"

SingleOscProbCalcer "$CONFIG" $PARAMS

echo "===== TASK $SLURM_ARRAY_TASK_ID END ====="
date
