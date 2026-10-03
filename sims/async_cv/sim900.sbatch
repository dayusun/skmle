#!/bin/bash
#SBATCH -J acv_n900
#SBATCH -A r00898
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 00:30:00
#SBATCH --array=0-399
#SBATCH -o logs/sim900_%A_%a.out
set -euo pipefail
cd /N/slate/dayusun/skmle_asynccv
module load r/4.4.1
export R_LIBS=$PWD/lib:$HOME/r-packages R_LIBS_USER=$PWD/lib OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# 2 links x 200 chunks of 10 replicates = 2000 replicates per link at n = 900.
L=(identity logistic)
cell=$(( SLURM_ARRAY_TASK_ID / 200 )); chunk=$(( SLURM_ARRAY_TASK_ID % 200 ))
link=${L[$cell]}; first=$(( chunk * 10 + 1 ))
Rscript sim_fixed_lower.R 900 "$link" "$first" 10 "results900/${link}_n900_c${chunk}.rds"
