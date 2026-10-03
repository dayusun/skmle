#!/bin/bash
#SBATCH -J acv_loss
#SBATCH -A r00898
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 00:30:00
#SBATCH --array=0-599
#SBATCH -o logs/sim4_%A_%a.out
set -euo pipefail
cd /N/slate/dayusun/skmle_asynccv
module load r/4.4.1
export R_LIBS=$PWD/lib:$HOME/r-packages R_LIBS_USER=$PWD/lib OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# 6 cells (n x link) x 100 chunks of 20 replicates = 2000 replicates per cell.
N=(100 200 400 100 200 400); L=(identity identity identity logistic logistic logistic)
cell=$(( SLURM_ARRAY_TASK_ID / 100 )); chunk=$(( SLURM_ARRAY_TASK_ID % 100 ))
n=${N[$cell]}; link=${L[$cell]}; first=$(( chunk * 20 + 1 ))
Rscript sim_acv4.R "$n" "$link" "$first" 20 "results4/${link}_n${n}_c${chunk}.rds"
