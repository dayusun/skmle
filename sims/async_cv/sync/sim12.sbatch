#!/bin/bash
#SBATCH -J log12
#SBATCH -A r00209
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 04:00:00
#SBATCH --array=0-161
#SBATCH -o logs/s12_%A_%a.out
set -euo pipefail
cd /N/slate/dayusun/skmle_asynccv/cov1000
module load r/4.4.1
export R_LIBS=$PWD/lib:$HOME/r-packages R_LIBS_USER=$PWD/lib OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
read gen link n first last < <(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" tasks12.txt)
Rscript run12.R "$gen" "$link" "$n" "$first" "$last"
