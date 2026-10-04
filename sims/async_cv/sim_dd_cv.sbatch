#!/bin/bash
#SBATCH -J acv_ddcv
#SBATCH -A r00898
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 01:00:00
#SBATCH --array=1-840
#SBATCH -o logs/ddcv_%A_%a.out
set -euo pipefail
cd /N/slate/dayusun/skmle_asynccv
module load r/4.4.1
export R_LIBS=$PWD/lib:$HOME/r-packages R_LIBS_USER=$PWD/lib OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# One line of tasks_dd.txt per task: n link first_rep n_reps.
read -r n link first nrep < <(sed -n "${SLURM_ARRAY_TASK_ID}p" tasks_dd.txt)
Rscript sim_dd_cv.R "$n" "$link" "$first" "$nrep" "results_ddcv/${link}_n${n}_r${first}.rds"
