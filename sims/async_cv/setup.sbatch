#!/bin/bash
#SBATCH -J acv_setup
#SBATCH -A r00898
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 00:30:00
#SBATCH -o logs/setup_%j.out
set -euo pipefail
cd /N/slate/dayusun/skmle_asynccv
module load r/4.4.1
mkdir -p lib
export R_LIBS=$PWD/lib:$HOME/r-packages R_LIBS_USER=$PWD/lib MAKEFLAGS=-j4
R CMD INSTALL -l lib skmle_0.1.0.tar.gz
Rscript -e 'library(skmle); cat("skmle OK; locf+eeq:", all(c("locf","eeq") %in% eval(formals(kee_async_cv)$criterion)), "\n")'
