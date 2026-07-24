#!/bin/bash

#SBATCH --partition=norm
#SBATCH --cpus-per-task=8
#SBATCH --mem=24g
#SBATCH --gres=lscratch:10
#SBATCH --time=23:00:00

module load bwtool/1.0
TMP_DIR=/lscratch/$SLURM_JOBID

# set variable

INPUT1=/path/AUC/Ascl1_F_narrowpeaks.bed
INPUT2=/path/merged_bigwig/Ascl1-input_F.bw
INPUT3=/path/AUC/Ascl1_M_narrowpeaks.bed
INPUT4=/path/merged_bigwig/Ascl1-input_M.bw
OUTPUT=/path/AUC

bwtool summary ${INPUT1} ${INPUT2} ${OUTPUT}/Ascl1_F_narrowpeaks_q0.001_AUC.txt -header -with-sum

bwtool summary ${INPUT3} ${INPUT4} ${OUTPUT}/Ascl1_M_narrowpeaks_q0.001_AUC.txt -header -with-sum

