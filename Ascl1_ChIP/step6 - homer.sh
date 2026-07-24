#!/bin/bash

#SBATCH --partition=norm
#SBATCH --cpus-per-task=12
#SBATCH --time=48:00:00
#SBATCH --mem=40g
#SBATCH --gres=lscratch:20

#Motif analysis using 100bp as default size and using a random background
module load homer/5.1

path1="/path/Homer/Ascl1_F_MACS_q0.001"
path2="/path/Homer/Ascl1_M_MACS_q0.001"

# for Ascl1_F
findMotifsGenome.pl ${path1}/Ascl1_F_narrowpeaks.bed mm10 ${path1} -size 100 -useNewBg

# for Ascl1_M
findMotifsGenome.pl ${path2}/Ascl1_M_narrowpeaks.bed mm10 ${path2} -size 100 -useNewBg




