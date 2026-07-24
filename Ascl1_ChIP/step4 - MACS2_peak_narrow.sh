#!/bin/bash

#SBATCH --partition=norm
#SBATCH --cpus-per-task=12
#SBATCH --time=48:00:00
#SBATCH --mem=100g
#SBATCH --gres=lscratch:20

module load macs/2.2.7.1

projpath=/path/alignment/bam
OUTPUT=/path/macs2_peakcalling_narrow_q0.001

macs2 callpeak -t $projpath/Ascl1_M_merged.bam -c $projpath/Input_M_merged.bam -n Ascl1_M_MACS --outdir ${OUTPUT} -f BAMPE -g 1.87e9 -B -q 0.001 --keep-dup all

macs2 callpeak -t $projpath/Ascl1_F_merged.bam -c $projpath/Input_F_merged.bam -n Ascl1_F_MACS --outdir ${OUTPUT} -f BAMPE -g 1.87e9 -B -q 0.001 --keep-dup all
