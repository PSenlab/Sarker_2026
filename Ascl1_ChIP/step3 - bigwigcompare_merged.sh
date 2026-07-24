#!/bin/bash

#SBATCH --partition=norm
#SBATCH --cpus-per-task=12
#SBATCH --time=48:00:00
#SBATCH --mem=100g

module load deeptools/3.5.4
projpath="/path"

#Ascl1_M
subtract input from Ascl1_M
bigwigCompare -b1 $projpath/merged_bigwig/Ascl1_M.merged.RPKM.bw -b2 $projpath/merged_bigwig/Input_M.merged.RPKM.bw --operation subtract -o $projpath/merged_bigwig/Ascl1-input_M.bw -of bigwig

#Ascl1_F
subtract input from Ascl1_F
bigwigCompare -b1 $projpath/merged_bigwig/Ascl1_F.merged.RPKM.bw -b2 $projpath/merged_bigwig/Input_F.merged.RPKM.bw --operation subtract -o $projpath/merged_bigwig/Ascl1-input_F.bw -of bigwig
