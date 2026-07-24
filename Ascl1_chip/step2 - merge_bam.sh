#!/bin/bash

#SBATCH --partition=norm
#SBATCH --cpus-per-task=16
#SBATCH --mem=24g
#SBATCH --gres=lscratch:20
#SBATCH --time=8:00:00

module load samtools/1.23
module load deeptools/3.5.6

TMP_DIR=/lscratch/$SLURM_JOBID
INPUT=/path/alignment/bam
INPUT1=/path/alignment/bam/*merged.bam
mkdir -p /path/merged_bigwig
OUTPUT=/path/merged_bigwig

samtools merge ${INPUT}/Ascl1_F_merged.bam ${INPUT}/Ascl1_F*.noB_unique.bam
samtools merge ${INPUT}/Ascl1_M_merged.bam ${INPUT}/Ascl1_M*.noB_unique.bam
samtools merge ${INPUT}/Input_F_merged.bam ${INPUT}/Input_F*.noB_unique.bam
samtools merge ${INPUT}/Input_M_merged.bam ${INPUT}/Input_M*.noB_unique.bam

for i in ${INPUT1}
    do
    base_name=$(basename $i)
    sample_name=${base_name%%_merged*}
    samtools index ${i}
    #Generating bigwig files
    bamCoverage -b ${i} -o ${OUTPUT}/${sample_name}.merged.RPKM.bw -of bigwig --normalizeUsing RPKM -p $SLURM_CPUS_PER_TASK
    done

  
    































