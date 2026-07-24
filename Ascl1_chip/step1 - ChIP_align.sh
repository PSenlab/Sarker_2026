#!/bin/bash

#SBATCH --partition=norm
#SBATCH --cpus-per-task=24
#SBATCH --time=72:00:00
#SBATCH --mem=200g
#SBATCH --gres=lscratch:100

module load bowtie/2
module load fastqc/0.12.1
module load multiqc/1.34
module load samtools/1.23
module load picard/3.1.0
module load sambamba/1.0.1
module load gcc/11.3.0
module load trimgalore/2.2.0
module load bedtools/2.31.1
module load deeptools/3.5.6

INPUT=/path/fastq
INPUT1=/path/fastq/*R1.fastq.gz
OUTPUT1=/path/adapter_trimmed
INPUT2=/path/adapter_trimmed/*gz
OUTPUT2=/path/fastqc
INPUT3=/path/adapter_trimmed/*_*_R1*gz
mkdir -p /path/alignment/sam
OUTPUT3=/path/alignment/sam
mkdir -p /path/alignment/sam/bowtie_summary
mkdir -p /path/alignment/bam
OUTPUT4=/path/alignment/sam/bowtie_summary
OUTPUT5=/path/alignment/bam
mkdir -p /path/alignment/picard_summary
OUTPUT6=/path/alignment/picard_summary
OUTPUT7=/path/bigwig
BLACKR=/path/mm10-blacklist.bed.gz
#Directing Bowtie2 to the genome build (here mm10). Bowtie2 indices are available as part of the igenomes package
export BOWTIE2_INDEXES=/fdb/igenomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/

#Trim adaptors with trim_galore
for i in ${INPUT1}
do
base_name=$(basename $i)
sample_name=${base_name%%_R1*}
trim_galore --paired ${i} ${INPUT}/${sample_name}_R2.fastq.gz -o $OUTPUT1
done

#Performing fastqc on the adapter trimmed files
for i in ${INPUT2}
do
fastqc -f fastq -o ${OUTPUT2} ${i} 
done

#multiqc on fastqc files
multiqc -i Ascl1_ChIP_trial2_multiqc -o /path/multiqc /path/fastqc

#Performing alignment (bowtie2), generating filtered bam files (samtools)
#PCR duplicates removal with picard
#Removing regions that intersect with ENCODE blacklisted regions
#index bam
#generate bigwig
for i in ${INPUT3}
do
base_name=$(basename $i)
sample_name=${base_name%%_R*}
bowtie2 -p 24 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 --threads=$SLURM_CPUS_PER_TASK -I 10 -X 700 -x genome -1 ${i} \
-2 ${OUTPUT1}/${sample_name}_R2_val_2.fq.gz  -S ${OUTPUT3}/${sample_name}.sam &> ${OUTPUT4}/${sample_name}_bowtie2.txt
    
#Generating filtered bam files (samtools)
samtools view -@ 24 -h -F 4 -q 10 -bS ${OUTPUT3}/${sample_name}.sam > ${OUTPUT5}/${sample_name}_filtered.bam
    
#PCR duplicates removal with picard
java -jar $PICARDJARPATH/picard.jar SortSam -I ${OUTPUT5}/${sample_name}_filtered.bam -O ${OUTPUT5}/${sample_name}_sorted.bam -SORT_ORDER coordinate
sambamba view -h -t 2 -f bam -F "[XS] == null and not unmapped and not duplicate" ${OUTPUT5}/${sample_name}_sorted.bam > ${OUTPUT5}/${sample_name}_sorted_unique.bam
java -jar $PICARDJARPATH/picard.jar MarkDuplicates I=${OUTPUT5}/${sample_name}_sorted_unique.bam O=${OUTPUT5}/${sample_name}_noDUP_unique.bam \
REMOVE_DUPLICATES=true METRICS_FILE=${OUTPUT6}/${sample_name}_noDUP.txt
   
#Removing regions that intersect with ENCODE blacklisted regions
bedtools intersect -a ${OUTPUT5}/${sample_name}_noDUP_unique.bam -b $BLACKR -v > ${OUTPUT5}/${sample_name}.noB_unique.bam
    
#Indexing new Bam files
samtools index ${OUTPUT5}/${sample_name}.noB_unique.bam
    
#Generating bigwig files
bamCoverage -b ${OUTPUT5}/${sample_name}.noB_unique.bam -o ${OUTPUT7}/${sample_name}.noB_unique.bw -of bigwig --normalizeUsing RPKM -p $SLURM_CPUS_PER_TASK
    
#Remove intermediate bam files
rm ${OUTPUT5}/${sample_name}_filtered.bam
rm ${OUTPUT5}/${sample_name}_sorted.bam 
rm ${OUTPUT5}/${sample_name}_sorted_unique.bam
rm ${OUTPUT5}/${sample_name}_noDUP_unique.bam    
done
