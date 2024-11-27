#!/bin/bash
#BSUB -o extractmhc.log
#BSUB -e extractmhc.err
#BSUB -q short
#BSUB -W 4:00
#BSUB -R rusage[mem=4000]
#BSUB -R span[hosts=1]
#BSUB -n 16

module load samtools/1.9
regions=/project/uma_lisa_komoroske/Jamie/chmy_RNAseq_MHC/MHC_classI_regions.bed
 
for file in ./*.bam
do
sample=$(basename $file .bam)
samtools view -bh $file -L $regions  > {$sample}_mhc.bam
done

