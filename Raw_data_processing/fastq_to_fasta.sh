#!/bin/bash
#BSUB -o ./logs/fastq_to_fasta.log
#BSUB -e ./logs/fastq_to_fasta.err
#BSUB -q long
#BSUB -W 24:00
#BSUB -R rusage[mem=4000]
#BSUB -R span[hosts=1]
#BSUB -n 16

module load fastx_toolkit/0.0.14



for file in /project/uma_lisa_komoroske/Jamie/chmy_RNAseq_MHC/trimmed_fastq_files/HI_samps/Sample*_R1_SS.fastq

do
sample=$(basename $file _R1_SS.fastq)
echo $sample
#gunzip $file
 
fastq_to_fasta -z -i /project/uma_lisa_komoroske/Jamie/chmy_RNAseq_MHC/trimmed_fastq_files/HI_samps/${sample}_R1_SS.fastq -o /project/uma_lisa_komoroske/Jamie/chmy_RNAseq_MHC/trimmed_fastq_files/HI_samps/${sample}_R1_SS.fasta.gz

done

