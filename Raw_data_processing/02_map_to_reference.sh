#!/bin/bash
#Map reads to indexed reference genome

#BSUB -q long
#BSUB -W 40:00
#BSUB -R rusage[mem=32000]
#BSUB -n 8
#BSUB -R span[hosts=1]
#BSUB -e logs/02_map_to_reference.err
#BSUB -oo logs/02_map_to_reference.log

#Load required modules
module load star/2.7.0e
module load gcc/8.1.0


####Step 0: check if files created in the script below already exist. If they do, delete previous versions.
if test -f ./metrics/line_counts.txt; then
rm ./metrics/line_counts.txt
fi

if test -f ./metrics/unmapped_reads.txt; then
rm ./metrics/unmapped_reads.txt
fi

####Step 3: Map each sample to indexed reference genome using STAR
touch ./metrics/unmapped_reads.txt

for f in /project/uma_lisa_komoroske/Jamie/chmy_RNAseq_MHC/trimmed_fastq_files/HI_samps/*R1_SS.fastq.gz
do

sample2=${f%%R1_SS.fastq.gz} #Store path and basename of file minus read info
out_name2=${sample2##*/} #Remove path, but keep basename top be used as a prefix to output

STAR --runThreadN 8 --genomeDir 01_genome_index --readFilesIn ${sample2}R1_SS.fastq.gz ${sample2}R2_SS.fastq.gz --readFilesCommand zcat --quantMode TranscriptomeSAM --outFileNamePrefix ./02_mapped_reads/$out_name2

wait

echo "sample $out_name2" >> ./metrics/unmapped_reads.txt
grep "unmapped" ./02_mapped_reads/${out_name2}Log.final.out >> ./metrics/unmapped_reads.txt

done
