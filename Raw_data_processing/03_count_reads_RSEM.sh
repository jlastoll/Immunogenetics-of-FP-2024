#!/bin/bash
#Calculate transcript abundance

#BSUB -q long
#BSUB -W 40:00
#BSUB -R rusage[mem=16000]
#BSUB -n 12
#BSUB -R span[hosts=1]
#BSUB -e logs/03_count_reads2.err
#BSUB -oo logs/03_count_reads2.log

#Load required modules
module load RSEM/1.3.0

for file  in 02_mapped_reads/Sample*Aligned.toTranscriptome.out.bam 

do

sample=${file%%_Aligned.toTranscriptome.out.bam} #Store path and basename of file minus read info

out_name=${sample##*/} #Remove path, but keep basename top be used as a prefix to output
echo $sample
echo $out_name

#Arguments are (after paired-end): input file is bam file, don't output bam file, use 12 threads, specify paired end sequences, stranded data, location/prefix of rsem index files, location/prefix of rsem calulate expression output files (from this run)

rsem-calculate-expression --bam --no-bam-output -p 12 --paired-end --strandedness reverse $file 01_genome_index/GCF_015237465.2_rCheMyd1.pri.v2_genomic 03_transcript_quant/$out_name 

done

#Once the above loop is finished, combine results from all samples into one data matrix for R. This matrix is of RAW COUNTS at gene level, not normalized
rsem-generate-data-matrix 03_transcript_quant/*.genes.results > all.counts.matrix
