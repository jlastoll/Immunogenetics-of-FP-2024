#!/bin/bash
# sickle scythe combo
#run sickle and scythe on all samples
#for MHC data which is not large, can run on interactive queue in less than 5 min

#SBATCH -e sicklescythe.err
#SBATCH -o sicklescythe.log


date

#load modules and set variables

ADAPTOR="/nese/meclab/Jamie/chmy_RNAseq_MHC/scripts/MHC_adaptors.fasta"
CONTAM="0.1"

module load conda
export CONDA_ENVS_PATH=/nese/meclab/Jamie/.conda/envs/:/home/jstoll_umass_edu/.conda/envs/
conda activate test

#adaptor trim
for file in *.fastq
do
echo $file
#date
#gunzip $file
sample=$(basename $file .fastq ) #adjust so that $sample includes name+ R1 or R2 etc.
echo $sample
scythe -a $ADAPTOR -p $CONTAM -o ./"$sample"_S.fastq "$sample".fastq
date


done

wait
for R1 in *_R1_S.fastq
do

echo $R1
R2=${R1/_R1/_R2}
echo $R2
sample=$(basename $R1 _R1_S.fastq)
echo $sample
sickle pe -f $R1 -r $R2 -t sanger -o "$sample"_R1_SS.fastq -p "$sample"_R2_SS.fastq -s "$sample"_orphans.fastq.gz -q 20
date
#rm $R1
#rm $R2
#gzip "$sample"_R1_SS.fastq
#gzip "$sample"_R2_SS.fastq
done

conda deactivate
