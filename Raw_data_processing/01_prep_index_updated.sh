#!/bin/bash
#Index reference genome 

#SBATCH -J "rsem_index"
#SBATCH -o %j%x.log
#SBATCH -e %j%x.err
#SBATCH -c 8
#SBATCH -p cpu
#SBATCH --time=4:00:00
#SBATCH --mem=80G


#Load required modules
module load star/2.7.9a
module load gcc/10.2.0
module load rsem/1.3.3


#Arguments in order are: annotation file (--gtf flag), reference genome file, [output directory/file_prefix], use STAR to prepare the index
rsem-prepare-reference --star --gtf /nese/meclab/Shared/reference_genomes/ST_reference_genomes/rCheMyd1_202105/Annotation/GCF_015237465.2_rCheMyd1.pri.v2_genomic.gtf /nese/meclab/Shared/reference_genomes/ST_reference_genomes/rCheMyd1_202105/Annotation/GCF_015237465.2_rCheMyd1.pri.v2_genomic.fna ./01_genome_index/GCF_015237465.2_rCheMyd1.pri.v2_genomic


