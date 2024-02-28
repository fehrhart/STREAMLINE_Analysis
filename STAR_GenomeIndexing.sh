#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=5:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address


# loading STAR software
module load STAR

# genome indexing

STAR --runThreadN 16 \ 
     --runMode genomeGenerate \
     --genomeDir ./genome_index3_output/ \
     --genomeFastaFiles ./GRCh37.p13.genome.fa \
     --sjdbGTFfile ./gencode.v19.annotation.gff3 \
     --sjdbOverhang 150\
