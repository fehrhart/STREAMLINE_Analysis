#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=5:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address

export PATH=$PATH:root/persistent/software/STAR-2.7.11b/source
# loading STAR software
module load STAR

# genome indexing

STAR --runThreadN 50 \ 
     --runMode genomeGenerate \
     --genomeDir ./STAR_OutPut/ \
     --genomeFastaFiles ./S \
     --sjdbGTFfile ./gencode.v45.annotation.gtf \
     --sjdbOverhang 150\
