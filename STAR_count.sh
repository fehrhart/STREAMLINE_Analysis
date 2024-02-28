#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=24:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address


# Loading STAR software in the server
module load STAR

for i in ./*R1_001_fastp.fastq.gz
do

    R1=$i
    R2=${i/R1/R2}
    sample_name=$(echo $i| cut -d'/' -f 7)
    sample_name=$(echo ${sample_name%R1_001_fastp.fastq.gz})
    out_dir=./STAR_OutPut/${sample_name}
    temp_dir=./STAR_temp/${sample_name}

    STAR --runThreadN 16 \
         --quantMode GeneCounts \
         --genomeDir ./genome_index3_output/ \
         --outSAMtype None \
         --outFileNamePrefix ${out_dir}. \
         --readFilesIn $R1 $R2 \
         --readFilesCommand zcat \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMultimapNmax 20 \
         --genomeLoad LoadAndKeep \
         --limitBAMsortRAM 50000000000 \
         --outFilterScoreMinOverLread 0.5 \
         --outFilterMatchNminOverLread 0.5 \
         --outTmpDir $temp_dir \

done

