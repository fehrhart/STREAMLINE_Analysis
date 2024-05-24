for i in ~/persistent/S/*R1_001.fastq.gz
do

    R1=$i
    R2=${i/R1/R2}
    sample_name=$(echo $i| cut -d'/' -f 7)
    sample_name=$(echo ${sample_name%R1_001.fastq.gz})
    out_dir=~/persistent/Results/${sample_name}
    temp_dir=~/persistent/STAR_temp/${sample_name}

    STAR --runThreadN 120 \
         --quantMode GeneCounts \
         --genomeDir ~/persistent/STAR_OutPut/ \
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

