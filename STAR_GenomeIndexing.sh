export PATH=$PATH:root/persistent/software/STAR-2.7.11b/source
# loading STAR software
module load STAR

# genome indexing

STAR --runThreadN 50 \ 
     --runMode genomeGenerate \
     --genomeDir ./STAR_OutPut/ \
     --genomeFastaFiles ./S/ \
     --sjdbGTFfile ./gencode.v45.annotation.gtf \
     --sjdbOverhang 150\
