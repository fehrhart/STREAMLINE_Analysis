export PATH=$PATH:~/persistent/software/STAR-2.7.11b/source

./STAR --runThreadN 50 \ 
     --runMode genomeGenerate \
     --genomeDir ~/persistent/STAR_OutPut/ \
     --genomeFastaFiles ~/persistent/GRCh37.p13.genome.fa \
     --sjdbGTFfile ~/persistent/gencode.v45.annotation.gtf \
     --sjdbOverhang 150\
