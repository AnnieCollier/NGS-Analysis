###Specify path to bam files
p1="/path/to/bamfile/"
#name your sample(s)
b1="$p1/filename1.bam"
b2="$p1/filename2.bam"

ml python/3.6.1
#index your bam file
samtools index $b2
#count how many reads map to each feature
htseq-count -s no --additional-attr=gene_name --additional-attr=gene_type -f bam $b1 $b2  /path/to/gtf/gtfile.gtf  > output_count_table.txt 
