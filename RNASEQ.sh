#!/bin/bash                                                                                                                               

#SBATCH --job-name=RNASEQ --output=RNASEQ.out --error=RNASEQ.err --time=4:00:00 --qos=normal --nodes=1 --mem-per-cpu=4000 --ntasks-per-node=8

PATH=$PATH:/home/users/gautams/bioinformatics_tools/samtools-0.1.19/
PATH=$PATH:/share/PI/oro/software/ghostscript-9.18-linux-x86_64/
PATH=$PATH:/share/PI/oro/software/weblogo/
PATH=$PATH:/share/PI/oro/software/
PATH=$PATH:/share/PI/oro/software/homer/bin/
PATH=$PATH:/share/PI/oro/software/bowtie-1.1.2/
ml biology
ml tophat/2.1.1
ml bowtie2/2.3.4.1

#unzip and rename your fastq.gz files
zcat ABCG_124_R1.fastq.gz > WT1_R1.fastq

#map your fastq files to the hg38 genome
tophat -p 10 --library-type fr-firststrand -r 100 --mate-std-dev 100 -o WT1hg38 /oak/stanford/groups/oro/reference/latest_refGenome/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome WT1_R1.fastq WT1_R.fastq

samtools view -bS -F 4 WT1.sam > WT1_norm.bam

samtools sort WT1_norm.bam WT1_norm_sorted

samtools rmdup -s WT1_norm_sorted.bam WT1_norm_sorted_rmdup.bam

samtools flagstat WT1_norm_sorted_rmdup.bam > WT1_norm_sorted_rmdup_flagstat.txt

rm -rf WT1_norm.bam

makeTagDirectory Tags_WT1/ WT1_norm_sorted_rmdup.bam

makeUCSCfile Tags_WT1/ -o auto

#make raw and rpkm union tables with all samples you want to compare with DESEQ1
#specify the path to each file
d1=WT1/WT1
d2=WT2/WT2
d3=KO1/KO1
d4=KO2/KO2

#use the RAW file for DESEQ2 differential expression analysis
analyzeRepeats.pl rna hg38 -strand both -count exons -condenseGenes -d $d1 $d2 $d3 $d4 -rpkm > WT_AHDC1KO_iKC_rpkm.txt
#use the RPKM file to get look at expression levels of different genes
analyzeRepeats.pl rna hg38 -strand both -count exons -condenseGenes -d $d1 $d2 $d3 $d4 -noadj > WT_AHDC1KO_iKC__raw.txt

