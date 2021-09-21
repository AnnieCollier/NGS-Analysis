#!/bin/bash                                                                                                                               

#SBATCH --job-name=processing --output=processing.out --error=processing.err --time=4:00:00 --qos=normal --nodes=1 --mem-per-cpu=4000 --ntasks-per-node=8

PATH=$PATH:/home/users/gautams/bioinformatics_tools/samtools-0.1.19/
PATH=$PATH:/share/PI/oro/software/ghostscript-9.18-linux-x86_64/
PATH=$PATH:/share/PI/oro/software/weblogo/
PATH=$PATH:/share/PI/oro/software/
PATH=$PATH:/share/PI/oro/software/homer/bin/
PATH=$PATH:/share/PI/oro/software/bowtie-1.1.2/

samtools view -bS -F 4 TFAP2CBASU1.sam > TFAP2CBASU1_norm.bam
samtools view -bS -F 4 TFAP2CBASU2.sam > TFAP2CBASU2_norm.bam
samtools view -bS -F 4 INPUT.sam > INPUT_norm.bam
samtools view -bS -F 4 GRHL2D7A.sam > GRHL2D7A_norm.bam
samtools view -bS -F 4 GRHL2D7B.sam > GRHL2D7B_norm.bam

samtools sort GATA1_norm.bam GATA1_norm_sorted
samtools sort TFAP2CBASU2_norm.bam TFAP2CBASU2_norm_sorted
samtools sort INPUT_norm.bam INPUT_norm_sorted
samtools sort GRHL2D7A_norm.bam GRHL2D7A_norm_sorted
samtools sort GRHL2D7B_norm.bam GRHL2D7B_norm_sorted

samtools rmdup -s TFAP2CBASU1_norm_sorted.bam TFAP2CBASU1_norm_sorted_rmdup.bam
samtools rmdup -s TFAP2CBASU2_norm_sorted.bam TFAP2CBASU2_norm_sorted_rmdup.bam
samtools rmdup -s INPUT_norm_sorted.bam INPUT_norm_sorted_rmdup.bam
samtools rmdup -s GRHL2D7A_norm_sorted.bam GRHL2D7A_norm_sorted_rmdup.bam
samtools rmdup -s GRHL2D7B_norm_sorted.bam GRHL2D7B_norm_sorted_rmdup.bam

samtools flagstat TFAP2CBASU1_norm_sorted_rmdup.bam > TFAP2CBASU1_norm_sorted_rmdup_flagstat.txt
samtools flagstat TFAP2CBASU2_norm_sorted_rmdup.bam > TFAP2CBASU2_norm_sorted_rmdup_flagstat.txt
samtools flagstat INPUT_norm_sorted_rmdup.bam > INPUT_norm_sorted_rmdup_flagstat.txt
samtools flagstat GRHL2D7A_norm_sorted_rmdup.bam > GRHL2D7A_norm_sorted_rmdup_flagstat.txt
samtools flagstat GRHL2D7B_norm_sorted_rmdup.bam > GRHL2D7B_norm_sorted_rmdup_flagstat.txt

rm -rf TFAP2CBASU1_norm.bam
rm -rf TFAP2CBASU2_norm.bam
rm -rf INPUT_norm.bam
rm -rf GRHL2D7A_norm.bam
rm -rf GRHL2D7B_norm.bam

makeTagDirectory Tags_TFAP2CBASU1/ TFAP2CBASU1_norm_sorted_rmdup.bam
makeTagDirectory Tags_TFAP2CBASU2/ TFAP2CBASU2_norm_sorted_rmdup.bam
makeTagDirectory Tags_INPUT/ INPUT_norm_sorted_rmdup.bam
makeTagDirectory Tags_GRHL2D7A/ GRHL2D7A_norm_sorted_rmdup.bam
makeTagDirectory Tags_GRHL2D7B/ GRHL2D7B_norm_sorted_rmdup.bam

makeUCSCfile Tags_TFAP2CBASU1/ -o auto
makeUCSCfile Tags_TFAP2CBASU2/ -o auto
makeUCSCfile Tags_INPUT/ -o auto
makeUCSCfile Tags_GRHL2D7A/ -o auto
makeUCSCfile Tags_GRHL2D7B/ -o auto