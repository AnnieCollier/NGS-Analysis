#!/bin/bash                                                                                                                               

#SBATCH --job-name=homer --output=homer.out --error=homer.err --time=24:00:00 --qos=normal --nodes=1 --mem-per-cpu=4000 --ntasks-per-node=8

PATH=$PATH:/home/gautams/bioinformatics_tools/samtools-0.1.19/
PATH=$PATH:/share/PI/oro/software/ghostscript-9.18-linux-x86_64/
PATH=$PATH:/share/PI/oro/software/weblogo/
PATH=$PATH:/share/PI/oro/software/
PATH=$PATH:/share/PI/oro/software/homer/bin/

#Call peaks with Homer
findPeaks Tags_WT1bowtie2/ -o WT1botwie2peakfile -minDist 150 -region

#Find motifs in a bed file
findMotifsGenome.pl KO_only_10KB_overlap_0.01_together.bed hg38 KO_only_homer/ -size 200 -preparsedDir /oak/stanford/groups/oro/anncoll/temp_homer 

#find the enrichment of peaks in a bed file
annotatePeaks.pl WT_only_10kb.bed hg38 -size 40000 -hist 500 -d /oak/stanford/groups/oro/anncoll/NextSeqChIP/new_GATA3_AHDC1KO/Mapping/WT1GATA3 > GATA3_in_WT_only.txt

