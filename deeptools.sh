#!/bin/bash                                                                                                                               

#SBATCH --job-name=deeptools --output=deeptools.out --error=deeptools.err --time=24:00:00 --qos=normal --nodes=1 --mem-per-cpu=4000 --ntasks-per-node=8

ml biology
ml py-deeptools

#Convert individual BAMs to bigwig using deeptools
bamCoverage -b WT2hg38_srt_nochrM_rmdup.bam -o deeptools_WT1.bw -bs=1 -p=max

#Convert two BAMs to bigwig together so that they are normalized
bamCompare -b1 WT2hg38_srt_nochrM_rmdup.bam -b2 AHDC1hg38_srt_nochrM_rmdup.bam -o WTAHDC1compare.bw

#Make a matrix file of read counts centered in the middle of chosen peaks with 1000 bp on each end
#For 1 treatment group
computeMatrix reference-point -S WTAHDC1compare.bw -R GATA3_hg38_0.05IDR.txt.npk.bed --referencePoint center -b 1000 -a 1000 -bs=1 -p=max -out deeptools_WT_AHDC1_compare.matrix.gz
#For 2 treatment groups
computeMatrix reference-point -S deeptools_WT1.bw deeptools_GATA1.bw -R GATA3_hg38_0.05IDR.txt.npk.bed --referencePoint center -b 1000 -a 1000 -bs=1 -p=max -out deeptools_WT_GATA1_ATAC.matrix.gz

#Plotting
#can change the file format, coordinates, labels etc.
#Generate a heatmap PDF
plotHeatmap -m deeptoolsAHDC1chi.matrix.gz -out deeptoolsAHDC1GATA.plot.svg --colorMap=Blues --plotFileFormat svg
#Generate a profile
plotProfile -m deeptoolschips.matrix.gz -out COMPAREchips.png
