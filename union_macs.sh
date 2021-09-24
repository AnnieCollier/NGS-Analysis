# Manually enter this into Sherlock

file0=0
file1=1
file2=2
file3=3
file4=4
file5=5
file6=6
file7=7
file8=8
file9=9
file10=10
file11=11
file12=12
file13=13
file14=14
file15=15

for n in $file0 $file1 $file2 $file3 $file4 $file5 $file6 $file7 $file8 $file9 $file10 $file11 $file12 $file13 $file14 $file15; do
	echo '#! /bin/bash' > multicov_${n}.sh; 
	echo "#SBATCH --job-name=multicovATAC --output=multicovATAC.out --error=multicovATAC.err --time=24:00:00 --qos=normal --nodes=1 --mem-per-cpu=4000 --ntasks-per-node=8" >> multicov_${n}.sh;
	echo a='/oak/stanford/groups/oro/sandra/bams/ATACseq/ATACseq_TN5_H9_p63WT_E8_day0_merge_rep1_norm_sorted_rmdup_nochrM_rmdup.bam' >> multicov_${n}.sh;
	echo b='/oak/stanford/groups/oro/sandra/bams/ATACseq/ATACseq_TN5_H9_p63WT_E8_day0_merge_rep2_norm_sorted_rmdup_nochrM_rmdup.bam' >> multicov_${n}.sh;
	echo c='/oak/stanford/groups/oro/sandra/bams/ATACseq/ATACseq_TN5_H9_p63WT_RA_BMP4_day7_rep1_norm_sorted_rmdup_nochrM.bam' >> multicov_${n}.sh;
	echo d='/oak/stanford/groups/oro/sandra/bams/ATACseq/ATACseq_TN5_H9_p63WT_RA_BMP4_day7_rep2_norm_sorted_rmdup_nochrM.bam' >> multicov_${n}.sh;
	echo e='/oak/stanford/groups/oro/jillian/MiSeq/102616_2AKO_d7_ATAC/ATACseq_TN5_2AKO_d7_rep1_norm_sorted_rmdup_nochrM.bam' >> multicov_${n}.sh;
	echo f='/oak/stanford/groups/oro/jillian/MiSeq/102616_2AKO_d7_ATAC/ATACseq_TN5_2AKO_d7_rep2_norm_sorted_rmdup_nochrM.bam' >> multicov_${n}.sh;
	echo s='/scratch/users/jillianp/LL_data/ATACseq/ATACseq_TN5_AP2Cd7off_rep1_norm_sorted_rmdup_nochrM.bam' >> multicov_${n}.sh;
	echo t='/scratch/users/jillianp/LL_data/ATACseq/ATACseq_TN5_AP2Cd7off_rep2_norm_sorted_rmdup_nochrM.bam' >> multicov_${n}.sh;
	echo q='/scratch/users/jillianp/LL_data/ATACseq/ATACseq_TN5_AP2Cd7on_rep1_norm_sorted_rmdup_nochrM.bam' >> multicov_${n}.sh;
	echo r='/scratch/users/jillianp/LL_data/ATACseq/ATACseq_TN5_AP2Cd7on_rep2_norm_sorted_rmdup_nochrM.bam' >> multicov_${n}.sh;
	echo PATH='$PATH:/home/users/gautams/bioinformatics_tools/samtools-1.2/' >> multicov_${n}.sh; 
	echo PATH='$PATH:/home/users/gautams/bioinformatics_tools/bedtools2/bin/' >> multicov_${n}.sh; 
	echo PATH='$PATH:/home/users/gautams/bioinformatics_tools/picard-tools-1.129/' >> multicov_${n}.sh; 
	echo PATH='$PATH:/home/users/gautams/bioinformatics_tools/ucsc_tools/executables/' >> multicov_${n}.sh; 
	echo PATH='$PATH:/home/users/gautams/bioinformatics_tools/bowtie2-2.2.5/' >> multicov_${n}.sh; 
	echo bedtools multicov -bams '$a $b $c $d $e $f $s $t $q $r' -bed union_${n} '>multicov_union_'${n}.txt >> multicov_${n}.sh;
done

sbatch multicov_0.sh
sbatch multicov_1.sh
sbatch multicov_2.sh
sbatch multicov_3.sh
sbatch multicov_4.sh
sbatch multicov_5.sh
sbatch multicov_6.sh
sbatch multicov_7.sh
sbatch multicov_8.sh
sbatch multicov_9.sh
sbatch multicov_10.sh
sbatch multicov_11.sh
sbatch multicov_12.sh
sbatch multicov_13.sh
sbatch multicov_14.sh
sbatch multicov_15.sh
