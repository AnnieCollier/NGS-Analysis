# NGS-Analysis
This repository contains various code for mapping and analyzing NGS data, including differential expression and motif finding.

# RNA-seq 
Use the RNASEQ.sh file to map and process raw fastq files from RNA-seq experiments. This method uses TOPHAT but other mapping programs are probably recommended. 
Use the RNA_DESEQ2.R script to run differential expression analysis on the raw count file. 

# ATAC or ChIP-seq Mapping:
1. First, unzip and rename your fastq files
2. Next, edit the config_file text file with the correct fastq file prefix, output prefix, species, p-value, and peak type for your sample.
3. Save the config file with a name unique to that sample (i.e. "config_file_WT1").
4. Make sure the config file, run_seq, and seq_pipeline files are all in the same directory. Do not change anything in the other two scripts.
5. Finally, run the "run_seq.sh" as a batch script (sbatch run_seq.sh). This will run the seq_pipeline.sh script using the information provided in each config file. 
6. The pipeline will eventually output Mapping, peakCalling, Motif and tracks output directories, which you can use for downstream analysis.

# ChIP-seq 
To make a final bed file with reproducible peaks, perform Irreproducible Discovery Rate ,using the IDR.sh script
For differential analysis, use the DiffBind.R script and bam files

# ATAC-seq 
For differential analysis, use the union_macs script to make a union table and then the ATAC_DESEQ2.R script to perform differential analysis. 

# Illumina Methylation Arrays
Follow the Methylation_Array.R script. 

# ChromHMM
ChromHMMM maps the enrichment of different chromatin states relative to a set of coordinates (bed file). First you need to establish the chromatin states using histone mark ChIP-seq or ATAC-seq data. 
1. Obtain bamtobed ("srt_rmdup.bed") files from ChIP or ATAC mapping
2. fill out a cellmarktable text file with the files you want to use to make chromatin states
3. Binarize the bed files using the chromHMM_binarizebed.sh script
4. Specify the names and paths of your region of interest bed file in the chromHMM python script and run. 
5. Program will output images and text files for your region of interest.

![image](https://user-images.githubusercontent.com/90862478/134739895-391de034-1452-4c5b-ba5b-ac29dae2696d.png)



# Other useful tooks and Plotting
Useful HOMER scripts can be found in the homer.sh file. These scripts are useful for:
- motif finding
- peak enrichment in a set of coordinates

ROSE is used to rank enhancers by H3K27Ac signal. Change the paths and names in the ROSE.py script and run as a batch command. 

Deeptools can be found in the deeptools.sh script. This can be used to:
- make heatmaps of ChIP/ATAC data
- make profile plots of ChIP/ATAC data

![image](https://user-images.githubusercontent.com/90862478/134739913-1e962143-c604-498a-a282-f64783eee403.png)


Use MinDist.py script to find distances between ChIP-seq peaks (this is Sam's)

![image](https://user-images.githubusercontent.com/90862478/134739935-ff22b847-886c-41f8-b9d0-69ed03a89c82.png)


Use Virtual4C script to make virtual 4C plots to visualize raw cohesin hiChIP data. For this, you will need normalization numbers and coordinates of interest.

![GREM2H9d7_5000 chr4_240770000 matrix](https://user-images.githubusercontent.com/90862478/134740008-7a453b1c-8ba9-4e01-95e6-207e38c1122f.png)

