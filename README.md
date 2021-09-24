# NGS-Analysis
Code for mapping and analyzing NGS data

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
For differential analysis, use the DiffBind.R script and bam files

# ATAC-seq 
For differential analysis, use the union_macs script to make a union table and then the ATAC_DESEQ2.R script to perform differential analysis. 

# Illumina Methylation Arrays
Follow the Methylation_Array.R script. 

Useful HOMER scripts can be found in the homer.sh file. These scripts are useful for:
- motif finding
- 
