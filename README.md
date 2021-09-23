# NGS-Analysis
Code for mapping and analyzing NGS data

Use the RNASEQ.sh file to map and process raw fastq files from RNA-seq experiments. This method uses TOPHAT but other mapping programs are probably recommended.

For ChIP-seq use the following workflow:
1. First, unzip and rename your fastq files
2. Next, edit the config_file text file with the correct fastq file prefix, output prefix, species, p-value, and peak type for your sample.
3. Save the config file with a name unique to that sample (i.e. "config_file_WT1").
4. Make sure the config file, run_chipseq, and ChIPseq_pipeline files are all in the same directory. Do not change anything in the other two scripts.
5. Finally, run the "run_chipseq.sh" as a batch script (sbatch run_chipseq.sh). This will run the ChIPseq_pipeline.sh script using the information provided in each config file. 
6. The pipeline will eventually output Mapping, peakCalling, and tracks output directories, which you can use for downstream analysis.
