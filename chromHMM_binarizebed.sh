#!/bin/bash

#SBATCH --job-name=chromHMM --output=chromHMM.out --error=chromHMM.err --time=12:00:00 --qos=normal --nodes=1 --mem-per-cpu=4G --ntasks-per-node=8

#bed files are actually bamtobed files ("srt_rmdup.bed"). For replicates, list both in the cellmark cellmarktable text file.

#In terminal. Make the binarized bed files. Make sure your bamtobed and cellmarktable files are in the directory
java -jar /home/groups/oro/software/ChromHMM/ChromHMM.jar BinarizeBed /home/groups/oro/software/ChromHMM/CHROMSIZES/hg38.txt 5_lines_bed_files cellmarktable.txt states
#In terminal. Make the states
java -jar /home/groups/oro/software/ChromHMM/ChromHMM.jar LearnModel -printposterior states states/10 10 hg38 #change number if you want more or less states

#can now reference these states when you run chromHMM
