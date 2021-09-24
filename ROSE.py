#!/bin/bash

#SBATCH --job-name=Homer --output=Homer.out --error=Homer.err --time=1:00:00 --qos=normal --nodes=1 --mem-per-cpu=4000 --ntasks-per-node=8

python ROSE_main.py -g hg19 -i H3K27Ac.gff -r H3K27Ac.bam -o /scratch/users/anncoll