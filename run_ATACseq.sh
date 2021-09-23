#!/bin/bash
PATH=$PATH:"/home/groups/oro/software/scripts/ATACseq/"

for FILE in config_file_WT1; do
    echo ${FILE}
    sbatch  -n 1 --time=48:00:00 --job-name=ATAC_$FILE --output=ATAC_$FILE.out --error=ATAC_$FILE.err --nodes=6 --ntasks=6 -c 2 --mem-per-cpu=4000 --ntasks-per-node=1 --wrap="bash ATACseq_pipeline.sh -C $FILE"
    sleep 1
done
