#!/bin/bash
PATH=$PATH:"/home/groups/oro/software/scripts/ChIPseq/"

for FILE in config_file*; do
    echo ${FILE}
    sbatch -n 1 --time=24:00:00 --job-name=CTCFChIP_$FILE --output=CTCFChIP_$FILE.out --error=CTCFChIP_$FILE.err  --nodes=6 --ntasks=6 -c 2 --mem-per-cpu=4000 --ntasks-per-node=1 --wrap="bash ChIPseq_pipeline.sh -C $FILE"
    sleep 1
done


#for FILE in Input_config*; do
#    echo ${FILE}
#    sbatch -n 1 --time=10:00:00 --job-name=ChIP_$FILE --output=ChIP_$FILE.out --error=ChIP_$FILE.err --mem-per-cpu=4000 --ntasks-per-node=8 --wrap="bash ChIPseq_pipeline_input.sh -C $FILE"
#    sleep 1
#done

