###THIS SCRIPT WAS CREATED BY SADHANA GADDAM, STANFORD UNIVERSITY
#!/bin/bash                                                                                                                                                        

ml biology
ml python/2.7.13
ml py-cutadapt/1.18_py27
ml samtools/1.8
ml bowtie2/2.3.4.1
ml trim_galore/0.5.0
ml py-macs2
ml bedtools
ml R/3.5.1
PATH=$PATH:/share/PI/oro/software/homer/bin/
PATH=$PATH:/share/PI/oro/software/ucsc_tools/executables/


usage(){
cat << EOF

Options:
        -C  ConfigFile          Name of the configuration file.
EOF
}

while getopts "C:" opt;
do
        case "$opt" in
                C) ConfigFile=$OPTARG;;
                \?) usage
                        echo "error: unrecognized option -$OPTARG";
                        exit 1
                        ;;
        esac
done

###parameters$###

fw=""
rw=""
samplename=""
controlname=""
refGen="hs" ##mm/hs
type="TF"
pvalue=""
out_map="Mapping"
out_tracks="tracks"
out_peaks="peakCalling"
out_motif="Motif"

# read the configuration file and store various parameters

echo -e "\n ================ Parsing input configuration file ================="

# separator used in the config file
IFS="="
while read -r name value
do
	param=$name
	paramval=${value//\"/}
	if [[ -n $param ]]; then
		if [[ $param != \#* ]]; then
			# if there are multiple parameter values (separated by # - old values are kept)
			# then the following operation selects the current one
			paramval=$(echo "$paramval" | awk -F['#\t'] '{print $1}' | tr -d '[:space:]');
			echo -e "Content of $param is $paramval"
			if [ $param == "fowardFastq" ]; then
			    fw=$paramval
			fi
			if [ $param == "reverseFastq" ]; then
				if [[ ! -z $paramval ]]; then
					rw=$paramval
				fi
			fi
			if [ $param == "samplename" ]; then
				samplename=$paramval
			fi
			if [ $param == "controlname" ]; then
				if [[ ! -z $paramval ]]; then
					controlname=$paramval
				fi
			fi
			if [ $param == "referenceGenome" ]; then
				refGen=$paramval
			fi
			if [ $param == "ChIPseq_type" ]; then
				type=$paramval
			fi
			if [ $param == "pvalue" ]; then
				pvalue=$paramval
			fi
		fi
	fi
done < $ConfigFile

## assign lib information single or paired
if [ ! -z $rw ] 
then 
 lib="paired"
else
 lib="single"
fi

##assigning index and genome sizes 
echo "trimming adapters"
echo $(date)
if [[ $refGen == "mm" ]]; then
 ref_size="/home/groups/oro/software/ChromHMM/CHROMSIZES/mm10.txt"
 blacklist="/oak/stanford/groups/oro/reference/latest_refGenome/Mus_musculus/UCSC/mm10/Annotation/mm10.blacklist.bed"
 ref_index="/oak/stanford/groups/oro/reference/latest_refGenome/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome"
else
 ref_size="/home/groups/oro/software/ChromHMM/CHROMSIZES/hg38.txt"
 blacklist="/oak/stanford/groups/oro/reference/latest_refGenome/Homo_sapiens/UCSC/hg38/Annotation/hg38.blacklist.bed"
 ref_index="/oak/stanford/groups/oro/reference/latest_refGenome/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome"
fi
echo $(date)

##create directories
if [ ! -d $out_map ]; then
 mkdir $out_map
fi
if [ ! -d $out_tracks ]; then
 mkdir $out_tracks
fi
if [ ! -d $out_peaks ]; then
 mkdir $out_peaks
fi

echo "genome size file:" $ref_size
echo "path to blacklist:" $blacklist
echo "reference index path:" $ref_index
echo "library type:" $lib

# 1. remove adapters
echo "trimming adapters"
echo $(date)

if [[ $lib == "single" ]]; then
 trim_galore  -q 10 --dont_gzip --fastqc $fw.fastq
else
 trim_galore -q 10 --dont_gzip --fastqc --paired ${fw}.fastq ${rw}.fastq
fi

echo $(date)

## 2. alignment
echo "Mapping to reference genome"
echo $(date)

if [[ $lib == "single" ]]; then
 bowtie2 -p 4 --very-sensitive -x ${ref_index} ${fw}_trimmed.fq -S $out_map/${samplename}.sam
 samtools view -S -b -F 4 -q 10 $out_map/${samplename}.sam | samtools sort -o $out_map/${samplename}_srt.bam
 rm ${fw}_trimmed.fq 
else
 bowtie2 -p 4 --very-sensitive -x ${ref_index} -1 ${fw}_val_1.fq -2 ${rw}_val_2.fq -S $out_map/${samplename}.sam
 samtools view -S -b -f 0x2 -q 10 $out_map/${samplename}.sam | samtools sort -o $out_map/${samplename}_srt.bam
 rm ${fw}_val_1.fq
 rm ${rw}_val_2.fq
fi

echo $(date)


## 3. remove PCR duplicates
echo "remove duplicates and create Index"
echo $(date)

samtools rmdup $out_map/${samplename}_srt.bam $out_map/${samplename}_srt_rmdup.bam

###Create Index file
samtools index $out_map/${samplename}_srt_rmdup.bam

#count remaining reads
samtools flagstat $out_map/${samplename}_srt_rmdup.bam > $out_map/mapping_stats_${samplename}_srt_rmdup

echo $(date)


## 4. Bam to bed conversion
echo "BamtoBed"
echo $(date)
bedtools bamtobed -i $out_map/${samplename}_srt_rmdup.bam > $out_map/${samplename}_srt_rmdup.bed
echo $(date)

## 5. Bed graph and bigwig
echo "Create tracks"
echo $(date)
echo "make bedGraph"

genomeCoverageBed -bg -split -i $out_map/${samplename}_srt_rmdup.bed  -g $ref_size > $out_tracks/${samplename}.bedGraph

echo "make normalized bedGraph"
perl /home/groups/oro/software/scripts/bedGraph_normalization.pl $out_tracks/${samplename}.bedGraph $out_tracks/${samplename}_norm.bedGraph &>  $out_tracks/${samplename}_norm.bedGraph.log

echo "create bigwig"
bedGraphToBigWig $out_tracks/${samplename}_norm.bedGraph ${ref_size} $out_tracks/${samplename}_norm.bw 

echo $(date)


## 6. Peak calling and black list removal
echo "peak calling using macs2"
echo $(date)

##remove black list regions
if [[ type=="TF" ]]; then
 macs2 callpeak -c $out_map/${controlname} -t $out_map/${samplename}_srt_rmdup.bed -f BED -g ${refGen} -n $out_peaks/${samplename} -p ${pvalue} -B
 bedtools intersect -v -a $out_peaks/${samplename}_peaks.narrowPeak -b ${blacklist} > $out_peaks/${samplename}_peaks.filterBL.bed
else
 macs2 callpeak -c $out_map/${controlname}  -t $out_map/${samplename}_srt_rmdup.bed -f BED -g ${refGen} -n $out_peaks/${samplename} -p ${pvalue} -B --broad --broad-cutoff ${pvalue}
 bedtools intersect -v -a $out_peaks/${samplename}_peaks.broadPeak -b ${blacklist} > $out_peaks/${samplename}_Peaks.filterBL.bed
fi

echo $(date)

###QC
echo "FRiP score"
# total reads
total_reads=$(samtools view -c $out_map/${samplename}_srt_rmdup.bam)

# reads in peaks
reads_in_peaks=$(bedtools sort -i $out_peaks/${samplename}_peaks.filterBL.bed | bedtools merge -i stdin | bedtools intersect -u -nonamecheck -a $out_map/${samplename}_srt_rmdup.bam -b stdin -ubam | samtools view -c)

# FRiP score
FRiP=$(awk "BEGIN {print "${reads_in_peaks}"/"${total_reads}"}")

echo ${samplename} "FRiP_score:" $FRiP >> $out_qc/FRiP_score.txt 


## 7. make tag directory
makeTagDirectory $out_map/${samplename}/ $out_map/${samplename}_srt_rmdup.bam

#make visualization file bigwig file
makeUCSCfile $out_map/${samplename}/ -o auto

## 8.Motif Analysis
##to find primary Motif, primary and co-enriched motif
echo "Motif analysis"
echo $(date)

if [[ $refGen == "mm" ]]; then
 findMotifsGenome.pl $out_peaks/${samplename}_peaks.filterBL.bed mm10 $out_motif/${samplename}_200size/ -size 200 -mask
else
 findMotifsGenome.pl $out_peaks/${samplename}_peaks.filterBL.bed hg38 $out_motif/${samplename}_200size/ -size 200 -mask
fi

#echo $(date)
