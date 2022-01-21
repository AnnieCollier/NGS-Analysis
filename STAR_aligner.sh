ml biology
ml python/3.6.1
ml py-cutadapt/1.18_py36
ml samtools/1.8
ml bowtie2/2.3.4.1
ml trim_galore/0.5.0
ml bedtools
ml star
ml bowtie
ml fastx_toolkit/0.0.14

###parameters to set

path="/path/to/fastq/"

r1="file_prefix"
r2="" ### leave empty for single end 

sample="Sample_Name"

refGen="hs"  ### mm or hs

mkdir $sample

sequenceType="single" ### single or paired


if [[ $refGen == "mm" ]]; then
 ref_index="/path/to/mouse/STAR/index/"
 gtf="/path/to/mouse/STAR/index/genes.gtf"
else
 ref_index="/path/to/human/STAR/index/"
 gtf="/path/to/human/STAR/index/genes.gtf"
fi
#unzip your fastq.gz file, run QC on the file, and run genome alignment
if [[ $sequenceType == "single" ]]; then
 trim_galore -q 10 --fastqc $path/${r1}.fastq.gz
 gunzip ${r1}_trimmed.fq.gz
 STAR --runThreadN 12 --genomeDir $ref_index --sjdbGTFfile $gtf  --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --readFilesIn ${r1}_trimmed.fq --quantMode GeneCounts --outFileNamePrefix ${sample}/
else
 trim_galore -q 10 --fastqc --paired $path/${r1}.fastq.gz $path/${r2}.fastq.gz
 gunzip ${r1}_val_1.fq.gz
 gunzip ${r2}_val_2.fq.gz
 STAR --runThreadN 12 --genomeDir $ref_index --sjdbGTFfile $gtf --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif  --outFilterIntronMotifs RemoveNoncanonical --readFilesIn ${r1}_val_1.fq ${r2}_val_2.fq --quantMode GeneCounts --outFileNamePrefix ${sample}/
fi

#make the name prettier
mv $sample/Aligned.sortedByCoord.out.bam $sample/${sample}.bam

#make sure you have homer in your path
#convert to  sam file if desired
samtools view $sample/${sample}.bam > $sample/${sample}.sam
#make a tag directory if desired
makeTagDirectory ${sample}/${sample}/ $sample/${sample}.sam -format sam
#make a UCSC file for visualization if desired
makeUCSCfile ${sample}/${sample}/ -fragLength given -o auto
