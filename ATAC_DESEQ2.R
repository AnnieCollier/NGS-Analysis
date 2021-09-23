library('DESeq2')
library(limma)
library("edgeR")
library("geneplotter")
library("ggpubr")
library(dplyr)
library(plyr)
library(reshape)
library(ggplot2)
library(tidyr)

#use the union peak file made with union_macs.sh script
targets <- read.delim("AC_pipeline_AHDC1_GATA3_hg38_union_peak_count.txt", stringsAsFactors=FALSE, header=FALSE)
head(targets)
targets$names<-paste0(targets$V1,":",targets$V2,"-",targets$V3)
targets=rename(targets, c("V1"="chr","V2"="Start","V3"="end","V4"="WT_1","V5"="WT_2","V6"="AHDC1_1","V7"="AHDC1_2","V8"="GATA3_1","V9"="GATA3_2"))
dat1=targets[,c(4,5,8,9)]
head(dat1)
dat1<-as.matrix(dat1)
row.names(dat1)<-targets$names
storage.mode(dat1)="integer"
head(dat1)
sample1 <- data.frame(row.names =colnames(dat1[,c(1,2,3,4,5,6,7,8,9,10,11,12)]),condition=factor(c("WT","WT","AHDC1","AHDC1","GATA3","GATA3","p63","p63","GRHL2","GRHL2","TFAP2A","TFAP2A")))
sample2<-data.frame(row.names =colnames(dat1[,c(1,2,3,4)]),condition=c("WT","WT","GATA3","GATA3"))

dds <- DESeqDataSetFromMatrix(countData = dat1, colData=sample2, design=~condition)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
dds <- DESeq(dds)
dds_res <- results(dds)
dds_res <- dds_res[complete.cases(dds_res),] 

rld <- rlog( dds_AHDC1 )
plotPCA( rld, intgroup = c( "condition") )
plotCounts(TFAP2A, intgroup = c("batch","condition"))

dds_res = as.data.frame(dds_res)
#res_n <- merge(names,dds_res, by=0)
res_p = dds_res[ dds_res$padj < 0.05,]
dim(res_p)
#write.csv(res_p,file="CTCF.csv")
res_p_up=res_p[res_p$log2FoldChange > 1,]
res_p_up = as.data.frame(res_p_up)
dim(res_p_up)
head(res_p_up)
#write.csv(res_p_up,file="AHDC1_hg38_down_AC.csv")
res_p_up <- tibble::rownames_to_column(res_p_up, "ID")
res_p_up <- separate(data = res_p_up, col = ID, into = c("chr", "Start"), sep = ":")
res_p_up <- separate(data = res_p_up, col = Start, into = c("start", "stop"), sep = "-")
res_p_up <- tibble::rowid_to_column(res_p_up, "ID")
bed_up<-res_p_up[,1:4]
bed_up<-bed_up[c(2, 3, 4, 1)]
write.table(bed_up, file="GATA3_ATAC_down_hg38_newBL_ACpipeline_0.05.bed", sep="\t", row.names = FALSE, col.names = FALSE, quote= FALSE)

AHDC1 <-results(dds, contrast=c("condition","GATA3","WT"))
#AHDC1 <- AHDC1[complete.cases(AHDC1),] 
AHDC1<-as.data.frame(AHDC1)
#AHDC1 <- AHDC1[complete.cases(AHDC1),]
gene=AHDC1

ggplot(data=gene)+
  geom_point(data=gene,mapping=aes(x=log(baseMean), y=log2FoldChange),color="grey65",alpha=0.1)+
  geom_point(data=gene[gene$padj<0.05 & gene$log2FoldChange>1,],mapping=aes(x=log(baseMean), y=log2FoldChange),color="dodgerblue",alpha=0.5, size=3)+
  geom_point(data=gene[gene$padj<0.05 & gene$log2FoldChange<= -1,],mapping=aes(x=log(baseMean), y=log2FoldChange),color="firebrick2",alpha=0.5, size =3)+
  labs(title="WT vs AHDC1KO")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_hline(mapping= aes(yintercept=0),color = "red")


#run through limma removebatcheffect if you need batch correction
design<- model.matrix(~sample1$condition)
ncounts<-counts(dds)
a<- log2(ncounts[,1:4])
#limma remove batch effect
library(limma)
b<-removeBatchEffect(a, batch=c("1","2","1","2"), design = design)
c<- 2**b[,1:4]
d<- round(c[,1:4])
dds2<-DESeqDataSetFromMatrix(countData=d,colData=sample1,design=~condition)
dds2<-DESeq(dds2)
rld2 <- rlog( dds2 )
plotPCA( rld2, intgroup = c( "condition") )   

