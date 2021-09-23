library('DESeq2')
library(limma)
library("edgeR")
library("geneplotter")
library("ggpubr")
library(dplyr)
library(plyr)
library(reshape)
library(ggplot2)

#input your files and format the table so that it will work properly with DESEQ2
#MAKE SURE TO USE THE RAW FILE, NOT RPKM
files<- file.path("RNAseq_AHDC1KO_condenseGenes_raw.txt")
dat=read.table(files , quote="", fill=T, header=T, sep="\t")
sampleTable<- data.frame(condition=factor(c("KO", "KO","WT","WT")))
names<- dat[,0:8]
d <- dat[, 9:12]
names(d)<- c("KO_1", "KO_2","WT_1", "WT_2")
d<-as.matrix(d)
storage.mode(d)="integer"
d<-na.omit(d)
#make a DESEQ object
dds <- DESeqDataSetFromMatrix(countData = d, colData = sampleTable, design = ~condition)
dds <- DESeq(dds)
dds_res <- results(dds)
#PCA plot of all samples
rld <- rlog( dds )
plotPCA( rld, intgroup = c( "condition") )

#make a data frame of your results, subset into different gene lists, and make an MA plot
AHDC1 <-results(dds, contrast=c("condition","AHDC1","WT"))
AHDC1<-as.data.frame(AHDC1)
AHDC1<-merge(AHDC1, names, by = 0)
AHDC1$gene_name<-sub("\\|..*", "", AHDC1$Annotation.Divergence)
AHDC1 <- AHDC1[complete.cases(AHDC1),]
AHDC1 <- AHDC1[order(AHDC1$padj),] 
AHDC1_sig = AHDC1[AHDC1$padj < 0.05,]
AHDC1_up = AHDC1_sig[AHDC1_sig$log2FoldChange > 1,]
AHDC1_down = AHDC1_sig[AHDC1_sig$log2FoldChange <= -1,]
write.csv(AHDC1_down,file="AHDC1_down.csv")
write.csv(AHDC1_up,file="AHDC1_up.csv")

gene=AHDC1
ggplot(data=gene)+
  geom_point(data=gene,mapping=aes(x=log(baseMean), y=log2FoldChange),color="grey65",alpha=0.1)+
  geom_point(data=gene[gene$padj<0.05 & gene$log2FoldChange>1,],mapping=aes(x=log(baseMean), y=log2FoldChange),color="dodgerblue",alpha=0.5, size=3)+
  geom_point(data=gene[gene$padj<0.05 & gene$log2FoldChange<= -1,],mapping=aes(x=log(baseMean), y=log2FoldChange),color="firebrick2",alpha=0.5, size =3)+
  labs(title="WT vs AHDC1KO")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_hline(mapping= aes(yintercept=0),color = "red")
 
 
#IF YOU NEED TO REMOVE BATCH EFFECT YOU CAN ADD THIS CODE:
design<- model.matrix(~sampleTable$condition)
ncounts<-counts(dds)
a<- log2(ncounts[,1:10])
#limma remove batch effect
library(limma)
b<-removeBatchEffect(a, batch=c("1","1","2","2","2","2","2","2","2","2"), design = design)
c<- 2**b[,1:10]
d<- round(c[,1:10])
