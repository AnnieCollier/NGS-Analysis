library("DiffBind")
library("tibble")
setwd("~/DiffBind")

#####
#Differential Binding Analysis
#####

#For DiffBind, you need to make a sample sheet to tell the program where your files are stored
samplesheet <- read.csv("H3K9me3_SampleSheet.csv", comment.char="#")
samplesheet
CTCF <- dba(sampleSheet=samplesheet)
#set peak size on each side of summit
CTCFcount <- dba.count(CTCF)
#set condition on which to test, change replicates since default is 3, add blocking factor to fix batch effects if needed
#CTCFcontrast <- dba.contrast(CTCFcount, categories=DBA_CONDITION, minMembers = 2, block=DBA_REPLICATE)
CTCFcontrast <- dba.contrast(CTCFcount, categories=DBA_CONDITION, minMembers = 2)
#CTCFcontrast2 <- dba.contrast(CTCFcount, group1=CTCFcount$masks$KO, group2=CTCFcount$masks$WT, name1="KO", name2="WT", block=DBA_REPLICATE)
CTCFanalyze<- dba.analyze(CTCFcontrast)
CTCFanalyze #tells you how many differentially expressed sites for each method used

#get significantly changed results, include all peaks, flip contrast groups if needed
CTCF.DB <- dba.report(CTCFanalyze, th=1, bFlip=TRUE, method=DBA_DESEQ2_BLOCK)
CTCFout<-as.data.frame(CTCF.DB)

vsd <- vst(CTCFanalyze)
plotPCA(vsd, "condition")

res_p = CTCFout[ CTCFout$FDR < 0.05,]
dim(res_p)
write.csv(res_p,file="CTCF.csv")
res_p_up=res_p[res_p$Fold > 1,]
res_p_up = as.data.frame(res_p_up)
dim(res_p_up)
head(res_p_up)
#convert data frame to bed file with unique identifiers
res_p_up <- tibble::rownames_to_column(res_p_up, "ID")
bed_up<-res_p_up[,1:4]
bed_up<-bed_up[c(2, 3, 4, 1)]
write.table(bed_up, file="CTCF_up.bed", sep="\t", row.names = FALSE, col.names = FALSE, quote= FALSE)
res_p_dn=res_p[res_p$Fold < -1,]
res_p_dn = as.data.frame(res_p_dn)
dim(res_p_dn)
res_p_dn <- tibble::rownames_to_column(res_p_dn, "ID")
bed_dn<-res_p_dn[,1:4]
bed_dn<-bed_dn[c(2, 3, 4, 1)]
write.table(bed_dn, file="CTCF_down.bed", sep="\t", row.names = FALSE, col.names = FALSE, quote= FALSE)



#####
#PLOTS
#####

gene=CTCFout

ggplot(data=gene)+
  geom_point(data=gene,mapping=aes(x=Conc, y=Fold),color="grey65",alpha=0.1)+
  geom_point(data=gene[gene$FDR<0.05 & gene$Fold>1,],mapping=aes(x=Conc, y=Fold),color="dodgerblue",alpha=0.5, size=3)+
  geom_point(data=gene[gene$FDR<0.05 & gene$Fold<= -1,],mapping=aes(x=Conc, y=Fold),color="firebrick2",alpha=0.5, size =3)+
  labs(title="WT vs AHDC1KO")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_hline(mapping= aes(yintercept=0),color = "red")

dba.plotMA(CTCFanalyze, method=DBA_DESEQ2_BLOCK, mapping=aes(x=Conc, y=Fold),color="grey65",alpha=0.1)

dba.plotMA(CTCFanalyze)
dba.plotPCA(CTCF)
pvals<-dba.plotBox(CTCF) #shows enrichement changes and separates them by up/down
pvals #returns p values from box plot
dba.plotMA(CTCFanalyze, method=DBA_DESEQ2_BLOCK)
dba.plotHeatmap(CTCF, method=DBA_DESEQ2_BLOCK)
dba.plotPCA(CTCF, method=DBA_DESEQ2_BLOCK, blog=TRUE, blind=FALSE)
dba.plotVolcano(CTCFanalyze, method=DBA_DESEQ2_BLOCK)
pvals<-dba.plotBox(CTCFanalyze, method=DBA_DESEQ2_BLOCK)
