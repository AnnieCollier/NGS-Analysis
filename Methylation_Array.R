####################
# Illumina EPIC Methylation Array analysis using Minfi, Limma, and Quantro
# Author: Annie Collier
####################


############
### Step 1 : Install packages and load IDAT files
############

#set your working directory where you have the IDAT folder
setwd("~/Box/Ann_Collier/Methylation/METHYLATION")
#load all packages. Install first if necessary
library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(quantro)
library(doParallel)
library(ggplot2)

#load the hg19 annotation file (no hg38 exists yet, will need to do a liftover at the end if needed)
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
head(anno)
#specify the path to the iDat files and filled out sample sheet
#sample sheet template can be found at https://support.illumina.com/downloads/infinium-methylationepic-sample-sheet.html
#fill out samples sheet and put in folder with IDAT files
dataDirectory <- file.path("IDAT")
list.files(dataDirectory, recursive = TRUE)
targets <- read.metharray.sheet(dataDirectory, pattern="SampleSheet.csv")
targets
#Subset if you only want to look at specific samples
#targets <- targets[1:4,]
rgSet <- read.metharray.exp(targets=targets)
rgSet
targets$ID <- paste(targets$Sample_Name,sep=".")
sampleNames(rgSet) <- targets$ID
rgSet

############
### Step 2 : QC
############

phenoData = pData(rgSet)
sampleID = phenoData$Sample_Name

MSet = preprocessRaw(rgSet)
ratioSet = ratioConvert(MSet, what="both", keepCN=TRUE)
gset = mapToGenome(ratioSet)

#Calculate detection (signal-to-background) p-values
#plot the means of each sample to find any samples that are too high (> 0.05)
detP <- detectionP(rgSet)
head(detP)
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2, 
        cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal,
       bg="white")
barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2, 
        cex.names=0.8, ylim=c(0,0.002), ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal, 
       bg="white")
#generate a full QC report
qcReport(rgSet, sampNames=targets$ID, sampGroups=targets$Sample_Group, 
         pdf="qcReport.pdf")
#plot QC of all samples
qc = getQC(MSet)
rownames(qc) = sampleID
fit = lm(qc$uMed~qc$mMed)
plotQC(qc)

#remove any samples above a certain p-value cutoff (IF NEEDED)
keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]
rgSet


############
### Step 3 : Normalization, M/Beta Values, and filtering of low quality probes
############

#SWAN NORMALIZATION (use if quantro p-value is < 0.05 AND/OR if you expect global changes between samples...like different cell types or KO cells)

MSet.swan = preprocessSWAN(rgSet)
ratioSet.swan = ratioConvert(MSet.swan, what="both", keepCN=TRUE)
gset.swan = mapToGenome(ratioSet.swan)
#remove low quality probes
detP <- detP[match(featureNames(gset.swan),rownames(detP)),] 
keep <- rowSums(detP < 0.05) == ncol(gset.swan) 
table(keep)
gset.swan <- gset.swan[keep,]
gset.swan
plotMDS(getM(gset.swan), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)])
#calculate and plot M and Beta Values
mVals <- getM(gset.swan)
bVals <- getBeta(gset.swan)
par(mfrow=c(1,2))
densityPlot(bVals, sampGroups=targets$Sample_Group, main="Beta values", 
            legend=FALSE, xlab="Beta values")
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(mVals, sampGroups=targets$Sample_Group, main="M-values", 
            legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))
#filter out SNPs
snps.swan = getSnpInfo(gset.swan)
gset.swan = addSnpInfo(gset.swan)
head(granges(gset.swan))
gset.swan = dropLociWithSnps(gset.swan, snps=c("SBE","CpG","Probe"), maf=0)
getAnnotationObject(gset.swan)
beta.swan = getBeta(gset.swan)
colnames(beta.swan) = sampleID
#plot Beta values before and after normalization to see if the normalization helped
par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=targets$Sample_Group,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(gset.swan), sampGroups=targets$Sample_Group,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))
#bean plot
par(mfcol=c(1,1),mar=c(2,10,1,1), oma=c(1,1,0,0),lwd=1)
densityBeanPlot(MSet.swan, sampGroups = phenoData$Sample_Group, main="", sampNames = phenoData$Sample_ID, pal=hmcol)
mtext("Beta values", side=1, line=3.5, cex=1.5)


#########
# DIFFERENTIAL METHYLATION
#########

#probe-wise differential analysis (DMP's)
#specify factor of interest (cell type, patient, etc.)
#in this example the factors are genotypes, specified by names from sample sheet
cellType <- factor(targets$Sample_Group, levels = c("WT_D7","AHDC1KO_D7","GATA3KO_D7","WT_D0"))
#if you have a second parameter, like patient, specify as a seprate vector
design <- model.matrix(~0+cellType)
#to add a second parameter do this:
#design <- model.matrix(~0+cellType+individual)
colnames(design) <- c(levels(cellType))
fit <- lmFit(mVals, design)
#decide which contrasts you want to make
contMatrix <- makeContrasts(WT_D7vsAHDC1KO_D7=WT_D7-AHDC1KO_D7, WT_D7vsGATA3KO_D7=WT_D7-GATA3KO_D7, WT_D7vsWT_D0=WT_D7-WT_D0, levels=design)
contMatrix
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
summary(decideTests(fit2))
top<-topTable(fit2,coef=1)
#pull the results table for each contrast by changing the coef number
DMP_AHDC1 <- topTable(fit2, num=Inf, coef=1, genelist=anno)
DMP_GATA3 <- topTable(fit2, num=Inf, coef=2, genelist=anno)
DMP_D0 <- topTable(fit2, num=Inf, coef=3, genelist=anno)
#export results to excel files
write.table(DMP_AHDC1, file="DMP_AHDC1_SWAN.csv", sep=",",row.names=FALSE)
write.table(DMP_GATA3, file="DMP_GATA3_SWAN.csv", sep=",",row.names=FALSE)
write.table(DMP_D0, file="DMP_D0_SWAN.csv", sep=",",row.names=FALSE)
write.table(bVals, file="bVals.csv", sep=",",row.names=FALSE)
#volcano plot of all probes for the first contrast
ggplot(data=DMP_AHDC1)+
  geom_point(data=DMP_AHDC1,mapping=aes(x=logFC, y=-log10(P.Value)),color="grey65",alpha=0.1)+
  geom_point(data=DMP_AHDC1[DMP_AHDC1$P.Value<0.01 & DMP_AHDC1$logFC>1,],mapping=aes(x=logFC, y=-log10(P.Value)),color="dodgerblue",alpha=0.1)+
  geom_point(data=DMP_AHDC1[DMP_AHDC1$P.Value<0.01 & DMP_AHDC1$logFC<= -1,],mapping=aes(x=logFC, y=-log10(P.Value)),color="firebrick2",alpha=0.1)+
  labs(title="WT vs AHDC1KO")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))
#plot the top 4 most significantly differentially methylated CpGs 
par(mfrow=c(2,2))
sapply(rownames(DMPs)[1:4], function(cpg){
  plotCpg(bVals, cpg=cpg, pheno=targets$Sample_Group, ylab = "Beta values")})

###Differentially regulated regions (DMR's)
myAnnotation <- cpg.annotate(object = mVals, datatype = "array", what = "M", 
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix, 
                             coef = "WT_D7vsAHDC1KO_D7", arraytype = "EPIC")
str(myAnnotation)
#combine the CPGs into regions
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs)
results.ranges
write.table(results.ranges, file="DMRs_SWAN.csv", sep=",", row.names=FALSE)
#plot results
groups <- pal[1:length(unique(targets$Sample_Group))]
names(groups) <- levels(factor(targets$Sample_Group))
cols <- groups[as.character(factor(targets$Sample_Group))]
par(mfrow=c(1,1))
DMR.plot(ranges = results.ranges, dmr = 6765, CpGs = bVals, phen.col = cols, 
         what = "Beta", arraytype = "EPIC", genome = "hg19")
# indicate which genome is being used
gen <- "hg19"
# the index of the DMR that we will plot 
dmrIndex <- 1
# extract chromosome number and location from DMR results 
chrom <- as.character(seqnames(results.ranges[dmrIndex]))
start <- as.numeric(start(results.ranges[dmrIndex]))
end <- as.numeric(end(results.ranges[dmrIndex]))
# add 25% extra space to plot
minbase <- start - (0.25*(end-start))
maxbase <- end + (0.25*(end-start))

enrichment_GO <- goregion(results.ranges[1:100], all.cpg = rownames(targets$Sample_Group),
                          collection = "GO", array.type = "EPIC")
enrichment_GO <- enrichment_GO[order(enrichment_GO$P.DE),]
head(as.matrix(enrichment_GO), 10)


# Get the significant CpG sites at less than 5% FDR
sigCpGs <- DMPs$Name[DMPs$adj.P.Val<0.05]

# First 10 significant CpGs
sigCpGs[1:10]
# Total number of significant CpGs at 5% FDR
length(sigCpGs)
# Get all the CpG sites used in the analysis to form the background
all <- DMPs$Name
# Total number of CpG sites tested
length(all)
par(mfrow=c(1,1))
gst <- gometh(sig.cpg=sigCpGs, all.cpg=all, plot.bias=TRUE)

#########
# ONTOLOGY ANALYSIS
#########

# Top 10 GO categories
topGSA(gst, number=10)
# load Broad human curated (C2) gene sets
load(paste(dataDirectory,"human_c2_v5p2.rdata",sep="/"))
# perform the gene set test(s)
gsa <- gsameth(sig.cpg=sigCpGs, all.cpg=all, collection=Hs.c2)
# top 10 gene sets
topGSA(gsa, number=10)
