## miR-71 data, E. multilocularis (Matias Perez, Mara Rosenzvit)


##########
## LOAD ##
##########

library(gridExtra)
#library(ggsci)

library(ggplot2)
library(DESeq2)
#library(pheatmap)

#library(pcaExplorer)
#library(FactoMineR)
#library(factoextra)



##########
## MAIN ##
##########


## load gene count matrix
load("miR71_gene_count.RData")
head(counts)
dim(counts)
#[1] 10663    18

## load phenotype
pheno <- read.csv("Phenodata.csv", header=FALSE)
rownames(pheno) <- pheno[,1]
colnames(pheno) <- c("sample", "treatment", "type")


## check columns in both table
all(rownames(pheno) %in% colnames(counts))
#[1] TRUE
## reorder columns in counts, given rows in phenodata
counts <- counts[, rownames(pheno)]
## check identity
all(rownames(pheno) == colnames(counts))
#[1] TRUE



## build dds
dds <- DESeqDataSetFromMatrix(countData=counts, colData=pheno, design= ~ treatment)
dds
#class: DESeqDataSet 
#dim: 10663 18 
#metadata(1): version
#assays(1): counts
#rownames(10663): EmuJ_001059300 EmuJ_000212400 ... EmuJ_002178900
#  EmuJ_000425400
#rowData names(0):
#colnames(18): sLNA1 sLNA2 ... sCSE2 sCSE3
#colData names(2): sample treatment  type


## Pre-Filtering
keep <- rowSums(counts(dds)) >=10
sum(keep)
#[1] 9703
dds <- dds[keep,]
dds
#class: DESeqDataSet 
#dim: 9703 18 
#...



## Differential expression analysis ##
dds <- DESeq(dds)
dds
#class: DESeqDataSet 
#dim: 9703 18 
#metadata(1): version
#assays(3): counts mu cooks
#rownames(9703): EmuJ_001059300 EmuJ_000784100 ... EmuJ_000690600
#  EmuJ_002178900
#rowData names(37): baseMean baseVar ... deviance maxCooks
#colnames(18): sLNA1 sLNA2 ... sCSE2 sCSE3
#colData names(3): sample treatment type sizeFactor

summary(rowSums(counts(dds)))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#     10    1756   15748   33966   40820 3384235 


dds$sizeFactor
#    sLNA1     sLNA2     sLNA3    sLNAS1    sLNAS2    sLNAS3    s2OME1    s2OME2 
#1.0014575 1.1874231 1.0520471 0.9579565 1.0672567 0.9901200 0.9691030 0.8853602 
#   s2OME3   s2OMES1   s2OMES2   s2OMES3      sCE1      sCE2      sCE3     sCSE1 
#0.8910831 0.8387084 1.0168616 1.1039658 0.8996085 1.0918279 1.1131470 1.0637019 
#    sCSE2     sCSE3 
#0.9972495 0.9553987 


resultsNames(dds)
#[1] "Intercept"                  "treatment_2OMESi_vs_2OMESc"
#[3] "treatment_CEle_vs_2OMESc"   "treatment_CSele_vs_2OMESc" 
#[5] "treatment_LNASc_vs_2OMESc"  "treatment_LNASi_vs_2OMESc" 




#######################
## PCA - ALL SAMPLES ##
#######################

## transforming data (blind dispersion = TRUE / default)
rld <- rlog(dds)
p.rld <- meanSdPlot(assay(rld))
## this gives log2(n + 1)
ntd <- normTransform(dds)
p.ntd <- meanSdPlot(assay(ntd))

pdf("meanSdPlot_Emultilocularis.pdf", width=10, height=5)
grid.arrange(p.ntd$gg  + ggtitle("normTransform"), p.rld$gg + ggtitle("rlog"), nrow=1)
dev.off()


#library(pcaExplorer)
#pcaExplorer(dds=dds, rlt=rld)


## plotPCA
plotPCA(rld, intgroup=c("treatment"), ntop=500)

my.colors <- c("coral1","coral1","lightsteelblue3","darkgrey","darkcyan","darkcyan")
names(my.colors) <- c("2OMESc","2OMESi","CEle","CSele","LNASc","LNASi")

my.shapes <- c("square","triangle","circle")
names(my.shapes) <- c("electro","scrambled","silenced")

pcaData <- plotPCA(rld, intgroup=c("treatment","type"), ntop=500, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p.pca <- ggplot(pcaData, aes(PC1, PC2, color=treatment, shape=type)) +
  geom_point(size=5) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + theme_bw() + scale_colour_manual(values=my.colors, labels=c("20MESc","20MESi","CEle","CSele","LNASc","LNASi")) + scale_shape_manual(values=my.shapes, labels=c("electro","scrambled","silenced"))
p.pca


## PCA (FactoMineR, factoextra)
##get the 500 leading genes (with highest row variance)
myvar <- sort(apply(assay(rld),1,var), decreasing=TRUE)
names(myvar[1:500])
res.pca <- PCA(t(assay(rld)[names(myvar[1:500]), ]), scale.unit=FALSE, ncp=5, graph=FALSE)

get_eigenvalue(res.pca)
#       eigenvalue variance.percent cumulative.variance.percent
#Dim.1  23.2750662       25.9256467                    25.92565
#Dim.2  15.9077160       17.7192976                    43.64494
#Dim.3   9.9257518       11.0561032                    54.70105
#Dim.4   7.5076660        8.3626442                    63.06369


pca.ell <- fviz_pca_ind(res.pca, geom.ind="point", # show points only (nbut not "text")
pointsize=6, 
col.ind = pheno$treatment, mean.point = FALSE, # color by groups
palette= c("coral1","coral1","lightsteelblue3","darkgrey","darkcyan","darkcyan"),
addEllipses = TRUE, ellipse.type="confidence", ellipse.level=0.9, ellipse.alpha=0.1,# Concentration ellipses
legend.title = "treatment", title='') + scale_shape_manual(values=c(19,17,15,15,19,17)) + xlab("PC1: 25.93% variance") + ylab("PC2: 17.72% variance") + theme_bw()
print(pca.ell)


##check sample id
pca.elln <- fviz_pca_ind(res.pca, geom.ind=c("point","text"), 
pointsize=6, 
col.ind = pheno$treatment, mean.point = FALSE, # color by groups
palette= c("coral1","coral1","lightsteelblue3","darkgrey","darkcyan","darkcyan"),
addEllipses = TRUE, ellipse.type="confidence", ellipse.level=0.9, ellipse.alpha=0.1,# Concentration ellipses
legend.title = "treatment", title='') + scale_shape_manual(values=c(19,17,15,15,19,17)) + xlab("PC1: 25.93% variance") + ylab("PC2: 17.72% variance") + theme_bw()
print(pca.elln)




######################
## DE: LNASivsLNASc ##
######################

## LNASivsLNASc
res.LNA <- results(dds, contrast=c("treatment","LNASi","LNASc"))

res.LNA

summary(res.LNA$padj)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.0000  0.9998  0.9998  0.9984  0.9998  1.0000      48

summary(res.LNA$log2FoldChange)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-31.91247  -0.10327  -0.00589  -0.00197   0.08372  35.78164 

summary(res.LNA)
#out of 9703 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 8, 0.082% 
#LFC < 0 (down)   : 4, 0.041% 
#outliers [1]     : 48, 0.49% 
#low counts [2]   : 0, 0% 
#(mean count < 0)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

## have a look on a combination of padj and log2FoldChange
sum(res.LNA$padj <=0.1 & abs(res.LNA$log2FoldChange) >=1, na.rm=TRUE)
#[1] 12
sum(res.LNA$padj <=0.1, na.rm=TRUE)
#[1] 12

genesDE.LNA <- rownames(res.LNA[res.LNA$padj<=0.1 & !is.na(res.LNA$padj),])
length(genesDE.LNA)
#[1] 12


## Log fold change shrinkage for visualization and ranking 
res.LNA.normal <- lfcShrink(dds, contrast=c("treatment","LNASi","LNASc"), type="normal")

## plotMA with both FoldChange treatments (raw y shrinkage "normal")
par(mfrow=c(1,2), mar=c(4,4,2,1))
xlim <- c(1,1e6) ; ylim <-c(-10,10)
plotMA(res.LNA, xlim=xlim, ylim=ylim, main="raw")
plotMA(res.LNA.normal, xlim=xlim, ylim=ylim, main="normal")

## log2FC summary
summary(res.LNA[genesDE.LNA,]$log2FoldChange)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -27.99  -24.06   23.47   11.13   33.29   35.78
## shrunken log2FC summary (normal)
summary(res.LNA.normal[genesDE.LNA,]$log2FoldChange)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-0.009633 -0.006551  0.012841  0.007743  0.015822  0.024963


## write table LNA
res.out <- as.data.frame(res.LNA)
res.out$normal <- as.vector(res.LNA.normal$log2FoldChange)
res.out$normal.pa <- as.vector(res.LNA.normal$padj)

write.table(res.out,"genes_LNA_log2FC_and_normal-shrinkage.csv", quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE, dec=",")



########################
## DE: 2OMESivs2OMESc ##
########################


##2OMESivs2OMESc
res.2OME <- results(dds, contrast=c("treatment","2OMESi","2OMESc"))
res.2OME

summary(res.2OME$padj)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.0000  0.3285  0.6864  0.5991  0.8989  0.9999      48

summary(res.2OME$log2FoldChange)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-27.31269  -0.11402   0.01022   0.19684   0.20597  27.79717

summary(res.2OME)
#adjusted p-value < 0.1
#LFC > 0 (up)     : 744, 7.7% 
#LFC < 0 (down)   : 548, 5.6% 
#outliers [1]     : 48, 0.49% 
#low counts [2]   : 0, 0% 
#(mean count < 0)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

## have a look on a combination of padj and log2FoldChange
sum(res.2OME$padj <=0.1 & abs(res.2OME$log2FoldChange) >=1, na.rm=TRUE)
#[1] 309
sum(res.2OME$padj <=0.1, na.rm=TRUE)
#[1] 1292

genesDE.2OME <- rownames(res.2OME[res.2OME$padj<=0.1 & !is.na(res.2OME$padj),])
length(genesDE.2OME)
#[1] 1292


## Log fold change shrinkage for visualization and ranking 
res.2OME.normal <- lfcShrink(dds, contrast=c("treatment","2OMESi","2OMESc"), type="normal")

## plotMA with both FoldChange treatments (raw y shrinkage "normal")
par(mfrow=c(1,2), mar=c(4,4,2,1))
xlim <- c(1,1e6) ; ylim <-c(-10,10)
plotMA(res.2OME, xlim=xlim, ylim=ylim, main="raw")
plotMA(res.2OME.normal, xlim=xlim, ylim=ylim, main="normal")


## log2FC summary
summary(res.2OME[genesDE.2OME,]$log2FoldChange)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-27.3127  -0.2636   0.2991   0.4840   0.8538  25.4147
## shrunken log2FC summary (normal)
summary(res.2OME.normal[genesDE.2OME,]$log2FoldChange)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-1.0039 -0.2471  0.2666  0.2069  0.5262  2.2796


## write table DE for 2OME
res.out <- as.data.frame(res.2OME)
res.out$normal <- as.vector(res.2OME.normal$log2FoldChange)
res.out$normal.pa <- as.vector(res.2OME.normal$padj)

write.table(res.out,"genes_2OME_log2FC_and_normal-shrinkage.csv", quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE, dec=",")
## see excel genes_2OME_log2FC_and_normal-shrinkage.xlsx


## VOLCANO PLOT

topT <- as.data.frame(res.2OME.normal)

pdf("volcano_shrinkage-normal_2OME.pdf", width=10, height=10)
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)

with(topT, plot(log2FoldChange, -log10(padj), pch=20, cex=1.0, xlab=bquote(~log[2]~(fold~change)), ylab=bquote(~-log[10]~(adjusted~p~value))))
with(subset(topT, padj<0.1 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.9))

#Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.1
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-1, col="black", lty=4, lwd=1.5)
abline(v=1, col="black", lty=4, lwd=1.5)
abline(h=-log10(0.1), col="black", lty=4, lwd=1.5)
dev.off()

dim(subset(topT, padj<0.1 & abs(log2FoldChange)>1))
#[1] 94  6
#94 son los genes rojos del volcano plot



## VOLCANO PLOT NEW FoldChange
## log2(1,4)= 0.485  ==>> 0.48

pdf("volcano_shrinkage-normal_2OME_NEW.pdf", width=10, height=10)
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)

with(topT, plot(log2FoldChange, -log10(padj), pch=20, cex=1.0, xlab=bquote(~log[2]~(fold~change)), ylab=bquote(~-log[10]~(adjusted~p~value))))
with(subset(topT, padj<0.1 & abs(log2FoldChange)>0.48), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.9))

#Add lines for log2(1,4)= 0.485  ==>> 0.48 and P-value cut-off at FDR 0.1
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-0.48, col="black", lty=4, lwd=1.5)
abline(v=0.48, col="black", lty=4, lwd=1.5)
abline(h=-log10(0.1), col="black", lty=4, lwd=1.5)
dev.off()

dim(subset(topT, padj<0.1 & abs(log2FoldChange)>0.48))
#[1] 416  6




#####################
## DE: CElevsCSele ##
#####################


##CElevsCSele
res.CE <- results(dds, contrast=c("treatment","CEle","CSele"))
res.CE


summary(res.CE$padj)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.0000  0.0664  0.3557  0.4051  0.7178  0.9993     425

summary(res.CE$log2FoldChange)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-30.2990  -0.1999  -0.0035   0.0516   0.2078  36.5538 

summary(res.CE)
#out of 9703 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 1283, 13% 
#LFC < 0 (down)   : 1412, 15% 
#outliers [1]     : 48, 0.49% 
#low counts [2]   : 377, 3.9% 
#(mean count < 2)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

## have a look on a combination of padj and log2FoldChange
sum(res.CE$padj <=0.1 & abs(res.CE$log2FoldChange) >=1, na.rm=TRUE)
#[1] 374
sum(res.CE$padj <=0.1, na.rm=TRUE)
#[1] 2695

genesDE.CE <- rownames(res.CE[res.CE$padj<=0.1 & !is.na(res.CE$padj),])
length(genesDE.CE)
#[1] 2695



## Log fold change shrinkage for visualization and ranking 
res.CE.normal <- lfcShrink(dds, contrast=c("treatment","CEle","CSele"), type="normal")

## plotMA with both FoldChange treatments (raw y shrinkage "normal")
par(mfrow=c(1,2), mar=c(4,4,2,1))
xlim <- c(1,1e6) ; ylim <-c(-10,10)
plotMA(res.CE, xlim=xlim, ylim=ylim, main="raw")
plotMA(res.CE.normal, xlim=xlim, ylim=ylim, main="normal")


## log2FC summary
summary(res.CE[genesDE.CE,]$log2FoldChange)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-25.2344  -0.4245  -0.1889   0.1464   0.3910  36.5538
## shrunken log2FC summary (normal)
summary(res.CE.normal[genesDE.CE,]$log2FoldChange)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-1.52250 -0.36883 -0.17649 -0.02941  0.31941  1.91292


## table DE with CE treatment
res.out <- as.data.frame(res.CE)
res.out$normal <- as.vector(res.CE.normal$log2FoldChange)
res.out$normal.pa <- as.vector(res.CE.normal$padj)

write.table(res.out,"genes_CE_log2FC_and_normal-shrinkage.csv", quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE, dec=",")




