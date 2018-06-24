# run this if it's your first time using it to install
# Make sure you the Plots window is maximised, otherwise errors will occur --> figures can be exported to an external file.
# install.packages(c("RColorBrewer", "mixOmics", "VennDiagram", "gplots", "stringi", "rstudioapi"))
# source("https://bioconductor.org/biocLite.R") # install bioconductor if not already installed
# install.packages(c("Rcpp", "httpuv", "shiny"))
# biocLite("HTSFilter") 
# biocLite("DBI")
library(limma)
library(edgeR) # for DGEList
library(RColorBrewer)
library(MASS)
library(lattice)
library(ggplot2)
library(mixOmics)
library(gplots)
library(DBI)
library(HTSFilter)
library(mixOmics)
#############Create a function set_wd that sets the working directory############## 
set_wd <- function() {
  # load rstudioapi # make sure you have it installed
  library(rstudioapi)
  # the following line is for getting the path of your current open file
  current_path <- getActiveDocumentContext()$path
  # The next line set the working directory to the relevant one:
  setwd(dirname(current_path ))
  # you can make sure you are in the right directory
  print( getwd() )
}

# Execute set_wd
set_wd()

################################Data#######################################
# Importing the files
# data <- read.delim(file.choose(), header=T)
# Make sure the files are on the scripts location.
WCFS1_anno <- read.delim("WCFS1_anno.txt", header=T, sep="\t")
RNA_Seq_counts <- read.delim("RNA_Seq_counts.txt", header=T, sep="\t")

# Have a look at the count data:
head(RNA_Seq_counts)
nrow(RNA_Seq_counts)

# merge the data
merged_data = merge(RNA_Seq_counts, WCFS1_anno, by.x="ID", by.y="ORF")

# remove NA values
merged_data <- merged_data[1:(length(merged_data)-4)]

#
row.names(merged_data) <- merged_data[,"ID"]

#############################Starting from count table##############################
# create DGEList
exp_WCFS1 <- c("WCFS1.glc","WCFS1.glc","WCFS1.rib","WCFS1.rib")
exp_NC8 <- c("NC8.glc", "NC8.glc", "NC8.rib", "NC8.rib")
group_WCFS1 <- factor(exp_WCFS1)
group_NC8 <-factor(exp_NC8)
y_WCFS1 <- DGEList(counts=merged_data[,2:5],group=group_WCFS1)
y_NC8 <- DGEList(counts=merged_data[,6:9],group=group_NC8)

#######################Data exploration and quality assessment######################
# Extract RNA_Seq_Counts
RNASeqCounts_WCFS1 <- log2(y_WCFS1$counts+1)
RNASeqCounts_NC8 <- log2(y_NC8$counts+1)
head(RNASeqCounts_WCFS1)
head(RNASeqCounts_NC8)

# Set the graphical parameters
par(mfrow=c(2,2))
# Create histogram for Counts # lijkt ongeveer normaal verdeeld, counts tegenover sample
# Note that these are plots before normalisation, for data exploration and quality assesment.
hist(RNASeqCounts_WCFS1[,"WCFS1.glc.1"])
hist(RNASeqCounts_WCFS1[,"WCFS1.glc.2"])
hist(RNASeqCounts_WCFS1[,"WCFS1.rib.1"])
hist(RNASeqCounts_WCFS1[,"WCFS1.rib.2"])
hist(RNASeqCounts_NC8[,"NC8.glc.1"])
hist(RNASeqCounts_NC8[,"NC8.glc.2"])
hist(RNASeqCounts_NC8[,"NC8.rib.1"])
hist(RNASeqCounts_NC8[,"NC8.rib.2"])

# Create boxplot for Counts
# Note that these are plots before normalisation, for data exploration and quality assesment.
boxplot(RNASeqCounts_WCFS1, col="gray", las=3, main="RNASeqCounts_WCFS1")
boxplot(RNASeqCounts_NC8, col="gray", las=3, main="RNASeqCounts_NC8")

# Create MA plot for Counts

# M values
## WCFS1.glc.1 vs WCFS1.glc.2
M_values_WCFS1_glc <- (RNASeqCounts_WCFS1[,1] - RNASeqCounts_WCFS1[,2])
## WCFS1.rib.1 vs WCFS1.rib.2
M_values_WCFS1_rib <- (RNASeqCounts_WCFS1[,3] - RNASeqCounts_WCFS1[,4])
## NC8.glc.1 vs NC8.glc.2
M_values_NC8_glc <- (RNASeqCounts_NC8[,1] - RNASeqCounts_NC8[,2])
## NC8.rib.1 vs NC8.rib.2
M_values_NC8_rib <- (RNASeqCounts_NC8[,3] - RNASeqCounts_NC8[,4])

# A values
## WCFS1.glc.1 vs WCFS1.glc.2
A_values_WCFS1_glc <- (RNASeqCounts_WCFS1[,1] + RNASeqCounts_WCFS1[,2])/2
## WCFS1.rib.1 vs WCFS1.rib.2
A_values_WCFS1_rib <- (RNASeqCounts_WCFS1[,3] + RNASeqCounts_WCFS1[,4])/2
## NC8.glc.1 vs NC8.glc.2
A_values_NC8_glc <- (RNASeqCounts_NC8[,1] + RNASeqCounts_NC8[,2])/2
## NC8.rib.1 vs NC8.rib.2
A_values_NC8_rib <- (RNASeqCounts_NC8[,3] + RNASeqCounts_NC8[,4])/2

# plot MA
# Note that these are plots before normalisation, for data exploration and quality assesment.
plot(A_values_WCFS1_glc, M_values_WCFS1_glc, xlab="A", ylab="M", pch=19, main="WCFS1_glc")
abline(h=0, col="red")
plot(A_values_NC8_glc, M_values_NC8_glc, xlab="A", ylab="M", pch=19, main="NC8_glc")
abline(h=0, col="red")
plot(A_values_WCFS1_rib, M_values_WCFS1_rib, xlab="A", ylab="M", pch=19, main="WCFS1_rib")
abline(h=0, col="red")
plot(A_values_NC8_rib, M_values_NC8_rib, xlab="A", ylab="M", pch=19, main="NC8_rib")
abline(h=0, col="red")

# MDS for Counts (using limma package) # soort PCA
# Note that these are plots before normalisation, for data exploration and quality assesment.
par(mfrow=c(1,2))
plotMDS(RNASeqCounts_WCFS1, main="RNASeqCounts_WCFS1")
plotMDS(RNASeqCounts_NC8, main="RNASeqCounts_NC8")

# heatmap for Counts (using mixOmics package)
# Note that these are plots before normalisation, for data exploration and quality assesment.
sampleDists_WCFS1 <- as.matrix(dist(t(RNASeqCounts_WCFS1)))
sampleDists_NC8 <- as.matrix(dist(t(RNASeqCounts_NC8)))
Color <- colorRampPalette(c("Yellow", "orange")) #set color
heatmap.2(sampleDists_WCFS1, col=Color, main="RNASeqCounts WCFS1",trace="none", margins = c(10,12), cexRow=0.5, cexCol=0.5)
heatmap.2(sampleDists_NC8, col=Color, main="RNASeqCounts NC8",trace="none", margins = c(10,12), cexRow=0.5, cexCol=0.5)

#################Differential expression analysis########################
# remove genes with zero counts for all samples
dgeFull_WCFS1 <- DGEList(y_WCFS1$counts[apply(y_WCFS1$counts, 1, sum) != 0, ],
                   group=y_WCFS1$samples$group)
dgeFull_NC8 <- DGEList(y_NC8$counts[apply(y_NC8$counts, 1, sum) != 0, ],
                   group=y_NC8$samples$group)

# filter data
keep.genes_WCFS1 <- rowSums(cpm(y_WCFS1)>50) >= 2
keep.genes_NC8 <- rowSums(cpm(y_NC8)>50) >= 2
y_NC8 <- y_NC8[keep.genes_NC8,]
y_WCFS1 <- y_WCFS1[keep.genes_WCFS1,]
y_NC8$samples$lib.size <- colSums(y_NC8$counts)
y_WCFS1$samples$lib.size <- colSums(y_WCFS1$counts)

# estimate the normalization factors
y_WCFS1 <- calcNormFactors(dgeFull_WCFS1, method="TMM") # 2:9
y_NC8 <- calcNormFactors(dgeFull_NC8, method="TMM") # 2:9
y_WCFS1$samples
y_NC8$samples
head(y_WCFS1$counts)
head(y_NC8$counts)

# From the normalization factors and the original count table, find the normalized counts and use the log2-transformation to inspect them with boxplots and a MDS. 
# Normalized counts can be extracted from dgeFull using the function cpm:
eff.lib.size <- y_WCFS1$samples$lib.size*y_WCFS1$samples$norm.factors
normCounts_WCFS1 <- cpm(dgeFull_WCFS1)
eff.lib.size <- y_NC8$samples$lib.size*y_NC8$samples$norm.factors
normCounts_NC8 <- cpm(dgeFull_NC8)

# Plot
RNASeqNormCounts_WCFS1 <- log2(normCounts_WCFS1 + 1)
RNASeqNormCounts_NC8 <- log2(normCounts_NC8 + 1)
boxplot(RNASeqNormCounts_WCFS1, col="gray", las=3, main="RNASeqNormCounts_WCFS1 (log2)")
boxplot(RNASeqCounts_WCFS1, col="gray", las=3, main="RNASeqNormCounts_NC8 (log2)")
plotMDS(RNASeqNormCounts_WCFS1, main="RNASeqNormCounts_WCFS1 (log2)")
plotMDS(RNASeqNormCounts_NC8, main="RNASeqNormCounts_NC8 (log2)")

sampleDists_WCFS1 <- as.matrix(dist(t(RNASeqNormCounts_WCFS1)))
sampleDists_NC8 <- as.matrix(dist(t(RNASeqNormCounts_NC8)))
Color <- colorRampPalette(c("Yellow", "orange")) #set color
heatmap.2(sampleDists_WCFS1, col=Color, main="sampleDists_WCFS1",trace="none", margins = c(10,12), cexRow=0.5, cexCol=0.5)
heatmap.2(sampleDists_NC8, col=Color, main="sampleDists_NC8",trace="none", margins = c(10,12), cexRow=0.5, cexCol=0.5)

# Estimate common and tagwise dispersion
dgeFull_WCFS1 <- estimateCommonDisp(dgeFull_WCFS1)
dgeFull_WCFS1 <- estimateTagwiseDisp(dgeFull_WCFS1)
dgeTest_WCFS1 <- exactTest(dgeFull_WCFS1)
filtData_WCFS1 <- HTSFilter(dgeFull_WCFS1)$filteredData

# Perform an exact test for the difference in expression between the conditions
dgeTestFilt_WCFS1 <- exactTest(filtData_WCFS1)
dgeTestFilt_WCFS1

# plot an histogram of unadjusted p-values after filtering
hist(dgeTest_WCFS1$table[,"PValue"], breaks=50, main="dgeTest_WCFS1")
hist(dgeTestFilt_WCFS1$table[,"PValue"], breaks=50, main="dgeTestFilt_WCFS1")

# extract a summary of the differential expression statistics
resNoFilt <- topTags(dgeTest_WCFS1, n=nrow(dgeTest_WCFS1$table))
head(resNoFilt$table)
resFilt <- topTags(dgeTestFilt_WCFS1, n=nrow(dgeTest_WCFS1$table))
head(resFilt$table)

# compare the number of differentially expressed genes with and without filtering (risk: 1%)
# before independent filtering
sum(resNoFilt$table$FDR < 0.01)
# after independent filtering
sum(resFilt$table$FDR < 0.01)

# extract and sort differentially expressed genes
sigDownReg <- resFilt$table[resFilt$table$FDR<0.01,]
sigDownReg <- sigDownReg[order(sigDownReg$logFC),]

sigUpReg <- sigDownReg[order(sigDownReg$logFC, decreasing=TRUE),]
head(sigUpReg)

# write the results in csv files
write.csv(sigDownReg, file="sigDownReg_WCFS1.csv")
write.csv(sigUpReg, file="sigUpReg_WCFS1.csv")

# create a MA plot with 1% differentially expressed genes
plotSmear(dgeTestFilt_WCFS1,
          de.tags = rownames(resFilt$table)[which(resFilt$table$FDR<0.01)], main="dgeTestFilt_WCFS1 1% diff exp genes")

# create a Volcano plot
volcanoData <- cbind(resFilt$table$logFC, -log10(resFilt$table$FDR))
colnames(volcanoData) <- c("logFC", "negLogPval")
head(volcanoData)
plot(volcanoData, pch=19, main="Volcanoplot WCFS1")

# transform the normalized counts in log-counts-per-million
y <- cpm(dgeFull_WCFS1, log=TRUE, prior.count = 1)
head(y)

# select 1% differentially expressed genes and produce a heatmap
selY <- y[rownames(resFilt$table)[resFilt$table$FDR<0.01 & 
                                    abs(resFilt$table$logFC)>1.5],]
head(selY)

# If you are interesting in the result of the gene clustering, it can be obtained from the previous command. 
# More precisely, the result of HAC is stored into $ddc.
cimColor <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)[255:1]
finalHM <- cim(t(selY), color=cimColor, symkey=FALSE)

plot(finalHM$ddc, leaflab="none", main="Heatmap WCFS1 (1% diff expressed genes)")
abline(h=10, lwd=2, col="pink")
# Heatmap results in a lot of locus tags, less data or a larger figure solves this.


# Using this dendrogram, we might want to cut the tree at level h=10 (for instance), 
# which can be performed using the function cutree, which will provide a cluster membership for each gene.
geneClust <- cutree(as.hclust(finalHM$ddc), h=10)
head(geneClust)

# For instance, the number of clusters is equal to
length(unique(geneClust))

# and the genes in cluster 1 are:
names(which(geneClust==1))
