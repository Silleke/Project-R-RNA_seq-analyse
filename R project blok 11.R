# run this if it's your first time using it to install
# install.packages(c("RColorBrewer", "mixOmics", "VennDiagram", "gplots", "stringi", "rstudioapi"))
source("https://bioconductor.org/biocLite.R") # install bioconductor if not already installed
library(limma)
library(edgeR) # for DGEList
library(RColorBrewer)
library(MASS)
library(lattice)
library(ggplot2)
library(mixOmics)
library(gplots)
biocLite("HTSFilter") 
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
WCFS1_anno <- read.delim("WCFS1_anno.txt", header=T)
RNA_Seq_counts <- read.delim("RNA_Seq_counts.txt", header=T)

# Have a look at the count data:
# head(RNA_Seq_counts)
# nrow(RNA_Seq_counts)

# merge the data
merged_data = merge(RNA_Seq_counts, WCFS1_anno, by.x="ID", by.y="ORF")

# remove NA values
merged_data <- merged_data[1:(length(merged_data)-4)]

#############################Starting from count table##############################
# create DGEList
exp_WCFS1 <- c("WCFS1.glc.1","WCFS1.glc.2","WCFS1.rib.1","WCFS1.rib.2")
exp_NC8 <- c("NC8.glc.1", "NC8.glc.2", "NC8.rib.1", "NC8.rib.2")
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
# Create histogram for Counts 
hist(RNASeqCounts_WCFS1[,"WCFS1.glc.1"])
hist(RNASeqCounts_WCFS1[,"WCFS1.glc.2"])
hist(RNASeqCounts_WCFS1[,"WCFS1.rib.1"])
hist(RNASeqCounts_WCFS1[,"WCFS1.rib.2"])
hist(RNASeqCounts_NC8[,"NC8.glc.1"])
hist(RNASeqCounts_NC8[,"NC8.glc.2"])
hist(RNASeqCounts_NC8[,"NC8.rib.1"])
hist(RNASeqCounts_NC8[,"NC8.rib.2"])

# Create boxplot for Counts
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
plot(A_values_WCFS1_glc, M_values_WCFS1_glc, xlab="A", ylab="M", pch=19, main="WCFS1_glc treated")
abline(h=0, col="red")
plot(A_values_NC8_glc, M_values_NC8_glc, xlab="A", ylab="M", pch=19, main="NC8_glc treated")
abline(h=0, col="red")
plot(A_values_WCFS1_rib, M_values_WCFS1_rib, xlab="A", ylab="M", pch=19, main="WCFS1_rib treated")
abline(h=0, col="red")
plot(A_values_NC8_rib, M_values_NC8_rib, xlab="A", ylab="M", pch=19, main="NC8_rib treated")
abline(h=0, col="red")

# MDS for Counts (using limma package) # soort PCA
par(mfrow=c(1,2))
plotMDS(RNASeqCounts_WCFS1, main="RNASeqCounts_WCFS1")
plotMDS(RNASeqCounts_NC8, main="RNASeqCounts_NC8")

# heatmap for Counts (using mixOmics package)
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

# estimate the normalization factors
y_WCFS1 <- calcNormFactors(dgeFull_WCFS1, method="TMM") # 2:9
y_NC8 <- calcNormFactors(dgeFull_NC8, method="TMM") # 2:9
y_WCFS1$samples
y_NC8$samples
head(y_WCFS1$counts)
head(y_NC8$counts)

eff.lib.size <- y_WCFS1$samples$lib.size*y_WCFS1$samples$norm.factors
normCounts_WCFS1 <- cpm(dgeFull_WCFS1)
eff.lib.size <- y_NC8$samples$lib.size*y_NC8$samples$norm.factors
normCounts_NC8 <- cpm(dgeFull_NC8)

RNASeqNormCounts_WCFS1 <- log2(normCounts_WCFS1 + 1)
RNASeqNormCounts_NC8 <- log2(normCounts_NC8 + 1)
boxplot(RNASeqCounts_WCFS1, col="gray", las=3)
boxplot(RNASeqCounts_NC8, col="gray", las=3)
plotMDS(RNASeqNormCounts_WCFS1)
plotMDS(RNASeqNormCounts_NC8)

#make PCA plot
PCA <- prcomp(normCounts, retx=F, center = TRUE, scale = TRUE)
res <- as.data.frame(PCA$x)
res$label - rownames(res)
ggplot(res, aes(PC1, PC2)) + geom_point() + scale_color_manual(values-c("red", "black"))
# text --> text toevoegen aan pca plot