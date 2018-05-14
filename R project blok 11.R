# install.packages("rstudioapi") # run this if it's your first time using it to install
# install.packages("stringi")
# install.packages(c("RColorBrewer", "mixOmics", "VennDiagram"))
# install.packages("gplots")
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
#############create a function set_wd that sets the working directory############## 
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

# execute set_wd
set_wd()

################################data#######################################
#Importing the files
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
exp <- c("WCFS1.glc","WCFS1.glc","WCFS1.rib","WCFS1.rib")
group <- factor(exp)
y <- DGEList(counts=RNA_Seq_counts[,2:5],group=group)


#######################Data exploration and quality assessment######################
# Extract RNA_Seq_Counts
RNA_Seq_Counts <- log2(y$counts+1)
head(RNA_Seq_Counts)

# Create histogram for pseudo-counts 
hist(RNA_Seq_Counts[,"WCFS1.glc.1"])
hist(RNA_Seq_Counts[,"WCFS1.glc.2"])
hist(RNA_Seq_Counts[,"WCFS1.rib.1"])
hist(RNA_Seq_Counts[,"WCFS1.rib.2"])

# Create boxplot for pseudo-counts
boxplot(RNA_Seq_Counts, col="gray", las=3)

# Create MA plot for pseudo-counts
par(mfrow=c(1,2))

## WT1 vs WT2
# M values
M_values <- (RNA_Seq_Counts[,1] - RNA_Seq_Counts[,2])

# A values
A_values <- (RNA_Seq_Counts[,1] + RNA_Seq_Counts[,2])/2

plot(A_values, M_values, xlab="A", ylab="M", pch=19, main="treated")
abline(h=0, col="red")

## Mt1 vs Mt2
# M values
# Nog even kijken of de waardes van M_values / A_values kloppen
M_values <- (RNA_Seq_Counts[,2] - RNA_Seq_Counts[,3])
A_values <- (RNA_Seq_Counts[,2] + RNA_Seq_Counts[,3])/2
plot(A_values, M_values, xlab="A", ylab="M", pch=19, main="control")
abline(h=0, col="red")

# MDS for pseudo-counts (using limma package)
plotMDS(RNA_Seq_Counts)

# heatmap for pseudo-counts (using mixOmics package)
sampleDists <- as.matrix(dist(t(RNA_Seq_Counts)))
Color <- colorRampPalette(c("Yellow", "orange")) #set color
heatmap.2(sampleDists, col=Color, main="Sample",trace="none", margins = c(10,12), cexRow=0.5, cexCol=0.5)

# werkt nog niet nog naar kijken
#cimColor <- colorRampPalette(c(brewer.pal(9, "Blues")))(16)
#par(mar=rep(2,4))
#heatmap <- cim(sampleDists, color=cimColor, symkey=FALSE)
#nba.m <- melt(nba)
#nba.m <- ddply(nba.m, .(variable), transform,
#               +     rescale = rescale(value))



#################Differential expression analysis########################
# remove genes with zero counts for all samples
dgeFull <- DGEList(y$counts[apply(y$counts, 1, sum) != 0, ],
                   group=y$samples$group)

# estimate the normalization factors
y <- calcNormFactors(dgeFull, method="TMM")
y$samples
head(y$counts)

eff.lib.size <- y$samples$lib.size*y$samples$norm.factors
normCounts <- cpm(dgeFull)
RNA_Seq_NormCounts <- log2(normCounts + 1)
boxplot(RNA_Seq_NormCounts, col="gray", las=3)
plotMDS(RNA_Seq_NormCounts)

# make PCA plot
# PCA <- prcomp(x.norm, retx=F, center = TRUE, scale = TRUE)
# text --> text toevoegen aan pca plot