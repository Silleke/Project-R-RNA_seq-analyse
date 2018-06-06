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


WCFS1_anno <- read.delim(file="C:/Users/xx_xx/Downloads/WCFS1_anno1.txt", header=TRUE, sep="\t")
RNA <- read.delim(file="C:/Users/xx_xx/Downloads/RNA-Seq-counts.txt", header=TRUE, sep="\t")



merged_data <- merge(RNA, WCFS1_anno, by.x=("ID"), by.y=("ORF"))

merged_data <- merged_data[1:(length(merged_data)-4)]


#counts
fDir <-  "C:/Users/xx_xx/Downloads/"
fName <- "WCFS1_anno.txt"

cnts <- read.delim(paste0(fDir,fName))
row.names(merged_data) <- merged_data[,"ID"]

#DGE list
exp_WCFS1 <- c("WCFS1.glc","WCFS1.glc","WCFS1.rib","WCFS1.rib")
exp_NC8 <- c("NC8.glc", "NC8.glc", "NC8.rib", "NC8.rib")
group_WCFS1 <- factor(exp_WCFS1)
group_NC8 <-factor(exp_NC8)

y_WCFS1 <- DGEList(counts=merged_data[,2:5],group=group_WCFS1)
y_NC8 <- DGEList(counts=merged_data[,6:9],group=group_NC8)

#filter data
keep.genes_WCFS1 <- rowSums(cpm(y_WCFS1)>50) >= 2
keep.genes_NC8 <- rowSums(cpm(y_NC8)>50) >= 2
y_NC8 <- y_NC8[keep.genes_NC8,]
y_WCFS1 <- y_WCFS1[keep.genes_WCFS1,]
y_NC8$samples$lib.size <- colSums(y_NC8$counts)
y_WCFS1$samples$lib.size <- colSums(y_WCFS1$counts)

#normaliseer
y_NC8 <- calcNormFactors(y_NC8, method="TMM" )
y_WCFS1 <- calcNormFactors(y_WCFS1, method="TMM" )

#check statistieken
print("Count statistics NC8")
print(summary(y_NC8$counts))
print(y_NC8$samples)

print("Count statistics WCFS1")
print(summary(y_WCFS1$counts))
print(y_WCFS1$samples)

#creeer design matrix
design_NC8 <- model.matrix(~0+group, data=y_NC8$samples)
colnames(design_NC8) <- levels(y_NC8$samples$group)
print(design_NC8)

design_WCFS1 <- model.matrix(~0+group, data=y_WCFS1$samples)
colnames(design_WCFS1) <- levels(y_WCFS1$samples$group)
print(design_WCFS1)

#dispersie bepalen
y_NC8 <- estimateGLMCommonDisp(y_NC8,design_NC8)
y_NC8 <- estimateGLMTrendedDisp(y_NC8,design_NC8, method="power")
y_NC8 <- estimateGLMTagwiseDisp(y_NC8,design_NC8)

y_WCFS1 <- estimateGLMCommonDisp(y_WCFS1,design_WCFS1)
y_WCFS1 <- estimateGLMTrendedDisp(y_WCFS1,design_WCFS1, method="power")
y_WCFS1 <- estimateGLMTagwiseDisp(y_WCFS1,design_WCFS1)

#resultaten plotten
pdf(paste0(fDir,"LP_edgeR.pdf"))
plotMDS(y_NC8)
plotBCV(y_NC8)
plotMDS(y_WCFS1)
plotBCV(y_WCFS1)
dev.off()

#data fitten
fit_NC8 <- glmFit(y_NC8,design_NC8)
fit_WCFS1 <- glmFit(y_WCFS1,design_WCFS1)

#fold changes

mc_NC8 <- makeContrasts(exp.r=NC8.glc-NC8.rib, levels=design_NC8)
mc_WCFS1 <- makeContrasts(exp.r=WCFS1.glc-WCFS1.rib, levels=design_WCFS1)

fit_NC8 <- glmLRT(fit_NC8, contrast=mc_NC8)
fit_WCFS1 <- glmLRT(fit_WCFS1, contrast=mc_WCFS1)

#print top Tags
res_NC8<-topTags(fit_NC8)
print(res_NC8)
res_WCFS1<-topTags(fit_WCFS1)
print(res_WCFS1)

write.csv(res_NC8, file = "res_NC8.csv")
read.csv("res_NC8.csv", row.names = 1)
write.csv(res_WCFS1, file = "res_WCFS1.csv")
read.csv("res_WCFS1.csv", row.names = 1)


write.csv(fit_NC8[["table"]], file = "all_NC8.csv")
#read.csv("all_NC8.csv", row.names = 1)
write.csv(fit_WCFS1[["table"]], file = "all_WCFS1.csv")
#read.csv("all_WCFS1.csv", row.names = 1)


#fold changes rib-glc

mc_NC8 <- makeContrasts(exp.r=NC8.rib-NC8.glc, levels=design_NC8)
mc_WCFS1 <- makeContrasts(exp.r=WCFS1.rib-WCFS1.glc, levels=design_WCFS1)

fit_NC8 <- glmLRT(fit_NC8, contrast=mc_NC8)
fit_WCFS1 <- glmLRT(fit_WCFS1, contrast=mc_WCFS1)

#print top Tags
res_NC8<-topTags(fit_NC8)
print(res_NC8)
res_WCFS1<-topTags(fit_WCFS1)
print(res_WCFS1)

write.csv(res_NC8, file = "res_NC8rib_glc.csv")
read.csv("res_NC8.csv", row.names = 1)
write.csv(res_WCFS1, file = "res_WCFS1rib_glc.csv")
read.csv("res_WCFS1.csv", row.names = 1)


write.csv(fit_NC8[["table"]], file = "all_NC8rib_glc.csv")
#read.csv("all_NC8rib-glc.csv", row.names = 1)
write.csv(fit_WCFS1[["table"]], file = "all_WCFS1rib_glc.csv")
#read.csv("all_WCFS1rib-glc.csv", row.names = 1)