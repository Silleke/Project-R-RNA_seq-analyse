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
library(edgeR)


exp <- c("WCFS1.glc","WCFS1.glc","WCFS1.rib","WCFS1.rib")

group <- factor(exp)
y <- DGEList(counts=merged_data[,2:5],group=group)
#filter data
keep.genes <- rowSums(cpm(y)>50) >= 2
y <- y[keep.genes,]

y$samples$lib.size <- colSums(y$counts)

#normaliseer
y <- calcNormFactors(y, method="TMM" )

#check statistieken
print("Count statistics")
print(summary(y$counts))
print(y$samples)

#creeer design matrix
design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(y$samples$group)
print(design)

#dispersie bepalen
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design, method="power")
y <- estimateGLMTagwiseDisp(y,design)

#resultaten plotten
pdf(paste0(fDir,"LP_edgeR.pdf"))
plotMDS(y)
plotBCV(y)
dev.off()

#data fitten
fit <- glmFit(y,design)

#fold changes

mc <- makeContrasts(exp.r=WCFS1.glc-WCFS1.rib, levels=design)

fit <- glmLRT(fit, contrast=mc)

#print top Tags
res<-topTags(fit)
print(res)
