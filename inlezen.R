library(edgeR)

#Laden van de data
data = read.delim("data/RNA-Seq-counts.txt", header = TRUE, sep = "\t", comment.char = "#")
annotatie = read.delim("data/WCFS1_anno.txt", header = TRUE, sep = "\t", comment.char = "#")

#Mergen van data
combine = cbind(data, annotatie[rownames(data),])
combined_data = merge(data, annotatie, by.x="ID", by.y="ORF")
row.names(combined_data) <- combined_data[,"ID"]

#Verkrijgen benodigde data
needed_data <- c("WCFS1.glc","WCFS1.glc","WCFS1.rib","WCFS1.rib")
group_data <- factor(needed_data)
new_data <- DGEList(counts=combined_data[,2:5],group=group_data)

#Filteren data
keep.genes <- rowSums(cpm(new_data)>50) >= 2
new_data<- new_data[keep.genes,]
new_data$samples$lib.size <- colSums(new_data$counts)

#Normaliseren data
new_data<- calcNormFactors(new_data, method="TMM" )

#Design matrix
design <- model.matrix(~0+group_data, data=new_data$samples)
colnames(design) <- levels(new_data$samples$group_data)
print(design)

#Estimate dispersion
new_data <- estimateGLMCommonDisp(new_data,design)
new_data <- estimateGLMTrendedDisp(new_data,design, method="power")
new_data <- estimateGLMTagwiseDisp(new_data,design)

#Plot van genormaliseerde data
pdf("Results.pdf")
plotMDS(new_data)
plotBCV(new_data)
dev.off()

#Bepalen differential expressed genes
fit <- glmFit(new_data,design)
mc <- makeContrasts(needed_data=WCFS1.glc-WCFS1.rib, levels=design)
fit <- glmLRT(fit, contrast=mc)
res <- topTags(fit)
print(res)