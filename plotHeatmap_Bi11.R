# Todt (c) HAN

library(gplots)
library(RColorBrewer)

fDir <- "C:/Users/xx_xx/Downloads/"

### read counts or normalized CPM values
cnts <- read.table(paste0(fDir,"CPM_glc_rib.txt"), header=TRUE, sep="\t", quote="", fill=F,comment.char="", stringsAsFactor=F)
cnts <- fit[["AveLogCPM"]]
### annotation
genes <- read.table(paste0(fDir,"LacplantDB_v11_may2010.txt"), header=TRUE, sep="\t", quote="", fill=F, stringsAsFactor=F)
genes <- genes[,c("ORF","name","function.","class","subclass")]

### fold changes
fc <- read.table(paste0(fDir,"FC_glc_rib.txt"), header=TRUE, sep="\t", quote="", fill=F,comment.char="", stringsAsFactor=F)

sig <- merge(cnts, genes, all.x=T, by.x="ID",by.y="ORF")
sig <- merge(sig, fc, by.x="ID",by.y="ID")

print(head(sig))

#########################################################################

pdf(paste0(fDir,"resheatmap.pdf"),width=10,height=10)

sig <- sig[,c("WCFS1.glc.1","WCFS1.glc.2","WCFS1.rib.1","WCFS1.rib.2","WCFS1.glc.over.WCFS1.rib","fdr.WCFS1.glc.over.WCFS1.rib","class","name","function.")]
sig <- sig[ sig[,"fdr.WCFS1.glc.over.WCFS1.rib"]<0.05, ]
sig <- sig[ abs(sig[,"WCFS1.glc.over.WCFS1.rib"])>4, ]
print(nrow(sig))
m <- as.matrix(sig[,c(1,2,3,4)])
colnames(m) <- c("WCFS1.glc.1","WCFS1.glc.2","WCFS1.rib.1","WCFS1.rib.2")
rownames(m) <- paste(sig[,"name"],sig[,"function."])
        
print(head(m))        

my_palette <- colorRampPalette(c("blue","green"))(n = 200)
  
heatmap.2(m, 
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,32),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier 
  dendrogram="row",     # only draw a row dendrogram
  Colv="NA")            # turn off column clustering
	
rownames(m) <- paste(sig[,"name"],sig[,"class"])
plot(hclust(dist(m, method="eu")))

dev.off()
