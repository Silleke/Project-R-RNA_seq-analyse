# Todt (c) HAN
#load file
fDir <- "C:/Users/xx_xx/Downloads/"

cnts <- read.delim(paste0(fDir,"RNA-Seq-counts.txt"))
#prints cnts
print(head(cnts))

cnts <- cnts[,-1]
#prints cnts without ID
m <- as.matrix(cnts)
print(head(m))  
#prints sum of m per column
print(head(colSums(m)))      

#
lm <- log2(m+1)

lm <- na.omit(lm)
m <- na.omit(m)
#make pdf file with boxplot and cluster results
pdf(paste0(fDir,"resPCA.pdf"),width=10,height=10)
 
boxplot(lm)

plot(hclust(dist(t(m), method="eu")))
#PCA plot
pr <- prcomp(m, centre=TRUE, scale.=TRUE)
plot(pr)
summary(pr)
biplot(pr)
dev.off()
