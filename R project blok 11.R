# install.packages("rstudioapi") # run this if it's your first time using it to install
# source("https://bioconductor.org/biocLite.R") install bioconductor if not already installed
# library(edgeR) # for DGEList
# library(limma)
# library(RColorBrewer)
# library(mixOmics)
biocLite("HTSFilter") # 
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

################################preprocessing data#######################################
#Importing the files
# data <- read.delim(file.choose(), header=T)
WCFS1_anno <- read.delim("WCFS1_anno.txt", header=T)
RNA_Seq_counts <- read.delim("RNA_Seq_counts.txt", header=T)

# Have a look at the count data:
# head(RNA_Seq_counts) 
# nrow(RNA_Seq_counts)  

# merge the data
merged_data = merge(merged_data_Seq_counts, WCFS1_anno, by.x="ID", by.y="ORF")

# remove NA values
merged_data <- merged_data[1:(length(merged_data)-4)]

#############################Starting from count table##############################
# create DGEList
exp <- c("WCFS1.glc","WCFS1.glc","WCFS1.rib","WCFS1.rib")
y <- DGEList(counts=RNA_Seq_counts[,2:5],group=group)


#######################Data exploration and quality assessment######################
# Extract pseudo-counts
pseudoCounts <- log2(y$counts+1)
head(pseudoCounts)

# Create histogram for pseudo-counts 
hist(pseudoCounts[,"WCFS1.glc.1"])
hist(pseudoCounts[,"WCFS1.glc.2"])
hist(pseudoCounts[,"WCFS1.rib.1"])
hist(pseudoCounts[,"WCFS1.rib.2"])

# Create boxplot for pseudo-counts
boxplot(pseudoCounts, col="gray", las=3)

# Create MA plot for pseudo-counts
par(mfrow=c(1,2))

## WT1 vs WT2
# M values
M_values <- (pseudoCounts[,1] - pseudoCounts[,2])

# A values
A_values <- (pseudoCounts[,1] + pseudoCounts[,2])/2

plot(A_values, M_values, xlab="A", ylab="M", pch=19, main="treated")
abline(h=0, col="red")

## Mt1 vs Mt2
# M values
# Nog even kijken of de waardes van M_values / A_values kloppen
M_values <- (pseudoCounts[,2] - pseudoCounts[,3])
A_values <- (pseudoCounts[,2] + pseudoCounts[,3])/2
plot(A_values, M_values, xlab="A", ylab="M", pch=19, main="control")
abline(h=0, col="red")

