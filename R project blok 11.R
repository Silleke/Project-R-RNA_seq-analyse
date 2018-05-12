# install.packages("rstudioapi") # run this if it's your first time using it to install

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

set_wd()

# Get files
# data <- read.delim(file.choose(), header=T)
WCFS1_anno <- read.delim("WCFS1_anno.txt", header=T)
RNA_Seq_counts <- read.delim("RNA_Seq_counts.txt", header=T)

merged_data = merge(RNA_Seq_counts, WCFS1_anno, by.x="ID", by.y="ORF")
