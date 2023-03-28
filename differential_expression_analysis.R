
# setting the correct working directory
setwd('D:/Classes/Applied Statistics/Project')

# DEseq2 will be used to perform the differential expression analysis

# Install the latest version of DEseq2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2", version = "3.16")


# loading the libraries
library(DESeq2)

# Reading the dataset
rawCounts <- read.csv("./dataset/kaggle_dataset/fpkm_table_unnormalized.csv")
head(rawCounts)

# Reading the sample data
sampleData <- read.csv("./dataset/kaggle_dataset/columns-samples.csv")
head(sampleData)









