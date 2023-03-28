
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

# saving the sample data for later use
sampleData_v2 <- sampleData 

# Convert count data to a matrix of appropriate form that DEseq2 can read
geneID <- rawCounts$gene_id...rnaseq_profile_id
sampleIndex <- grepl("X\\d+", colnames(rawCounts))
rawCounts <- as.matrix(rawCounts[,sampleIndex])
rownames(rawCounts) <- geneID
head(rawCounts)

# importing more information for each sample
info_samples <- read.csv("./dataset/kaggle_dataset/info-samples.csv")

# joining info-samples and sampledata
sampleDataFull = merge(x = sampleData, y = info_samples, by = "rnaseq_profile_id")

# Convert sample variable mappings to an appropriate form that DESeq2 can read
rownames(sampleDataFull) <- sampleDataFull$rnaseq_profile_id
keep <- c("sex", "act_demented")
sampleDataFull <- sampleDataFull[,keep]
colnames(sampleDataFull) <- c("sex", "dementia")
sampleDataFull$sex <- factor(sampleDataFull$sex)
head(sampleDataFull)

# remove the additional X in the coumn names in row counts
colnames(rawCounts) <- substring(colnames(rawCounts), 2)

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
rawCounts <- rawCounts[,unique(rownames(sampleDataFull))]
all(colnames(rawCounts) == rownames(sampleDataFull))











