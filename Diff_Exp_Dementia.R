# Loading Packages
library(DESeq2)
library(ggplot2)

# Set path to input files
path <- "C:/Users/Desktop/MSAM-Northeastern/MATH7343/Project/Data/"

# Importing data
counts <- read.table(paste0(path, "fpkm_table_normalized.csv"), sep = ",", header = TRUE, row.names = 1)
coldata <- read.csv(paste0(path, "info-samples.csv"), header=TRUE, row.names=1)

# Filter out low-expressed genes
keep <- rowSums(counts) >= 10
counts_filt <- counts[keep,]

# Creating DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(counts_filt), colData = coldata, design =~ act_demented)

# Normalizing data
dds <- DESeq(dds)

# Obtaining differential expression results
res <- results(dds)
res <- res[order(res$pvalue), ]

# Plot data relations to define filtering criteria
ggplot(data = as.data.frame(res), aes(x = pvalue)) + 
  geom_histogram(bins = 100)

# Filtering significant genes
sig_genes <- subset(res, pvalue < 0.05) 
sig_genes_df <- as.data.frame(sig_genes)
sig_genes_with_id <- rownames_to_column(sig_genes_df, var = "gene_id")

# Save Diffentially Expressed Genes for further analysis
write.csv(sig_genes_with_id, paste0(path,"differential_expression_dementia.csv"))