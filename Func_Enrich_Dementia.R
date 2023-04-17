
# Set path to input files
path <- "C:/Users/Desktop/MSAM-Northeastern/MATH7343/Project/Data/"

# Importing data
res <- read.csv(paste0(path, "differential_expression_dementia.csv"), header = TRUE)
countData <- read.csv(paste0(path, "fpkm_table_normalized.csv"), header = TRUE)
row_genes <- read.csv(paste0(path,"rows-genes.csv"))

# rank genes based on log-fold change
rankedGenes <- res[order(res$log2FoldChange, decreasing = TRUE),]

# extract upregulated genes
upGenes <- subset(rankedGenes, log2FoldChange > 1)
nrow(upGenes)

# extract downregulated genes
downGenes <- subset(rankedGenes, log2FoldChange < -1)
nrow(downGenes)

# calculate adjusted p-value (FDR) for each gene
resSig <- res[res$pvalue < 0.05,]
resSig

# filter DEG list to retain only significant genes with LFC or FC cutoff
sigGenes <- resSig[abs(resSig$log2FoldChange) > 1,]
sigGenes
# generate heatmap of significant DEGs
sig_count <- countData[row.names(sigGenes),]
heatmap(as.matrix(sig_count), scale = "col")

# generate volcano plot of significant DEGs
with(resSig, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot of Significant DEGs", xlim=c(-3,3)))
with(subset(resSig, padj < 0.05), points(log2FoldChange, -log10(padj), col="red", pch=20))


