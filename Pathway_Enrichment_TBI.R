# Loading packages
library(clusterProfiler)

# Set path to input files
path <- "C:/Users/Desktop/MSAM-Northeastern/MATH7343/Project/Data/"

# Importing data
DE_Results <- read.table(paste0(path, "differential_expression_tbi.csv"), sep = ",", header = TRUE, row.names = 1)
row_genes <- read.csv(paste0(path,"rows-genes.csv"))

DEResults_id <- merge(DE_Results, row_genes, by = "gene_id")
DEResults_id

# Extract gene Ids from DEseq2 results
DEgenes_id <- DEResults_id[["gene_entrez_id"]]

kegg_enrichment <- enrichKEGG(gene = DEgenes_id,
                              organism = 'hsa',
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.05)


kegg_enrichment
barplot(kegg_enrichment, showCategory = 10)


