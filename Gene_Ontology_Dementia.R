# Loading packages
library(clusterProfiler)

# Set path to input files
path <- "C:/Users/Desktop/MSAM-Northeastern/MATH7343/Project/Data/"

# Importing data
DE_Results <- read.table(paste0(path, "differential_expression_dementia.csv"), sep = ",", header = TRUE, row.names = 1)
row_genes <- read.csv(paste0(path,"rows-genes.csv"))

DEResults_id <- merge(DE_Results, row_genes, by = "gene_id")
DEResults_id

# Extract gene symbols from DEseq2 results
DEgenes_symbol <- DEResults_id[["gene_symbol"]]
DEgenes_symbol

# Enrichment analysis
go_enrichment <- enrichGO(gene = DEgenes_symbol,
                          OrgDb = org.Hs.eg.db,
                          keyType = "SYMBOL",
                          ont = "BP",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)

go_enrichment
dotplot(go_enrichment, showCategory = 15)
