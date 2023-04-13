library(clusterProfiler)
library(org.Hs.eg.db)

# Set path to input files
path <- "C:/Users/abhil/OneDrive/Desktop/git clone/MSAM-Northeastern/MATH7343/Project/Data/"

row_genes <- read.csv(paste0(path,"rows-genes.csv"))

# Load DE results data
DEresults <- read.csv(paste0(path,"differential_expression_results.csv"), header = TRUE)
DEresults

merged_df <- merge(DEresults, row_genes, by = "gene_id")
merged_df

# Extract gene names from DEseq2 results
DEgenes_symbol <- merged_df[["gene_symbol"]]
DEgenes_symbol

DEgenes_id <- merged_df[["gene_entrez_id"]]
DEgenes_id

# Perform Gene Ontology analysis
go_enrichment <- enrichGO(gene = DEgenes_symbol,
                          OrgDb = org.Hs.eg.db, # Specify the organism database
                          keyType = "SYMBOL", # Specify the gene ID type
                          ont = "BP", # Specify the ontology (e.g., biological process, molecular function, cellular component)
                          pvalueCutoff = 0.05, # Specify the p-value cutoff for significance
                          qvalueCutoff = 0.05) # Specify the q-value cutoff for significance

# Print the GO enrichment results
print(go_enrichment)

# Perform KEGG pathway analysis
kegg_enrichment <- enrichKEGG(gene = DEgenes_id,
                               organism = "hsa", # Specify the organism
                               pvalueCutoff = 0.05, # Specify the p-value cutoff for significance
                               qvalueCutoff = 0.0 5) # Specify the q-value cutoff for significance

# Print the KEGG 
print(kegg_enrichment)

# Create a bar plot of GO enrichment results
barplot(go_enrichment, showCategory = 10)

# Create a bar plot of KEGG results
barplot(kegg_enrichment, showCategory = 10)