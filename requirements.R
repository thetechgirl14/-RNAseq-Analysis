# Run this script to install all the necessary packages required for this project.

# Install required packages
install.packages(c("BiocManager", "AnnotationDbi"))

# Install Bioconductor packages
BiocManager::install(c("DESeq2","clusterProfiler","org.Hs.eg.db"))

# Install additional packages for data visualization and analysis
install.packages(c("ggplot2", "dplyr", "viridis","tidyverse", "reshape2"))
