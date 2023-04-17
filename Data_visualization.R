# Load required packages
library(ggplot2)
library(dplyr)
library(tidyverse)
library(reshape2)

# Set path to input files
path <- "C:/Users/Desktop/MSAM-Northeastern/MATH7343/Project/Data/"

# Load FPKM normalized data
fpkm_table_normalized <- read.csv(paste0(path,"fpkm_table_normalized.csv"), header = TRUE)
colnames(fpkm_table_normalized)[1] <- "gene_id"

# Load sample information data
sample_info <- read.csv(paste0(path,"info-samples.csv"), header = TRUE)

# transpose fpkm matrix to convert sample names to row names
fpkm_table_normalized_t <- t(fpkm_table_normalized)

# remove the additional X in the row names
rownames(fpkm_table_normalized_t) <- substring(rownames(fpkm_table_normalized_t), 2)

#Shift gene_id to column indexing
colnames(fpkm_table_normalized_t) <- head(fpkm_table_normalized_t,1)

#rename first column variable
fpkm_table_normalized_t <- rownames_to_column(as.data.frame(fpkm_table_normalized_t), var = "rnaseq_profile_id")

#remove extra first row
fpkm_table_normalized_t <- fpkm_table_normalized_t[-1,]
fpkm_table_normalized_t <- data.frame(fpkm_table_normalized_t)

# Merge sample information data with FPKM data
merged_df <- merge(sample_info, fpkm_table_normalized_t, by = "rnaseq_profile_id")

# print the merged dataframe
print(merged_df)
merged_df[1:4,1:4]
colnames(merged_df[1:35])

# Subset data for Dementia Presence and Subtypes
da_dsm_data <- merged_df %>%
  filter(dsm_iv_clinical_diagnosis != "No_Dementia") %>%
  select(starts_with("X"), age, dsm_iv_clinical_diagnosis)


# Calculate average expression values for each gene by diagnosis group
da_dsm_mean <- da_dsm_data %>%
  group_by(dsm_iv_clinical_diagnosis) %>%
  summarize(across(starts_with("X"), mean, na.rm = TRUE))


# Reshape data for ggplot2
da_dsm_mean_melted <- melt(da_dsm_mean, id.vars = "dsm_iv_clinical_diagnosis", variable.name = "Gene", value.name = "Average FPKM")


# Create boxplot of average gene expression by diagnosis group
ggplot(data = da_dsm_mean_melted, aes(x = dsm_iv_clinical_diagnosis, y = `Average FPKM`)) +
  geom_boxplot() +
  ggtitle("Average Gene Expression by Diagnosis Group") +
  ylab("Average FPKM") +
  xlab("Classification found in the DSM of Mental Disorders") +
  theme_bw()

# Plot gene expression data for Dementia samples
dementia_data <- merged_df %>% filter(act_demented == "Dementia") %>% select(starts_with("X"), age, sex, act_demented)

# Calculate average expression values for each gene by dementia
dementia_mean <- dementia_data %>%
  group_by(act_demented) %>%
  summarize(across(starts_with("X"), mean, na.rm = TRUE))


# Reshape data for ggplot2
dementia_mean_melted <- melt(dementia_mean, id.vars = "act_demented", variable.name = "Gene", value.name = "Average FPKM")


ggplot(data = dementia_mean_melted, aes(x = act_demented, y = `Average FPKM`, color = Gene)) + 
  geom_point() + 
  ggtitle("Average Gene Expression by Dementia") +
  ylab("log2(FPKM + 1)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Plot gene expression data for Aging samples
aging_data <- merged_df %>% filter(age >= 65) %>% select(starts_with("X"),gene_id, age)

# Calculate average expression values for each gene by age
aging_mean <- aging_data %>%
  group_by(age) %>%
  summarize(across(starts_with("X"), mean, na.rm = TRUE))


# Reshape data for ggplot2
aging_mean_melted <- melt(aging_mean, id.vars = "age", variable.name = "Gene", value.name = "Average FPKM")


ggplot(data = aging_mean_melted, aes(x = age, y = `Average FPKM`, color = gene_id)) + 
  geom_point() +  
  ggtitle("Average Gene Expression by Age") +
  ylab("log2(FPKM + 1)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Plot gene expression data for TBI samples
tbi_data <- merged_df %>% filter(num_tbi_w_loc >= 1) %>% select(starts_with("X"),gene_id, num_tbi_w_loc)

# Calculate average expression values for each gene by TBI
tbi_mean <- tbi_data %>%
  group_by(num_tbi_w_loc) %>%
  summarize(across(starts_with("X"), mean, na.rm = TRUE))


# Reshape data for ggplot2
tbi_mean_melted <- melt(tbi_mean, id.vars = "num_tbi_w_loc", variable.name = "Gene", value.name = "Average FPKM")


ggplot(data = tbi_mean_melted, aes(x = num_tbi_w_loc, y = `Average FPKM`, color = gene_id)) + 
  geom_point() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Plot correlation between FPKM values and age
aging_data <- merged_df %>% filter(!is.na(age)) %>% select(starts_with("X"), age)

# Calculate average expression values for each gene by age
aging_mean <- aging_data %>%
  group_by(age) %>%
  summarize(across(starts_with("X"), mean, na.rm = TRUE))


# Reshape data for ggplot2
aging_mean_melted <- melt(aging_mean, id.vars = "age", variable.name = "Gene", value.name = "Average FPKM")


ggplot(data = aging_mean_melted, aes(x = age, y = `Average FPKM`)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)

# Plot correlation between FPKM values and CERAD score
cerad_data <- merged_df %>% filter(!is.na(cerad)) %>% select(starts_with("X"), cerad)

# Calculate average expression values for each gene by CERAD score
cerad_mean <- cerad_data %>%
  group_by(cerad) %>%
  summarize(across(starts_with("X"), mean, na.rm = TRUE))


# Reshape data for ggplot2
cerad_mean_melted <- melt(cerad_mean, id.vars = "cerad", variable.name = "Gene", value.name = "Average FPKM")


ggplot(data = cerad_mean_melted, aes(x = cerad, y = `Average FPKM`)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)




