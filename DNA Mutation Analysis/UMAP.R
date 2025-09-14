# Load required packages
library(maftools)
library(ggplot2)
library(ggrepel)
library(umap)


# Load MAF files
maf_na <- read.maf(maf = "filtered_NA_final.maf")
maf_ca <- read.maf(maf = "filtered_CA_final.maf")

# Create binary mutation matrices (rows: genes, columns: samples)
matrix_na <- mutCountMatrix(maf_na)
matrix_ca <- mutCountMatrix(maf_ca)

# Merge and preprocess
mutation_matrix <- merge(as.data.frame(matrix_na), as.data.frame(matrix_ca), by = "row.names", all = TRUE)
rownames(mutation_matrix) <- mutation_matrix$Row.names
mutation_matrix <- mutation_matrix[, -1]
mutation_matrix[mutation_matrix > 0] <- 1
mutation_matrix[is.na(mutation_matrix)] <- 0

# Construct sample metadata
sample_metadata <- data.frame(
  Sample = colnames(mutation_matrix),
  Race = ifelse(grepl("TCGA", colnames(mutation_matrix)), "Caucasian", "Native American")
)

# UMAP before filtering
umap_matrix <- t(mutation_matrix)
umap_result <- umap(umap_matrix, n_neighbors = 15, min_dist = 0.1, metric = "cosine")
umap_data <- as.data.frame(umap_result)
colnames(umap_data) <- c("UMAP1", "UMAP2")
umap_data$Sample <- rownames(umap_matrix)
umap_data <- merge(umap_data, sample_metadata, by = "Sample")

# UMAP plot
ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = Race)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(title = "UMAP of Breast Cancer Mutations", x = "UMAP1", y = "UMAP2") +
  scale_color_manual(values = c("Native American" = "blue", "Caucasian" = "red"))

