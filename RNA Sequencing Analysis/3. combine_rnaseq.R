library(data.table)

# Load expression matrix for Native American breast cancer patients
na_expr <- fread("data/na_rna_expression_matrix.csv")
na_expr <- as.data.frame(na_expr)

# Extract sample columns (assumes first column is gene symbol)
na_sample_ids <- setdiff(colnames(na_expr), "hgnc_symbol")
na_breast <- na_expr[, na_sample_ids, drop = FALSE]
rownames(na_breast) <- na_expr$hgnc_symbol

# Load expression matrix for TCGA Caucasian breast cancer patients
tcga_breast <- readRDS("data/BRCA_Caucasian_RNASeq.rds")

# Ensure common gene set
common_genes <- intersect(rownames(na_breast), rownames(tcga_breast))
na_breast_common <- na_breast[common_genes, ]
tcga_breast_common <- tcga_breast[common_genes, ]
stopifnot(all(rownames(na_breast_common) == rownames(tcga_breast_common)))

# Combine Native and TCGA breast cancer samples
combined_breast <- cbind(na_breast_common, tcga_breast_common)

# Save merged expression matrix (rows = genes, columns = samples)
write.csv(combined_breast, file = "combined_breast_rnaseq.csv", row.names = TRUE)
