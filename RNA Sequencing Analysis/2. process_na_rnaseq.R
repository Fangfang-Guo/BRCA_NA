library(data.table)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v75)
library(readr)
library(dplyr)

# Load RNA expression data from three batches
rna1 <- read_csv("data/RNA_Expression1/Data/Group_Level_Molecular/normalized_rna.csv")
rna2 <- read_csv("data/RNA_Expression2/Data/Group_Level_Molecular/normalized_rna.csv")
rna3 <- read_csv("data/RNA_Expression3/Data/Group_Level_Molecular/normalized_rna.csv")

# Combine data into a single data.table
data <- rbindlist(list(rna1, rna2, rna3))
setDT(data)

# Extract required columns
selected_data <- data[, .(partner_patient_id, ensembl_gene, gene_raw)]

# Get unique Ensembl gene IDs
unique_genes <- unique(selected_data$ensembl_gene)

# Map Ensembl gene IDs to HGNC gene symbols using GRCh37 (v75)
org_symbols <- mapIds(EnsDb.Hsapiens.v75,
                      keys = unique_genes,
                      column = "SYMBOL",
                      keytype = "GENEID",
                      multiVals = "first")

# Convert mapping to data.frame
org_symbols_df <- data.frame(ensembl_gene_id = names(org_symbols),
                             hgnc_symbol = org_symbols,
                             stringsAsFactors = FALSE)

# Check unmapped gene symbols
missing_symbols <- org_symbols_df$ensembl_gene_id[is.na(org_symbols_df$hgnc_symbol)]
cat("Number of unmapped gene symbols:", length(missing_symbols), "\n")

# Check duplicated gene symbols
symbol_counts <- table(org_symbols_df$hgnc_symbol)
duplicate_symbols <- symbol_counts[symbol_counts > 1]
cat("Number of duplicated gene symbols (after merge):", length(duplicate_symbols), "\n")
print("Duplicated symbols and counts:")
print(duplicate_symbols)

# Merge annotation with original data
merged_data <- merge(selected_data, org_symbols_df,
                     by.x = "ensembl_gene", by.y = "ensembl_gene_id", all.x = TRUE)

# Aggregate expression values by sample and gene symbol
summed_data <- merged_data[, .(gene_expression_sum = sum(gene_raw, na.rm = TRUE)),
                           by = .(partner_patient_id, hgnc_symbol)]

# Cast long-format to wide-format expression matrix
expression_matrix <- dcast(summed_data, hgnc_symbol ~ partner_patient_id,
                           value.var = "gene_expression_sum")

# Save final expression matrix
write.csv(expression_matrix,
          file = "na_rna_expression_matrix.csv",
          row.names = FALSE)
