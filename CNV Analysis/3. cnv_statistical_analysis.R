library(dplyr)
library(readxl)

# -------------------------------------------------------------------
# Step 1: Load gene list
# Expecting an Excel file with gene names in the first column
# -------------------------------------------------------------------
gene_list <- read_excel("Tempus_Genes_List.xlsx")[[1]]
gene_list <- unique(gene_list)

# Initialize empty CNV matrix (genes x samples)
cnv_data <- matrix(
  NA, nrow = length(gene_list), ncol = 0,
  dimnames = list(gene_list, NULL)
)

# -------------------------------------------------------------------
# Step 2: Define input folders
# Each folder should contain per-sample CNV .tsv files
# -------------------------------------------------------------------
na_folder <- "data/NA_samples"         # Native American cohort
ca_folder <- "data/Caucasian_samples"  # Caucasian cohort

na_files <- list.files(na_folder, pattern = "\\.tsv$", full.names = TRUE)
ca_files <- list.files(ca_folder, pattern = "\\.tsv$", full.names = TRUE)

# -------------------------------------------------------------------
# Step 3: Function to extract CNV calls per gene
# -------------------------------------------------------------------
fill_cnv_data <- function(file, sample_name, cnv_data) {
  tsv_data <- read.delim(file, header = TRUE, stringsAsFactors = FALSE)

  if (!all(c("gene_name", "amplification") %in% colnames(tsv_data))) {
    warning(paste("File", file, "is missing required columns. Skipping..."))
    return(cnv_data)
  }

  matched_values <- tsv_data[match(rownames(cnv_data), tsv_data$gene_name), "amplification"]
  cnv_data <- cbind(cnv_data, matched_values)
  colnames(cnv_data)[ncol(cnv_data)] <- sample_name
  return(cnv_data)
}

# -------------------------------------------------------------------
# Step 4: Add Native American samples
# -------------------------------------------------------------------
for (i in seq_along(na_files)) {
  sample_name <- paste0("NA_", sprintf("%03d", i))
  cnv_data <- fill_cnv_data(na_files[i], sample_name, cnv_data)
}

# -------------------------------------------------------------------
# Step 5: Add Caucasian samples
# -------------------------------------------------------------------
for (i in seq_along(ca_files)) {
  sample_name <- paste0("CA_", sprintf("%03d", i))
  cnv_data <- fill_cnv_data(ca_files[i], sample_name, cnv_data)
}

# Replace blanks with NA
cnv_data[cnv_data == ""] <- NA
cnv_data <- as.data.frame(cnv_data)

# Save the final CNV matrix
write.table(
  cnv_data,
  file = "data/final_cnv_matrix.tsv",
  sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA
)

# -------------------------------------------------------------------
# Step 6: Fisher's exact test per gene
# -------------------------------------------------------------------
results <- data.frame(
  Gene = rownames(cnv_data),
  Gain_P_Value = NA,
  Gain_Odds_Ratio = NA,
  Loss_P_Value = NA,
  Loss_Odds_Ratio = NA,
  Change_P_Value = NA,
  Change_Odds_Ratio = NA,
  Gain_Padj = NA,
  Loss_Padj = NA,
  Change_Padj = NA,
  Gain_Higher_In = NA,
  Loss_Higher_In = NA,
  CNV_Change_Higher_In = NA
)

native_samples <- grep("^NA_", colnames(cnv_data))
caucasian_samples <- grep("^CA_", colnames(cnv_data))

for (i in 1:nrow(cnv_data)) {
  gene_data <- cnv_data[i, ]

  na_total <- sum(!is.na(gene_data[native_samples]))
  ca_total <- sum(!is.na(gene_data[caucasian_samples]))
  if (na_total == 0 || ca_total == 0) next

  # Gain
  na_gain <- sum(gene_data[native_samples] == "gain", na.rm = TRUE)
  ca_gain <- sum(gene_data[caucasian_samples] == "gain", na.rm = TRUE)
  gain_table <- matrix(c(na_gain, na_total - na_gain, ca_gain, ca_total - ca_gain), nrow = 2)
  if (all(gain_table >= 0)) {
    gain_test <- fisher.test(gain_table)
    results$Gain_P_Value[i] <- gain_test$p.value
    results$Gain_Odds_Ratio[i] <- gain_test$estimate
  }

  # Loss
  na_loss <- sum(gene_data[native_samples] == "loss", na.rm = TRUE)
  ca_loss <- sum(gene_data[caucasian_samples] == "loss", na.rm = TRUE)
  loss_table <- matrix(c(na_loss, na_total - na_loss, ca_loss, ca_total - ca_loss), nrow = 2)
  if (all(loss_table >= 0)) {
    loss_test <- fisher.test(loss_table)
    results$Loss_P_Value[i] <- loss_test$p.value
    results$Loss_Odds_Ratio[i] <- loss_test$estimate
  }

  # CNV Change
  na_change <- sum(gene_data[native_samples] %in% c("gain", "loss"), na.rm = TRUE)
  ca_change <- sum(gene_data[caucasian_samples] %in% c("gain", "loss"), na.rm = TRUE)
  change_table <- matrix(c(na_change, na_total - na_change, ca_change, ca_total - ca_change), nrow = 2)
  if (all(change_table >= 0)) {
    change_test <- fisher.test(change_table)
    results$Change_P_Value[i] <- change_test$p.value
    results$Change_Odds_Ratio[i] <- change_test$estimate
  }

  # Frequencies
  na_gain_freq <- na_gain / na_total
  ca_gain_freq <- ca_gain / ca_total
  na_loss_freq <- na_loss / na_total
  ca_loss_freq <- ca_loss / ca_total
  na_change_freq <- na_change / na_total
  ca_change_freq <- ca_change / ca_total

  results$Gain_Higher_In[i] <- ifelse(!is.na(results$Gain_P_Value[i]) && results$Gain_P_Value[i] < 0.05,
                                      ifelse(na_gain_freq > ca_gain_freq, "Native American", "Caucasian"), NA)
  results$Loss_Higher_In[i] <- ifelse(!is.na(results$Loss_P_Value[i]) && results$Loss_P_Value[i] < 0.05,
                                      ifelse(na_loss_freq > ca_loss_freq, "Native American", "Caucasian"), NA)
  results$CNV_Change_Higher_In[i] <- ifelse(!is.na(results$Change_P_Value[i]) && results$Change_P_Value[i] < 0.05,
                                            ifelse(na_change_freq > ca_change_freq, "Native American", "Caucasian"), NA)
}

# Multiple testing correction
results$Gain_Padj <- p.adjust(results$Gain_P_Value, method = "BH")
results$Loss_Padj <- p.adjust(results$Loss_P_Value, method = "BH")
results$Change_Padj <- p.adjust(results$Change_P_Value, method = "BH")

# Save results
write.table(results, file = "results/cnv_fisher_results.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
