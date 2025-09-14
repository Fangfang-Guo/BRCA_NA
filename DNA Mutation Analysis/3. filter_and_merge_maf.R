library(readxl)
library(dplyr)
library(maftools)
library(openxlsx)

# Read gene list
gene_list <- read_excel("Tempus_Genes_List.xlsx")[[1]]

# Load NA cohort MAF files
maf_files_na <- list.files("path/to/NA_maf", pattern = "\\.maf$", full.names = TRUE)
maf_na_list <- lapply(maf_files_na, read.maf)
combined_na <- merge_mafs(maf_na_list)
combined_na@clinical.data$Race <- "Native American"
filtered_na <- subsetMaf(combined_na, gene = gene_list)

# Load CA (TCGA) cohort MAF files
maf_files_ca <- list.files("path/to/TCGA_maf", pattern = "\\.maf\\.gz$", full.names = TRUE, recursive = TRUE)
maf_ca_list <- list()

for (f in maf_files_ca) {
  temp <- tryCatch(read.maf(f, vc_nonSyn = NULL), error = function(e) NULL)
  if (!is.null(temp)) {
    non_syn <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site",
                 "Translation_Start_Site", "Nonsense_Mutation", "Nonstop_Mutation",
                 "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation")
    if (sum(temp@data$Variant_Classification %in% non_syn) > 0) {
      maf_ca_list[[f]] <- temp
    }
  }
}
combined_ca <- merge_mafs(maf_ca_list)
combined_ca@clinical.data$Race <- "Caucasian"

# Keep one MAF per sample
barcodes <- as.character(combined_ca@variants.per.sample$Tumor_Sample_Barcode)
sample_prefix <- sapply(strsplit(barcodes, "-"), function(x) paste(x[1:3], collapse = "-"))
samples_to_keep <- barcodes[!duplicated(sample_prefix)]
filtered_ca <- subsetMaf(combined_ca, tsb = samples_to_keep)
filtered_ca <- subsetMaf(filtered_ca, gene = gene_list)

# Define filtering function
filter_query_gnomAD <- function(prefix) {
  paste0("(", paste0("is.na(", prefix, "_AF) | ", prefix, "_AF < 0.0005", collapse = ") & ("), ")")
}

# Apply frequency filtering
filtered_na_final <- subsetMaf(filtered_na, query = filter_query_gnomAD("gnomADe"))
filtered_ca_final <- subsetMaf(filtered_ca, query = filter_query_gnomAD("gnomAD"))

# Merge MAFs and plot oncoplot
final_maf <- merge_mafs(list(filtered_na_final, filtered_ca_final))
oncoplot(final_maf, clinicalFeatures = "Race", sortByAnnotation = TRUE, draw_titv = TRUE, top = 20)

# Save results
write.table(filtered_na_final@data, "filtered_NA_final.maf", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(filtered_ca_final@data, "filtered_CA_final.maf", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(final_maf@data, "final_combined_maf.maf", sep = "\t", quote = FALSE, row.names = FALSE)
