library(dplyr)
library(readxl)

# Define input/output directories
base_dir <- "data/gene_level_cnv"
output_dir <- "results/cnv_filtered_by_panel_genes"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load gene panel list
panel_gene_list <- read_excel("data/panel_gene_list.xlsx")[[1]]

# Identify all .tsv files in input directory and subfolders
all_files <- list.files(base_dir, pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE)

# Process each CNV file
for (file in all_files) {
  tryCatch({
    df <- read.delim(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

    if ("gene_name" %in% colnames(df)) {
      filtered_df <- df %>% filter(gene_name %in% panel_gene_list)

      if ("copy_number" %in% colnames(filtered_df)) {
        filtered_df <- filtered_df %>%
          mutate(amplification = ifelse(copy_number > 2.2, "gain",
                                 ifelse(copy_number < 1.8, "loss", "neutral")))
      } else {
        warning(paste("Missing column 'copy_number' in:", file))
      }

      # Save filtered result
      filtered_file_name <- paste0(tools::file_path_sans_ext(basename(file)), "_filtered.tsv")
      filtered_file_path <- file.path(output_dir, filtered_file_name)
      write.table(filtered_df, file = filtered_file_path, sep = "\t", row.names = FALSE, quote = FALSE)
      message("Filtered file saved: ", filtered_file_path)

    } else {
      warning(paste("Missing column 'gene_name' in:", file))
    }
  }, error = function(e) {
    message("Error processing ", file, ": ", e$message)
  })
}
