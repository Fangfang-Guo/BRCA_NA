rm(list = ls())
set.seed(7)

library(dplyr)
library(readr)
library(stringr)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)
library(GenomicFeatures)
library(clusterProfiler)
library(biomaRt)
library(org.Hs.eg.db)
library(ReactomePA)
library(readxl)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v75)

# -----------------------------
# Step 1: Convert VCF-derived CSV to BED
# -----------------------------
vcf_dir <- "data/vcf_input"
output_bed_dir <- "data/bed_files"
dir.create(output_bed_dir, showWarnings = FALSE)

vcf_files <- list.files(vcf_dir, pattern = "\\.freebayes\\.vcf$", full.names = TRUE)

for (vcf_file in vcf_files) {
  file_name <- basename(vcf_file)
  file_prefix <- tools::file_path_sans_ext(file_name)

  # Example: assume a CNV CSV is generated elsewhere, here just a placeholder
  csv_file <- paste0(vcf_dir, "/", file_prefix, "_CNV.csv")
  if (!file.exists(csv_file)) {
    message("CNV CSV not found for: ", file_prefix)
    next
  }

  csv_data <- read.csv(csv_file)
  if (!all(c("chrom", "start", "stop", "amplification") %in% colnames(csv_data))) {
    message("Required columns not found in: ", csv_file)
    next
  }

  bed_data <- csv_data %>%
    dplyr::select(chrom, start, stop, amplification) %>%
    mutate(start = start - 1)

  bed_file <- file.path(output_bed_dir, paste0(file_prefix, ".bed"))
  write.table(bed_data, file = bed_file,
              quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
  message("Generated BED: ", bed_file)
}

# -----------------------------
# Step 2: Annotate BED with gene information
# -----------------------------
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
folder_path <- "data/bed_files"
output_filter_folder <- "data/bed_files_annotated"
dir.create(output_filter_folder, showWarnings = FALSE)

panel_genes <- read_excel("data/panel_gene_list.xlsx")
panel_gene_list <- panel_genes$Gene

bed_files <- list.files(folder_path, pattern = "\\.bed$", full.names = TRUE)

for (bed_file in bed_files) {
  peaks.bed <- read.table(bed_file, sep = "\t", header = FALSE, blank.lines.skip = TRUE)
  peaks.bed[, 1] <- paste0("chr", peaks.bed[, 1])

  peaks.gr <- GRanges(seqnames = peaks.bed[, 1],
                      ranges = IRanges(peaks.bed[, 2], peaks.bed[, 3]),
                      strand = "*",
                      status = peaks.bed[, 4])

  overlaps <- findOverlaps(peaks.gr, genes(txdb, single.strand.genes.only = FALSE))

  segment_gene_map <- data.frame(
    peak = queryHits(overlaps),
    gene_id = subjectHits(overlaps),
    status = peaks.gr$status[queryHits(overlaps)]
  )

  segment_gene_map$SYMBOL <- mapIds(org.Hs.eg.db,
                                    keys = as.character(segment_gene_map$gene_id),
                                    column = "SYMBOL",
                                    keytype = "ENTREZID",
                                    multiVals = "first")

  filtered_genes <- segment_gene_map %>%
    dplyr::filter(SYMBOL %in% panel_gene_list)

  colnames(filtered_genes)[colnames(filtered_genes) == "status"] <- "amplification"
  colnames(filtered_genes)[colnames(filtered_genes) == "SYMBOL"] <- "gene_name"

  output_file <- file.path(output_filter_folder,
                           paste0(tools::file_path_sans_ext(basename(bed_file)),
                                  "_filtered_gene_segments.tsv"))

  write.table(filtered_genes, file = output_file,
              sep = "\t", row.names = FALSE, quote = FALSE)
  message("Annotated and filtered BED saved: ", output_file)
}

message("All files processed. Results in: ", output_filter_folder)
