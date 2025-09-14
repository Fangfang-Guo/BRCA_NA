# Load required packages
library(TCGAbiolinks)
library(EnsDb.Hsapiens.v75)
library(AnnotationDbi)
library(xml2)
library(SummarizedExperiment)

# Define project
project <- "TCGA-BRCA"

# Query RNA-seq count data
query <- GDCquery(
  project = project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# Download and prepare RNA-seq data
GDCdownload(query)
data <- GDCprepare(query)

# Query clinical data to filter by race
query_clinical <- GDCquery(
  project = project,
  data.category = "Clinical",
  data.type = "Clinical Supplement"
)

GDCdownload(query_clinical)

# Locate XML files
clinical_path <- file.path("GDCdata", project, "Clinical", "Clinical_Supplement")
xml_files <- list.files(clinical_path, pattern = "\\.xml$", full.names = TRUE, recursive = TRUE)

# Parse clinical XML files for submitter ID and race
clinical_info <- data.frame(submitter_id = character(), race = character(), stringsAsFactors = FALSE)
for (file in xml_files) {
  doc <- read_xml(file)
  ns <- xml_ns(doc)
  sid <- xml_text(xml_find_first(doc, ".//shared:bcr_patient_barcode", ns))
  race <- xml_text(xml_find_first(doc, ".//clin_shared:race", ns))
  if (!is.na(sid) && sid != "") {
    clinical_info <- rbind(clinical_info, data.frame(submitter_id = sid, race = race))
  }
}

# Select only White (Caucasian) patients
caucasian_barcodes <- clinical_info$submitter_id[tolower(clinical_info$race) == "white"]

# Extract expression matrix
counts_data <- assay(data)
sample_barcodes <- substr(colnames(counts_data), 1, 12)
caucasian_data <- counts_data[, sample_barcodes %in% caucasian_barcodes, drop = FALSE]
colnames(caucasian_data) <- sample_barcodes[sample_barcodes %in% caucasian_barcodes]

# Retain only one sample per patient
unique_samples <- !duplicated(colnames(caucasian_data))
caucasian_data <- caucasian_data[, unique_samples]
caucasian_data <- as.data.frame(caucasian_data)

# Remove Ensembl version suffix from gene IDs
unique_genes <- gsub("\\.[0-9]+_PAR_Y$", "", rownames(caucasian_data))
unique_genes <- gsub("\\.[0-9]+$", "", unique_genes)
caucasian_data$unique_genes <- unique_genes

# Aggregate counts by unique gene ID
caucasian_data_aggregated <- aggregate(. ~ unique_genes, data = caucasian_data, FUN = mean)
rownames(caucasian_data_aggregated) <- caucasian_data_aggregated$unique_genes
caucasian_data_aggregated <- caucasian_data_aggregated[, -which(names(caucasian_data_aggregated) == "unique_genes")]

# Map Ensembl gene IDs to gene symbols
caucasian_data_aggregated$ensembl_gene_id <- rownames(caucasian_data_aggregated)
org_symbols <- mapIds(EnsDb.Hsapiens.v75,
                      keys = caucasian_data_aggregated$ensembl_gene_id,
                      column = "SYMBOL",
                      keytype = "GENEID",
                      multiVals = "first")
caucasian_data_aggregated$symbol <- org_symbols[rownames(caucasian_data_aggregated)]

# Remove genes without mapped gene symbols
caucasian_data_aggregated <- caucasian_data_aggregated[!is.na(caucasian_data_aggregated$symbol), ]
caucasian_data_aggregated <- caucasian_data_aggregated[, !(names(caucasian_data_aggregated) == "ensembl_gene_id")]

# Collapse rows with duplicate gene symbols by averaging
numeric_columns <- sapply(caucasian_data_aggregated, is.numeric)
caucasian_data_aggregated <- aggregate(caucasian_data_aggregated[, numeric_columns],
                                       by = list(symbol = caucasian_data_aggregated$symbol),
                                       FUN = mean)

# Set gene symbols as rownames
rownames(caucasian_data_aggregated) <- caucasian_data_aggregated$symbol
caucasian_data_aggregated <- caucasian_data_aggregated[, !(names(caucasian_data_aggregated) == "symbol")]

# Save final expression matrix
saveRDS(caucasian_data_aggregated, file = "BRCA_Caucasian_RNASeq.rds")
