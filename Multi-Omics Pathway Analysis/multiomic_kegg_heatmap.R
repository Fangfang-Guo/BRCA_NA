# ===============================
# Libraries
# ===============================
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(dplyr)
library(readr)
library(pheatmap)
library(grid)

# ===============================
# Function: KEGG enrichment for a gene set
# ===============================
get_kegg_padj_and_genes <- function(gene_symbols) {
  gene_df <- tryCatch({
    bitr(gene_symbols,
         fromType = "SYMBOL",
         toType = "ENTREZID",
         OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    warning("Gene symbol conversion failed.")
    return(NULL)
  })
  
  if (is.null(gene_df) || nrow(gene_df) == 0) {
    empty_df <- data.frame(Description = NA, p.adjust = NA, geneID = NA)
    return(list(padj_result = empty_df[, 1:2], geneid_result = empty_df[, c(1,3)]))
  }
  
  entrez_ids <- unique(na.omit(gene_df$ENTREZID))
  m_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
  
  ego <- enricher(entrez_ids,
                  TERM2GENE = m_df[, c("gs_name", "entrez_gene")],
                  pvalueCutoff = 1,
                  pAdjustMethod = "BH")
  
  ego <- setReadable(ego, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  result <- as.data.frame(ego)[, c("Description", "p.adjust", "geneID")]
  
  all_paths <- unique(m_df$gs_name)
  
  padj_result <- data.frame(Description = all_paths) %>%
    left_join(result[, c("Description", "p.adjust")], by = "Description")
  
  geneid_result <- data.frame(Description = all_paths) %>%
    left_join(result[, c("Description", "geneID")], by = "Description")
  
  return(list(padj_result = padj_result, geneid_result = geneid_result))
}

# ===============================
# Load input data
# ===============================
deg_info <- read.csv("DEG_INFO_breast.csv", header = TRUE)
rownames(deg_info) <- deg_info$X

# Genes upregulated in each ancestry group
deg_up_na_genes <- deg_info[deg_info$Upregulated_group == "Native American", "X"]
deg_up_ca_genes <- deg_info[deg_info$Upregulated_group == "Caucasian", "X"]

# Significant CNV genes (Gain/Loss)
CNV_Gain_Info <- read_csv("Significant_CNV_Gain_Genes.csv")
cnv_gain_na_genes <- CNV_Gain_Info[CNV_Gain_Info$Gain_Higher_In == "Native American", "Gene"]
cnv_gain_ca_genes <- CNV_Gain_Info[CNV_Gain_Info$Gain_Higher_In == "Caucasian", "Gene"]

CNV_Loss_Info <- read_csv("Significant_CNV_Loss_Genes.csv")
cnv_loss_na_genes <- CNV_Loss_Info[CNV_Loss_Info$Loss_Higher_In == "Native American", "Gene"]
cnv_loss_ca_genes <- CNV_Loss_Info[CNV_Loss_Info$Loss_Higher_In == "Caucasian", "Gene"]

# ===============================
# Define gene sets from multiple omics
# ===============================
gene_sets <- list(
  "Mutation_NA"   = c("ARID1B", "BCOR", "DNMT3A", "ERCC5", "FANCL", 
                      "FOXO1", "HLA-DRB1", "HLA-DRB5", "INPP4B", "NOTCH4", "POLE"),
  "DEG_Up_NA"     = deg_up_na_genes,
  "DEG_Up_CA"     = deg_up_ca_genes,
  "CNV_Gain_NA"   = cnv_gain_na_genes$Gene,
  "CNV_Loss_NA"   = cnv_loss_na_genes$Gene,
  "CNV_Gain_CA"   = cnv_gain_ca_genes$Gene,
  "CNV_Loss_CA"   = cnv_loss_ca_genes$Gene
)

# ===============================
# Run KEGG enrichment for each gene set
# ===============================
all_kegg_results <- lapply(gene_sets, get_kegg_padj_and_genes)
all_kegg_padj_results <- lapply(all_kegg_results, function(res) res$padj_result)
all_kegg_geneid_results <- lapply(all_kegg_results, function(res) res$geneid_result)

# Build padj matrix (pathway × gene set)
named_padj_list <- lapply(all_kegg_padj_results, function(df) {
  setNames(df$p.adjust, df$Description)
})

padj_df <- bind_rows(named_padj_list, .id = "GeneSet") %>%
  tibble::column_to_rownames("GeneSet") %>%
  t() %>%
  as.data.frame()

padj_mat <- apply(padj_df, 2, as.numeric)
rownames(padj_mat) <- rownames(padj_df)

# ===============================
# Prepare heatmap data
# ===============================
padj_log <- -log10(padj_mat)
padj_log_clean <- padj_log[rowSums(is.na(padj_log)) < ncol(padj_log),
                           colSums(is.na(padj_log)) < nrow(padj_log)]

# Mark significant pathways (padj < 0.05)
star_mat <- ifelse(padj_mat < 0.05, "*", NA)
star_mat_clean <- star_mat[rownames(padj_log_clean), colnames(padj_log_clean)]
star_mat_clean[is.na(star_mat_clean)] <- ""

# Sort rows by number of significant hits
sig_count <- rowSums(padj_mat < 0.05, na.rm = TRUE)
padj_log_clean <- padj_log_clean[order(-sig_count[rownames(padj_log_clean)]), ]
star_mat_clean <- star_mat_clean[rownames(padj_log_clean), ]

# Remove KEGG_ prefix for display
rownames(padj_log_clean) <- gsub("^KEGG_", "", rownames(padj_log_clean))
rownames(star_mat_clean) <- rownames(padj_log_clean)

# ===============================
# Plot heatmap
# ===============================
p <- pheatmap(
  padj_log_clean,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  angle_col = 0,
  na_col = "white",
  main = "Pathway Enrichment across Multi-omic Alterations",
  display_numbers = star_mat_clean,
  number_color = "black",
  fontsize_row = 8,
  fontsize_col = 7,
  legend = TRUE
)

grid.text(expression("-log"[10] * "(padj)"),
          x = 0.96, y = 0.98, gp = gpar(fontsize = 9))

# ===============================
# Extract significant genes per pathway
# ===============================
named_geneid_list <- lapply(all_kegg_geneid_results, function(df) {
  setNames(df$geneID, df$Description)
})

geneid_df <- bind_rows(named_geneid_list, .id = "GeneSet") %>%
  tibble::column_to_rownames("GeneSet") %>%
  t() %>%
  as.data.frame()

geneid_clean <- geneid_df[rowSums(is.na(geneid_df)) < ncol(geneid_df),
                          colSums(is.na(geneid_df)) < nrow(geneid_df)]
rownames(geneid_clean) <- gsub("^KEGG_", "", rownames(geneid_clean))
geneid_clean <- geneid_clean[rownames(star_mat_clean), , drop = FALSE]

# Mask non-significant entries
sig_geneid_clean <- geneid_clean
sig_geneid_clean[as.data.frame(star_mat_clean) == ""] <- NA

# Save if needed
# write.csv(sig_geneid_clean, "sig_pathway_genesubset_table.csv", row.names = TRUE)
