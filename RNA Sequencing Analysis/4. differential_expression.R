library(data.table)
library(EnhancedVolcano)

# Load combined RNA-seq count matrix (genes x samples)
counts <- read.csv("data/combined_breast_rnaseq.csv", row.names = 1)

# Define group sizes
n1 <- 17  # Native American samples (columns 1 to 17)
n2 <- ncol(counts) - n1  # Caucasian samples (columns 18+)

# Filter genes with sufficient expression in both groups
C <- 10
mean_na <- rowMeans(counts[, 1:n1])
mean_ca <- rowMeans(counts[, (n1 + 1):(n1 + n2)])
filtered_counts <- counts[mean_na >= C & mean_ca >= C, ]

# Rank-based transformation per sample
ranked_counts <- apply(filtered_counts, 2, function(x) rank(x, ties.method = "average"))

# Mann-Whitney U test for each gene
p_values <- apply(ranked_counts, 1, function(x) {
  test_result <- wilcox.test(x[1:n1], x[(n1 + 1):(n1 + n2)])
  return(test_result$p.value)
})

# Adjust p-values (Benjamini-Hochberg)
padj_values <- p.adjust(p_values, method = "BH")
sig_ind <- padj_values < 0.01

# Compute log2 fold change (using scaled values)
normalized_counts <- scale(filtered_counts, center = FALSE, scale = TRUE)
log2fc <- log2(rowMeans(normalized_counts[, 1:n1]) /
                 rowMeans(normalized_counts[, (n1 + 1):(n1 + n2)]))

# Identify significant DEGs with large fold change
log2_fc_cutoff <- 2
log2_fc_keep <- abs(log2fc) >= log2_fc_cutoff
deg_ind <- sig_ind & log2_fc_keep

# Result table
res <- data.frame(
  log2FoldChange = log2fc,
  P_Value = p_values,
  Padj = padj_values,
  Upregulated_group = ifelse(log2fc > 0, "Native American", "Caucasian"),
  row.names = rownames(filtered_counts)
)

# Subset: DEGs only
res_deg <- res[deg_ind, ]

# Save result tables
write.csv(res, "results/DEG_results_all_genes.csv")
write.csv(res_deg, "results/DEG_results_significant.csv")

# Volcano plot
png("figures/Volcano_Plot_Breast_DEG.png", width = 3000, height = 2400, res = 300)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'Padj',
                pCutoff = 0.01,
                FCcutoff = 2,
                title = 'Volcano Plot: Breast Cancer DEGs',
                xlab = 'Log2 Fold Change',
                ylab = '-Log10 Adjusted p-value',
                labSize = 4.0,
                col = c('grey', 'red', 'blue', 'green'),
                pointSize = 1.5,
                drawConnectors = TRUE,
                max.overlaps = 18,
                selectLab = c("GNG4", "CBLN2", "OLFM4", "MSLN", "CHGB",
                              "CRISP3", "SLC30A8", "ZPLD1", "POTEF", "POTEC",
                              "POTEJ", "CTAGE4", "APOB", "ANKRD30A", "RGPD8"))
dev.off()
