library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(EnhancedVolcano)
library(pheatmap)
library(enrichplot)
library(ggplot2)

# Load DEG gene symbols
deg_genes <- read.csv("results/DEG_results_significant.csv", header = TRUE)
gene_symbols <- deg_genes$X  # SYMBOLs are in column "X"

# Map gene symbols to Entrez IDs
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = gene_symbols,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")
entrez_ids <- na.omit(entrez_ids)

# Prepare KEGG gene sets from MSigDB (C2:CP:KEGG)
msigdb_kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")

# KEGG pathway enrichment
ego_kegg <- enricher(entrez_ids,
                     TERM2GENE = msigdb_kegg[, c("gs_name", "entrez_gene")],
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05)
ego_kegg <- setReadable(ego_kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# GO enrichment (Biological Process only)
ego_go <- enrichGO(gene = entrez_ids,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   readable = TRUE)

# Dotplot: GO
p_go <- dotplot(ego_go, showCategory = 10, title = "GO Enrichment Analysis (Biological Process)") +
  theme(
    text = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    plot.title = element_text(size = 16, hjust = 0.5)
  )

# Dotplot: KEGG
p_kegg <- dotplot(ego_kegg, showCategory = 10, title = "KEGG Enrichment Analysis") +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    plot.title = element_text(size = 16, hjust = 0.5)
  )

# Save dotplots as PDF
ggsave("figures/go_enrichment_dotplot.pdf", plot = p_go, width = 8, height = 8)
ggsave("figures/kegg_enrichment_dotplot.pdf", plot = p_kegg, width = 8, height = 8)

# Save enrichment results
write.csv(as.data.frame(ego_go), "results/go_enrichment_results.csv", row.names = FALSE)
write.csv(as.data.frame(ego_kegg), "results/kegg_enrichment_results.csv", row.names = FALSE)
