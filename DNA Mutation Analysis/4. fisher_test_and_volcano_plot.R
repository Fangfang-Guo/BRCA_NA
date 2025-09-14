library(maftools)

# Load MAF files
na_maf <- read.maf(maf = "filtered_NA_final.maf")
ca_maf <- read.maf(maf = "filtered_CA_final.maf")

# Assign group labels
na_maf@clinical.data$Group <- "Native American"
ca_maf@clinical.data$Group <- "Caucasian"

# Gene-level mutation counts
na_summary <- getGeneSummary(na_maf)[, c("Hugo_Symbol", "MutatedSamples")]
ca_summary <- getGeneSummary(ca_maf)[, c("Hugo_Symbol", "MutatedSamples")]

# Merge and align gene tables
mutation_table <- merge(
  na_summary,
  ca_summary,
  by = "Hugo_Symbol",
  suffixes = c("_NA", "_CA"),
  all = TRUE
)
mutation_table[is.na(mutation_table)] <- 0

# Total sample numbers
total_na <- 17      
total_ca <- 689      

# Non-mutated sample counts
mutation_table$NonMutated_NativeAmerican <- total_na - mutation_table$MutatedSamples_NativeAmerican
mutation_table$NonMutated_Caucasian <- total_ca - mutation_table$MutatedSamples_Caucasian

# Fisher's exact test
p_values <- numeric(nrow(mutation_table))
odds_ratios <- numeric(nrow(mutation_table))

for (i in 1:nrow(mutation_table)) {
  table_2x2 <- matrix(
    c(
      mutation_table$MutatedSamples_NativeAmerican[i], mutation_table$NonMutated_NativeAmerican[i],
      mutation_table$MutatedSamples_Caucasian[i], mutation_table$NonMutated_Caucasian[i]
    ),
    nrow = 2,
    byrow = TRUE
  )
  fisher_res <- fisher.test(table_2x2)
  p_values[i] <- fisher_res$p.value
  odds_ratios[i] <- fisher_res$estimate
}

# Append results
mutation_table$p_value <- p_values
mutation_table$odds_ratio <- odds_ratios
mutation_table$adjusted_p_value <- p.adjust(p_values, method = "BH")

# Annotate significance and direction
mutation_table$significant <- mutation_table$adjusted_p_value < 0.05
mutation_table$mutation_tendency <- ifelse(
  mutation_table$significant,
  ifelse(mutation_table$odds_ratio > 1, "More in NA", "More in CA"),
  "No Significant Difference"
)

# Frequency difference
mutation_table$freq_diff <- mutation_table$MutatedSamples_NativeAmerican / total_na * 100 -
                             mutation_table$MutatedSamples_Caucasian / total_ca * 100

# -log10(FDR)
mutation_table$minus_log10_FDR <- -log10(mutation_table$adjusted_p_value)

# Save result table
write.csv(mutation_table, "mutation_frequency_comparison_results.csv", row.names = FALSE)

# ============================================
# Volcano plot
# ============================================

library(ggplot2)
library(ggrepel)

ggplot(mutation_table, aes(x = minus_log10_FDR, y = freq_diff, label = Hugo_Symbol)) + 
  geom_point(aes(color = significant), size = 3, alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.8, color = "black") +
  geom_vline(xintercept = -log10(0.1), linetype = "dashed", color = "blue", size = 0.8) +
  scale_color_manual(values = c("gray", "red")) +
  geom_text_repel(aes(label = ifelse(significant, Hugo_Symbol, "")), size = 4, max.overlaps = 20) +
  scale_x_continuous(
    limits = c(min(mutation_table$minus_log10_FDR) - 0.2, max(mutation_table$minus_log10_FDR) + 0.2),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(min(mutation_table$freq_diff) - 2, max(mutation_table$freq_diff) + 2),
    expand = c(0, 0)
  ) +
  labs(
    x = "-log10(FDR)",
    y = "Difference in Mutation Frequency (%)",
    title = "Mutation Frequency Comparison"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.line = element_line(size = 0.8, color = "black"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, hjust = 0.5),
    panel.grid = element_blank()
  )

# Save plot
ggsave("mutation_frequency_volcano_plot.png", width = 8, height = 6, dpi = 300)
