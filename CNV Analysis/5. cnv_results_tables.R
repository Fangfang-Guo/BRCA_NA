# Load necessary libraries
library(dplyr)

# Read CNV summary data
df <- read.csv("results/cnv_fisher_results.tsv", sep = "\t")

# Filter significantly different genes based on adjusted p-value < 0.05
change_sig <- df %>% filter(Change_Padj < 0.05)
gain_sig <- df %>% filter(Gain_Padj < 0.05)
loss_sig <- df %>% filter(Loss_Padj < 0.05)

# Select relevant columns for each category
change_table <- change_sig %>%
  select(Gene, Change_Odds_Ratio, Change_P_Value, Change_Padj, CNV_Change_Higher_In)

gain_table <- gain_sig %>%
  select(Gene, Gain_Odds_Ratio, Gain_P_Value, Gain_Padj, Gain_Higher_In)

loss_table <- loss_sig %>%
  select(Gene, Loss_Odds_Ratio, Loss_P_Value, Loss_Padj, Loss_Higher_In)

# Remove entries with missing directionality information
change_table <- change_table %>% filter(!is.na(CNV_Change_Higher_In))
gain_table <- gain_table %>% filter(!is.na(Gain_Higher_In))
loss_table <- loss_table %>% filter(!is.na(Loss_Higher_In))

# Save results to CSV files
write.csv(change_table, "Significant_CNV_Change_Genes.csv", row.names = FALSE)
write.csv(gain_table, "Significant_CNV_Gain_Genes.csv", row.names = FALSE)
write.csv(loss_table, "Significant_CNV_Loss_Genes.csv", row.names = FALSE)

# Preview results
head(change_table)
head(gain_table)
head(loss_table)
