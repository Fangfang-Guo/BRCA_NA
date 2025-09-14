library(ggplot2)
library(dplyr)
library(ggrepel)
library(readr)

df <- read_tsv("results/cnv_fisher_results.tsv")

# ------------------- Change -------------------
df <- df %>%
  mutate(
    log2_Odds = log2(Change_Odds_Ratio),
    neg_log10_Padj = -log10(Change_Padj)
  )

x_min <- floor(min(df$log2_Odds[is.finite(df$log2_Odds)], na.rm = TRUE)) - 1
x_max <- ceiling(max(df$log2_Odds[is.finite(df$log2_Odds)], na.rm = TRUE)) + 1
y_max <- ceiling(quantile(df$neg_log10_Padj[is.finite(df$neg_log10_Padj)], 0.99, na.rm = TRUE)) + 1

df <- df %>%
  mutate(
    log2_Odds = case_when(
      is.infinite(log2_Odds) & Change_Odds_Ratio == 0 ~ x_min,
      is.infinite(log2_Odds) & Change_Odds_Ratio == Inf ~ x_max,
      TRUE ~ log2_Odds
    ),
    neg_log10_Padj = ifelse(is.infinite(neg_log10_Padj), y_max, neg_log10_Padj)
  )

sig_genes <- df %>% filter(Change_Padj < 0.05)
left_genes <- sig_genes %>% filter(log2_Odds < 0) %>% arrange(-neg_log10_Padj) %>% slice(1:16)
right_genes <- sig_genes %>% filter(log2_Odds > 0)

ggplot(df, aes(x = log2_Odds, y = neg_log10_Padj)) +
  geom_point(aes(color = Change_Padj < 0.05), alpha = 0.6) +
  scale_color_manual(name = "Significance", values = c("gray", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = left_genes, aes(label = Gene), size = 3.5, box.padding = 0.6, force = 20) +
  geom_text_repel(data = right_genes, aes(label = Gene), size = 3.5, box.padding = 0.6) +
  labs(title = "Volcano Plot: CNV Change", x = "log2(Odds Ratio)", y = "-log10(Padj)") +
  theme_minimal() +
  xlim(x_min, x_max) +
  ylim(0, y_max)

# ------------------- Gain -------------------
df_gain <- df %>%
  mutate(
    log2_Gain_Odds = log2(Gain_Odds_Ratio),
    neg_log10_Gain_Padj = -log10(Gain_Padj)
  )

x_min_gain <- floor(min(df_gain$log2_Gain_Odds[is.finite(df_gain$log2_Gain_Odds)], na.rm = TRUE)) - 1
x_max_gain <- ceiling(max(df_gain$log2_Gain_Odds[is.finite(df_gain$log2_Gain_Odds)], na.rm = TRUE)) + 1
y_max_gain <- ceiling(quantile(df_gain$neg_log10_Gain_Padj[is.finite(df_gain$neg_log10_Gain_Padj)], 0.99, na.rm = TRUE)) + 1

df_gain <- df_gain %>%
  mutate(
    log2_Gain_Odds = case_when(
      is.infinite(log2_Gain_Odds) & Gain_Odds_Ratio == 0 ~ x_min_gain,
      is.infinite(log2_Gain_Odds) & Gain_Odds_Ratio == Inf ~ x_max_gain,
      TRUE ~ log2_Gain_Odds
    ),
    neg_log10_Gain_Padj = ifelse(is.infinite(neg_log10_Gain_Padj), y_max_gain, neg_log10_Gain_Padj)
  )

sig_gain_genes <- df_gain %>% filter(Gain_Padj < 0.05)
left_gain_genes <- sig_gain_genes %>% filter(log2_Gain_Odds < 0) %>% arrange(Gain_Padj) %>% slice(1:15)
right_gain_genes <- sig_gain_genes %>% filter(log2_Gain_Odds > 0) %>% arrange(Gain_Padj)

ggplot(df_gain, aes(x = log2_Gain_Odds, y = neg_log10_Gain_Padj)) +
  geom_point(aes(color = Gain_Padj < 0.05), alpha = 0.6) +
  scale_color_manual(name = "Significance", values = c("gray", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = left_gain_genes, aes(label = Gene), size = 3.2, box.padding = 0.6, force = 20) +
  geom_text_repel(data = right_gain_genes, aes(label = Gene), size = 3.2, box.padding = 0.6) +
  labs(title = "Volcano Plot: CNV Gain", x = "log2(Odds Ratio)", y = "-log10(Padj)") +
  theme_minimal() +
  xlim(x_min_gain, x_max_gain) +
  ylim(0, y_max_gain)

# ------------------- Loss -------------------
df_loss <- df %>%
  mutate(
    log2_Loss_Odds = log2(Loss_Odds_Ratio),
    neg_log10_Loss_Padj = -log10(Loss_Padj)
  )

x_min_loss <- floor(min(df_loss$log2_Loss_Odds[is.finite(df_loss$log2_Loss_Odds)], na.rm = TRUE)) - 1
x_max_loss <- ceiling(max(df_loss$log2_Loss_Odds[is.finite(df_loss$log2_Loss_Odds)], na.rm = TRUE)) + 1
y_max_loss <- ceiling(quantile(df_loss$neg_log10_Loss_Padj[is.finite(df_loss$neg_log10_Loss_Padj)], 0.99, na.rm = TRUE)) + 1

df_loss <- df_loss %>%
  mutate(
    log2_Loss_Odds = case_when(
      is.infinite(log2_Loss_Odds) & Loss_Odds_Ratio == 0 ~ x_min_loss,
      is.infinite(log2_Loss_Odds) & Loss_Odds_Ratio == Inf ~ x_max_loss,
      TRUE ~ log2_Loss_Odds
    ),
    neg_log10_Loss_Padj = ifelse(is.infinite(neg_log10_Loss_Padj), y_max_loss, neg_log10_Loss_Padj)
  )

sig_loss_genes <- df_loss %>% filter(Loss_Padj < 0.05)
right_loss_genes <- sig_loss_genes %>% filter(log2_Loss_Odds > 0) %>% arrange(Loss_Padj) %>% slice(1:16)
left_loss_genes <- sig_loss_genes %>% filter(log2_Loss_Odds < 0) %>% arrange(Loss_Padj) %>% slice(1:20)

ggplot(df_loss, aes(x = log2_Loss_Odds, y = neg_log10_Loss_Padj)) +
  geom_point(aes(color = Loss_Padj < 0.05), alpha = 0.6) +
  scale_color_manual(name = "Significance", values = c("gray", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = left_loss_genes, aes(label = Gene), size = 3) +
  geom_text_repel(data = right_loss_genes, aes(label = Gene), size = 3) +
  labs(title = "Volcano Plot: CNV Loss", x = "log2(Odds Ratio)", y = "-log10(Padj)") +
  theme_minimal() +
  xlim(x_min_loss, x_max_loss) +
  ylim(0, y_max_loss)
