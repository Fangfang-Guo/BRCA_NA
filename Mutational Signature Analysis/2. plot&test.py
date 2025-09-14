import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats

# ===============================
# Load signature activities (SBS96, 3-signature solution)
# ===============================
file_path = "mutational_signature_results/SBS96/All_Solutions/SBS96_3_Signatures/Activities/SBS96_S3_NMF_Activities.txt"
df = pd.read_csv(file_path, sep="\t")

# ===============================
# Add ancestry labels
# ===============================
# Number of samples after filtering
n_caucasian <- 683   
n_native    <- 16 
df["ancestry"] = ["Caucasian"] * n_caucasian + ["Native American"] * n_native

# ===============================
# Rename signature columns
# ===============================
df = df.rename(columns={
    "SBS96A": "MMR signature",
    "SBS96B": "AID/POLE signature",
    "SBS96C": "AZA signature"
})
signatures = ["MMR signature", "AID/POLE signature", "AZA signature"]

# ===============================
# Log-transform exposures for plotting
# ===============================
for sig in signatures:
    df[f"log_{sig}"] = np.log10(df[sig] + 1)

# ===============================
# Visualization: violin + strip plots
# ===============================
sns.set(style="whitegrid", font_scale=1.5)
fig, axes = plt.subplots(1, len(signatures), figsize=(18, 6))

for i, sig in enumerate(signatures):
    log_sig = f"log_{sig}"
    ax = axes[i]

    # Violin plot
    sns.violinplot(
        data=df, x="ancestry", y=log_sig,
        palette=["blue", "red"], ax=ax, inner=None, alpha=0.3
    )

    # Strip plot
    sns.stripplot(
        data=df, x="ancestry", y=log_sig, hue="ancestry",
        palette=["blue", "red"], jitter=0.2, alpha=0.7, ax=ax
    )

    ax.set_title(sig, fontsize=20)
    ax.set_ylabel("log10(1 + exposure)", fontsize=16)
    ax.set_xlabel("")
    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", labelsize=14)

    # Remove duplicate legends
    legend = ax.get_legend()
    if legend:
        legend.remove()

# Unified legend
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, ["Caucasian", "Native American"],
           loc="upper center", ncol=2, frameon=True, fontsize=14)

plt.suptitle("Breast Cancer Signatures", y=0.95, fontsize=18)
plt.tight_layout()
plt.savefig("breast_cancer_signatures_violin_strip.pdf", dpi=300, bbox_inches="tight")
plt.show()

# ===============================
# Statistical comparison (Mann-Whitney U test, normalized proportions)
# ===============================

# Normalize so each row sums to 1 (relative contribution of signatures)
df_norm = df.copy()
df_norm[signatures] = df_norm[signatures].div(df_norm[signatures].sum(axis=1), axis=0)

# Split groups
group1 = df_norm[df_norm["ancestry"] == "Caucasian"]
group2 = df_norm[df_norm["ancestry"] == "Native American"]

# Compare each signature
alpha = 0.05
for sig in signatures:
    mean1, median1 = group1[sig].mean(), group1[sig].median()
    mean2, median2 = group2[sig].mean(), group2[sig].median()

    print(f"Signature: {sig}")
    print(f"  Caucasian - Mean: {mean1:.4f}, Median: {median1:.4f}")
    print(f"  Native American - Mean: {mean2:.4f}, Median: {median2:.4f}")

    # Mann-Whitney U test
    stat, p_value = stats.mannwhitneyu(group1[sig], group2[sig], alternative="two-sided")
    print(f"  U statistic: {stat}, p-value: {p_value:.4g}")

    if p_value < alpha:
        print("  Significant difference between groups.")
    else:
        print("  No significant difference between groups.")

    # Report which ancestry has higher average contribution
    if mean1 > mean2:
        print(f"  {sig} is higher in Caucasian group.\n")
    else:
        print(f"  {sig} is higher in Native American group.\n")
