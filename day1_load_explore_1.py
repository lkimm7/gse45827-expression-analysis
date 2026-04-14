# ============================================================
# Gene Expression Analysis — Day 1: Load & Explore
# Dataset: GSE45827 (breast cancer gene expression, GEO)
# Download from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45827
# After downloading, rename the file to: gse45827.csv
# ============================================================

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Saves plots as files instead of opening windows

# ── 1. LOAD THE DATA ────────────────────────────────────────
# GEO series matrix files have metadata lines starting with "!"
# We skip those and read only the actual expression data
import gzip

# Data starts at line 69 (after all the "!" metadata lines)
print("Loading data from line 69...")
df = pd.read_csv("gse45827.csv", skiprows=68, index_col=0, sep="\t")

print("=== Dataset Overview ===")
print(f"Number of genes:   {df.shape[0]}")
print(f"Number of samples: {df.shape[1]}")
print()
print("First 5 rows:")
print(df.head())

# ── 2. BASIC STATISTICS ─────────────────────────────────────
print("\n=== Basic Statistics (first 5 genes) ===")
print(df.iloc[:5].T.describe())  # .T flips rows/cols so stats are per gene

# ── 3. CHECK FOR MISSING VALUES ─────────────────────────────
missing = df.isnull().sum().sum()
print(f"\nTotal missing values: {missing}")

# Drop any rows (genes) that have missing data
df_clean = df.dropna()
print(f"Genes remaining after cleaning: {df_clean.shape[0]}")

# ── 4. PLOT: Distribution of expression values ──────────────
# This shows us what the "typical" expression level looks like
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle("Gene Expression Distribution", fontsize=14, fontweight='bold')

# Left plot: all values flattened into one histogram
all_values = df_clean.values.flatten()
axes[0].hist(all_values, bins=60, color='steelblue', edgecolor='white', linewidth=0.3)
axes[0].set_title("All expression values")
axes[0].set_xlabel("Expression level")
axes[0].set_ylabel("Frequency")

# Right plot: mean expression per gene
gene_means = df_clean.mean(axis=1)
axes[1].hist(gene_means, bins=60, color='mediumseagreen', edgecolor='white', linewidth=0.3)
axes[1].set_title("Mean expression per gene")
axes[1].set_xlabel("Mean expression level")
axes[1].set_ylabel("Number of genes")

plt.tight_layout()
plt.savefig("day1_distribution.png", dpi=150)
print("\nPlot saved: day1_distribution.png")

# ── 5. SAVE CLEANED DATA FOR DAY 2 ──────────────────────────
df_clean.to_csv("gse45827_clean.csv")
print("Cleaned data saved: gse45827_clean.csv")