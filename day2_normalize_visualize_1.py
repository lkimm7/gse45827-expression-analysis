# ============================================================
# Gene Expression Analysis — Day 2: Normalize & Visualize
# Run this after day1_load_explore.py
# ============================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import seaborn as sns

# ── 1. LOAD CLEANED DATA ────────────────────────────────────
df = pd.read_csv("gse45827_clean.csv", index_col=0)
print(f"Loaded cleaned data: {df.shape[0]} genes x {df.shape[1]} samples")

# ── 2. NORMALIZE: Z-score per gene ──────────────────────────
# Z-score normalization centers each gene at 0 with std=1
# This lets us compare genes that have very different baseline levels
df_norm = df.subtract(df.mean(axis=1), axis=0).divide(df.std(axis=1), axis=0)

print("\nNormalization complete.")
print(f"Mean after normalization (should be ~0): {df_norm.values.mean():.4f}")
print(f"Std after normalization  (should be ~1): {df_norm.values.std():.4f}")

# ── 3. FILTER: Keep only high-variance genes ─────────────────
# High-variance genes are the most biologically interesting —
# they differ the most across samples
gene_variance = df_norm.var(axis=1)
top_genes = gene_variance.nlargest(50).index   # top 50 most variable genes
df_top = df_norm.loc[top_genes]

print(f"\nTop 50 high-variance genes selected for heatmap.")

# ── 4. PLOT: Heatmap of top 50 genes ────────────────────────
fig, ax = plt.subplots(figsize=(14, 10))

heatmap = ax.imshow(
    df_top.values,
    aspect='auto',
    cmap='RdBu_r',       # Red = high expression, Blue = low expression
    vmin=-2, vmax=2      # Clamp color scale to ±2 std deviations
)

# Labels
ax.set_title("Top 50 High-Variance Genes Across Samples\n(Z-score normalized)", 
             fontsize=13, fontweight='bold', pad=12)
ax.set_xlabel("Samples", fontsize=11)
ax.set_ylabel("Genes", fontsize=11)

# Gene names on y-axis (shortened to fit)
ax.set_yticks(range(len(df_top.index)))
ax.set_yticklabels([g[:20] for g in df_top.index], fontsize=6)

# Sample numbers on x-axis
ax.set_xticks(range(0, df_top.shape[1], max(1, df_top.shape[1]//10)))
ax.set_xticklabels(
    [f"S{i+1}" for i in range(0, df_top.shape[1], max(1, df_top.shape[1]//10))],
    fontsize=8
)

# Color bar
cbar = plt.colorbar(heatmap, ax=ax, fraction=0.02, pad=0.02)
cbar.set_label("Z-score (expression level)", fontsize=9)

plt.tight_layout()
plt.savefig("day2_heatmap.png", dpi=150)
print("Plot saved: day2_heatmap.png")

# ── 5. PLOT: Top 10 most variable genes — bar chart ─────────
top10 = gene_variance.nlargest(10)

fig2, ax2 = plt.subplots(figsize=(10, 5))
bars = ax2.barh(
    [g[:25] for g in top10.index[::-1]],
    top10.values[::-1],
    color='steelblue', edgecolor='white', linewidth=0.4
)
ax2.set_title("Top 10 Most Variable Genes", fontsize=13, fontweight='bold')
ax2.set_xlabel("Variance (after normalization)", fontsize=11)
ax2.set_ylabel("Gene", fontsize=11)
plt.tight_layout()
plt.savefig("day2_top10_variance.png", dpi=150)
print("Plot saved: day2_top10_variance.png")

# ── 6. SAVE NORMALIZED DATA FOR DAY 3 ───────────────────────
df_norm.to_csv("gse45827_normalized.csv")
print("Normalized data saved: gse45827_normalized.csv")
