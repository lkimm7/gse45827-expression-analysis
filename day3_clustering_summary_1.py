# ============================================================
# Gene Expression Analysis — Day 3: Clustering & Summary
# Run this after day2_normalize_visualize.py
# ============================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

# ── 1. LOAD NORMALIZED DATA ─────────────────────────────────
df = pd.read_csv("gse45827_normalized.csv", index_col=0)
print(f"Loaded normalized data: {df.shape[0]} genes x {df.shape[1]} samples")

# We work sample-wise from here: transpose so rows = samples
df_T = df.T   # shape: samples x genes

# ── 2. PCA: Reduce to 2 dimensions ──────────────────────────
# PCA finds the directions of maximum variation in the data.
# Plotting the first 2 components lets us SEE how samples cluster.
pca = PCA(n_components=2, random_state=42)
pca_coords = pca.fit_transform(df_T.fillna(0))

variance_explained = pca.explained_variance_ratio_ * 100
print(f"\nPCA complete.")
print(f"PC1 explains {variance_explained[0]:.1f}% of variance")
print(f"PC2 explains {variance_explained[1]:.1f}% of variance")

# ── 3. K-MEANS CLUSTERING ────────────────────────────────────
# Group samples into 3 clusters (a reasonable starting point for cancer subtypes)
N_CLUSTERS = 3
kmeans = KMeans(n_clusters=N_CLUSTERS, random_state=42, n_init=10)
cluster_labels = kmeans.fit_predict(pca_coords)
print(f"\nK-means clustering: {N_CLUSTERS} clusters assigned to {len(cluster_labels)} samples")

# ── 4. PLOT: PCA scatter colored by cluster ──────────────────
colors = ['#E8593C', '#3B8BD4', '#2BAE7E']   # coral, blue, teal
cluster_names = [f"Cluster {i+1}" for i in range(N_CLUSTERS)]

fig, ax = plt.subplots(figsize=(9, 7))
for i in range(N_CLUSTERS):
    mask = cluster_labels == i
    ax.scatter(
        pca_coords[mask, 0], pca_coords[mask, 1],
        c=colors[i], label=cluster_names[i],
        s=60, alpha=0.8, edgecolors='white', linewidths=0.5
    )

ax.set_title("PCA of Gene Expression Samples\n(colored by K-means cluster)",
             fontsize=13, fontweight='bold', pad=12)
ax.set_xlabel(f"PC1 ({variance_explained[0]:.1f}% variance)", fontsize=11)
ax.set_ylabel(f"PC2 ({variance_explained[1]:.1f}% variance)", fontsize=11)
ax.legend(fontsize=10, framealpha=0.7)
ax.axhline(0, color='gray', linewidth=0.5, linestyle='--', alpha=0.4)
ax.axvline(0, color='gray', linewidth=0.5, linestyle='--', alpha=0.4)

plt.tight_layout()
plt.savefig("day3_pca_clusters.png", dpi=150)
print("Plot saved: day3_pca_clusters.png")

# ── 5. PLOT: Summary figure (3 panels) ───────────────────────
# This is the "hero" figure you'd show Dr. Park as a work sample
fig2, axes = plt.subplots(1, 3, figsize=(16, 5))
fig2.suptitle("Gene Expression Analysis — Summary", fontsize=14, fontweight='bold')

# Panel A: Expression distribution
all_vals = df.values.flatten()
axes[0].hist(all_vals[~np.isnan(all_vals)], bins=60,
             color='steelblue', edgecolor='white', linewidth=0.3)
axes[0].set_title("A. Expression distribution", fontsize=11)
axes[0].set_xlabel("Z-score")
axes[0].set_ylabel("Frequency")

# Panel B: Top 10 variance genes
gene_variance = df.var(axis=1)
top10 = gene_variance.nlargest(10)
axes[1].barh(
    [g[:20] for g in top10.index[::-1]], top10.values[::-1],
    color='mediumseagreen', edgecolor='white', linewidth=0.3
)
axes[1].set_title("B. Top 10 variable genes", fontsize=11)
axes[1].set_xlabel("Variance")

# Panel C: PCA clusters
for i in range(N_CLUSTERS):
    mask = cluster_labels == i
    axes[2].scatter(
        pca_coords[mask, 0], pca_coords[mask, 1],
        c=colors[i], label=cluster_names[i],
        s=40, alpha=0.8, edgecolors='white', linewidths=0.4
    )
axes[2].set_title("C. PCA — sample clusters", fontsize=11)
axes[2].set_xlabel(f"PC1 ({variance_explained[0]:.1f}%)")
axes[2].set_ylabel(f"PC2 ({variance_explained[1]:.1f}%)")
axes[2].legend(fontsize=8)

plt.tight_layout()
plt.savefig("day3_summary_figure.png", dpi=150)
print("Plot saved: day3_summary_figure.png  ← use this as your work sample!")

print("\n=== Project complete! Files to include in your GitHub repo: ===")
print("  day1_load_explore.py")
print("  day2_normalize_visualize.py")
print("  day3_clustering_summary.py")
print("  day1_distribution.png")
print("  day2_heatmap.png")
print("  day2_top10_variance.png")
print("  day3_pca_clusters.png")
print("  day3_summary_figure.png  ← attach this one to your email to Dr. Park")
