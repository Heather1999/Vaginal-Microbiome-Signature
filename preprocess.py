"""
preprocess.py
-----------------------------------------------------------------------------
Bridges the raw sequence-processing outputs (Metaxa2 + Lotus3, see
01_sequence_processing.sh) to the R statistics script (permanova_maaslin2.R)
and the Python elastic net scripts (prediction_archaea.py, prediction_birth.py).

Corresponds to Methods Sections 3.2-3.3 (compositional transformation, CLR).
PERMANOVA/PERMDISP/MaAsLin2/the primary GLM (Sections 3.4, 4.1-4.3) are
implemented separately in permanova_maaslin2.R, as noted inline below.

-----------------------------------------------------------------------------
"""
import os
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from skbio.stats.composition import clr
from scipy.spatial.distance import pdist, squareform
from skbio.stats.distance import DistanceMatrix

OUTDIR = "permanova"
os.makedirs(OUTDIR, exist_ok=True)

# ----------------------------- 1) Load inputs --------------------------------
# counts: pandas DataFrame, rows=samples, cols=features (OTU/ASV), produced by
# 01_sequence_processing.sh -> aggregate_to_family.py
us_counts = pd.read_csv("us_counts_after_preprocess.tsv", sep="\t", index_col=0)
de_counts = pd.read_csv("de_counts_after_preprocess.tsv", sep="\t", index_col=0)
us_meta = pd.read_csv("us_meta_after_preprocess.tsv", sep="\t", index_col=0)
de_meta = pd.read_csv("de_meta_after_preprocess.tsv", sep="\t", index_col=0)

# archaea.txt: list of archaea-positive sample IDs from Metaxa2 screening
# (Section 2). The Danish cohort's "Archaea" label is expected to already be
# present in de_meta (produced upstream, alongside the 8-sample shotgun
# confirmation in Section 6); this only needs to inject the label for the
# US cohort.
metaxa = pd.read_csv("archaea.txt", header=None)
us_meta["archaea"] = us_meta.index.isin(metaxa.iloc[:, 0]).astype(int)

# ----------------------------- 2) CLR transformation -------------------------
# Section 3.3: pseudocount of 0.5, selected below the global minimum non-zero
# count in the combined dataset (see manuscript for the 0.1/0.5/1.0
# sensitivity analysis).
common_taxa = us_counts.index.intersection(de_counts.index)  # keep taxa present in both cohorts
us_counts = us_counts.loc[common_taxa].T
de_counts = de_counts.loc[common_taxa].T
pca_input = pd.concat([us_counts, de_counts], axis=0)

non_zero_min = pca_input[pca_input != 0].min().min()
assert non_zero_min > 0.5, (
    f"Global minimum non-zero count is {non_zero_min}, which is not > 0.5 -- "
    "the paper's rationale for pseudocount=0.5 assumes it sits below this "
    "minimum. Re-check before proceeding."
)
pca_input = pca_input + 0.5

X_clr = pd.DataFrame(clr(pca_input.values), index=pca_input.index, columns=pca_input.columns)

# ----------------------------- 3) Exploratory PCA (Fig 1) --------------------
labels = ["US"] * us_counts.shape[0] + ["DE"] * de_counts.shape[0]
pca = PCA(n_components=2)
pc = pca.fit_transform(X_clr)
pca_df = pd.DataFrame(pc, columns=["PC1", "PC2"], index=X_clr.index)
pca_df["Group"] = labels

plt.figure(figsize=(6, 5))
for g, color in zip(["US", "DE"], ["tab:blue", "tab:orange"]):
    subset = pca_df[pca_df["Group"] == g]
    plt.scatter(subset["PC1"], subset["PC2"], label=g, alpha=0.7, color=color)
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
plt.legend()
plt.title("Exploratory PCA of CLR-transformed Family-level microbiome profiles")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "pca_us_de.pdf"), bbox_inches="tight")
plt.close()

# ----------------------------- 4) Distance matrix ----------------------------
# Euclidean distance on CLR-transformed data == Aitchison distance (Section 3.3).
# Beta-diversity significance testing (PERMANOVA/PERMDISP, Sections 3.4, 4.2)
# and the mixed-effects/GLM modeling (Sections 4.1, 4.3) are done in
# permanova_maaslin2.R, which reads R_meta.tsv and concated_clr.tsv below.
concat_clr = X_clr  # sample x taxa, matches vegan::vegdist expectations in R
D = squareform(pdist(concat_clr.values, metric="euclidean"))
dist = DistanceMatrix(D, ids=concat_clr.index)
df_dist = pd.DataFrame(dist.data, index=dist.ids, columns=dist.ids)

# ----------------------------- 5) Harmonized metadata ------------------------
de_sub_meta = de_meta[["PatID", "Archaea", "meta_Time", "Preterm"]]
us_sub_meta = us_meta[["host_subject_id", "archaea", "gest_day_collection", "term_vs_preterm_delivery"]]
us_sub_meta.columns = ["PatID", "Archaea", "meta_Time", "Preterm"]
archaea_meta = pd.concat([us_sub_meta, de_sub_meta], axis=0)

archaea_meta["cohort"] = np.where(
    archaea_meta.index.astype(str).str.startswith("S"), "US", "DK"
)

# ----------------------------- 6) Write outputs ------------------------------
archaea_meta.to_csv(os.path.join(OUTDIR, "R_meta.tsv"), sep="\t", index=True)
df_dist.to_csv(os.path.join(OUTDIR, "dist.tsv"), sep="\t", index=True)
concat_clr.to_csv(os.path.join(OUTDIR, "concated_clr.tsv"), sep="\t", index=True)

print(f"Wrote {OUTDIR}/R_meta.tsv, {OUTDIR}/dist.tsv, {OUTDIR}/concated_clr.tsv, "
      f"and {OUTDIR}/pca_us_de.pdf")
