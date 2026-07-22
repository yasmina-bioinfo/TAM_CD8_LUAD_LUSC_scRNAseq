#!/usr/bin/env python3
# ============================================================
# GSE243013 — Bloc5_Annotation_TME, Script 02: CellTypist annotation (LUSC)
# Reads the manually-exported log-normalized matrix (MatrixMarket format,
# see 02_LUSC_ExportForCellTypist.R), reconstructs it as an AnnData object,
# and runs CellTypist (Immune_All_Low model) with majority_voting enabled
# for a more robust per-cluster consensus label, less sensitive to
# individual-cell noise than a raw per-cell prediction alone.
# Input:  Objects/celltypist_export/LUSC_matrix.mtx, barcodes.tsv, features.tsv
# Output: Objects/celltypist_export/LUSC_celltypist_predictions.csv
# ============================================================
import scanpy as sc
import pandas as pd

IN_DIR = "Objects/celltypist_export"

# -----------------------------
# 1) Read the exported matrix and reconstruct an AnnData object
# -----------------------------
print("Reading exported matrix...")
adata = sc.read_mtx(f"{IN_DIR}/LUSC_matrix.mtx").T  # transpose: cells x genes

barcodes = pd.read_csv(f"{IN_DIR}/LUSC_barcodes.tsv", header=None)[0].values
features = pd.read_csv(f"{IN_DIR}/LUSC_features.tsv", header=None)[0].values

adata.obs_names = barcodes
adata.var_names = features

print(f"Cells: {adata.n_obs} | Genes: {adata.n_vars}")

# Sanity check: confirm the data is log-normalized, matching the R-side check
print(f"Max value in matrix (sanity check, expect a small float, ~9.02 like R): {adata.X.max():.2f}")

# -----------------------------
# 2) Run CellTypist (Immune_All_Low model, majority_voting enabled)
# -----------------------------
import celltypist
from celltypist import models

print("Running CellTypist annotation...")
predictions = celltypist.annotate(
    adata,
    model="Immune_All_Low.pkl",
    majority_voting=True
)

print("CellTypist annotation done.")
print(predictions.predicted_labels.head(10))

# -----------------------------
# 3) Save predictions to CSV, for import back into R
# -----------------------------
OUT_DIR = "Objects/celltypist_export"

predictions.predicted_labels.to_csv(f"{OUT_DIR}/LUSC_celltypist_predictions.csv")
print(f"Saved: {OUT_DIR}/LUSC_celltypist_predictions.csv")
