#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc5_Annotation_TME, Script 02: Export for CellTypist (LUAD)
# Replaces the previous SingleR-based cross-check (MonacoImmuneData),
# found insufficiently robust to distinguish closely related T cell
# states (Th1/Th17/Tfh/MAIT) — see CAHIER_METHODOLOGIQUE.md. Manual
# export (MatrixMarket format) chosen over an automated Seurat->AnnData
# converter (e.g. SeuratDisk), to keep every step inspectable rather than
# relying on a black-box conversion — SeuratDisk has known compatibility
# issues with Seurat v5 objects in the community.
# Exports the log-normalized "data" layer (already produced by
# NormalizeData, scale.factor=1e4), matching the normalization CellTypist
# models expect (log1p of counts-per-10k).
# Input:  Objects/01_LUAD_seu_clustered.rds
# Output: Objects/celltypist_export/LUAD_matrix.mtx, barcodes.tsv, features.tsv
# ============================================================
suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(Matrix)
})

IN_OBJ  <- here("Objects", "01_LUAD_seu_clustered.rds")
OUT_DIR <- here("Objects", "celltypist_export")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

message("Loading clustered Seurat object...")
seu <- readRDS(IN_OBJ)
seu <- JoinLayers(seu)
message("Cells: ", ncol(seu), " | Genes: ", nrow(seu))

# -----------------------------
# 1) Extract log-normalized data (genes x cells, sparse)
# -----------------------------
data_matrix <- GetAssayData(seu, layer = "data")
message("Matrix class: ", class(data_matrix)[1])
message("Matrix dims: ", nrow(data_matrix), " genes x ", ncol(data_matrix), " cells")

# Sanity check before export: confirm this is log-normalized (not raw
# counts) — max value should be a modest float, not large integers
message("Max value in matrix (sanity check, expect a small float, not raw counts): ",
        round(max(data_matrix), 2))

# -----------------------------
# 2) Write MatrixMarket format (genes x cells) + barcode/feature lists
# -----------------------------
writeMM(data_matrix, file.path(OUT_DIR, "LUAD_matrix.mtx"))
write.table(colnames(data_matrix), file.path(OUT_DIR, "LUAD_barcodes.tsv"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(rownames(data_matrix), file.path(OUT_DIR, "LUAD_features.tsv"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

message("Saved: LUAD_matrix.mtx, LUAD_barcodes.tsv, LUAD_features.tsv")
message("DONE Script 02 — LUAD export for CellTypist")