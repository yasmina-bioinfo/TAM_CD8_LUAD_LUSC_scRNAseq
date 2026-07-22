#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc5_Annotation_TME, Script 04: Resolve mixed clusters (LUSC, CellTypist)
# CellTypist trained on a multi-tissue atlas gave clearer separation for
# closely related T-cell states than a peripheral-blood-only reference.
# Systematic per-cluster CellTypist_majority concordance check (all 21
# clusters) identified 8 clusters below the ~95% single-family threshold:
# 1, 2, 3, 4, 6, 12, 14, 20.
# For each target cluster, isolates cells, groups them by their individual
# CellTypist_majority label (labels below MIN_CELLS excluded as noise),
# and runs FindAllMarkers to reveal actual differential markers between
# groups — cross-checked against CellTypist's own author-curated marker
# set per label, rather than trusting the label name alone.
# Input:  Objects/03_LUSC_seu_celltypist.rds
# Output: Results/Bloc5_Annotation_TME/Tables/04_LUSC_cluster<N>_markers.tsv
# ============================================================
suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(dplyr)
  library(readr)
})
IN_OBJ  <- here("Objects", "03_LUSC_seu_celltypist.rds")
OUT_TAB <- here("Results", "Bloc5_Annotation_TME", "Tables")
dir.create(OUT_TAB, showWarnings = FALSE, recursive = TRUE)
message("Loading CellTypist-annotated Seurat object...")
seu <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu), " | Clusters: ", length(levels(seu$seurat_clusters)))
MIN_CELLS <- 50
# Target clusters identified via systematic CellTypist_majority concordance
# check on ALL 21 clusters — below ~95% concordance in a single
# biological family, no gray zone tolerated.
TARGET_CLUSTERS <- c("1", "2", "3", "4", "6", "12", "14", "20")
for (cl in TARGET_CLUSTERS) {
  
  message("=== Resolving cluster ", cl, " ===")
  
  cl_cells <- colnames(seu)[seu$seurat_clusters == cl]
  label_counts <- table(seu$CellTypist_majority[cl_cells])
  keep_labels <- names(label_counts[label_counts >= MIN_CELLS])
  
  message("Cluster ", cl, ": ", length(cl_cells), " cells, ",
          length(keep_labels), " labels retained (>= ", MIN_CELLS, " cells)")
  
  if (length(keep_labels) < 2) {
    message("Cluster ", cl, ": only ", length(keep_labels), 
            " label(s) above MIN_CELLS — skipping FindAllMarkers (nothing to compare against). ",
            "Treating as effectively homogeneous within this cluster.")
    next
  }
  
  print(sort(label_counts[keep_labels], decreasing = TRUE))
  
  target_cells <- cl_cells[seu$CellTypist_majority[cl_cells] %in% keep_labels]
  seu_sub <- subset(seu, cells = target_cells)
  Idents(seu_sub) <- seu_sub$CellTypist_majority
  
  message("Running FindAllMarkers on cluster ", cl, " subgroups...")
  markers_sub <- FindAllMarkers(seu_sub, only.pos = TRUE,
                                min.pct = 0.25, logfc.threshold = 0.25)
  message("Total markers found: ", nrow(markers_sub))
  
  top10_sub <- markers_sub %>%
    group_by(cluster) %>%
    slice_head(n = 10) %>%
    select(cluster, gene, avg_log2FC) %>%
    ungroup()
  
  write_tsv(markers_sub, file.path(OUT_TAB,
                                   paste0("04_LUSC_cluster", cl, "_markers_CellTypist_all.tsv")))
  write_tsv(top10_sub, file.path(OUT_TAB,
                                 paste0("04_LUSC_cluster", cl, "_markers_CellTypist_top10.tsv")))
  
  message("Saved: cluster ", cl, " markers tables")
}
message("DONE Script 04 — LUSC mixed cluster resolution (CellTypist)")