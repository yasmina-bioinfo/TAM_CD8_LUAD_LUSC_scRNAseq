#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc5_Annotation_TME, Script 03: Resolve mixed clusters (LUAD)
# Extended after a systematic SingleR_fine concordance check across ALL
# 17 clusters (not just the ones that looked suspicious from the top10
# marker list alone) revealed 5 additional mixed clusters beyond the
# originally identified 6 and 11: clusters 1, 2, 3, 8 fell below the
# ~95% single-family concordance threshold (64%-84%), and cluster 12
# showed a notable secondary "Progenitor cells" component requiring the
# same CD34/CD38 verification already used for cluster 6's ILC3 subgroup.
# No gray zone tolerated: any cluster below ~95% concordance in a single
# biological family is resolved at the cell level, regardless of how
# close to the threshold it appears.
# For each target cluster, isolates cells, groups them by their
# individual SingleR_fine label (labels below MIN_CELLS excluded as
# noise, not merged into larger groups by assumption), and runs
# FindAllMarkers to reveal actual differential markers between groups —
# unbiased by prior assumption about what those markers should be.
# Input:  Objects/02_LUAD_seu_singler.rds
# Output: Results/Bloc5_Annotation_TME/Tables/03_LUAD_cluster<N>_markers.tsv
# ============================================================
suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(dplyr)
  library(readr)
})

IN_OBJ  <- here("Objects", "02_LUAD_seu_singler.rds")
OUT_TAB <- here("Results", "Bloc5_Annotation_TME", "Tables")
dir.create(OUT_TAB, showWarnings = FALSE, recursive = TRUE)

message("Loading annotated Seurat object...")
seu <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu), " | Clusters: ", length(levels(seu$seurat_clusters)))

MIN_CELLS <- 50

# Full extended target list: 6 and 11 (originally identified), plus
# 1, 2, 3, 8, 12 (found on systematic re-check of all 17 clusters'
# SingleR_fine distribution — see CAHIER_METHODOLOGIQUE.md for the
# concordance table that led to this correction). No threshold exception
# made for clusters close to 95% (e.g. cluster 3 at 83.9%).
TARGET_CLUSTERS <- c("1", "2", "3", "6", "8", "11", "12")

for (cl in TARGET_CLUSTERS) {
  
  message("=== Resolving cluster ", cl, " ===")
  
  cl_cells <- colnames(seu)[seu$seurat_clusters == cl]
  label_counts <- table(seu$SingleR_fine[cl_cells])
  keep_labels <- names(label_counts[label_counts >= MIN_CELLS])
  
  message("Cluster ", cl, ": ", length(cl_cells), " cells, ",
          length(keep_labels), " labels retained (>= ", MIN_CELLS, " cells)")
  print(sort(label_counts[keep_labels], decreasing = TRUE))
  
  target_cells <- cl_cells[seu$SingleR_fine[cl_cells] %in% keep_labels]
  seu_sub <- subset(seu, cells = target_cells)
  Idents(seu_sub) <- seu_sub$SingleR_fine
  
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
                                   paste0("03_LUAD_cluster", cl, "_markers_all.tsv")))
  write_tsv(top10_sub, file.path(OUT_TAB,
                                 paste0("03_LUAD_cluster", cl, "_markers_top10.tsv")))
  
  message("Saved: cluster ", cl, " markers tables")
}

message("DONE Script 03 — LUAD mixed cluster resolution (extended, no gray zone)")