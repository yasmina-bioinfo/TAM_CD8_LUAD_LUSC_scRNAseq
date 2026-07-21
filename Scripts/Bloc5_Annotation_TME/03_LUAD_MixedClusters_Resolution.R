#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc5_Annotation_TME, Script 03: Resolve mixed clusters (LUAD)
# Generalized version: works for any mixed cluster identified after
# SingleR cross-check (clusters 6 and 11 on this run), not limited to
# CD8 subtypes. For each target cluster, isolates cells, groups them by
# their individual SingleR_fine label (labels below MIN_CELLS excluded
# as noise, not merged into larger groups by assumption), and runs
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

# Minimum cells per SingleR_fine label to be treated as its own group;
# below this threshold, a label likely reflects noise or misassignment
# rather than a real subpopulation worth its own marker analysis.
# Raised from 20 to 50 after observing that labels in the 20-30 cell
# range on cluster 6 produced marker profiles inconsistent with their
# own SingleR label (e.g. a 20-cell "Myeloid dendritic cells" group
# showing no myeloid markers at all) — insufficient statistical power
# for a reliable FindAllMarkers result at that size.
MIN_CELLS <- 50


# Target clusters requiring resolution, identified via SingleR cross-check
# (see Bloc5 Script02 concordance table)
TARGET_CLUSTERS <- c("6", "11")

for (cl in TARGET_CLUSTERS) {
  
  message("=== Resolving cluster ", cl, " ===")
  
  cl_cells <- colnames(seu)[seu$seurat_clusters == cl]
  label_counts <- table(seu$SingleR_fine[cl_cells])
  keep_labels <- names(label_counts[label_counts >= MIN_CELLS])
  
  message("Cluster ", cl, ": ", length(cl_cells), " cells, ",
          length(keep_labels), " labels retained (>= ", MIN_CELLS, " cells)")
  print(label_counts[keep_labels])
  
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
  print(top10_sub, n = 50)
  
  write_tsv(markers_sub, file.path(OUT_TAB,
                                   paste0("03_LUAD_cluster", cl, "_markers_all.tsv")))
  write_tsv(top10_sub, file.path(OUT_TAB,
                                 paste0("03_LUAD_cluster", cl, "_markers_top10.tsv")))
  
  message("Saved: cluster ", cl, " markers tables")
}

message("DONE Script 03 — LUAD mixed cluster resolution (clusters 6, 11)")