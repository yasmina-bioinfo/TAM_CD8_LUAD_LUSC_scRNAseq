#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc5_Annotation_TME, Script 04: Consensus annotation (LUAD)
# Combines manual canonical-marker annotation (17 homogeneous clusters)
# with the resolved CD8 subtypes from cluster 11 (Script03), based on
# targeted FindAllMarkers rather than an assumed mapping to existing
# categories (avoids confirmation bias, see CAHIER_METHODOLOGIQUE.md).
# Input:  Objects/02_LUAD_seu_singler.rds
# Output: Objects/04_LUAD_seu_annotated.rds
#         Results/Bloc5_Annotation_TME/Tables/04_LUAD_consensus_annotation.csv
# ============================================================
suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(data.table)
  library(dplyr)
})

IN_OBJ  <- here("Objects", "02_LUAD_seu_singler.rds")
OUT_TAB <- here("Results", "Bloc5_Annotation_TME", "Tables")

message("Loading SingleR-annotated Seurat object...")
seu <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu), " | Clusters: ", length(levels(seu$seurat_clusters)))

# -----------------------------
# Manual annotation, 17 homogeneous clusters (cluster 11 excluded here)
# -----------------------------
manual_lookup <- c(
  "0" = "B cell",
  "1" = "CD4 T (naive/memory)",
  "2" = "CD8 T (GZMK+ effector/memory)",
  "3" = "CD8 T (LAG3+ exhausted/TRM)",
  "4" = "Macrophage/TAM (C1Q+/lipid-associated)",
  "5" = "Treg",
  "6" = "NK cell",
  "7" = "CD4 T (CXCL13+ exhausted/Tfh-like)",
  "8" = "Monocyte (classical, profile A)",
  "9" = "Plasma cell",
  "10" = "Mast cell",
  "12" = "Neutrophil",
  "13" = "DC (mature/migratory, mregDC)",
  "14" = "ILC3",
  "15" = "pDC",
  "16" = "GC B cell",
  "17" = "Monocyte (classical, profile B)"
)

# -----------------------------
# Cluster 11: resolved CD8 subtypes (Script03), no assumed mapping
# -----------------------------
cluster11_map <- c(
  "Effector memory CD8 T cells"   = "CD8 T (PDCD1+ exhausted, cycling)",
  "Central memory CD8 T cells"    = "CD8 T (TCF7+ naive/memory, cycling)",
  "Terminal effector CD8 T cells" = "CD8 T (KLRC2+ NK-like terminal, cycling)"
)

# -----------------------------
# Build final cell-level annotation
# -----------------------------
seu$annotation_consensus <- ifelse(
  seu$seurat_clusters == "11",
  ifelse(seu$SingleR_fine %in% names(cluster11_map),
         cluster11_map[seu$SingleR_fine],
         seu$SingleR_fine),  # non-CD8 cells in cluster 11 keep raw SingleR_fine label
  manual_lookup[as.character(seu$seurat_clusters)]
)

message("Cells with NA annotation: ", sum(is.na(seu$annotation_consensus)))
message("Final annotation distribution:")
print(sort(table(seu$annotation_consensus, useNA = "ifany"), decreasing = TRUE))

# -----------------------------
# Cluster-level summary table (for the 17 homogeneous clusters)
# -----------------------------
consensus_summary <- data.frame(
  cluster = as.integer(names(manual_lookup)),
  manual_annotation = unname(manual_lookup)
)
consensus_summary$note <- "homogeneous cluster, manual annotation retained"

cluster11_summary <- data.frame(
  cluster = 11,
  manual_annotation = "MIXED - resolved at cell level (3 CD8 subtypes + other lineages)",
  note = "see 03_LUAD_cluster11_CD8_markers_top10.tsv for resolution detail"
)

consensus_summary <- bind_rows(consensus_summary, cluster11_summary) %>%
  arrange(cluster)

dir.create(OUT_TAB, showWarnings = FALSE, recursive = TRUE)
fwrite(consensus_summary, file.path(OUT_TAB, "04_LUAD_consensus_annotation.csv"))

saveRDS(seu, here("Objects", "04_LUAD_seu_annotated.rds"))
message("DONE Script 04 — LUAD consensus annotation")