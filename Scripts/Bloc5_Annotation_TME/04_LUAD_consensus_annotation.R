#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc5_Annotation_TME, Script 04: Consensus annotation (LUAD)
# Combines manual canonical-marker annotation (15 homogeneous/confirmed
# clusters, including 2 confirmed via SingleR cross-check despite a
# cell-cycle-dominated top10) with cell-level resolution for the 2 truly
# mixed clusters (6 and 11), based on targeted FindAllMarkers per
# SingleR_fine subgroup rather than an assumed mapping to existing
# categories (avoids confirmation bias).
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
dir.create(OUT_TAB, showWarnings = FALSE, recursive = TRUE)

message("Loading SingleR-annotated Seurat object...")
seu <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu), " | Clusters: ", length(levels(seu$seurat_clusters)))

# -----------------------------
# Manual annotation, 15 homogeneous/confirmed clusters (6 and 11 excluded
# here — resolved at cell level below). Clusters 13 and 16 are confirmed
# via SingleR family-level cross-check (>=95% concordance within one
# biological family) despite a cell-cycle-dominated top10 marker set.
# -----------------------------
manual_lookup <- c(
  "0"  = "B cell",
  "1"  = "CD4 T (naive/memory)",
  "2"  = "CD8 T (GZMK+ effector/memory)",
  "3"  = "CD8 T (LAG3+ exhausted/TRM)",
  "4"  = "Treg",
  "5"  = "Macrophage/TAM (C1Q+/lipid-associated)",
  "7"  = "Monocyte classique (S100A8/9+/FCN1+)",
  "8"  = "CD4 T (CXCL13+ exhausted/Tfh-like)",
  "9"  = "DC (cDC2, CD1A+/FCER1A+)",
  "10" = "Plasma cell",
  "12" = "Mast cell",
  "13" = "Macrophage/TAM (lipid-associated, cycling)",
  "14" = "pDC",
  "15" = "GC B cell",
  "16" = "Plasma cell (cycling/plasmablast)"
)

# -----------------------------
# Cluster 6: mixed mature-lineage cluster, resolved cell-by-cell via
# SingleR_fine (no cell-cycle signal — distinct in nature from cluster 11)
# -----------------------------
cluster6_map <- c(
  "Natural killer cells"          = "NK cell",
  "Terminal effector CD8 T cells" = "CD8 T (terminal effector, cluster 6)",
  "Non-Vd2 gd T cells"            = "T γδ (Vδ1/non-Vδ2)",
  "Effector memory CD8 T cells"   = "CD8 T (effector memory, cluster 6 variant)",
  "Vd2 gd T cells"                = "T γδ (Vδ2)",
  "MAIT cells"                    = "MAIT cell",
  "Terminal effector CD4 T cells" = "CD4 T (terminal effector)",
  "Progenitor cells"              = "ILC3"
)

# -----------------------------
# Cluster 11: mixed cycling cluster, resolved cell-by-cell via SingleR_fine.
# Two SingleR labels ("Myeloid dendritic cells", "Intermediate monocytes")
# both map to the same true identity (TAM, lipid-associated) — a known
# MonacoImmuneData labeling artifact for this signature (no dedicated
# tissue-macrophage category in a peripheral-blood-derived reference).
# -----------------------------
cluster11_map <- c(
  "Effector memory CD8 T cells"   = "CD8 T (cytotoxic effector/memory, cycling)",
  "Myeloid dendritic cells"       = "Macrophage/TAM (lipid-associated, cycling)",
  "Intermediate monocytes"        = "Macrophage/TAM (lipid-associated, cycling)",
  "T regulatory cells"            = "Treg (cycling)",
  "Vd2 gd T cells"                = "T γδ (Vδ2, cycling)",
  "Non-Vd2 gd T cells"            = "T γδ (Vδ1/non-Vδ2, cycling)",
  "MAIT cells"                    = "MAIT cell (cycling)",
  "Switched memory B cells"       = "B cell (cycling)",
  "Exhausted B cells"             = "B cell (cycling)",
  "Natural killer cells"          = "NK cell (cycling)",
  "Th1 cells"                     = "CD4 T (Th1-like dysfunctional, cycling)",
  "Terminal effector CD8 T cells" = "CD8 T (KLRC2+ NK-like terminal, cycling)"
)

# -----------------------------
# Build final cell-level annotation
# -----------------------------
# Cells in clusters 6/11 whose SingleR_fine label fell below the MIN_CELLS
# threshold used in Script03 (excluded from FindAllMarkers, no
# marker-verified identity) are assigned a residual label rather than
# their raw, unverified SingleR_fine tag.
OTHER_LABEL <- "Other (mixed, minor)"

seu$annotation_consensus <- unname(case_when(
  seu$seurat_clusters == "6" ~
    ifelse(seu$SingleR_fine %in% names(cluster6_map),
           cluster6_map[seu$SingleR_fine], OTHER_LABEL),
  seu$seurat_clusters == "11" ~
    ifelse(seu$SingleR_fine %in% names(cluster11_map),
           cluster11_map[seu$SingleR_fine], OTHER_LABEL),
  TRUE ~ manual_lookup[as.character(seu$seurat_clusters)]
))

message("Cells with NA annotation: ", sum(is.na(seu$annotation_consensus)))
message("Final annotation distribution:")
print(sort(table(seu$annotation_consensus, useNA = "ifany"), decreasing = TRUE))

# -----------------------------
# Cluster-level summary table
# -----------------------------
consensus_summary <- data.frame(
  cluster = as.integer(names(manual_lookup)),
  manual_annotation = unname(manual_lookup)
)
consensus_summary$note <- "homogeneous or SingleR-confirmed cluster, manual annotation retained"

mixed_summary <- data.frame(
  cluster = c(6, 11),
  manual_annotation = c(
    "MIXED - resolved at cell level (8 lineages, no cell-cycle signal)",
    "MIXED - resolved at cell level (cycling, multiple lineages)"
  ),
  note = c(
    "see 03_LUAD_cluster6_markers_top10.tsv for resolution detail",
    "see 03_LUAD_cluster11_markers_top10.tsv for resolution detail"
  )
)

consensus_summary <- bind_rows(consensus_summary, mixed_summary) %>%
  arrange(cluster)

fwrite(consensus_summary, file.path(OUT_TAB, "04_LUAD_consensus_annotation.csv"))
saveRDS(seu, here("Objects", "04_LUAD_seu_annotated.rds"))
message("DONE Script 04 — LUAD consensus annotation")