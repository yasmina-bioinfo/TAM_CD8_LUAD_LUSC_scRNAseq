#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc5_Annotation_TME, Script 05: Final annotation (LUSC, CellTypist)
# Final annotation built entirely from CellTypist_majority — the initial
# manual canonical-marker pass (Script01) is retained as an exploratory/
# sanity-check step only, not merged here, since CellTypist cross-checked
# against two source papers (Mullan et al. 2023; Domínguez Conde et al.
# 2022) provides a more robust and directly verifiable basis on its own.
# For clusters confirmed homogeneous (>=95% concordance in a single
# biological family), all cells in the cluster receive that cluster's
# dominant CellTypist label. For clusters resolved at the cell level
# (mixed clusters), each cell keeps its individual CellTypist_majority
# label if confirmed by official markers, otherwise an explicit
# "Unresolved" or "Other/contamination" tag.
# Input:  Objects/03_LUSC_seu_celltypist.rds
# Output: Objects/05_LUSC_seu_annotated.rds
#         Results/Bloc5_Annotation_TME/Tables/05_LUSC_final_annotation.csv
# ============================================================
suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(dplyr)
  library(data.table)
})

IN_OBJ  <- here("Objects", "03_LUSC_seu_celltypist.rds")
OUT_TAB <- here("Results", "Bloc5_Annotation_TME", "Tables")
dir.create(OUT_TAB, showWarnings = FALSE, recursive = TRUE)

message("Loading CellTypist-annotated Seurat object...")
seu <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu), " | Clusters: ", length(levels(seu$seurat_clusters)))

# -----------------------------
# 1) Clusters confirmed homogeneous directly (>=95% concordance) —
# dominant CellTypist label for the entire cluster
# -----------------------------
homogeneous_lookup <- c(
  "0"  = "Memory B cells",
  "5"  = "Regulatory T cells",
  "7"  = "Classical monocytes",
  "8"  = "Alveolar macrophages",
  "9"  = "Plasma cells",
  "10" = "DC2",
  "11" = "Mast cells",
  "13" = "Classical monocytes",
  "15" = "Germinal center B cells",
  "16" = "pDC",
  "17" = "Macrophages",
  "18" = "Alveolar macrophages",
  "19" = "Classical monocytes"
)

# -----------------------------
# 2) Mixed clusters, cell-level resolution
# -----------------------------
resolved_labels <- list(
  "1" = c(
    "Tem/Effector helper T cells" = "Tem/Effector helper T cells",
    "Tem/Trm cytotoxic T cells" = "Tem/Trm cytotoxic T cells",
    "MAIT cells" = "MAIT cells",
    "CRTAM+ gamma-delta T cells" = "CRTAM+ gamma-delta T cells",
    "Tem/Temra cytotoxic T cells" = "Tem/Temra cytotoxic T cells",
    "CD16- NK cells" = "CD16- NK cells",
    "Type 17 helper T cells" = "Type 17 helper T cells",
    "Type 1 helper T cells" = "Unresolved (cluster 1, CellTypist label: Type 1 helper T cells)",
    "Regulatory T cells" = "Unresolved (cluster 1, CellTypist label: Regulatory T cells)",
    "Trm cytotoxic T cells" = "Unresolved (cluster 1, CellTypist label: Trm cytotoxic T cells)",
    "gamma-delta T cells" = "Unresolved (cluster 1, CellTypist label: gamma-delta T cells)",
    "NK cells" = "Unresolved (cluster 1, CellTypist label: NK cells)"
  ),
  "2" = c(
    "CRTAM+ gamma-delta T cells" = "CRTAM+ gamma-delta T cells",
    "Tem/Trm cytotoxic T cells" = "Tem/Trm cytotoxic T cells",
    "MAIT cells" = "MAIT cells",
    "Regulatory T cells" = "Regulatory T cells",
    "CD16- NK cells" = "CD16- NK cells",
    "Trm cytotoxic T cells" = "Unresolved (cluster 2, CellTypist label: Trm cytotoxic T cells)",
    "Tem/Effector helper T cells" = "Other/contamination (cluster 2, CellTypist label: Tem/Effector helper T cells)",
    "Tem/Temra cytotoxic T cells" = "Unresolved (cluster 2, CellTypist label: Tem/Temra cytotoxic T cells)",
    "Type 1 helper T cells" = "Unresolved (cluster 2, CellTypist label: Type 1 helper T cells)",
    "gamma-delta T cells" = "Unresolved (cluster 2, CellTypist label: gamma-delta T cells)",
    "NK cells" = "Unresolved (cluster 2, CellTypist label: NK cells)"
  ),
  "3" = c(
    "Tem/Trm cytotoxic T cells" = "Tem/Trm cytotoxic T cells",
    "Trm cytotoxic T cells" = "Trm cytotoxic T cells",
    "MAIT cells" = "MAIT cells",
    "ILC3" = "ILC3",
    "Regulatory T cells" = "Regulatory T cells",
    "Tcm/Naive helper T cells" = "Tcm/Naive helper T cells",
    "Type 17 helper T cells" = "Type 17 helper T cells",
    "CD16- NK cells" = "CD16- NK cells",
    "Tem/Effector helper T cells" = "Unresolved (cluster 3, CellTypist label: Tem/Effector helper T cells)",
    "Type 1 helper T cells" = "Unresolved (cluster 3, CellTypist label: Type 1 helper T cells)",
    "Follicular helper T cells" = "Unresolved (cluster 3, CellTypist label: Follicular helper T cells)"
  ),
  "4" = c(
    "CD16+ NK cells" = "CD16+ NK cells",
    "Tem/Temra cytotoxic T cells" = "Tem/Temra cytotoxic T cells",
    "Tem/Trm cytotoxic T cells" = "Tem/Trm cytotoxic T cells",
    "CD16- NK cells" = "CD16- NK cells",
    "CRTAM+ gamma-delta T cells" = "CRTAM+ gamma-delta T cells",
    "MAIT cells" = "MAIT cells",
    "ILC3" = "ILC3",
    "Regulatory T cells" = "Regulatory T cells",
    "gamma-delta T cells" = "gamma-delta T cells",
    "Trm cytotoxic T cells" = "Unresolved (cluster 4, CellTypist label: Trm cytotoxic T cells)",
    "NK cells" = "Unresolved (cluster 4, CellTypist label: NK cells)",
    "Tem/Effector helper T cells" = "Other/contamination (cluster 4, CellTypist label: Tem/Effector helper T cells)"
  ),
  "6" = c(
    "Regulatory T cells" = "Regulatory T cells",
    "Tcm/Naive helper T cells" = "Tcm/Naive helper T cells",
    "Type 17 helper T cells" = "Type 17 helper T cells",
    "Tem/Trm cytotoxic T cells" = "Unresolved (cluster 6, CellTypist label: Tem/Trm cytotoxic T cells)",
    "Type 1 helper T cells" = "Unresolved (cluster 6, CellTypist label: Type 1 helper T cells)",
    "Follicular helper T cells" = "Unresolved (cluster 6, CellTypist label: Follicular helper T cells)",
    "Tem/Effector helper T cells" = "Other/contamination (cluster 6, CellTypist label: Tem/Effector helper T cells)"
  ),
  "12" = c(
    "CD16+ NK cells" = "CD16+ NK cells",
    "Regulatory T cells" = "Regulatory T cells",
    "Tem/Trm cytotoxic T cells" = "Tem/Trm cytotoxic T cells",
    "NK cells" = "NK cells",
    "CRTAM+ gamma-delta T cells" = "CRTAM+ gamma-delta T cells",
    "CD16- NK cells" = "CD16- NK cells",
    "Plasmablasts" = "Unresolved (cluster 12, CellTypist label: Plasmablasts)",
    "Epithelial cells" = "Other/contamination (cluster 12, CellTypist label: Epithelial cells)"
  ),
  "14" = c(
    "ILC3" = "ILC3",
    "CD16- NK cells" = "CD16- NK cells",
    "Tem/Effector helper T cells" = "Unresolved (cluster 14, CellTypist label: Tem/Effector helper T cells)",
    "Tem/Trm cytotoxic T cells" = "Unresolved (cluster 14, CellTypist label: Tem/Trm cytotoxic T cells)"
  ),
  "20" = c(
    "Tem/Trm cytotoxic T cells" = "Tem/Trm cytotoxic T cells"
  )
)

OTHER_LABEL <- "Other (mixed, minor)"

# -----------------------------
# 3) Build final cell-level annotation (vectorized, not cell-by-cell loop)
# -----------------------------
seu$annotation_final <- NA_character_

is_homogeneous <- as.character(seu$seurat_clusters) %in% names(homogeneous_lookup)
seu$annotation_final[is_homogeneous] <- unname(
  homogeneous_lookup[as.character(seu$seurat_clusters[is_homogeneous])]
)

for (cl in names(resolved_labels)) {
  in_cluster <- as.character(seu$seurat_clusters) == cl
  ct_labels <- as.character(seu$CellTypist_majority[in_cluster])
  map <- resolved_labels[[cl]]
  seu$annotation_final[in_cluster] <- ifelse(
    ct_labels %in% names(map), unname(map[ct_labels]), OTHER_LABEL
  )
}

message("Cells with NA annotation: ", sum(is.na(seu$annotation_final)))
message("Final annotation distribution:")
print(sort(table(seu$annotation_final, useNA = "ifany"), decreasing = TRUE))

# -----------------------------
# 4) Technical lineage flags for downstream Bloc6 (CD8) / Bloc7 (TAM)
# extraction — confirmed labels only, never "Unresolved"/"Other"
# -----------------------------
CD8_confirmed_labels <- c(
  "Tem/Trm cytotoxic T cells",
  "Tem/Temra cytotoxic T cells",
  "Trm cytotoxic T cells"
)
TAM_confirmed_labels <- c("Macrophages", "Alveolar macrophages")
TAM_related_labels <- c("Macrophages", "Alveolar macrophages", "Classical monocytes")

seu$lineage_CD8 <- seu$annotation_final %in% CD8_confirmed_labels
seu$lineage_TAM <- seu$annotation_final %in% TAM_confirmed_labels
seu$lineage_TAM_related <- seu$annotation_final %in% TAM_related_labels

message("CD8 lineage-confirmed cells: ", sum(seu$lineage_CD8))
message("TAM lineage-confirmed cells: ", sum(seu$lineage_TAM))
message("TAM-related cells (macrophages + classical monocytes, for Bloc7 reclustering): ", sum(seu$lineage_TAM_related))

saveRDS(seu, here("Objects", "05_LUSC_seu_annotated.rds"))

final_summary <- as.data.frame(table(seu$annotation_final))
colnames(final_summary) <- c("annotation_final", "n_cells")
fwrite(final_summary, file.path(OUT_TAB, "05_LUSC_final_annotation.csv"))

message("Saved: 05_LUSC_seu_annotated.rds, 05_LUSC_final_annotation.csv")
message("DONE Script 05 — LUSC final annotation (CellTypist)")