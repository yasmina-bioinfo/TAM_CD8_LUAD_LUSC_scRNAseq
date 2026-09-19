#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc5_Annotation_TME, Script 05: Final annotation (LUSC, CellTypist)
# Corrected cohort (173 patients), 20 clusters. Every CellTypist_majority
# subgroup, in every cluster, was verified individually against its
# official marker — no homogeneity threshold applied at the cluster
# level. Each cell keeps its individual CellTypist_majority label if
# confirmed by official markers, otherwise an explicit "Unresolved" or
# "Other/contamination" tag. Labels below the minimum group size in
# their cluster (never tested) fall to "Other (mixed, minor)".
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
# Per-cluster resolution, all 20 clusters
# -----------------------------
resolved_labels <- list(
  "0" = c(
    "Memory B cells" = "Memory B cells",
    "Naive B cells" = "Unresolved (cluster 0, CellTypist label: Naive B cells)",
    "Age-associated B cells" = "Unresolved (cluster 0, CellTypist label: Age-associated B cells)",
    "Germinal center B cells" = "Unresolved (cluster 0, CellTypist label: Germinal center B cells)"
  ),
  "1" = c(
    "Tem/Trm cytotoxic T cells" = "Tem/Trm cytotoxic T cells",
    "Regulatory T cells" = "Regulatory T cells",
    "Tcm/Naive helper T cells" = "Tcm/Naive helper T cells",
    "Type 1 helper T cells" = "Type 1 helper T cells",
    "MAIT cells" = "MAIT cells",
    "Type 17 helper T cells" = "Type 17 helper T cells",
    "ILC3" = "ILC3",
    "CRTAM+ gamma-delta T cells" = "CRTAM+ gamma-delta T cells",
    "Tem/Effector helper T cells" = "Unresolved (cluster 1, CellTypist label: Tem/Effector helper T cells)",
    "Follicular helper T cells" = "Unresolved (cluster 1, CellTypist label: Follicular helper T cells)",
    "Trm cytotoxic T cells" = "Unresolved (cluster 1, CellTypist label: Trm cytotoxic T cells)",
    "CD16+ NK cells" = "Unresolved (cluster 1, CellTypist label: CD16+ NK cells)",
    "CD16- NK cells" = "Unresolved (cluster 1, CellTypist label: CD16- NK cells)"
  ),
  "2" = c(
    "Tem/Trm cytotoxic T cells" = "Tem/Trm cytotoxic T cells",
    "Tem/Effector helper T cells" = "Tem/Effector helper T cells",
    "MAIT cells" = "MAIT cells",
    "CRTAM+ gamma-delta T cells" = "CRTAM+ gamma-delta T cells",
    "Tem/Temra cytotoxic T cells" = "Tem/Temra cytotoxic T cells",
    "Regulatory T cells" = "Regulatory T cells",
    "Trm cytotoxic T cells" = "Unresolved (cluster 2, CellTypist label: Trm cytotoxic T cells)",
    "Type 1 helper T cells" = "Unresolved (cluster 2, CellTypist label: Type 1 helper T cells)",
    "CD16- NK cells" = "Unresolved (cluster 2, CellTypist label: CD16- NK cells)",
    "CD16+ NK cells" = "Unresolved (cluster 2, CellTypist label: CD16+ NK cells)"
  ),
  "3" = c(
    "Tem/Trm cytotoxic T cells" = "Tem/Trm cytotoxic T cells",
    "CRTAM+ gamma-delta T cells" = "CRTAM+ gamma-delta T cells",
    "MAIT cells" = "MAIT cells",
    "NK cells" = "NK cells",
    "Trm cytotoxic T cells" = "Unresolved (cluster 3, CellTypist label: Trm cytotoxic T cells)",
    "CD16- NK cells" = "Unresolved (cluster 3, CellTypist label: CD16- NK cells)",
    "Tem/Effector helper T cells" = "Unresolved (cluster 3, CellTypist label: Tem/Effector helper T cells)",
    "Tem/Temra cytotoxic T cells" = "Unresolved (cluster 3, CellTypist label: Tem/Temra cytotoxic T cells)",
    "Type 1 helper T cells" = "Unresolved (cluster 3, CellTypist label: Type 1 helper T cells)",
    "Regulatory T cells" = "Unresolved (cluster 3, CellTypist label: Regulatory T cells)",
    "gamma-delta T cells" = "Unresolved (cluster 3, CellTypist label: gamma-delta T cells)"
  ),
  "4" = c(
    "CD16+ NK cells" = "CD16+ NK cells",
    "Tem/Trm cytotoxic T cells" = "Tem/Trm cytotoxic T cells",
    "MAIT cells" = "MAIT cells",
    "Regulatory T cells" = "Regulatory T cells",
    "Tem/Temra cytotoxic T cells" = "Unresolved (cluster 4, CellTypist label: Tem/Temra cytotoxic T cells)",
    "CD16- NK cells" = "Unresolved (cluster 4, CellTypist label: CD16- NK cells)",
    "CRTAM+ gamma-delta T cells" = "Unresolved (cluster 4, CellTypist label: CRTAM+ gamma-delta T cells)",
    "Trm cytotoxic T cells" = "Unresolved (cluster 4, CellTypist label: Trm cytotoxic T cells)",
    "Tem/Effector helper T cells" = "Unresolved (cluster 4, CellTypist label: Tem/Effector helper T cells)",
    "NK cells" = "Unresolved (cluster 4, CellTypist label: NK cells)",
    "ILC3" = "Unresolved (cluster 4, CellTypist label: ILC3)",
    "gamma-delta T cells" = "Unresolved (cluster 4, CellTypist label: gamma-delta T cells)"
  ),
  "5" = c(
    "Regulatory T cells" = "Regulatory T cells",
    "Tem/Trm cytotoxic T cells" = "Tem/Trm cytotoxic T cells",
    "Type 1 helper T cells" = "Type 1 helper T cells",
    "Tem/Effector helper T cells" = "Unresolved (cluster 5, CellTypist label: Tem/Effector helper T cells)",
    "Follicular helper T cells" = "Unresolved (cluster 5, CellTypist label: Follicular helper T cells)",
    "Tcm/Naive helper T cells" = "Unresolved (cluster 5, CellTypist label: Tcm/Naive helper T cells)"
  ),
  "6" = c(
    "Classical monocytes" = "Classical monocytes",
    "Non-classical monocytes" = "Non-classical monocytes",
    "DC2" = "DC2",
    "Macrophages" = "Macrophages",
    "Intermediate macrophages" = "Unresolved (cluster 6, CellTypist label: Intermediate macrophages)",
    "Alveolar macrophages" = "Unresolved (cluster 6, CellTypist label: Alveolar macrophages)"
  ),
  "7" = c(
    "Regulatory T cells" = "Regulatory T cells",
    "Tem/Trm cytotoxic T cells" = "Tem/Trm cytotoxic T cells",
    "Tcm/Naive helper T cells" = "Tcm/Naive helper T cells",
    "Type 1 helper T cells" = "Type 1 helper T cells",
    "Tem/Effector helper T cells" = "Unresolved (cluster 7, CellTypist label: Tem/Effector helper T cells)",
    "Follicular helper T cells" = "Unresolved (cluster 7, CellTypist label: Follicular helper T cells)"
  ),
  "8" = c(
    "DC2" = "DC2",
    "Classical monocytes" = "Classical monocytes",
    "Macrophages" = "Macrophages",
    "Alveolar macrophages" = "Unresolved (cluster 8, CellTypist label: Alveolar macrophages)",
    "Intestinal macrophages" = "Unresolved (cluster 8, CellTypist label: Intestinal macrophages)",
    "Intermediate macrophages" = "Unresolved (cluster 8, CellTypist label: Intermediate macrophages)",
    "Non-classical monocytes" = "Unresolved (cluster 8, CellTypist label: Non-classical monocytes)"
  ),
  "9" = c(
    "Plasma cells" = "Plasma cells"
  ),
  "10" = c(
    "DC2" = "DC2",
    "DC1" = "DC1",
    "Migratory DCs" = "Migratory DCs",
    "Classical monocytes" = "Classical monocytes",
    "Macrophages" = "Macrophages",
    "Intermediate macrophages" = "Unresolved (cluster 10, CellTypist label: Intermediate macrophages)",
    "Alveolar macrophages" = "Unresolved (cluster 10, CellTypist label: Alveolar macrophages)"
  ),
  "11" = c(
    "Mast cells" = "Mast cells"
  ),
  "12" = c(
    "Tem/Trm cytotoxic T cells" = "Tem/Trm cytotoxic T cells",
    "Regulatory T cells" = "Regulatory T cells",
    "NK cells" = "NK cells",
    "Plasmablasts" = "Unresolved (cluster 12, CellTypist label: Plasmablasts)",
    "Epithelial cells" = "Other/contamination (cluster 12, CellTypist label: Epithelial cells)"
  ),
  "13" = c(
    "ILC3" = "ILC3",
    "Tem/Trm cytotoxic T cells" = "Tem/Trm cytotoxic T cells",
    "Tem/Effector helper T cells" = "Unresolved (cluster 13, CellTypist label: Tem/Effector helper T cells)",
    "CD16- NK cells" = "Unresolved (cluster 13, CellTypist label: CD16- NK cells)"
  ),
  "14" = c(
    "pDC" = "pDC"
  ),
  "15" = c(
    "Plasma cells" = "Plasma cells",
    "Plasmablasts" = "Plasmablasts",
    "Proliferative germinal center B cells" = "Unresolved (cluster 15, CellTypist label: Proliferative germinal center B cells)",
    "Germinal center B cells" = "Unresolved (cluster 15, CellTypist label: Germinal center B cells)",
    "Memory B cells" = "Unresolved (cluster 15, CellTypist label: Memory B cells)"
  ),
  "16" = c(
    "Macrophages" = "Macrophages",
    "DC2" = "DC2",
    "DC1" = "DC1",
    "Alveolar macrophages" = "Unresolved (cluster 16, CellTypist label: Alveolar macrophages)"
  ),
  "17" = c(
    "Classical monocytes" = "Classical monocytes",
    "DC2" = "Unresolved (cluster 17, CellTypist label: DC2)",
    "Non-classical monocytes" = "Unresolved (cluster 17, CellTypist label: Non-classical monocytes)"
  ),
  "18" = c(
    "Classical monocytes" = "Classical monocytes",
    "Non-classical monocytes" = "Non-classical monocytes"
  ),
  "19" = c(
    "Tem/Trm cytotoxic T cells" = "Tem/Trm cytotoxic T cells"
  )
)

OTHER_LABEL <- "Other (mixed, minor)"

# -----------------------------
# Build final cell-level annotation (vectorized, per cluster)
# -----------------------------
seu$annotation_final <- NA_character_

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
# Technical lineage flags, built from raw CellTypist_majority labels
# (not from annotation_final's confirmed/Unresolved status), so that
# cells with a transitional or non-canonical CellTypist label remain
# available under these flags rather than being dropped.
# -----------------------------
CD8_confirmed_labels <- c(
  "Tem/Trm cytotoxic T cells",
  "Tem/Temra cytotoxic T cells"
)
TAM_celltypist_labels <- c(
  "Macrophages", "Intermediate macrophages", "Alveolar macrophages"
)
Monocyte_celltypist_labels <- c(
  "Classical monocytes", "Non-classical monocytes"
)

seu$lineage_CD8 <- seu$annotation_final %in% CD8_confirmed_labels
seu$lineage_TAM <- seu$CellTypist_majority %in% TAM_celltypist_labels
seu$lineage_TAM_related <- seu$CellTypist_majority %in% Monocyte_celltypist_labels

message("CD8 lineage-confirmed cells: ", sum(seu$lineage_CD8))
message("TAM lineage cells (macrophage continuum, incl. transitional): ", sum(seu$lineage_TAM))
message("TAM-related cells (monocytes): ", sum(seu$lineage_TAM_related))

saveRDS(seu, here("Objects", "05_LUSC_seu_annotated.rds"))

final_summary <- as.data.frame(table(seu$annotation_final))
colnames(final_summary) <- c("annotation_final", "n_cells")
fwrite(final_summary, file.path(OUT_TAB, "05_LUSC_final_annotation.csv"))

message("Saved: 05_LUSC_seu_annotated.rds, 05_LUSC_final_annotation.csv")
message("DONE Script 05 — LUSC final annotation (CellTypist)")