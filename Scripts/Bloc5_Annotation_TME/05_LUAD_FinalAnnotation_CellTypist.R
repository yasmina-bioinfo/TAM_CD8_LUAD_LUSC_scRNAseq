#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc5_Annotation_TME, Script 05: Final annotation (LUAD, CellTypist)
# Corrected cohort (61 patients), 20 clusters. Every CellTypist_majority
# subgroup, in every cluster, was verified individually against its
# official marker — no homogeneity threshold applied at the cluster
# level. Each cell keeps its individual CellTypist_majority label if
# confirmed by official markers, otherwise an explicit "Unresolved" or
# "Other/contamination" tag.
# Input:  Objects/03_LUAD_seu_celltypist.rds
# Output: Objects/05_LUAD_seu_annotated.rds
#         Results/Bloc5_Annotation_TME/Tables/05_LUAD_final_annotation.csv
# ============================================================
suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(dplyr)
  library(data.table)
})

IN_OBJ  <- here("Objects", "03_LUAD_seu_celltypist.rds")
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
    "Age-associated B cells" = "Age-associated B cells",
    "Memory B cells" = "Unresolved (cluster 0, CellTypist label: Memory B cells)",
    "Naive B cells" = "Unresolved (cluster 0, CellTypist label: Naive B cells)",
    "Germinal center B cells" = "Unresolved (cluster 0, CellTypist label: Germinal center B cells)"
  ),
  "1" = c(
    "Regulatory T cells" = "Regulatory T cells",
    "Tem/Effector helper T cells" = "Tem/Effector helper T cells",
    "Tem/Trm cytotoxic T cells" = "Tem/Trm cytotoxic T cells",
    "MAIT cells" = "MAIT cells",
    "CRTAM+ gamma-delta T cells" = "CRTAM+ gamma-delta T cells",
    "Tem/Temra cytotoxic T cells" = "Tem/Temra cytotoxic T cells",
    "Tcm/Naive helper T cells" = "Tcm/Naive helper T cells",
    "Trm cytotoxic T cells" = "Unresolved (cluster 1, CellTypist label: Trm cytotoxic T cells)",
    "Type 1 helper T cells" = "Unresolved (cluster 1, CellTypist label: Type 1 helper T cells)",
    "NK cells" = "Unresolved (cluster 1, CellTypist label: NK cells)",
    "CD16- NK cells" = "Unresolved (cluster 1, CellTypist label: CD16- NK cells)",
    "CD16+ NK cells" = "Unresolved (cluster 1, CellTypist label: CD16+ NK cells)"
  ),
  "2" = c(
    "Regulatory T cells" = "Regulatory T cells",
    "Trm cytotoxic T cells" = "Unresolved (cluster 2, CellTypist label: Trm cytotoxic T cells)",
    "Tem/Trm cytotoxic T cells" = "Tem/Trm cytotoxic T cells",
    "CRTAM+ gamma-delta T cells" = "CRTAM+ gamma-delta T cells",
    "CD16- NK cells" = "Unresolved (cluster 2, CellTypist label: CD16- NK cells)",
    "Tem/Effector helper T cells" = "Tem/Effector helper T cells",
    "NK cells" = "Unresolved (cluster 2, CellTypist label: NK cells)",
    "MAIT cells" = "MAIT cells",
    "Tem/Temra cytotoxic T cells" = "Unresolved (cluster 2, CellTypist label: Tem/Temra cytotoxic T cells)",
    "gamma-delta T cells" = "gamma-delta T cells"
  ),
  "3" = c(
    "Tem/Effector helper T cells" = "Unresolved (cluster 3, CellTypist label: Tem/Effector helper T cells)",
    "MAIT cells" = "MAIT cells",
    "Regulatory T cells" = "Regulatory T cells",
    "Tem/Trm cytotoxic T cells" = "Tem/Trm cytotoxic T cells",
    "Tcm/Naive helper T cells" = "Tcm/Naive helper T cells",
    "Follicular helper T cells" = "Unresolved (cluster 3, CellTypist label: Follicular helper T cells)",
    "Type 1 helper T cells" = "Unresolved (cluster 3, CellTypist label: Type 1 helper T cells)",
    "Trm cytotoxic T cells" = "Unresolved (cluster 3, CellTypist label: Trm cytotoxic T cells)",
    "ILC3" = "Unresolved (cluster 3, CellTypist label: ILC3)",
    "Epithelial cells" = "Other/contamination (cluster 3, CellTypist label: Epithelial cells)"
  ),
  "4" = c(
    "Macrophages" = "Macrophages",
    "Alveolar macrophages" = "Unresolved (cluster 4, CellTypist label: Alveolar macrophages)",
    "Intermediate macrophages" = "Unresolved (cluster 4, CellTypist label: Intermediate macrophages)",
    "DC2" = "DC2",
    "Classical monocytes" = "Classical monocytes",
    "Non-classical monocytes" = "Unresolved (cluster 4, CellTypist label: Non-classical monocytes)"
  ),
  "5" = c(
    "Regulatory T cells" = "Regulatory T cells",
    "Tem/Effector helper T cells" = "Unresolved (cluster 5, CellTypist label: Tem/Effector helper T cells)",
    "Tem/Trm cytotoxic T cells" = "Tem/Trm cytotoxic T cells",
    "Tcm/Naive helper T cells" = "Unresolved (cluster 5, CellTypist label: Tcm/Naive helper T cells)",
    "Type 1 helper T cells" = "Unresolved (cluster 5, CellTypist label: Type 1 helper T cells)"
  ),
  "6" = c(
    "Regulatory T cells" = "Unresolved (cluster 6, CellTypist label: Regulatory T cells)",
    "Tem/Effector helper T cells" = "Unresolved (cluster 6, CellTypist label: Tem/Effector helper T cells)",
    "Tem/Trm cytotoxic T cells" = "Tem/Trm cytotoxic T cells",
    "Follicular helper T cells" = "Unresolved (cluster 6, CellTypist label: Follicular helper T cells)",
    "Type 1 helper T cells" = "Unresolved (cluster 6, CellTypist label: Type 1 helper T cells)",
    "Tcm/Naive helper T cells" = "Unresolved (cluster 6, CellTypist label: Tcm/Naive helper T cells)"
  ),
  "7" = c(
    "CD16+ NK cells" = "CD16+ NK cells",
    "CD16- NK cells" = "Unresolved (cluster 7, CellTypist label: CD16- NK cells)",
    "Trm cytotoxic T cells" = "Unresolved (cluster 7, CellTypist label: Trm cytotoxic T cells)",
    "Tem/Trm cytotoxic T cells" = "Tem/Trm cytotoxic T cells",
    "Tem/Temra cytotoxic T cells" = "Unresolved (cluster 7, CellTypist label: Tem/Temra cytotoxic T cells)",
    "Regulatory T cells" = "Regulatory T cells",
    "NK cells" = "Unresolved (cluster 7, CellTypist label: NK cells)",
    "gamma-delta T cells" = "Unresolved (cluster 7, CellTypist label: gamma-delta T cells)",
    "CRTAM+ gamma-delta T cells" = "CRTAM+ gamma-delta T cells",
    "ILC3" = "Unresolved (cluster 7, CellTypist label: ILC3)",
    "MAIT cells" = "MAIT cells"
  ),
  "8" = c(
    "Classical monocytes" = "Classical monocytes",
    "Intermediate macrophages" = "Unresolved (cluster 8, CellTypist label: Intermediate macrophages)",
    "Macrophages" = "Macrophages",
    "Non-classical monocytes" = "Unresolved (cluster 8, CellTypist label: Non-classical monocytes)",
    "DC2" = "Unresolved (cluster 8, CellTypist label: DC2)",
    "Alveolar macrophages" = "Unresolved (cluster 8, CellTypist label: Alveolar macrophages)",
    "Regulatory T cells" = "Unresolved (cluster 8, CellTypist label: Regulatory T cells)"
  ),
  "9" = c(
    "Plasma cells" = "Plasma cells",
    "Memory B cells" = "Unresolved (cluster 9, CellTypist label: Memory B cells)"
  ),
  "10" = c(
    "DC2" = "DC2",
    "Migratory DCs" = "Migratory DCs",
    "DC" = "Unresolved (cluster 10, CellTypist label: DC)",
    "Macrophages" = "Macrophages",
    "Intermediate macrophages" = "Unresolved (cluster 10, CellTypist label: Intermediate macrophages)",
    "DC1" = "DC1",
    "Alveolar macrophages" = "Unresolved (cluster 10, CellTypist label: Alveolar macrophages)",
    "Classical monocytes" = "Classical monocytes"
  ),
  "11" = c(
    "Tem/Trm cytotoxic T cells" = "Tem/Trm cytotoxic T cells",
    "Regulatory T cells" = "Regulatory T cells",
    "Plasmablasts" = "Unresolved (cluster 11, CellTypist label: Plasmablasts)",
    "Memory B cells" = "Memory B cells",
    "Macrophages" = "Macrophages",
    "NK cells" = "NK cells",
    "CD16+ NK cells" = "CD16+ NK cells",
    "CD16- NK cells" = "Unresolved (cluster 11, CellTypist label: CD16- NK cells)",
    "Proliferative germinal center B cells" = "Unresolved (cluster 11, CellTypist label: Proliferative germinal center B cells)",
    "Age-associated B cells" = "Unresolved (cluster 11, CellTypist label: Age-associated B cells)",
    "CRTAM+ gamma-delta T cells" = "CRTAM+ gamma-delta T cells",
    "Germinal center B cells" = "Unresolved (cluster 11, CellTypist label: Germinal center B cells)",
    "gamma-delta T cells" = "gamma-delta T cells",
    "Epithelial cells" = "Other/contamination (cluster 11, CellTypist label: Epithelial cells)"
  ),
  "12" = c(
    "Low-density basophils" = "Unresolved (cluster 12, CellTypist label: Low-density basophils)",
    "Myeloid dendritic cells" = "Unresolved (cluster 12, CellTypist label: Myeloid dendritic cells)",
    "Progenitor cells" = "Unresolved (cluster 12, CellTypist label: Progenitor cells)"
  ),
  "13" = c(
    "Macrophages" = "Macrophages",
    "DC" = "Unresolved (cluster 13, CellTypist label: DC)",
    "DC2" = "DC2",
    "Intermediate macrophages" = "Unresolved (cluster 13, CellTypist label: Intermediate macrophages)",
    "Alveolar macrophages" = "Unresolved (cluster 13, CellTypist label: Alveolar macrophages)",
    "Plasmablasts" = "Unresolved (cluster 13, CellTypist label: Plasmablasts)",
    "DC1" = "DC1"
  ),
  "14" = c(
    "ILC3" = "Unresolved (cluster 14, CellTypist label: ILC3)",
    "Tem/Effector helper T cells" = "Unresolved (cluster 14, CellTypist label: Tem/Effector helper T cells)",
    "gamma-delta T cells" = "Unresolved (cluster 14, CellTypist label: gamma-delta T cells)"
  ),
  "15" = c(
    "Tem/Trm cytotoxic T cells" = "Tem/Trm cytotoxic T cells",
    "Regulatory T cells" = "Unresolved (cluster 15, CellTypist label: Regulatory T cells)",
    "Tem/Effector helper T cells" = "Unresolved (cluster 15, CellTypist label: Tem/Effector helper T cells)"
  ),
  "16" = c(
    "pDC" = "Unresolved (cluster 16, CellTypist label: pDC)",
    "DC2" = "Unresolved (cluster 16, CellTypist label: DC2)"
  ),
  "17" = c(
    "Germinal center B cells" = "Unresolved (cluster 17, CellTypist label: Germinal center B cells)",
    "Naive B cells" = "Unresolved (cluster 17, CellTypist label: Naive B cells)"
  ),
  "18" = c(
    "Plasma cells" = "Unresolved (cluster 18, CellTypist label: Plasma cells)",
    "Plasmablasts" = "Plasmablasts"
  ),
  "19" = c(
    "Non-classical monocytes" = "Non-classical monocytes",
    "Classical monocytes" = "Classical monocytes"
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

saveRDS(seu, here("Objects", "05_LUAD_seu_annotated.rds"))

final_summary <- as.data.frame(table(seu$annotation_final))
colnames(final_summary) <- c("annotation_final", "n_cells")
fwrite(final_summary, file.path(OUT_TAB, "05_LUAD_final_annotation.csv"))

message("Saved: 05_LUAD_seu_annotated.rds, 05_LUAD_final_annotation.csv")
message("DONE Script 05 — LUAD final annotation (CellTypist)")