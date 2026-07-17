#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc5_Annotation_TME, Script 03: Resolve cluster 11 CD8 subtypes (LUAD)
# Cluster 11 (cycling, mixed lineage) contains CD8 T cells labeled by
# SingleR as three distinct states (effector memory, central memory,
# terminal effector). Rather than assuming these map onto the existing
# manual CD8 categories (clusters 2/3) or presupposing which markers to
# check (confirmation bias), this script lets FindAllMarkers reveal the
# actual differential markers between these three SingleR-defined
# subgroups, unbiased by prior assumption.
# Input:  Objects/03_LUAD_seu_annotated.rds
# Output: Results/Bloc5_Annotation_TME/Tables/03_LUAD_cluster11_CD8_markers.tsv
# ============================================================
suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(dplyr)
  library(readr)
})

IN_OBJ <- here("Objects", "02_LUAD_seu_singler.rds")
OUT_TAB <- here("Results", "Bloc5_Annotation_TME", "Tables")
dir.create(OUT_TAB, showWarnings = FALSE, recursive = TRUE)

message("Loading annotated Seurat object...")
seu <- readRDS(IN_OBJ)

# -----------------------------
# 1) Isolate CD8 cells within cluster 11, by SingleR_fine label
# -----------------------------
cd8_labels_c11 <- c("Effector memory CD8 T cells",
                    "Central memory CD8 T cells",
                    "Terminal effector CD8 T cells")

cd8_c11_cells <- colnames(seu)[seu$seurat_clusters == "11" &
                                 seu$SingleR_fine %in% cd8_labels_c11]
message("CD8 cells found in cluster 11: ", length(cd8_c11_cells))
print(table(seu$SingleR_fine[cd8_c11_cells]))

# -----------------------------
# 2) Subset and re-identify by SingleR label, not manual assumption
# -----------------------------
seu_cd8_c11 <- subset(seu, cells = cd8_c11_cells)
Idents(seu_cd8_c11) <- seu_cd8_c11$SingleR_fine

# -----------------------------
# 3) FindAllMarkers, no prior marker assumption
# -----------------------------
message("Running FindAllMarkers on cluster 11 CD8 subgroups...")
markers_c11_cd8 <- FindAllMarkers(seu_cd8_c11, only.pos = TRUE,
                                  min.pct = 0.25, logfc.threshold = 0.25)
message("Total markers found: ", nrow(markers_c11_cd8))

top10_c11_cd8 <- markers_c11_cd8 %>%
  group_by(cluster) %>%
  slice_head(n = 10) %>%
  select(cluster, gene, avg_log2FC) %>%
  ungroup()

print(top10_c11_cd8, n = 30)

write_tsv(markers_c11_cd8, file.path(OUT_TAB, "03_LUAD_cluster11_CD8_markers_all.tsv"))
write_tsv(top10_c11_cd8, file.path(OUT_TAB, "03_LUAD_cluster11_CD8_markers_top10.tsv"))
message("Saved: cluster11 CD8 markers tables")
message("DONE Script 03 — LUAD cluster 11 CD8 subtype resolution")