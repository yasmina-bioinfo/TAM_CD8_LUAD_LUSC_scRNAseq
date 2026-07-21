#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc5_Annotation_TME, Script 02: SingleR annotation (LUAD)
# Reference: MonacoImmuneData (immune-only, appropriate for this
# CD45+ sorted dataset). label.fine is the primary cross-check tool
# against canonical marker annotation.
# ============================================================
suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(SingleR)
  library(celldex)
  library(ggplot2)
  library(dplyr)
  library(data.table)
})

# -----------------------------
# Paths
# -----------------------------
IN_OBJ  <- here("Objects", "01_LUAD_seu_clustered.rds")
OUT_FIG <- here("Results", "Bloc5_Annotation_TME", "Figures")
OUT_TAB <- here("Results", "Bloc5_Annotation_TME", "Tables")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_TAB, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1) Load clustered object
# -----------------------------
message("Loading clustered Seurat object...")
seu <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu), " | Clusters: ", length(levels(seu$seurat_clusters)))

# -----------------------------
# 2) Load reference: MonacoImmuneData
# -----------------------------
# Immune-specific reference, appropriate for this CD45+ sorted dataset
# (no epithelial/stromal compartment to justify a broader reference like
# HumanPrimaryCellAtlas). See CAHIER_METHODOLOGIQUE.md for full rationale.
message("Loading MonacoImmuneData reference...")
ref <- MonacoImmuneData()

# -----------------------------
# 3) Extract normalized expression matrix
# -----------------------------
# NOTE: JoinLayers() required for Seurat v5 multi-layer objects after merge()
# Object is already a standard dgCMatrix (converted in Bloc1 Script02),
# no BPCells-specific conversion needed here.
seu <- JoinLayers(seu)
mat <- GetAssayData(seu, layer = "data")

# -----------------------------
# 4) Run SingleR: label.main (broad cell types)
# -----------------------------
message("Running SingleR (label.main)...")
pred_main <- SingleR(
  test    = mat,
  ref     = ref,
  labels  = ref$label.main,
  BPPARAM = BiocParallel::SnowParam(1)
)
message("SingleR main annotation done.")

# -----------------------------
# 5) Run SingleR: label.fine (fine-grained subtypes)
# -----------------------------
# label.fine is the primary cross-check tool against canonical marker
# annotation, per project methodology (see ROADMAP.md)
message("Running SingleR (label.fine)...")
pred_fine <- SingleR(
  test    = mat,
  ref     = ref,
  labels  = ref$label.fine,
  BPPARAM = BiocParallel::SnowParam(1)
)
message("SingleR fine annotation done.")

# -----------------------------
# 6) Add annotations to Seurat object
# -----------------------------
seu$SingleR_main <- pred_main$labels
seu$SingleR_fine <- pred_fine$labels

message("Generating summary tables...")

summary_main <- as.data.frame(table(
  Cluster = seu$seurat_clusters,
  SingleR = seu$SingleR_main
))
summary_main <- summary_main[summary_main$Freq > 0, ]
fwrite(summary_main, file.path(OUT_TAB, "02_LUAD_SingleR_main_per_cluster.csv"))

summary_fine <- as.data.frame(table(
  Cluster = seu$seurat_clusters,
  SingleR = seu$SingleR_fine
))
summary_fine <- summary_fine[summary_fine$Freq > 0, ]
fwrite(summary_fine, file.path(OUT_TAB, "02_LUAD_SingleR_fine_per_cluster.csv"))

message("Saved: SingleR summary tables")

# -----------------------------
# 7) Heatmaps
# -----------------------------
message("Generating heatmaps...")

png(file.path(OUT_FIG, "02_LUAD_SingleR_heatmap_main.png"),
    width = 12, height = 8, units = "in", res = 300)
plotScoreHeatmap(pred_main)
dev.off()

png(file.path(OUT_FIG, "02_LUAD_SingleR_heatmap_fine.png"),
    width = 16, height = 10, units = "in", res = 300)
plotScoreHeatmap(pred_fine)
dev.off()
message("Saved: heatmaps")

# -----------------------------
# 8) UMAP figures
# -----------------------------
message("Generating UMAP figures...")

png(file.path(OUT_FIG, "02_LUAD_UMAP_SingleR_main.png"),
    width = 14, height = 8, units = "in", res = 300)
print(DimPlot(seu, reduction = "umap", group.by = "SingleR_main",
              label = FALSE, raster = FALSE) +
        theme_bw() + labs(title = "SingleR annotation (main) — GSE243013 LUAD"))
dev.off()

png(file.path(OUT_FIG, "02_LUAD_UMAP_SingleR_main_labeled.png"),
    width = 14, height = 8, units = "in", res = 300)
print(DimPlot(seu, reduction = "umap", group.by = "SingleR_main",
              label = FALSE, raster = FALSE) +
        geom_text(data = data.frame(
          Embeddings(seu, "umap"),
          cluster = seu$seurat_clusters) %>%
            group_by(cluster) %>%
            summarise(umap_1 = median(umap_1), umap_2 = median(umap_2)),
          aes(x = umap_1, y = umap_2, label = cluster),
          size = 3, fontface = "bold", color = "black") +
        theme_bw() +
        labs(title = "SingleR annotation (main) with cluster numbers — GSE243013 LUAD"))
dev.off()

png(file.path(OUT_FIG, "02_LUAD_UMAP_SingleR_fine.png"),
    width = 18, height = 10, units = "in", res = 300)
print(DimPlot(seu, reduction = "umap", group.by = "SingleR_fine",
              label = FALSE, raster = FALSE) +
        theme_bw() + labs(title = "SingleR annotation (fine) — GSE243013 LUAD") +
        theme(legend.text = element_text(size = 6)))
dev.off()
message("Saved: UMAP figures")

# -----------------------------
# 9) Save updated Seurat object
# -----------------------------
message("Saving updated Seurat object...")
saveRDS(seu, here("Objects", "02_LUAD_seu_singler.rds"))
message("DONE Script 02 — LUAD SingleR annotation")