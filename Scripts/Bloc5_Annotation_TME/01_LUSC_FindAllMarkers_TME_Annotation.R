#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc5_Annotation_TME Script 01: FindAllMarkers (LUSC)
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(here)
  library(presto)
})

IN_OBJ   <- here("Objects", "01_LUSC_seu_clustered.rds")
OUT_DIR  <- here("Results", "Bloc5_Annotation_TME", "Tables")

message("Loading clustered Seurat object...")
seu <- readRDS(IN_OBJ)

# Remove scale.data: dense matrix, only needed for RunPCA() at Bloc2,
# unused downstream (FindAllMarkers works on "data", not "scale.data").
# Removing it avoids unnecessary RAM load on this large object.
seu[["RNA"]]$scale.data <- NULL
gc()

if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
  stop("seurat_clusters absent. Run FindNeighbors + FindClusters first.")
}
Idents(seu) <- "seurat_clusters"
DefaultAssay(seu) <- "RNA"

seu <- JoinLayers(seu)
gc()
message("Layers joined.")
message("Cells: ", ncol(seu), 
        " | Genes: ", nrow(seu), 
        " | Clusters: ", length(levels(Idents(seu))))

message("Running FindAllMarkers (this may take a while)...")
markers <- FindAllMarkers(
  object          = seu,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  test.use        = "wilcox"
)

if (!"gene" %in% colnames(markers)) {
  markers <- markers %>% tibble::rownames_to_column("gene")
}
message("Total markers found: ", nrow(markers))

write_tsv(markers, file.path(OUT_DIR, "01_LUSC_markers_all_clusters.tsv"))
message("Saved: 01_LUSC_markers_all_clusters.tsv")

top50 <- markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 50) %>%
  ungroup()
write_tsv(top50, file.path(OUT_DIR, "01_LUSC_markers_top50_per_cluster.tsv"))
message("Saved: 01_LUSC_markers_top50_per_cluster.tsv")

message("DONE Bloc5 Script 01_LUSC")