#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc2_Norm_HVG_PCA, Script 01: Normalization + HVG + PCA (LUSC)
# ============================================================
suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(ggplot2)
})

IN_OBJ  <- here("Objects", "02_LUSC_seu_qc.rds")
OUT_OBJ <- here("Objects")
OUT_FIG <- here("Results", "Bloc2_Norm_HVG_PCA", "Figures")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)

message("Loading QC Seurat object...")
seu <- readRDS(IN_OBJ)
DefaultAssay(seu) <- "RNA"
message("Cells: ", ncol(seu), " | Genes: ", nrow(seu))

message("Normalizing...")
seu <- NormalizeData(seu, normalization.method = "LogNormalize",
                     scale.factor = 1e4, verbose = FALSE)

# nfeatures = 3000 (vs Seurat default 2000): dataset is CD45+ sorted
# (confirmed via major_cell_type: B cell, Myeloid cell, T/NK cell only,
# no epithelial population). Increased HVG count chosen to better capture
# fine-grained heterogeneity within immune subpopulations (e.g. TAM
# subtypes, CD8 states), which a default threshold tuned for general
# tissue diversity could under-represent.
message("Selecting highly variable genes (nfeatures = 3000)...")
seu <- FindVariableFeatures(seu, selection.method = "vst",
                            nfeatures = 3000, verbose = FALSE)
message("HVGs selected: ", length(VariableFeatures(seu)))

# Cell cycle scoring, to regress out proliferation signal before PCA.
# Added after Bloc5 annotation revealed a cluster (cluster 11) mixing
# multiple cell lineages, driven by cell-cycle genes rather than lineage
# identity (see CAHIER_METHODOLOGIQUE.md, Bloc5 Script 01).
message("Scoring cell cycle phase...")
s.genes   <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
seu <- CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes,
                        set.ident = FALSE)

message("Scaling HVGs...")
seu <- ScaleData(seu, features = VariableFeatures(seu),
                 vars.to.regress = c("S.Score", "G2M.Score"),
                 verbose = FALSE)

message("Running PCA...")
seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 50, verbose = FALSE)

p_elbow <- ElbowPlot(seu, ndims = 50) + theme_bw() +
  labs(title = "PCA Elbow Plot — GSE243013 LUSC",
       x = "Principal Component", y = "Standard Deviation") +
  theme(plot.title = element_text(size = 12, face = "bold"))

print(p_elbow)

ggsave(file.path(OUT_FIG, "01_LUSC_ElbowPlot.png"), p_elbow,
       width = 6, height = 4, dpi = 300, bg = "white")

saveRDS(seu, file.path(OUT_OBJ, "01_LUSC_seu_pca.rds"))
message("DONE Script 01 — LUSC Normalization/HVG/PCA")