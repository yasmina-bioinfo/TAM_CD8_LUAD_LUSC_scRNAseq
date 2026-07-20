#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc3_Harmony, Script 01: Batch correction (LUSC)
# ============================================================
suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(harmony)
})

IN_OBJ  <- here("Objects", "01_LUSC_seu_pca.rds")
OUT_OBJ <- here("Objects")

message("Loading PCA Seurat object...")
seu <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu), " | Genes: ", nrow(seu))

# dims.use = 1:20: number of PCs determined and justified in Bloc2
# (variance explained per PC, elbow ~PC17-18, margin retained for rare
# subpopulations)
message("Running Harmony batch correction...")
seu <- RunHarmony(seu,
                  group.by.vars = "sampleID",
                  dims.use      = 1:20,
                  project.dim   = FALSE,
                  verbose       = FALSE)

saveRDS(seu, file.path(OUT_OBJ, "01_LUSC_seu_harmony.rds"))
message("DONE Script 01 — LUSC Harmony")