# ============================================================
# Bloc6_CD8 — Script01_scGate_LUSC.R
# Purpose: purify CD8 T cells from the full LUSC TME object using
# scGate (marker-based, hierarchical gating), independent of any
# prior annotation (CellTypist / annotation_final / lineage_CD8).
# Reference: Andreatta & Carmona, Bioinformatics, 2022.
# Model: scGate_models_DB$human$TME_HiRes
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(scGate)
  library(dplyr)
  library(readr)
})

IN_OBJ  <- here("Objects", "06_LUSC_seu_final.rds")
OUT_OBJ <- here("Objects", "01_LUSC_seu_scGate.rds")
OUT_TAB <- here("Results", "Bloc6_CD8", "Tables")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

message("Loading LUSC final object... ", Sys.time())
seu <- readRDS(IN_OBJ)
message("Cells loaded: ", ncol(seu))

stopifnot(all(c("annotation_final", "lineage_CD8") %in% colnames(seu@meta.data)))
stopifnot("harmony" %in% names(seu@reductions))

message("Loading scGate TME_HiRes model database...")
models.TME <- get_scGateDB()

# ---- Lighten object before scGate: RNA assay only, no scale.data,
#      but KEEP the harmony reduction to avoid PCA recomputation at
#      every hierarchical level of the model ----
message("Creating lightweight object for scGate (RNA + harmony only)...")
seu_light <- DietSeurat(seu, assays = "RNA", dimreducs = "harmony")

message("Running scGate (reduction = 'harmony', ncores = 4)... ", Sys.time())
seu_light <- scGate(
  seu_light,
  model     = models.TME$human$TME_HiRes,
  reduction = "harmony",
  ncores    = 4
)
message("scGate done. ", Sys.time())

# ---- Transfer result back to full object ----
seu$scGate_multi <- seu_light$scGate_multi
rm(seu_light); gc()

message("scGate_multi distribution:")
print(table(seu$scGate_multi, useNA = "ifany"))

crosscheck <- table(
  scGate_CD8T  = seu$scGate_multi == "CD8T",
  lineage_CD8  = seu$lineage_CD8
)
print(crosscheck)

write_csv(
  as.data.frame.matrix(crosscheck) %>% tibble::rownames_to_column("scGate_CD8T"),
  file.path(OUT_TAB, "01_LUSC_scGate_vs_lineageCD8_crosscheck.csv")
)

n_scGate_cd8   <- sum(seu$scGate_multi == "CD8T", na.rm = TRUE)
n_lineage_cd8  <- sum(seu$lineage_CD8, na.rm = TRUE)
n_intersection <- sum(seu$scGate_multi == "CD8T" & seu$lineage_CD8, na.rm = TRUE)

message("scGate CD8T: ", n_scGate_cd8)
message("lineage_CD8 (CellTypist): ", n_lineage_cd8)
message("Intersection: ", n_intersection)

saveRDS(seu, OUT_OBJ)
message("Script01 (scGate, LUSC) complete. ", Sys.time())