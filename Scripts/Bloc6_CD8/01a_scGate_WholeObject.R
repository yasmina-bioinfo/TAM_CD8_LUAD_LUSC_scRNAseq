# ============================================================
# Bloc6_CD8 — 01a_LUAD_scGate_WholeObject.R
# Purpose: run scGate on the WHOLE LUAD object (no batching).
# This run is NOT used for the results of the pipeline. It is the
# reference against which the batched run (01_LUAD_scGate.R) is
# compared in 01b_LUAD_scGate_BatchVerification.R.
# It is only feasible for LUAD: the whole LUSC object needs about
# 96 GB of RAM (32 GB available).
# Output: Objects/01_LUAD_seu_scGate_wholeobject.rds (read by 01b).
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

IN_OBJ  <- here("Objects", "06_LUAD_seu_final.rds")
OUT_OBJ <- here("Objects", "01_LUAD_seu_scGate_wholeobject.rds")
OUT_TAB <- here("Results", "Bloc6_CD8", "Tables")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

message("Loading LUAD final object... ", Sys.time())
seu <- readRDS(IN_OBJ)
message("Cells loaded: ", ncol(seu))

stopifnot(all(c("annotation_final", "lineage_CD8") %in% colnames(seu@meta.data)))
stopifnot("harmony" %in% names(seu@reductions))

# ---- Cohort guard. The object holds MPR and non-MPR patients only
#      (pCR already excluded): LUAD 50 patients (10 MPR + 40 non-MPR) ----
print(
  seu@meta.data %>%
    distinct(sampleID, pathological_response) %>%
    count(pathological_response)
)
print(table(seu$pathological_response))
stopifnot(n_distinct(seu$sampleID) == 50)

message("Loading scGate TME_HiRes model database...")
models.TME <- get_scGateDB()

writeLines(c(
  "Script : 01a_LUAD_scGate_WholeObject.R",
  paste("Execution date :", Sys.time()),
  paste("scGate :", as.character(packageVersion("scGate"))),
  paste("Seurat :", as.character(packageVersion("Seurat"))),
  "Model : TME_HiRes, downloaded via get_scGateDB()"
), file.path(OUT_TAB, "01a_LUAD_scGate_wholeobject_run_info.txt"))

# ---- Lighten object: RNA assay only, KEEP harmony ----
seu_light <- DietSeurat(seu, assays = "RNA", dimreducs = "harmony")

message("Running scGate on the whole object... ", Sys.time())
seu_light <- scGate(
  seu_light,
  model     = models.TME$human$TME_HiRes,
  reduction = "harmony",
  ncores    = 1 
)
message("scGate done. ", Sys.time())

seu$scGate_multi <- seu_light$scGate_multi
rm(seu_light); gc()

message("scGate_multi distribution:")
print(table(seu$scGate_multi, useNA = "ifany"))

# ---- Cross-check with CellTypist lineage (NA-safe) ----
scGate_CD8T <- seu$scGate_multi %in% "CD8T"
crosscheck <- table(
  scGate_CD8T = scGate_CD8T,
  lineage_CD8 = seu$lineage_CD8,
  useNA = "ifany"
)
print(crosscheck)
stopifnot(sum(crosscheck) == ncol(seu))

message("scGate CD8T: ", sum(scGate_CD8T))
message("lineage_CD8 (CellTypist): ", sum(seu$lineage_CD8, na.rm = TRUE))
message("Intersection: ", sum(scGate_CD8T & seu$lineage_CD8 %in% TRUE))

saveRDS(seu, OUT_OBJ)
message("Script 01a (scGate, LUAD, whole object) complete. ", Sys.time())