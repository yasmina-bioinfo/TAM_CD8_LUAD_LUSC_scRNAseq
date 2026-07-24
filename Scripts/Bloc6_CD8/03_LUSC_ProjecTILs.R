# ============================================================
# Bloc6_CD8 — Script03_ProjecTILs_LUSC.R
# Purpose: project the purified, cycling/ISG-filtered LUSC CD8 population
# (Script02 output) onto the official human CD8 TIL reference atlas
# (Andreatta et al., 2021, Nat Commun; atlas construction: Carmona lab,
# CD8_human_TIL_atlas_construction), to obtain fine-grained functional
# cluster labels and cross-validate against CellTypist/scGate.
# ============================================================

# ---- Installation (run once) ----
# if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
# remotes::install_github("carmonalab/STACAS")
# remotes::install_github("carmonalab/ProjecTILs")

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(ProjecTILs)
  library(dplyr)
  library(readr)
})

IN_OBJ    <- here("Objects", "02_LUSC_seu_ucell.rds")
OUT_OBJ   <- here("Objects", "03_LUSC_seu_projectils.rds")
OUT_TAB   <- here("Results", "Bloc6_CD8", "Tables")
REF_DIR   <- here("References")
REF_FILE  <- file.path(REF_DIR, "CD8T_human_ref_v1.rds")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(REF_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Download reference atlas (once; cached locally afterward) ----
if (!file.exists(REF_FILE)) {
  message("Downloading human CD8 TIL reference atlas...")
  options(timeout = max(900, getOption("timeout")))
  download.file(
    "https://figshare.com/ndownloader/files/41414556",
    destfile = REF_FILE, mode = "wb"
  )
}

message("Loading reference atlas... ", Sys.time())
ref.cd8 <- load.reference.map(REF_FILE)

# ---- Load query (LUSC, post scGate + UCell) ----
message("Loading LUSC query object...")
seu <- readRDS(IN_OBJ)
message("Query cells: ", ncol(seu))
stopifnot("sampleID" %in% colnames(seu@meta.data))

# ---- Run ProjecTILs ----
# filter.cells = FALSE: query already purified (scGate TME_HiRes CD8T +
# UCell cycling/ISG+ removal), per official guidance for pre-filtered data.
# split.by = "sampleID": project each patient separately for optimal
# batch-effect correction, as recommended in the official tutorial.
message("Running Run.ProjecTILs... ", Sys.time())
query.projected <- Run.ProjecTILs(
  seu,
  ref = ref.cd8,
  filter.cells = FALSE,
  split.by = "sampleID",
  ncores = 4
)
message("Projection done. ", Sys.time())

# ---- Distribution of predicted functional clusters ----
message("functional.cluster distribution:")
print(table(query.projected$functional.cluster, useNA = "ifany"))

message("functional.cluster.conf summary:")
print(summary(query.projected$functional.cluster.conf))

# ---- Cross-check against CellTypist annotation_final (CD8 substates) ----
crosscheck <- table(
  functional.cluster = query.projected$functional.cluster,
  annotation_final    = query.projected$annotation_final,
  useNA = "ifany"
)
write.csv(
  as.data.frame.matrix(crosscheck),
  file.path(OUT_TAB, "03_LUSC_ProjecTILs_vs_annotationFinal_crosscheck.csv")
)

# ---- Save ----
saveRDS(query.projected, OUT_OBJ)
message("Script03 (ProjecTILs, LUSC) complete. ", Sys.time())