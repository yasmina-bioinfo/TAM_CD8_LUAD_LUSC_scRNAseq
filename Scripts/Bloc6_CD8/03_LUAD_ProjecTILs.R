# ============================================================
# Bloc6_CD8 — 03_LUAD_ProjecTILs.R
# Purpose: project the purified, cycling/ISG-filtered LUAD CD8 population
# (Script 02 output) onto the official human CD8 TIL reference atlas
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

IN_OBJ    <- here("Objects", "02_LUAD_seu_ucell.rds")
OUT_OBJ   <- here("Objects", "03_LUAD_seu_projectils.rds")
OUT_TAB   <- here("Results", "Bloc6_CD8", "Tables")
REF_DIR   <- here("References")
REF_FILE  <- file.path(REF_DIR, "CD8T_human_ref_v1.rds")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(REF_DIR, recursive = TRUE, showWarnings = FALSE)

# Cells kept after the UCell filter (Script 02) and patients per group
EXPECTED_CELLS    <- 47996                        
EXPECTED_PATIENTS <- c("MPR" = 10, "non-MPR" = 40)

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

# ---- Load query (LUAD, post scGate + UCell) ----
message("Loading LUAD query object...")
seu <- readRDS(IN_OBJ)
message("Query cells: ", ncol(seu))
stopifnot("sampleID" %in% colnames(seu@meta.data))

# ---- Guard: input must be the object validated at Script 02 ----
stopifnot(ncol(seu) == EXPECTED_CELLS)
obs_pat <- seu@meta.data %>%
  distinct(sampleID, pathological_response) %>%
  count(pathological_response)
print(obs_pat)
obs_vec <- setNames(obs_pat$n, obs_pat$pathological_response)
stopifnot(all(obs_vec[names(EXPECTED_PATIENTS)] == EXPECTED_PATIENTS))

# ---- Traceability: package versions and reference atlas fingerprint ----
writeLines(c(
  "Script : 03_LUAD_ProjecTILs.R",
  paste("Execution date :", Sys.time()),
  paste("ProjecTILs :", as.character(packageVersion("ProjecTILs"))),
  paste("STACAS :", as.character(packageVersion("STACAS"))),
  paste("Seurat :", as.character(packageVersion("Seurat"))),
  paste("Reference atlas file :", basename(REF_FILE),
        "| md5 :", unname(tools::md5sum(REF_FILE))),
  "Reference source : figshare file 41414556 (human CD8 TIL atlas)"
), file.path(OUT_TAB, "03_LUAD_ProjecTILs_run_info.txt"))

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

# ---- No cell should be lost (filter.cells = FALSE) ----
message("Cells projected: ", ncol(query.projected), " of ", ncol(seu))
stopifnot(ncol(query.projected) == ncol(seu))

# ---- Distribution of predicted functional clusters ----
message("functional.cluster distribution:")
print(table(query.projected$functional.cluster, useNA = "ifany"))

message("functional.cluster.conf summary:")
print(summary(query.projected$functional.cluster.conf))

# ---- Cells and patients per functional state and response group ----
state_tab <- query.projected@meta.data %>%
  group_by(functional.cluster, pathological_response) %>%
  summarise(n_cells = n(), n_patients = n_distinct(sampleID), .groups = "drop")
print(state_tab, n = Inf)
write_csv(state_tab,
          file.path(OUT_TAB, "03_LUAD_ProjecTILs_states_by_response.csv"))

# ---- Cross-check against CellTypist annotation_final (CD8 substates) ----
crosscheck <- table(
  functional.cluster = query.projected$functional.cluster,
  annotation_final    = query.projected$annotation_final,
  useNA = "ifany"
)
write.csv(
  as.data.frame.matrix(crosscheck),
  file.path(OUT_TAB, "03_LUAD_ProjecTILs_vs_annotationFinal_crosscheck.csv")
)

# ---- Save ----
saveRDS(query.projected, OUT_OBJ)
message("Script 03 (ProjecTILs, LUAD) complete. ", Sys.time())