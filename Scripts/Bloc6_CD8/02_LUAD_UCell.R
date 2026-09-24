# ============================================================
# Bloc6_CD8 — 02_LUAD_UCell.R
# Purpose: remove cycling and IFN-high (ISG+) cells from the pure
# CD8T population isolated by scGate (Script 01), using the same
# UCell signatures and thresholds as the official CD8 human TIL
# atlas construction pipeline (carmonalab/CD8_human_TIL_atlas_construction,
# script 1-QC_CD8_atlas.Rmd, step 5).
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(SignatuR)
  library(UCell)
  library(dplyr)
  library(readr)
})

data(SignatuR)

IN_OBJ  <- here("Objects", "01_LUAD_seu_scGate.rds")
OUT_OBJ <- here("Objects", "02_LUAD_seu_ucell.rds")
OUT_TAB <- here("Results", "Bloc6_CD8", "Tables")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

# CD8T cells kept by scGate in Script 01 (guard against a wrong input file)
EXPECTED_CD8T <- 56347     

message("Loading LUAD scGate object... ", Sys.time())
seu <- readRDS(IN_OBJ)
message("Cells loaded: ", ncol(seu))

stopifnot("scGate_multi" %in% colnames(seu@meta.data))

# ---- Subset to pure CD8T cells (scGate): base population for this step ----
# %in% drops the NA labels (cells not classified by scGate) explicitly.
message("Subsetting to scGate_multi == 'CD8T'...")
cd8_cells <- colnames(seu)[seu$scGate_multi %in% "CD8T"]
stopifnot(length(cd8_cells) == EXPECTED_CD8T)
seu_cd8 <- subset(seu, cells = cd8_cells)
rm(seu); gc()
message("CD8T cells: ", ncol(seu_cd8))

# ---- Lighten object before UCell: RNA assay only ----
# scale.data and reductions are not needed by UCell (rank-based scoring).
seu_cd8 <- DietSeurat(seu_cd8, assays = "RNA", dimreducs = NULL)

# ---- Signatures: identical to the official atlas construction pipeline ----
# cellCycle.G1S / cellCycle.G2M from SignatuR (carmonalab); IFN signature
# replaced with the authors' manual 4-gene set, documented as performing
# better with UCell than the full SignatuR IFN signature.
filter.signatures <- GetSignature(SignatuR$Hs)
filter.signatures <- filter.signatures[c("cellCycle.G1S", "cellCycle.G2M", "IFN")]
filter.signatures[["IFN"]] <- c("ISG15", "IFI6", "IFI44L", "MX1")

# ---- Compute UCell scores ----
message("Computing UCell scores... ", Sys.time())
seu_cd8 <- AddModuleScore_UCell(
  seu_cd8,
  features = filter.signatures,
  ncores   = 4,
  assay    = "RNA",
  name     = ""
)
message("UCell scores done. ", Sys.time())

# ---- Distribution check before filtering ----
message("Summary of UCell scores (pre-filtering):")
print(summary(seu_cd8@meta.data[, c("cellCycle.G1S", "cellCycle.G2M", "IFN")]))

# ---- Apply official thresholds ----
# Keep: cellCycle.G1S < 0.1 & cellCycle.G2M < 0.1 & IFN < 0.25
n_before <- ncol(seu_cd8)

keep <- seu_cd8$cellCycle.G1S < 0.1 &
  seu_cd8$cellCycle.G2M < 0.1 &
  seu_cd8$IFN < 0.25
stopifnot(!anyNA(keep))
seu_cd8_filtered <- subset(seu_cd8, cells = colnames(seu_cd8)[keep])

n_after   <- ncol(seu_cd8_filtered)
n_removed <- n_before - n_after

message("Cells before UCell filtering: ", n_before)
message("Cells after UCell filtering: ", n_after)
message("Cells removed (cycling/IFN-high): ", n_removed,
        " (", round(100 * n_removed / n_before, 1), "%)")

# ---- Breakdown of removal reason (non-exclusive, for reporting) ----
removal_summary <- data.frame(
  criterion = c("cellCycle.G1S >= 0.1", "cellCycle.G2M >= 0.1", "IFN >= 0.25"),
  n_cells_flagged = c(sum(seu_cd8$cellCycle.G1S >= 0.1),
                      sum(seu_cd8$cellCycle.G2M >= 0.1),
                      sum(seu_cd8$IFN >= 0.25))
)
print(removal_summary)
write_csv(removal_summary,
          file.path(OUT_TAB, "02_LUAD_UCell_removal_breakdown.csv"))

# ---- Cells and patients per response group, before and after ----
# Needed for the patient-level tests of the following scripts, and to
# check that the removal does not differ between MPR and non-MPR.
by_group <- function(obj, step) {
  obj@meta.data %>%
    group_by(pathological_response) %>%
    summarise(n_cells = n(), n_patients = n_distinct(sampleID), .groups = "drop") %>%
    mutate(step = step, .before = 1)
}
counts_by_group <- bind_rows(by_group(seu_cd8, "scGate_CD8T"),
                             by_group(seu_cd8_filtered, "post_UCell_filter"))
print(counts_by_group)
write_csv(counts_by_group,
          file.path(OUT_TAB, "02_LUAD_UCell_counts_by_response.csv"))

write_csv(
  data.frame(step = c("scGate_CD8T", "post_UCell_filter"),
             n_cells = c(n_before, n_after)),
  file.path(OUT_TAB, "02_LUAD_UCell_cell_counts.csv")
)

# ---- Save filtered CD8 object ----
saveRDS(seu_cd8_filtered, OUT_OBJ)
message("Script 02 (UCell, LUAD) complete. ", Sys.time())