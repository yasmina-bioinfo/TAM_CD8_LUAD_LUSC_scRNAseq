# ============================================================
# Bloc6_CD8 — Script02_UCell_LUSC.R
# Purpose: remove cycling and IFN-high (ISG+) cells from the pure
# CD8T population isolated by scGate (Script01), using the same
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

IN_OBJ  <- here("Objects", "01_LUSC_seu_scGate.rds")
OUT_OBJ <- here("Objects", "02_LUSC_seu_ucell.rds")
OUT_TAB <- here("Results", "Bloc6_CD8", "Tables")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

message("Loading LUSC scGate object... ", Sys.time())
seu <- readRDS(IN_OBJ)
message("Cells loaded: ", ncol(seu))

stopifnot("scGate_multi" %in% colnames(seu@meta.data))

# ---- Subset to pure CD8T cells (scGate), base population for this step ----
message("Subsetting to scGate_multi == 'CD8T'...")
seu_cd8 <- subset(seu, subset = scGate_multi == "CD8T")
message("CD8T cells: ", ncol(seu_cd8))

# ---- Lighten object before UCell: RNA assay only ----
# UCell does not perform dimensionality reduction/kNN smoothing (unlike
# scGate), so no reduction needs to be preserved here — only scale.data
# and reductions are dropped to avoid the same layer-conflict slowdown
# seen in Script01.
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
# cellCycle.G1S < 0.1 & cellCycle.G2M < 0.1 & IFN < 0.25
n_before <- ncol(seu_cd8)

seu_cd8_filtered <- subset(
  seu_cd8,
  subset = cellCycle.G1S < 0.1 & cellCycle.G2M < 0.1 & IFN < 0.25
)

n_after   <- ncol(seu_cd8_filtered)
n_removed <- n_before - n_after

message("Cells before UCell filtering: ", n_before)
message("Cells after UCell filtering: ", n_after)
message("Cells removed (cycling/IFN-high): ", n_removed,
        " (", round(100 * n_removed / n_before, 1), "%)")

# ---- Breakdown of removal reason (non-exclusive, for reporting) ----
removed_cycling_G1S <- sum(seu_cd8$cellCycle.G1S >= 0.1)
removed_cycling_G2M <- sum(seu_cd8$cellCycle.G2M >= 0.1)
removed_IFN         <- sum(seu_cd8$IFN >= 0.25)

removal_summary <- data.frame(
  criterion = c("cellCycle.G1S >= 0.1", "cellCycle.G2M >= 0.1", "IFN >= 0.25"),
  n_cells_flagged = c(removed_cycling_G1S, removed_cycling_G2M, removed_IFN)
)
print(removal_summary)

write_csv(
  removal_summary,
  file.path(OUT_TAB, "02_LUSC_UCell_removal_breakdown.csv")
)

write_csv(
  data.frame(
    step = c("scGate_CD8T", "post_UCell_filter"),
    n_cells = c(n_before, n_after)
  ),
  file.path(OUT_TAB, "02_LUSC_UCell_cell_counts.csv")
)

# ---- Save filtered CD8 object ----
saveRDS(seu_cd8_filtered, OUT_OBJ)
message("Script02 (UCell, LUSC) complete. ", Sys.time())