#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc0_DataPrep, Script 01: Load MTX matrix + filter LUAD cells
# Uses BPCells to avoid RAM overload. Filters using the reference
# ID list produced by Script 00, not by re-filtering cancer_type.
# Input:  Data/ (MTX + barcodes + genes)
#         Objects/00_LUAD_reference_ids.rds
# Output: Objects/01_LUAD_matrix_BPCells/
# ============================================================
suppressPackageStartupMessages({
  library(here)
  library(BPCells)
  library(data.table)
})

# -----------------------------
# Paths
# -----------------------------
MTX_FILE      <- here("Data", "GSE243013_NSCLC_immune_scRNA_counts.mtx.gz")
BARCODES_FILE <- here("Data", "GSE243013_barcodes.csv.gz")
GENES_FILE    <- here("Data", "GSE243013_genes.csv.gz")
IDS_FILE      <- here("Objects", "00_LUAD_reference_ids.rds")
OUT_DIR       <- here("Objects", "01_LUAD_matrix_BPCells")

# -----------------------------
# 1) Load reference IDs from Script 00
# -----------------------------
message("Loading LUAD reference IDs...")
luad_ids <- readRDS(IDS_FILE)
message("Expected LUAD cells: ", length(luad_ids$cellID))

# -----------------------------
# 2) Open MTX matrix with BPCells
# -----------------------------
message("Opening MTX matrix with BPCells...")
barcode_names <- fread(BARCODES_FILE, header = TRUE)$barcode
gene_names    <- fread(GENES_FILE, header = TRUE)$geneSymbol

mat <- import_matrix_market(
  mtx_path  = MTX_FILE,
  row_names = barcode_names,
  col_names = gene_names
)
message("Full matrix dim (cells x genes): ", paste(dim(mat), collapse = " x "))

# -----------------------------
# 3) Transpose to genes x cells, convert type
# -----------------------------
mat <- t(mat)
mat <- convert_matrix_type(mat, type = "uint32_t")
message("Transposed matrix dim (genes x cells): ", paste(dim(mat), collapse = " x "))

# -----------------------------
# 4) Filter to LUAD cells using Script 00 reference list
# -----------------------------
message("Filtering to LUAD cells...")
mat_luad <- mat[, luad_ids$cellID]
message("LUAD matrix dim (genes x cells): ", paste(dim(mat_luad), collapse = " x "))

# -----------------------------
# 5) Verify filtered cell count matches Script 00
# -----------------------------
if (ncol(mat_luad) != length(luad_ids$cellID)) {
  stop("Filtered matrix cell count (", ncol(mat_luad), ") does not match ",
       "Script 00 reference count (", length(luad_ids$cellID), "). ",
       "Check for barcode mismatch between metadata and MTX files.")
}
message("Filtered cell count matches Script 00 reference. OK.")

# -----------------------------
# 6) Save LUAD matrix to disk (BPCells format)
# -----------------------------
message("Saving LUAD matrix to disk...")

# Remove existing output directory if present: write_matrix_dir() refuses
# to overwrite an existing folder (BPCells limitation), so this ensures
# the script can be safely re-run without manual cleanup.
unlink(OUT_DIR, recursive = TRUE, force = TRUE)

write_matrix_dir(mat_luad, OUT_DIR)
message("Saved: Objects/01_LUAD_matrix_BPCells/")
message("DONE Script 01 — LUAD matrix loading")