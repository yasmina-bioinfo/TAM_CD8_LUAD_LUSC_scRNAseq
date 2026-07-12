#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc0_DataPrep, Script 01: Load MTX matrix + filter LUSC cells
# ============================================================
suppressPackageStartupMessages({
  library(here)
  library(BPCells)
  library(data.table)
})

MTX_FILE      <- here("Data", "GSE243013_NSCLC_immune_scRNA_counts.mtx.gz")
BARCODES_FILE <- here("Data", "GSE243013_barcodes.csv.gz")
GENES_FILE    <- here("Data", "GSE243013_genes.csv.gz")
IDS_FILE      <- here("Objects", "00_LUSC_reference_ids.rds")
OUT_DIR       <- here("Objects", "01_LUSC_matrix_BPCells")

message("Loading LUSC reference IDs...")
lusc_ids <- readRDS(IDS_FILE)
message("Expected LUSC cells: ", length(lusc_ids$cellID))

message("Opening MTX matrix with BPCells...")
barcode_names <- fread(BARCODES_FILE, header = TRUE)$barcode
gene_names    <- fread(GENES_FILE, header = TRUE)$geneSymbol

mat <- import_matrix_market(
  mtx_path  = MTX_FILE,
  row_names = barcode_names,
  col_names = gene_names
)
message("Full matrix dim (cells x genes): ", paste(dim(mat), collapse = " x "))

mat <- t(mat)
mat <- convert_matrix_type(mat, type = "uint32_t")
message("Transposed matrix dim (genes x cells): ", paste(dim(mat), collapse = " x "))

message("Filtering to LUSC cells...")
mat_lusc <- mat[, lusc_ids$cellID]
message("LUSC matrix dim (genes x cells): ", paste(dim(mat_lusc), collapse = " x "))

if (ncol(mat_lusc) != length(lusc_ids$cellID)) {
  stop("Filtered matrix cell count (", ncol(mat_lusc), ") does not match ",
       "Script 00 reference count (", length(lusc_ids$cellID), "). ",
       "Check for barcode mismatch between metadata and MTX files.")
}
message("Filtered cell count matches Script 00 reference. OK.")

message("Saving LUSC matrix to disk...")
write_matrix_dir(mat_lusc, OUT_DIR)
message("Saved: Objects/01_LUSC_matrix_BPCells/")
message("DONE Script 01 — LUSC matrix loading")