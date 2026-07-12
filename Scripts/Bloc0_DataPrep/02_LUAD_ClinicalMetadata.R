#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc0_DataPrep, Script 02: Clinical metadata (LUAD)
# ============================================================
suppressPackageStartupMessages({
  library(here)
  library(data.table)
})

META_FILE  <- here("Data", "GSE243013_NSCLC_immune_scRNA_metadata.csv.gz")
IDS_FILE   <- here("Objects", "00_LUAD_reference_ids.rds")
TABLES_DIR <- here("Results", "Bloc0_DataPrep", "Tables")

luad_ids <- readRDS(IDS_FILE)
meta <- fread(META_FILE)

cols_patient <- c("sampleID", "pathological_response", "cancer_type",
                  "gender", "age", "anti-PD1_therapy", "radiological_response")
missing_cols <- setdiff(cols_patient, colnames(meta))
if (length(missing_cols) > 0) stop("Missing columns: ", paste(missing_cols, collapse = ", "))

meta_luad <- meta[sampleID %in% luad_ids$sampleID]

clinical_luad <- unique(meta_luad[, ..cols_patient], by = "sampleID")

if (nrow(clinical_luad) != length(luad_ids$sampleID)) {
  stop("Patient count mismatch with Script 00: got ", nrow(clinical_luad),
       ", expected ", length(luad_ids$sampleID))
}
message("Patient count matches Script 00. OK.")

# NOTE: 'unknowm' is a typo in the original dataset, kept as-is for traceability.
message("pathological_response: ", paste(unique(clinical_luad$pathological_response), collapse = ", "))
message("NA check:"); print(colSums(is.na(clinical_luad)))

dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
fwrite(clinical_luad, file.path(TABLES_DIR, "02_LUAD_clinical_metadata.csv"))
message("DONE Script 02 — LUAD clinical metadata")