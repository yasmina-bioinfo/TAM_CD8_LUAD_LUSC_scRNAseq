#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc0_DataPrep, Script 00: LUAD/LUSC separation
# Single source of truth for the LUAD/LUSC split used across the
# entire pipeline. All downstream scripts (Bloc0-Bloc8) reference
# the reference ID lists produced here rather than re-filtering
# cancer_type independently.
# Input:  Data/GSE243013_NSCLC_immune_scRNA_metadata.csv.gz
# Output: Objects/00_LUAD_reference_ids.rds
#         Objects/00_LUSC_reference_ids.rds
#         Results/Bloc0_DataPrep/Tables/00_separation_summary.csv
# ============================================================
suppressPackageStartupMessages({
  library(here)
  library(data.table)
})

# -----------------------------
# Paths (relative, via here())
# -----------------------------
META_FILE   <- here("Data", "GSE243013_NSCLC_immune_scRNA_metadata.csv.gz")
OBJECTS_DIR <- here("Objects")
TABLES_DIR  <- here("Results", "Bloc0_DataPrep", "Tables")

# -----------------------------
# 1) Read metadata
# -----------------------------
message("Reading metadata...")
meta <- fread(META_FILE)

# -----------------------------
# 2) Verify expected columns exist before doing anything else
# -----------------------------
expected_cols <- c("sampleID", "cellID", "cancer_type", "pathological_response", "anti-PD1_therapy")
missing_cols  <- setdiff(expected_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop("Missing expected columns: ", paste(missing_cols, collapse = ", "))
}

# -----------------------------
# 3) Verify cancer_type only contains expected values
# -----------------------------
observed_types <- unique(meta$cancer_type)
unexpected_types <- setdiff(observed_types, c("LUAD", "LUSC"))
if (length(unexpected_types) > 0) {
  stop("Unexpected cancer_type value(s) found: ", paste(unexpected_types, collapse = ", "),
       ". Investigate before proceeding; do not silently drop or ignore.")
}

# -----------------------------
# 3b) Exclude chemotherapy-only patients (no anti-PD-1)
# -----------------------------
# Liu et al. 2025 (Cell): of 343 patients followed, 243 underwent
# scRNA/TCR-seq — 234 with anti-PD-1 + chemotherapy, 9 with
# chemotherapy alone. This study concerns anti-PD-1 response only.
message("Identifying chemotherapy-only patients (no anti-PD-1)...")

chemo_only_ids <- unique(meta[`anti-PD1_therapy` == "No", sampleID])

if (length(chemo_only_ids) != 9) {
  stop("Expected exactly 9 chemotherapy-only patients, found ", length(chemo_only_ids),
       ". Investigate before proceeding.")
}

# Log the excluded patients for traceability before dropping them
exclusion_log <- unique(meta[sampleID %in% chemo_only_ids,
                             .(sampleID, cancer_type, pathological_response)])
fwrite(exclusion_log, file.path(TABLES_DIR, "00_excluded_chemo_only_patients.csv"))
message("Saved: Results/Bloc0_DataPrep/Tables/00_excluded_chemo_only_patients.csv")

message("Excluding ", length(chemo_only_ids), " chemotherapy-only patients (2 LUAD, 7 LUSC).")
meta <- meta[!sampleID %in% chemo_only_ids]



# -----------------------------
# 4) Build reference ID lists per histology
# -----------------------------
message("Building reference ID lists...")

luad_ids <- list(
  cellID   = meta[cancer_type == "LUAD", cellID],
  sampleID = unique(meta[cancer_type == "LUAD", sampleID])
)

lusc_ids <- list(
  cellID   = meta[cancer_type == "LUSC", cellID],
  sampleID = unique(meta[cancer_type == "LUSC", sampleID])
)

# -----------------------------
# 5) Verification against previously validated counts
# -----------------------------
expected_n_patients <- c(LUAD = 61, LUSC = 173)
observed_n_patients <- c(LUAD = length(luad_ids$sampleID), LUSC = length(lusc_ids$sampleID))

message("--- Separation check ---")
message("Expected patients: LUAD = ", expected_n_patients["LUAD"], ", LUSC = ", expected_n_patients["LUSC"])
message("Observed patients:  LUAD = ", observed_n_patients["LUAD"], ", LUSC = ", observed_n_patients["LUSC"])

if (!all(as.numeric(observed_n_patients) == as.numeric(expected_n_patients))) {
  stop("Patient counts after separation do not match previously verified values. ",
       "Do not proceed until this is understood; check for duplicated sampleID, ",
       "unexpected NA in cancer_type, or a metadata file mismatch.")
}
message("Separation matches previously verified patient counts. OK.")

# -----------------------------
# 6) Save reference ID lists
# -----------------------------
dir.create(OBJECTS_DIR, recursive = TRUE, showWarnings = FALSE)

saveRDS(luad_ids, file.path(OBJECTS_DIR, "00_LUAD_reference_ids.rds"))
saveRDS(lusc_ids, file.path(OBJECTS_DIR, "00_LUSC_reference_ids.rds"))
message("Saved: Objects/00_LUAD_reference_ids.rds")
message("Saved: Objects/00_LUSC_reference_ids.rds")

# -----------------------------
# 7) Save summary table
# -----------------------------
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)

summary_table <- data.table(
  cancer_type  = c("LUAD", "LUSC"),
  n_patients   = c(observed_n_patients["LUAD"], observed_n_patients["LUSC"]),
  n_cells      = c(length(luad_ids$cellID), length(lusc_ids$cellID))
)

fwrite(summary_table, file.path(TABLES_DIR, "00_separation_summary.csv"))
message("Saved: Results/Bloc0_DataPrep/Tables/00_separation_summary.csv")

message("DONE Script 00 — LUAD/LUSC separation")