#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc1_QC, Script 03: Filter MPR/non-MPR (LUAD)
# ============================================================
suppressPackageStartupMessages({
  library(here)
  library(Seurat)
})

IN_OBJ  <- here("Objects", "02_LUAD_seu_qc.rds")
OUT_OBJ <- here("Objects", "03_LUAD_seu_resp.rds")

message("Loading QC-filtered Seurat object...")
seu <- readRDS(IN_OBJ)

message("Cells before response filtering: ", ncol(seu))
message("Distribution before filtering:")
print(table(seu$pathological_response, useNA = "ifany"))

# --------------------------------------------------------------
# Remove pCR (and any group outside MPR/non-MPR): these patients
# are not used anywhere downstream in the main pipeline (Bloc6/7/8).
# Only a descriptive pCR/total ratio per histology is reported,
# computed separately from clinical metadata (Bloc0 Script02).
# --------------------------------------------------------------
seu_resp <- subset(seu, subset = pathological_response %in% c("MPR", "non-MPR"))

message("Cells after response filtering: ", ncol(seu_resp))
message("Distribution after filtering:")
print(table(seu_resp$pathological_response, useNA = "ifany"))

# Strict check: only MPR and non-MPR should remain
remaining_groups <- unique(seu_resp$pathological_response)
stopifnot(setequal(remaining_groups, c("MPR", "non-MPR")))

saveRDS(seu_resp, OUT_OBJ)
message("Saved: 03_LUAD_seu_resp.rds")
message("DONE Bloc1 Script 03_LUAD")