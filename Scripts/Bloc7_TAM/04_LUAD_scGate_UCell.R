# ============================================================
# ** NOT RUN — SUPERSEDED BY Script04b **
#
# Reason for replacement: two methodology decisions changed after this
# script was written. (1) scGate is now the pipeline's identity-deciding
# method for the TAM compartment (Macrophage_pure / Monocyte_pure /
# scGate_unresolved), not a validation step alongside CellTypist — no
# crosscheck between the two is performed, consistent with Bloc6 (CD8):
# CellTypist is an exploratory first-pass triage, scGate adjudicates on
# a curated reference atlas; comparing the two re-tests a decision
# already made, not a new result (see Bloc6_CD8/04_LUAD_CrossValidation_
# Summary.R for the same reasoning applied to CD8). (2) UCell scoring
# moved to Script05, restricted to Macrophage_pure only (Cheng et al.
# signatures describe macrophage states, not monocyte states).
#
# See Script04b_scGate_Classification_LUAD.R for the current logic.
# Outputs below are kept in Results/Bloc7_TAM/Tables/ for the record
# (methodological traceability, possible reuse) — do not regenerate them
# by running this script; do not treat them as current pipeline output.
# ============================================================

stop("Script04 is superseded by Script04b — see header above. Do not run.")

# ============================================================
# Bloc7_TAM — Script04_scGate_UCell_LUAD.R
# Purpose: cross-validate the CellTypist macrophage annotation with a
# marker-based scGate model, then score four functional axes with UCell.
#
# Input : Objects/05_LUAD_seu_annotated.rds
# Output: Objects/07_LUAD_seu_TAM_scored.rds
#         Results/Bloc7_TAM/Tables/04_LUAD_scGate_crosscheck.csv
#         Results/Bloc7_TAM/Tables/04_LUAD_TAM_scores_percell.csv
#
# Design notes:
#   - scGate is run on the FULL object, not on the lineage_TAM subset:
#     it must remain independent of the CellTypist annotation it is meant
#     to validate (same rationale as Bloc6 for CD8).
#   - Starting population for scoring: lineage_TAM (macrophages only).
#     Classical monocytes are excluded here; the Cheng signatures describe
#     macrophage states, and scGate's Macrophage model itself excludes
#     monocytes at level 4. Monocytes will be included in Bloc8.
#   - NO cycling / ISG filtering, unlike Bloc6: Prolif-TAM and IFN-TAM are
#     recognised TAM states (Ma et al., 2022), not artefacts. The Bloc6
#     filtering existed to match the ProjecTILs atlas, which excludes such
#     cells by construction; no reference atlas is used here.
#   - Cheng et al. (2021) Table S5 uses two non-HGNC symbols, IL23 and
#     FASL, corrected to IL23A and FASLG (verified present in both
#     histologies, Script01).
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(scGate)
  library(UCell)
  library(dplyr)
  library(readr)
  library(tibble)
})

set.seed(42)

HISTOLOGY <- "LUAD"
IN_OBJ    <- here("Objects", "05_LUAD_seu_annotated.rds")
OUT_OBJ   <- here("Objects", "07_LUAD_seu_TAM_scored.rds")
OUT_TAB   <- here("Results", "Bloc7_TAM", "Tables")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Load ------------------------------------------------------------
message("Loading ", HISTOLOGY, " annotated object... ", Sys.time())
seu <- readRDS(IN_OBJ)
message("Total cells: ", ncol(seu))

# ---- 2. scGate purification (independent of CellTypist) -----------------
message("Running scGate with TME_HiRes/Macrophage model... ", Sys.time())

models <- get_scGateDB()

seu_diet <- DietSeurat(seu, dimreducs = "harmony")
seu_diet <- scGate(seu_diet,
                   model  = models$human$TME_HiRes$Macrophage,
                   ncores = 4)

seu$scGate_Macrophage <- seu_diet$is.pure[colnames(seu)]
rm(seu_diet); gc()

message("scGate result:")
print(table(seu$scGate_Macrophage, useNA = "ifany"))

# ---- 3. Cross-check scGate against CellTypist ---------------------------
crosscheck <- table(annotation_final = seu$annotation_final,
                    scGate_pure      = seu$scGate_Macrophage,
                    useNA = "ifany")

cc_df <- as.data.frame.matrix(crosscheck) %>%
  rownames_to_column("annotation_final")
write_csv(cc_df, file.path(OUT_TAB, "04_LUAD_scGate_crosscheck.csv"))

cat("\n=== retention rate within lineage_TAM ===\n")
tam_cells <- colnames(seu)[seu$lineage_TAM]
rate <- mean(seu$scGate_Macrophage[tam_cells] == "Pure", na.rm = TRUE)
cat(sprintf("lineage_TAM cells: %d | scGate-pure: %d (%.1f%%)\n",
            length(tam_cells),
            sum(seu$scGate_Macrophage[tam_cells] == "Pure", na.rm = TRUE),
            100 * rate))

cat("\n=== retention rate per CellTypist label (myeloid only) ===\n")
myeloid_labels <- names(which(table(seu$annotation_final[
  seu$lineage_TAM | seu$lineage_TAM_related]) > 0))
for (lab in myeloid_labels) {
  idx <- seu$annotation_final == lab
  cat(sprintf("%-25s %5d cells | pure: %.1f%%\n", lab, sum(idx),
              100 * mean(seu$scGate_Macrophage[idx] == "Pure", na.rm = TRUE)))
}

# ---- 4. Subset to the TAM compartment for scoring -----------------------
seu_tam <- subset(seu, subset = lineage_TAM)
message("\nTAM cells retained for scoring: ", ncol(seu_tam))

# ---- 5. UCell functional scoring ----------------------------------------
sig <- list(
  M1  = c("IL23A","TNF","CXCL9","CXCL10","CXCL11","CD86","IL1A","IL1B","IL6",
          "CCL5","IRF5","IRF1","CD40","IDO1","KYNU","CCR7"),
  M2  = c("IL4R","CCL4","CCL13","CCL20","CCL17","CCL18","CCL22","CCL24",
          "LYVE1","VEGFA","VEGFB","VEGFC","VEGFD","EGF","CTSA","CTSB","CTSC",
          "CTSD","TGFB1","TGFB2","TGFB3","MMP14","MMP19","MMP9","CLEC7A",
          "WNT7B","FASLG","TNFSF12","TNFSF8","CD276","VTCN1","MSR1","FN1","IRF4"),
  Angiogenesis = c("CCND2","CCNE1","CD44","CXCR4","EDN1","E2F3","EZH2","FGF18",
                   "FGFR1","FYN","HEY1","ITGAV","JAG1","JAG2","MMP9","NOTCH1",
                   "PDGFA","PTK2","SPP1","STC1","TNFAIP6","TYMP","VAV2","VCAN",
                   "VEGFA"),
  Phagocytosis = c("MRC1","CD163","MERTK","C1QB")
)

message("Scoring functional axes with UCell... ", Sys.time())
seu_tam <- AddModuleScore_UCell(seu_tam, features = sig, ncores = 4)

score_cols <- paste0(names(sig), "_UCell")
cat("\n=== score distributions ===\n")
print(summary(seu_tam@meta.data[, score_cols]))

# ---- 6. Export per-cell scores ------------------------------------------
scores <- seu_tam@meta.data %>%
  rownames_to_column("cellID") %>%
  select(cellID, sampleID, pathological_response, annotation_final,
         scGate_Macrophage, all_of(score_cols))

write_csv(scores, file.path(OUT_TAB, "04_LUAD_TAM_scores_percell.csv"))

# ---- 7. Save ------------------------------------------------------------
seu_tam@assays$RNA@layers$scale.data <- NULL   # drop if present
saveRDS(seu_tam, OUT_OBJ)

message("\nScript04 (scGate + UCell, ", HISTOLOGY, ") complete. ", Sys.time())