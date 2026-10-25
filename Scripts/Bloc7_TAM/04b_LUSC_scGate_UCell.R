# ============================================================
# Bloc7_TAM — Script04b_scGate_Classification_LUSC.R
# Purpose: classify the TAM/TAM_related compartment (CellTypist,
# Bloc5) into scGate-derived Macrophage_pure / Monocyte_pure /
# scGate_unresolved states, then cross it against the CellTypist
# Confirmed/Unresolved/Other status established in Script01.
#
# Input : Objects/05_LUSC_seu_annotated.rds
# Output: Results/Bloc7_TAM/Tables/04b_LUSC_scGate_TAM_class.csv
#         Results/Bloc7_TAM/Tables/04b_LUSC_crosstab_Bloc5_vs_scGate.csv
#
# Design notes (decided after the whole-object run proved
# uninterpretable, see cahier de notes) :
#   - scGate is run ONLY on lineage_TAM | lineage_TAM_related, not on
#     the full object. Rationale: for TAM there is no independent
#     downstream tool (no ProjecTILs equivalent) that requires scGate
#     to remain blind to the CellTypist myeloid gate itself, unlike
#     CD8 (Bloc6) where scGate had to run on the full object because
#     ProjecTILs was applied afterward. Restricting the input here
#     avoids diluting "unresolved" with cells that are simply not
#     myeloid at all (T, B, NK, etc.), which is what happened in the
#     first, whole-object attempt.
#   - No batching: the myeloid subset is small enough (tens of
#     thousands of cells, not hundreds of thousands) to run in one
#     pass on this machine.
#   - Both TME_HiRes models are run (Macrophage and Monocyte), since
#     Script03 confirmed they are mutually exclusive by construction
#     (each excludes the other's core signature at its final level).
#   - This script does NOT attempt to characterise scGate_unresolved
#     with additional marker panels (e.g. Zilionis alveolar/recruited
#     markers). That path was tried before and abandoned: no
#     consensus reference atlas exists for TAM states the way
#     ProjecTILs exists for CD8, so further ad hoc marker checks would
#     just add another unarbitrated opinion on top of CellTypist and
#     scGate. This is recorded as a methodological limitation, not a
#     gap to fill here.
#   - Variable naming: Confirmed_Myeloid_Bloc5 / Unresolved_Bloc5 is
#     used verbatim (not "CellTypist_status" or similar) so it is
#     unambiguous which pipeline stage the variable comes from.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(scGate)
  library(dplyr)
  library(readr)
  library(tibble)
})

set.seed(42)

HISTOLOGY <- "LUSC"
IN_OBJ    <- here("Objects", "05_LUSC_seu_annotated.rds")
OUT_TAB   <- here("Results", "Bloc7_TAM", "Tables")
OUT_OBJ   <- here("Objects", "04b_LUSC_seu_TAM_classified.rds")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

# Expected size of the lineage_TAM | lineage_TAM_related subset,
# verified against Script01 output (README, Bloc5 methodology table).
EXPECTED_TAM_SUBSET_CELLS <- 46390

# ---- 1. Load and subset ---------------------------------------------------
message("Loading ", HISTOLOGY, " annotated object... ", Sys.time())
seu <- readRDS(IN_OBJ)
message("Total cells in object: ", ncol(seu))

seu_tam <- subset(seu, subset = lineage_TAM | lineage_TAM_related)
message("TAM | TAM_related subset: ", ncol(seu_tam), " cells")
stopifnot(ncol(seu_tam) == EXPECTED_TAM_SUBSET_CELLS)

rm(seu); gc()

# ---- 2. Confirmed_Myeloid_Bloc5 / Unresolved_Bloc5 -----------------------
# Intermediate macrophages is never "confirmed" here: Bloc5 had no
# validated marker set to confirm this state directly, so cells
# labelled Intermediate macrophages are treated as unresolved and
# go through scGate below like any other unresolved myeloid cell.
confirmed_labels <- c("Macrophages", "Classical monocytes", "Non-classical monocytes")

seu_tam$Confirmed_Myeloid_Bloc5 <- case_when(
  grepl("Intermediate macrophages", seu_tam$annotation_final)          ~ "Unresolved_Bloc5",
  seu_tam$annotation_final %in% confirmed_labels                       ~ "Confirmed_Myeloid_Bloc5",
  grepl("^Unresolved", seu_tam$annotation_final)                       ~ "Unresolved_Bloc5",
  seu_tam$annotation_final == "Other (mixed, minor)"                   ~ "Other",
  TRUE ~ NA_character_
)
stopifnot(!anyNA(seu_tam$Confirmed_Myeloid_Bloc5))

cat("\n=== Confirmed_Myeloid_Bloc5 / Unresolved_Bloc5 status ===\n")
print(table(seu_tam$Confirmed_Myeloid_Bloc5))

# ---- 3. scGate, restricted to this subset, no batching -------------------
# DietSeurat with dimreducs = "harmony" is mandatory before scGate
# (established in Bloc6: passing NULL breaks downstream steps).
message("\nRunning scGate (Macrophage + Monocyte models) on the TAM subset... ", Sys.time())

models <- get_scGateDB()
seu_diet <- DietSeurat(seu_tam, dimreducs = "harmony")

seu_diet <- scGate(seu_diet,
                   model  = models$human$TME_HiRes$Macrophage,
                   ncores = 4)
seu_tam$scGate_Macrophage <- seu_diet$is.pure[colnames(seu_tam)]

seu_diet <- scGate(seu_diet,
                   model  = models$human$TME_HiRes$Monocyte,
                   ncores = 4)
seu_tam$scGate_Monocyte <- seu_diet$is.pure[colnames(seu_tam)]

rm(seu_diet); gc()

cat("\n=== scGate raw calls ===\n")
cat("Macrophage model:\n"); print(table(seu_tam$scGate_Macrophage, useNA = "ifany"))
cat("Monocyte model:\n");   print(table(seu_tam$scGate_Monocyte,   useNA = "ifany"))

# ---- 4. 3-group scGate_TAM_class ------------------------------------------
# Macrophage_pure / Monocyte_pure / scGate_unresolved. A fourth,
# Ambiguous_double_pure, is kept separate rather than folded into
# scGate_unresolved, since it is a different failure mode (both models
# fire) worth inspecting on its own if it occurs.
seu_tam$scGate_TAM_class <- case_when(
  seu_tam$scGate_Macrophage == "Pure" & seu_tam$scGate_Monocyte == "Pure" ~ "Ambiguous_double_pure",
  seu_tam$scGate_Macrophage == "Pure"                                     ~ "Macrophage_pure",
  seu_tam$scGate_Monocyte   == "Pure"                                     ~ "Monocyte_pure",
  TRUE                                                                    ~ "scGate_unresolved"
)

cat("\n=== scGate_TAM_class ===\n")
print(table(seu_tam$scGate_TAM_class))

# ---- 5. Formal crosstab: Confirmed_Myeloid_Bloc5 x scGate_TAM_class -------
crosstab <- table(Confirmed_Myeloid_Bloc5 = seu_tam$Confirmed_Myeloid_Bloc5,
                  scGate_TAM_class       = seu_tam$scGate_TAM_class)

cat("\n=== crosstab: Confirmed_Myeloid_Bloc5 x scGate_TAM_class ===\n")
print(crosstab)

crosstab_df <- as.data.frame.matrix(crosstab) %>%
  rownames_to_column("Confirmed_Myeloid_Bloc5")
write_csv(crosstab_df, file.path(OUT_TAB, "04b_LUSC_crosstab_Bloc5_vs_scGate.csv"))

# ---- 6. Per-cell export ----------------------------------------------------
percell <- seu_tam@meta.data %>%
  rownames_to_column("cellID") %>%
  select(cellID, sampleID, pathological_response, annotation_final,
         Confirmed_Myeloid_Bloc5, scGate_Macrophage, scGate_Monocyte, scGate_TAM_class)

write_csv(percell, file.path(OUT_TAB, "04b_LUSC_scGate_TAM_class.csv"))

# ---- 7. Save classified object ------------------------------------------
saveRDS(seu_tam, OUT_OBJ)

message("\nScript04b (scGate classification, ", HISTOLOGY, ") complete. ", Sys.time())