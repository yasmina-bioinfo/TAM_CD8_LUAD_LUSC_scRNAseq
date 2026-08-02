# ==============================================================================
# Bloc7_TAM — Sensitivity analysis, preflight checks
#
# Purpose: decide the perimeter for re-scoring the four Cheng axes.
#   Main analysis  = scGate-confirmed macrophages only ("Pure")
#   Sensitivity    = full lineage_TAM (Pure + Impure)
#
# This script decides nothing on its own. It is descriptive: it verifies that
# UCell scores can be reused as-is, and reports which patients enter each
# perimeter. No statistical test is run here.
# ==============================================================================

library(Seurat)
library(dplyr)
library(here)

# Histology for this run. Change to "LUSC" and save as the LUSC counterpart.
HIST <- "LUAD"

# Minimum cells per patient. Same threshold as the main analysis (muscat
# default, Crowell et al. 2020), applied here to a different denominator.
MIN_CELLS <- 10

# Names of the four Cheng functional axes as stored by UCell.
AXES <- c("M1_UCell", "M2_UCell", "Angiogenesis_UCell", "Phagocytosis_UCell")


# ---- 1. Load both objects -----------------------------------------------------

# Full lineage_TAM, scored before scGate filtering. This is the sensitivity
# perimeter: it contains transitional cells that the main analysis excluded.
seu <- readRDS(here("Objects", paste0("07_", HIST, "_seu_TAM_scored.rds")))

# scGate-confirmed subset, i.e. the population used in the main analysis.
# Loaded only to verify score identity below.
pure <- readRDS(here("Objects", paste0("07_", HIST, "_seu_TAM_pure_scored.rds")))


# ---- 2. Object integrity ------------------------------------------------------

# Expected cell counts: LUAD 19,238 (14,708 Pure / 4,530 Impure)
#                       LUSC 16,920 (12,175 Pure / 4,745 Impure)
# A mismatch means the wrong object was loaded, or upstream files changed.
cat("\n=== Object dimensions ===\n")
print(dim(seu))

# The purity flag must be present, otherwise the two perimeters cannot be told
# apart and the whole comparison is impossible.
stopifnot("scGate_Macrophage" %in% colnames(seu@meta.data))

# Distribution of the purity flag. useNA catches unassigned cells, which would
# silently break the counts below if left undetected.
cat("\n=== scGate purity flag ===\n")
print(table(seu$scGate_Macrophage, useNA = "ifany"))

# All four axes must already be scored on the full lineage, otherwise the
# sensitivity analysis would require re-running UCell.
stopifnot(all(AXES %in% colnames(seu@meta.data)))


# ---- 3. UCell scores are cell-intrinsic --------------------------------------

# UCell ranks genes within each cell, so a cell's score does not depend on which
# other cells were present when scoring ran. If that holds, scores computed on
# the full lineage are valid for the pure subset too, and no re-scoring is
# needed — any change in results comes from the added cells alone.
# This check turns that theoretical property into a verified fact on these data.
common <- intersect(colnames(pure), colnames(seu))
cat("\n=== Score identity check (cells in common:", length(common), ") ===\n")

for (a in AXES) {
  ok <- identical(pure@meta.data[common, a], seu@meta.data[common, a])
  cat(a, ":", ok, "\n")
  # Stop immediately on failure: if scores differ, scoring was re-run after
  # filtering and the reasoning above collapses.
  stopifnot(ok)
}

# The pure object is no longer needed; freeing it avoids holding two large
# Seurat objects in memory.
rm(pure); gc()


# ---- 4. Patient-level coverage ------------------------------------------------

# One row per patient: total cells, confirmed cells, and purity.
# pct_pure is the key column — a patient entering only through the wider
# perimeter is by construction a low-purity patient, whose median will be
# computed mostly on transitional cells.
cov <- seu@meta.data %>%
  as_tibble() %>%
  group_by(sampleID, pathological_response) %>%
  summarise(
    n_total  = n(),
    n_pure   = sum(scGate_Macrophage == "Pure", na.rm = TRUE),
    pct_pure = round(100 * n_pure / n_total, 1),
    .groups  = "drop"
  ) %>%
  mutate(
    # Qualified in the main analysis: enough scGate-confirmed cells.
    in_main        = n_pure  >= MIN_CELLS,
    # Qualifies in the sensitivity analysis: enough cells overall.
    in_sensitivity = n_total >= MIN_CELLS,
    status = case_when(
      in_main & in_sensitivity  ~ "kept",
      !in_main & in_sensitivity ~ "NEW",        # returns via transitional cells
      TRUE                      ~ "still out"
    )
  ) %>%
  arrange(status, pathological_response, pct_pure)


# ---- 5. Report ----------------------------------------------------------------

# Patient counts per group under each perimeter. The main-analysis counts must
# reproduce those already documented (LUAD 9/37, LUSC 32/62); if they do not,
# something changed upstream and nothing below should be trusted.
cat("\n=== Patients — main analysis (n_pure >=", MIN_CELLS, ") ===\n")
print(table(cov$pathological_response[cov$in_main]))

cat("\n=== Patients — sensitivity (n_total >=", MIN_CELLS, ") ===\n")
print(table(cov$pathological_response[cov$in_sensitivity]))

# Returning patients, with their purity. This is the decision-relevant output:
# how many patients are added, to which response group, and how transitional
# their macrophage compartment is.
cat("\n=== Returning patients ===\n")
print(cov %>% filter(status == "NEW") %>%
        select(sampleID, pathological_response, n_total, n_pure, pct_pure))

# Patients excluded under both perimeters, reported for completeness.
cat("\n=== Still excluded ===\n")
print(cov %>% filter(status == "still out") %>%
        select(sampleID, pathological_response, n_total, n_pure))

# Cohort-level purity medians per response group, for context: a returning
# patient should be read against these, not in isolation.
cat("\n=== Median purity per response group ===\n")
print(cov %>% group_by(pathological_response) %>%
        summarise(median_pct_pure = median(pct_pure), n = n(), .groups = "drop"))


# ---- 6. Save ------------------------------------------------------------------

# Full per-patient table saved as .csv (project convention), so the perimeter
# decision can be re-examined later without re-loading the Seurat object.
write.csv(
  cov,
  here("Results", "Bloc7_TAM", "Tables",
       paste0("07_", HIST, "_TAM_sensitivity_coverage.csv")),
  row.names = FALSE
)

cat("\nDone —", HIST, "\n")