# ============================================================
# Bloc7_TAM — Script03_scGateModelInspection.R
# Purpose: inspect the scGate marker models before applying them, to
# confirm that the Macrophage model targets the intended population and
# does not overlap with the Monocyte model.
#
# Rationale: scGate is used here as a marker-based cross-validation of
# the CellTypist annotation, as was done for CD8 in Bloc6. Applying a
# model without first reading its definition would make it impossible to
# interpret which cells are retained or discarded, and why.
#
# Input : none (queries the scGate model database)
# Output: console only, plus optional export of the two model definitions
#
# Histology-independent: run once, applies to both LUAD and LUSC.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(scGate)
  library(dplyr)
  library(readr)
})

OUT_TAB <- here("Results", "Bloc7_TAM", "Tables")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

models <- get_scGateDB()

# ---- 1. Which model collections and cell types are available? ----------
cat("=== available model collections (human) ===\n")
print(names(models$human))

cat("\n=== cell types in TME_HiRes ===\n")
print(names(models$human$TME_HiRes))

# ---- 2. Definition of the Macrophage model -----------------------------
# Read to confirm: (a) which markers define the population, (b) which
# populations are actively excluded, (c) at which hierarchical level.
cat("\n=== TME_HiRes / Macrophage ===\n")
print(models$human$TME_HiRes$Macrophage)

# ---- 3. Definition of the Monocyte model -------------------------------
# Checked to confirm the two models do not retain the same cells.
cat("\n=== TME_HiRes / Monocyte ===\n")
print(models$human$TME_HiRes$Monocyte)

# ---- 4. Export both definitions for the record -------------------------
write_csv(models$human$TME_HiRes$Macrophage,
          file.path(OUT_TAB, "03_scGate_model_Macrophage.csv"))
write_csv(models$human$TME_HiRes$Monocyte,
          file.path(OUT_TAB, "03_scGate_model_Monocyte.csv"))

cat("\nScript03 (scGate model inspection) complete.\n")