# ============================================================
# Bloc7_TAM — Script06_FunctionalScoring_UCell_LUAD.R
# Purpose: compare UCell functional axis scores between MPR and non-MPR
# tumors at the patient level, within scGate-confirmed macrophages.
#
# Input : Objects/05_LUAD_seu_TAM_Macrophage_pure_scored.rds  (scored in Script05)
# Output: Results/Bloc7_TAM/Tables/06_LUAD_TAM_FunctionalScoring_UCell_patient_level.csv
#         Results/Bloc7_TAM/Tables/06_LUAD_TAM_FunctionalScoring_UCell_tests.csv
#
# Design:
#   - Per-patient summary = MEDIAN of the UCell score (as in Bloc6)
#   - Patients contributing <10 cells excluded (muscat default)
#   - Wilcoxon rank-sum, BH correction across the four Cheng axes
#   - EXPLORATORY: no directional prediction specified in advance
#
# Note: the four axes are not fully independent — M2 and Angiogenesis
# share VEGFA and MMP9 (2 genes out of 34 and 25), a property of the
# published signatures, identical in both histologies. Concordance
# between these two axes must not be presented as two independent
# results. The empirical correlation between axes is measured in
# section 3 and may differ between histologies.
#
# Note on scope: this script tests an exploratory characterization
# framework for the TAM compartment (Cheng et al. 2021 functional
# signatures, applied via UCell independently of the clustering
# pipeline the source study used to derive them). No consensual,
# established method exists for classifying TAM activation states or
# scoring TAM function at single-cell resolution in this literature
# (Murray et al., 2014, Immunity). For this reason, results from this
# script are documented in the README (Bloc7_TAM, "Functional and
# subtype scoring frameworks (exploratory)") but are not reported in
# the preprint, which relies on compositional analysis and CollecTRI
# transcription factor activity for TAM characterization instead.
# ============================================================

suppressPackageStartupMessages({
  library(here); library(Seurat); library(dplyr)
  library(readr); library(tibble); library(tidyr)
})

HISTOLOGY <- "LUAD"
IN_OBJ    <- here("Objects", "05_LUAD_seu_TAM_Macrophage_pure_scored.rds")
OUT_TAB   <- here("Results", "Bloc7_TAM", "Tables")

MIN_CELLS <- 10
AXES <- c("M1_UCell", "M2_UCell", "Angiogenesis_UCell", "Phagocytosis_UCell")

seu <- readRDS(IN_OBJ)

# ---- 1. Per-patient median score ---------------------------------------
res <- seu@meta.data %>%
  as_tibble() %>%
  group_by(sampleID, pathological_response) %>%
  summarise(n_tam = n(),
            across(all_of(AXES), median, .names = "{.col}"),
            .groups = "drop") %>%
  filter(n_tam >= MIN_CELLS) %>%
  arrange(pathological_response, sampleID)

cat("=== patients retained ===\n")
print(table(res$pathological_response))

print(res, n = Inf)

# ---- 2. Wilcoxon tests, BH across the four axes ------------------------
run_test <- function(df, var) {
  f <- as.formula(paste0("`", var, "` ~ pathological_response"))
  wilcox.test(f, data = df)$p.value
}

med_by_group <- function(df, var) {
  c(MPR    = median(df[[var]][df$pathological_response == "MPR"]),
    nonMPR = median(df[[var]][df$pathological_response == "non-MPR"]))
}

tests <- tibble(
  axis          = AXES,
  median_MPR    = sapply(AXES, function(v) med_by_group(res, v)["MPR"]),
  median_nonMPR = sapply(AXES, function(v) med_by_group(res, v)["nonMPR"]),
  p_value       = sapply(AXES, function(v) run_test(res, v))
) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

cat("\n=== PATIENT-LEVEL TESTS (BH across 4 axes) ===\n")
print(tests)

# ---- 3. Correlation between axes ---------------------------------------
# M2 and Angiogenesis share two genes; check how far the axes track
# each other before interpreting concordance as independent evidence.
cat("\n=== Spearman correlation between axes (patient-level) ===\n")
print(round(cor(res[, AXES], method = "spearman"), 3))

# ---- 4. Save ------------------------------------------------------------
write_csv(res,   file.path(OUT_TAB, "06_LUAD_TAM_FunctionalScoring_UCell_patient_level.csv"))
write_csv(tests, file.path(OUT_TAB, "06_LUAD_TAM_FunctionalScoring_UCell_tests.csv"))

cat("\nScript06 (patient-level comparison,", HISTOLOGY, ") complete.\n")