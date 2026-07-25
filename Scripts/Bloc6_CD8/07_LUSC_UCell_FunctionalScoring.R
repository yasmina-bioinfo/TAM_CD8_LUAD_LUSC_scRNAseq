# ============================================================
# Bloc6_CD8 — Script07_UCell_FunctionalScoring_PatientLevel_LUSC.R
# Purpose: patient-level re-analysis of CD8.TEX/CD8.TPEX functional
# scores (cytotoxicity, exhaustion, stemness). CD8.TEX receives a full
# patient-level Wilcoxon-BH test (sufficient cells per patient). CD8.TPEX
# is reported descriptively only (insufficient cells per patient for a
# reliable test in either histology; see Methods).
# Patient inclusion consistent with Script05 (≥50 CD8 cells/patient).
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(scales)
})

IN_SCORES <- here("Objects", "06_LUSC_seu_functional_scores.rds")
IN_PATIENTS <- here("Results", "Bloc6_CD8", "Tables",
                    "05_LUSC_CD8_patient_level_proportions.csv")
OUT_TAB <- here("Results", "Bloc6_CD8", "Tables")

message("Loading LUSC functional-score object... ", Sys.time())
seu_ex <- readRDS(IN_SCORES)

# ---- Restrict to patients already retained in Script05 (≥50 CD8 cells) ----
retained_patients <- read_csv(IN_PATIENTS, show_col_types = FALSE) %>%
  distinct(sampleID) %>%
  pull(sampleID)

meta <- seu_ex@meta.data %>%
  filter(sampleID %in% retained_patients)
message("Cells retained (patients passing 50-CD8-cell threshold): ", nrow(meta))

score_cols <- c("Tcell.cytotoxicity", "Tcell.exhaustion", "Tcell.stemness")

# ---- Per-patient median score, per state and signature ----
patient_scores <- meta %>%
  group_by(sampleID, functional.cluster, pathological_response) %>%
  summarise(across(all_of(score_cols), median, .names = "median_{.col}"),
            n_cells = n(), .groups = "drop")

write_csv(patient_scores, file.path(OUT_TAB, "07_LUSC_patient_level_scores_all.csv"))

# ---- CD8.TEX: full patient-level Wilcoxon-BH test ----
tex_scores <- patient_scores %>% filter(functional.cluster == "CD8.TEX")

tex_results <- lapply(score_cols, function(sig) {
  col <- paste0("median_", sig)
  mpr_vals    <- tex_scores[[col]][tex_scores$pathological_response == "MPR"]
  nonmpr_vals <- tex_scores[[col]][tex_scores$pathological_response == "non-MPR"]
  test <- wilcox.test(nonmpr_vals, mpr_vals)
  data.frame(
    signature = sig,
    n_MPR = length(mpr_vals), n_nonMPR = length(nonmpr_vals),
    median_MPR = median(mpr_vals), median_nonMPR = median(nonmpr_vals),
    p_value = test$p.value
  )
})
tex_results <- bind_rows(tex_results)
tex_results$p_adj_BH <- p.adjust(tex_results$p_value, method = "BH")
tex_results$p_display <- scales::pvalue(tex_results$p_adj_BH, accuracy = 0.001, add_p = TRUE)

message("CD8.TEX patient-level results:")
print(tex_results)
write_csv(tex_results, file.path(OUT_TAB, "07_LUSC_TEX_patientlevel_wilcoxon.csv"))

# ---- CD8.TPEX: descriptive only, no formal test ----
tpex_scores <- patient_scores %>% filter(functional.cluster == "CD8.TPEX")

message("CD8.TPEX cell count per patient — distribution check:")
print(summary(tpex_scores$n_cells))
message("Patients with <10 CD8.TPEX cells: ",
        sum(tpex_scores$n_cells < 10), " / ", nrow(tpex_scores))

tpex_descriptive <- tpex_scores %>%
  group_by(pathological_response) %>%
  summarise(across(all_of(paste0("median_", score_cols)),
                   list(median = median, IQR = IQR), .names = "{.col}_{.fn}"))

message("CD8.TPEX descriptive summary (no formal test — see Methods):")
print(tpex_descriptive)
write_csv(tpex_scores, file.path(OUT_TAB, "07_LUSC_TPEX_descriptive_only.csv"))

message("Script07 (patient-level functional scoring, LUSC) complete. ", Sys.time())