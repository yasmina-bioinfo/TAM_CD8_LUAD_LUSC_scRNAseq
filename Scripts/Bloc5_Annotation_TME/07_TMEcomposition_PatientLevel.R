# ============================================================
# Bloc5_Annotation — Script07_TMEcomposition_PatientLevel.R
# Purpose: re-test TME composition at the patient level, replacing the
# cell-level chi-squared and prop.test analyses.
#
# Rationale: cell-level tests treat cells from the same patient as
# independent observations, which inflates significance. This was
# identified while evaluating propeller on these same data. All other
# comparisons in this project are patient-level; this restores
# consistency across the pipeline.
#
# Populations tested: lineage_TAM and lineage_CD8_revised, as proportions
# of all annotated cells per patient. Selected before inspecting results,
# being the two compartments characterised in Blocs 6 and 7.
#
# Input : Objects/05_{HIST}_seu_annotated.rds
# Output: Results/Bloc5_Annotation/Tables/TME_composition_patient_level.csv
# ============================================================

suppressPackageStartupMessages({
  library(here); library(Seurat); library(dplyr)
  library(tibble); library(readr)
})

MIN_CELLS <- 50   # total annotated cells per patient

test_composition <- function(obj_path, histology) {
  
  seu <- readRDS(obj_path)
  
  cd8_col <- if ("lineage_CD8_revised" %in% colnames(seu@meta.data))
    "lineage_CD8_revised" else "lineage_CD8"
  
  res <- seu@meta.data %>%
    as_tibble() %>%
    group_by(sampleID, pathological_response) %>%
    summarise(n_total = n(),
              pct_TAM = 100 * mean(lineage_TAM),
              pct_CD8 = 100 * mean(.data[[cd8_col]]),
              .groups = "drop") %>%
    filter(n_total >= MIN_CELLS) %>%
    mutate(histology = histology)
  
  cat("\n=== ", histology, " — patients retained ===\n", sep = "")
  print(table(res$pathological_response))
  
  cat("\n=== medians (%) ===\n")
  print(res %>%
          group_by(pathological_response) %>%
          summarise(TAM = round(median(pct_TAM), 2),
                    CD8 = round(median(pct_CD8), 2),
                    .groups = "drop"))
  
  p_tam <- wilcox.test(pct_TAM ~ pathological_response, data = res)$p.value
  p_cd8 <- wilcox.test(pct_CD8 ~ pathological_response, data = res)$p.value
  
  tests <- tibble(
    histology = histology,
    population = c("TAM", "CD8"),
    median_MPR = c(median(res$pct_TAM[res$pathological_response == "MPR"]),
                   median(res$pct_CD8[res$pathological_response == "MPR"])),
    median_nonMPR = c(median(res$pct_TAM[res$pathological_response == "non-MPR"]),
                      median(res$pct_CD8[res$pathological_response == "non-MPR"])),
    p_value = c(p_tam, p_cd8)
  ) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH"))
  
  cat("\n=== tests (Wilcoxon, BH across the two populations) ===\n")
  print(tests)
  
  rm(seu); gc()
  list(res = res, tests = tests)
}

luad <- test_composition(here("Objects", "05_LUAD_seu_annotated.rds"), "LUAD")
lusc <- test_composition(here("Objects", "05_LUSC_seu_annotated.rds"), "LUSC")

OUT_TAB <- here("Results", "Bloc5_Annotation_TME", "Tables")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

write_csv(bind_rows(luad$res, lusc$res),
          file.path(OUT_TAB, "TME_composition_patient_level.csv"))
write_csv(bind_rows(luad$tests, lusc$tests),
          file.path(OUT_TAB, "TME_composition_patient_level_tests.csv"))

cat("\nDone.\n")