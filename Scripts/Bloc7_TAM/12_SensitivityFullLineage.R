# ============================================================
# Bloc7_TAM — Script12_SensitivityFullLineage.R
# Purpose: assess whether including transitional cells changes the
# functional scores and their comparison between response groups.
#
# Question: the Bloc7 scoring was restricted to scGate-confirmed
# macrophages, on the grounds that the Cheng signatures describe
# macrophage states. Transitional cells within lineage_TAM were excluded.
# This script asks what including them changes, in both histologies.
#
# Nothing is recomputed: UCell scores were already calculated on the full
# lineage_TAM object at Script04. Only the tests are redone.
#
# Input : Objects/07_{HIST}_seu_TAM_scored.rds
# Output: Results/Bloc7_TAM/Tables/11_{HIST}_sensitivity_full_lineage.csv
# ============================================================

suppressPackageStartupMessages({
  library(here); library(Seurat); library(dplyr)
  library(tibble); library(readr)
})

AXES      <- c("M1_UCell", "M2_UCell", "Angiogenesis_UCell", "Phagocytosis_UCell")
MIN_CELLS <- 10

run_sensitivity <- function(histology) {
  
  seu <- readRDS(here("Objects", paste0("07_", histology, "_seu_TAM_scored.rds")))
  md  <- seu@meta.data %>% as_tibble()
  rm(seu); gc()
  
  # ---- 1. Do transitional cells score differently from pure ones? ------
  cat("\n\n########## ", histology, " ##########\n")
  cat("\n=== median score by scGate status (all cells pooled) ===\n")
  print(md %>%
          group_by(scGate_Macrophage) %>%
          summarise(n = n(),
                    across(all_of(AXES), ~ round(median(.x), 4)),
                    .groups = "drop"))
  
  # ---- 2. Per-patient medians on the FULL lineage ----------------------
  res_full <- md %>%
    group_by(sampleID, pathological_response) %>%
    summarise(n_cells = n(),
              across(all_of(AXES), median),
              .groups = "drop") %>%
    filter(n_cells >= MIN_CELLS)
  
  # ---- 3. Same, restricted to pure cells (Bloc7 perimeter) -------------
  res_pure <- md %>%
    filter(scGate_Macrophage == "Pure") %>%
    group_by(sampleID, pathological_response) %>%
    summarise(n_cells = n(),
              across(all_of(AXES), median),
              .groups = "drop") %>%
    filter(n_cells >= MIN_CELLS)
  
  cat("\n=== patients retained ===\n")
  cat("full lineage: ")
  print(table(res_full$pathological_response))
  cat("pure only:    ")
  print(table(res_pure$pathological_response))
  
  # ---- 4. Tests, both perimeters ---------------------------------------
  test_one <- function(df, perimeter) {
    tibble(
      histology = histology,
      perimeter = perimeter,
      axis = AXES,
      median_MPR = sapply(AXES, function(a)
        median(df[[a]][df$pathological_response == "MPR"])),
      median_nonMPR = sapply(AXES, function(a)
        median(df[[a]][df$pathological_response == "non-MPR"])),
      p_value = sapply(AXES, function(a)
        wilcox.test(df[[a]] ~ df$pathological_response)$p.value)
    ) %>%
      mutate(p_adj = p.adjust(p_value, method = "BH"))
  }
  
  tests <- bind_rows(test_one(res_full, "full_lineage"),
                     test_one(res_pure, "pure_only"))
  
  cat("\n=== tests (Wilcoxon patient-level, BH across 4 axes) ===\n")
  print(tests %>%
          mutate(across(c(median_MPR, median_nonMPR), ~ round(.x, 4)),
                 across(c(p_value, p_adj), ~ round(.x, 4))) %>%
          arrange(axis, perimeter), n = Inf)
  
  list(tests = tests, res_full = res_full)
}

luad <- run_sensitivity("LUAD")
lusc <- run_sensitivity("LUSC")

OUT_TAB <- here("Results", "Bloc7_TAM", "Tables")
write_csv(bind_rows(luad$tests, lusc$tests),
          file.path(OUT_TAB, "11_sensitivity_full_lineage_tests.csv"))

cat("\nScript12 complete.\n")