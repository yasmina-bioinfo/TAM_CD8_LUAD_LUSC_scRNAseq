# ============================================================
# Bloc7_TAM — Script09_WangMarkers_Test_LUSC.R
# Purpose: test whether the per-patient proportion of each Wang marker
# group differs between MPR and non-MPR tumors.
#
# EXPLORATORY, POST HOC: this analysis was decided after reading Wang
# et al. (2022) and observing the marker distribution. Reported as such.
#
# Critical design point: patients with zero cells assigned to a given
# marker group must be counted as 0 percent, NOT excluded. Excluding
# them would drop the most extreme values and bias the comparison —
# the same pitfall encountered in Bloc5 with count(.drop = FALSE).
#
# Input : Objects/07_LUSC_seu_TAM_pure_scored.rds
# Output: Results/Bloc7_TAM/Tables/09_LUSC_WangMarkers_tests.csv
# ============================================================

suppressPackageStartupMessages({
  library(here); library(Seurat); library(dplyr)
  library(tidyr); library(readr); library(tibble)
})

HISTOLOGY <- "LUSC"
seu <- readRDS(here("Objects", "07_LUSC_seu_TAM_pure_scored.rds"))
OUT_TAB <- here("Results", "Bloc7_TAM", "Tables")

MK        <- c("FABP4","FCN1","SPP1","SELENOP")
MIN_CELLS <- 10   # per patient, assigned cells only

# ---- 1. Assign each cell to its highest-expressing marker --------------
expr <- as.matrix(GetAssayData(seu, assay = "RNA", layer = "data")[MK, ])
seu$wang_top <- apply(expr, 2, function(x)
  if (max(x) == 0) NA_character_ else MK[which.max(x)])

# ---- 2. Per-patient proportions, with explicit zeros ------------------
md <- seu@meta.data %>%
  as_tibble() %>%
  filter(!is.na(wang_top)) %>%
  mutate(wang_top = factor(wang_top, levels = MK))

res <- md %>%
  count(sampleID, pathological_response, wang_top, .drop = FALSE) %>%
  # .drop = FALSE keeps zero combinations because wang_top is a factor
  group_by(sampleID) %>%
  mutate(n_assigned = sum(n),
         pct        = 100 * n / n_assigned) %>%
  ungroup() %>%
  # response can be lost on zero rows: refill from the patient table
  select(-pathological_response) %>%
  left_join(distinct(md, sampleID, pathological_response), by = "sampleID") %>%
  filter(n_assigned >= MIN_CELLS)

cat("=== patients retained (identical across marker groups) ===\n")
print(res %>% distinct(sampleID, pathological_response) %>%
        count(pathological_response))

cat("\n=== per-patient proportions ===\n")
print(res %>% select(sampleID, pathological_response, wang_top, n, pct) %>%
        arrange(wang_top, desc(pct)), n = Inf)

# ---- 3. Wilcoxon per marker group, BH across the four -----------------
tests <- res %>%
  group_by(wang_top) %>%
  summarise(
    n_MPR         = sum(pathological_response == "MPR"),
    n_nonMPR      = sum(pathological_response == "non-MPR"),
    median_MPR    = median(pct[pathological_response == "MPR"]),
    median_nonMPR = median(pct[pathological_response == "non-MPR"]),
    p_value       = wilcox.test(pct ~ pathological_response)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

cat("\n=== TESTS (Wilcoxon patient-level, BH across 4 marker groups) ===\n")
print(tests)

write_csv(res,   file.path(OUT_TAB, "09_LUSC_WangMarkers_patient_level.csv"))
write_csv(tests, file.path(OUT_TAB, "09_LUSC_WangMarkers_tests.csv"))

cat("\nScript09 (Wang marker tests,", HISTOLOGY, ") complete.\n")