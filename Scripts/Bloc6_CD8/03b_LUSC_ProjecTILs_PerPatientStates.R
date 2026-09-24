# ============================================================
# Bloc6_CD8 — 03b_LUSC_ProjecTILs_PerPatientStates.R
# Purpose: compute, per patient, the number of cells in each
# ProjecTILs functional state, and derive the per-patient testability
# of CD8.TEX and CD8.TPEX (share of patients with >=1 cell, and among
# those, share below the 10-cell threshold used for UCell functional
# scoring).
# Run AFTER 03_LUSC_ProjecTILs.R. Does not modify the object.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
})

IN_OBJ  <- here("Objects", "03_LUSC_seu_projectils.rds")
OUT_TAB <- here("Results", "Bloc6_CD8", "Tables")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

# Cells and patients expected from Script 03 (guard against a wrong input file)
EXPECTED_CELLS    <- 93584                          
EXPECTED_PATIENTS <- c("MPR" = 32, "non-MPR" = 68)  
MIN_CELLS_STATE   <- 10   # threshold used for UCell functional scoring within a state

message("Loading LUSC ProjecTILs object... ", Sys.time())
seu <- readRDS(IN_OBJ)
md  <- seu@meta.data
rm(seu); gc()
message("Cells loaded: ", nrow(md))
stopifnot(nrow(md) == EXPECTED_CELLS)

obs_pat <- md %>%
  distinct(sampleID, pathological_response) %>%
  count(pathological_response)
obs_vec <- setNames(obs_pat$n, obs_pat$pathological_response)
stopifnot(all(obs_vec[names(EXPECTED_PATIENTS)] == EXPECTED_PATIENTS))
n_patients_total <- sum(EXPECTED_PATIENTS)

# ---- Cells per patient and functional state ----
per_patient_state <- md %>%
  group_by(sampleID, pathological_response, functional.cluster) %>%
  summarise(n_cells = n(), .groups = "drop")
write_csv(per_patient_state,
          file.path(OUT_TAB, "03b_LUSC_cells_per_patient_per_state.csv"))

# ---- Testability summary for every state (all patients pooled) ----
# For each state: % of patients with >=1 cell, and among those,
# % below MIN_CELLS_STATE cells (not enough for functional scoring).
all_patients <- md %>% distinct(sampleID)
states <- sort(unique(md$functional.cluster))

testability <- lapply(states, function(s) {
  n_by_patient <- per_patient_state %>%
    filter(functional.cluster == s) %>%
    right_join(all_patients, by = "sampleID") %>%
    mutate(n_cells = tidyr::replace_na(n_cells, 0))
  has_cell <- n_by_patient$n_cells >= 1
  below_min_among_present <- n_by_patient$n_cells[has_cell] < MIN_CELLS_STATE
  data.frame(
    state = s,
    pct_patients_with_ge1_cell = round(100 * mean(has_cell), 1),
    pct_of_those_below_min     = round(100 * mean(below_min_among_present), 1)
  )
}) %>% bind_rows()
print(testability)
write_csv(testability, file.path(OUT_TAB, "03b_LUSC_state_testability.csv"))
# ---- CD8.TEX / CD8.TPEX zero-cell patients: linked to response? ----
patient_response <- distinct(md, sampleID, pathological_response)

for (s in c("CD8.TEX", "CD8.TPEX")) {
  n_by_patient <- per_patient_state %>%
    filter(functional.cluster == s) %>%
    select(sampleID, n_cells) %>%
    right_join(patient_response, by = "sampleID") %>%
    mutate(n_cells = tidyr::replace_na(n_cells, 0), has_cell = n_cells >= 1)
  zero_by_group <- n_by_patient %>%
    group_by(pathological_response) %>%
    summarise(n_patients = n(), n_zero = sum(!has_cell),
              pct_zero = round(100 * mean(!has_cell), 1), .groups = "drop")
  message(s, " — patients with 0 cells, by response group:")
  print(zero_by_group)
  write_csv(zero_by_group,
            file.path(OUT_TAB, paste0("03b_LUAD_", gsub("\\.", "", s), "_zero_by_response.csv")))
}

message("Script 03b (per-patient state testability, LUSC) complete. ", Sys.time())