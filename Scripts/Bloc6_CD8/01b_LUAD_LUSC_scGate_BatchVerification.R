# ============================================================
# Bloc6_CD8 — 01b_LUAD_scGate_BatchVerification.R
# Purpose: verify that running scGate in batches of whole patients
# (32 GB RAM limit) did not alter the CD8T population, and did not
# create any difference between response groups (MPR / non-MPR).
#
# Run AFTER 01_LUAD_scGate.R (batched run) and, for the comparison,
# 01a_LUAD_scGate_WholeObject.R (whole-object run, LUAD only).
#
# Structure:
#   Part A — Exploration (LUAD only): how different are the two runs?
#            Measurements only, no threshold.
#   Part B — Rules: PASS / REVIEW checks.
#            Check 1 (totals): both histologies.
#            Checks 2 and 3 (comparison with the whole-object run): LUAD only.
#
# The thresholds of checks 2 and 3 are practical decision rules, set
# after seeing the LUAD measurements of Part A, with a margin. They
# are not published standards: report the measured values.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
})

IN_OBJ    <- here("Objects", "01_LUAD_seu_scGate.rds")               # batched run
WHOLE_OBJ <- here("Objects", "01_LUAD_seu_scGate_wholeobject.rds")   # whole-object run (LUAD only)
OUT_TAB   <- here("Results", "Bloc6_CD8", "Tables")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

# Expected number of patients per response group
EXPECTED_PATIENTS <- c("MPR" = 10, "non-MPR" = 40)    # LUSC: c("MPR" = 32, "non-MPR" = 68)

# ------------------------------------------------------------
# SETUP — load only the cell metadata (much lighter than the object)
# ------------------------------------------------------------
seu <- readRDS(IN_OBJ)
md  <- seu@meta.data
rm(seu); gc()
stopifnot(all(c("sampleID", "pathological_response",
                "scGate_multi", "scGate_batch") %in% colnames(md)))
md$pathological_response <- as.character(md$pathological_response)

# Cells without a scGate label are "Unclassified"
norm <- function(x) ifelse(is.na(x), "Unclassified", as.character(x))
md$new_label <- norm(md$scGate_multi)          # label in the batched run
md$is_new    <- md$new_label == "CD8T"

# Jaccard index = cells in common / cells in either run (1 = identical)
jaccard <- function(a, b) sum(a & b) / sum(a | b)

# Whole-object run available? (LUAD yes, LUSC no: does not fit in 32 GB RAM)
has_whole <- file.exists(WHOLE_OBJ)
if (has_whole) {
  w   <- readRDS(WHOLE_OBJ)
  ref <- setNames(as.character(w$scGate_multi), colnames(w))
  rm(w); gc()
  stopifnot(all(rownames(md) %in% names(ref)))
  md$ref_label <- norm(ref[rownames(md)])       # label in the whole-object run
  md$is_ref    <- md$ref_label == "CD8T"
}

# ============================================================
# PART A — EXPLORATION (LUAD only)
# Measurements without threshold: how different are the two runs?
# These values are the basis of the rules used in Part B.
# ============================================================
if (!has_whole) {
  message("Part A SKIPPED: no whole-object run available ",
          "(whole object does not fit in 32 GB RAM).")
} else {
  
  # A1 — All cells
  # OBJECTIVE : how many cells get the same label in both runs?
  # HOW TO READ: high = the runs agree. Most differences are cells near
  #   the scGate threshold moving to or from "Unclassified".
  message("A1 same label in both runs (all cells): ",
          round(100 * mean(md$new_label == md$ref_label), 1), " %")
  
  # A2 — CD8T population
  # OBJECTIVE : do both runs keep the same CD8T cells?
  # HOW TO READ: retained = share of whole-run CD8T still CD8T in the
  #   batched run; Jaccard = shared CD8T / CD8T in either run (1 = identical).
  retained <- sum(md$is_new & md$is_ref) / sum(md$is_ref)
  jac_all  <- jaccard(md$is_new, md$is_ref)
  message("A2 CD8T whole run: ", sum(md$is_ref), " | batched: ", sum(md$is_new),
          " | shared: ", sum(md$is_new & md$is_ref),
          " | retained: ", round(retained, 3),
          " | Jaccard: ", round(jac_all, 3))
  
  # A3 — By response group
  # OBJECTIVE : is the agreement the same in MPR and non-MPR?
  # HOW TO READ: similar values = the batching has no response-related effect.
  jac_grp <- sapply(names(EXPECTED_PATIENTS), function(g) {
    i <- md$pathological_response == g
    jaccard(md$is_new[i], md$is_ref[i])
  })
  message("A3 Jaccard CD8T by response group: ",
          paste(names(jac_grp), round(jac_grp, 3), collapse = " / "))
  
  # A4 — By batch
  # OBJECTIVE : does one batch behave differently from the others?
  # HOW TO READ: look for a batch clearly below the others, then check
  #   whether it also affects CD8T (jaccard_CD8T column).
  batch_cmp <- md %>%
    mutate(same = new_label == ref_label) %>%
    group_by(scGate_batch) %>%
    summarise(n_cells = n(),
              same_label   = round(mean(same), 3),
              jaccard_CD8T = round(sum(is_new & is_ref) / sum(is_new | is_ref), 3),
              .groups = "drop")
  print(batch_cmp, n = Inf)
  write_csv(batch_cmp, file.path(OUT_TAB, "01b_LUAD_batched_vs_whole_by_batch.csv"))
  
  # A5 — By patient
  # OBJECTIVE : is each patient's CD8T count preserved? (downstream
  #   analyses compare patients with each other)
  # HOW TO READ: Spearman = same ranking of patients by CD8T proportion
  #   (1 = identical ranking); diff_pct = change in CD8T count per patient
  #   (0 = no change).
  pat <- md %>%
    group_by(sampleID, pathological_response) %>%
    summarise(n_cells = n(), cd8_batched = sum(is_new), cd8_whole = sum(is_ref),
              .groups = "drop") %>%
    mutate(prop_batched = cd8_batched / n_cells,
           prop_whole   = cd8_whole / n_cells,
           diff_pct     = 100 * (cd8_batched - cd8_whole) / pmax(cd8_whole, 1))
  write_csv(pat, file.path(OUT_TAB, "01b_LUAD_per_patient_CD8T_batched_vs_whole.csv"))
  
  rho <- cor(pat$prop_batched, pat$prop_whole, method = "spearman")
  message("A5 Spearman (per-patient CD8T proportion): ", round(rho, 4))
  message("A5 per-patient change in CD8T count (%): min ",
          round(min(pat$diff_pct), 1), ", median ", round(median(pat$diff_pct), 2),
          ", max ", round(max(pat$diff_pct), 1))
}

# ============================================================
# PART B — RULES (PASS / REVIEW)
# ============================================================
results <- list()
add <- function(id, what, value, rule, ok) {
  results[[length(results) + 1]] <<- data.frame(
    check = id, what = what, value = value, rule = rule,
    status = ifelse(ok, "PASS", "REVIEW"), stringsAsFactors = FALSE)
}

# ------------------------------------------------------------
# CHECK 1 — Totals are consistent (both histologies)
# OBJECTIVE : make sure no patient or cell was lost or duplicated
#             when the object was split into batches.
# WHAT IT DOES : counts patients per response group, checks that each
#             patient sits in exactly one batch, and that every cell
#             has a batch and a unique name.
# HOW TO READ : all three must be PASS. A REVIEW means the batching or
#             the merge of results is wrong: fix before going on.
# ------------------------------------------------------------
obs_pat <- md %>%
  distinct(sampleID, pathological_response) %>%
  count(pathological_response)
obs_vec <- setNames(obs_pat$n, obs_pat$pathological_response)
add("1a", "Patients per response group",
    paste(names(obs_vec), obs_vec, collapse = " / "),
    paste(names(EXPECTED_PATIENTS), EXPECTED_PATIENTS, collapse = " / "),
    identical(sort(names(obs_vec)), sort(names(EXPECTED_PATIENTS))) &&
      all(obs_vec[names(EXPECTED_PATIENTS)] == EXPECTED_PATIENTS))

pat_batches <- md %>% distinct(sampleID, scGate_batch) %>% count(sampleID)
add("1b", "Each patient is in exactly one batch",
    paste("max batches per patient =", max(pat_batches$n)),
    "= 1", all(pat_batches$n == 1))

add("1c", "Every cell has a batch and a unique name",
    paste("cells without batch =", sum(is.na(md$scGate_batch)),
          "; duplicated cells =", sum(duplicated(rownames(md)))),
    "0 ; 0", !anyNA(md$scGate_batch) && !anyDuplicated(rownames(md)))

# Information only: cells per batch and response group.
# A batch with no MPR cell is possible when MPR patients are few.
batch_tab <- table(batch = md$scGate_batch, response = md$pathological_response)
message("Cells per batch and response group (information):")
print(batch_tab)
write_csv(as.data.frame.matrix(batch_tab) %>% tibble::rownames_to_column("batch"),
          file.path(OUT_TAB, "01b_LUAD_batch_composition.csv"))

# ------------------------------------------------------------
# CHECKS 2 and 3 — Comparison with the whole-object run (LUAD only)
# They reuse the values measured in Part A.
#
# CHECK 2 — The batched run keeps the same CD8T cells as the whole run
#   OBJECTIVE : splitting into batches must not change which cells
#               scGate calls CD8T.
#   HOW TO READ: both values >= 0.95 = PASS (LUAD measured: 0.986, 0.970).
#               Below: the batching alters the CD8T population; discuss
#               before using the batched results.
#
# CHECK 3 — No effect of the batching linked to the response
#   OBJECTIVE : the batching must not change CD8T calls differently in
#               MPR and non-MPR (it would bias MPR vs non-MPR comparisons).
#   HOW TO READ: 3a: difference between groups <= 0.02
#               (LUAD: 0.976 vs 0.969). 3b: Spearman >= 0.99 (LUAD: 0.9988).
#               3c: no patient beyond +/- 10 % (LUAD: -8.7 % to +2.1 %).
#               A REVIEW on 3a means a response-related effect: stop.
# ------------------------------------------------------------
if (!has_whole) {
  message("Checks 2 and 3 SKIPPED: no whole-object run available ",
          "(whole object does not fit in 32 GB RAM).")
} else {
  add("2a", "Whole-run CD8T cells still CD8T in batched run",
      sprintf("%.3f", retained), ">= 0.95", retained >= 0.95)
  add("2b", "Jaccard index of CD8T (all cells)",
      sprintf("%.3f", jac_all), ">= 0.95", jac_all >= 0.95)
  add("3a", "Jaccard of CD8T within each response group",
      paste(names(jac_grp), sprintf("%.3f", jac_grp), collapse = " / "),
      "difference between groups <= 0.02", abs(diff(jac_grp)) <= 0.02)
  add("3b", "Spearman correlation of per-patient CD8T proportion",
      sprintf("%.4f", rho), ">= 0.99", rho >= 0.99)
  add("3c", "Largest per-patient change in CD8T count",
      sprintf("%.1f %% to %.1f %%", min(pat$diff_pct), max(pat$diff_pct)),
      "within +/- 10 %", all(abs(pat$diff_pct) <= 10))
}

# ============================================================
# SUMMARY
# ============================================================
summary_tab <- bind_rows(results)
print(summary_tab, row.names = FALSE)
write_csv(summary_tab, file.path(OUT_TAB, "01b_LUAD_scGate_batch_verification.csv"))

if (any(summary_tab$status == "REVIEW")) {
  message("REVIEW NEEDED for: ",
          paste(summary_tab$check[summary_tab$status == "REVIEW"], collapse = ", "))
} else {
  message("All checks passed.")
}