# ============================================================
# Bloc7_TAM — Script14_CollecTRI_Overlap_LUAD.R
# Purpose: identify which transcription factors are shared across the
# three myeloid groups and which are specific to one, from the output of
# Script13.
#
# Rationale: Script13 tests each group separately along the
# monocyte-to-macrophage continuum. The question this script answers is
# whether the significant factors form a common set, which would point to
# a shared response to the microenvironment, or differ by differentiation
# stage. Factors specific to a single group are listed separately, since
# they are what distinguishes one stage from the others.
#
# No recomputation: this reads the exported table and cross-tabulates it.
#
# Input : Results/Bloc7_TAM/Tables/13_LUAD_CollecTRI_myeloid.csv
# Output: console only (tables are short enough to read directly)
# ============================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(readr); library(tibble)
})

HISTOLOGY <- "LUAD"
GROUPS    <- c("Macrophage_pure", "Transitional", "Monocyte")

acts <- read_csv(here("Results", "Bloc7_TAM", "Tables",
                      paste0("13_", HISTOLOGY, "_CollecTRI_myeloid.csv")),
                 show_col_types = FALSE)

sig <- acts %>% filter(p_adj < 0.05)

# ---- 1. In how many groups is each factor significant? -----------------
# The direction column records whether the factor is more active in MPR
# or non-MPR; a factor significant in several groups could in principle
# differ in direction between them, which this makes visible.
shared <- sig %>%
  group_by(source) %>%
  summarise(n_groups  = n_distinct(group),
            groups    = paste(sort(unique(group)), collapse = " + "),
            direction = paste(unique(ifelse(score > 0, "MPR", "nonMPR")),
                              collapse = "/"),
            .groups = "drop")

cat("=== factors by number of groups in which they are significant ===\n")
print(shared %>% dplyr::count(n_groups, groups, sort = TRUE), n = Inf)

# ---- 2. The set shared across all three groups -------------------------
core <- shared %>% filter(n_groups == 3) %>% pull(source)
cat("\n=== significant in ALL THREE groups (", length(core), ") ===\n", sep = "")
print(sort(core))

cat("\ndirections within the shared set:\n")
print(shared %>% filter(n_groups == 3) %>% dplyr::count(direction))

# ---- 3. What is specific to each group ---------------------------------
for (g in GROUPS) {
  only <- shared %>% filter(n_groups == 1, groups == g) %>% pull(source)
  cat("\n=== specific to ", g, " (", length(only), ") ===\n", sep = "")
  print(sort(only))
}

# ---- 4. Factors more active in responding tumors -----------------------
# Reported separately because they are few: the asymmetry between response
# groups is itself part of the result.
cat("\n=== more active in MPR ===\n")
sig %>%
  filter(score > 0) %>%
  select(group, source, score, p_adj) %>%
  arrange(group, p_adj) %>%
  print(n = Inf)

cat("\n=== counts by direction and group ===\n")
print(sig %>% dplyr::count(group, direction = ifelse(score > 0, "MPR", "nonMPR")))

cat("\nScript14 (overlap, ", HISTOLOGY, ") complete.\n", sep = "")