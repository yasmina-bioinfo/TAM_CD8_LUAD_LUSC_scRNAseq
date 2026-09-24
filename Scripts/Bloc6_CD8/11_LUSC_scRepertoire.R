# ============================================================
# Bloc6_CD8 — 11_LUSC_TCR_ClonalAnalysis.R
# TCR clonal analysis of the CD8.TEX compartment
#
# Input   : Objects/06_LUSC_seu_functional_scores.rds
#           Data/GSE243013_T_with_TCR_annotation.csv.gz
# Output  : patient-level tables + tests (Results/)
#
# Expansion threshold : >=3 cells per clonotype (Liu et al. 2025, Cell,
#                       STAR Methods — original authors of GSE243013)
# Patient threshold   : >=10 CD8.TEX cells WITH a recovered TCR
# Denominator         : TEX cells with recovered TCR (never all TEX cells)
#
# Pre-registered prediction (27/07/2026, before any result):
#   MPR -> LOWER % of cells in expanded clones, HIGHER Shannon diversity
#   Falsification: no difference on either metric
#
# PRIMARY analysis   : unrarefied, within-TEX, % expanded + Shannon
# SENSITIVITY checks : Gini-Simpson, rarefaction, global clone definition
#
# NOTE: this is the first run of this exact denominator (TCR-bearing
# TEX cells) on the corrected cohort. Patient counts and per-patient
# TCR-library coverage have not been verified elsewhere in this
# pipeline — read them from this run's own output, not from memory.
#
# All count() calls are namespaced dplyr::count() — matrixStats also
# exports a count() that masks dplyr's when loaded first in the session.
# ============================================================

library(here); library(Seurat); library(dplyr); library(readr); library(tibble)

set.seed(42)

HISTOLOGY <- "LUSC"
IN_OBJ    <- here("Objects", paste0("06_", HISTOLOGY, "_seu_functional_scores.rds"))
IN_TCR    <- here("Data", "GSE243013_T_with_TCR_annotation.csv.gz")
OUT_TAB   <- here("Results", "Bloc6_CD8", "Tables")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

MIN_CLONE <- 3    # Liu et al. 2025
MIN_CELLS <- 10   # per patient, TCR-bearing TEX cells only

# CD8.TEX cells expected from Script 06 (guard against a wrong input file)
EXPECTED_TEX_CELLS <- 21571   

# ---- Helper functions --------------------------------------------------

shannon_H <- function(clono) {
  p <- table(clono)
  p <- p / sum(p)
  as.numeric(-sum(p * log(p)))
}

gini_simpson <- function(clono) {
  p <- table(clono)
  p <- p / sum(p)
  as.numeric(1 - sum(p^2))
}

pct_expanded <- function(clono, min_size = MIN_CLONE) {
  tab <- table(clono)
  as.numeric(sum(tab[tab >= min_size]) / sum(tab) * 100)
}

run_test <- function(df, var) {
  f <- as.formula(paste(var, "~ pathological_response"))
  wilcox.test(f, data = df)$p.value
}

med_by_group <- function(df, var) {
  c(MPR    = median(df[[var]][df$pathological_response == "MPR"]),
    nonMPR = median(df[[var]][df$pathological_response == "non-MPR"]))
}

# ---- 1. Load -----------------------------------------------------------

seu <- readRDS(IN_OBJ)
message("CD8.TEX cells loaded: ", ncol(seu))
stopifnot(ncol(seu) == EXPECTED_TEX_CELLS)

tcr <- read_csv(IN_TCR, show_col_types = FALSE)
stopifnot(!any(duplicated(tcr$cellID)))

# ---- 2. TEX metadata + clonotype join ----------------------------------

md <- seu@meta.data %>%
  rownames_to_column("cellID") %>%
  filter(functional.cluster == "CD8.TEX") %>%
  select(cellID, sampleID, pathological_response) %>%
  left_join(tcr %>% select(cellID, clonotype), by = "cellID")

dat <- md %>% filter(!is.na(clonotype))
message("TEX cells with a recovered TCR: ", nrow(dat), " of ", nrow(md))

# ---- 3. Patient threshold ----------------------------------------------

keep <- dat %>% dplyr::count(sampleID) %>% filter(n >= MIN_CELLS) %>% pull(sampleID)
dat  <- dat %>% filter(sampleID %in% keep)

cat("Patients retained (>=", MIN_CELLS, "TCR-bearing CD8.TEX cells):\n")
print(table(distinct(dat, sampleID, pathological_response)$pathological_response))

# ---- 4. Patient-level metrics ------------------------------------------

size_global <- tcr %>% dplyr::count(clonotype, name = "size_global")

res <- dat %>%
  left_join(size_global, by = "clonotype") %>%
  group_by(sampleID, pathological_response) %>%
  summarise(
    n_cells             = n(),
    n_clonotypes        = n_distinct(clonotype),
    pct_expanded        = pct_expanded(clonotype),
    shannon             = shannon_H(clonotype),
    gini_simpson        = gini_simpson(clonotype),
    pct_expanded_global = 100 * mean(size_global >= MIN_CLONE),
    .groups = "drop"
  ) %>%
  arrange(pathological_response, sampleID)

print(res, n = Inf)

# ---- 5. PRIMARY tests (Wilcoxon, BH across the 2 pre-specified metrics)

p_raw <- c(pct_expanded = run_test(res, "pct_expanded"),
           shannon      = run_test(res, "shannon"))

tests <- tibble(
  metric        = names(p_raw),
  median_MPR    = c(med_by_group(res, "pct_expanded")["MPR"],
                    med_by_group(res, "shannon")["MPR"]),
  median_nonMPR = c(med_by_group(res, "pct_expanded")["nonMPR"],
                    med_by_group(res, "shannon")["nonMPR"]),
  p_value       = as.numeric(p_raw),
  p_adj         = p.adjust(as.numeric(p_raw), method = "BH")
)

cat("\n=== PRIMARY tests ===\n"); print(tests)

# ---- 6. Sensitivity A: Gini-Simpson ------------------------------------

gs_med <- med_by_group(res, "gini_simpson")
cat("\n=== Gini-Simpson (sensitivity) ===\n")
cat("median MPR =", round(gs_med["MPR"], 4),
    "| median non-MPR =", round(gs_med["nonMPR"], 4),
    "| p =", signif(run_test(res, "gini_simpson"), 3), "\n")

# ---- 7. Sensitivity B: residual link between Shannon and depth ---------

sp <- suppressWarnings(cor.test(res$shannon, res$n_cells, method = "spearman"))
cat("\nSpearman Shannon ~ n_cells: rho =", round(sp$estimate, 3),
    "| p =", signif(sp$p.value, 3), "\n")

# ---- 8. Sensitivity C: rarefaction to common depth ---------------------

DEPTH    <- 30
N_DRAWS  <- 100
keep_rar <- res %>% filter(n_cells >= DEPTH) %>% pull(sampleID)

rar <- dat %>%
  filter(sampleID %in% keep_rar) %>%
  group_by(sampleID, pathological_response) %>%
  summarise(
    shannon_rar      = mean(replicate(N_DRAWS, shannon_H(sample(clonotype, DEPTH)))),
    pct_expanded_rar = mean(replicate(N_DRAWS, pct_expanded(sample(clonotype, DEPTH)))),
    gini_rar         = mean(replicate(N_DRAWS, gini_simpson(sample(clonotype, DEPTH)))),
    .groups = "drop"
  )

cat("\nRarefied to", DEPTH, "cells:\n"); print(table(rar$pathological_response))

rar_medians <- rar %>%
  group_by(pathological_response) %>%
  summarise(median_pct_expanded = median(pct_expanded_rar),
            median_shannon      = median(shannon_rar),
            median_gini         = median(gini_rar),
            n_patients          = n(), .groups = "drop")

cat("\n=== Rarefied medians ===\n"); print(rar_medians)

p_rar <- c(pct_expanded_rar = run_test(rar, "pct_expanded_rar"),
           shannon_rar      = run_test(rar, "shannon_rar"))

cat("\n=== Rarefaction tests ===\n")
print(tibble(metric  = names(p_rar),
             p_value = as.numeric(p_rar),
             p_adj   = p.adjust(as.numeric(p_rar), method = "BH")))
cat("Gini-Simpson rarefied: p =", signif(run_test(rar, "gini_rar"), 3), "\n")

# ---- 9. Sensitivity D: expansion defined over all T cells --------------

cat("\n=== Expansion defined over all T cells ===\n")
print(wilcox.test(pct_expanded_global ~ pathological_response, data = res))

# ---- 10. Save ----------------------------------------------------------

write_csv(res,         file.path(OUT_TAB, paste0("11_", HISTOLOGY, "_TEX_clonal_patient_level.csv")))
write_csv(tests,       file.path(OUT_TAB, paste0("11_", HISTOLOGY, "_TEX_clonal_tests.csv")))
write_csv(rar,         file.path(OUT_TAB, paste0("11_", HISTOLOGY, "_TEX_clonal_rarefied.csv")))
write_csv(rar_medians, file.path(OUT_TAB, paste0("11_", HISTOLOGY, "_TEX_clonal_rarefied_medians.csv")))

cat("\nDone —", HISTOLOGY, "\n")