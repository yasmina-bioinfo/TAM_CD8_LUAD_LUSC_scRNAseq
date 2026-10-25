# ============================================================
# Bloc7_TAM — Script10_HallmarkMetabolic_LUSC.R
# Purpose: test whether metabolic programme activity differs between MPR
# and non-MPR macrophages, using four MSigDB Hallmark gene sets.
#
# Status: EXPLORATORY, and the THIRD characterisation framework tested on
# this compartment (after the four Cheng axes and the four Wang subtype
# markers). Declared as such.
#
# Perimeter: full lineage_TAM, comprising gating-confirmed macrophages
# and transitional cells along the monocyte-to-macrophage continuum,
# consistent with the perimeter adopted for the functional scoring and
# for the communication analysis in Bloc8.
#
# Rationale for selecting these four sets, fixed BEFORE any computation:
# Yang et al. (2021, Front Immunol) show that the two immunosuppressive
# ============================================================
# Bloc7_TAM — Script10_HallmarkMetabolic_LUSC.R
# Purpose: test whether metabolic programme activity differs between MPR
# and non-MPR macrophages, using four MSigDB Hallmark gene sets.
#
# Status: EXPLORATORY, and the THIRD characterisation framework tested on
# this compartment (after the four Cheng axes and the four Wang subtype
# markers). Declared as such.
#
# Perimeter: Macrophage_pure (scGate-confirmed macrophages), consistent
# with Script08/09 and with the classification established in Script04b.
#
# Rationale for selecting these four sets: Yang et al. (2021, Front
# Immunol) report that two immunosuppressive TAM subtypes they identify
# in NSCLC differ by metabolism rather than phenotype (SPP1+ macrophages
# glycolytic/angiogenic, CCL18+ macrophages fatty-acid-oxidative). This
# motivated testing glycolysis, hypoxia, OXPHOS and fatty acid metabolism
# here, since this cohort's macrophage compartment is SPP1-dominant
# (Script08). No directional prediction is made: whether this metabolic
# axis tracks response status in a post anti-PD-1, treatment-exposed
# cohort is not established by Yang et al.'s treatment-naive data, and
# asserting a direction without a validated cohort-specific basis would
# overstate what is known. This is a purely exploratory, two-sided test.
#
# Input : Objects/05_LUSC_seu_TAM_Macrophage_pure_scored.rds  (scored in Script05)
# Output: Results/Bloc7_TAM/Tables/10_LUSC_Hallmark_metabolic_*.csv
#
# Note on scope: this script tests an exploratory characterization
# framework for the TAM compartment (MSigDB Hallmark metabolic gene sets,
# applied via UCell independently of any established TAM metabolic
# classification). No consensual, established method exists for scoring
# TAM metabolic state at single-cell resolution in this literature
# (Murray et al., 2014, Immunity). For this reason, results from this
# script are documented in the README (Bloc7_TAM, "Functional and
# subtype scoring frameworks (exploratory)") but are not reported in
# the preprint, which relies on compositional analysis and CollecTRI
# transcription factor activity for TAM characterization instead.
# ============================================================

suppressPackageStartupMessages({
  library(here); library(Seurat); library(UCell); library(msigdbr)
  library(dplyr); library(readr); library(tibble)
})

set.seed(42)

HISTOLOGY <- "LUSC"
seu <- readRDS(here("Objects", "05_LUSC_seu_TAM_Macrophage_pure_scored.rds"))
OUT_TAB <- here("Results", "Bloc7_TAM", "Tables")

MIN_CELLS <- 10
SETS <- c("HALLMARK_GLYCOLYSIS",
          "HALLMARK_HYPOXIA",
          "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
          "HALLMARK_FATTY_ACID_METABOLISM")

# ---- 1. Retrieve the four Hallmark sets --------------------------------
hm <- msigdbr(species = "Homo sapiens", collection = "H")
cat("msigdbr version:", as.character(packageVersion("msigdbr")), "\n")

sig <- hm %>%
  filter(gs_name %in% SETS) %>%
  split(.$gs_name) %>%
  lapply(function(x) unique(x$gene_symbol))

cat("\n=== gene set sizes and coverage in this object ===\n")
genes <- rownames(seu)
for (s in names(sig)) {
  present <- sum(sig[[s]] %in% genes)
  cat(sprintf("%-38s %3d genes | %3d present (%.0f%%)\n",
              s, length(sig[[s]]), present,
              100 * present / length(sig[[s]])))
}

# Restrict each set to genes present, so scores are not diluted by absence
sig <- lapply(sig, function(x) intersect(x, genes))

# ---- 2. UCell scoring --------------------------------------------------
seu <- AddModuleScore_UCell(seu, features = sig, ncores = 4)
cols <- paste0(names(sig), "_UCell")

cat("\n=== score distributions ===\n")
print(summary(seu@meta.data[, cols]))

cat("\n=== % of cells with a zero score ===\n")
print(round(sapply(cols, function(v) 100 * mean(seu@meta.data[[v]] == 0)), 1))

# ---- 3. Per-patient medians and Wilcoxon, BH across the four sets ------
res <- seu@meta.data %>%
  as_tibble() %>%
  group_by(sampleID, pathological_response) %>%
  summarise(n_tam = n(), across(all_of(cols), median), .groups = "drop") %>%
  filter(n_tam >= MIN_CELLS)

cat("\n=== patients retained ===\n")
print(table(res$pathological_response))

run_test <- function(df, var)
  wilcox.test(as.formula(paste0("`", var, "` ~ pathological_response")),
              data = df)$p.value

med <- function(df, var, grp)
  median(df[[var]][df$pathological_response == grp])

tests <- tibble(
  hallmark      = cols,
  median_MPR    = sapply(cols, function(v) med(res, v, "MPR")),
  median_nonMPR = sapply(cols, function(v) med(res, v, "non-MPR")),
  p_value       = sapply(cols, function(v) run_test(res, v))
) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

cat("\n=== TESTS (Wilcoxon patient-level, BH across 4 Hallmark sets) ===\n")
print(tests)

# ---- 4. Correlation between the four programs ------------------------
cat("\n=== Spearman correlation between programs (patient-level) ===\n")
print(round(cor(res[, cols], method = "spearman"), 3))

write_csv(res,   file.path(OUT_TAB, "10_LUSC_Hallmark_metabolic_fullLineage_patient_level.csv"))
write_csv(tests, file.path(OUT_TAB, "10_LUSC_Hallmark_metabolic_fullLineage_tests.csv"))

cat("\nScript10 (Hallmark metabolic,", HISTOLOGY, ") complete.\n")