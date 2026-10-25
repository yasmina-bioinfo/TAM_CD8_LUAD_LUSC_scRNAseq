# ============================================================
# Bloc7_TAM — Script07_M1Decomposition_LUSC.R
# Purpose: identify which functional register within the Cheng M1
# signature drives the lower M1 score observed in responders.
#
# EXPLORATORY, POST HOC: this decomposition was decided after observing
# the M1 result. Reported descriptively, not as a formal test of the
# same standing as the four pre-specified axes.
#
# Rationale: the Cheng M1 signature mixes registers —
#   - inflammatory cytokines (TNF, IL1A, IL1B, IL6, IL23A)
#   - IFN-induced chemokines (CXCL9/10/11, CCL5), which recruit T cells
#   - co-stimulation / antigen presentation (CD86, CD40, CCR7, IRF1, IRF5)
#   - tryptophan catabolism (IDO1, KYNU), functionally immunosuppressive
#     but IFN-inducible, hence their presence in an "M1" signature
# A lower score could reflect reduced inflammation OR reduced IFN-gamma
# signalling — two different biological readings.
#
# Input : Objects/05_LUSC_seu_TAM_Macrophage_pure_scored.rds  (scored in Script05)
# Output: Results/Bloc7_TAM/Tables/07_LUSC_M1_decomposition.csv
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
  library(here); library(Seurat); library(UCell)
  library(dplyr); library(readr); library(tibble)
})

set.seed(42)

HISTOLOGY <- "LUSC"
seu <- readRDS(here("Objects", "05_LUSC_seu_TAM_Macrophage_pure_scored.rds"))
OUT_TAB <- here("Results", "Bloc7_TAM", "Tables")

sub_sig <- list(
  M1_cytokines = c("TNF","IL1A","IL1B","IL6","IL23A"),
  M1_IFNchemo  = c("CXCL9","CXCL10","CXCL11","CCL5"),
  M1_costim    = c("CD86","CD40","CCR7","IRF1","IRF5"),
  M1_trp       = c("IDO1","KYNU")
)

seu <- AddModuleScore_UCell(seu, features = sub_sig, ncores = 4)
cols <- paste0(names(sub_sig), "_UCell")

res <- seu@meta.data %>%
  as_tibble() %>%
  group_by(sampleID, pathological_response) %>%
  summarise(n_tam = n(),
            across(all_of(cols), median),
            .groups = "drop") %>%
  filter(n_tam >= 10)

out <- tibble(
  register      = cols,
  median_MPR    = sapply(cols, function(v)
    median(res[[v]][res$pathological_response == "MPR"])),
  median_nonMPR = sapply(cols, function(v)
    median(res[[v]][res$pathological_response == "non-MPR"])),
  p_value       = sapply(cols, function(v)
    wilcox.test(as.formula(paste0("`", v, "` ~ pathological_response")),
                data = res)$p.value)
) %>% arrange(p_value)

cat("\n=== M1 decomposition (exploratory, post hoc) ===\n")
print(out)

write_csv(res, file.path(OUT_TAB, "07_LUSC_M1_decomposition_patient_level.csv"))
write_csv(out, file.path(OUT_TAB, "07_LUSC_M1_decomposition.csv"))

cat("\nScript07 (M1 decomposition,", HISTOLOGY, ") complete.\n")