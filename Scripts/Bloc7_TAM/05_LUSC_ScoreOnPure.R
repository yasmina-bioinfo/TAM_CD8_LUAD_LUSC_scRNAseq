# ============================================================
# Bloc7_TAM — Script05_ScoreOnPure_LUSC.R
# Purpose: restrict the TAM compartment to scGate-confirmed macrophages
# and compute functional scores on that population.
#
# Rationale: Macrophage_pure is the population scGate confidently
# assigns as macrophage identity (see Script04b). Monocyte_pure,
# scGate_unresolved and Ambiguous_double_pure cells are excluded here:
# the Cheng et al. signatures describe macrophage states, and no
# consensus reference exists to characterise the excluded groups
# further (see Bloc7 Methods).
#
# Input : Objects/04b_LUSC_seu_TAM_classified.rds  (Script04b output,
#         lineage_TAM | lineage_TAM_related, scGate-classified)
# Output: Objects/05_LUSC_seu_TAM_Macrophage_pure_scored.rds
#         Results/Bloc7_TAM/Tables/05_LUSC_TAM_pure_scores_percell.csv
# ============================================================

suppressPackageStartupMessages({
  library(here); library(Seurat); library(UCell)
  library(dplyr); library(readr); library(tibble)
})

set.seed(42)

HISTOLOGY <- "LUSC"
IN_OBJ    <- here("Objects", "04b_LUSC_seu_TAM_classified.rds")
OUT_OBJ   <- here("Objects", "05_LUSC_seu_TAM_Macrophage_pure_scored.rds")
OUT_TAB   <- here("Results", "Bloc7_TAM", "Tables")

seu <- readRDS(IN_OBJ)
message("Cells before filtering: ", ncol(seu))

# ---- Restrict to scGate-confirmed macrophages --------------------------
seu <- subset(seu, subset = scGate_TAM_class == "Macrophage_pure")
message("Cells after scGate filtering: ", ncol(seu))

cat("\n=== cells per response group ===\n")
print(table(seu$pathological_response))

# ---- Compute scores on the restricted population ------------------------
sig <- list(
  M1  = c("IL23A","TNF","CXCL9","CXCL10","CXCL11","CD86","IL1A","IL1B","IL6",
          "CCL5","IRF5","IRF1","CD40","IDO1","KYNU","CCR7"),
  M2  = c("IL4R","CCL4","CCL13","CCL20","CCL17","CCL18","CCL22","CCL24",
          "LYVE1","VEGFA","VEGFB","VEGFC","VEGFD","EGF","CTSA","CTSB","CTSC",
          "CTSD","TGFB1","TGFB2","TGFB3","MMP14","MMP19","MMP9","CLEC7A",
          "WNT7B","FASLG","TNFSF12","TNFSF8","CD276","VTCN1","MSR1","FN1","IRF4"),
  Angiogenesis = c("CCND2","CCNE1","CD44","CXCR4","EDN1","E2F3","EZH2","FGF18",
                   "FGFR1","FYN","HEY1","ITGAV","JAG1","JAG2","MMP9","NOTCH1",
                   "PDGFA","PTK2","SPP1","STC1","TNFAIP6","TYMP","VAV2","VCAN",
                   "VEGFA"),
  Phagocytosis = c("MRC1","CD163","MERTK","C1QB")
)

seu <- AddModuleScore_UCell(seu, features = sig, ncores = 4)
score_cols <- paste0(names(sig), "_UCell")

cat("\n=== score distributions (Macrophage_pure only) ===\n")
print(summary(seu@meta.data[, score_cols]))

# ---- Cells per patient, at the analysis threshold ----------------------
cov <- seu@meta.data %>%
  as_tibble() %>%
  count(sampleID, pathological_response, name = "n_tam") %>%
  arrange(n_tam)

cat("\n=== patients retained at threshold 10 ===\n")
print(cov %>% filter(n_tam >= 10) %>% count(pathological_response))
cat("\n=== patients below threshold ===\n")
print(cov %>% filter(n_tam < 10))

# ---- Export ------------------------------------------------------------
scores <- seu@meta.data %>%
  rownames_to_column("cellID") %>%
  select(cellID, sampleID, pathological_response, annotation_final,
         all_of(score_cols))

write_csv(scores, file.path(OUT_TAB, "05_LUSC_TAM_pure_scores_percell.csv"))
saveRDS(seu, OUT_OBJ)

message("\nScript05 (rescore on pure, ", HISTOLOGY, ") complete. ", Sys.time())