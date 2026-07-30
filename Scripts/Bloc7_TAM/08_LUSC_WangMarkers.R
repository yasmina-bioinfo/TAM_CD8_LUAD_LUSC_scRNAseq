# ============================================================
# Bloc7_TAM — Script08_WangMarkers_LUSC.R
# Purpose: test whether the four macrophage subtype markers reported by
# Wang et al. (2022, Signal Transduct Target Ther) structure the
# macrophage compartment in this treatment-exposed cohort.
#
# Context: Wang et al. clustered 11,016 macrophages from 19 TREATMENT-NAIVE
# NSCLC patients into four subtypes named by marker gene: FABP4, FCN1,
# SPP1, SELENOP. FABP4-Mphi dominated LUSC (35.6% vs 7.3%), SPP1-Mphi
# dominated LUSC (44.6% vs 39.9%).
#
# Caveats to keep in mind:
#   - Their subtypes come from unsupervised clustering, not from a
#     threshold on four genes. This is therefore NOT a reproduction of
#     their method, only a test of whether these markers structure the
#     compartment here.
#   - They state that marker expression was not exclusive to a single
#     cluster, so separation is imperfect even in their own data.
#   - Their patients are untreated; this cohort is post anti-PD-1.
#
# Input : Objects/07_LUSC_seu_TAM_pure_scored.rds
# Output: console + Results/Bloc7_TAM/Tables/08_LUSC_WangMarkers.csv
# ============================================================

suppressPackageStartupMessages({
  library(here); library(Seurat); library(dplyr)
  library(readr); library(tibble)
})

HISTOLOGY <- "LUSC"
seu <- readRDS(here("Objects", "07_LUSC_seu_TAM_pure_scored.rds"))
OUT_TAB <- here("Results", "Bloc7_TAM", "Tables")

MK <- c("FABP4","FCN1","SPP1","SELENOP")
stopifnot(all(MK %in% rownames(seu)))

# ---- 1. Overall expression level of the four markers -------------------
cat("=== mean expression, all pure macrophages (", HISTOLOGY, ") ===\n")
print(round(AverageExpression(seu, features = MK, assays = "RNA",
                              layer = "data")$RNA, 2))

# ---- 2. Fraction of cells expressing each marker ----------------------
# A cell is counted as positive if the gene has any non-zero count.
expr <- GetAssayData(seu, assay = "RNA", layer = "data")[MK, ]
pos  <- as.matrix(expr > 0)

cat("\n=== % of macrophages positive for each marker ===\n")
print(round(100 * rowMeans(pos), 1))

# ---- 3. Do the markers separate cells, or overlap? ---------------------
# Assign each cell to the marker with the highest expression, but only
# if at least one marker is detected.
top_mk <- apply(as.matrix(expr), 2, function(x) {
  if (max(x) == 0) NA_character_ else MK[which.max(x)]
})
seu$wang_top <- top_mk

cat("\n=== cells assigned by highest-expressing marker ===\n")
print(table(seu$wang_top, useNA = "ifany"))
cat("\nproportions (%):\n")
print(round(100 * prop.table(table(seu$wang_top)), 1))

cat("\n=== co-expression: number of markers detected per cell ===\n")
print(table(colSums(pos)))

# ---- 4. Per-patient assignment, by response group ---------------------
res <- seu@meta.data %>%
  as_tibble() %>%
  filter(!is.na(wang_top)) %>%
  count(sampleID, pathological_response, wang_top) %>%
  group_by(sampleID, pathological_response) %>%
  mutate(pct = 100 * n / sum(n), n_assigned = sum(n)) %>%
  ungroup()

cat("\n=== median % per marker group, by response ===\n")
res %>%
  filter(n_assigned >= 10) %>%
  group_by(pathological_response, wang_top) %>%
  summarise(median_pct = round(median(pct), 1),
            n_patients = n(), .groups = "drop") %>%
  print(n = Inf)

write_csv(res, file.path(OUT_TAB, "08_LUSC_WangMarkers.csv"))
cat("\nScript08 (Wang markers,", HISTOLOGY, ") complete.\n")