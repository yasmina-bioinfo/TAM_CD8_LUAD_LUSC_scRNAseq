# ============================================================
# Bloc8_Communication — Script05_ExpressionChecks.R
# DEPRECATED: not part of the final analysis, do not run.
# Written to interpret mis-oriented pairs and the TGFB1 axis in the
# 200-interaction literature check, an approach later replaced by the
# search for known axes in the full interaction tables (done outside the scripts, 
# from the full tables of Script04).
# Kept for traceability only. Output was console-only; no result in
# the README or the preprint depends on it.
# ============================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr)
})

stop("Script05_ExpressionChecks.R is deprecated and is not part of the final analysis. Do not run.")

# ---- helper ------------------------------------------------------------
get_freq <- function(histology, genes, by_group = FALSE) {
  obj <- readRDS(here("Objects",
                      paste0("08_", histology, "_sce_multinichenet.rds")))
  grouping <- if (by_group) c("celltype", "gene", "group") else c("celltype", "gene")
  res <- obj$frq_list$frq_df %>%
    filter(gene %in% genes) %>%
    group_by(across(all_of(grouping))) %>%
    summarise(pct = round(100 * mean(fraction_sample), 1), .groups = "drop") %>%
    mutate(histology = histology)
  rm(obj); gc()
  res
}

# ---- 1. Myeloid receptors from mis-oriented pairs ----------------------
RECEPTORS <- c("SIGLEC7", "SIGLEC9", "SIGLEC10", "SIRPA", "LILRB2",
               "NECTIN2", "CD47", "TIGIT", "PTPRO", "IL4R")

res1 <- bind_rows(get_freq("LUAD", RECEPTORS),
                  get_freq("LUSC", RECEPTORS))

cat("\n=== myeloid receptors: mean % of cells expressing, by histology ===\n")
print(as.data.frame(
  res1 %>%
    pivot_wider(names_from = histology, values_from = pct) %>%
    arrange(celltype, gene)
))

# ---- 2. TGFB1 and its receptors, by response group ---------------------
TGFB <- c("TGFB1", "TGFBR1", "TGFBR2", "ENG", "ACVRL1",
          "ITGAV", "ITGB1", "ITGB8", "SDC2", "APP")

res2 <- bind_rows(get_freq("LUAD", TGFB, by_group = TRUE),
                  get_freq("LUSC", TGFB, by_group = TRUE))

cat("\n=== TGFB1 axis: % of cells expressing, by response group ===\n")
print(as.data.frame(
  res2 %>%
    pivot_wider(names_from = group, values_from = pct) %>%
    arrange(histology, gene, celltype)
))

cat("\nScript05 (expression checks) complete.\n")