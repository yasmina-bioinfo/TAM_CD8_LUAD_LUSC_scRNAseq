# ============================================================
# Bloc8_Communication — Script05_ExpressionChecks.R
# Purpose: two expression checks used to interpret the prioritised
# interactions.
#
# 1. Myeloid receptors appearing in mis-oriented pairs: are they more
#    expressed in LUSC than in LUAD? This tests whether their higher
#    frequency in the LUSC ranking follows from expression levels.
# 2. TGFB1 and its receptors by response group: is the ligand restricted
#    to the compartment the inference assigns as sender?
#
# Input : Objects/08_{HIST}_sce_multinichenet.rds
# Output: console only
# ============================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr)
})

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