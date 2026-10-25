# ============================================================
# Bloc8_Communication — Script04_FilterInteractions_LUAD.R
# Purpose: restriction to CD8.TEX–myeloid interactions
#
# Rationale: sender and receiver being the same cell type cannot be
# resolved from dissociated data.
#
# Note: an earlier version also required at least two sources per
# ligand-receptor pair in the reference network. That filter was dropped:
# all checkpoint pairs, including CD274-PDCD1 and CD80-CTLA4, carry a
# single source, so the count reflects how many databases list a pair
# rather than how well established it is.
# ============================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(readr); library(multinichenetr); library(tidyr)
})

HISTOLOGY <- "LUAD"

OUT_TAB <- here("Results", "Bloc8_Communication", "Tables")

# ---- 1. Load the MultiNicheNet output ---------------------------------
f_out <- here("Objects", paste0("08_", HISTOLOGY, "_multinichenet_output.rds"))
stopifnot(file.exists(f_out))
output <- readRDS(f_out)

stopifnot(all(c("group", "top_group", "sender", "receiver", "ligand", "receptor",
                "id", "prioritization_score",
                "fraction_expressing_ligand_receptor") %in%
                colnames(output$prioritization_tables$group_prioritization_tbl)),
          setequal(unique(output$prioritization_tables$group_prioritization_tbl$group),
                   c("MPR", "non.MPR")))

# ---- 2. Select and filter ----------------------------------------------
all_pairs <- get_top_n_lr_pairs(output$prioritization_tables, top_n = 10000)

# Keep only CD8.TEX <-> myeloid interactions (research question)
filtered <- all_pairs %>% ungroup() %>%
  filter(sender != receiver,
         sender == "CD8.TEX" | receiver == "CD8.TEX")

full <- filtered %>%
  group_by(group) %>%
  arrange(desc(prioritization_score), .by_group = TRUE) %>%
  mutate(rank_non_autocrine = row_number()) %>%
  ungroup()

cat("all prioritised      :", nrow(all_pairs), "\n")
cat("CD8.TEX-myeloid      :", nrow(full), "\n")

# ---- 3. Direction of the dialogue -------------------------------------
cat("\n=== interactions by direction ===\n")
print(as.data.frame(
  full %>%
    dplyr::count(group, sender, receiver) %>%
    tidyr::pivot_wider(names_from = group, values_from = n, values_fill = 0) %>%
    arrange(desc(non.MPR))
))

# ---- 4. Export ----------------------------------------------------------
write_csv(full, file.path(OUT_TAB,
                          paste0("04_", HISTOLOGY, "_filtered_interactions_FULL.csv")))