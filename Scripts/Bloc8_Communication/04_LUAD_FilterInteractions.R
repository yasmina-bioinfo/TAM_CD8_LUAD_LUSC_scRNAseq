# ============================================================
# Bloc8_Communication — Script04_FilterInteractions_LUAD.R
# Purpose: restrict the prioritised interactions to those that can be
# interpreted, by removing autocrine pairs.
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
N_TOP     <- 50   # interactions retained per group

OUT_TAB <- here("Results", "Bloc8_Communication", "Tables")

# ---- 1. Load the MultiNicheNet output ---------------------------------
output <- readRDS(here("Objects",
                       paste0("08_", HISTOLOGY, "_multinichenet_output.rds")))

prior <- output$prioritization_tables$group_prioritization_tbl

# ---- 2. Filter -------------------------------------------------------
filtered <- prior %>%
  filter(sender != receiver)

cat("all prioritised      :", nrow(prior), "\n")
cat("non-autocrine        :", sum(prior$sender != prior$receiver), "\n")

top <- filtered %>%
  group_by(group) %>%
  arrange(desc(prioritization_score)) %>%
  slice_head(n = N_TOP) %>%
  ungroup()

cat("\n=== retained per group ===\n")
print(table(top$group))

cat("\n=== top interactions, non-MPR ===\n")
print(as.data.frame(top %>% filter(group == "non.MPR") %>%
                      select(sender, receiver, ligand, receptor,
                             prioritization_score) %>%
                      slice_head(n = 50)))

cat("\n=== top interactions, MPR ===\n")
print(as.data.frame(top %>% filter(group == "MPR") %>%
                      select(sender, receiver, ligand, receptor,
                             prioritization_score) %>%
                      slice_head(n = 50)))

# ---- 3. Direction of the dialogue -------------------------------------
cat("\n=== interactions by direction ===\n")
print(as.data.frame(
  top %>%
    dplyr::count(group, sender, receiver) %>%
    tidyr::pivot_wider(names_from = group, values_from = n,
                       values_fill = 0) %>%
    arrange(desc(non.MPR))
))


write_csv(top, file.path(OUT_TAB,
                         paste0("04_", HISTOLOGY, "_filtered_interactions.csv")))