# ============================================================
# Bloc8_Communication — Script01_DownloadReferenceNetworks.R
# Purpose: download the NicheNet v2 prior knowledge networks required
# by MultiNicheNet and save them to References/.
#
# Source (two Zenodo records, version v1 of each; the NicheNet v2 model
# is the one named in their titles):
#   - doi:10.5281/zenodo.7074291: ligand-target matrix
#     (ligand_target_matrix_nsga2r_final.rds)
#   - doi:10.5281/zenodo.10229222: annotated ligand-receptor network
#     with OmniPath information (lr_network_human_allInfo_30112033.rds)
#
# Code adapted from the official MultiNicheNet tutorial
# (basic_analysis_steps_MISC, section "Load NicheNet's ligand-receptor
# network and ligand-target matrix"), human branch only. Gene symbols
# are converted from aliases to official symbols and made syntactically
# valid with make.names() in both objects. The ligand-receptor network
# is then restricted to ligands present in the matrix, and the matrix
# to those ligands.
#
# Output: References/nichenet_lr_network_all_human.rds (annotated, not filtered)
#         References/nichenet_lr_network_human.rds (ligand, receptor, filtered)
#         References/nichenet_ligand_target_matrix_human.rds (filtered)
# ============================================================

suppressPackageStartupMessages({
  library(here); library(nichenetr); library(multinichenetr); library(dplyr)
})

REF_DIR <- here("References")
dir.create(REF_DIR, recursive = TRUE, showWarnings = FALSE)

organism <- "human"
options(timeout = 3600)

# ---- Tutorial block (human branch only, mouse branch removed) ----
lr_network_all = readRDS(url(
  "https://zenodo.org/record/10229222/files/lr_network_human_allInfo_30112033.rds"
)) %>%
  mutate(
    ligand = convert_alias_to_symbols(ligand, organism = organism),
    receptor = convert_alias_to_symbols(receptor, organism = organism))

lr_network_all = lr_network_all %>%
  mutate(ligand = make.names(ligand), receptor = make.names(receptor))

lr_network = lr_network_all %>%
  distinct(ligand, receptor)

ligand_target_matrix = readRDS(url(
  "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"
))

colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>%
  convert_alias_to_symbols(organism = organism) %>% make.names()
rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>%
  convert_alias_to_symbols(organism = organism) %>% make.names()

lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]

# ---- End of tutorial block ----

saveRDS(lr_network_all, file.path(REF_DIR, "nichenet_lr_network_all_human.rds"))
saveRDS(lr_network, file.path(REF_DIR, "nichenet_lr_network_human.rds"))
saveRDS(ligand_target_matrix, file.path(REF_DIR, "nichenet_ligand_target_matrix_human.rds"))

cat("ligand-receptor pairs:", nrow(lr_network), "\n")
cat("unique ligands:", n_distinct(lr_network$ligand),
    "| unique receptors:", n_distinct(lr_network$receptor), "\n")
cat("ligand-target matrix:", nrow(ligand_target_matrix), "targets x",
    ncol(ligand_target_matrix), "ligands\n")

message("\nScript01 (reference networks) complete. ", Sys.time())