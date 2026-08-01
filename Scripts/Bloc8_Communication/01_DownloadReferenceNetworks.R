# ============================================================
# Bloc8_Communication — Script01_DownloadReferenceNetworks.R
# Purpose: download and cache the NicheNet v2 prior knowledge networks
# required by MultiNicheNet.
#
# Source: Zenodo, doi:10.5281/zenodo.7074291 (networks and matrices)
#         and record 10229222 for the annotated LR network used by
#         MultiNicheNet specifically.
#
# Note: MultiNicheNet uses lr_network_human_allInfo, which carries
# annotation columns absent from the plain lr_network file used by
# NicheNet alone. Gene names are made syntactically valid with
# make.names() as done in the official vignette, to avoid losing genes
# in downstream visualisations.
#
# Output: References/nichenet_lr_network_all_human.rds
#         References/nichenet_lr_network_human.rds
#         References/nichenet_ligand_target_matrix_human.rds
# ============================================================

suppressPackageStartupMessages({
  library(here); library(nichenetr); library(multinichenetr); library(dplyr)
})

REF_DIR <- here("References")
dir.create(REF_DIR, recursive = TRUE, showWarnings = FALSE)

options(timeout = max(3600, getOption("timeout")))

# ---- 1. Ligand-receptor network, annotated version -----------------------
f_lr_all <- file.path(REF_DIR, "nichenet_lr_network_all_human.rds")

if (!file.exists(f_lr_all)) {
  message("Downloading annotated LR network... ", Sys.time())
  lr_network_all <- readRDS(url(
    "https://zenodo.org/record/10229222/files/lr_network_human_allInfo_30112033.rds"
  )) %>%
    mutate(ligand   = convert_alias_to_symbols(ligand,   organism = "human"),
           receptor = convert_alias_to_symbols(receptor, organism = "human")) %>%
    mutate(ligand   = make.names(ligand),
           receptor = make.names(receptor))
  saveRDS(lr_network_all, f_lr_all)
} else {
  lr_network_all <- readRDS(f_lr_all)
}

lr_network <- lr_network_all %>% distinct(ligand, receptor)
saveRDS(lr_network, file.path(REF_DIR, "nichenet_lr_network_human.rds"))

cat("ligand-receptor pairs:", nrow(lr_network), "\n")
cat("unique ligands:", n_distinct(lr_network$ligand),
    "| unique receptors:", n_distinct(lr_network$receptor), "\n")

# ---- 2. Ligand-target matrix --------------------------------------------
f_lt <- file.path(REF_DIR, "nichenet_ligand_target_matrix_human.rds")

if (!file.exists(f_lt)) {
  message("Downloading ligand-target matrix... ", Sys.time())
  ligand_target_matrix <- readRDS(url(
    "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"
  ))
  colnames(ligand_target_matrix) <- make.names(colnames(ligand_target_matrix))
  rownames(ligand_target_matrix) <- make.names(rownames(ligand_target_matrix))
  saveRDS(ligand_target_matrix, f_lt)
} else {
  ligand_target_matrix <- readRDS(f_lt)
}

cat("ligand-target matrix:", nrow(ligand_target_matrix), "targets x",
    ncol(ligand_target_matrix), "ligands\n")

message("\nScript01 (reference networks) complete. ", Sys.time())