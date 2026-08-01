# ============================================================
# Bloc8_Communication — Script03_LUAD_MultiNicheNet_Analysis.R
# Purpose: run the MultiNicheNet core analysis on the object prepared in
# Script02: pseudobulk differential expression, geneset ratio check,
# ligand activity prediction, multi-criteria prioritisation, and
# ligand-receptor / target expression correlation.
#
# Input : Objects/08_LUAD_sce_multinichenet.rds
#         References/nichenet_lr_network_all_human.rds
#         References/nichenet_ligand_target_matrix_human.rds
# Output: Objects/08_LUAD_multinichenet_output.rds
#         Results/Bloc8_Communication/Tables/03_LUAD_*.csv
#
# p_val_adj = FALSE here. The package authors apply the p-value cutoff to
# raw p-values by default, since most multi-sample single-cell datasets
# have few samples per group, and recommend adjusted p-values only when
# the smallest group has at least 20 samples. LUAD has 9 MPR patients,
# so raw p-values are used. LUSC, with 29 to 34, uses TRUE.
# ============================================================

suppressPackageStartupMessages({
  library(here); library(SingleCellExperiment)
  library(nichenetr); library(multinichenetr)
  library(dplyr); library(tibble); library(readr)
})

set.seed(42)

HISTOLOGY <- "LUAD"
IN_OBJ    <- here("Objects", "08_LUAD_sce_multinichenet.rds")
OUT_OBJ   <- here("Objects", "08_LUAD_multinichenet_output.rds")
OUT_TAB   <- here("Results", "Bloc8_Communication", "Tables")
OUT_FIG <- here("Results", "Bloc8_Communication", "Figures")
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

# ---- Parameters ---------------------------------------------------------
LOGFC_THRESHOLD <- 0.50
P_VAL_THRESHOLD <- 0.05
P_VAL_ADJ       <- FALSE     # LUAD: smallest group = 9 patients
TOP_N_TARGET    <- 250
SCENARIO        <- "regular"
LIGAND_ACT_DOWN <- FALSE
N_CORES         <- 4

# ---- 1. Load prepared object and reference networks --------------------
message("Loading prepared object... ", Sys.time())
prep <- readRDS(IN_OBJ)
sce            <- prep$sce
frq_list       <- prep$frq_list
abundance_info <- prep$abundance_info
st             <- prep$settings

sample_id   <- st$sample_id
group_id    <- st$group_id
celltype_id <- st$celltype_id
contrasts_oi <- st$contrasts_oi
contrast_tbl <- st$contrast_tbl
MIN_CELLS       <- st$min_cells
FRACTION_CUTOFF <- st$fraction_cutoff

batches    <- NA
covariates <- NA

senders_oi   <- SummarizedExperiment::colData(sce)[, celltype_id] %>% unique()
receivers_oi <- senders_oi

lr_network_all <- readRDS(here("References", "nichenet_lr_network_all_human.rds"))
lr_network     <- readRDS(here("References", "nichenet_lr_network_human.rds"))
ligand_target_matrix <- readRDS(here("References",
                                     "nichenet_ligand_target_matrix_human.rds"))

# Restrict the network to ligands present in the matrix, as in the vignette
lr_network <- lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
ligand_target_matrix <- ligand_target_matrix[, unique(lr_network$ligand)]

cat("cell groups:", paste(senders_oi, collapse = ", "), "\n")
cat("ligand-receptor pairs after restriction:", nrow(lr_network), "\n")

# ---- 2. Pseudobulk expression info -------------------------------------
message("\nComputing pseudobulk expression info... ", Sys.time())
abundance_expression_info <- process_abundance_expression_info(
  sce = sce,
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id,
  min_cells = MIN_CELLS,
  senders_oi = senders_oi, receivers_oi = receivers_oi,
  lr_network = lr_network,
  batches = batches,
  frq_list = frq_list,
  abundance_info = abundance_info
)

# ---- 3. Differential expression: pseudobulk + EdgeR --------------------
message("\nRunning differential expression... ", Sys.time())
DE_info <- get_DE_info(
  sce = sce,
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id,
  batches = batches, covariates = covariates,
  contrasts_oi = contrasts_oi,
  min_cells = MIN_CELLS,
  expressed_df = frq_list$expressed_df
)

celltype_de <- DE_info$celltype_de$de_output_tidy

cat("\n=== top DE genes by logFC ===\n")
print(celltype_de %>% arrange(-logFC) %>% head(10))

# p-value distributions: should be roughly uniform with a peak below 0.05.
# An irregular shape would call for empirical p-values instead.
ggplot2::ggsave(
  file.path(OUT_FIG, "03_LUAD_pvalue_histograms.png"),
  DE_info$hist_pvals, width = 10, height = 6, dpi = 300, bg = "white"
)

# ---- 4. Geneset / background ratio check -------------------------------
# The ligand activity step is an enrichment test: the ratio of DE genes to
# expressed genes should lie between 1/200 and 1/10. Outside that range,
# thresholds must be adjusted and the adjustment documented.
message("\nChecking geneset/background ratios... ")

geneset_assessment <- contrast_tbl$contrast %>%
  lapply(process_geneset_data,
         celltype_de, LOGFC_THRESHOLD, P_VAL_ADJ, P_VAL_THRESHOLD) %>%
  bind_rows()

cat("\n=== geneset assessment (p_val_adj =", P_VAL_ADJ, ") ===\n")
print(geneset_assessment, n = Inf)
write_csv(geneset_assessment,
          file.path(OUT_TAB, "03_LUAD_geneset_assessment.csv"))

# For comparison only, not used unless the above is out of range
geneset_assessment_adj <- contrast_tbl$contrast %>%
  lapply(process_geneset_data,
         celltype_de, LOGFC_THRESHOLD, TRUE, P_VAL_THRESHOLD) %>%
  bind_rows()
cat("\n=== same, with adjusted p-values (for comparison) ===\n")
print(geneset_assessment_adj, n = Inf)

# ---- 5. Combine sender and receiver DE ---------------------------------
sender_receiver_de <- combine_sender_receiver_de(
  sender_de = celltype_de, receiver_de = celltype_de,
  senders_oi = senders_oi, receivers_oi = receivers_oi,
  lr_network = lr_network
)

# ---- 6. Ligand activity prediction -------------------------------------
message("\nPredicting ligand activities... ", Sys.time())
ligand_activities_targets_DEgenes <- suppressMessages(suppressWarnings(
  get_ligand_activities_targets_DEgenes(
    receiver_de = celltype_de,
    receivers_oi = intersect(receivers_oi, unique(celltype_de$cluster_id)),
    ligand_target_matrix = ligand_target_matrix,
    logFC_threshold = LOGFC_THRESHOLD,
    p_val_threshold = P_VAL_THRESHOLD,
    p_val_adj = P_VAL_ADJ,
    top_n_target = TOP_N_TARGET,
    verbose = TRUE,
    n.cores = N_CORES
  )
))

# ---- 7. Prioritisation --------------------------------------------------
message("\nGenerating prioritisation tables... ", Sys.time())
sender_receiver_tbl <- sender_receiver_de %>% distinct(sender, receiver)

metadata_combined <- SummarizedExperiment::colData(sce) %>% as_tibble()
grouping_tbl <- metadata_combined[, c(sample_id, group_id)] %>%
  as_tibble() %>% distinct()
colnames(grouping_tbl) <- c("sample", "group")

prioritization_tables <- suppressMessages(generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  scenario = SCENARIO,
  fraction_cutoff = FRACTION_CUTOFF,
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender,
  ligand_activity_down = LIGAND_ACT_DOWN
))

cat("\n=== top 20 prioritised interactions ===\n")
print(prioritization_tables$group_prioritization_tbl %>%
        select(group, sender, receiver, ligand, receptor,
               prioritization_score) %>% head(20))

cat("\n=== number of interactions per group ===\n")
print(table(prioritization_tables$group_prioritization_tbl$group))

# ---- 8. Across-sample LR / target expression correlation ---------------
message("\nComputing LR-target expression correlation... ", Sys.time())
lr_target_prior_cor <- lr_target_prior_cor_inference(
  receivers_oi = unique(prioritization_tables$group_prioritization_tbl$receiver),
  abundance_expression_info = abundance_expression_info,
  celltype_de = celltype_de,
  grouping_tbl = grouping_tbl,
  prioritization_tables = prioritization_tables,
  ligand_target_matrix = ligand_target_matrix,
  logFC_threshold = LOGFC_THRESHOLD,
  p_val_threshold = P_VAL_THRESHOLD,
  p_val_adj = P_VAL_ADJ
)

# ---- 9. Save ------------------------------------------------------------
multinichenet_output <- list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = lr_target_prior_cor
)
multinichenet_output <- make_lite_output(multinichenet_output)

saveRDS(multinichenet_output, OUT_OBJ)

write_csv(
  prioritization_tables$group_prioritization_tbl %>%
    select(where(~ !is.list(.))),
  file.path(OUT_TAB, "03_LUAD_group_prioritization.csv")
)

message("\nScript03 (analysis, ", HISTOLOGY, ") complete. ", Sys.time())