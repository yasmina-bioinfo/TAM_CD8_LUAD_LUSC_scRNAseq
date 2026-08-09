# ============================================================
# Bloc8_Communication — Script02_LUAD_MultiNicheNet_Prepare.R
# Purpose: build the SingleCellExperiment input for MultiNicheNet and run
# the cell-type abundance and gene filtering diagnostics, before any
# differential analysis.
#
# Scope: three cell groups, chosen because they are the only compartments
# whose biology has been characterised in this cohort. Widening to NK,
# Tregs, DC or B cells would produce interactions that could not be
# interpreted against anything established here.
#   - Macrophages   : full lineage_TAM, comprising scGate-confirmed cells
#                     and transitional cells along the monocyte-to-
#                     macrophage continuum. The restriction to
#                     scGate-confirmed cells applied in Bloc7 followed
#                     from the Cheng signatures describing macrophage
#                     states; that constraint does not apply here, and
#                     transitional cells express ligands and receptors
#                     like any other
#   - CD8.TEX       : ProjecTILs functional state, Bloc6
#   - Monocytes     : classical monocytes, a distinct population rather
#                     than a stage of the above continuum. Excluded from
#                     Bloc7 scoring since the signatures do not describe
#                     undifferentiated cells, but retained here as a real
#                     source and target of signalling
#
# Parameters: defaults from the official vignette are kept unless stated.
#   min_cells = 10        : the package default, and the threshold used
#                           throughout this project (muscat convention)
#   fraction_cutoff = 0.05: package default; authors advise against >0.10
#   min_sample_prop = 0.50: package default
#
# Input : Objects/07_LUAD_seu_TAM_scored.rds       (macrophages)
#         Objects/06_LUAD_seu_functional_scores.rds (CD8.TEX)
#         Objects/05_LUAD_seu_annotated.rds         (monocytes)
# Output: Objects/08_LUAD_sce_multinichenet.rds
# ============================================================

suppressPackageStartupMessages({
  library(here); library(Seurat); library(SingleCellExperiment)
  library(nichenetr); library(multinichenetr)
  library(dplyr); library(tibble); library(ggplot2)
})

set.seed(42)

HISTOLOGY <- "LUAD"
OUT_OBJ   <- here("Objects", "08_LUAD_sce_multinichenet.rds")
OUT_FIG   <- here("Results", "Bloc8_Communication", "Figures")
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

MIN_CELLS       <- 10
FRACTION_CUTOFF <- 0.05
MIN_SAMPLE_PROP <- 0.50

# ---- 1. Assemble the three cell groups ---------------------------------
message("Loading macrophages (full lineage_TAM)... ", Sys.time())
seu_mac <- readRDS(here("Objects", "07_LUAD_seu_TAM_scored.rds"))
mac_cells <- colnames(seu_mac)
rm(seu_mac); gc()

message("Loading CD8 object... ")
seu_cd8 <- readRDS(here("Objects", "06_LUAD_seu_functional_scores.rds"))
tex_cells <- colnames(seu_cd8)[seu_cd8$functional.cluster == "CD8.TEX"]
rm(seu_cd8); gc()

message("Loading full annotated object... ")
seu <- readRDS(here("Objects", "05_LUAD_seu_annotated.rds"))
mono_cells <- colnames(seu)[seu$lineage_TAM_related & !seu$lineage_TAM]

cat("\n=== cells per group ===\n")
cat("Macrophages (full lineage_TAM):", length(mac_cells), "\n")
cat("CD8.TEX:", length(tex_cells), "\n")
cat("Classical monocytes:", length(mono_cells), "\n")

# Build the cell type column
seu$mnn_celltype <- NA_character_
seu$mnn_celltype[colnames(seu) %in% mac_cells]  <- "Macrophage"
seu$mnn_celltype[colnames(seu) %in% tex_cells]  <- "CD8_TEX"
seu$mnn_celltype[colnames(seu) %in% mono_cells] <- "Monocyte"

seu_sub <- subset(seu, cells = colnames(seu)[!is.na(seu$mnn_celltype)])
rm(seu); gc()

cat("\n=== final composition ===\n")
print(table(seu_sub$mnn_celltype, seu_sub$pathological_response))

# ---- 2. Syntactically valid names (required by MultiNicheNet) ----------
seu_sub$mnn_sample   <- make.names(seu_sub$sampleID)
seu_sub$mnn_group    <- make.names(seu_sub$pathological_response)
seu_sub$mnn_celltype <- make.names(seu_sub$mnn_celltype)

cat("\n=== group labels after make.names() ===\n")
print(table(seu_sub$mnn_group))

# ---- 3. Convert to SingleCellExperiment --------------------------------
message("\nConverting to SingleCellExperiment... ", Sys.time())
sce <- Seurat::as.SingleCellExperiment(seu_sub, assay = "RNA")
rm(seu_sub); gc()

# Update gene symbols to match the NicheNet v2 networks, then make valid
sce <- alias_to_symbol_SCE(sce, "human") %>% makenames_SCE()
cat("genes after symbol update:", nrow(sce), "\n")

# ---- 4. Analysis settings ----------------------------------------------
sample_id   <- "mnn_sample"
group_id    <- "mnn_group"
celltype_id <- "mnn_celltype"
covariates  <- NA
batches     <- NA

senders_oi   <- SummarizedExperiment::colData(sce)[, celltype_id] %>% unique()
receivers_oi <- senders_oi

# Both directions are declared explicitly: MultiNicheNet does not test the
# reverse contrast automatically, so a single direction would only reveal
# interactions specific to one group.
contrasts_oi <- c("'MPR-non.MPR','non.MPR-MPR'")
contrast_tbl <- tibble(contrast = c("MPR-non.MPR", "non.MPR-MPR"),
                       group    = c("MPR", "non.MPR"))

cat("\n=== check that contrast labels match the group labels ===\n")
print(unique(SummarizedExperiment::colData(sce)[, group_id]))

# ---- 5. Cell type abundance diagnostics --------------------------------
message("\nComputing abundance info... ", Sys.time())
abundance_info <- get_abundance_info(
  sce = sce,
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id,
  min_cells = MIN_CELLS,
  senders_oi = senders_oi, receivers_oi = receivers_oi,
  batches = batches
)

p_abund <- abundance_info$abund_plot_sample
ggsave(file.path(OUT_FIG, "02_LUAD_abundance_per_sample.png"),
       p_abund, width = 10, height = 6, dpi = 300, bg = "white")

# Which cell types are absent or condition-specific?
abundance_df_summarized <- abundance_info$abundance_data %>%
  mutate(keep = as.logical(keep)) %>%
  group_by(group_id, celltype_id) %>%
  summarise(samples_present = sum(keep), .groups = "drop")

cat("\n=== samples passing min_cells, per group and cell type ===\n")
print(abundance_df_summarized, n = Inf)

celltypes_absent_one_condition <- abundance_df_summarized %>%
  filter(samples_present == 0) %>% pull(celltype_id) %>% unique()
celltypes_present_one_condition <- abundance_df_summarized %>%
  filter(samples_present >= 2) %>% pull(celltype_id) %>% unique()

condition_specific_celltypes <- intersect(celltypes_absent_one_condition,
                                          celltypes_present_one_condition)

cat("\ncondition-specific cell types:", 
    ifelse(length(condition_specific_celltypes) == 0, "none",
           paste(condition_specific_celltypes, collapse = ", ")), "\n")

# ---- 6. Gene filtering --------------------------------------------------
message("\nComputing expressed genes... ", Sys.time())
frq_list <- get_frac_exprs(
  sce = sce,
  sample_id = sample_id, celltype_id = celltype_id, group_id = group_id,
  batches = batches,
  min_cells = MIN_CELLS,
  fraction_cutoff = FRACTION_CUTOFF,
  min_sample_prop = MIN_SAMPLE_PROP
)

genes_oi <- frq_list$expressed_df %>%
  filter(expressed == TRUE) %>% pull(gene) %>% unique()
sce <- sce[genes_oi, ]

cat("\ngenes retained:", length(genes_oi), "\n")

# ---- 7. Save ------------------------------------------------------------
saveRDS(list(sce = sce,
             frq_list = frq_list,
             abundance_info = abundance_info,
             settings = list(sample_id = sample_id, group_id = group_id,
                             celltype_id = celltype_id,
                             contrasts_oi = contrasts_oi,
                             contrast_tbl = contrast_tbl,
                             min_cells = MIN_CELLS,
                             fraction_cutoff = FRACTION_CUTOFF,
                             min_sample_prop = MIN_SAMPLE_PROP)),
        OUT_OBJ)

message("\nScript02 (prepare, ", HISTOLOGY, ") complete. ", Sys.time())