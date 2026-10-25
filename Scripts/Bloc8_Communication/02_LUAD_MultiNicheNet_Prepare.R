# ============================================================
# Bloc8_Communication, Script02_LUAD_MultiNicheNet_Prepare.R
# Purpose: build the SingleCellExperiment input for MultiNicheNet and run
# the cell type abundance and gene filtering steps of the official
# vignette (basic_analysis_steps_MISC), before any differential analysis.
#
# Cell groups (three, kept separate):
#   - CD8.TEX         : ProjecTILs functional.cluster == "CD8.TEX" (Bloc6)
#   - Macrophage_pure : scGate_TAM_class == "Macrophage_pure" (Bloc7)
#   - Monocyte_pure   : scGate_TAM_class == "Monocyte_pure" (Bloc7)
# Ambiguous_double_pure and scGate_unresolved cells are not used.
# A cell present in both the CD8 and the TAM object is excluded.
#
# Genes: the union of the genes of the two objects is kept (in LUAD the
# CD8 object holds 34115 genes, the TAM object 26636). Genes absent from
# the TAM object are set to zero in TAM cells. The gene filtering step of
# MultiNicheNet then removes genes that are not expressed in a cell type.
#
# Parameters: defaults from the official vignette are kept.
#   min_cells = 10, fraction_cutoff = 0.05, min_sample_prop = 0.50
#
# Input : Objects/03_LUAD_seu_projectils.rds
#         Objects/04b_LUAD_seu_TAM_classified.rds
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
f_tex <- here("Objects", "03_LUAD_seu_projectils.rds")
f_tam <- here("Objects", "04b_LUAD_seu_TAM_classified.rds")
stopifnot(file.exists(f_tex), file.exists(f_tam))

message("Loading CD8 object... ", Sys.time())
seu_tex <- readRDS(f_tex)
message("Loading TAM object... ", Sys.time())
seu_tam <- readRDS(f_tam)

stopifnot(all(c("sampleID", "pathological_response", "functional.cluster") %in%
                colnames(seu_tex@meta.data)))
stopifnot(all(c("sampleID", "pathological_response", "scGate_TAM_class") %in%
                colnames(seu_tam@meta.data)))

# A cell found in both objects cannot carry two labels: it is excluded
shared_cells <- intersect(colnames(seu_tex), colnames(seu_tam))
cat("\ncells present in both objects (excluded):", length(shared_cells), "\n")

tex_cells  <- setdiff(colnames(seu_tex)[seu_tex$functional.cluster %in% "CD8.TEX"],
                      shared_cells)
mac_cells  <- setdiff(colnames(seu_tam)[seu_tam$scGate_TAM_class %in% "Macrophage_pure"],
                      shared_cells)
mono_cells <- setdiff(colnames(seu_tam)[seu_tam$scGate_TAM_class %in% "Monocyte_pure"],
                      shared_cells)

cat("CD8.TEX:", length(tex_cells), "\n")
cat("Macrophage_pure:", length(mac_cells), "\n")
cat("Monocyte_pure:", length(mono_cells), "\n")

# Raw counts of the RNA assay
counts_tex <- GetAssayData(seu_tex, assay = "RNA", layer = "counts")[, tex_cells]
counts_tam <- GetAssayData(seu_tam, assay = "RNA", layer = "counts")[, c(mac_cells, mono_cells)]

# Union of genes: pad the TAM matrix with zero rows for the genes it lacks
genes_all <- rownames(counts_tex)
stopifnot(all(rownames(counts_tam) %in% genes_all))

s <- Matrix::summary(counts_tam)
counts_tam_full <- Matrix::sparseMatrix(
  i = match(rownames(counts_tam), genes_all)[s$i], j = s$j, x = s$x,
  dims = c(length(genes_all), ncol(counts_tam)),
  dimnames = list(genes_all, colnames(counts_tam))
)
counts_all <- cbind(counts_tex, counts_tam_full)

# Cell metadata: one label per cell
mk_meta <- function(seu, cells, label) {
  data.frame(sampleID = as.character(seu@meta.data[cells, "sampleID"]),
             pathological_response = as.character(seu@meta.data[cells, "pathological_response"]),
             mnn_celltype = label,
             row.names = cells, stringsAsFactors = FALSE)
}
meta <- rbind(mk_meta(seu_tex, tex_cells,  "CD8.TEX"),
              mk_meta(seu_tam, mac_cells,  "Macrophage_pure"),
              mk_meta(seu_tam, mono_cells, "Monocyte_pure"))
meta <- meta[colnames(counts_all), ]

# Guards
stopifnot(!anyDuplicated(colnames(counts_all)),
          nrow(meta) == ncol(counts_all),
          !anyNA(meta$sampleID),
          all(meta$pathological_response %in% c("MPR", "non-MPR")))

rm(seu_tex, seu_tam, counts_tex, counts_tam, counts_tam_full); gc()

# Seurat object with the same normalisation as the rest of the project
seu_sub <- CreateSeuratObject(counts = counts_all, meta.data = meta, assay = "RNA")
seu_sub <- NormalizeData(seu_sub, normalization.method = "LogNormalize",
                         scale.factor = 10000, verbose = FALSE)
rm(counts_all); gc()

cat("\n=== final composition (cells) ===\n")
print(table(seu_sub$mnn_celltype, seu_sub$pathological_response))
cat("genes:", nrow(seu_sub), "| patients:", length(unique(seu_sub$sampleID)), "\n")

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

# Both contrast directions are declared explicitly, as in the vignette.
contrasts_oi <- c("'MPR-non.MPR','non.MPR-MPR'")
contrast_tbl <- tibble(contrast = c("MPR-non.MPR", "non.MPR-MPR"),
                       group    = c("MPR", "non.MPR"))

cat("\n=== check that contrast labels match the group labels ===\n")
print(unique(SummarizedExperiment::colData(sce)[, group_id]))

# ---- 5. Cell type filtering (vignette step 1) --------------------------
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

total_nr_conditions <- SummarizedExperiment::colData(sce)[, group_id] %>%
  unique() %>% length()

absent_celltypes <- abundance_df_summarized %>%
  filter(samples_present < 2) %>%
  group_by(celltype_id) %>%
  count() %>%
  filter(n == total_nr_conditions) %>%
  pull(celltype_id)

cat("\ncondition-specific cell types:",
    ifelse(length(condition_specific_celltypes) == 0, "none",
           paste(condition_specific_celltypes, collapse = ", ")), "\n")
cat("absent cell types:",
    ifelse(length(absent_celltypes) == 0, "none",
           paste(absent_celltypes, collapse = ", ")), "\n")

# Vignette default: condition-specific cell types are not analysed
analyse_condition_specific_celltypes <- FALSE

if (analyse_condition_specific_celltypes == TRUE) {
  senders_oi   <- senders_oi   %>% setdiff(absent_celltypes)
  receivers_oi <- receivers_oi %>% setdiff(absent_celltypes)
} else {
  senders_oi   <- senders_oi   %>% setdiff(union(absent_celltypes, condition_specific_celltypes))
  receivers_oi <- receivers_oi %>% setdiff(union(absent_celltypes, condition_specific_celltypes))
}

sce <- sce[, SummarizedExperiment::colData(sce)[, celltype_id] %in%
             c(senders_oi, receivers_oi)]

# Guard: the three cell groups must all be kept
stopifnot(setequal(senders_oi, c("CD8.TEX", "Macrophage_pure", "Monocyte_pure")))
cat("\ncell types kept:", paste(senders_oi, collapse = ", "), "\n")

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