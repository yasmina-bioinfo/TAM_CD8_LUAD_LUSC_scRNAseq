# ============================================================
# Bloc7_TAM — Script13_CollecTRI_Myeloid_LUSC.R
# Purpose: infer transcription factor activity in the myeloid compartment,
# comparing MPR and non-MPR within each of two groups confirmed by scGate: 
# Macrophage_pure and Monocyte_pure.
#
# Status: EXPLORATORY, no directional prediction. This is not a further
# phenotypic framework but a different level of analysis: CollecTRI
# infers regulatory activity from target gene expression, independently
# of any cell type signature.
#
# Gene filtering: edgeR's filterByExpr (Chen et al., 2016), the
# recommended strategy for pseudobulk. It keeps genes whose counts per
# million exceed a threshold in a sufficient number of samples, the
# threshold being set from min.count and library sizes and the number of
# samples from the design. Using counts per million rather than raw
# counts avoids favouring patients contributing more cells. This replaces
# the earlier project-specific filters, which differed between Bloc6 and
# Bloc7 and were therefore not directly comparable.
#
# Input : Objects/05_LUSC_seu_annotated.rds
#         Objects/04b_LUSC_seu_TAM_classified.rds  (for the scGate column)
#         Data/interactions_collectri_9606.tsv.gz
# Output: Results/Bloc7_TAM/Tables/13_LUSC_CollecTRI_myeloid.csv
# ============================================================

suppressPackageStartupMessages({
  library(here); library(Seurat); library(dplyr); library(tibble)
  library(readr); library(stringr); library(DESeq2); library(decoupleR)
  library(edgeR)
})

set.seed(42)

HISTOLOGY <- "LUSC"
OUT_TAB   <- here("Results", "Bloc7_TAM", "Tables")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

MIN_CELLS    <- 10   # per patient, within each group
MIN_PATIENTS <- 4    # per response group, to run the contrast

# ---- 1. Build the two-group column -------------------------------------

message("Loading TAM-classified object (Bloc7)... ", Sys.time())
seu <- readRDS(here("Objects", paste0("04b_", HISTOLOGY, "_seu_TAM_classified.rds")))

seu$myeloid_group <- NA_character_
seu$myeloid_group[seu$scGate_TAM_class == "Macrophage_pure"] <- "Macrophage_pure"
seu$myeloid_group[seu$scGate_TAM_class == "Monocyte_pure"]   <- "Monocyte"

cat("\n=== cells per group and response ===\n")
print(table(seu$myeloid_group, seu$pathological_response, useNA = "no"))

# ---- 2. CollecTRI network ----------------------------------------------
collectri_raw <- read_tsv(here("Data", "interactions_collectri_9606.tsv.gz"),
                          show_col_types = FALSE)

net <- collectri_raw %>%
  filter(!str_detect(source, "COMPLEX")) %>%
  filter(!(is_stimulation == 1 & is_inhibition == 1)) %>%
  transmute(source = source_genesymbol,
            target = target_genesymbol,
            mor    = ifelse(is_inhibition == 1, -1, 1)) %>%
  distinct()

cat("\nCollecTRI: ", n_distinct(net$source), " TFs, ",
    nrow(net), " interactions\n", sep = "")

# ---- 3. One analysis per group -----------------------------------------
run_group <- function(group_name) {
  
  cat("\n\n########## ", group_name, " ##########\n")
  
  cells <- colnames(seu)[!is.na(seu$myeloid_group) &
                           seu$myeloid_group == group_name]
  md <- seu@meta.data[cells, c("sampleID", "pathological_response")]
  
  keep <- md %>%
    rownames_to_column("cell") %>%
    dplyr::count(sampleID, pathological_response) %>%
    filter(n >= MIN_CELLS)
  
  cat("patients retained: ")
  print(table(keep$pathological_response))
  
  if (min(table(keep$pathological_response)) < MIN_PATIENTS) {
    cat("!! fewer than ", MIN_PATIENTS,
        " patients in one group, contrast not run\n", sep = "")
    return(NULL)
  }
  
  cells <- cells[md$sampleID %in% keep$sampleID]
  
  # pseudobulk: raw counts summed per patient
  counts <- GetAssayData(seu, assay = "RNA", layer = "counts")[, cells]
  pb <- sapply(keep$sampleID, function(s) {
    idx <- which(seu$sampleID[match(cells, colnames(seu))] == s)
    Matrix::rowSums(counts[, idx, drop = FALSE])
  })
  colnames(pb) <- keep$sampleID
  
  condition <- factor(keep$pathological_response,
                      levels = c("non-MPR", "MPR"))
  coldata <- data.frame(row.names = keep$sampleID, condition = condition)
  
  # ---- gene filtering with filterByExpr ----
  # A DGEList is built solely to run the filter; the retained genes are
  # then passed to DESeq2, which performs the test itself.
  dge  <- DGEList(counts = pb, group = condition)
  keep_genes <- filterByExpr(dge, group = condition)
  pb <- pb[keep_genes, ]
  cat("genes retained (filterByExpr): ", nrow(pb), "\n", sep = "")
  
  dds <- DESeqDataSetFromMatrix(round(pb), coldata, ~ condition)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", "MPR", "non-MPR"))
  
  stat_mat <- as.matrix(data.frame(stat = res$stat,
                                   row.names = rownames(res)))
  stat_mat <- stat_mat[!is.na(stat_mat[, 1]), , drop = FALSE]
  
  acts <- run_ulm(mat = stat_mat, net = net,
                  .source = "source", .target = "target", .mor = "mor",
                  minsize = 5)
  
  acts <- acts %>%
    mutate(group = group_name, histology = HISTOLOGY,
           p_adj = p.adjust(p_value, method = "BH")) %>%
    arrange(p_adj)
  
  cat("\nTFs with adjusted P < 0.05: ", sum(acts$p_adj < 0.05), "\n", sep = "")
  cat("  higher in MPR:     ", sum(acts$p_adj < 0.05 & acts$score > 0), "\n", sep = "")
  cat("  higher in non-MPR: ", sum(acts$p_adj < 0.05 & acts$score < 0), "\n", sep = "")
  
  print(acts %>% filter(p_adj < 0.05) %>% head(15))
  
  acts
}

results_all <- bind_rows(
  run_group("Macrophage_pure"),
  run_group("Monocyte")
)

write_csv(results_all,
          file.path(OUT_TAB, paste0("13_", HISTOLOGY, "_CollecTRI_myeloid.csv")))

cat("\nScript13 (CollecTRI myeloid, ", HISTOLOGY, ") complete.\n", sep = "")