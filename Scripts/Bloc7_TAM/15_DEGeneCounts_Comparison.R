# ============================================================
# Bloc7_TAM — Script15_DEGeneCounts_Comparison.R
# Purpose: compare the number of differentially expressed genes between
# LUAD and LUSC in each myeloid group, to test whether the smaller number
# of significant transcription factors in LUSC follows from a smaller
# differential signal upstream.
#
# Rationale: decoupleR infers activity from the per-gene DESeq2 statistic.
# If fewer genes differ between response groups in LUSC, fewer factors
# will reach significance regardless of patient numbers.
# ============================================================

suppressPackageStartupMessages({
  library(here); library(Seurat); library(dplyr); library(tibble)
  library(DESeq2); library(edgeR)
})

MIN_CELLS    <- 10
MIN_PATIENTS <- 4
GROUPS       <- c("Macrophage_pure", "Transitional", "Monocyte")

count_de <- function(histology) {
  
  seu <- readRDS(here("Objects", paste0("05_", histology, "_seu_annotated.rds")))
  seu_tam <- readRDS(here("Objects", paste0("07_", histology, "_seu_TAM_scored.rds")))
  pure_cells <- colnames(seu_tam)[seu_tam$scGate_Macrophage == "Pure"]
  rm(seu_tam); gc()
  
  seu$myeloid_group <- NA_character_
  seu$myeloid_group[seu$lineage_TAM & colnames(seu) %in% pure_cells]    <- "Macrophage_pure"
  seu$myeloid_group[seu$lineage_TAM & !(colnames(seu) %in% pure_cells)] <- "Transitional"
  seu$myeloid_group[seu$lineage_TAM_related & !seu$lineage_TAM]         <- "Monocyte"
  
  out <- lapply(GROUPS, function(g) {
    
    cells <- colnames(seu)[!is.na(seu$myeloid_group) & seu$myeloid_group == g]
    md <- seu@meta.data[cells, c("sampleID", "pathological_response")]
    
    keep <- md %>% rownames_to_column("cell") %>%
      dplyr::count(sampleID, pathological_response) %>%
      filter(n >= MIN_CELLS)
    
    if (min(table(keep$pathological_response)) < MIN_PATIENTS) return(NULL)
    
    cells <- cells[md$sampleID %in% keep$sampleID]
    counts <- GetAssayData(seu, assay = "RNA", layer = "counts")[, cells]
    
    pb <- sapply(keep$sampleID, function(s) {
      idx <- which(seu$sampleID[match(cells, colnames(seu))] == s)
      Matrix::rowSums(counts[, idx, drop = FALSE])
    })
    colnames(pb) <- keep$sampleID
    
    coldata <- data.frame(row.names = keep$sampleID,
                          condition = factor(keep$pathological_response,
                                             levels = c("non-MPR", "MPR")))
    condition <- factor(keep$pathological_response, levels = c("non-MPR", "MPR"))
    dge <- edgeR::DGEList(counts = pb, group = condition)
    pb <- pb[edgeR::filterByExpr(dge, group = condition), ]
    
    dds <- DESeqDataSetFromMatrix(round(pb), coldata, ~ condition)
    dds <- DESeq(dds, quiet = TRUE)
    res <- results(dds, contrast = c("condition", "MPR", "non-MPR"))
    res <- res[!is.na(res$padj), ]
    
    tibble(histology = histology, group = g,
           n_patients = ncol(pb),
           n_genes_tested = nrow(res),
           n_padj_05  = sum(res$padj < 0.05),
           n_padj_01  = sum(res$padj < 0.01),
           n_p_raw_05 = sum(res$pvalue < 0.05, na.rm = TRUE),
           median_abs_stat = round(median(abs(res$stat), na.rm = TRUE), 3))
  })
  
  rm(seu); gc()
  bind_rows(out)
}

res_all <- bind_rows(count_de("LUAD"), count_de("LUSC"))

cat("\n=== differentially expressed genes per group and histology ===\n")
print(res_all, width = Inf)