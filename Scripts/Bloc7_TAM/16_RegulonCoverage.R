# ============================================================
# Bloc7_TAM — Script16_RegulonCoverage.R
# Purpose: test whether the smaller number of significant transcription
# factors in LUSC follows from fewer regulons being evaluable, rather
# than from a weaker biological signal.
#
# Rationale: decoupleR requires a minimum of five targets present in the
# tested gene space. filterByExpr is more demanding as cohort size grows,
# since the number of required samples follows the smallest group. Fewer
# genes tested may therefore mean fewer evaluable regulons.
# ============================================================

suppressPackageStartupMessages({
  library(here); library(Seurat); library(dplyr); library(tibble)
  library(readr); library(stringr); library(edgeR)
})

MIN_CELLS   <- 10
MIN_TARGETS <- 5
GROUPS      <- c("Macrophage_pure", "Monocyte")

# ---- CollecTRI network, filtered as in Script13 ------------------------
net <- read_tsv(here("Data", "interactions_collectri_9606.tsv.gz"),
                show_col_types = FALSE) %>%
  filter(!str_detect(source, "COMPLEX")) %>%
  filter(!(is_stimulation == 1 & is_inhibition == 1)) %>%
  transmute(source = source_genesymbol, target = target_genesymbol) %>%
  distinct()

cat("CollecTRI: ", n_distinct(net$source), " factors\n", sep = "")

# ---- Genes retained per group and histology ---------------------------
coverage <- function(histology) {
  
  seu <- readRDS(here("Objects", paste0("04b_", histology, "_seu_TAM_classified.rds")))
  
  seu$myeloid_group <- NA_character_
  seu$myeloid_group[seu$scGate_TAM_class == "Macrophage_pure"] <- "Macrophage_pure"
  seu$myeloid_group[seu$scGate_TAM_class == "Monocyte_pure"]   <- "Monocyte"      
  
  out <- lapply(GROUPS, function(g) {
    
    cells <- colnames(seu)[!is.na(seu$myeloid_group) & seu$myeloid_group == g]
    md <- seu@meta.data[cells, c("sampleID", "pathological_response")]
    
    # patient-level condition table, kept for the filterByExpr design
    keep <- md %>%
      rownames_to_column("cell") %>%
      dplyr::count(sampleID, pathological_response) %>%
      filter(n >= MIN_CELLS)
    
    cells  <- cells[md$sampleID %in% keep$sampleID]
    counts <- GetAssayData(seu, assay = "RNA", layer = "counts")[, cells]
    
    pb <- sapply(keep$sampleID, function(s) {
      idx <- which(seu$sampleID[match(cells, colnames(seu))] == s)
      Matrix::rowSums(counts[, idx, drop = FALSE])
    })
    colnames(pb) <- keep$sampleID
    
    cond  <- factor(keep$pathological_response, levels = c("non-MPR", "MPR"))
    dge   <- DGEList(counts = pb, group = cond)
    genes <- rownames(pb)[filterByExpr(dge, group = cond)]
    
    # how many regulons keep at least MIN_TARGETS in this gene space
    sizes <- net %>%
      filter(target %in% genes) %>%
      dplyr::count(source) %>%
      filter(n >= MIN_TARGETS)
    
    tibble(histology = histology, group = g,
           n_patients = ncol(pb),
           n_genes = length(genes),
           n_regulons_evaluable = nrow(sizes),
           median_regulon_size = median(sizes$n))
  })
  
  rm(seu); gc()
  bind_rows(out)
}

res <- bind_rows(coverage("LUAD"), coverage("LUSC"))

cat("\n=== regulon coverage per group and histology ===\n")
print(res, width = Inf)