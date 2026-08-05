# ============================================================
# Bloc6_CD8 — Script9b_DEGenes_RegulonCoverage.R
# Purpose: apply to CD8.TEX the two checks performed on the myeloid
# compartment, testing whether the absence of significant factors in
# LUSC follows from a weaker differential signal or from regulon
# coverage rather than from statistical power.
#
# Same gene filter as the CollecTRI script (filterByExpr), so the
# numbers correspond to what run_ulm actually received.
# ============================================================

suppressPackageStartupMessages({
  library(here); library(Seurat); library(DESeq2); library(edgeR)
  library(dplyr); library(tibble); library(readr); library(stringr)
})

MIN_CELLS   <- 10
MIN_TARGETS <- 5

net <- read_tsv(here("Data", "interactions_collectri_9606.tsv.gz"),
                show_col_types = FALSE) %>%
  filter(!str_detect(source, "COMPLEX")) %>%
  filter(!(is_stimulation == 1 & is_inhibition == 1)) %>%
  transmute(source = source_genesymbol, target = target_genesymbol) %>%
  distinct()

check_tex <- function(histology) {
  
  seu <- readRDS(here("Objects", paste0("03_", histology, "_seu_projectils.rds")))
  
  retained <- seu@meta.data %>%
    filter(lineage_CD8_revised == TRUE, functional.cluster == "CD8.TEX") %>%
    dplyr::count(sampleID) %>%
    filter(n >= MIN_CELLS) %>%
    pull(sampleID)
  
  seu_tex <- subset(seu, subset = lineage_CD8_revised == TRUE &
                      functional.cluster == "CD8.TEX" &
                      sampleID %in% retained)
  rm(seu); gc()
  
  pb <- AggregateExpression(seu_tex, group.by = "sampleID", assays = "RNA",
                            slot = "counts", return.seurat = FALSE)$RNA
  
  meta <- seu_tex@meta.data %>% as_tibble() %>%
    distinct(sampleID, pathological_response) %>%
    filter(sampleID %in% colnames(pb)) %>%
    arrange(match(sampleID, colnames(pb))) %>%
    column_to_rownames("sampleID")
  
  condition <- factor(meta$pathological_response,
                      levels = c("non-MPR", "MPR"))
  dge <- DGEList(counts = pb, group = condition)
  pb_filt <- pb[filterByExpr(dge, group = condition), ]
  
  dds <- DESeqDataSetFromMatrix(round(pb_filt), meta, ~ pathological_response)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("pathological_response", "MPR", "non-MPR"))
  res <- res[!is.na(res$padj), ]
  
  sizes <- net %>%
    filter(target %in% rownames(pb_filt)) %>%
    dplyr::count(source) %>%
    filter(n >= MIN_TARGETS)
  
  rm(seu_tex); gc()
  
  tibble(histology = histology,
         n_patients = ncol(pb_filt),
         n_genes_kept = nrow(pb_filt),
         n_genes_tested = nrow(res),
         n_padj_05 = sum(res$padj < 0.05),
         n_padj_01 = sum(res$padj < 0.01),
         median_abs_stat = round(median(abs(res$stat), na.rm = TRUE), 3),
         n_regulons_evaluable = nrow(sizes),
         median_regulon_size = median(sizes$n))
}

res_all <- bind_rows(check_tex("LUAD"), check_tex("LUSC"))

cat("\n=== CD8.TEX — differential signal and regulon coverage ===\n")
print(res_all, width = Inf)