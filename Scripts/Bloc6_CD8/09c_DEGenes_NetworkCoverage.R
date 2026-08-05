# ============================================================
# Bloc6_CD8 — Script09c_DEGenes_NetworkCoverage.R
# Purpose: identify which genes differ between response groups in
# CD8.TEX, and assess how well the CollecTRI network covers them.
#
# Rationale: LUSC yields more differentially expressed genes than LUAD
# yet no significant transcription factor. If the genes driving the LUSC
# signal have few or no documented regulators in the network, the
# discrepancy follows from network coverage rather than from biology.
#
# Input : Objects/03_{HIST}_seu_projectils.rds
#         Data/interactions_collectri_9606.tsv.gz
# ============================================================

suppressPackageStartupMessages({
  library(here); library(Seurat); library(DESeq2); library(edgeR)
  library(dplyr); library(tibble); library(readr); library(stringr)
})

HISTOLOGY <- "LUAD"   # <-- or LUSC
MIN_CELLS <- 10

# ---- 1. Differential expression on CD8.TEX pseudobulk ------------------
seu <- readRDS(here("Objects", paste0("03_", HISTOLOGY, "_seu_projectils.rds")))

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

meta <- seu_tex@meta.data %>%
  as_tibble() %>%
  distinct(sampleID, pathological_response) %>%
  filter(sampleID %in% colnames(pb)) %>%
  arrange(match(sampleID, colnames(pb))) %>%
  column_to_rownames("sampleID")
rm(seu_tex); gc()

condition <- factor(meta$pathological_response, levels = c("non-MPR", "MPR"))
dge <- DGEList(counts = pb, group = condition)
pb  <- pb[filterByExpr(dge, group = condition), ]

dds <- DESeqDataSetFromMatrix(round(pb), meta, ~ pathological_response)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("pathological_response", "MPR", "non-MPR"))

deg <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  filter(!is.na(padj), padj < 0.05) %>%
  arrange(padj)

cat("\n##########", HISTOLOGY, "##########\n")
cat("patients:", ncol(pb), "| genes tested:", nrow(res), "\n")
cat("differential genes (padj < 0.05):", nrow(deg), "\n")
cat("  higher in MPR    :", sum(deg$log2FoldChange > 0), "\n")
cat("  higher in non-MPR:", sum(deg$log2FoldChange < 0), "\n")

cat("\n=== top 30, higher in non-MPR ===\n")
print(as.data.frame(deg %>% filter(log2FoldChange < 0) %>%
                      slice_head(n = 30) %>%
                      select(gene, log2FoldChange, padj)))

cat("\n=== top 30, higher in MPR ===\n")
print(as.data.frame(deg %>% filter(log2FoldChange > 0) %>%
                      slice_head(n = 30) %>%
                      select(gene, log2FoldChange, padj)))

# ---- 2. Network coverage of these genes -------------------------------
net <- read_tsv(here("Data", "interactions_collectri_9606.tsv.gz"),
                show_col_types = FALSE) %>%
  filter(!str_detect(source, "COMPLEX")) %>%
  filter(!(is_stimulation == 1 & is_inhibition == 1)) %>%
  transmute(source = source_genesymbol, target = target_genesymbol) %>%
  distinct()

coverage <- deg %>%
  mutate(n_regulators = sapply(gene, function(g)
    sum(net$target == g))) %>%
  select(gene, log2FoldChange, padj, n_regulators)

cat("\n=== network coverage of the differential genes ===\n")
cat("with no regulator in CollecTRI:",
    sum(coverage$n_regulators == 0), "of", nrow(coverage),
    paste0("(", round(100 * mean(coverage$n_regulators == 0), 1), " %)"), "\n")
cat("median number of regulators:", median(coverage$n_regulators), "\n")

cat("\n=== 20 most significant genes, with their coverage ===\n")
print(as.data.frame(coverage %>% slice_head(n = 20)))