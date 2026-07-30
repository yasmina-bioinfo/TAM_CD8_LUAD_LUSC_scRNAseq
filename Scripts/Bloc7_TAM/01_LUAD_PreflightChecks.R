# ============================================================
# Bloc7_TAM — Script01_PreflightChecks_LUAD.R
# Purpose: verify the composition, per-patient coverage and signature
# gene availability of the TAM compartment before any scoring.
# Run once per histology (see Script01_PreflightChecks_LUSC.R).
#
# Input : Objects/05_LUAD_seu_annotated.rds
# Output: console only — no object or table written
#
# Checks performed:
#   1. Which CellTypist labels carry lineage_TAM / lineage_TAM_related
#   2. TAM cells per patient, and impact of candidate inclusion thresholds
#   3. Presence of Cheng et al. (2021) Table S5 signature genes
#   4. Gene overlap between signatures (non-independence)
#
# Note: Cheng et al. Table S5 uses two non-standard HGNC symbols,
# IL23 and FASL, corrected here to IL23A and FASLG.
# ============================================================

library(here); library(Seurat); library(dplyr); library(tibble)

HISTOLOGY <- "LUAD"
seu <- readRDS(here("Objects", "05_LUAD_seu_annotated.rds"))

# ---- 1. Composition of the TAM lineage flags ---------------------------
cat("=== annotation_final labels flagged as TAM ===\n")
print(table(seu$annotation_final[seu$lineage_TAM], useNA = "ifany"))

cat("\n=== labels flagged as TAM_related (not TAM) ===\n")
print(table(seu$annotation_final[seu$lineage_TAM_related & !seu$lineage_TAM],
            useNA = "ifany"))

# ---- 2. TAM cells per patient ------------------------------------------
md <- seu@meta.data %>%
  rownames_to_column("cellID") %>%
  filter(lineage_TAM) %>%
  select(cellID, sampleID, pathological_response)

cov <- md %>%
  count(sampleID, pathological_response, name = "n_tam") %>%
  as_tibble() %>%
  arrange(n_tam)
print(cov, n = Inf)

cat("\n=== patients retained at candidate thresholds ===\n")
for (thr in c(10, 20, 30, 50, 100)) {
  kept <- cov %>% filter(n_tam >= thr) %>% count(pathological_response)
  cat("threshold", thr, ":", paste(kept$pathological_response, kept$n,
                                   collapse = " | "), "\n")
}

# ---- 3. Signature gene coverage ----------------------------------------
sig <- list(
  M1  = c("IL23A","TNF","CXCL9","CXCL10","CXCL11","CD86","IL1A","IL1B","IL6",
          "CCL5","IRF5","IRF1","CD40","IDO1","KYNU","CCR7"),
  M2  = c("IL4R","CCL4","CCL13","CCL20","CCL17","CCL18","CCL22","CCL24",
          "LYVE1","VEGFA","VEGFB","VEGFC","VEGFD","EGF","CTSA","CTSB","CTSC",
          "CTSD","TGFB1","TGFB2","TGFB3","MMP14","MMP19","MMP9","CLEC7A",
          "WNT7B","FASLG","TNFSF12","TNFSF8","CD276","VTCN1","MSR1","FN1","IRF4"),
  Angiogenesis = c("CCND2","CCNE1","CD44","CXCR4","EDN1","E2F3","EZH2","FGF18",
                   "FGFR1","FYN","HEY1","ITGAV","JAG1","JAG2","MMP9","NOTCH1",
                   "PDGFA","PTK2","SPP1","STC1","TNFAIP6","TYMP","VAV2","VCAN",
                   "VEGFA"),
  Phagocytosis = c("MRC1","CD163","MERTK","C1QB")
)

genes <- rownames(seu)

cat("\n=== signature gene coverage ===\n")
for (s in names(sig)) {
  missing <- setdiff(sig[[s]], genes)
  cat(sprintf("%-13s %2d/%2d present | missing: %s\n",
              s, length(sig[[s]]) - length(missing), length(sig[[s]]),
              if (length(missing)) paste(missing, collapse = ", ") else "none"))
}

# ---- 4. Signature overlap ----------------------------------------------
cat("\n=== genes shared between signatures ===\n")
for (i in 1:(length(sig)-1)) {
  for (j in (i+1):length(sig)) {
    shared <- intersect(sig[[i]], sig[[j]])
    if (length(shared))
      cat(sprintf("%s / %s : %s\n", names(sig)[i], names(sig)[j],
                  paste(shared, collapse = ", ")))
  }
}

cat("\nScript01 (pre-flight checks,", HISTOLOGY, ") complete.\n")