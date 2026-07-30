# ============================================================
# Bloc7_TAM — Script02_AlveolarLabelCheck_LUAD.R
# Purpose: test whether cells labelled "Alveolar macrophages" by
# CellTypist are genuine tissue-resident alveolar macrophages or
# recruited TAMs mislabelled by the reference.
#
# Rationale: the TAM lineage has a different label composition across
# histologies (LUAD: Macrophages only; LUSC: 96% Alveolar macrophages).
# Since alveolar macrophages are developmentally distinct from
# monocyte-derived TAMs, and since the Cheng et al. signatures describe
# recruited TAM states, this asymmetry would invalidate the cross-
# histology comparison if the labels were biologically accurate.
#
# Marker logic:
#   MARCO, FABP4, MCEMP1, PPARG -> tissue-resident alveolar
#   SPP1, C1QC, APOE, TREM2     -> recruited TAM states
#   CD68, CD14                  -> pan-macrophage control
#
# Input : Objects/05_LUAD_seu_annotated.rds
# Output: console only
# ============================================================

library(here); library(Seurat); library(dplyr)

HISTOLOGY <- "LUAD"
seu <- readRDS(here("Objects", "05_LUAD_seu_annotated.rds"))

cat("=== all myeloid-related labels present ===\n")
labs <- table(seu$annotation_final, useNA = "ifany")
print(labs[grepl("acroph|lveol|onocyt", names(labs))])

mk <- c("MARCO","FABP4","MCEMP1","PPARG","CD68","CD14","SPP1","C1QC",
        "APOE","TREM2")
mk <- mk[mk %in% rownames(seu)]

myeloid <- WhichCells(seu, expression = lineage_TAM | lineage_TAM_related)

cat("\n=== mean expression per label (myeloid cells only) ===\n")
avg <- AverageExpression(subset(seu, cells = myeloid),
                         features = mk, group.by = "annotation_final",
                         assays = "RNA", layer = "data")$RNA
print(round(avg, 2))

cat("\nScript02 (alveolar label check,", HISTOLOGY, ") complete.\n")