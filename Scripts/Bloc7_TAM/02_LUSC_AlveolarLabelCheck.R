# ============================================================
# ** NOT RUN — EXCLUDED FROM THE Bloc7_TAM PIPELINE **
#
# Original rationale: the TAM lineage includes a substantial fraction of
# cells labelled "Alveolar macrophages" / "Intermediate macrophages" by
# CellTypist (Bloc5). Alveolar macrophages are developmentally distinct,
# tissue-resident cells, whereas the Cheng et al. (2021) signatures used
# for functional scoring (Script04+) describe recruited TAM states. This
# script tested whether a marker-based check could distinguish resident
# from recruited cells within lineage_TAM.
#
# Reason for exclusion: TAM annotation has no consensus reference atlas
# (unlike ProjecTILs for CD8, Bloc6). Layering an additional, non-
# consensual marker-based test on top of CellTypist + scGate risks
# introducing bias rather than resolving it, with no arbiter to
# adjudicate disagreements between the three. See README / Cahier de
# notes for the full discussion.
#
# This script produces no output consumed by any downstream script
# (Script03/04 onward). Kept for the record only — do not source.
# ============================================================

# ============================================================
# Bloc7_TAM — Script02_AlveolarLabelCheck_LUSC.R

stop("Script02 is excluded from the Bloc7_TAM pipeline — see header above. Do not run.")

library(here); library(Seurat); library(dplyr)

HISTOLOGY <- "LUSC"
seu <- readRDS(here("Objects", "05_LUSC_seu_annotated.rds"))

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