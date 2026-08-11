# ============================================================
# Bloc6_CD8 — Script06_UCell_FunctionalScoring_LUAD.R
# Purpose: score the intensity of three functional programs
# (cytotoxicity, exhaustion, stemness) within the exhausted CD8
# compartment, using UCell with signatures from SignatuR.
#
# Scope: restricted to CD8.TEX and CD8.TPEX within the revised CD8
# lineage, i.e. after the cross-validation correction described in
# Script04. CD8.TPEX is scored here but excluded from the patient-level
# analysis (Script07) for insufficient cell numbers per patient.
#
# Signatures: Tcell.cytotoxicity, Tcell.exhaustion and Tcell.stemness
# from the human collection of SignatuR (Carmona et al., 2020), the
# curated signature database maintained by the same laboratory as
# scGate, UCell and the ProjecTILs reference atlas.
#
# Input : Objects/03_LUAD_seu_projectils.rds
# Output: Objects/06_LUAD_seu_functional_scores.rds
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(SignatuR)
  library(UCell)
  library(dplyr)
})

data(SignatuR)

HISTOLOGY <- "LUAD"
IN_OBJ  <- here("Objects", paste0("03_", HISTOLOGY, "_seu_projectils.rds"))
OUT_OBJ <- here("Objects", paste0("06_", HISTOLOGY, "_seu_functional_scores.rds"))

STATES <- c("CD8.TEX", "CD8.TPEX")

message("Loading ProjecTILs object... ", Sys.time())
seu <- readRDS(IN_OBJ)
message("Cells loaded: ", ncol(seu))

# ---- Restrict to the exhausted compartment ----------------------------
seu_ex <- subset(seu, subset = lineage_CD8_revised == TRUE &
                   functional.cluster %in% STATES)
rm(seu); gc()

message("Cells retained (", paste(STATES, collapse = ", "), "): ", ncol(seu_ex))
print(table(seu_ex$functional.cluster))

# ---- Signatures from the human collection of SignatuR ------------------
sig <- GetSignature(SignatuR$Hs)
sig <- sig[c("Tcell.cytotoxicity", "Tcell.exhaustion", "Tcell.stemness")]

message("Signature sizes:")
print(sapply(sig, length))

# ---- Compute UCell scores ---------------------------------------------
message("Computing UCell scores... ", Sys.time())
seu_ex <- AddModuleScore_UCell(
  seu_ex,
  features = sig,
  ncores   = 4,
  assay    = "RNA",
  name     = ""
)

message("Score distributions:")
print(summary(seu_ex@meta.data[, names(sig)]))

saveRDS(seu_ex, OUT_OBJ)
message("Script06 (functional scoring, ", HISTOLOGY, ") complete. ", Sys.time())