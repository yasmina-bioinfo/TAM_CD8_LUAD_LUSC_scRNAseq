# ============================================================
# Bloc6_CD8 — 06_LUAD_UCell_FunctionalScoring.R
# Purpose: score the intensity of three functional programs
# (cytotoxicity, exhaustion, stemness) within the exhausted CD8
# compartment (CD8.TEX), using UCell with signatures from SignatuR.
#
# Scope: restricted to CD8.TEX cells (purified population from
# scGate + UCell cycling/IFN-high removal + ProjecTILs, Scripts 01-03).
# CD8.TPEX is not scored here: it is excluded from all patient-level
# functional analyses for insufficient cell numbers per patient (see
# 03b_LUAD_state_testability — 86 % of patients have >=1 TPEX cell,
# but 39.5 % of those have fewer than 10).
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

STATE <- "CD8.TEX"
EXPECTED_TEX_CELLS <- 14865 

message("Loading ProjecTILs object... ", Sys.time())
seu <- readRDS(IN_OBJ)
message("Cells loaded: ", ncol(seu))
stopifnot("functional.cluster" %in% colnames(seu@meta.data))

# ---- Restrict to the exhausted compartment ----------------------------
seu_ex <- subset(seu, subset = functional.cluster == STATE)
rm(seu); gc()

message("Cells retained (", STATE, "): ", ncol(seu_ex))
stopifnot(ncol(seu_ex) == EXPECTED_TEX_CELLS)
print(table(seu_ex$pathological_response))

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
message("Script 06 (functional scoring, ", HISTOLOGY, ") complete. ", Sys.time())