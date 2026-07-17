#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc1_QC, Script 02: QC filtering (LUSC)
# Uniform thresholds across all patients and response groups.
# ============================================================
suppressPackageStartupMessages({
  library(here)
  library(Seurat)
})

IN_OBJ  <- here("Objects", "01_LUSC_seu_raw.rds")
OUT_OBJ <- here("Objects")
OUT_FIG <- here("Results", "Bloc1_QC", "Figures")

message("Loading raw Seurat object...")
seu <- readRDS(IN_OBJ)
message("Cells before filtering: ", ncol(seu))

# Uniform thresholds, based on 99th percentile of the full LUSC distribution
# Never differentiated by pathological_response, to avoid biasing the MPR vs non-MPR comparison.
MIN_FEATURES <- 200
MAX_FEATURES <- 4300
MIN_COUNTS   <- 1000
MAX_COUNTS   <- 26000
MAX_MT       <- 10

message("Distribution before filtering:")
print(table(seu$pathological_response, useNA = "ifany"))

seu_qc <- subset(seu,
                 subset = nFeature_RNA > MIN_FEATURES &
                   nFeature_RNA < MAX_FEATURES &
                   nCount_RNA   > MIN_COUNTS &
                   nCount_RNA   < MAX_COUNTS &
                   percent.mt   < MAX_MT)

message("Cells after filtering: ", ncol(seu_qc))
message("Cells removed: ", ncol(seu) - ncol(seu_qc),
        " (", round(100 * (ncol(seu) - ncol(seu_qc)) / ncol(seu), 1), "%)")

message("Distribution after filtering:")
print(table(seu_qc$pathological_response, useNA = "ifany"))

p_vln <- VlnPlot(seu_qc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                 group.by = "pathological_response", ncol = 3, pt.size = 0) & theme_bw()
ggsave(file.path(OUT_FIG, "02_LUSC_QC_violin_post_filter.png"), p_vln,
       width = 12, height = 4, dpi = 300, bg = "white")

# Convert from BPCells (disk-backed, path-dependent) to standard in-memory
# sparse matrix. Makes the object portable across OS/environments (WSL,
# Windows) from this point onward, avoiding BPCells path-breakage when
# switching environments.
seu_qc[["RNA"]]$counts <- as(seu_qc[["RNA"]]$counts, "dgCMatrix")

saveRDS(seu_qc, file.path(OUT_OBJ, "02_LUSC_seu_qc.rds"))
message("DONE Script 02 — LUSC QC filtering")