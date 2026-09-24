# ============================================================
# Bloc6_CD8 — 01_LUAD_scGate.R  (batched version)
# Purpose: purify CD8 T cells from the full TME object using
# scGate (marker-based, hierarchical gating), independent of any
# prior annotation (CellTypist / annotation_final / lineage_CD8).
# Reference: Andreatta & Carmona, Bioinformatics, 2022.
# Model: scGate_models_DB$human$TME_HiRes
# Batching: whole patients, ~25,000 cells per batch, balanced for
# MPR / non-MPR, to avoid dense-matrix coercion of the full object.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(scGate)
  library(dplyr)
  library(readr)
})

IN_OBJ    <- here("Objects", "06_LUAD_seu_final.rds")
OUT_OBJ   <- here("Objects", "01_LUAD_seu_scGate.rds")
WHOLE_OBJ <- here("Objects", "01_LUAD_seu_scGate_wholeobject.rds")
OUT_TAB   <- here("Results", "Bloc6_CD8", "Tables")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

CELLS_PER_BATCH <- 25000
SEED <- 123

message("Loading LUAD final object... ", Sys.time())
seu <- readRDS(IN_OBJ)
message("Cells loaded: ", ncol(seu))

stopifnot(all(c("annotation_final", "lineage_CD8") %in% colnames(seu@meta.data)))
stopifnot("harmony" %in% names(seu@reductions))

# ---- Cohort guard: 2 chemo-only patients must be absent ----
# Expected: LUAD 50 patients (10 MPR + 40 non-MPR)
print(
  seu@meta.data %>%
    distinct(sampleID, pathological_response) %>%
    count(pathological_response)
)
print(table(seu$pathological_response))
stopifnot(n_distinct(seu$sampleID) == 50)

message("Loading scGate TME_HiRes model database...")
models.TME <- get_scGateDB()

# ---- Traceability: exact model and run information ----
saveRDS(models.TME$human$TME_HiRes,
        here("Objects", "01_LUAD_scGate_TME_HiRes_model.rds"))
writeLines(c(
  "Script : 01_LUAD_scGate.R (batched)",
  paste("Execution date :", Sys.time()),
  paste("scGate :", as.character(packageVersion("scGate"))),
  paste("UCell :", as.character(packageVersion("UCell"))),
  paste("Seurat :", as.character(packageVersion("Seurat"))),
  paste("Cells per batch :", CELLS_PER_BATCH, "| seed :", SEED),
  "Model : TME_HiRes, downloaded via get_scGateDB()"
), file.path(OUT_TAB, "01_LUAD_scGate_run_info.txt"))

# ---- Batch design: whole patients, balanced for response group ----
set.seed(SEED)
n_batches <- ceiling(ncol(seu) / CELLS_PER_BATCH)

pat <- seu@meta.data %>%
  count(sampleID, pathological_response, name = "n_cells") %>%
  mutate(pathological_response = as.character(pathological_response)) %>%
  group_by(pathological_response) %>%
  slice_sample(prop = 1) %>%
  ungroup()

load_mat <- matrix(0, n_batches, 2,
                   dimnames = list(NULL, c("MPR", "non-MPR")))
pat$batch <- NA_integer_
for (i in seq_len(nrow(pat))) {
  g <- pat$pathological_response[i]
  b <- which.min(load_mat[, g] + rowSums(load_mat) * 1e-6)
  pat$batch[i] <- b
  load_mat[b, g] <- load_mat[b, g] + pat$n_cells[i]
}
seu$scGate_batch <- pat$batch[match(seu$sampleID, pat$sampleID)]
stopifnot(!anyNA(seu$scGate_batch))
message("Cells per batch and response group:")
print(table(batch = seu$scGate_batch, seu$pathological_response))

# ---- Lighten object: RNA assay only, KEEP harmony ----
seu_light <- DietSeurat(seu, assays = "RNA", dimreducs = "harmony")

# ---- Run scGate batch by batch ----
res <- vector("list", n_batches)
for (b in seq_len(n_batches)) {
  message("Batch ", b, "/", n_batches, " — ", Sys.time())
  sub <- subset(seu_light, cells = colnames(seu)[seu$scGate_batch == b])
  sub <- scGate(
    sub,
    model     = models.TME$human$TME_HiRes,
    reduction = "harmony",
    ncores    = 1
  )
  res[[b]] <- setNames(as.character(sub$scGate_multi), colnames(sub))
  rm(sub); gc()
}
rm(seu_light); gc()

lab <- unlist(res, use.names = TRUE)
stopifnot(length(lab) == ncol(seu), !anyDuplicated(names(lab)))
seu$scGate_multi <- unname(lab[colnames(seu)])

message("scGate_multi distribution:")
print(table(seu$scGate_multi, useNA = "ifany"))

# ---- Comparison with the whole-object run (if available) ----
if (file.exists(WHOLE_OBJ)) {
  w   <- readRDS(WHOLE_OBJ)
  ref <- setNames(as.character(w$scGate_multi), colnames(w)); rm(w); gc()
  lab_new <- ifelse(is.na(seu$scGate_multi), "Unclassified", seu$scGate_multi)
  lab_ref <- ref[colnames(seu)]
  lab_ref <- ifelse(is.na(lab_ref), "Unclassified", lab_ref)
  message("Same label as whole-object run: ",
          round(100 * mean(lab_new == lab_ref), 2), " %")
  message("CD8T Jaccard (batched vs whole): ",
          round(sum(lab_new == "CD8T" & lab_ref == "CD8T") /
                  sum(lab_new == "CD8T" | lab_ref == "CD8T"), 4))
  cmp <- table(whole = lab_ref, batched = lab_new)
  print(cmp)
  write_csv(as.data.frame.matrix(cmp) %>% tibble::rownames_to_column("whole"),
            file.path(OUT_TAB, "01_LUAD_scGate_batched_vs_whole.csv"))
}

# ---- Cross-check with CellTypist lineage (NA-safe) ----
scGate_CD8T <- seu$scGate_multi %in% "CD8T"
crosscheck <- table(
  scGate_CD8T = scGate_CD8T,
  lineage_CD8 = seu$lineage_CD8,
  useNA = "ifany"
)
print(crosscheck)
stopifnot(sum(crosscheck) == ncol(seu))

write_csv(
  as.data.frame.matrix(crosscheck) %>% tibble::rownames_to_column("scGate_CD8T"),
  file.path(OUT_TAB, "01_LUAD_scGate_vs_lineageCD8_crosscheck.csv")
)

message("scGate CD8T: ", sum(scGate_CD8T))
message("lineage_CD8 (CellTypist): ", sum(seu$lineage_CD8, na.rm = TRUE))
message("Intersection: ", sum(scGate_CD8T & seu$lineage_CD8 %in% TRUE))

saveRDS(seu, OUT_OBJ)
message("Script 01 (scGate, LUAD, batched) complete. ", Sys.time())