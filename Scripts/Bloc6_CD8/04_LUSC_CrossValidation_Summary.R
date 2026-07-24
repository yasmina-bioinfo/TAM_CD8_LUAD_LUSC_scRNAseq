# ============================================================
# Bloc6_CD8 — Script04_CrossValidation_Summary_LUSC.R
# Purpose: categorize the ProjecTILs vs CellTypist crosscheck (Script03)
# into interpretable groups, and quantify specific cross-validation
# patterns (MAIT confirmation, Treg/NK overlap) before deciding on any
# downstream filtering.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(dplyr)
  library(readr)
})

IN_OBJ  <- here("Objects", "03_LUSC_seu_projectils.rds")
OUT_TAB <- here("Results", "Bloc6_CD8", "Tables")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

message("Loading LUSC ProjecTILs object... ", Sys.time())
seu <- readRDS(IN_OBJ)
message("Cells loaded: ", ncol(seu))

meta <- seu@meta.data

# ---- Categorize annotation_final (CellTypist, Bloc5) into 3 groups ----
# Group 1: CD8 lineage confirmed by CellTypist marker verification (Bloc5)
confirmed_cd8_labels <- c(
  "Tem/Trm cytotoxic T cells",
  "Tem/Temra cytotoxic T cells",
  "Trm cytotoxic T cells",
  "MAIT cells"
)

meta$crosscheck_category <- case_when(
  meta$annotation_final %in% confirmed_cd8_labels ~ "CD8_confirmed_CellTypist",
  grepl("^Unresolved", meta$annotation_final)      ~ "Unresolved_CellTypist",
  grepl("^Other/contamination", meta$annotation_final) ~ "Contamination_flagged_Bloc5",
  TRUE ~ "Non_CD8_other_lineage"
)

# ---- Overall category counts ----
message("Overall crosscheck category counts:")
cat_summary <- table(meta$crosscheck_category)
print(cat_summary)
print(round(100 * prop.table(cat_summary), 1))

write_csv(
  as.data.frame(cat_summary) %>% rename(category = Var1, n_cells = Freq),
  file.path(OUT_TAB, "04_LUSC_crosscheck_category_summary.csv")
)

# ---- Category x functional.cluster breakdown ----
cat_by_cluster <- table(
  crosscheck_category = meta$crosscheck_category,
  functional.cluster   = meta$functional.cluster
)
write.csv(
  as.data.frame.matrix(cat_by_cluster),
  file.path(OUT_TAB, "04_LUSC_category_by_functionalCluster.csv")
)

# ---- Specific pattern 1: MAIT cross-confirmation rate ----
mait_celltypist   <- sum(meta$annotation_final == "MAIT cells")
mait_projectils   <- sum(meta$annotation_final == "MAIT cells" & meta$functional.cluster == "CD8.MAIT")
message("MAIT cells (CellTypist): ", mait_celltypist,
        " | also CD8.MAIT (ProjecTILs): ", mait_projectils,
        " (", round(100 * mait_projectils / mait_celltypist, 1), "%)")

# ---- Specific pattern 2: Regulatory T cells overlap (largest anomaly) ----
treg_total     <- sum(meta$annotation_final == "Regulatory T cells")
treg_by_ptils  <- table(meta$functional.cluster[meta$annotation_final == "Regulatory T cells"])
message("Regulatory T cells (CellTypist): ", treg_total, " — breakdown by ProjecTILs label:")
print(treg_by_ptils)

# ---- Specific pattern 3: NK cells overlap (CD16+/CD16- vs CD8.TEMRA) ----
nk_labels <- c("CD16+ NK cells", "CD16- NK cells")
nk_total  <- sum(meta$annotation_final %in% nk_labels)
nk_temra  <- sum(meta$annotation_final %in% nk_labels & meta$functional.cluster == "CD8.TEMRA")
message("NK cells (CellTypist, CD16+/-): ", nk_total,
        " | of which CD8.TEMRA (ProjecTILs): ", nk_temra,
        " (", round(100 * nk_temra / nk_total, 1), "%)")

message("Script04 (cross-validation summary, LUSC) complete. ", Sys.time())