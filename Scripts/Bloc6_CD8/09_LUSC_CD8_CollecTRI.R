# ============================================================
# Bloc6_CD8 — Script09_CollecTRI_Pseudobulk_LUSC.R
# Purpose: identify TF whose activity differs between MPR and non-MPR
# in CD8.TEX, using a pseudobulk + DESeq2 differential signature (one
# statistic per gene, patient-level comparison already resolved) as
# input to a single run_ulm call
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(DESeq2)
  library(decoupleR)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

IN_OBJ  <- here("Objects", "03_LUSC_seu_projectils.rds")
IN_TCRI <- here("Data", "interactions_collectri_9606.tsv.gz")
OUT_TAB <- here("Results", "Bloc6_CD8", "Tables")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

# ---- Step 1: load CD8 object, restrict to CD8.TEX cells of patients
# with >=10 TEX cells (identical threshold and patient set as Script07,
# for consistency between the two functional analyses) ----
message("Loading LUSC CD8 object... ", Sys.time())
seu <- readRDS(IN_OBJ)

retained_patients <- seu@meta.data %>%
  filter(lineage_CD8_revised == TRUE, functional.cluster == "CD8.TEX") %>%
  count(sampleID) %>%
  filter(n >= 10) %>%
  pull(sampleID)

seu_tex <- subset(seu, subset = lineage_CD8_revised == TRUE &
                    functional.cluster == "CD8.TEX" &
                    sampleID %in% retained_patients)
message("Patients retained: ", length(retained_patients),
        " | CD8.TEX cells retained: ", ncol(seu_tex))

# ---- Step 2: pseudobulk aggregation -- sum RAW counts per patient ----
# AggregateExpression() sums raw counts (not normalized data) across all
# cells belonging to the same sampleID; this is exactly what DESeq2's
# statistical model requires as input (see CAHIER_METHODOLOGIQUE for the
# citation on why raw, un-normalized counts are mandatory here).
message("Aggregating raw counts per patient (pseudobulk)...")
pseudobulk <- AggregateExpression(
  seu_tex, group.by = "sampleID", assays = "RNA", slot = "counts",
  return.seurat = FALSE
)$RNA
message("Pseudobulk matrix: ", nrow(pseudobulk), " genes x ", ncol(pseudobulk), " patients")

# ---- Step 3: DESeq2 differential analysis, MPR vs non-MPR ----
# Build per-patient metadata (one row per pseudobulk column)
patient_meta <- seu_tex@meta.data %>%
  as_tibble() %>%
  distinct(sampleID, pathological_response) %>%
  filter(sampleID %in% colnames(pseudobulk)) %>%
  arrange(match(sampleID, colnames(pseudobulk))) %>%
  column_to_rownames("sampleID")

# Standard practice: remove genes with very low counts across all
# patients before fitting -- improves dispersion estimation and reduces
# the multiple-testing burden without affecting the biology of interest.
keep_genes <- rowSums(pseudobulk) >= 10
pseudobulk_filt <- pseudobulk[keep_genes, ]
message("Genes retained after low-count filtering: ", nrow(pseudobulk_filt),
        " / ", nrow(pseudobulk))

dds <- DESeqDataSetFromMatrix(
  countData = round(pseudobulk_filt),
  colData   = patient_meta,
  design    = ~ pathological_response
)
dds <- DESeq(dds)

# Wald statistic per gene: direction + significance in one number,
# exactly the kind of input the decoupleR authors recommend for run_ulm
# (equivalent to limma's t-value in their bulk RNA-seq vignette).
res <- results(dds, contrast = c("pathological_response", "non-MPR", "MPR"))
message("DESeq2 differential analysis done. ", Sys.time())

deg_stat <- as.data.frame(res) %>%
  filter(!is.na(stat)) %>%
  tibble::rownames_to_column("gene") %>%
  select(gene, stat)

message("Genes with a valid differential statistic: ", nrow(deg_stat))

# ---- Step 4: CollecTRI network (identical to Script09, reused as-is) ----
collectri_raw <- read_tsv(IN_TCRI, show_col_types = FALSE)

net <- collectri_raw %>%
  filter(!str_detect(source, "COMPLEX")) %>%
  filter(!(is_stimulation == 1 & is_inhibition == 1)) %>%
  transmute(
    source = source_genesymbol,
    target = target_genesymbol,
    mor    = if_else(is_inhibition == 1, -1, 1)
  ) %>%
  distinct()

# ---- Step 5: run_ulm ONCE on the differential statistic vector ----
# mat here is a single-column matrix (genes x 1 "sample", which is really
# the MPR-vs-non-MPR contrast already resolved by DESeq2) -- not a
# cell-by-gene matrix as in the earlier per-cell approach.
stat_mat <- deg_stat %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()
colnames(stat_mat) <- "nonMPR_vs_MPR"

tf_results <- run_ulm(mat = stat_mat, net = net,
                      .source = 'source', .target = 'target', .mor = 'mor',
                      minsize = 5) %>%
  arrange(p_value) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

message("TFs tested (passed minsize filter): ", nrow(tf_results))
message("Top 10 TF by p-value:")
print(head(tf_results, 10))

write_csv(tf_results, file.path(OUT_TAB, "09_LUSC_TEX_CollecTRI_pseudobulk_DESeq2.csv"))

message("Script09b (CollecTRI pseudobulk, LUSC) complete. ", Sys.time())