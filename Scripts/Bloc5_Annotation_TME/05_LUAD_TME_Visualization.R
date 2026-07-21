#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc5_Annotation_TME, Script 05: TME visualization (LUAD)
# UMAP + barplot with simplified display categories (categories with
# >= 500 cells shown individually, smaller ones grouped into "Other").
# Full annotation_consensus (all 33 labels) preserved in the saved object
# for downstream Bloc6/7 analysis. Statistics on MPR vs non-MPR — pCR
# already excluded from the object since end of Bloc1 QC, no filtering
# needed at this stage.
# Input:  Objects/04_LUAD_seu_annotated.rds
# ============================================================
suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(scales)
})

IN_OBJ  <- here("Objects", "04_LUAD_seu_annotated.rds")
OUT_FIG <- here("Results", "Bloc5_Annotation_TME", "Figures")
OUT_TAB <- here("Results", "Bloc5_Annotation_TME", "Tables")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)

message("Loading annotated Seurat object...")
seu <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu))

# -----------------------------
# 1) Simplified display category (visualization only; full
#    annotation_consensus untouched for downstream use)
# -----------------------------
# Categories with >= 500 cells shown individually; smaller categories
# grouped into "Other" for figure readability only (annotation_consensus,
# used for downstream Bloc6/7 extraction, is untouched by this threshold)
DISPLAY_MIN_CELLS <- 500

counts_consensus <- table(seu$annotation_consensus)
main_categories  <- names(counts_consensus)[counts_consensus >= DISPLAY_MIN_CELLS]

seu$annotation_display <- ifelse(
  seu$annotation_consensus %in% main_categories,
  seu$annotation_consensus,
  "Other (< 500 cells)"
)
message("Display categories:")
print(sort(table(seu$annotation_display), decreasing = TRUE))

# -----------------------------
# 2) Color palette, organized by lineage
# -----------------------------
tme_colors <- c(
  "B cell"                                       = "#33a02c",
  "GC B cell"                                    = "#74c476",
  "Plasma cell"                                  = "#006d2c",
  "CD4 T (naive/memory)"                         = "#a6cee3",
  "CD4 T (CXCL13+ exhausted/Tfh-like)"           = "#6baed6",
  "Treg"                                         = "#08519c",
  "Treg (cycling)"                               = "#4292c6",
  "CD8 T (GZMK+ effector/memory)"                = "#fdbb84",
  "CD8 T (LAG3+ exhausted/TRM)"                  = "#D73027",
  "CD8 T (cytotoxic effector/memory, cycling)"   = "#7f0000",
  "CD8 T (terminal effector, cluster 6)"         = "#B2182B",
  "CD8 T (effector memory, cluster 6 variant)"   = "#fee08b",
  "NK cell"                                      = "#984ea3",
  "T γδ (Vδ1/non-Vδ2)"                           = "#f781bf",
  "T γδ (Vδ2)"                                   = "#c51b8a",
  "Macrophage/TAM (C1Q+/lipid-associated)"       = "#a63603",
  "Macrophage/TAM (lipid-associated, cycling)"   = "#e6550d",
  "Monocyte classique (S100A8/9+/FCN1+)"         = "#807dba",
  "DC (cDC2, CD1A+/FCER1A+)"                     = "#BCBD22",
  "pDC"                                          = "#17BECF",
  "Mast cell"                                    = "#c994c7",
  "Other (< 500 cells)"                          = "#d9d9d9"
)

seu$annotation_display <- factor(seu$annotation_display, levels = names(tme_colors))

# -----------------------------
# 3) UMAP: final annotation
# -----------------------------
message("Generating UMAP...")
png(file.path(OUT_FIG, "05_LUAD_UMAP_final_annotation.png"),
    width = 12, height = 8, units = "in", res = 300)
print(DimPlot(seu, reduction = "umap", group.by = "annotation_display",
              cols = tme_colors, label = FALSE, pt.size = 0.1, raster = FALSE) +
        theme_classic() +
        labs(title = "TME global annotation — GSE243013 LUAD") +
        theme(plot.title = element_text(size = 14, face = "bold"),
              legend.text = element_text(size = 8),
              legend.key.size = unit(0.4, "cm"),
              legend.title = element_blank()))
dev.off()

# -----------------------------
# 4) UMAP split by response (MPR vs non-MPR — only groups present in
# the object, pCR already excluded upstream at Bloc1 QC)
# -----------------------------
png(file.path(OUT_FIG, "05_LUAD_UMAP_split_response.png"),
    width = 14, height = 6, units = "in", res = 300)
print(DimPlot(seu, reduction = "umap", group.by = "annotation_display",
              split.by = "pathological_response", cols = tme_colors,
              label = FALSE, pt.size = 0.1, raster = FALSE) +
        theme_classic() +
        labs(title = "TME annotation — MPR vs non-MPR") +
        theme(plot.title = element_text(size = 13, face = "bold"),
              strip.text = element_text(size = 12, face = "bold"),
              legend.text = element_text(size = 8),
              legend.key.size = unit(0.4, "cm"),
              legend.title = element_blank()))
dev.off()

# -----------------------------
# 5) Barplot: proportions, MPR vs non-MPR (only groups present in the object)
# -----------------------------
df_prop <- seu@meta.data %>%
  filter(!is.na(annotation_display), !is.na(pathological_response)) %>%
  count(pathological_response, annotation_display, name = "n") %>%
  group_by(pathological_response) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p_bar <- ggplot(df_prop, aes(x = pathological_response, y = prop, fill = annotation_display)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = tme_colors, drop = TRUE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.05)) +
  ylab("Cell-type proportion") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_blank(),
        plot.caption = element_text(size = 8, face = "italic", hjust = 0)) +
  guides(fill = guide_legend(ncol = 1))

# Format p-value for display (avoid showing literal "0")
p_display <- ifelse(chisq_res$p.value < 2.2e-16, "p < 2.2e-16", 
                    formatC(chisq_res$p.value, format = "e", digits = 2))

p_bar <- p_bar +
  labs(caption = paste0("Chi-squared test: X² = ", round(chisq_res$statistic, 0),
                        ", df = ", chisq_res$parameter, ", ", p_display))

ggsave(file.path(OUT_FIG, "05_LUAD_Barplot_proportions.png"), p_bar,
       width = 9, height = 7, dpi = 300, bg = "white")

message("Saved: barplot")

# -----------------------------
# 6) Statistical tests (MPR vs non-MPR — only groups present in the object)
# -----------------------------
message("Running statistical tests (MPR vs non-MPR)...")
cont_table <- table(seu$annotation_display, seu$pathological_response)

chisq_res <- chisq.test(cont_table)
message("Chi-2 test p-value: ", chisq_res$p.value)

# R's floating-point precision floors extremely small p-values to exactly
# 0 on export — not a literal p=0, but below representable precision.
# Report using the standard convention already used in the previous
# pipeline run (p < 2.2e-16, R's machine epsilon).
p_report <- ifelse(chisq_res$p.value == 0, "< 2.2e-16", 
                   formatC(chisq_res$p.value, format = "e", digits = 3))

chisq_summary <- data.frame(
  statistic = unname(chisq_res$statistic),
  df        = unname(chisq_res$parameter),
  p.value   = p_report
)

dir.create(OUT_TAB, showWarnings = FALSE, recursive = TRUE)
fwrite(chisq_summary, file.path(OUT_TAB, "05_LUAD_chisq_test.csv"))
message("Saved: chi-2 test")

# -----------------------------
# 7) Targeted proportion tests, MPR vs non-MPR
# Categories selected before seeing results, to avoid selectively testing
# what looks visually striking after the fact:
# - B cell: sampling-artifact check (MPR has few patients)
# - CD8 and TAM: central narrative categories of this study's hypothesis
# Combined categories (all CD8 / all TAM states pooled) computed from
# annotation_consensus (full granularity), for comparability with the
# single-category narrative framing used in the previous pipeline run.
# -----------------------------
message("Running targeted proportion tests...")

prop_test_category <- function(seu, category_name, group_col, label = category_name) {
  cat_table <- table(seu@meta.data[[group_col]] == category_name, seu$pathological_response)
  res <- prop.test(cat_table["TRUE", ], colSums(cat_table))
  data.frame(
    category    = label,
    prop_MPR    = unname(res$estimate[1]),
    prop_nonMPR = unname(res$estimate[2]),
    p.value     = unname(res$p.value)
  )
}

prop_test_combined <- function(seu, flag_col, label) {
  cat_table <- table(seu@meta.data[[flag_col]], seu$pathological_response)
  res <- prop.test(cat_table["TRUE", ], colSums(cat_table))
  data.frame(
    category    = label,
    prop_MPR    = unname(res$estimate[1]),
    prop_nonMPR = unname(res$estimate[2]),
    p.value     = unname(res$p.value)
  )
}

results_individual <- dplyr::bind_rows(
  prop_test_category(seu, "B cell", "annotation_display"),
  prop_test_category(seu, "CD8 T (GZMK+ effector/memory)", "annotation_display"),
  prop_test_category(seu, "CD8 T (LAG3+ exhausted/TRM)", "annotation_display"),
  prop_test_category(seu, "Macrophage/TAM (C1Q+/lipid-associated)", "annotation_display")
)

seu$is_CD8_any <- grepl("^CD8 T", seu$annotation_consensus)
seu$is_TAM_any <- grepl("^Macrophage/TAM", seu$annotation_consensus)

results_combined <- dplyr::bind_rows(
  prop_test_combined(seu, "is_CD8_any", "CD8 T (all states combined)"),
  prop_test_combined(seu, "is_TAM_any", "Macrophage/TAM (all states combined)")
)

prop_results <- dplyr::bind_rows(results_individual, results_combined)
print(prop_results)

fwrite(prop_results, file.path(OUT_TAB, "05_LUAD_proportion_tests.csv"))
message("Saved: targeted proportion tests")

# -----------------------------
# 8) Save final object
# -----------------------------
saveRDS(seu, here("Objects", "05_LUAD_seu_final.rds"))
message("DONE Script 05 — LUAD TME visualization")