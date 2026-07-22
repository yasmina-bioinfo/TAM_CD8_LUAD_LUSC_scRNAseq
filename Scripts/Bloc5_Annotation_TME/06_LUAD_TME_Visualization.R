#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc5_Annotation_TME, Script 06: TME visualization (LUAD)
# UMAP + barplot from annotation_final (CellTypist). Display categories
# with >= 500 cells shown individually, smaller ones grouped into
# "Other (< 500 cells)". Full annotation_final (all labels, including
# Unresolved/Other-contamination) preserved in the saved object for
# downstream Bloc6/7 use via lineage_CD8/lineage_TAM/lineage_TAM_related.
# Statistics on MPR vs non-MPR — pCR already excluded from the object
# since end of Bloc1 QC, no filtering needed at this stage.
# Input:  Objects/05_LUAD_seu_annotated.rds
# ============================================================
suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(scales)
})

IN_OBJ  <- here("Objects", "05_LUAD_seu_annotated.rds")
OUT_FIG <- here("Results", "Bloc5_Annotation_TME", "Figures")
OUT_TAB <- here("Results", "Bloc5_Annotation_TME", "Tables")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_TAB, showWarnings = FALSE, recursive = TRUE)

message("Loading final annotated Seurat object...")
seu <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu))

# -----------------------------
# 1) Simplified display category (visualization only)
# For figure readability, small or unconfirmed categories are grouped:
# non-immune contamination signals are grouped separately from all other
# minor categories (confirmed-but-small, or below the reliable-detection
# threshold for marker verification). The full annotation_final column,
# with every individually confirmed and unconfirmed label, remains
# unchanged and is used for all downstream analysis (see Methods).
# -----------------------------
DISPLAY_MIN_CELLS <- 500

labels_raw <- as.character(seu$annotation_final)

labels_grouped <- labels_raw
labels_grouped[grepl("^Unresolved", labels_raw)] <- "Minor/unresolved categories"
labels_grouped[grepl("^Other/contamination", labels_raw)] <- "Non-immune contamination"
labels_grouped[labels_raw == "Other (mixed, minor)"] <- "Minor/unresolved categories"

counts_display <- table(labels_grouped)
main_categories <- names(counts_display)[counts_display >= DISPLAY_MIN_CELLS]

labels_final <- labels_grouped
labels_final[!(labels_grouped %in% main_categories)] <- "Minor/unresolved categories"

seu$annotation_display <- labels_final

message("Final display categories:")
print(sort(table(seu$annotation_display), decreasing = TRUE))

# -----------------------------
# 2) Color palette — saturation reflects analytical priority, not just
# lineage grouping. Cytotoxic/effector CD8 states and Macrophage/TAM
# (this study's central hypothesis categories) use strong, saturated
# colors; other confirmed populations use moderate colors; any
# "Unresolved"/"Other-contamination" label uses flat greys, deliberately
# de-emphasized since these identities were not marker-confirmed.
# -----------------------------
tme_colors <- c(
  # --- Cells of interest: CD8 cytotoxic/effector (saturated) ---
  "Tem/Trm cytotoxic T cells"     = "#B2182B",
  "Tem/Temra cytotoxic T cells"   = "#67001F",
  
  # --- Cells of interest: Macrophage/TAM (saturated) ---
  "Macrophages"                   = "#D94701",
  
  # --- Other confirmed T/NK/ILC populations (moderate) ---
  "Regulatory T cells"            = "#6BAED6",
  "Tem/Effector helper T cells"   = "#9ECAE1",
  "CRTAM+ gamma-delta T cells"    = "#C994C7",
  "MAIT cells"                    = "#DE9ED6",
  "CD16+ NK cells"                = "#807DBA",
  "ILC3"                          = "#BCBD22",
  
  # --- Confirmed myeloid (non-TAM) and B/plasma (moderate) ---
  "Classical monocytes"           = "#FDBE85",
  "DC2"                           = "#E6AB02",
  "pDC"                           = "#17BECF",
  "Mast cells"                    = "#DDA0DD",
  "Memory B cells"                = "#41AB5D",
  "Plasma cells"                  = "#006D2C",
  "Germinal center B cells"       = "#A1D99B",
  
  # --- Non-immune contamination and minor/unresolved categories, grouped
  # for display readability only (two flat greys) ---
  "Non-immune contamination"    = "#8C8C8C",
  "Minor/unresolved categories" = "#D9D9D9"
)

seu$annotation_display <- factor(seu$annotation_display, levels = names(tme_colors))

# -----------------------------
# 3) UMAP: final annotation
# -----------------------------
message("Generating UMAP...")
png(file.path(OUT_FIG, "06_LUAD_UMAP_final_annotation.png"),
    width = 14, height = 9, units = "in", res = 300)
print(DimPlot(seu, reduction = "umap", group.by = "annotation_display",
              cols = tme_colors, label = FALSE, pt.size = 0.1, raster = FALSE) +
        theme_classic() +
        labs(title = "TME global annotation — GSE243013 LUAD") +
        theme(plot.title = element_text(size = 14, face = "bold"),
              legend.text = element_text(size = 9),
              legend.key.size = unit(0.35, "cm"),
              legend.title = element_blank()))
dev.off()

# -----------------------------
# 4) UMAP split by response (MPR vs non-MPR)
# -----------------------------
png(file.path(OUT_FIG, "06_LUAD_UMAP_split_response.png"),
    width = 16, height = 7, units = "in", res = 300)
print(DimPlot(seu, reduction = "umap", group.by = "annotation_display",
              split.by = "pathological_response", cols = tme_colors,
              label = FALSE, pt.size = 0.1, raster = FALSE) +
        theme_classic() +
        labs(title = "TME annotation — MPR vs non-MPR") +
        theme(plot.title = element_text(size = 13, face = "bold"),
              strip.text = element_text(size = 12, face = "bold"),
              legend.text = element_text(size = 9),
              legend.key.size = unit(0.35, "cm"),
              legend.title = element_blank()))
dev.off()

# -----------------------------
# 5) Statistical tests (MPR vs non-MPR)
# Computed before the barplot so the chi-squared result can be shown
# directly in the figure caption.
# -----------------------------
message("Running statistical tests (MPR vs non-MPR)...")
cont_table <- table(seu$annotation_display, seu$pathological_response)

chisq_res <- chisq.test(cont_table)
message("Chi-2 test p-value: ", chisq_res$p.value)

p_report <- ifelse(chisq_res$p.value == 0, "< 2.2e-16",
                   formatC(chisq_res$p.value, format = "e", digits = 3))

chisq_summary <- data.frame(
  statistic = unname(chisq_res$statistic),
  df        = unname(chisq_res$parameter),
  p.value   = p_report
)
fwrite(chisq_summary, file.path(OUT_TAB, "06_LUAD_chisq_test.csv"))
message("Saved: chi-2 test")

# -----------------------------
# 6) Barplot: proportions, MPR vs non-MPR
# -----------------------------
df_prop <- seu@meta.data %>%
  filter(!is.na(annotation_display), !is.na(pathological_response)) %>%
  count(pathological_response, annotation_display, name = "n") %>%
  group_by(pathological_response) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p_display <- ifelse(chisq_res$p.value < 2.2e-16, "p < 2.2e-16",
                    formatC(chisq_res$p.value, format = "e", digits = 2))

p_bar <- ggplot(df_prop, aes(x = pathological_response, y = prop, fill = annotation_display)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = tme_colors, drop = TRUE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.05)) +
  ylab("Cell-type proportion") +
  labs(caption = paste0("Chi-squared test: X\u00b2 = ", round(chisq_res$statistic, 0),
                        ", df = ", chisq_res$parameter, ", ", p_display)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.35, "cm"),
        legend.title = element_blank(),
        plot.caption = element_text(size = 8, face = "italic", hjust = 0)) +
  guides(fill = guide_legend(ncol = 1))

ggsave(file.path(OUT_FIG, "06_LUAD_Barplot_final_annotation_proportions.png"), p_bar,
       width = 10, height = 8, dpi = 300, bg = "white")
message("Saved: barplot")

# -----------------------------
# 7) Targeted proportion tests, MPR vs non-MPR
# Categories selected before seeing results: B cell (sampling-artifact
# check), CD8 and TAM (central narrative categories). Combined CD8/TAM
# use the lineage_CD8/lineage_TAM flags already built in Script05
# (marker-confirmed labels only), not a name-pattern match — CellTypist
# labels do not share a common prefix like the old "CD8 T (...)" naming.
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
  prop_test_category(seu, "Memory B cells", "annotation_display"),
  prop_test_category(seu, "Tem/Trm cytotoxic T cells", "annotation_display"),
  prop_test_category(seu, "Macrophages", "annotation_display")
)

results_combined <- dplyr::bind_rows(
  prop_test_combined(seu, "lineage_CD8", "CD8 T (all confirmed states combined)"),
  prop_test_combined(seu, "lineage_TAM", "Macrophage/TAM (confirmed only)"),
  prop_test_combined(seu, "lineage_TAM_related", "Macrophage/TAM + classical monocytes")
)

prop_results <- dplyr::bind_rows(results_individual, results_combined)
print(prop_results)

fwrite(prop_results, file.path(OUT_TAB, "06_LUAD_proportion_tests.csv"))
message("Saved: targeted proportion tests")

# -----------------------------
# 8) Save final object
# -----------------------------
saveRDS(seu, here("Objects", "06_LUAD_seu_final.rds"))
message("DONE Script 06 — LUAD TME visualization")