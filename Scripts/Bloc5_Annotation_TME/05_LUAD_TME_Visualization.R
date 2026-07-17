#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc5_Annotation_TME, Script 05: TME visualization (LUAD)
# UMAP + barplot with simplified display categories (20 main categories +
# "Other (cycling, minor)" for rare residual cluster-11 labels). Full
# annotation_consensus (all labels) preserved in the saved object for
# downstream analysis. Statistics restricted to MPR vs non-MPR; pCR
# retained through clustering/annotation for diversity, excluded from
# response-group comparisons onward.
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
main_categories <- c(
  "B cell", "CD4 T (naive/memory)", "CD8 T (GZMK+ effector/memory)",
  "CD8 T (LAG3+ exhausted/TRM)", "Macrophage/TAM (C1Q+/lipid-associated)",
  "Treg", "NK cell", "CD4 T (CXCL13+ exhausted/Tfh-like)",
  "Monocyte (classical, profile A)", "Plasma cell", "Mast cell",
  "Neutrophil", "DC (mature/migratory, mregDC)", "ILC3", "pDC", "GC B cell",
  "Monocyte (classical, profile B)",
  "CD8 T (PDCD1+ exhausted, cycling)", "CD8 T (TCF7+ naive/memory, cycling)",
  "CD8 T (KLRC2+ NK-like terminal, cycling)"
)

seu$annotation_display <- ifelse(
  seu$annotation_consensus %in% main_categories,
  seu$annotation_consensus,
  "Other (cycling, minor)"
)
message("Display categories:")
print(table(seu$annotation_display))

# -----------------------------
# 2) Color palette, organized by lineage
# -----------------------------
tme_colors <- c(
  "B cell"                                     = "#33a02c",
  "GC B cell"                                  = "#74c476",
  "Plasma cell"                                = "#006d2c",
  "CD4 T (naive/memory)"                       = "#a6cee3",
  "CD4 T (CXCL13+ exhausted/Tfh-like)"         = "#6baed6",
  "Treg"                                       = "#08519c",
  "CD8 T (GZMK+ effector/memory)"              = "#f4a582",
  "CD8 T (LAG3+ exhausted/TRM)"                = "#D73027",
  "CD8 T (PDCD1+ exhausted, cycling)"          = "#B2182B",
  "CD8 T (TCF7+ naive/memory, cycling)"        = "#FCAE91",
  "CD8 T (KLRC2+ NK-like terminal, cycling)"   = "#67001F",
  "NK cell"                                    = "#e6550d",
  "ILC3"                                       = "#fdae6b",
  "Macrophage/TAM (C1Q+/lipid-associated)"     = "#a63603",
  "Monocyte (classical, profile A)"            = "#fd8d3c",
  "Monocyte (classical, profile B)"            = "#FFD700",
  "DC (mature/migratory, mregDC)"              = "#BCBD22",
  "pDC"                                        = "#17BECF",
  "Mast cell"                                  = "#c994c7",
  "Neutrophil"                                 = "#969696",
  "Other (cycling, minor)"                     = "#d9d9d9"
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
# 4) UMAP split by response (MPR vs non-MPR only, pCR excluded)
# -----------------------------
seu_resp <- subset(seu, subset = pathological_response %in% c("MPR", "non-MPR"))

png(file.path(OUT_FIG, "05_LUAD_UMAP_split_response.png"),
    width = 14, height = 6, units = "in", res = 300)
print(DimPlot(seu_resp, reduction = "umap", group.by = "annotation_display",
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
message("Saved: UMAP figures")

# -----------------------------
# 5) Barplot: proportions, MPR vs non-MPR only
# -----------------------------
df_prop <- seu_resp@meta.data %>%
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

ggsave(file.path(OUT_FIG, "05_LUAD_Barplot_proportions.png"), p_bar,
       width = 9, height = 7, dpi = 300, bg = "white")
message("Saved: barplot")

# -----------------------------
# 6) Statistical tests, MPR vs non-MPR only
# -----------------------------
message("Running statistical tests (MPR vs non-MPR)...")
cont_table <- table(seu_resp$annotation_display, seu_resp$pathological_response)

chisq_res <- chisq.test(cont_table)
message("Chi-2 test p-value: ", chisq_res$p.value)

chisq_summary <- data.frame(
  statistic = unname(chisq_res$statistic),
  df        = unname(chisq_res$parameter),
  p.value   = unname(chisq_res$p.value)
)
dir.create(OUT_TAB, showWarnings = FALSE, recursive = TRUE)
fwrite(chisq_summary, file.path(OUT_TAB, "05_LUAD_chisq_test.csv"))
message("Saved: chi-2 test")

# -----------------------------
# 7) Save final object
# -----------------------------
saveRDS(seu, here("Objects", "05_LUAD_seu_final.rds"))
message("DONE Script 05 — LUAD TME visualization")