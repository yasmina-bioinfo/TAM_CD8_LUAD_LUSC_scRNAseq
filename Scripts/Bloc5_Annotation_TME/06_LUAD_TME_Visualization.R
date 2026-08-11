#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc5_Annotation_TME, Script 06: TME visualization (LUAD)
#
# REVISED 01/08/2026: cell-level tests replaced by patient-level tests.
# The chi-squared test and prop.test treated cells from the same patient
# as independent observations, which inflates significance — the same
# pseudoreplication corrected elsewhere in this project. The barplot is
# retained as a descriptive figure, without an accompanying test, and
# boxplots of per-patient proportions are added for the two compartments
# characterised downstream.
#
# Input:  Objects/05_LUAD_seu_annotated.rds
# ============================================================
suppressPackageStartupMessages({
  library(here); library(Seurat); library(ggplot2)
  library(dplyr); library(data.table); library(scales)
  library(tibble); library(patchwork)
})

IN_OBJ  <- here("Objects", "05_LUAD_seu_annotated.rds")
OUT_FIG <- here("Results", "Bloc5_Annotation_TME", "Figures")
OUT_TAB <- here("Results", "Bloc5_Annotation_TME", "Tables")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_TAB, showWarnings = FALSE, recursive = TRUE)

HISTOLOGY <- "LUAD"
MIN_CELLS <- 50   # total annotated cells per patient. A proportion
                  # computed over a very small denominator is unstable
                  # regardless of the population's abundance: one cell
                  # out of twelve already represents eight percent.

message("Loading final annotated Seurat object...")
seu <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu))

# -----------------------------
# 1) Simplified display category (visualization only) — unchanged
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
# 2) Color palette — unchanged
# -----------------------------
tme_colors <- c(
  "Tem/Trm cytotoxic T cells"     = "#B2182B",
  "Tem/Temra cytotoxic T cells"   = "#67001F",
  "Macrophages"                   = "#D94701",
  "Regulatory T cells"            = "#6BAED6",
  "Tem/Effector helper T cells"   = "#9ECAE1",
  "CRTAM+ gamma-delta T cells"    = "#C994C7",
  "MAIT cells"                    = "#DE9ED6",
  "CD16+ NK cells"                = "#807DBA",
  "ILC3"                          = "#BCBD22",
  "Classical monocytes"           = "#FDBE85",
  "DC2"                           = "#E6AB02",
  "pDC"                           = "#17BECF",
  "Mast cells"                    = "#DDA0DD",
  "Memory B cells"                = "#41AB5D",
  "Plasma cells"                  = "#006D2C",
  "Germinal center B cells"       = "#A1D99B",
  "Non-immune contamination"      = "#8C8C8C",
  "Minor/unresolved categories"   = "#D9D9D9"
)

seu$annotation_display <- factor(seu$annotation_display, levels = names(tme_colors))

# -----------------------------
# 3-4) UMAPs — unchanged
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
# 5) Panel A — descriptive barplot, NO test
# Cumulative proportions across all cells of each group. Shown for
# overview only: cells from one patient are not independent, so no
# statistic is attached to this display.
# -----------------------------
df_prop <- seu@meta.data %>%
  filter(!is.na(annotation_display), !is.na(pathological_response)) %>%
  count(pathological_response, annotation_display, name = "n") %>%
  group_by(pathological_response) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p_bar <- ggplot(df_prop, aes(x = pathological_response, y = prop,
                             fill = annotation_display)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = tme_colors, drop = TRUE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.05)) +
  ylab("Cell-type proportion") +
  labs(subtitle = "A. Cumulative composition (descriptive)") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_blank()) +
  guides(fill = guide_legend(ncol = 1))

# -----------------------------
# 6) Patient-level proportions and tests
# The two compartments characterised downstream (Bloc6, Bloc7).
# -----------------------------
message("Computing patient-level proportions...")

cd8_col <- if ("lineage_CD8_revised" %in% colnames(seu@meta.data))
  "lineage_CD8_revised" else "lineage_CD8"

res <- seu@meta.data %>%
  as_tibble() %>%
  group_by(sampleID, pathological_response) %>%
  summarise(n_total = n(),
            pct_CD8 = 100 * mean(.data[[cd8_col]]),
            pct_TAM = 100 * mean(lineage_TAM),
            .groups = "drop") %>%
  filter(n_total >= MIN_CELLS)

cat("\n=== patients retained (>=", MIN_CELLS, "annotated cells) ===\n")
print(table(res$pathological_response))

tests <- tibble(
  population = c("CD8 T cells", "Macrophages"),
  median_MPR = c(median(res$pct_CD8[res$pathological_response == "MPR"]),
                 median(res$pct_TAM[res$pathological_response == "MPR"])),
  median_nonMPR = c(median(res$pct_CD8[res$pathological_response == "non-MPR"]),
                    median(res$pct_TAM[res$pathological_response == "non-MPR"])),
  p_value = c(wilcox.test(pct_CD8 ~ pathological_response, data = res)$p.value,
              wilcox.test(pct_TAM ~ pathological_response, data = res)$p.value)
) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

cat("\n=== TESTS (Wilcoxon patient-level, BH across the two populations) ===\n")
print(tests)

fwrite(res,   file.path(OUT_TAB, "06_LUAD_TME_patient_level.csv"))
fwrite(tests, file.path(OUT_TAB, "06_LUAD_TME_patient_level_tests.csv"))

# -----------------------------
# 7) Panels B and C — per-patient boxplots with points
# -----------------------------
make_box <- function(varname, popname, colour, panel_letter, p_adj_value) {
  ggplot(res, aes(x = pathological_response, y = .data[[varname]])) +
    geom_boxplot(outlier.shape = NA, width = 0.55, fill = colour, alpha = 0.35) +
    geom_jitter(width = 0.15, size = 1.6, alpha = 0.7, colour = colour) +
    ylab(paste0(popname, " (% of annotated cells)")) +
    labs(subtitle = paste0(panel_letter, ". ", popname, " per patient"),
         caption = paste0("Wilcoxon, BH-adjusted P = ",
                          formatC(p_adj_value, format = "f", digits = 2))) +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 12, face = "bold"),
          plot.caption = element_text(size = 8, face = "italic", hjust = 0))
}

p_cd8 <- make_box("pct_CD8", "CD8 T cells", "#B2182B", "B", tests$p_adj[1])
p_tam <- make_box("pct_TAM", "Macrophages", "#D94701", "C", tests$p_adj[2])

p_final <- (p_bar | (p_cd8 / p_tam)) +
  plot_layout(widths = c(1, 1), guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(OUT_FIG, "06_LUAD_TME_composition.png"), p_final,
       width = 14, height = 9, dpi = 300, bg = "white")
message("Saved: composition figure (3 panels)")

# -----------------------------
# 8) Save final object
# -----------------------------
saveRDS(seu, here("Objects", "06_LUAD_seu_final.rds"))
message("DONE Script 06 — LUAD TME visualization")