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
main_categories <- names(counts_display)[
  counts_display >= DISPLAY_MIN_CELLS | names(counts_display) == "Non-immune contamination"
]

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
  "Regulatory T cells"            = "#6BAED6",
  "Tem/Effector helper T cells"   = "#9ECAE1",
  "Tcm/Naive helper T cells"      = "#4292C6",
  "CRTAM+ gamma-delta T cells"    = "#C994C7",
  "MAIT cells"                    = "#DE9ED6",
  "CD16+ NK cells"                = "#807DBA",
  "Macrophages"                   = "#D94701",
  "Classical monocytes"           = "#FDBE85",
  "DC2"                           = "#E6AB02",
  "DC1"                           = "#B15928",
  "Migratory DCs"                 = "#FDD0A2",
  "Age-associated B cells"        = "#238B45",
  "Plasma cells"                  = "#006D2C",
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
# 5) Panel A — descriptive barplots, NO test
# Two versions saved: (1) full composition including Unresolved/
# contamination bins, kept for internal reference and traceability;
# (2) composition among marker-confirmed cells only, used in the
# manuscript for readability — the Confirmed/Unresolved/Contamination
# percentages are reported as text alongside this second figure, not
# encoded in it.
# Cells from one patient are not independent, so no statistic is
# attached to either display.
# -----------------------------

# --- Version 1: full composition (all cells, all bins) ---
df_prop_full <- seu@meta.data %>%
  filter(!is.na(annotation_display), !is.na(pathological_response)) %>%
  count(pathological_response, annotation_display, name = "n") %>%
  group_by(pathological_response) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p_bar_full <- ggplot(df_prop_full, aes(x = pathological_response, y = prop,
                                       fill = annotation_display)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = tme_colors, drop = TRUE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.05)) +
  ylab("Cell-type proportion") +
  labs(subtitle = "A. Cumulative composition, all cells (descriptive)") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_blank()) +
  guides(fill = guide_legend(ncol = 1))

# --- Version 2: marker-confirmed cells only (manuscript figure) ---
df_prop_confirmed <- seu@meta.data %>%
  filter(!is.na(annotation_display), !is.na(pathological_response),
         !(annotation_display %in% c("Non-immune contamination", "Minor/unresolved categories"))) %>%
  count(pathological_response, annotation_display, name = "n") %>%
  group_by(pathological_response) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p_bar_confirmed <- ggplot(df_prop_confirmed, aes(x = pathological_response, y = prop,
                                                 fill = annotation_display)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = tme_colors, drop = TRUE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.05)) +
  ylab("Cell-type proportion (marker-confirmed cells only)") +
  labs(subtitle = "A. Composition among marker-confirmed cells") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_blank()) +
  guides(fill = guide_legend(ncol = 1))

ggsave(file.path(OUT_FIG, paste0("06_", HISTOLOGY, "_TME_composition_barplot_full.png")), p_bar_full,
       width = 8, height = 9, dpi = 300, bg = "white")
ggsave(file.path(OUT_FIG, paste0("06_", HISTOLOGY, "_TME_composition_barplot_confirmed.png")), p_bar_confirmed,
       width = 8, height = 9, dpi = 300, bg = "white")
message("Saved: both composition barplots (full and confirmed-only)")

# -----------------------------
# 6) Patient-level proportions and tests — every category shown in the
# barplot (annotation_display), tested identically. CD8 and TAM are not
# singled out here: each is computed and tested in its own dedicated
# block (Bloc6, Bloc7), on its own more rigorously defined population.
# -----------------------------
message("Computing patient-level proportions for all displayed TME categories...")

test_groups <- setdiff(levels(seu$annotation_display),
                       c("Non-immune contamination", "Minor/unresolved categories"))

patient_meta <- seu@meta.data %>%
  as_tibble() %>%
  group_by(sampleID, pathological_response) %>%
  summarise(n_total = n(), .groups = "drop") %>%
  filter(n_total >= MIN_CELLS)

res_all <- purrr::map_dfr(test_groups, function(grp) {
  seu@meta.data %>%
    as_tibble() %>%
    filter(sampleID %in% patient_meta$sampleID) %>%
    group_by(sampleID, pathological_response) %>%
    summarise(pct = 100 * mean(annotation_display == grp), .groups = "drop") %>%
    mutate(population = grp)
})

tests_all <- res_all %>%
  group_by(population) %>%
  summarise(
    median_MPR = median(pct[pathological_response == "MPR"]),
    median_nonMPR = median(pct[pathological_response == "non-MPR"]),
    p_value = tryCatch(
      wilcox.test(pct ~ pathological_response)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

cat("\n=== TESTS (Wilcoxon patient-level, BH across", length(test_groups), "displayed TME categories) ===\n")
print(tests_all, n = Inf)

fwrite(res_all,   file.path(OUT_TAB, paste0("06_", HISTOLOGY, "_TME_patient_level_all.csv")))
fwrite(tests_all, file.path(OUT_TAB, paste0("06_", HISTOLOGY, "_TME_patient_level_tests_all.csv")))

# -----------------------------
# 7) Save
# -----------------------------

saveRDS(seu, here("Objects", paste0("06_", HISTOLOGY, "_seu_final.rds")))
message("DONE Script 06 — ", HISTOLOGY, " TME visualization")