# ============================================================
# Bloc6_CD8 — Script08_UCell_FunctionalScoring_PatientLevel_Plot_LUSC.R
# Purpose: patient-level boxplot with individual patient points (median
# score per patient) for CD8.TEX and CD8.TPEX functional scores,
# MPR vs non-MPR. Replaces the earlier cell-level violin plot, since
# a violin is misleading with the small number of patient-level
# observations involved (10-38 points), and the annotated p-value must
# match the statistical unit actually shown (patient, not cell).
# CD8.TEX: BH-adjusted patient-level Wilcoxon p-value annotated.
# CD8.TPEX: no p-value shown (descriptive only, see Methods).
# Colors: Okabe & Ito (2008) colorblind-safe palette.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
})

IN_SCORES <- here("Results", "Bloc6_CD8", "Tables", "07_LUSC_patient_level_scores_all.csv")
IN_TEX_STATS <- here("Results", "Bloc6_CD8", "Tables", "07_LUSC_TEX_patientlevel_wilcoxon.csv")
OUT_FIG <- here("Results", "Bloc6_CD8", "Figures")
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

patient_scores <- read_csv(IN_SCORES, show_col_types = FALSE)
tex_stats <- read_csv(IN_TEX_STATS, show_col_types = FALSE)

score_cols <- c("median_Tcell.cytotoxicity", "median_Tcell.exhaustion", "median_Tcell.stemness")
signature_labels <- c(
  median_Tcell.cytotoxicity = "Tcell.cytotoxicity",
  median_Tcell.exhaustion = "Tcell.exhaustion",
  median_Tcell.stemness = "Tcell.stemness"
)

# ---- Long format, only CD8.TEX and CD8.TPEX ----
plot_df <- patient_scores %>%
  filter(functional.cluster %in% c("CD8.TEX", "CD8.TPEX")) %>%
  pivot_longer(cols = all_of(score_cols), names_to = "signature", values_to = "score") %>%
  mutate(signature = recode(signature, !!!signature_labels))

# ---- Annotation: BH-adjusted p-value, CD8.TEX only ----
annot_df <- tex_stats %>%
  mutate(
    functional.cluster = "CD8.TEX",
    p_display = scales::pvalue(p_adj_BH, accuracy = 0.001, add_p = TRUE),
    x = 1.5,
    y = 1.05
  )

# ---- Boxplot with individual patient points ----
p <- ggplot(plot_df, aes(x = pathological_response, y = score, fill = pathological_response)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.08, size = 1.5, alpha = 0.7, color = "black") +
  facet_grid(rows = vars(signature), cols = vars(functional.cluster)) +
  geom_text(data = annot_df, aes(x = x, y = y, label = p_display),
            inherit.aes = FALSE, size = 3.2, fontface = "italic") +
  scale_fill_manual(values = c("MPR" = "#0072B2", "non-MPR" = "#D55E00")) +
  ylim(0, 1.1) +
  labs(title = "LUSC — CD8.TEX / CD8.TPEX patient-level functional scores by pathological response",
       subtitle = "CD8.TPEX shown descriptively (no formal patient-level test; see Methods)",
       x = NULL, y = "Median UCell score per patient") +
  theme_classic() +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 9),
        axis.text.x = element_text(size = 9),
        plot.subtitle = element_text(size = 8, face = "italic"))

ggsave(file.path(OUT_FIG, "08_LUSC_patientlevel_scores_boxplot.png"),
       plot = p, width = 8, height = 9, dpi = 300, bg = "white")

message("Script08 (patient-level boxplot, LUSC) complete.")