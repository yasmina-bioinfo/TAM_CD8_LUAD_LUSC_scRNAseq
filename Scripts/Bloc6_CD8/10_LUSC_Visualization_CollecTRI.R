# ============================================================
# Bloc6_CD8 — Script10_CollecTRI_DotPlot_LUSC.R
# Purpose: dot plot of differential TF activity (CD8.TEX, non-MPR vs
# MPR, pseudobulk + DESeq2 + run_ulm), following the visual convention
# used in the CollecTRI publication's own case study: one dot per TF,
# colored by activity score (direction/magnitude), sized inversely to
# p-value (larger = more significant). No pre-selection of TFs shown.
# Colors: Okabe & Ito (2008) colorblind-safe palette.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
})

IN_RESULTS <- here("Results", "Bloc6_CD8", "Tables",
                   "09_LUSC_TEX_CollecTRI_pseudobulk_DESeq2.csv")
OUT_FIG <- here("Results", "Bloc6_CD8", "Figures")
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

tf_results <- read_csv(IN_RESULTS, show_col_types = FALSE)

# ---- Flag significant TFs (for labeling only, not for filtering the
# plot -- all tested TFs remain shown) ----
plot_df <- tf_results %>%
  mutate(
    significant = p_adj_BH < 0.05,
    neg_log10_p = -log10(p_value)
  )

n_sig <- sum(plot_df$significant)
message("TFs shown: ", nrow(plot_df), " | significant (p_adj < 0.05): ", n_sig)

# ---- Threshold line on the RAW p-value scale that corresponds to
# p_adj_BH = 0.05 -- NOT a fixed -log10(0.05) line, which would mark raw
# significance rather than the BH-corrected significance used everywhere
# else in this project. Found as the largest raw p-value among TFs
# flagged significant after correction. ----
threshold_y <- -log10(max(plot_df$p_value[plot_df$significant], na.rm = TRUE))

# ---- Label top 5 by p-value on each side (non-MPR and MPR), rather
# than top 10 overall -- avoids an asymmetric label set dominated by
# the stronger non-MPR-skewed signal, and keeps the figure legible on
# its own with minimal caption text ----
label_df <- bind_rows(
  plot_df %>% filter(score > 0) %>% arrange(p_value) %>% slice_head(n = 5),
  plot_df %>% filter(score < 0) %>% arrange(p_value) %>% slice_head(n = 5)
)

p <- ggplot(plot_df, aes(x = score, y = neg_log10_p)) +
  geom_point(aes(size = neg_log10_p, color = score), alpha = 0.7) +
  geom_text_repel(data = label_df, aes(label = source),
                  size = 3, max.overlaps = 20, seed = 1) +
  scale_color_distiller(palette = "PRGn", name = "Activity\nscore") +
  scale_size_continuous(range = c(1, 5), guide = "none") +
  geom_hline(yintercept = threshold_y, linetype = "dashed", color = "grey50") +
  labs(title = "LUSC — CD8.TEX differential TF activity, non-MPR vs MPR",
       subtitle = paste0("Dashed line = BH-adjusted p = 0.05 threshold (", n_sig,
                         " significant TFs of ", nrow(plot_df), " tested)"),
       x = "Transcription factor activity (positive = higher in non-MPR, negative = higher in MPR)",
       y = expression(-log[10](p~value))) +
  theme_classic() +
  theme(legend.position = "right",
        plot.subtitle = element_text(size = 9, face = "italic"))

ggsave(file.path(OUT_FIG, "10_LUSC_CollecTRI_dotplot.png"),
       plot = p, width = 8, height = 6, dpi = 300, bg = "white")

message("Script10 (CollecTRI dot plot, LUSC) complete.")