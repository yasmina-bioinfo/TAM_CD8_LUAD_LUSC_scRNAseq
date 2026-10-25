# ============================================================
# Bloc6_CD8 — 10_CD8TEX_Visualization_CollecTRI.R
# Purpose: dot plot of differential TF activity in CD8.TEX for both
# histologies, following the visual convention of the CollecTRI
# publication: one dot per TF, coloured by activity score, sized
# inversely to the p-value. No pre-selection of TFs shown.
#
# Contrast is MPR vs non-MPR, so a positive score means higher in
# responding tumors. Color scale uses the same MPR/non-MPR convention
# as the rest of Bloc6 (blue #0072B2 / orange #D55E00), applied here
# as a continuous diverging gradient rather than two discrete groups.
# No title baked in the figure — used in the preprint text.
# ============================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(readr)
  library(ggplot2); library(ggrepel); library(patchwork)
})

HISTOLOGIES <- c("LUAD", "LUSC")
OUT_FIG <- here("Results", "Bloc6_CD8", "Figures")
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

# TFs tested expected from Script 09 (guard against a wrong input file)
EXPECTED_TF <- c(LUAD = 477, LUSC = 459)

make_panel <- function(hist) {
  
  df <- read_csv(here("Results", "Bloc6_CD8", "Tables",
                      paste0("09_", hist, "_TEX_CollecTRI_pseudobulk_DESeq2.csv")),
                 show_col_types = FALSE) %>%
    mutate(significant = p_adj_BH < 0.05,
           neg_log10_p = -log10(p_value))
  
  message(hist, " — TFs tested: ", nrow(df))
  stopifnot(nrow(df) == EXPECTED_TF[[hist]])
  
  n_sig <- sum(df$significant)
  
  # Threshold on the raw p-value scale corresponding to BH-adjusted 0.05,
  # not a fixed -log10(0.05) line. NA when nothing is significant.
  threshold_y <- if (n_sig > 0)
    -log10(max(df$p_value[df$significant], na.rm = TRUE)) else NA
  
  # Top 5 by p-value on each side, so the label set is not dominated by
  # the non-responder-skewed signal
  label_df <- bind_rows(
    df %>% filter(score > 0) %>% arrange(p_value) %>% slice_head(n = 5),
    df %>% filter(score < 0) %>% arrange(p_value) %>% slice_head(n = 5)
  )
  
  subtitle <- if (n_sig > 0)
    paste0(hist, " : ", n_sig, " significant of ", nrow(df), " tested")
  else
    paste0(hist, " : none significant of ", nrow(df), " tested")
  
  ggplot(df, aes(x = score, y = neg_log10_p)) +
    geom_point(aes(size = neg_log10_p, color = score), alpha = 0.7) +
    geom_text_repel(data = label_df, aes(label = source),
                    size = 3.5, max.overlaps = 20, seed = 1) +
    scale_color_gradient2(low = "#D55E00", mid = "grey85", high = "#0072B2",
                          midpoint = 0, name = "Activity\nscore",
                          limits = c(-max(abs(df$score)), max(abs(df$score)))) +
    scale_size_continuous(range = c(1, 5), guide = "none") +
    { if (!is.na(threshold_y))
      geom_hline(yintercept = threshold_y, linetype = "dashed",
                 color = "grey50") } +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
    labs(subtitle = subtitle, x = NULL,
         y = expression(-log[10](p~value))) +
    theme_classic() +
    theme(plot.subtitle = element_text(size = 10, face = "bold"))
}

panels <- lapply(HISTOLOGIES, make_panel)

p_final <- (panels[[1]] | panels[[2]]) +
  plot_layout(guides = "collect") +
  plot_annotation(
    caption = "Dashed line = BH-adjusted p = 0.05. Positive score = higher in responding tumors.",
    theme = theme(plot.caption = element_text(size = 10, face = "italic", hjust = 0)))

ggsave(file.path(OUT_FIG, "10_CD8TEX_CollecTRI_dotplot.png"),
       plot = p_final, width = 12, height = 6, dpi = 300, bg = "white")

message("Script 10 (CollecTRI dot plot, CD8.TEX) complete.")