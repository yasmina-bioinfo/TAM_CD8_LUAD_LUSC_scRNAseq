# ============================================================
# Bloc7_TAM — Script17_CollecTRI_DotPlot_LUAD.R
# Purpose: dot plot of differential TF activity across the three myeloid
# groups, following the same visual convention as the CD8.TEX figure:
# one dot per TF, colored by activity score, sized inversely to p-value.
# No pre-selection of TFs shown.
#
# ⚠ Sign convention differs from Bloc6: here the DESeq2 contrast is
# MPR vs non-MPR, so a positive score means higher in responders. The
# axis label states this explicitly.
# ============================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(readr)
  library(ggplot2); library(ggrepel); library(patchwork)
})

HISTOLOGY <- "LUAD"
GROUPS    <- c("Macrophage_pure", "Transitional", "Monocyte")
LABELS    <- c("Gating-confirmed macrophages", "Transitional cells",
               "Classical monocytes")

IN_RESULTS <- here("Results", "Bloc7_TAM", "Tables",
                   paste0("13_", HISTOLOGY, "_CollecTRI_myeloid.csv"))
OUT_FIG <- here("Results", "Bloc7_TAM", "Figures")
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

acts <- read_csv(IN_RESULTS, show_col_types = FALSE)

make_panel <- function(g, label) {
  
  df <- acts %>%
    filter(group == g) %>%
    mutate(significant = p_adj < 0.05,
           neg_log10_p = -log10(p_value))
  
  n_sig <- sum(df$significant)
  
  # Threshold on the raw p-value scale corresponding to BH-adjusted 0.05,
  # not a fixed -log10(0.05) line
  threshold_y <- if (n_sig > 0)
    -log10(max(df$p_value[df$significant], na.rm = TRUE)) else NA
  
  # Top 5 by p-value on each side, to avoid a label set dominated by the
  # non-responder-skewed signal
  label_df <- bind_rows(
    df %>% filter(score > 0) %>% arrange(p_value) %>% slice_head(n = 5),
    df %>% filter(score < 0) %>% arrange(p_value) %>% slice_head(n = 5)
  )
  
  ggplot(df, aes(x = score, y = neg_log10_p)) +
    geom_point(aes(size = neg_log10_p, color = score), alpha = 0.7) +
    geom_text_repel(data = label_df, aes(label = source),
                    size = 3.5, max.overlaps = 20, seed = 1) +
    scale_color_distiller(palette = "PRGn", name = "Activity\nscore",
                          limits = c(-max(abs(df$score)), max(abs(df$score)))) +
    scale_size_continuous(range = c(1, 5), guide = "none") +
    { if (!is.na(threshold_y))
      geom_hline(yintercept = threshold_y, linetype = "dashed",
                 color = "grey50") } +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
    labs(subtitle = paste0(label, " — ", n_sig, " significant of ",
                           nrow(df), " tested"),
         x = NULL, y = expression(-log[10](p~value))) +
    theme_classic() +
    theme(plot.subtitle = element_text(size = 9, face = "bold"),
          legend.position = "none")
}

panels <- Map(make_panel, GROUPS, LABELS)

p_final <- (panels[[1]] | panels[[2]] | panels[[3]]) +
  plot_annotation(
    title = paste0(HISTOLOGY, " — myeloid differential TF activity, MPR vs non-MPR"),
    caption = "Dashed line = BH-adjusted p = 0.05. Positive score = higher in responding tumors.",
    theme = theme(plot.title = element_text(size = 12, face = "bold"),
                  plot.caption = element_text(size = 8, face = "italic", hjust = 0)))

ggsave(file.path(OUT_FIG, paste0("17_", HISTOLOGY, "_CollecTRI_dotplot.png")),
       plot = p_final, width = 18, height = 6.5, dpi = 300, bg = "white")

message("Script17 (CollecTRI dot plot, ", HISTOLOGY, ") complete.")