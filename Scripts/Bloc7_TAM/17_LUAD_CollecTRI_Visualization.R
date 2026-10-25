# ============================================================
# Bloc7_TAM — Script17_CollecTRI_DotPlot_LUAD.R
# Purpose: dot plot of differential TF activity in the two myeloid
# populations (Macrophage_pure, Monocyte), following the same visual
# convention as the CD8.TEX figure: one dot per TF, colored by activity
# score, sized by -log10(p). No pre-selection of TFs shown.
#
# Sign convention: the DESeq2 contrast is MPR vs non-MPR, so a positive
# score means higher inferred activity in responders (MPR, blue) and a
# negative score higher in non-responders (non-MPR, orange).
# Colors: Okabe-Ito, as for all MPR vs non-MPR figures of the project.
# ============================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(readr)
  library(ggplot2); library(ggrepel); library(patchwork)
})

HISTOLOGY <- "LUAD"
GROUPS <- c("Macrophage_pure", "Monocyte")   
LABELS <- c("Macrophage_pure", "Monocyte_pure")

COL_MPR    <- "#0072B2"
COL_NONMPR <- "#D55E00"

IN_RESULTS <- here("Results", "Bloc7_TAM", "Tables",
                   paste0("13_", HISTOLOGY, "_CollecTRI_myeloid.csv"))
OUT_FIG <- here("Results", "Bloc7_TAM", "Figures")
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

acts <- read_csv(IN_RESULTS, show_col_types = FALSE)

# Input guard: right file, right histology, expected groups only
stopifnot(all(acts$histology == HISTOLOGY),
          all(GROUPS %in% acts$group),
          !"Transitional" %in% acts$group)

make_panel <- function(g, label) {
  
  df <- acts %>%
    filter(group == g) %>%
    mutate(significant = p_adj < 0.05,
           neg_log10_p = -log10(p_value))
  
  n_sig <- sum(df$significant)
  
  # Threshold on the raw p-value scale corresponding to BH-adjusted 0.05
  threshold_y <- if (n_sig > 0)
    -log10(max(df$p_value[df$significant], na.rm = TRUE)) else NA
  
  # Top 5 by p-value on each side
  label_df <- bind_rows(
    df %>% filter(score > 0) %>% arrange(p_value) %>% slice_head(n = 5),
    df %>% filter(score < 0) %>% arrange(p_value) %>% slice_head(n = 5)
  )
  
  lim <- max(abs(df$score))
  
  ggplot(df, aes(x = score, y = neg_log10_p)) +
    geom_point(aes(size = neg_log10_p, color = score), alpha = 0.7) +
    geom_text_repel(data = label_df, aes(label = source),
                    size = 3.5, max.overlaps = 20, seed = 1) +
    scale_color_gradient2(low = COL_NONMPR, mid = "grey85", high = COL_MPR,
                          midpoint = 0, limits = c(-lim, lim),
                          name = "Activity\nscore") +
    scale_size_continuous(range = c(1, 5), guide = "none") +
    { if (!is.na(threshold_y))
      geom_hline(yintercept = threshold_y, linetype = "dashed",
                 color = "grey50") } +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
    labs(subtitle = paste0(label, ": ", n_sig, " significant of ",
                           nrow(df), " tested"),
         x = NULL, y = expression(-log[10](p~value))) +
    theme_classic() +
    theme(plot.subtitle = element_text(size = 11, face = "bold"))
}

panels <- Map(make_panel, GROUPS, LABELS)

# No title/caption in the figure: they are written in the preprint legend
p_final <- (panels[[1]] | panels[[2]]) +
  plot_layout(guides = "collect")

ggsave(file.path(OUT_FIG, paste0("17_", HISTOLOGY, "_CollecTRI_dotplot.png")),
       plot = p_final, width = 12, height = 6.5, dpi = 300, bg = "white")

message("Script17 (CollecTRI dot plot, ", HISTOLOGY, ") complete.")