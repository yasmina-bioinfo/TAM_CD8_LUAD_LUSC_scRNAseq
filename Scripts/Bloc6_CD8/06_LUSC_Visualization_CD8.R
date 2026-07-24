# ============================================================
# Bloc6_CD8 — Script05_Composition_MPR_LUSC.R
# Purpose: statistical comparison and visualization of ProjecTILs CD8
# functional-state composition (MPR vs non-MPR), restricted to the
# revised CD8 lineage (scGate + UCell + ProjecTILs cross-validation,
# including the Treg->CD8 reclassification from Script04).
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(ProjecTILs)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
})

IN_OBJ   <- here("Objects", "03_LUSC_seu_projectils.rds")
REF_FILE <- here("References", "CD8T_human_ref_v1.rds")
OUT_TAB  <- here("Results", "Bloc6_CD8", "Tables")
OUT_FIG  <- here("Results", "Bloc6_CD8", "Figures")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

# ---- Load reference atlas (for official color palette) ----
message("Loading reference atlas (for official color palette)...")
ref.cd8 <- load.reference.map(REF_FILE)
atlas_palette <- ref.cd8@misc$atlas.palette

# ---- Load LUSC ProjecTILs object ----
message("Loading LUSC ProjecTILs object... ", Sys.time())
seu <- readRDS(IN_OBJ)
stopifnot(all(c("lineage_CD8_revised", "pathological_response", "functional.cluster") %in%
                colnames(seu@meta.data)))

# ---- Restrict to revised CD8 lineage ----
cd8 <- subset(seu, subset = lineage_CD8_revised == TRUE)
message("CD8 cells (revised lineage): ", ncol(cd8))
message("pathological_response levels:")
print(table(cd8$pathological_response, useNA = "ifany"))

# ---- Composition table: functional.cluster x pathological_response ----
comp_table <- table(cd8$functional.cluster, cd8$pathological_response)
message("Composition table:")
print(comp_table)

# ---- Global chi-squared test ----
chisq_res <- chisq.test(comp_table)
message("Global chi-squared test (functional.cluster x pathological_response):")
print(chisq_res)

# Standard p-value display format (consistent with Bloc5 TME barplot,
# to reuse for every barplot going forward)
p_display <- scales::pvalue(chisq_res$p.value, accuracy = 0.001, add_p = TRUE)

# ---- Proportions per condition, for plotting ----
df_prop <- as.data.frame(comp_table)
colnames(df_prop) <- c("functional_cluster", "pathological_response", "n_cells")

df_prop <- df_prop %>%
  group_by(pathological_response) %>%
  mutate(prop = n_cells / sum(n_cells)) %>%
  ungroup()

write_csv(df_prop, file.path(OUT_TAB, "05_LUSC_CD8_composition_by_MPR.csv"))

# ---- Targeted prop.test per functional state ----
prop_results <- lapply(levels(factor(df_prop$functional_cluster)), function(state) {
  sub <- comp_table[state, ]
  total_by_cond <- colSums(comp_table)
  test <- prop.test(x = sub, n = total_by_cond)
  data.frame(
    functional_cluster = state,
    p_value = test$p.value,
    estimate_diff = diff(test$estimate)
  )
})
prop_results <- bind_rows(prop_results)
message("Targeted prop.test results (per functional state):")
print(prop_results)
write_csv(prop_results, file.path(OUT_TAB, "05_LUSC_CD8_proptest_by_state.csv"))

# ---- Barplot (same visual model as Bloc5 TME barplot) ----
p_bar <- ggplot(df_prop, aes(x = pathological_response, y = prop, fill = functional_cluster)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = atlas_palette, drop = TRUE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.05)) +
  ylab("CD8 functional-state proportion") +
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

ggsave(file.path(OUT_FIG, "05_LUSC_CD8_composition_barplot.png"), p_bar,
       width = 8, height = 6, dpi = 300, bg = "white")

message("Saved: barplot")
message("Script05 (CD8 composition by MPR, LUSC) complete. ", Sys.time())