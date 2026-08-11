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
  library(tidyr)
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

# ---- Global chi-squared test (cell-level; kept for traceability, not
# used for reported conclusions — see patient-level test below) ----
chisq_res <- chisq.test(comp_table)
message("Global chi-squared test (cell-level, archival):")
print(chisq_res)

p_display <- scales::pvalue(chisq_res$p.value, accuracy = 0.001, add_p = TRUE)

# ---- Proportions per condition, for plotting ----
df_prop <- as.data.frame(comp_table)
colnames(df_prop) <- c("functional_cluster", "pathological_response", "n_cells")

df_prop <- df_prop %>%
  group_by(pathological_response) %>%
  mutate(prop = n_cells / sum(n_cells)) %>%
  ungroup()

write_csv(df_prop, file.path(OUT_TAB, "05_LUSC_CD8_composition_by_MPR.csv"))

# ---- Targeted prop.test per functional state (cell-level; archival) ----
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
message("Targeted prop.test results (cell-level, archival):")
print(prop_results)
write_csv(prop_results, file.path(OUT_TAB, "05_LUSC_CD8_proptest_by_state.csv"))

# ============================================================
# PATIENT-LEVEL TEST (Wilcoxon rank-sum, BH-corrected)
# Rationale: cell-level tests above treat each cell as an independent
# observation despite within-patient correlation (pseudoreplication).
# Following the same approach as Hu et al., 2023, Genome Medicine
# (comparable NSCLC neoadjuvant anti-PD-1 scRNA-seq study): compute
# per-patient state proportions, then Wilcoxon rank-sum test, BH-
# corrected across states. This is the test used for reported
# conclusions; cell-level results above are kept for traceability only.
# ============================================================

MIN_CD8_PER_PATIENT <- 50  # documented threshold: below this, a rare
# state (~1-3% of CD8, e.g. CD8.TPEX) has a high chance of a false zero
# purely from under-sampling, not biology

message("Computing per-patient proportions (min ", MIN_CD8_PER_PATIENT,
        " CD8 cells per patient)... ", Sys.time())

patient_totals <- cd8@meta.data %>%
  count(sampleID, name = "n_total_cd8")

patient_counts <- cd8@meta.data %>%
  mutate(functional.cluster = factor(functional.cluster)) %>%
  count(sampleID, functional.cluster, .drop = FALSE, name = "n_state") %>%
  tidyr::complete(sampleID, functional.cluster, fill = list(n_state = 0))

patient_prop <- patient_counts %>%
  left_join(patient_totals, by = "sampleID") %>%
  filter(n_total_cd8 >= MIN_CD8_PER_PATIENT) %>%
  mutate(prop = n_state / n_total_cd8) %>%
  left_join(
    cd8@meta.data %>% distinct(sampleID, pathological_response),
    by = "sampleID"
  )

message("Patients retained: ", n_distinct(patient_prop$sampleID))
print(patient_prop %>% distinct(sampleID, pathological_response) %>%
        count(pathological_response))

write_csv(patient_prop, file.path(OUT_TAB, "05_LUSC_CD8_patient_level_proportions.csv"))

wilcoxon_results <- lapply(levels(factor(patient_prop$functional.cluster)), function(state) {
  sub <- patient_prop %>% filter(functional.cluster == state)
  mpr_vals    <- sub$prop[sub$pathological_response == "MPR"]
  nonmpr_vals <- sub$prop[sub$pathological_response == "non-MPR"]
  test <- wilcox.test(nonmpr_vals, mpr_vals)
  data.frame(
    functional_cluster = state,
    n_MPR         = length(mpr_vals),
    n_nonMPR      = length(nonmpr_vals),
    median_MPR    = median(mpr_vals),
    median_nonMPR = median(nonmpr_vals),
    p_value       = test$p.value
  )
})
wilcoxon_results <- bind_rows(wilcoxon_results)
wilcoxon_results$p_adj_BH <- p.adjust(wilcoxon_results$p_value, method = "BH")
wilcoxon_results$p_display <- scales::pvalue(wilcoxon_results$p_adj_BH, accuracy = 0.001, add_p = TRUE)

message("Patient-level Wilcoxon results (BH-corrected, reported conclusions):")
print(wilcoxon_results)
write_csv(wilcoxon_results, file.path(OUT_TAB, "05_LUSC_CD8_patientlevel_wilcoxon.csv"))

# ---- Barplot: legend marks states significant at patient level (no
# single caption statistic — patient-level test has no cell-level-style
# global omnibus equivalent; see FICHE_TECHNIQUE_OUTILS) ----
sig_states <- wilcoxon_results$functional_cluster[wilcoxon_results$p_adj_BH < 0.05]

legend_labels <- setNames(
  ifelse(names(atlas_palette) %in% sig_states,
         paste0(names(atlas_palette), " *"),
         names(atlas_palette)),
  names(atlas_palette)
)

p_bar <- ggplot(df_prop, aes(x = pathological_response, y = prop, fill = functional_cluster)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = atlas_palette, labels = legend_labels, drop = TRUE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.05)) +
  ylab("CD8 functional-state proportion") +
  labs(caption = "* patient-level Wilcoxon test, BH-adjusted P < 0.05") +
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