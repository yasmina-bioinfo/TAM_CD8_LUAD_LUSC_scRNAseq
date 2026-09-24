# ============================================================
# Bloc6_CD8 — 05_LUAD_CD8_Composition.R
# Purpose: statistical comparison and visualization of ProjecTILs CD8
# functional-state composition (MPR vs non-MPR), on the purified CD8
# population (scGate + UCell, Scripts 01-02; ProjecTILs, Script 03).
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

IN_OBJ   <- here("Objects", "03_LUAD_seu_projectils.rds")
REF_FILE <- here("References", "CD8T_human_ref_v1.rds")
OUT_TAB  <- here("Results", "Bloc6_CD8", "Tables")
OUT_FIG  <- here("Results", "Bloc6_CD8", "Figures")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

# Cells and patients expected from Script 03 (guard against a wrong input file)
EXPECTED_CELLS    <- 47996                          # LUSC: 93584
EXPECTED_PATIENTS <- c("MPR" = 10, "non-MPR" = 40)  # LUSC: c("MPR" = 32, "non-MPR" = 68)
MIN_CD8_PER_PATIENT <- 50  # documented threshold: below this, a rare
# state (~1-3% of CD8, e.g. CD8.TPEX) has a high chance of a false zero
# purely from under-sampling, not biology.

# ---- Load reference atlas (for official color palette) ----
message("Loading reference atlas (for official color palette)...")
ref.cd8 <- load.reference.map(REF_FILE)
atlas_palette <- ref.cd8@misc$atlas.palette

# ---- Load LUAD ProjecTILs object ----
message("Loading LUAD ProjecTILs object... ", Sys.time())
seu <- readRDS(IN_OBJ)
stopifnot(all(c("pathological_response", "functional.cluster", "sampleID") %in%
                colnames(seu@meta.data)))

# ---- Population: already purified by scGate + UCell (Scripts 01-02) ----
cd8 <- seu
message("CD8 cells: ", ncol(cd8))
stopifnot(ncol(cd8) == EXPECTED_CELLS)

obs_pat <- cd8@meta.data %>%
  distinct(sampleID, pathological_response) %>%
  count(pathological_response)
print(obs_pat)
obs_vec <- setNames(obs_pat$n, obs_pat$pathological_response)
stopifnot(all(obs_vec[names(EXPECTED_PATIENTS)] == EXPECTED_PATIENTS))

# ============================================================
# PATIENT-LEVEL TEST (Wilcoxon rank-sum, BH-corrected)
# Rationale: treating each cell as an independent observation would
# ignore within-patient correlation (pseudoreplication). Following the
# same approach as Hu et al., 2023, Genome Medicine (comparable NSCLC
# neoadjuvant anti-PD-1 scRNA-seq study): compute per-patient state
# proportions, then Wilcoxon rank-sum test, BH-corrected across states.
# This is the test used for reported conclusions.
# ============================================================

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

n_patients_retained <- n_distinct(patient_prop$sampleID)
message("Patients retained (>= ", MIN_CD8_PER_PATIENT, " CD8 cells): ", n_patients_retained)
print(patient_prop %>% distinct(sampleID, pathological_response) %>%
        count(pathological_response))

write_csv(patient_prop, file.path(OUT_TAB, "05_LUAD_CD8_patient_level_proportions.csv"))

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
write_csv(wilcoxon_results, file.path(OUT_TAB, "05_LUAD_CD8_patientlevel_wilcoxon.csv"))

# ---- Composition table: functional.cluster x pathological_response ----
# Restricted to the same patients as the patient-level test (>= 50 CD8
# cells), so the figure and the statistics describe the same population.
comp_table <- patient_counts %>%
  filter(sampleID %in% patient_prop$sampleID) %>%
  left_join(cd8@meta.data %>% distinct(sampleID, pathological_response), by = "sampleID") %>%
  group_by(functional.cluster, pathological_response) %>%
  summarise(n_cells = sum(n_state), .groups = "drop")
message("Composition table (cells, patients with >= ", MIN_CD8_PER_PATIENT, " CD8 cells only):")
print(comp_table)

# ---- Proportions per condition, for plotting ----
df_prop <- comp_table %>%
  rename(functional_cluster = functional.cluster) %>%
  group_by(pathological_response) %>%
  mutate(prop = n_cells / sum(n_cells)) %>%
  ungroup()

write_csv(df_prop, file.path(OUT_TAB, "05_LUAD_CD8_composition_by_response.csv"))

# ---- Figure A: composition barplot (legend marks patient-level
# significant states; no single caption statistic, since the
# patient-level test has no cell-level-style global omnibus equivalent) ----
sig_states <- wilcoxon_results$functional_cluster[wilcoxon_results$p_adj_BH < 0.05]

legend_labels <- setNames(
  ifelse(names(atlas_palette) %in% sig_states,
         paste0(names(atlas_palette), " *"),
         names(atlas_palette)),
  names(atlas_palette)
)

caption_text <- if (length(sig_states) > 0) {
  "* patient-level Wilcoxon test, BH-adjusted P < 0.05"
} else {
  "No functional state significantly different (patient-level Wilcoxon test, BH-adjusted P \u2265 0.05)"
}

p_bar <- ggplot(df_prop, aes(x = pathological_response, y = prop, fill = functional_cluster)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = atlas_palette, labels = legend_labels, drop = TRUE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.05)) +
  ylab("CD8 functional-state proportion") +
  labs(caption = caption_text) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.35, "cm"),
        legend.title = element_blank(),
        plot.caption = element_text(size = 8, face = "italic", hjust = 0)) +
  guides(fill = guide_legend(ncol = 1))

ggsave(file.path(OUT_FIG, "05_LUAD_CD8_composition_barplot.png"), p_bar,
       width = 8, height = 9, dpi = 300, bg = "white")
message("Saved: composition barplot")

# ---- Figure B: CD8.TEX per patient (boxplot + points), the focal
# state of this study. CD8.TPEX is not tested at the patient level
# (see 03b_LUAD_state_testability) and is not plotted here. ----
tex_prop <- patient_prop %>% filter(functional.cluster == "CD8.TEX")
tex_p_display <- wilcoxon_results$p_display[wilcoxon_results$functional_cluster == "CD8.TEX"]

p_tex <- ggplot(tex_prop, aes(x = pathological_response, y = prop, fill = pathological_response)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.7) +
  geom_jitter(width = 0.08, size = 1.5, alpha = 0.7, color = "black") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("MPR" = "#0072B2", "non-MPR" = "#D55E00")) +
  ylab("CD8.TEX (% of CD8 cells)") +
  labs(caption = paste0("Wilcoxon, BH-adjusted ", tex_p_display)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 13, face = "bold"),
        legend.position = "none",
        plot.caption = element_text(size = 8, face = "italic", hjust = 0))

ggsave(file.path(OUT_FIG, "05_LUAD_CD8_TEX_boxplot.png"), p_tex,
       width = 4, height = 5, dpi = 300, bg = "white")
message("Saved: CD8.TEX boxplot")

message("Script 05 (CD8 composition by response, LUAD) complete. ", Sys.time())