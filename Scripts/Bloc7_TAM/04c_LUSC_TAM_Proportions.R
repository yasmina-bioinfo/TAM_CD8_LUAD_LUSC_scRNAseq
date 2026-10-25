# ============================================================
# Bloc7_TAM — Script4c_TAM_Proportions_LUSC.R
# Purpose: per-patient proportion of Macrophage_pure among the scGate-pure
# myeloid cells (Macrophage_pure + Monocyte_pure), MPR vs non-MPR.
# Same logic as Bloc6/05_LUSC_CD8_Composition.R: denominator = the patient's
# cells in the purified population; patient-level Wilcoxon rank-sum test.
# The Monocyte_pure proportion is its complement, so a single test is run
# per histology (no multiplicity correction needed). scGate_unresolved and
# Ambiguous_double_pure cells are excluded, as in all downstream analyses.
# No direct comparison between histologies is made.
#
# Input : Objects/07_LUSC_seu_TAM_scored.rds  (check object and column names)
# Output: Results/Bloc7_TAM/Tables/18_LUSC_TAM_proportions_*.csv
#         Results/Bloc7_TAM/Figures/18_LUSC_TAM_proportions_boxplot.png
# ============================================================

suppressPackageStartupMessages({
  library(here); library(Seurat); library(dplyr); library(tidyr)
  library(readr); library(ggplot2); library(scales)
})

HISTOLOGY <- "LUSC"
MIN_PURE_PER_PATIENT <- 50   # same threshold as the CD8 composition analysis
# Macrophage_pure + Monocyte_pure cells (Supplementary Table 3)
EXPECTED_PURE <- c(LUAD = 27550, LUSC = 33083)
COL_MPR <- "#0072B2"; COL_NONMPR <- "#D55E00"

IN_OBJ  <- here("Objects", paste0("04b_", HISTOLOGY, "_seu_TAM_classified.rds"))
OUT_TAB <- here("Results", "Bloc7_TAM", "Tables")
OUT_FIG <- here("Results", "Bloc7_TAM", "Figures")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Metadata and input guard ---------------------------------------
seu <- readRDS(IN_OBJ)
md  <- seu@meta.data
rm(seu); gc()
stopifnot(all(c("scGate_TAM_class", "sampleID", "pathological_response") %in%
                colnames(md)))

# Descriptive check (console only, not tested): share of myeloid cells excluded
# from the pure population (scGate_unresolved + Ambiguous_double_pure), per
# patient, then median by response group
md %>%
  filter(pathological_response %in% c("MPR", "non-MPR")) %>%
  group_by(sampleID, pathological_response) %>%
  summarise(prop_excluded = mean(!scGate_TAM_class %in%
                                   c("Macrophage_pure", "Monocyte_pure")),
            .groups = "drop") %>%
  group_by(pathological_response) %>%
  summarise(n_patients = n(), median_prop_excluded = median(prop_excluded)) %>%
  print()

md <- md %>% filter(scGate_TAM_class %in% c("Macrophage_pure", "Monocyte_pure"))
stopifnot(nrow(md) == EXPECTED_PURE[[HISTOLOGY]])

md <- md %>% filter(pathological_response %in% c("MPR", "non-MPR"))

# ---- 2. Per-patient proportions ----------------------------------------
patient_counts <- md %>%
  count(sampleID, pathological_response, scGate_TAM_class) %>%
  tidyr::complete(tidyr::nesting(sampleID, pathological_response),
                  scGate_TAM_class = c("Macrophage_pure", "Monocyte_pure"),
                  fill = list(n = 0))

patient_tot <- patient_counts %>%
  group_by(sampleID, pathological_response) %>%
  summarise(n_pure  = sum(n),
            n_macro = n[scGate_TAM_class == "Macrophage_pure"],
            .groups = "drop") %>%
  mutate(prop_macro = n_macro / n_pure)

# Impact of candidate thresholds (decide on patient numbers, not on results)
for (thr in c(10, 25, 50)) {
  cat("\nPatients with >=", thr, "pure myeloid cells:\n")
  print(table(patient_tot$pathological_response[patient_tot$n_pure >= thr]))
}

patient_prop <- patient_tot %>% filter(n_pure >= MIN_PURE_PER_PATIENT)
cat("\nPatients retained (>=", MIN_PURE_PER_PATIENT, "pure cells):\n")
print(table(patient_prop$pathological_response))

write_csv(patient_prop,
          file.path(OUT_TAB, paste0("18_", HISTOLOGY, "_TAM_proportions_patient_level.csv")))

# ---- 3. Patient-level Wilcoxon test ------------------------------------
mpr_vals    <- patient_prop$prop_macro[patient_prop$pathological_response == "MPR"]
nonmpr_vals <- patient_prop$prop_macro[patient_prop$pathological_response == "non-MPR"]
test <- wilcox.test(nonmpr_vals, mpr_vals)
test_mono <- wilcox.test(1 - nonmpr_vals, 1 - mpr_vals)
cat("p Macrophage_pure:", test$p.value, " | p Monocyte_pure:", test_mono$p.value, "\n")

res <- data.frame(
  histology     = HISTOLOGY,
  n_MPR         = length(mpr_vals),
  n_nonMPR      = length(nonmpr_vals),
  median_MPR    = median(mpr_vals),
  median_nonMPR = median(nonmpr_vals),
  p_value       = test$p.value,
  p_display     = scales::pvalue(test$p.value, accuracy = 0.001, add_p = TRUE)
)
print(res)
write_csv(res, file.path(OUT_TAB, paste0("18_", HISTOLOGY, "_TAM_proportions_wilcoxon.csv")))

# ---- 4. Boxplot, both populations --------------------------------------
# Monocyte_pure is the complement of Macrophage_pure (same cells, same test):
# the same Wilcoxon result is annotated in each panel.
plot_df <- patient_prop %>%
  mutate(Macrophage_pure = prop_macro,
         Monocyte_pure   = 1 - prop_macro) %>%
  pivot_longer(c(Macrophage_pure, Monocyte_pure),
               names_to = "population", values_to = "prop") %>%
  mutate(pathological_response = factor(pathological_response,
                                        levels = c("MPR", "non-MPR")),
         population = factor(population,
                             levels = c("Macrophage_pure", "Monocyte_pure")))

lab_df <- data.frame(
  population = factor(c("Macrophage_pure", "Monocyte_pure"),
                      levels = c("Macrophage_pure", "Monocyte_pure")),
  label = paste0("Wilcoxon, ", res$p_display)
)

p <- ggplot(plot_df, aes(x = pathological_response, y = prop,
                         fill = pathological_response)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.7) +
  geom_jitter(width = 0.08, size = 1.5, alpha = 0.7, color = "black") +
  geom_text(data = lab_df, aes(x = 1.5, y = 1.07, label = label),
            inherit.aes = FALSE, size = 3.3, fontface = "italic") +
  facet_wrap(~ population) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     breaks = seq(0, 1, 0.25), limits = c(0, 1.1)) +
  scale_fill_manual(values = c("MPR" = COL_MPR, "non-MPR" = COL_NONMPR)) +
  ylab("% of pure myeloid cells") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 13, face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 11),
        legend.position = "none")

ggsave(file.path(OUT_FIG, paste0("18_", HISTOLOGY, "_TAM_proportions_boxplot.png")),
       p, width = 7, height = 5, dpi = 300, bg = "white")

message("Script18 (TAM proportions, ", HISTOLOGY, ") complete.")