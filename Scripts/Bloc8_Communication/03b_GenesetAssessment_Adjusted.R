#-------------------------------------------------------------------------------
#Script03b_GenesetAssessment_Adjusted.R
# Geneset-to-background ratio with raw vs adjusted p-values, from the saved
# Script03 output. No DE is rerun and no existing file is overwritten.
# The background (n_background) is taken from the table saved by Script03,
# because it does not depend on the p-value choice.
#-------------------------------------------------------------------------------

library(here); library(dplyr); library(readr); library(multinichenetr)

HISTOLOGY <- "LUAD"   # replace by "LUSC" for the second run
OUT_TAB   <- here("Results", "Bloc8_Communication", "Tables")

f_out <- here("Objects", paste0("08_", HISTOLOGY, "_multinichenet_output.rds"))
f_ref <- file.path(OUT_TAB, paste0("03_", HISTOLOGY, "_geneset_assessment.csv"))
stopifnot(file.exists(f_out), file.exists(f_ref))

out <- readRDS(f_out)
stopifnot("celltype_de is missing from the saved output" = "celltype_de" %in% names(out))
celltype_de <- out$celltype_de
ref <- read_csv(f_ref, show_col_types = FALSE)

contrasts <- c("MPR-non.MPR", "non.MPR-MPR")
assess <- function(p_val_adj) {
  contrasts %>%
    lapply(process_geneset_data, celltype_de, 0.50, p_val_adj, 0.05) %>%
    bind_rows()
}

bg <- ref %>% select(contrast, cluster_id, n_background_ref = n_background) %>% distinct()

res <- bind_rows(assess(FALSE), assess(TRUE)) %>%
  left_join(bg, by = c("contrast", "cluster_id")) %>%
  mutate(n_background       = n_background_ref,
         prop_geneset_up    = n_geneset_up   / n_background,
         prop_geneset_down  = n_geneset_down / n_background,
         in_range_up        = prop_geneset_up   >= 1/200 & prop_geneset_up   <= 1/10,
         in_range_down      = prop_geneset_down >= 1/200 & prop_geneset_down <= 1/10) %>%
  select(-n_background_ref)

# Control: the raw version must reproduce the table saved by Script03
chk <- res %>% filter(!adjusted) %>%
  inner_join(ref, by = c("contrast", "cluster_id"), suffix = c("", "_ref"))
stopifnot("raw rows missing" = nrow(chk) == 6,
          "gene counts differ from Script03" =
            all(chk$n_geneset_up == chk$n_geneset_up_ref) &&
            all(chk$n_geneset_down == chk$n_geneset_down_ref),
          "ratios differ from Script03" =
            isTRUE(all.equal(chk$prop_geneset_up, chk$prop_geneset_up_ref, tolerance = 1e-6)) &&
            isTRUE(all.equal(chk$prop_geneset_down, chk$prop_geneset_down_ref, tolerance = 1e-6)),
          "in_range flags differ from Script03" =
            all(chk$in_range_up == chk$in_range_up_ref) &&
            all(chk$in_range_down == chk$in_range_down_ref))

write_csv(res, file.path(OUT_TAB,
                         paste0("03b_", HISTOLOGY, "_geneset_assessment_raw_vs_adjusted.csv")))

print(res %>% select(adjusted, contrast, cluster_id, n_background, n_geneset_up,
                     n_geneset_down, prop_geneset_up, prop_geneset_down,
                     in_range_up, in_range_down), n = Inf, width = Inf)
cat("Script03b (", HISTOLOGY, ") complete\n")