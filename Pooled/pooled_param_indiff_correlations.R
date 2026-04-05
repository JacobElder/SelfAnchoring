#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# Pooled parameter × individual-difference correlations
# Mirrors the pooled CCA preprocessing:
# 1. Load study-level individual-differences tables
# 2. Keep shared parameter and personality variables
# 3. Z-score parameters and scales within study
# 4. Compute pooled Pearson correlations across the stacked standardized data
# 5. Apply FDR correction separately for mechanism vs. model-preference families

scale_vars <- c("DS", "Proto", "SCC", "SI", "RSE", "NTB", "NFC", "SING.Ind", "SING.Inter")
mech_params <- c("m", "bias", "lambda", "subject_mcr")
elpd_params <- "delta_elpd"

s1 <- read.csv(here("Results", "ind_diffs_s1_full.csv")) %>%
  mutate(
    study = "Study 1 (Minimal Groups)",
    condition = "Minimal Groups"
  )

s2_file <- here("Results", "ind_diffs_s2_full_enriched.csv")
s2 <- read.csv(if (file.exists(s2_file)) s2_file else here("Results", "ind_diffs_s2_full.csv")) %>%
  mutate(
    study = "Study 2 (University Groups)",
    condition = case_when(
      outgroup == "UCLA" ~ "High-Status",
      outgroup == "Not UCR" ~ "Negation",
      outgroup == "CSU LA" ~ "Low-Status",
      TRUE ~ as.character(outgroup)
    )
  )

s3_file <- here("Results", "ind_diffs_s3_full_enriched.csv")
s3 <- read.csv(if (file.exists(s3_file)) s3_file else here("Results", "ind_diffs_s3_full.csv")) %>%
  mutate(
    study = "Study 3 (Racial Groups)",
    condition = as.character(condition)
  )

keep_cols <- c(
  "subID", "study", "condition",
  intersect(c(mech_params, elpd_params), union(names(s1), union(names(s2), names(s3)))),
  intersect(scale_vars, union(names(s1), union(names(s2), names(s3))))
)

df <- bind_rows(
  s1 %>% select(any_of(keep_cols)),
  s2 %>% select(any_of(keep_cols)),
  s3 %>% select(any_of(keep_cols))
)

avail_mech_params <- intersect(mech_params, names(df))
avail_elpd_params <- intersect(elpd_params, names(df))
avail_scales <- intersect(scale_vars, names(df))

z_cols <- c(avail_mech_params, avail_elpd_params, avail_scales)

df_z <- df %>%
  group_by(study) %>%
  mutate(across(all_of(z_cols), ~ as.numeric(scale(.)))) %>%
  ungroup()

run_cors <- function(param_set, dat, scale_cols, family_label) {
  rows <- list()
  for (p in param_set) {
    for (s in scale_cols) {
      cc <- complete.cases(dat[, c(p, s)])
      ct <- cor.test(dat[[p]][cc], dat[[s]][cc], method = "pearson")
      rows[[length(rows) + 1]] <- tibble(
        family = family_label,
        param = p,
        scale = s,
        r = unname(ct$estimate),
        p_raw = ct$p.value,
        n = sum(cc)
      )
    }
  }
  bind_rows(rows)
}

res_mech <- run_cors(avail_mech_params, df_z, avail_scales, "mechanism")
res_elpd <- run_cors(avail_elpd_params, df_z, avail_scales, "model_preference")

if (nrow(res_mech) > 0) {
  res_mech <- res_mech %>% mutate(p_fdr = p.adjust(p_raw, method = "BH"))
}
if (nrow(res_elpd) > 0) {
  res_elpd <- res_elpd %>% mutate(p_fdr = p.adjust(p_raw, method = "BH"))
}

res <- bind_rows(res_mech, res_elpd) %>%
  arrange(p_raw) %>%
  mutate(
    r = round(r, 3),
    p_raw = round(p_raw, 4),
    p_fdr = round(p_fdr, 4)
  )

write.csv(df_z, here("Results", "ind_diffs_pooled_standardized.csv"), row.names = FALSE)
write.csv(res, here("Results", "correlations_pooled_standardized.csv"), row.names = FALSE)

message("Saved: Results/ind_diffs_pooled_standardized.csv")
message("Saved: Results/correlations_pooled_standardized.csv")
