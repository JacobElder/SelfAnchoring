# RESCUE: Study 2 — minimal version using only saved LOO RDS files
# Does NOT need chain CSVs. Handles: model comparison, pareto_k, diagnostics.
# params_ind/summary for asym_lambda will need a re-run (secondary priority).

library(loo)
library(tidyverse)
library(here)

message("Loading saved LOO objects...")
loo_bias        <- readRDS(here("Fits", "loo_s2_bias.rds"))
loo_symmetric   <- readRDS(here("Fits", "loo_s2_symmetric.rds"))
loo_sym_lambda  <- readRDS(here("Fits", "loo_s2_sym_lambda.rds"))
loo_asym_lambda <- readRDS(here("Fits", "loo_s2_asym_lambda.rds"))

# 1. Model comparison
comp <- loo_compare(list(
  bias        = loo_bias,
  symmetric   = loo_symmetric,
  sym_lambda  = loo_sym_lambda,
  asym_lambda = loo_asym_lambda
))
print(comp)
write.csv(as.data.frame(comp), here("Results", "model_comparison_s2_results.csv"))
message("model_comparison_s2_results.csv written.")

# 2. Per-subject Pareto k for asym_lambda (from saved LOO)
fulldf  <- read.csv(here("Study 2/Cleaning/output/fullTest.csv")) %>% filter(!is.na(ingChoiceN))
traindf <- read.csv(here("Study 2/Cleaning/output/fullTrain_fixed.csv")) %>% filter(!is.na(selfResp))
common_ids  <- intersect(unique(fulldf$subID), unique(traindf$subID))
fulldf      <- fulldf %>% filter(subID %in% common_ids)
uIds        <- sort(common_ids)
maxTrials   <- max(fulldf$trialTotalT2)
nTrials_vec <- sapply(uIds, function(id) nrow(filter(fulldf, subID == id)))

compute_subj_pareto_k <- function(l, study_label, model_name, uIds, nTrials_vec, maxTrials) {
  pk_vec <- l$diagnostics$pareto_k
  rows <- lapply(seq_along(uIds), function(i) {
    idx <- ((i - 1) * maxTrials + 1):((i - 1) * maxTrials + nTrials_vec[i])
    k_i <- pk_vec[idx]
    data.frame(
      study        = study_label,
      model        = model_name,
      subID        = uIds[i],
      subj_idx     = i,
      n_trials     = nTrials_vec[i],
      k_mean       = round(mean(k_i), 4),
      k_max        = round(max(k_i), 4),
      n_good       = sum(k_i < 0.5),
      n_ok         = sum(k_i >= 0.5 & k_i < 0.7),
      n_bad        = sum(k_i >= 0.7 & k_i < 1.0),
      n_verybad    = sum(k_i >= 1.0),
      pct_reliable = round(100 * mean(k_i < 0.7), 1),
      concern      = ifelse(any(k_i >= 1.0), "HIGH",
                    ifelse(any(k_i >= 0.7), "MODERATE",
                    ifelse(any(k_i >= 0.5), "LOW", "NONE")))
    )
  })
  do.call(rbind, rows)
}

subj_pk <- compute_subj_pareto_k(loo_asym_lambda, "S2", "asym_lambda",
                                  uIds, nTrials_vec, maxTrials)
write.csv(subj_pk, here("Results", "pareto_k_subj_s2_asym_lambda.csv"), row.names = FALSE)
message("pareto_k_subj_s2_asym_lambda.csv written.")

# 3. Print LOO summary for all 4 models
message("\n=== LOO Summary ===")
for (nm in c("bias", "symmetric", "sym_lambda", "asym_lambda")) {
  l <- get(paste0("loo_", nm))
  pk <- l$diagnostics$pareto_k
  cat(sprintf("%-15s ELPD=%.1f (SE=%.1f) | pk_ok=%.1f%%\n",
    nm, l$estimates["elpd_loo","Estimate"], l$estimates["elpd_loo","SE"],
    100 * mean(pk < 0.7)))
}

message("\nMinimal rescue complete.")
message("NOTE: params_ind_s2_asym_lambda.csv and summary_s2_asym_lambda.csv")
message("      will need a re-run (chain CSVs had 0 data rows).")
