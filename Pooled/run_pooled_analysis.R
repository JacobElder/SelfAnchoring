# Pooled Model Comparison Analysis
# Optimized for Speed (cmdstanr) and Portability (CSV Export)
# Memory-Efficient Version: Saves each model's results to disk and reloads at the end.

library(cmdstanr)
library(tidyverse)
library(loo)
library(here)
library(parallel)

# 1. Setup
options(mc.cores = parallel::detectCores())
stan_data <- readRDS(here("Fits", "pooled_stan_data.rds"))

# Helper: per-subject Pareto k summary
compute_subj_pareto_k <- function(l, study_label, model_name, uIds, nTrials_vec, maxTrials) {
  pk <- l$diagnostics$pareto_k
  rows <- lapply(seq_along(uIds), function(i) {
    idx <- ((i - 1) * maxTrials + 1):((i - 1) * maxTrials + nTrials_vec[i])
    k_i <- pk[idx]
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

# Recover subject indices
pooled_uIds <- if (!is.null(stan_data$subIDs)) stan_data$subIDs else seq_len(stan_data$nSubjects)
pooled_maxTrials <- max(stan_data$nTrials)

# 2. Define Models
models <- list(
  bias = "S_Pooled_Bias.stan",
  symmetric = "S_Pooled_Symmetric.stan",
  sym_lambda = "S_Pooled_Sym_Lambda.stan",
  asym_lambda = "S_Pooled_Asym_Lambda.stan"
)

# 3. Fit function (now saves to disk and returns nothing to save memory)
fit_and_save_pooled <- function(model_name) {
  message(paste("Fitting Pooled model:", model_name))
  mod <- cmdstan_model(here("Computational Models", models[[model_name]]))
  
  fit <- mod$sample(
    data = stan_data, 
    seed = 1234, 
    chains = 4, 
    parallel_chains = 2,
    iter_warmup = 1000,
    iter_sampling = 2000,
    adapt_delta = 0.99,
    max_treedepth = 12,
    init = 0,
    refresh = 100
  )
  
  # A. Save LOO results
  l <- fit$loo()
  saveRDS(l, here("Fits", paste0("loo_pooled_", model_name, ".rds")))
  
  # B. Save Summary
  all_vars    <- fit$metadata()$stan_variables
  struct_vars <- all_vars[!all_vars %in% c("log_lik", "p_pred", "mcr")]
  sum_fit <- fit$summary(variables = struct_vars)
  write.csv(sum_fit, here("Results", paste0("summary_pooled_", model_name, ".csv")), row.names = FALSE)
  
  # C. Save Individual Parameters
  ind_params <- sum_fit %>%
    filter(str_detect(variable, "\\[")) %>%
    select(variable, median, rhat, ess_bulk)
  write.csv(ind_params, here("Results", paste0("params_ind_pooled_", model_name, ".csv")), row.names = FALSE)
  
  # D. Save Diagnostics
  diag <- fit$diagnostic_summary()
  pk   <- l$diagnostics$pareto_k
  diag_df <- data.frame(
    model           = model_name,
    max_rhat        = max(sum_fit$rhat, na.rm = TRUE),
    num_divergent   = sum(diag$num_divergent),
    pct_divergent   = round(100 * sum(diag$num_divergent) / (4 * 2000), 3),
    min_ess_bulk    = min(sum_fit$ess_bulk, na.rm = TRUE),
    converged       = (max(sum_fit$rhat, na.rm = TRUE) < 1.01 & sum(diag$num_divergent) == 0),
    pk_good         = sum(pk < 0.5),
    pk_ok           = sum(pk >= 0.5 & pk < 0.7),
    pk_bad          = sum(pk >= 0.7 & pk < 1.0),
    pk_verybad      = sum(pk >= 1.0),
    pk_pct_ok       = round(100 * mean(pk < 0.7), 1),
    pk_max          = round(max(pk), 3),
    loo_reliable    = mean(pk < 0.7) > 0.9
  )
  write.csv(diag_df, here("Results", paste0("diag_pooled_", model_name, ".csv")), row.names = FALSE)

  # E. Save Per-subject Pareto k
  subj_pk <- compute_subj_pareto_k(l, "Pooled", model_name, pooled_uIds, stan_data$nTrials, pooled_maxTrials)
  write.csv(subj_pk, here("Results", paste0("pareto_k_subj_pooled_", model_name, ".csv")), row.names = FALSE)

  # Explicitly clear memory
  rm(fit, l, sum_fit, ind_params, diag, pk, diag_df, subj_pk)
  gc()
  
  return(NULL)
}

# 4. Run Analysis (Loop and save)
for(m in names(models)) {
  fit_and_save_pooled(m)
}

# 5. Model Comparison (Load from disk)
message("All models fit. Comparing results from saved files.")
loo_files <- list.files(path = here("Fits"), pattern = "loo_pooled_.*\\.rds", full.names = TRUE)
diag_files <- list.files(path = here("Results"), pattern = "diag_pooled_.*\\.csv", full.names = TRUE)

# Check if all files are present
if (length(loo_files) == length(models) && length(diag_files) == length(models)) {
  loo_list <- lapply(loo_files, readRDS)
  names(loo_list) <- str_remove(str_remove(basename(loo_files), "loo_pooled_"), ".rds")
  
  diag_list <- lapply(diag_files, read.csv)
  
  # Perform comparison and save
  comp <- loo_compare(loo_list)
  write.csv(as.data.frame(comp), here("Results", "pooled_model_comparison_results.csv"), row.names = FALSE)
  write.csv(do.call(rbind, diag_list), here("Results", "pooled_model_diagnostics.csv"), row.names = FALSE)
  
  message("Pooled model comparison complete.")
} else {
  warning("Could not find all necessary loo/diagnostic files for comparison. Please check Fits/ and Results/ directories.")
}
