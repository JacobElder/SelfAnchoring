# Single-model pooled fitter — called via command line:
#   Rscript Pooled/run_pooled_one_model.R sym_lambda
#   Rscript Pooled/run_pooled_one_model.R asym_lambda

library(cmdstanr)
library(tidyverse)
library(loo)
library(here)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
model_name <- args[1]

stan_files <- list(
  bias        = "S_Pooled_Bias.stan",
  symmetric   = "S_Pooled_Symmetric.stan",
  sym_lambda  = "S_Pooled_Sym_Lambda.stan",
  asym_lambda = "S_Pooled_Asym_Lambda.stan"
)

if (!model_name %in% names(stan_files))
  stop("model_name must be one of: ", paste(names(stan_files), collapse = ", "))

options(mc.cores = parallel::detectCores())
stan_data      <- readRDS(here("Fits", "pooled_stan_data.rds"))
pooled_uIds    <- if (!is.null(stan_data$subIDs)) stan_data$subIDs else seq_len(stan_data$nSubjects)
pooled_maxTrials <- max(stan_data$nTrials)

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

message(paste("Fitting Pooled model:", model_name))
mod <- cmdstan_model(here("Computational Models", stan_files[[model_name]]))

fit <- mod$sample(
  data            = stan_data,
  seed            = 1234,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 800,
  iter_sampling   = 1000,
  adapt_delta     = 0.95,
  max_treedepth   = 10,
  init            = 0,
  refresh         = 100
)

gc()
l <- fit$loo()
saveRDS(l, here("Fits", paste0("loo_pooled_", model_name, ".rds")))

all_vars    <- fit$metadata()$stan_variables
struct_vars <- all_vars[!all_vars %in% c("log_lik", "p_pred", "mcr")]
sum_fit     <- fit$summary(variables = struct_vars)
write.csv(sum_fit, here("Results", paste0("summary_pooled_", model_name, ".csv")), row.names = FALSE)

ind_params <- sum_fit %>%
  filter(str_detect(variable, "\\[")) %>%
  select(variable, median, rhat, ess_bulk)
write.csv(ind_params, here("Results", paste0("params_ind_pooled_", model_name, ".csv")), row.names = FALSE)

diag <- fit$diagnostic_summary()
pk   <- l$diagnostics$pareto_k
diag_df <- data.frame(
  model           = model_name,
  max_rhat        = max(sum_fit$rhat, na.rm = TRUE),
  num_divergent   = sum(diag$num_divergent),
  pct_divergent   = round(100 * sum(diag$num_divergent) / (4 * 1000), 3),
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

subj_pk <- compute_subj_pareto_k(l, "Pooled", model_name, pooled_uIds, stan_data$nTrials, pooled_maxTrials)
write.csv(subj_pk, here("Results", paste0("pareto_k_subj_pooled_", model_name, ".csv")), row.names = FALSE)

message(paste("Done:", model_name))
