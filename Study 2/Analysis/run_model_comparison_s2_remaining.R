# AD HOC: Study 2 — sym_lambda and asym_lambda only
# bias and symmetric LOO already saved (Mar 7-8). Delete this script when done.

library(cmdstanr)
library(tidyverse)
library(loo)
library(here)
library(igraph)
library(parallel)

options(mc.cores = parallel::detectCores())

# 1. Load Data
fulldf  <- read.csv(here("Study 2/Cleaning/output/fullTest.csv")) %>% filter(!is.na(ingChoiceN))
traindf <- read.csv(here("Study 2/Cleaning/output/fullTrain_fixed.csv")) %>% filter(!is.na(selfResp))

common_ids <- intersect(unique(fulldf$subID), unique(traindf$subID))
fulldf  <- fulldf  %>% filter(subID %in% common_ids)
traindf <- traindf %>% filter(subID %in% common_ids)

posDf    <- read.csv(here("Pooled/input/adjacencyMatrix_p.csv"))
posGraph <- graph_from_adjacency_matrix(as.matrix(posDf), mode = "max")
simMat   <- similarity(posGraph, method = "dice")

uIds     <- sort(common_ids)
maxSubjs <- length(uIds)
maxTrials <- max(fulldf$trialTotalT2)
maxTrain  <- 91

stan_data <- list(
  nSubjects = maxSubjs,
  maxTrials = maxTrials,
  maxTrain  = maxTrain,
  nTrain    = integer(maxSubjs),
  nTrials   = integer(maxSubjs),
  groupChoice = array(0, c(maxSubjs, maxTrials)),
  prevSim     = array(0, c(maxSubjs, maxTrials, maxTrain)),
  prevSelf    = array(0, c(maxSubjs, maxTrain))
)

for (i in 1:length(uIds)) {
  s_df    <- filter(fulldf,  subID == uIds[i])
  s_train <- filter(traindf, subID == uIds[i])
  t_count  <- nrow(s_df);  tr_count <- nrow(s_train)
  stan_data$nTrials[i] <- t_count
  stan_data$nTrain[i]  <- tr_count
  if (t_count > 0) {
    stan_data$groupChoice[i, 1:t_count] <- s_df$ingChoiceN + 1
    if (tr_count > 0)
      stan_data$prevSim[i, 1:t_count, 1:tr_count] <- simMat[s_df$Idx, s_train$Idx]
  }
  if (tr_count > 0)
    stan_data$prevSelf[i, 1:tr_count] <- s_train$selfResp
}

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

# 2. Remaining models only
models <- list(
  sym_lambda  = "S_Sym_Lambda.stan",
  asym_lambda = "S_Asym_Lambda.stan"
)

fit_and_save <- function(model_name) {
  message(paste("Starting Study 2 model:", model_name))

  mod <- cmdstan_model(here("Computational Models", models[[model_name]]))

  fit <- mod$sample(
    data            = stan_data,
    seed            = 456,
    chains          = 4,
    parallel_chains = 4,
    iter_warmup     = 1000,
    iter_sampling   = 2000,
    adapt_delta     = 0.99,
    max_treedepth   = 12,
    init            = 0,
    refresh         = 100
  )

  # A. Save LOO
  l <- fit$loo()
  saveRDS(l, here("Fits", paste0("loo_s2_", model_name, ".rds")))

  # B. Summary — exclude trial-level GQ arrays (OOM), keep subject_mcr (small)
  struct_vars <- grep("^(log_lik|p_pred\\[|mcr\\[)",
                      fit$metadata()$stan_variables, value = TRUE, invert = TRUE)
  sum_fit <- fit$summary(variables = struct_vars)
  write.csv(sum_fit, here("Results", paste0("summary_s2_", model_name, ".csv")))

  # C. Individual parameters
  ind_params <- sum_fit %>%
    filter(str_detect(variable, "\\[")) %>%
    select(variable, median, rhat, ess_bulk)
  write.csv(ind_params, here("Results", paste0("params_ind_s2_", model_name, ".csv")))

  # D. Diagnostics
  diag <- fit$diagnostic_summary()
  pk   <- l$diagnostics$pareto_k
  diag_df <- data.frame(
    model         = model_name,
    max_rhat      = max(sum_fit$rhat, na.rm = TRUE),
    num_divergent = sum(diag$num_divergent),
    pct_divergent = round(100 * sum(diag$num_divergent) / (4 * 2000), 3),
    min_ess_bulk  = min(sum_fit$ess_bulk, na.rm = TRUE),
    converged     = (max(sum_fit$rhat, na.rm = TRUE) < 1.01 & sum(diag$num_divergent) == 0),
    pk_good       = sum(pk < 0.5),
    pk_ok         = sum(pk >= 0.5 & pk < 0.7),
    pk_bad        = sum(pk >= 0.7 & pk < 1.0),
    pk_verybad    = sum(pk >= 1.0),
    pk_pct_ok     = round(100 * mean(pk < 0.7), 1),
    pk_max        = round(max(pk), 3),
    loo_reliable  = mean(pk < 0.7) > 0.9
  )

  # E. Per-subject Pareto k
  subj_pk <- compute_subj_pareto_k(l, "S2", model_name, uIds, stan_data$nTrials, maxTrials)
  write.csv(subj_pk, here("Results", paste0("pareto_k_subj_s2_", model_name, ".csv")), row.names = FALSE)

  return(list(loo = l, diag = diag_df))
}

# 3. Fit remaining models
results_list <- list()
for (m in names(models)) {
  results_list[[m]] <- fit_and_save(m)
}

# 4. Load previously saved LOO for bias and symmetric, then run full comparison
loo_bias      <- readRDS(here("Fits", "loo_s2_bias.rds"))
loo_symmetric <- readRDS(here("Fits", "loo_s2_symmetric.rds"))

loo_list <- list(
  bias       = loo_bias,
  symmetric  = loo_symmetric,
  sym_lambda  = results_list$sym_lambda$loo,
  asym_lambda = results_list$asym_lambda$loo
)

comp <- loo_compare(loo_list)
write.csv(as.data.frame(comp), here("Results", "model_comparison_s2_results.csv"))

diag_list <- lapply(results_list, function(x) x$diag)
write.csv(do.call(rbind, diag_list), here("Results", "model_diagnostics_s2.csv"))

message("Study 2 comparison complete.")
