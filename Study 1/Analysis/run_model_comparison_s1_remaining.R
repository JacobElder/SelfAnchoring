# AD HOC: Study 1 — sym_lambda only (re-run to capture subject_mcr in GQ)
# Previous run used old struct_vars pattern that excluded subject_mcr.
# This re-runs sym_lambda with corrected struct_vars, then rebuilds full comparison
# using previously saved LOO files for bias, symmetric, asym_lambda.
# Delete this script after winning-model MCR is confirmed in params_ind_s1_sym_lambda.csv.

library(cmdstanr)
library(tidyverse)
library(loo)
library(here)
library(igraph)
library(parallel)

options(mc.cores = parallel::detectCores())

# 1. Load Data
fulldf  <- read.csv(here("Study 1/Cleaning/output/fullTest.csv")) %>% filter(!is.na(ingChoiceN))
traindf <- read.csv(here("Study 1/Cleaning/output/fullTrain.csv")) %>% filter(!is.na(selfResp))

common_ids <- intersect(unique(fulldf$subID), unique(traindf$subID))
fulldf  <- fulldf  %>% filter(subID %in% common_ids)
traindf <- traindf %>% filter(subID %in% common_ids)

posDf    <- read.csv(here("Pooled/input/adjacencyMatrix_p.csv"))
posGraph <- graph_from_adjacency_matrix(as.matrix(posDf), mode = "max")
simMat   <- similarity(posGraph, method = "dice")

uIds      <- sort(common_ids)
maxSubjs  <- length(uIds)
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

# 2. Run sym_lambda only
message("Starting Study 1 model: sym_lambda (re-run for subject_mcr)")

mod <- cmdstan_model(here("Computational Models", "S_Sym_Lambda.stan"))

fit <- mod$sample(
  data            = stan_data,
  seed            = 123,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1000,
  iter_sampling   = 2000,
  adapt_delta     = 0.99,
  max_treedepth   = 12,
  init            = 0,
  refresh         = 100
)

# A. Save LOO (overwrite)
l <- fit$loo()
saveRDS(l, here("Fits", "loo_s1_sym_lambda.rds"))

# B. Summary — corrected struct_vars keeps subject_mcr[i] (N-length), excludes trial-level arrays
struct_vars <- grep("^(log_lik|p_pred\\[|mcr\\[)",
                    fit$metadata()$stan_variables, value = TRUE, invert = TRUE)
sum_fit <- fit$summary(variables = struct_vars)
write.csv(sum_fit, here("Results", "summary_s1_sym_lambda.csv"))

# C. Individual parameters (includes subject_mcr[i])
ind_params <- sum_fit %>%
  filter(str_detect(variable, "\\[")) %>%
  select(variable, median, rhat, ess_bulk)
write.csv(ind_params, here("Results", "params_ind_s1_sym_lambda.csv"))

# D. Diagnostics
diag <- fit$diagnostic_summary()
pk   <- l$diagnostics$pareto_k
diag_df <- data.frame(
  model         = "sym_lambda",
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

# E. Per-subject Pareto k (overwrite)
subj_pk <- compute_subj_pareto_k(l, "S1", "sym_lambda", uIds, stan_data$nTrials, maxTrials)
write.csv(subj_pk, here("Results", "pareto_k_subj_s1_sym_lambda.csv"), row.names = FALSE)

# 3. Load previously saved LOO for remaining models, then run full comparison
loo_bias      <- readRDS(here("Fits", "loo_s1_bias.rds"))
loo_symmetric <- readRDS(here("Fits", "loo_s1_symmetric.rds"))
loo_asym      <- readRDS(here("Fits", "loo_s1_asym_lambda.rds"))

loo_list <- list(
  bias        = loo_bias,
  symmetric   = loo_symmetric,
  sym_lambda  = l,
  asym_lambda = loo_asym
)

comp <- loo_compare(loo_list)
write.csv(as.data.frame(comp), here("Results", "model_comparison_s1_results.csv"))
write.csv(diag_df, here("Results", "model_diagnostics_s1_sym_lambda.csv"), row.names = FALSE)

message("Study 1 sym_lambda re-run complete. subject_mcr now in params_ind_s1_sym_lambda.csv")
