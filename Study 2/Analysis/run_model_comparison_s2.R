# Model Comparison Script for Study 2 (University Status)
# Optimized for Speed (cmdstanr) and Portability (CSV Export)
# Memory-Efficient Version: Saves each model's results to disk and reloads at the end.

library(cmdstanr)
library(tidyverse)
library(loo)
library(here)
library(igraph)
library(parallel)

# Setup
options(mc.cores = parallel::detectCores()) 

# 1. Load Data (unchanged)
fulldf <- read.csv(here("Study 2/Cleaning/output/fullTest.csv")) %>% filter(!is.na(ingChoiceN))
traindf <- read.csv(here("Study 2/Cleaning/output/fullTrain_fixed.csv")) %>% filter(!is.na(selfResp))
common_ids <- intersect(unique(fulldf$subID), unique(traindf$subID))
fulldf <- fulldf %>% filter(subID %in% common_ids)
traindf <- traindf %>% filter(subID %in% common_ids)
posDf <- read.csv(here("Pooled/input/adjacencyMatrix_p.csv"))
posMat <- as.matrix(posDf)
posGraph <- graph_from_adjacency_matrix(posMat, mode = "max")
simMat <- similarity(posGraph, method = "dice")
uIds <- sort(common_ids)
maxSubjs <- length(uIds)
maxTrials <- max(fulldf$trialTotalT2)
maxTrain <- 91

stan_data <- list(
  nSubjects = maxSubjs,
  maxTrials = maxTrials,
  maxTrain = maxTrain,
  nTrain = integer(maxSubjs),
  nTrials = integer(maxSubjs),
  groupChoice = array(0, c(maxSubjs, maxTrials)),
  prevSim = array(0, c(maxSubjs, maxTrials, maxTrain)),
  prevSelf = array(0, c(maxSubjs, maxTrain))
)
for(i in 1:length(uIds)) {
  s_df <- filter(fulldf, subID == uIds[i])
  s_train <- filter(traindf, subID == uIds[i])
  t_count <- nrow(s_df)
  tr_count <- nrow(s_train)
  stan_data$nTrials[i] <- t_count
  stan_data$nTrain[i] <- tr_count
  if (t_count > 0) {
    stan_data$groupChoice[i, 1:t_count] <- s_df$ingChoiceN + 1
    if (tr_count > 0) {
      stan_data$prevSim[i, 1:t_count, 1:tr_count] <- simMat[s_df$Idx, s_train$Idx]
    }
  }
  if (tr_count > 0) {
    stan_data$prevSelf[i, 1:tr_count] <- s_train$selfResp
  }
}

# Helper: per-subject Pareto k summary (unchanged)
compute_subj_pareto_k <- function(l, study_label, model_name, uIds, nTrials_vec, maxTrials) {
  pk <- l$diagnostics$pareto_k
  rows <- lapply(seq_along(uIds), function(i) {
    idx <- ((i - 1) * maxTrials + 1):((i - 1) * maxTrials + nTrials_vec[i])
    k_i <- pk[idx]
    data.frame(
      study        = study_label, model = model_name, subID = uIds[i], subj_idx = i,
      n_trials = nTrials_vec[i], k_mean = round(mean(k_i), 4), k_max = round(max(k_i), 4),
      n_good = sum(k_i < 0.5), n_ok = sum(k_i >= 0.5 & k_i < 0.7), n_bad = sum(k_i >= 0.7 & k_i < 1.0),
      n_verybad = sum(k_i >= 1.0), pct_reliable = round(100 * mean(k_i < 0.7), 1),
      concern = ifelse(any(k_i >= 1.0), "HIGH", ifelse(any(k_i >= 0.7), "MODERATE", ifelse(any(k_i >= 0.5), "LOW", "NONE")))
    )
  })
  do.call(rbind, rows)
}

# 2. Define Models
models <- list(
  bias = "S_Bias.stan",
  symmetric = "S_Symmetric.stan",
  sym_lambda = "S_Sym_Lambda.stan",
  asym_lambda = "S_Asym_Lambda.stan"
)

# 3. Fit function (memory-efficient)
fit_and_save <- function(model_name) {
  message(paste("Starting Study 2 model:", model_name))
  mod <- cmdstan_model(here("Computational Models", models[[model_name]]))
  
  fit <- mod$sample(
    data = stan_data, seed = 456, chains = 4, parallel_chains = 4,
    iter_warmup = 1000, iter_sampling = 2000, adapt_delta = 0.99,
    max_treedepth = 12, init = 0, refresh = 100
  )
  
  l <- fit$loo()
  saveRDS(l, here("Fits", paste0("loo_s2_", model_name, ".rds")))
  
  all_vars    <- fit$metadata()$stan_variables
  struct_vars <- all_vars[!all_vars %in% c("log_lik", "p_pred", "mcr")]
  sum_fit <- fit$summary(variables = struct_vars)
  write.csv(sum_fit, here("Results", paste0("summary_s2_", model_name, ".csv")), row.names = FALSE)
  
  ind_params <- sum_fit %>%
    filter(str_detect(variable, "\\[")) %>%
    select(variable, median, rhat, ess_bulk)
  write.csv(ind_params, here("Results", paste0("params_ind_s2_", model_name, ".csv")), row.names = FALSE)
  
  diag <- fit$diagnostic_summary()
  pk   <- l$diagnostics$pareto_k
  diag_df <- data.frame(
    model = model_name, max_rhat = max(sum_fit$rhat, na.rm = TRUE),
    num_divergent = sum(diag$num_divergent), pct_divergent = round(100 * sum(diag$num_divergent) / (4 * 2000), 3),
    min_ess_bulk = min(sum_fit$ess_bulk, na.rm = TRUE), converged = (max(sum_fit$rhat, na.rm = TRUE) < 1.01 & sum(diag$num_divergent) == 0),
    pk_good = sum(pk < 0.5), pk_ok = sum(pk >= 0.5 & pk < 0.7), pk_bad = sum(pk >= 0.7 & pk < 1.0),
    pk_verybad = sum(pk >= 1.0), pk_pct_ok = round(100 * mean(pk < 0.7), 1),
    pk_max = round(max(pk), 3), loo_reliable = mean(pk < 0.7) > 0.9
  )
  write.csv(diag_df, here("Results", paste0("diag_s2_", model_name, ".csv")), row.names = FALSE)

  subj_pk <- compute_subj_pareto_k(l, "S2", model_name, uIds, stan_data$nTrials, maxTrials)
  write.csv(subj_pk, here("Results", paste0("pareto_k_subj_s2_", model_name, ".csv")), row.names = FALSE)
  
  rm(fit, l, sum_fit, ind_params, diag, pk, diag_df, subj_pk); gc()
  return(NULL)
}

# 4. Fit Models
for(m in names(models)) {
  fit_and_save(m)
}

# 5. Save Comparison Results (Load from disk)
message("All models fit for Study 2. Comparing results from saved files.")
loo_files <- list.files(path = here("Fits"), pattern = "loo_s2_.*\\.rds", full.names = TRUE)
diag_files <- list.files(path = here("Results"), pattern = "diag_s2_.*\\.csv", full.names = TRUE)

if (length(loo_files) >= length(models) && length(diag_files) >= length(models)) {
  model_basenames <- paste0(names(models), ".rds")
  loo_files <- loo_files[grepl(paste(model_basenames, collapse="|"), basename(loo_files))]

  diag_model_basenames <- paste0(names(models), ".csv")
  diag_files <- diag_files[grepl(paste(diag_model_basenames, collapse="|"), basename(diag_files))]

  loo_list <- lapply(loo_files, readRDS)
  names(loo_list) <- str_remove(str_remove(basename(loo_files), "loo_s2_"), ".rds")
  
  diag_list <- lapply(diag_files, read.csv)
  
  comp <- loo_compare(loo_list)
  write.csv(as.data.frame(comp), here("Results", "model_comparison_s2_results.csv"), row.names = FALSE)
  write.csv(do.call(rbind, diag_list), here("Results", "model_diagnostics_s2.csv"), row.names = FALSE)
  
  message("Study 2 comparison complete.")
} else {
  warning("Could not find all necessary loo/diagnostic files for S2 comparison.")
}
