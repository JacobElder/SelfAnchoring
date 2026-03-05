# Model Comparison Script for Study 3 (Racial Group Size)
# Optimized for Speed using cmdstanr and parallel model fitting.

library(cmdstanr)
library(tidyverse)
library(loo)
library(here)
library(igraph)
library(parallel)

# Setup
options(mc.cores = 2) 

# 1. Load Data
fulldf <- read.csv(here("Study 3/Cleaning/output/fullTest_fixed.csv")) %>% filter(!is.na(ingChoiceN))
traindf <- read.csv(here("Study 3/Cleaning/output/fullTrain_fixed.csv")) %>% filter(!is.na(selfResp))

# Ensure subjects exist in BOTH training and test sets
common_ids <- intersect(unique(fulldf$subID), unique(traindf$subID))
fulldf <- fulldf %>% filter(subID %in% common_ids)
traindf <- traindf %>% filter(subID %in% common_ids)

# Load Network for Similarity Matrix
posDf <- read.csv(here("Combined/input/adjacencyMatrix_p.csv"))
posMat <- as.matrix(posDf)
posGraph <- graph_from_adjacency_matrix(posMat, mode = "max")
simMat <- similarity(posGraph, method = "dice")

uIds <- sort(common_ids)
maxSubjs <- length(uIds)
maxTrials <- max(fulldf$trialTotalT2)
maxTrain <- 91

# Prepare Stan Data List
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
  t_count <- nrow(s_df); tr_count <- nrow(s_train)
  stan_data$nTrials[i] <- t_count; stan_data$nTrain[i] <- tr_count
  if (t_count > 0) {
    stan_data$groupChoice[i, 1:t_count] <- s_df$ingChoiceN + 1
    if (tr_count > 0) stan_data$prevSim[i, 1:t_count, 1:tr_count] <- simMat[s_df$Idx, s_train$Idx]
  }
  if (tr_count > 0) stan_data$prevSelf[i, 1:tr_count] <- s_train$selfResp
}

# 2. Define Models
models <- list(
  bias = "S_Bias.stan",
  symmetric = "S_Symmetric.stan",
  sym_lambda = "S_Sym_Lambda.stan",
  asym_lambda = "S_Asym_Lambda.stan"
)

# 3. Fit in Parallel
message("Fitting Study 3 models in parallel...")
results_list <- mclapply(names(models), function(model_name) {
  message(paste("Starting model:", model_name))
  mod <- cmdstan_model(here("Computational Models", models[[model_name]]))
  fit <- mod$sample(
    data = stan_data, seed = 789, chains = 4, parallel_chains = 4,
    iter_warmup = 2000, iter_sampling = 2000, adapt_delta = 0.95, max_treedepth = 15, init = 0, refresh = 100
  )
  fit$save_object(here("Fits", paste0("fit_s3_", model_name, ".rds")))
  l <- fit$loo(); saveRDS(l, here("Fits", paste0("loo_s3_", model_name, ".rds")))
  sum_fit <- fit$summary(); write.csv(sum_fit, here("Results", paste0("summary_s3_", model_name, ".csv")))
  diag <- fit$diagnostic_summary()
  diag_df <- data.frame(
    model = model_name, max_rhat = max(sum_fit$rhat, na.rm = TRUE),
    num_divergent = sum(diag$num_divergent),
    converged = (max(sum_fit$rhat, na.rm = TRUE) < 1.01 & sum(diag$num_divergent) == 0)
  )
  return(list(loo = l, diag = diag_df))
}, mc.cores = 4)

# 4. Save Results
loo_list <- lapply(results_list, function(x) x$loo)
diag_list <- lapply(results_list, function(x) x$diag)
comp <- loo_compare(loo_list)
write.csv(as.data.frame(comp), here("Results", "model_comparison_s3_results.csv"))
write.csv(do.call(rbind, diag_list), here("Results", "model_diagnostics_s3.csv"))
message("Study 3 complete.")
