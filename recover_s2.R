# Recovery script for Study 2
# Fits only the models missing from the Fits/ folder
library(cmdstanr)
library(tidyverse)
library(loo)
library(here)
library(igraph)

# 1. Setup Data (Same as main script)
fulldf <- read.csv(here("Study 2/Cleaning/output/fullTest.csv")) %>% filter(!is.na(ingChoiceN))
traindf <- read.csv(here("Study 2/Cleaning/output/fullTrain_fixed.csv")) %>% filter(!is.na(selfResp))
common_ids <- intersect(unique(fulldf$subID), unique(traindf$subID))
fulldf <- fulldf %>% filter(subID %in% common_ids)
traindf <- traindf %>% filter(subID %in% common_ids)
posDf <- read.csv(here("Combined/input/adjacencyMatrix_p.csv"))
posMat <- as.matrix(posDf)
posGraph <- graph_from_adjacency_matrix(posMat, mode = "max")
simMat <- similarity(posGraph, method = "dice")
uIds <- sort(common_ids)
maxSubjs <- length(uIds)
maxTrials <- max(fulldf$trialTotalT2)
maxTrain <- 91
stan_data <- list(
  nSubjects = maxSubjs, maxTrials = maxTrials, maxTrain = maxTrain,
  nTrain = integer(maxSubjs), nTrials = integer(maxSubjs),
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

# 3. Recovery Loop
loo_results <- list()

for(m in names(models)) {
  loo_file <- here("Fits", paste0("loo_s2_", m, ".rds"))
  
  if (file.exists(loo_file)) {
    message(paste("Skipping", m, "- LOO already exists."))
    loo_results[[m]] <- readRDS(loo_file)
  } else {
    message(paste("Fitting missing model:", m))
    mod <- cmdstan_model(here("Computational Models", models[[m]]))
    fit <- mod$sample(
      data = stan_data, seed = 456, chains = 4, parallel_chains = 4,
      iter_warmup = 2000, iter_sampling = 2000, adapt_delta = 0.95, max_treedepth = 15, init = 0
    )
    # Save Results
    l <- fit$loo(); saveRDS(l, loo_file)
    sum_fit <- fit$summary(); write.csv(sum_fit, here("Results", paste0("summary_s2_", m, ".csv")))
    ind_params <- sum_fit %>% filter(str_detect(variable, "\\[")) %>% select(variable, median, rhat, ess_bulk)
    write.csv(ind_params, here("Results", paste0("params_ind_s2_", m, ".csv")))
    loo_results[[m]] <- l
  }
}

# 4. Final Comparison
comp <- loo_compare(loo_results)
write.csv(as.data.frame(comp), here("Results", "model_comparison_s2_results.csv"))
message("Study 2 recovery and comparison complete.")
