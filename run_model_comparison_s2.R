# Model Comparison Script for Study 2 (University Status)
# Optimized for Speed (cmdstanr) and Portability (CSV Export)

library(cmdstanr)
library(tidyverse)
library(loo)
library(here)
library(igraph)
library(parallel)

# Setup
options(mc.cores = parallel::detectCores()) 

# 1. Load Data
# Using the genuine CSV text files (converted from mislabeled Parquet)
fulldf <- read.csv(here("Study 2/Cleaning/output/fullTest.csv")) %>% filter(!is.na(ingChoiceN))
traindf <- read.csv(here("Study 2/Cleaning/output/fullTrain_fixed.csv")) %>% filter(!is.na(selfResp))

# Ensure subjects exist in BOTH training and test sets
common_ids <- intersect(unique(fulldf$subID), unique(traindf$subID))
fulldf <- fulldf %>% filter(subID %in% common_ids)
traindf <- traindf %>% filter(subID %in% common_ids)

# Load Network for Similarity Matrix
posDf <- read.csv(here("Combined/input/adjacencyMatrix_p.csv"))
posMat <- as.matrix(posDf)
posGraph <- graph_from_adjacency_matrix(posMat, mode = "undirected")
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

# 2. Define Models
models <- list(
  symmetric = "S_Logistic_1mOppose_Bias.stan",
  asymmetric = "S_Logistic_Asym.stan",
  asym_lambda = "S_Logistic_Asym_Lambda_NoW.stan",
  asym_lambda_w = "S_Logistic_Asym_Lambda.stan"
)

# 3. Fit function
fit_and_save <- function(model_name) {
  message(paste("Starting model:", model_name))
  
  mod <- cmdstan_model(here("Computational Models", models[[model_name]]))
  
  fit <- mod$sample(
    data = stan_data,
    seed = 456,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.95
  )
  
  # A. Save LOO results
  l <- fit$loo()
  saveRDS(l, paste0("loo_s2_", model_name, ".rds"))
  
  # B. Save Posterior Medians (Portability)
  sum_fit <- fit$summary()
  write.csv(sum_fit, paste0("summary_s2_", model_name, ".csv"))
  
  # C. Extract and save individual level parameters specifically
  # This makes loading individual medians much easier for correlations
  ind_params <- sum_fit %>%
    filter(str_detect(variable, "\\[")) %>%
    select(variable, median, rhat, ess_bulk)
  write.csv(ind_params, paste0("params_ind_s2_", model_name, ".csv"))
  
  return(l)
}

# 4. Sequential Fit (since user running externally)
loo_results <- list()
for(m in names(models)) {
  loo_results[[m]] <- fit_and_save(m)
}

comp <- loo_compare(loo_results)
write.csv(as.data.frame(comp), "model_comparison_s2_results.csv")
message("Study 2 comparison complete.")
