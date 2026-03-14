# Ad-Hoc Exploratory Script for Study 2
# Runs the Symmetric + Lambda + Certainty model on Study 2 data.

library(cmdstanr)
library(tidyverse)
library(loo)
library(here)
library(parallel)

# Setup
options(mc.cores = parallel::detectCores()) 

# 1. Load Data (copied from main S2 script)
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
  nSubjects = maxSubjs, maxTrials = maxTrials, maxTrain = maxTrain,
  nTrain = integer(maxSubjs), nTrials = integer(maxSubjs),
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

# 2. Fit the Certainty Model
model_name <- "S_Symmetric_Lambda_Certainty"
message(paste("Starting Study 2 exploratory model:", model_name))

mod <- cmdstan_model(here("Computational Models", paste0(model_name, ".stan")))

fit <- mod$sample(
  data = stan_data, seed = 456, chains = 4, parallel_chains = 4,
  iter_warmup = 1000, iter_sampling = 2000, adapt_delta = 0.99,
  max_treedepth = 12, init = 0, refresh = 100
)

# 3. Save key outputs
message("Model fit complete. Saving results...")
fit$save_object(here("Fits", paste0("fit_s2_certainty_exploratory.rds")))

l <- fit$loo()
saveRDS(l, here("Fits", paste0("loo_s2_certainty_exploratory.rds")))

sum_fit <- fit$summary()
write.csv(sum_fit, here("Results", paste0("summary_s2_certainty_exploratory.csv")), row.names = FALSE)

message("Exploratory run for Study 2 complete.")
