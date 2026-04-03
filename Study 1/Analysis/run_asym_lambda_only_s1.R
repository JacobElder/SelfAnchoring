# Refit ONLY asym_lambda for Study 1 to get subject_mcr
library(cmdstanr)
library(tidyverse)
library(loo)
library(here)
library(igraph)

# Setup
options(mc.cores = parallel::detectCores())

# 1. Load Data
fulldf <- read.csv(here("Study 1/Cleaning/output/fullTest.csv")) %>% filter(!is.na(ingChoiceN))
traindf <- read.csv(here("Study 1/Cleaning/output/fullTrain.csv")) %>% filter(!is.na(selfResp))
common_ids <- intersect(unique(fulldf$subID), unique(traindf$subID))
fulldf <- fulldf %>% filter(subID %in% common_ids)
traindf <- traindf %>% filter(subID %in% common_ids)

# Sim Matrix
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

# 2. Fit Asym_Lambda_NoW
model_name <- "asym_lambda"
message("Refitting Study 1 model: ", model_name)
mod <- cmdstan_model(here("Computational Models", "S_Asym_Lambda_NoW.stan"))

fit <- mod$sample(
  data = stan_data, seed = 123, chains = 4, parallel_chains = 4,
  iter_warmup = 1000, iter_sampling = 2000, adapt_delta = 0.99,
  max_treedepth = 12, init = 0, refresh = 100
)

# 3. Save Results
l <- fit$loo()
saveRDS(l, here("Fits", paste0("loo_s1_", model_name, ".rds")))

all_vars    <- fit$metadata()$stan_variables
struct_vars <- all_vars[!all_vars %in% c("log_lik", "p_pred", "mcr")]
sum_fit <- fit$summary(variables = struct_vars)
write.csv(sum_fit, here("Results", paste0("summary_s1_", model_name, ".csv")), row.names = FALSE)

ind_params <- sum_fit %>%
  filter(str_detect(variable, "\[")) %>%
  select(variable, median, rhat, ess_bulk)
write.csv(ind_params, here("Results", paste0("params_ind_s1_", model_name, ".csv")), row.names = FALSE)

message("Study 1 refit complete for asym_lambda.")
