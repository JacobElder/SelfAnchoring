# Pooled Model Comparison Analysis
library(cmdstanr)
library(tidyverse)
library(loo)
library(here)
library(parallel)

# 1. Setup
options(mc.cores = 4) 
stan_data <- readRDS(here("Fits", "pooled_stan_data.rds"))

# 2. Define Models
models <- list(
  bias = "S_Pooled_Bias.stan",
  symmetric = "S_Pooled_Symmetric.stan",
  sym_lambda = "S_Pooled_Sym_Lambda.stan",
  asym_lambda = "S_Pooled_Asym_Lambda.stan"
)

# 3. Fit function
fit_pooled <- function(model_name) {
  message(paste("Fitting Pooled model:", model_name))
  mod <- cmdstan_model(here("Computational Models", models[[model_name]]))
  
  fit <- mod$sample(
    data = stan_data, seed = 1234, chains = 4, parallel_chains = 4,
    iter_warmup = 2000, iter_sampling = 2000, adapt_delta = 0.95, max_treedepth = 15, init = 1.0, refresh = 100
  )
  
  # Save Fit and results
  fit$save_object(here("Fits", paste0("fit_pooled_", model_name, ".rds")))
  l <- fit$loo(); saveRDS(l, here("Fits", paste0("loo_pooled_", model_name, ".rds")))
  sum_fit <- fit$summary(); write.csv(sum_fit, here("Results", paste0("summary_pooled_", model_name, ".csv")))
  
  return(l)
}

# 4. Run Analysis
loo_results <- list()
for(m in names(models)) {
  loo_results[[m]] <- fit_pooled(m)
}

# 5. Model Comparison
comp <- loo_compare(loo_results)
write.csv(as.data.frame(comp), here("Results", "pooled_model_comparison_results.csv"))
message("Pooled model comparison complete.")
