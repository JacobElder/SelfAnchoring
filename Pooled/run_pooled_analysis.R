# Fit Pooled Model
library(cmdstanr)
library(tidyverse)

# 1. Load Data
library(here)
stan_data <- readRDS(here("Fits", "pooled_stan_data.rds"))

# 2. Compile Model
mod <- cmdstan_model(here("Computational Models", "S_Pooled_Asym_Lambda.stan"))

# 3. Fit
fit <- mod$sample(
  data = stan_data,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 2000,
  adapt_delta = 0.95,
  max_treedepth = 15,
  init = 0,
  refresh = 100
)

# 4. Save
fit$save_object(here("Fits", "fit_pooled_asym_lambda.rds"))
summary_pooled <- fit$summary()
write.csv(summary_pooled, here("Results", "summary_pooled_asym_lambda.csv"))
message("Pooled model fit complete.")
