# Helper script to produce Study 2 CSVs from finished fits
library(loo)
library(here)
library(tidyverse)

message("Generating Study 2 CSVs from existing fits in Fits/...")

# List models used in Study 2
models <- c("bias", "symmetric", "sym_lambda", "asym_lambda")

loo_results <- list()

for (m in models) {
  loo_file <- here("Fits", paste0("loo_s2_", m, ".rds"))
  if (file.exists(loo_file)) {
    message(paste("Loading LOO for:", m))
    loo_results[[m]] <- readRDS(loo_file)
  } else {
    message(paste("Warning: LOO file not found for", m))
  }
}

if (length(loo_results) > 0) {
  message("Comparing models...")
  comp <- loo_compare(loo_results)
  write.csv(as.data.frame(comp), here("Results", "model_comparison_s2_results.csv"))
  message("Saved Results/model_comparison_s2_results.csv")
} else {
  message("Error: No LOO results found to compare.")
}

message("Done.")
