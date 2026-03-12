# RESCUE: Study 2 asym_lambda — recover post-processing from chain CSVs
# Chains finished; LOO already saved. Process was killed before summary/params/comparison.
# Reads chain CSVs directly with base R — bypasses as_cmdstan_fit / read_cmdstan_csv.

library(posterior)
library(tidyverse)
library(loo)
library(here)

chain_files <- c(
  "/var/folders/4d/8dmpyclj3x10xhz5nx1ynz780000gp/T/RtmpYzxT13/S_Asym_Lambda-202603111433-1-1f07af.csv",
  "/var/folders/4d/8dmpyclj3x10xhz5nx1ynz780000gp/T/RtmpYzxT13/S_Asym_Lambda-202603111433-2-1f07af.csv",
  "/var/folders/4d/8dmpyclj3x10xhz5nx1ynz780000gp/T/RtmpYzxT13/S_Asym_Lambda-202603111433-3-1f07af.csv",
  "/var/folders/4d/8dmpyclj3x10xhz5nx1ynz780000gp/T/RtmpYzxT13/S_Asym_Lambda-202603111433-4-1f07af.csv"
)

stopifnot(all(file.exists(chain_files)))
message("Reading chain CSVs directly (base R)...")

read_stan_csv_simple <- function(file, chain_id) {
  lines <- readLines(file, warn = FALSE)
  # Keep only non-comment, non-empty lines
  keep <- !startsWith(lines, "#") & nchar(trimws(lines)) > 0
  data_lines <- lines[keep]
  df <- read.csv(text = paste(data_lines, collapse = "\n"), header = TRUE)
  df$.chain <- chain_id
  df$.iteration <- seq_len(nrow(df))
  df
}

chain_dfs <- lapply(seq_along(chain_files), function(i) {
  message(paste("  Reading chain", i, "..."))
  read_stan_csv_simple(chain_files[[i]], i)
})
all_df <- bind_rows(chain_dfs)
message(paste("Total rows:", nrow(all_df), "| Cols:", ncol(all_df)))

# Divergences
num_divergent <- sum(all_df$divergent__)
message(paste("Divergences:", num_divergent))

# Build posterior draws_df — exclude sampler diagnostics and GQ flat vectors
sampler_cols <- c("lp__", "accept_stat__", "stepsize__", "treedepth__",
                  "n_leapfrog__", "divergent__", "energy__")
param_df <- all_df[, !names(all_df) %in% sampler_cols]
param_df  <- param_df[, !grepl("^(log_lik|p_pred|mcr)", names(param_df))]
# posterior needs .chain and .iteration
draws <- posterior::as_draws_df(param_df)
message(paste("Variables in draws:", posterior::nvariables(draws)))

# B. Summary
message("Computing summary (may take a minute)...")
sum_fit <- posterior::summarise_draws(draws,
  mean, median, posterior::sd,
  ~quantile(.x, probs = c(0.025, 0.975)),
  posterior::default_convergence_measures()
)
names(sum_fit)[names(sum_fit) == "2.5%"]  <- "q5"
names(sum_fit)[names(sum_fit) == "97.5%"] <- "q95"
write.csv(sum_fit, here("Results", "summary_s2_asym_lambda.csv"), row.names = FALSE)
message("summary_s2_asym_lambda.csv written.")

# C. Individual parameters
ind_params <- sum_fit %>%
  filter(str_detect(variable, "\\[")) %>%
  select(variable, median, rhat, ess_bulk)
write.csv(ind_params, here("Results", "params_ind_s2_asym_lambda.csv"), row.names = FALSE)
message("params_ind_s2_asym_lambda.csv written.")

# D. Diagnostics — use saved LOO for Pareto k
l  <- readRDS(here("Fits", "loo_s2_asym_lambda.rds"))
pk <- l$diagnostics$pareto_k

diag_df <- data.frame(
  model         = "asym_lambda",
  max_rhat      = max(sum_fit$rhat, na.rm = TRUE),
  num_divergent = num_divergent,
  pct_divergent = round(100 * num_divergent / (4 * 2000), 3),
  min_ess_bulk  = min(sum_fit$ess_bulk, na.rm = TRUE),
  converged     = (max(sum_fit$rhat, na.rm = TRUE) < 1.01 & num_divergent == 0),
  pk_good       = sum(pk < 0.5),
  pk_ok         = sum(pk >= 0.5 & pk < 0.7),
  pk_bad        = sum(pk >= 0.7 & pk < 1.0),
  pk_verybad    = sum(pk >= 1.0),
  pk_pct_ok     = round(100 * mean(pk < 0.7), 1),
  pk_max        = round(max(pk), 3),
  loo_reliable  = mean(pk < 0.7) > 0.9
)
message(paste("Max Rhat:", round(max(sum_fit$rhat, na.rm=TRUE), 4),
              "| Min ESS:", round(min(sum_fit$ess_bulk, na.rm=TRUE), 0)))

# E. Per-subject Pareto k — from saved LOO
fulldf  <- read.csv(here("Study 2/Cleaning/output/fullTest.csv")) %>% filter(!is.na(ingChoiceN))
traindf <- read.csv(here("Study 2/Cleaning/output/fullTrain_fixed.csv")) %>% filter(!is.na(selfResp))
common_ids  <- intersect(unique(fulldf$subID), unique(traindf$subID))
fulldf      <- fulldf %>% filter(subID %in% common_ids)
uIds        <- sort(common_ids)
maxTrials   <- max(fulldf$trialTotalT2)
nTrials_vec <- sapply(uIds, function(id) nrow(dplyr::filter(fulldf, subID == id)))

compute_subj_pareto_k <- function(l, study_label, model_name, uIds, nTrials_vec, maxTrials) {
  pk_vec <- l$diagnostics$pareto_k
  rows <- lapply(seq_along(uIds), function(i) {
    idx <- ((i - 1) * maxTrials + 1):((i - 1) * maxTrials + nTrials_vec[i])
    k_i <- pk_vec[idx]
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

subj_pk <- compute_subj_pareto_k(l, "S2", "asym_lambda", uIds, nTrials_vec, maxTrials)
write.csv(subj_pk, here("Results", "pareto_k_subj_s2_asym_lambda.csv"), row.names = FALSE)
message("pareto_k_subj_s2_asym_lambda.csv written.")

# Merge diagnostics
prev_diag <- read.csv(here("Results", "model_diagnostics_s2.csv")) %>%
  dplyr::filter(model != "asym_lambda")
write.csv(bind_rows(prev_diag, diag_df),
          here("Results", "model_diagnostics_s2.csv"), row.names = FALSE)
message("model_diagnostics_s2.csv updated.")

# 4. Full model comparison
loo_bias       <- readRDS(here("Fits", "loo_s2_bias.rds"))
loo_symmetric  <- readRDS(here("Fits", "loo_s2_symmetric.rds"))
loo_sym_lambda <- readRDS(here("Fits", "loo_s2_sym_lambda.rds"))

comp <- loo_compare(list(
  bias        = loo_bias,
  symmetric   = loo_symmetric,
  sym_lambda  = loo_sym_lambda,
  asym_lambda = l
))
write.csv(as.data.frame(comp), here("Results", "model_comparison_s2_results.csv"))
message("model_comparison_s2_results.csv written.")
message("Study 2 rescue complete. No re-fitting required.")
