# Parameter Recovery (Empirical) — NoW Architecture
# Uses actual participant posterior medians as true parameters.
# Uses actual self-ratings (prevSelf_real) for both simulation and Stan input.
# This tests whether the specific parameter values driving the paper's results
# are recoverable — a stronger claim than random-distribution recovery.
#
# Two-process split to avoid macOS OOM-kill:
#   Process 1 (setup): load data + empirical params, simulate choices,
#                       save pr_empirical_cache.rds, exit.
#   Process 2 (fit):   load cache, fit Stan models, report recovery stats.
#
# Usage:
#   caffeinate -i sh -c 'Rscript "Parameter Recovery/run_parameter_recovery_empirical.R" setup && \
#                         Rscript "Parameter Recovery/run_parameter_recovery_empirical.R"'

library(here)

inv_logit <- function(x) 1 / (1 + exp(-x))

# ── Parse command-line arg ───────────────────────────────────────────────────
args       <- commandArgs(trailingOnly = TRUE)
run_models <- if (length(args) > 0) args[1] else "both"

# ── PROCESS 1: setup ─────────────────────────────────────────────────────────
if (run_models == "setup") {

  # 1. Load Study 1 data structure (identical to random recovery setup)
  message("Loading Study 1 data structure...")
  fulldf  <- read.csv(here("Study 1/Cleaning/output/fullTest.csv"))
  fulldf  <- fulldf[!is.na(fulldf$ingChoiceN), ]
  traindf <- read.csv(here("Study 1/Cleaning/output/fullTrain.csv"))
  traindf <- traindf[!is.na(traindf$selfResp), ]
  common_ids <- sort(intersect(unique(fulldf$subID), unique(traindf$subID)))
  fulldf  <- fulldf[fulldf$subID  %in% common_ids, ]
  traindf <- traindf[traindf$subID %in% common_ids, ]
  simMat  <- readRDS(here("Parameter Recovery/simMat_cached.rds"))   # 148×148
  uIds      <- sort(common_ids)
  nSubjects <- length(uIds)
  maxTrials <- max(fulldf$trialTotal)
  maxTrain  <- 91
  prevSim_real  <- array(0, c(nSubjects, maxTrials, maxTrain))
  prevSelf_real <- array(0, c(nSubjects, maxTrain))
  nTrials_vec   <- integer(nSubjects)
  nTrain_vec    <- integer(nSubjects)
  for (i in seq_along(uIds)) {
    s_df    <- fulldf[fulldf$subID   == uIds[i], ]
    s_train <- traindf[traindf$subID == uIds[i], ]
    t_count  <- nrow(s_df)
    tr_count <- nrow(s_train)
    nTrials_vec[i] <- t_count
    nTrain_vec[i]  <- tr_count
    if (t_count > 0 && tr_count > 0)
      prevSim_real[i, 1:t_count, 1:tr_count] <- simMat[s_df$Idx, s_train$Idx]
    if (tr_count > 0)
      prevSelf_real[i, 1:tr_count] <- s_train$selfResp
  }
  message(sprintf("Data loaded: %d subjects, maxTrials=%d, maxTrain=%d",
                  nSubjects, maxTrials, maxTrain))

  # 2. Load empirical true params — Sym_Lambda_NoW
  # ind_diffs_s1_full.csv: individual posterior medians on natural scale
  message("Loading empirical true params...")
  ind_diffs <- read.csv(here("Results/ind_diffs_s1_full.csv"))
  ind_diffs <- ind_diffs[order(ind_diffs$subj_idx), ]
  true_sym  <- data.frame(
    subj_idx    = ind_diffs$subj_idx,
    m           = ind_diffs$m,
    bias        = ind_diffs$bias,
    lambda      = ind_diffs$lambda,
    subject_mcr = ind_diffs$subject_mcr
  )
  stopifnot(nrow(true_sym) == nSubjects)

  # 3. Load empirical true params — Asym_Lambda_NoW
  # params_ind_s1_asym_lambda.csv: individual posterior medians, natural scale
  asym_raw <- read.csv(here("Results/params_ind_s1_asym_lambda.csv"))
  extract_param <- function(df, prefix) {
    rows <- df[grepl(paste0("^", prefix, "\\["), df$variable), ]
    rows$idx <- as.integer(regmatches(rows$variable, regexpr("[0-9]+", rows$variable)))
    rows[order(rows$idx), "median"]
  }
  true_asym <- data.frame(
    subj_idx = 1:nSubjects,
    m_in     = extract_param(asym_raw, "m_in"),
    m_out    = extract_param(asym_raw, "m_out"),
    bias     = extract_param(asym_raw, "bias"),
    lambda   = extract_param(asym_raw, "lambda")
  )
  stopifnot(nrow(true_asym) == nSubjects)

  # Compute true subject_mcr for asym from actual params + actual self-ratings.
  # subject_mcr = mean_t(simW_in_t / simW_out_t)
  compute_subject_mcr_asym <- function(params) {
    mcr_vec <- numeric(nSubjects)
    for (i in 1:nSubjects) {
      tp   <- params[i, ]
      nt   <- nTrials_vec[i]; ntr <- nTrain_vec[i]
      self_i   <- prevSelf_real[i, 1:ntr]
      GPin     <- inv_logit( tp$m_in  * (self_i - 4))
      GPout    <- inv_logit(-tp$m_out * (self_i - 4))
      PS_lam   <- prevSim_real[i, 1:nt, 1:ntr] ^ tp$lambda
      simW_in  <- as.vector(PS_lam %*% GPin)  + 1e-9
      simW_out <- as.vector(PS_lam %*% GPout) + 1e-9
      mcr_vec[i] <- mean(simW_in / simW_out)
    }
    mcr_vec
  }
  true_asym$subject_mcr <- compute_subject_mcr_asym(true_asym)

  cat("\nSym_Lambda_NoW empirical true param ranges:\n")
  for (p in c("m", "bias", "lambda", "subject_mcr"))
    cat(sprintf("  %s: [%.3f, %.3f] mean=%.3f\n",
                p, min(true_sym[[p]]), max(true_sym[[p]]), mean(true_sym[[p]])))

  cat("\nAsym_Lambda_NoW empirical true param ranges:\n")
  for (p in c("m_in", "m_out", "bias", "lambda", "subject_mcr"))
    cat(sprintf("  %s: [%.3f, %.3f] mean=%.3f\n",
                p, min(true_asym[[p]]), max(true_asym[[p]]), mean(true_asym[[p]])))

  # 4. Simulate choices using prevSelf_real (not GMRF synthetic ratings)
  simulate_choices_sym_now <- function(true_params) {
    choices <- array(0L, c(nSubjects, maxTrials))
    for (i in 1:nSubjects) {
      tp  <- true_params[i, ]
      nt  <- nTrials_vec[i]; ntr <- nTrain_vec[i]
      GP   <- inv_logit(tp$m * (prevSelf_real[i, 1:ntr] - 4))
      PS   <- prevSim_real[i, 1:nt, 1:ntr] ^ tp$lambda
      simW_in  <- as.vector(PS %*% GP) + 1e-9
      simW_out <- as.vector(PS %*% (1 - GP)) + 1e-9
      logit_p  <- log(tp$bias + 1e-9) - log(1 - tp$bias + 1e-9) +
                  log(simW_in) - log(simW_out)
      choices[i, 1:nt] <- rbinom(nt, 1, inv_logit(logit_p)) + 1L
    }
    choices
  }

  simulate_choices_asym_now <- function(true_params) {
    choices <- array(0L, c(nSubjects, maxTrials))
    for (i in 1:nSubjects) {
      tp  <- true_params[i, ]
      nt  <- nTrials_vec[i]; ntr <- nTrain_vec[i]
      self_i <- prevSelf_real[i, 1:ntr]
      GPin    <- inv_logit( tp$m_in  * (self_i - 4))
      GPout   <- inv_logit(-tp$m_out * (self_i - 4))
      PS      <- prevSim_real[i, 1:nt, 1:ntr] ^ tp$lambda
      simW_in  <- as.vector(PS %*% GPin)  + 1e-9
      simW_out <- as.vector(PS %*% GPout) + 1e-9
      logit_p  <- log(tp$bias + 1e-9) - log(1 - tp$bias + 1e-9) +
                  log(simW_in) - log(simW_out)
      choices[i, 1:nt] <- rbinom(nt, 1, inv_logit(logit_p)) + 1L
    }
    choices
  }

  set.seed(42)
  sim_sym  <- simulate_choices_sym_now(true_sym)
  set.seed(43)
  sim_asym <- simulate_choices_asym_now(true_asym)

  # 5. Save cache and exit
  cache_path <- here("Parameter Recovery/pr_empirical_cache.rds")
  saveRDS(list(
    nSubjects     = nSubjects,
    maxTrials     = maxTrials,
    maxTrain      = maxTrain,
    nTrials_vec   = nTrials_vec,
    nTrain_vec    = nTrain_vec,
    prevSim_real  = prevSim_real,
    prevSelf_real = prevSelf_real,
    true_sym      = true_sym,
    sim_sym       = sim_sym,
    true_asym     = true_asym,
    sim_asym      = sim_asym
  ), cache_path)
  message(sprintf("Cache saved to %s", cache_path))
  message("Setup complete. Run fit process next.")
  quit(save = "no", status = 0)
}

# ── PROCESS 2: fit ────────────────────────────────────────────────────────────
# run_models values:
#   "both"       — Sym_Lambda_NoW + Asym_Lambda_NoW (default)
#   "sym_lambda" — Sym_Lambda_NoW only
#   "asym_lambda"— Asym_Lambda_NoW only
run_sym  <- run_models %in% c("both", "sym_lambda")
run_asym <- run_models %in% c("both", "asym_lambda")

cache_path <- here("Parameter Recovery/pr_empirical_cache.rds")
if (!file.exists(cache_path))
  stop("pr_empirical_cache.rds not found. Run setup first: Rscript ... setup")
message("Loading pr_empirical_cache.rds...")
cache <- readRDS(cache_path)
list2env(cache, envir = environment())
rm(cache); gc()
message(sprintf("Cache loaded: %d subjects, maxTrials=%d, maxTrain=%d",
                nSubjects, maxTrials, maxTrain))

suppressPackageStartupMessages(library(cmdstanr))

# Stan data uses prevSelf_real (actual participant ratings, not GMRF)
make_stan_data_emp <- function(sim_choices) {
  list(
    nSubjects = nSubjects, maxTrials = maxTrials, maxTrain = maxTrain,
    nTrain = nTrain_vec, nTrials = nTrials_vec, groupChoice = sim_choices,
    prevSim = prevSim_real, prevSelf = prevSelf_real
  )
}

get_ind_param <- function(sum_pr, param_name) {
  rows <- sum_pr[grepl(paste0("^", param_name, "\\["), sum_pr$variable), ]
  rows$subj_idx <- as.integer(regmatches(rows$variable, regexpr("[0-9]+", rows$variable)))
  rows[order(rows$subj_idx), c("subj_idx", "median")]
}

run_recovery_emp <- function(model_name, true_params, param_names, sim_choices, seed_fit) {
  message(sprintf("\n--- Empirical Parameter Recovery: %s ---", model_name))

  stan_data_pr <- make_stan_data_emp(sim_choices)
  mod <- cmdstan_model(here("Computational Models", paste0("S_", model_name, ".stan")))

  fit_pr <- mod$sample(
    data = stan_data_pr, seed = seed_fit, chains = 4, parallel_chains = 4,
    iter_warmup = 1000, iter_sampling = 1000, adapt_delta = 0.99,
    max_treedepth = 12, init = 0, refresh = 100
  )

  sum_pr  <- fit_pr$summary(variables = param_names)
  diag_pr <- fit_pr$diagnostic_summary()

  max_rhat  <- max(sum_pr$rhat, na.rm = TRUE)
  n_div     <- sum(diag_pr$num_divergent)
  min_ess   <- min(sum_pr$ess_bulk, na.rm = TRUE)
  converged <- max_rhat < 1.01 && n_div == 0
  cat(sprintf("Convergence: max Rhat=%.3f, divergent transitions=%d, min ESS_bulk=%.0f\n",
              max_rhat, n_div, min_ess))
  if (!converged)
    warning(sprintf("%s empirical recovery fit has convergence issues.", model_name))

  recovery_df <- do.call(rbind, lapply(param_names, function(p) {
    rec <- get_ind_param(sum_pr, p)
    data.frame(
      model = model_name, param = p, subj_idx = true_params$subj_idx,
      true = true_params[[p]], recovered = rec$median
    )
  }))

  cat(sprintf("\n=== Empirical Recovery Summary: %s ===\n", model_name))
  cat(sprintf("%-14s %6s %8s %8s\n", "Param", "r", "RMSE", "Bias"))
  recovery_stats <- do.call(rbind, lapply(param_names, function(p) {
    d <- recovery_df[recovery_df$param == p, ]
    r    <- cor(d$true, d$recovered)
    rmse <- sqrt(mean((d$true - d$recovered)^2))
    bias <- mean(d$recovered - d$true)
    cat(sprintf("%-14s %6.3f %8.3f %8.3f\n", p, r, rmse, bias))
    data.frame(model = model_name, param = p, r = r, rmse = rmse, bias = bias)
  }))

  slug <- tolower(gsub("_", "", model_name))

  write.csv(recovery_df,
            here("Results", paste0("parameter_recovery_empirical_", slug, "_results.csv")),
            row.names = FALSE)
  write.csv(recovery_stats,
            here("Results", paste0("parameter_recovery_empirical_", slug, "_summary.csv")),
            row.names = FALSE)

  conv_df <- data.frame(
    model         = model_name,
    max_rhat      = max_rhat,
    n_divergent   = n_div,
    pct_divergent = round(100 * n_div / (4 * 1000), 3),
    min_ess_bulk  = min_ess,
    converged     = converged
  )
  write.csv(conv_df,
            here("Results", paste0("parameter_recovery_empirical_", slug, "_convergence.csv")),
            row.names = FALSE)

  message(sprintf(
    "Results saved to Results/parameter_recovery_empirical_%s_{results,summary,convergence}.csv",
    slug))

  rm(fit_pr, sum_pr, diag_pr, recovery_stats, conv_df); gc()
  invisible(recovery_df)
}

# ── Run recoveries ────────────────────────────────────────────────────────────
# Sym_Lambda_NoW: m, bias, lambda
# (S_Sym_Lambda_NoW.stan does not output subject_mcr; it is a symmetric model)
if (run_sym)
  run_recovery_emp("Sym_Lambda_NoW", true_sym,
                   c("m", "bias", "lambda"),
                   sim_sym, seed_fit = 996)

# Asym_Lambda_NoW: m_in, m_out, bias, lambda, subject_mcr
# subject_mcr computed from true asym params + actual self-ratings in setup
if (run_asym)
  run_recovery_emp("Asym_Lambda_NoW", true_asym,
                   c("m_in", "m_out", "bias", "lambda", "subject_mcr"),
                   sim_asym, seed_fit = 995)

message("\nEmpirical parameter recovery complete.")
