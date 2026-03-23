# Parameter Recovery — Current Model Architecture
# Two-process split to avoid macOS OOM-kill:
#   Process 1 (setup): Rscript ... setup
#     Base R + here only (~20 MB). Loads data, generates GMRF, draws params,
#     simulates choices, saves pr_cache.rds, exits (OS reclaims all pages).
#   Process 2 (fit):   Rscript ... [both|sym_lambda|asym_lambda]
#     Fresh R process loads cmdstanr + pr_cache.rds (~170 MB peak) and fits.
#
# Usage:
#   caffeinate -i sh -c 'Rscript "Parameter Recovery/run_parameter_recovery.R" setup && \
#                         Rscript "Parameter Recovery/run_parameter_recovery.R"'

library(here)

set.seed(42)

# ── Parse command-line arg ───────────────────────────────────────────────────
args       <- commandArgs(trailingOnly = TRUE)
run_models <- if (length(args) > 0) args[1] else "both"

inv_logit <- function(x) 1 / (1 + exp(-x))

# ── PROCESS 1: setup ─────────────────────────────────────────────────────────
if (run_models == "setup") {

  # 1. Load Study 1 data structure
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
  message(sprintf("Data loaded: %d subjects, maxTrials=%d, maxTrain=%d", nSubjects, maxTrials, maxTrain))

  # 2. GMRF self-ratings
  # Synthetic self-ratings: centered at 4 (avoids GP ceiling from self-enhancement
  # bias), smooth over trait network via GMRF precision = alpha*L + eps*I.
  # Each subject gets an independent GMRF draw; SD standardised to 1.13 (real data).
  message("Generating GMRF self-ratings...")
  n_traits   <- nrow(simMat)   # 148
  L_mat      <- readRDS(here("Parameter Recovery/L_mat_cached.rds"))
  gmrf_alpha <- 2.0   # higher = more network-aligned ratings
  gmrf_eps   <- 0.01  # regularises singular Laplacian
  Sigma_gmrf <- solve(gmrf_alpha * L_mat + gmrf_eps * diag(n_traits))

  set.seed(77)
  gmrf_chol  <- chol(Sigma_gmrf)
  gmrf_raw   <- matrix(rnorm(nSubjects * n_traits), nrow = nSubjects) %*% gmrf_chol
  gmrf_scaled <- t(apply(gmrf_raw, 1, function(x) {
    pmin(pmax(round((x - mean(x)) / sd(x) * 1.13 + 4), 1), 7)
  }))  # nSubjects × n_traits integer matrix

  prevSelf_gmrf <- array(0, c(nSubjects, maxTrain))
  for (i in seq_along(uIds)) {
    s_train <- traindf[traindf$subID == uIds[i], ]
    prevSelf_gmrf[i, 1:nrow(s_train)] <- gmrf_scaled[i, s_train$Idx]
  }

  gmrf_gp <- inv_logit(2.80 * (prevSelf_gmrf[prevSelf_gmrf != 0] - 4))
  message(sprintf(
    "GMRF GP diagnostics (m=2.80): mean=%.3f, pct>0.7=%.1f%%, pct<0.3=%.1f%%",
    mean(gmrf_gp), 100 * mean(gmrf_gp > 0.7), 100 * mean(gmrf_gp < 0.3)
  ))
  message("GMRF self-ratings generated.")

  # 3. Draw true parameters
  draw_params_sym <- function(n, seed = 1) {
    set.seed(seed)
    mu_pr <- c(m = -0.5744, bias = -0.3594, lambda = 0.5649, w = 0.1781)
    sigma <- c(m = 0.4000,  bias =  0.7091, lambda = 0.3343, w = 1.1808)
    pr    <- matrix(rnorm(n * 4), nrow = n)
    data.frame(
      subj_idx = 1:n,
      m        = pnorm(mu_pr["m"]      + sigma["m"]      * pr[,1]) * 10,
      bias     = pnorm(mu_pr["bias"]   + sigma["bias"]   * pr[,2]),
      lambda   = pnorm(mu_pr["lambda"] + sigma["lambda"] * pr[,3]) * 5,
      w        = pnorm(mu_pr["w"]      + sigma["w"]      * pr[,4])
    )
  }
  draw_params_asym <- function(n, seed = 2) {
    set.seed(seed)
    mu_pr <- c(m_in = -0.1106, m_out = -0.4780,
               bias = -0.3794, lambda =  0.5162, w = 0.2069)
    sigma <- c(m_in =  0.4000, m_out =  0.4000,
               bias =  0.7213, lambda =  0.3768, w = 1.1884)
    pr    <- matrix(rnorm(n * 5), nrow = n)
    data.frame(
      subj_idx = 1:n,
      m_in     = pnorm(mu_pr["m_in"]   + sigma["m_in"]   * pr[,1]) * 10,
      m_out    = pnorm(mu_pr["m_out"]  + sigma["m_out"]  * pr[,2]) * 10,
      bias     = pnorm(mu_pr["bias"]   + sigma["bias"]   * pr[,3]),
      lambda   = pnorm(mu_pr["lambda"] + sigma["lambda"] * pr[,4]) * 5,
      w        = pnorm(mu_pr["w"]      + sigma["w"]      * pr[,5])
    )
  }

  # 4. Simulate choices
  simulate_choices_sym <- function(true_params) {
    choices <- array(0L, c(nSubjects, maxTrials))
    for (i in 1:nSubjects) {
      tp  <- true_params[i, ]
      nt  <- nTrials_vec[i]; ntr <- nTrain_vec[i]
      GP   <- inv_logit(tp$m * (prevSelf_gmrf[i, 1:ntr] - 4))
      PS   <- prevSim_real[i, 1:nt, 1:ntr] ^ tp$lambda
      simW_in  <- as.vector(PS %*% GP) + 1e-9
      simW_out <- as.vector(PS %*% (1 - GP)) + 1e-9
      logit_p  <- log(tp$bias + 1e-9) - log(1 - tp$bias + 1e-9) + log(simW_in) - log(simW_out)
      p_in     <- tp$w * 0.5 + (1 - tp$w) * inv_logit(logit_p)
      choices[i, 1:nt] <- rbinom(nt, 1, p_in) + 1L
    }
    choices
  }
  simulate_choices_asym <- function(true_params) {
    choices <- array(0L, c(nSubjects, maxTrials))
    for (i in 1:nSubjects) {
      tp  <- true_params[i, ]
      nt  <- nTrials_vec[i]; ntr <- nTrain_vec[i]
      self_i <- prevSelf_gmrf[i, 1:ntr]
      GPin    <- inv_logit( tp$m_in  * (self_i - 4))
      GPout   <- inv_logit(-tp$m_out * (self_i - 4))
      PS      <- prevSim_real[i, 1:nt, 1:ntr] ^ tp$lambda
      simW_in  <- as.vector(PS %*% GPin)  + 1e-9
      simW_out <- as.vector(PS %*% GPout) + 1e-9
      logit_p  <- log(tp$bias + 1e-9) - log(1 - tp$bias + 1e-9) + log(simW_in) - log(simW_out)
      p_in     <- tp$w * 0.5 + (1 - tp$w) * inv_logit(logit_p)
      choices[i, 1:nt] <- rbinom(nt, 1, p_in) + 1L
    }
    choices
  }

  true_sym  <- draw_params_sym(nSubjects, seed = 1)
  sim_sym   <- simulate_choices_sym(true_sym)
  true_asym <- draw_params_asym(nSubjects, seed = 2)
  sim_asym  <- simulate_choices_asym(true_asym)

  cat("\nSym_Lambda true param ranges:\n")
  for (p in c("m","bias","lambda","w"))
    cat(sprintf("  %s: [%.2f, %.2f] mean=%.2f\n", p, min(true_sym[[p]]), max(true_sym[[p]]), mean(true_sym[[p]])))
  cat("\nAsym_Lambda true param ranges:\n")
  for (p in c("m_in","m_out","bias","lambda","w"))
    cat(sprintf("  %s: [%.2f, %.2f] mean=%.2f\n", p, min(true_asym[[p]]), max(true_asym[[p]]), mean(true_asym[[p]])))

  # 5. Save cache and exit
  cache_path <- here("Parameter Recovery/pr_cache.rds")
  saveRDS(list(
    nSubjects     = nSubjects,
    maxTrials     = maxTrials,
    maxTrain      = maxTrain,
    nTrials_vec   = nTrials_vec,
    nTrain_vec    = nTrain_vec,
    prevSim_real  = prevSim_real,
    prevSelf_gmrf = prevSelf_gmrf,
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
#   "both"         — Sym_Lambda + Asym_Lambda (with w)
#   "sym_lambda"   — Sym_Lambda only (with w)
#   "asym_lambda"  — Asym_Lambda only (with w)
#   "now"          — Sym_Lambda_NoW + Asym_Lambda_NoW (no lapse rate; recovery test)
run_sym  <- run_models %in% c("both", "sym_lambda")
run_asym <- run_models %in% c("both", "asym_lambda")
run_now  <- run_models == "now"

# Load cache written by setup process
cache_path <- here("Parameter Recovery/pr_cache.rds")
if (!file.exists(cache_path))
  stop("pr_cache.rds not found. Run setup first: Rscript ... setup")
message("Loading pr_cache.rds...")
cache <- readRDS(cache_path)
list2env(cache, envir = environment())
rm(cache); gc()
message(sprintf("Cache loaded: %d subjects, maxTrials=%d, maxTrain=%d", nSubjects, maxTrials, maxTrain))

# Load cmdstanr in the fresh process (genuinely minimal OS footprint before this point)
suppressPackageStartupMessages(library(cmdstanr))

# ── Helper: assemble stan_data ────────────────────────────────────────────────
make_stan_data <- function(sim_choices) {
  list(
    nSubjects = nSubjects, maxTrials = maxTrials, maxTrain = maxTrain,
    nTrain = nTrain_vec, nTrials = nTrials_vec, groupChoice = sim_choices,
    prevSim = prevSim_real, prevSelf = prevSelf_gmrf
  )
}

# ── Helper: extract recovered individual params ───────────────────────────────
get_ind_param <- function(sum_pr, param_name) {
  rows <- sum_pr[grepl(paste0("^", param_name, "\\["), sum_pr$variable), ]
  rows$subj_idx <- as.integer(regmatches(rows$variable, regexpr("[0-9]+", rows$variable)))
  rows[order(rows$subj_idx), c("subj_idx", "median")]
}

# ── Helper: fit Stan model and report recovery ────────────────────────────────
run_recovery <- function(model_name, true_params, param_names, sim_choices, seed_fit) {
  message(sprintf("\n--- Parameter Recovery: %s ---", model_name))

  stan_data_pr <- make_stan_data(sim_choices)
  mod <- cmdstan_model(here("Computational Models", paste0("S_", model_name, ".stan")))

  fit_pr <- mod$sample(
    data = stan_data_pr, seed = seed_fit, chains = 4, parallel_chains = 4,
    iter_warmup = 1000, iter_sampling = 1000, adapt_delta = 0.99,
    max_treedepth = 12, init = 0, refresh = 100
  )

  sum_pr  <- fit_pr$summary(variables = param_names)  # avoid loading log_lik/p_pred/mcr
  diag_pr <- fit_pr$diagnostic_summary()

  max_rhat  <- max(sum_pr$rhat, na.rm = TRUE)
  n_div     <- sum(diag_pr$num_divergent)
  min_ess   <- min(sum_pr$ess_bulk, na.rm = TRUE)
  converged <- max_rhat < 1.01 && n_div == 0
  cat(sprintf("Convergence: max Rhat=%.3f, divergent transitions=%d, min ESS_bulk=%.0f\n",
              max_rhat, n_div, min_ess))
  if (!converged)
    warning(sprintf("%s recovery fit has convergence issues.", model_name))

  recovery_df <- do.call(rbind, lapply(param_names, function(p) {
    rec <- get_ind_param(sum_pr, p)
    data.frame(
      model = model_name, param = p, subj_idx = true_params$subj_idx,
      true = true_params[[p]], recovered = rec$median
    )
  }))

  cat(sprintf("\n=== Recovery Summary: %s ===\n", model_name))
  cat(sprintf("%-8s %6s %8s %8s\n", "Param", "r", "RMSE", "Bias"))
  recovery_stats <- do.call(rbind, lapply(param_names, function(p) {
    d <- recovery_df[recovery_df$param == p, ]
    r    <- cor(d$true, d$recovered)
    rmse <- sqrt(mean((d$true - d$recovered)^2))
    bias <- mean(d$recovered - d$true)
    cat(sprintf("%-8s %6.3f %8.3f %8.3f\n", p, r, rmse, bias))
    data.frame(model = model_name, param = p, r = r, rmse = rmse, bias = bias)
  }))

  slug <- tolower(gsub("_", "", model_name))

  # Save recovery results
  out_file <- here("Results", paste0("parameter_recovery_", slug, "_results.csv"))
  write.csv(recovery_df, out_file, row.names = FALSE)

  # Save recovery summary (r, RMSE, bias per param)
  summary_file <- here("Results", paste0("parameter_recovery_", slug, "_summary.csv"))
  write.csv(recovery_stats, summary_file, row.names = FALSE)

  # Save convergence diagnostics
  conv_df <- data.frame(
    model         = model_name,
    max_rhat      = max_rhat,
    n_divergent   = n_div,
    pct_divergent = round(100 * n_div / (4 * 1000), 3),
    min_ess_bulk  = min_ess,
    converged     = converged
  )
  conv_file <- here("Results", paste0("parameter_recovery_", slug, "_convergence.csv"))
  write.csv(conv_df, conv_file, row.names = FALSE)

  message(sprintf("Results saved to Results/parameter_recovery_%s_{results,summary,convergence}.csv", slug))

  rm(fit_pr, sum_pr, diag_pr, recovery_stats, conv_df); gc()

  invisible(recovery_df)
}

# ── Run recoveries ────────────────────────────────────────────────────────────
if (run_sym)
  run_recovery("Sym_Lambda", true_sym, c("m","bias","lambda","w"), sim_sym, seed_fit = 999)

if (run_asym)
  run_recovery("Asym_Lambda", true_asym, c("m_in","m_out","bias","lambda","w"), sim_asym, seed_fit = 1000)

# ── No-lapse recovery (NoW) ───────────────────────────────────────────────────
# Simulate without lapse term to test whether w is responsible for poor m/lambda
# recovery. Uses same true params as above minus w; same prevSelf and prevSim.
if (run_now) {
  inv_logit <- function(x) 1 / (1 + exp(-x))

  # Draw params (no w)
  draw_params_sym_now <- function(n, seed = 1) {
    set.seed(seed)
    # Calibrated to Study 1 NoW sym_lambda posterior (probit-scale medians)
    mu_pr <- c(m = -1.8230, bias = -0.1208, lambda = 0.6478)
    sigma <- c(m =  0.7323, bias =  0.3433, lambda = 0.0760)
    pr    <- matrix(rnorm(n * 3), nrow = n)
    data.frame(
      subj_idx = 1:n,
      m        = pnorm(mu_pr["m"]      + sigma["m"]      * pr[,1]) * 10,
      bias     = pnorm(mu_pr["bias"]   + sigma["bias"]   * pr[,2]),
      lambda   = pnorm(mu_pr["lambda"] + sigma["lambda"] * pr[,3]) * 5
    )
  }
  draw_params_asym_now <- function(n, seed = 2) {
    set.seed(seed)
    # Best guess based on Study 1 NoW sym_lambda posterior; update when asym_lambda NoW fit completes
    mu_pr <- c(m_in = -1.8230, m_out = -1.8230, bias = -0.1208, lambda = 0.6478)
    sigma <- c(m_in =  0.7323, m_out =  0.7323, bias =  0.3433, lambda = 0.0760)
    pr    <- matrix(rnorm(n * 4), nrow = n)
    data.frame(
      subj_idx = 1:n,
      m_in     = pnorm(mu_pr["m_in"]   + sigma["m_in"]   * pr[,1]) * 10,
      m_out    = pnorm(mu_pr["m_out"]  + sigma["m_out"]  * pr[,2]) * 10,
      bias     = pnorm(mu_pr["bias"]   + sigma["bias"]   * pr[,3]),
      lambda   = pnorm(mu_pr["lambda"] + sigma["lambda"] * pr[,4]) * 5
    )
  }

  # Simulate choices (no lapse — p_in = inv_logit(logit_p) directly)
  simulate_choices_sym_now <- function(true_params) {
    choices <- array(0L, c(nSubjects, maxTrials))
    for (i in 1:nSubjects) {
      tp  <- true_params[i, ]
      nt  <- nTrials_vec[i]; ntr <- nTrain_vec[i]
      GP   <- inv_logit(tp$m * (prevSelf_gmrf[i, 1:ntr] - 4))
      PS   <- prevSim_real[i, 1:nt, 1:ntr] ^ tp$lambda
      simW_in  <- as.vector(PS %*% GP) + 1e-9
      simW_out <- as.vector(PS %*% (1 - GP)) + 1e-9
      logit_p  <- log(tp$bias + 1e-9) - log(1 - tp$bias + 1e-9) + log(simW_in) - log(simW_out)
      choices[i, 1:nt] <- rbinom(nt, 1, inv_logit(logit_p)) + 1L
    }
    choices
  }
  simulate_choices_asym_now <- function(true_params) {
    choices <- array(0L, c(nSubjects, maxTrials))
    for (i in 1:nSubjects) {
      tp  <- true_params[i, ]
      nt  <- nTrials_vec[i]; ntr <- nTrain_vec[i]
      self_i <- prevSelf_gmrf[i, 1:ntr]
      GPin    <- inv_logit( tp$m_in  * (self_i - 4))
      GPout   <- inv_logit(-tp$m_out * (self_i - 4))
      PS      <- prevSim_real[i, 1:nt, 1:ntr] ^ tp$lambda
      simW_in  <- as.vector(PS %*% GPin)  + 1e-9
      simW_out <- as.vector(PS %*% GPout) + 1e-9
      logit_p  <- log(tp$bias + 1e-9) - log(1 - tp$bias + 1e-9) + log(simW_in) - log(simW_out)
      choices[i, 1:nt] <- rbinom(nt, 1, inv_logit(logit_p)) + 1L
    }
    choices
  }

  true_sym_now  <- draw_params_sym_now(nSubjects, seed = 1)
  sim_sym_now   <- simulate_choices_sym_now(true_sym_now)
  true_asym_now <- draw_params_asym_now(nSubjects, seed = 2)
  sim_asym_now  <- simulate_choices_asym_now(true_asym_now)

  cat("\nSym_Lambda_NoW true param ranges:\n")
  for (p in c("m","bias","lambda"))
    cat(sprintf("  %s: [%.2f, %.2f] mean=%.2f\n", p, min(true_sym_now[[p]]), max(true_sym_now[[p]]), mean(true_sym_now[[p]])))

  cat("\nAsym_Lambda_NoW true param ranges:\n")
  for (p in c("m_in","m_out","bias","lambda"))
    cat(sprintf("  %s: [%.2f, %.2f] mean=%.2f\n", p, min(true_asym_now[[p]]), max(true_asym_now[[p]]), mean(true_asym_now[[p]])))

  run_recovery("Sym_Lambda_NoW",  true_sym_now,  c("m","bias","lambda"),          sim_sym_now,  seed_fit = 998)
  run_recovery("Asym_Lambda_NoW", true_asym_now, c("m_in","m_out","bias","lambda"), sim_asym_now, seed_fit = 997)
}

message("\nParameter recovery complete.")
