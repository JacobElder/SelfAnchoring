# Parameter Recovery — Current Model Architecture
# Memory-Efficient Version

library(cmdstanr)
library(tidyverse)
library(here)
library(igraph)

set.seed(42)

# ── Parse command-line arg ───────────────────────────────────────────────────
args       <- commandArgs(trailingOnly = TRUE)
run_models <- if (length(args) > 0) args[1] else "both"
run_sym    <- run_models %in% c("both", "sym_lambda")
run_asym   <- run_models %in% c("both", "asym_lambda")

inv_logit <- function(x) 1 / (1 + exp(-x))

# ── 1. Load Study 1 data structure ───────────────────────────────────────────
message("Loading Study 1 data structure...")
fulldf  <- read.csv(here("Study 1/Cleaning/output/fullTest.csv"))  |> filter(!is.na(ingChoiceN))
traindf <- read.csv(here("Study 1/Cleaning/output/fullTrain.csv")) |> filter(!is.na(selfResp))
common_ids <- sort(intersect(unique(fulldf$subID), unique(traindf$subID)))
fulldf  <- filter(fulldf,  subID %in% common_ids)
traindf <- filter(traindf, subID %in% common_ids)
posDf    <- read.csv(here("Pooled/input/adjacencyMatrix_p.csv"))
posGraph <- graph_from_adjacency_matrix(as.matrix(posDf), mode = "max")
simMat   <- similarity(posGraph, method = "dice")
uIds      <- sort(common_ids)
nSubjects <- length(uIds)
maxTrials <- max(fulldf$trialTotal)
maxTrain  <- 91
prevSim_real  <- array(0, c(nSubjects, maxTrials, maxTrain))
prevSelf_real <- array(0, c(nSubjects, maxTrain))
nTrials_vec   <- integer(nSubjects)
nTrain_vec    <- integer(nSubjects)
for (i in seq_along(uIds)) {
  s_df    <- filter(fulldf,  subID == uIds[i])
  s_train <- filter(traindf, subID == uIds[i])
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

# ── Helper: draw synthetic parameters ────────────────────────────────────────
draw_params_sym <- function(n, seed = 1) {
  set.seed(seed)
  mu_pr <- c(m = qnorm(2.80/10), bias = qnorm(0.36), lambda = qnorm(3.56/5), w = qnorm(0.56))
  sigma <- c(m = 0.5, bias = 0.4, lambda = 0.4, w = 0.4)
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
  mu_pr <- c(m_in = qnorm(3.5/10), m_out = qnorm(2.5/10),
             bias = qnorm(0.36), lambda = qnorm(3.56/5), w = qnorm(0.56))
  sigma <- c(m_in = 0.5, m_out = 0.5, bias = 0.4, lambda = 0.4, w = 0.4)
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

# ── Helper: simulate choices ──────────────────────────────────────────────────
simulate_choices_sym <- function(true_params) {
  choices <- array(0L, c(nSubjects, maxTrials))
  for (i in 1:nSubjects) {
    tp  <- true_params[i, ]
    nt  <- nTrials_vec[i]; ntr <- nTrain_vec[i]
    GP   <- inv_logit(tp$m * (prevSelf_real[i, 1:ntr] - 4))
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
    self_i <- prevSelf_real[i, 1:ntr]
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

# ── Helper: assemble stan_data ────────────────────────────────────────────────
make_stan_data <- function(sim_choices) {
  list(
    nSubjects = nSubjects, maxTrials = maxTrials, maxTrain = maxTrain,
    nTrain = nTrain_vec, nTrials = nTrials_vec, groupChoice = sim_choices,
    prevSim = prevSim_real, prevSelf = prevSelf_real
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

  sum_pr <- fit_pr$summary()
  diag_pr  <- fit_pr$diagnostic_summary()

  max_rhat <- max(sum_pr$rhat, na.rm = TRUE)
  n_div    <- sum(diag_pr$num_divergent)
  cat(sprintf("Convergence: max Rhat=%.3f, divergent transitions=%d\n", max_rhat, n_div))
  if (max_rhat > 1.01 || n_div > 0)
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
  for (p in param_names) {
    d <- filter(recovery_df, param == p)
    r    <- cor(d$true, d$recovered)
    rmse <- sqrt(mean((d$true - d$recovered)^2))
    bias <- mean(d$recovered - d$true)
    cat(sprintf("%-8s %6.3f %8.3f %8.3f\n", p, r, rmse, bias))
  }

  out_file <- here("Results", paste0("parameter_recovery_", tolower(gsub("_", "", model_name)), "_results.csv"))
  write.csv(recovery_df, out_file, row.names = FALSE)
  message(sprintf("Results saved to %s", out_file))

  # Explicitly clear memory
  rm(fit_pr, sum_pr, diag_pr); gc()

  invisible(recovery_df)
}

# ── 2. Run S_Sym_Lambda recovery ──────────────────────────────────────────────
if (run_sym) {
  true_sym <- draw_params_sym(nSubjects, seed = 1)
  cat("\nSym_Lambda true param ranges:\n")
  for (p in c("m","bias","lambda","w"))
    cat(sprintf("  %s: [%.2f, %.2f] mean=%.2f\n", p, min(true_sym[[p]]), max(true_sym[[p]]), mean(true_sym[[p]])))

  sim_sym <- simulate_choices_sym(true_sym)
  run_recovery("Sym_Lambda", true_sym, c("m","bias","lambda","w"), sim_sym, seed_fit = 999)
}

# ── 3. Run S_Asym_Lambda recovery ─────────────────────────────────────────────
if (run_asym) {
  true_asym <- draw_params_asym(nSubjects, seed = 2)
  cat("\nAsym_Lambda true param ranges:\n")
  for (p in c("m_in","m_out","bias","lambda","w"))
    cat(sprintf("  %s: [%.2f, %.2f] mean=%.2f\n", p, min(true_asym[[p]]), max(true_asym[[p]]), mean(true_asym[[p]])))

  sim_asym <- simulate_choices_asym(true_asym)
  run_recovery("Asym_Lambda", true_asym, c("m_in","m_out","bias","lambda","w"), sim_asym, seed_fit = 1000)
}

message("\nParameter recovery complete.")
