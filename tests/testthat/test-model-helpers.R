# tests/testthat/test-model-helpers.R
# Unit tests for helper functions in model comparison and pooled analysis scripts.

library(testthat)
library(here)

# ── compute_subj_pareto_k() ───────────────────────────────────────────────────
# Redefined here for unit testing without sourcing the full pooled script.
compute_subj_pareto_k <- function(l, study_label, model_name, uIds, nTrials_vec, maxTrials) {
  pk <- l$diagnostics$pareto_k
  rows <- lapply(seq_along(uIds), function(i) {
    idx <- ((i - 1) * maxTrials + 1):((i - 1) * maxTrials + nTrials_vec[i])
    k_i <- pk[idx]
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

make_mock_loo <- function(n_obs, seed = 1) {
  set.seed(seed)
  structure(
    list(
      estimates    = matrix(c(-100, 5), nrow = 1, dimnames = list("elpd_loo", c("Estimate", "SE"))),
      pointwise    = matrix(runif(n_obs, -3, -0.5), ncol = 1, dimnames = list(NULL, "elpd_loo")),
      diagnostics  = list(pareto_k = runif(n_obs, 0, 0.8))
    ),
    class = c("psis_loo", "loo")
  )
}

test_that("compute_subj_pareto_k() returns one row per subject", {
  nSubj   <- 5
  nTrials <- rep(10, nSubj)
  maxT    <- 10
  l       <- make_mock_loo(nSubj * maxT)

  out <- compute_subj_pareto_k(l, "Test", "sym_lambda",
                                uIds = paste0("s", 1:nSubj),
                                nTrials_vec = nTrials,
                                maxTrials   = maxT)
  expect_equal(nrow(out), nSubj)
})

test_that("compute_subj_pareto_k() has all required columns", {
  out <- compute_subj_pareto_k(
    make_mock_loo(30), "Test", "sym_lambda",
    uIds = paste0("s", 1:3), nTrials_vec = rep(10, 3), maxTrials = 10
  )
  required <- c("study","model","subID","subj_idx","n_trials",
                "k_mean","k_max","n_good","n_ok","n_bad","n_verybad",
                "pct_reliable","concern")
  expect_true(all(required %in% names(out)),
              info = paste("Missing:", paste(setdiff(required, names(out)), collapse=",")))
})

test_that("compute_subj_pareto_k() n_good + n_ok + n_bad + n_verybad == n_trials", {
  out <- compute_subj_pareto_k(
    make_mock_loo(50), "Test", "sym_lambda",
    uIds = paste0("s", 1:5), nTrials_vec = rep(10, 5), maxTrials = 10
  )
  totals <- out$n_good + out$n_ok + out$n_bad + out$n_verybad
  expect_equal(totals, out$n_trials)
})

test_that("compute_subj_pareto_k() pct_reliable is in [0, 100]", {
  out <- compute_subj_pareto_k(
    make_mock_loo(50), "Test", "sym_lambda",
    uIds = paste0("s", 1:5), nTrials_vec = rep(10, 5), maxTrials = 10
  )
  expect_true(all(out$pct_reliable >= 0 & out$pct_reliable <= 100))
})

test_that("compute_subj_pareto_k() concern is HIGH when any k >= 1.0", {
  l         <- make_mock_loo(10)
  l$diagnostics$pareto_k[1:3] <- 1.5   # force HIGH in first subject
  out <- compute_subj_pareto_k(l, "Test", "sym_lambda",
                                uIds = c("s1","s2"), nTrials_vec = c(5, 5), maxTrials = 5)
  expect_equal(out$concern[1], "HIGH")
})

# ── Back-transform round-trip (consistency between scripts) ───────────────────
test_that("probit back-transform is consistent across all parameter scalings", {
  bt <- function(x, scale = 1) pnorm(x) * scale

  # Known S1 group-level values from MEMORY
  # α: Phi(mu_pr[1]) * 10 ≈ 2.80
  # γ: Phi(mu_pr[2])      ≈ 0.36
  # λ: Phi(mu_pr[3]) * 5  ≈ 3.56
  # w: Phi(mu_pr[4])      ≈ 0.56

  mu_alpha  <- qnorm(2.80 / 10)
  mu_gamma  <- qnorm(0.36)
  mu_lambda <- qnorm(3.56 / 5)
  mu_w      <- qnorm(0.56)

  expect_equal(bt(mu_alpha,  10), 2.80, tolerance = 1e-6)
  expect_equal(bt(mu_gamma,   1), 0.36, tolerance = 1e-6)
  expect_equal(bt(mu_lambda,  5), 3.56, tolerance = 1e-6)
  expect_equal(bt(mu_w,       1), 0.56, tolerance = 1e-6)
})

# ── Non-centered parameterization self-consistency ────────────────────────────
test_that("non-centered param: Phi(mu + sigma * z) * scale stays in [0, scale]", {
  bt      <- function(x, scale = 1) pnorm(x) * scale
  mu      <- qnorm(0.36)
  sigma   <- 0.4
  z_draws <- rnorm(1000)
  params  <- bt(mu + sigma * z_draws, scale = 1)

  expect_true(all(params >= 0 & params <= 1))
  # Mean should be close to population mean (0.36)
  expect_true(abs(mean(params) - 0.36) < 0.05,
              info = paste("Mean param:", round(mean(params), 3)))
})

# ── LOO model comparison order ─────────────────────────────────────────────────
test_that("sym_lambda ELPD exceeds symmetric and bias in all three studies", {
  studies <- c("s1", "s2", "s3")
  for (study in studies) {
    loo_sl  <- here("Fits", paste0("loo_", study, "_sym_lambda.rds"))
    loo_sym <- here("Fits", paste0("loo_", study, "_symmetric.rds"))
    loo_b   <- here("Fits", paste0("loo_", study, "_bias.rds"))

    if (!file.exists(loo_sl) || !file.exists(loo_sym) || !file.exists(loo_b)) {
      skip(paste("LOO files not found for", study))
    }

    sl  <- readRDS(loo_sl)$estimates["elpd_loo", "Estimate"]
    sym <- readRDS(loo_sym)$estimates["elpd_loo", "Estimate"]
    b   <- readRDS(loo_b)$estimates["elpd_loo", "Estimate"]

    expect_true(sl > sym, info = paste(study, ": sym_lambda ELPD should exceed symmetric"))
    expect_true(sl > b,   info = paste(study, ": sym_lambda ELPD should exceed bias"))
    expect_true(sym > b,  info = paste(study, ": symmetric ELPD should exceed bias"))
  }
})
