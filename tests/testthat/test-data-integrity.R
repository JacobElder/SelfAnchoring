# tests/testthat/test-data-integrity.R
# Structural and range checks on key results CSV files.
# These tests verify that model runs completed successfully and produced
# outputs with expected columns, sample sizes, and parameter value ranges.
# Tests are skipped if files don't yet exist (i.e., models haven't been run).

library(testthat)
library(here)

skip_if_missing <- function(path) {
  if (!file.exists(path)) skip(paste("File not found:", basename(path)))
}

# ── Individual difference outputs ────────────────────────────────────────────
test_that("ind_diffs_s1_full.csv has expected structure", {
  path <- here("Results", "ind_diffs_s1_full.csv")
  skip_if_missing(path)
  df <- read.csv(path)

  expect_equal(nrow(df), 61)
  expected_cols <- c("subj_idx", "m", "bias", "lambda", "w", "subject_mcr", "SING.Ind")
  expect_true(all(expected_cols %in% names(df)),
              info = paste("Missing columns:", paste(setdiff(expected_cols, names(df)), collapse=", ")))

  expect_true(all(df$m      >= 0   & df$m      <= 10),  info = "α out of [0,10]")
  expect_true(all(df$bias   >= 0   & df$bias   <= 1),   info = "γ out of [0,1]")
  expect_true(all(df$lambda >= 0   & df$lambda <= 5),   info = "λ out of [0,5]")
  expect_true(all(df$w      >= 0   & df$w      <= 1),   info = "w out of [0,1]")
  expect_true(all(df$subject_mcr >= 0, na.rm = TRUE),   info = "MCR should be non-negative")
})

test_that("ind_diffs_s2_full.csv has expected structure", {
  path <- here("Results", "ind_diffs_s2_full.csv")
  skip_if_missing(path)
  df <- read.csv(path)

  expect_equal(nrow(df), 181)
  expect_true("outgroup" %in% names(df))
  expect_setequal(unique(df$outgroup), c("Not UCR", "UCLA", "CSU LA"))
  expect_true(all(c("m", "bias", "lambda", "w", "subject_mcr") %in% names(df)))

  expect_true(all(df$m      >= 0 & df$m      <= 10))
  expect_true(all(df$bias   >= 0 & df$bias   <= 1))
  expect_true(all(df$lambda >= 0 & df$lambda <= 5))
  expect_true(all(df$w      >= 0 & df$w      <= 1))
})

test_that("ind_diffs_s3_full.csv has expected structure", {
  path <- here("Results", "ind_diffs_s3_full.csv")
  skip_if_missing(path)
  df <- read.csv(path)

  expect_true(nrow(df) >= 260 & nrow(df) <= 270,
              info = paste("Unexpected S3 N:", nrow(df)))
  expect_true("condition" %in% names(df))
  expect_setequal(unique(df$condition), c("Minority", "Majority"))
  expect_true(all(c("m", "bias", "lambda", "w", "subject_mcr") %in% names(df)))
})

# ── Condition-level parameter summaries ───────────────────────────────────────
test_that("param_by_condition_s2_desc.csv has all conditions and params", {
  path <- here("Results", "param_by_condition_s2_desc.csv")
  skip_if_missing(path)
  df <- read.csv(path)

  expect_true(all(c("condition", "M", "SD", "n", "param") %in% names(df)))
  expect_setequal(unique(df$condition), c("Negation", "High-Status", "Low-Status"))
  expect_true("bias" %in% df$param)
  expect_true("m"    %in% df$param)
  expect_true("lambda" %in% df$param)

  # γ should be dramatically lower in High-Status
  bias_hs  <- df$M[df$condition == "High-Status"  & df$param == "bias"]
  bias_neg <- df$M[df$condition == "Negation"     & df$param == "bias"]
  expect_true(bias_hs < bias_neg - 0.10,
              info = paste("Expected γ_HighStatus < γ_Negation - 0.10; got", bias_hs, "vs", bias_neg))

  # α should be similar across conditions (within 0.5 units on [0,10] scale)
  m_vals <- df$M[df$param == "m"]
  expect_true(diff(range(m_vals)) < 0.5,
              info = paste("α varies too much across conditions:", diff(range(m_vals))))
})

test_that("param_by_condition_s3_desc.csv has both conditions", {
  path <- here("Results", "param_by_condition_s3_desc.csv")
  skip_if_missing(path)
  df <- read.csv(path)

  expect_setequal(unique(df$condition), c("Minority", "Majority"))
  expect_true("lambda" %in% df$param)

  # λ should be nearly identical between conditions (within 0.05)
  lam_vals <- df$M[df$param == "lambda"]
  expect_true(diff(range(lam_vals)) < 0.05,
              info = paste("Unexpected λ difference S3:", diff(range(lam_vals))))
})

# ── Model summary files ───────────────────────────────────────────────────────
test_that("summary_s1_sym_lambda.csv has mu_pr rows with valid values", {
  path <- here("Results", "summary_s1_sym_lambda.csv")
  skip_if_missing(path)
  df <- read.csv(path)

  expect_true(all(c("variable", "median", "q5", "q95", "rhat") %in% names(df)))

  mu_pr_rows <- df[grepl("mu_pr\\[", df$variable), ]
  expect_equal(nrow(mu_pr_rows), 4, info = "Expected 4 mu_pr rows for sym_lambda")

  # All rhat < 1.01 for hyperparameters
  expect_true(all(mu_pr_rows$rhat < 1.01, na.rm = TRUE),
              info = paste("mu_pr rhat > 1.01:", mu_pr_rows$rhat[mu_pr_rows$rhat >= 1.01]))
})

test_that("summary_s2_sym_lambda.csv converged", {
  path <- here("Results", "summary_s2_sym_lambda.csv")
  skip_if_missing(path)
  df <- read.csv(path)
  expect_true(max(df$rhat, na.rm = TRUE) < 1.05,
              info = paste("Max rhat:", max(df$rhat, na.rm = TRUE)))
})

test_that("summary_s3_sym_lambda.csv converged", {
  path <- here("Results", "summary_s3_sym_lambda.csv")
  skip_if_missing(path)
  df <- read.csv(path)
  expect_true(max(df$rhat, na.rm = TRUE) < 1.05,
              info = paste("Max rhat:", max(df$rhat, na.rm = TRUE)))
})

# ── LOO files (must be RDS, not CSV — just check existence and class) ─────────
test_that("LOO RDS files exist and are loo objects for S1", {
  models <- c("bias", "symmetric", "sym_lambda", "asym_lambda")
  for (m in models) {
    path <- here("Fits", paste0("loo_s1_", m, ".rds"))
    skip_if_missing(path)
    obj <- readRDS(path)
    expect_s3_class(obj, "loo", info = paste("S1", m, "LOO is not a loo object"))
    expect_true("pointwise" %in% names(obj),
                info = paste("S1", m, "LOO missing pointwise field"))
  }
})

test_that("sym_lambda beats bias model in S1 LOO", {
  loo_sl <- here("Fits", "loo_s1_sym_lambda.rds")
  loo_b  <- here("Fits", "loo_s1_bias.rds")
  skip_if_missing(loo_sl)
  skip_if_missing(loo_b)
  sl <- readRDS(loo_sl)
  b  <- readRDS(loo_b)
  expect_true(sl$estimates["elpd_loo", "Estimate"] > b$estimates["elpd_loo", "Estimate"],
              info = "sym_lambda should have higher ELPD than bias model in S1")
})

# ── Correlation outputs ───────────────────────────────────────────────────────
test_that("correlations_s1_full.csv has expected structure", {
  path <- here("Results", "correlations_s1_full.csv")
  skip_if_missing(path)
  df <- read.csv(path)

  expect_true(all(c("param", "scale", "r", "p_raw", "p_fdr") %in% names(df)))
  expect_true(all(abs(df$r) <= 1 + 1e-9), info = "r values must be in [-1, 1]")
  expect_true(all(df$p_raw >= 0 & df$p_raw <= 1), info = "p values must be in [0, 1]")
  expect_true(all(df$p_fdr >= 0 & df$p_fdr <= 1), info = "FDR p values must be in [0, 1]")
  expect_true(all(df$p_fdr >= df$p_raw - 1e-9),
              info = "FDR-corrected p should be >= raw p")

  # MCR × SING.Ind should be the strongest S1 correlate
  mcr_sing <- df[df$param == "subject_mcr" & df$scale == "SING.Ind", ]
  expect_true(nrow(mcr_sing) > 0, info = "Missing MCR × SING.Ind row")
  expect_true(mcr_sing$r > 0.40, info = paste("MCR × SING.Ind r too low:", mcr_sing$r))
})

# ── Pooled behavioral data ────────────────────────────────────────────────────
test_that("marginal_effects_pooled_ames.csv has expected conditions", {
  path <- here("Results", "marginal_effects_pooled_ames.csv")
  skip_if_missing(path)
  df <- read.csv(path)

  expect_true(all(c("estimate", "std.error", "p.value", "conf.low", "conf.high") %in% names(df)))
  expect_true(all(df$std.error > 0, na.rm = TRUE))
  # All confidence intervals should straddle estimate
  expect_true(all(df$conf.low <= df$estimate + 1e-9, na.rm = TRUE))
  expect_true(all(df$conf.high >= df$estimate - 1e-9, na.rm = TRUE))
})

# ── Parameter recovery outputs ────────────────────────────────────────────────
test_that("parameter_recovery_symlambda_results.csv has valid recovery correlations", {
  path <- here("Results", "parameter_recovery_symlambda_results.csv")
  skip_if_missing(path)
  df <- read.csv(path)

  expect_true(all(c("model", "param", "subj_idx", "true", "recovered") %in% names(df)))
  expect_setequal(unique(df$param), c("m", "bias", "lambda", "w"))

  # Recovery correlations should all be positive and reasonably high
  for (p in unique(df$param)) {
    d <- df[df$param == p, ]
    r <- cor(d$true, d$recovered)
    expect_true(r > 0.30, info = paste("Recovery r for", p, "is low:", round(r, 3)))
  }
})

test_that("parameter_recovery_asymlambda_results.csv has valid recovery", {
  path <- here("Results", "parameter_recovery_asymlambda_results.csv")
  skip_if_missing(path)
  df <- read.csv(path)

  expect_setequal(unique(df$param), c("m_in", "m_out", "bias", "lambda", "w"))
  for (p in unique(df$param)) {
    d <- df[df$param == p, ]
    r <- cor(d$true, d$recovered)
    expect_true(r > 0.20, info = paste("Recovery r for", p, ":", round(r, 3)))
  }
})
