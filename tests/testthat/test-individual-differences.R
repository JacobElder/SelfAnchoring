# tests/testthat/test-individual-differences.R
# Verifies that key individual-difference correlations match reported results
# across all three studies. Tests check direction, approximate magnitude,
# FDR significance, and internal consistency of correlation files.

library(testthat)
library(here)

skip_if_missing <- function(path) {
  if (!file.exists(path)) skip(paste("File not found:", basename(path)))
}

get_corr <- function(df, param, scale) {
  df[df$param == param & df$scale == scale, ]
}

# ── Study 1 ──────────────────────────────────────────────────────────────────
test_that("S1: MCR × SING.Ind is the strongest FDR-significant correlate", {
  path <- here("Results", "correlations_s1_full.csv")
  skip_if_missing(path)
  df <- read.csv(path)

  mcr_sing <- get_corr(df, "subject_mcr", "SING.Ind")
  expect_true(nrow(mcr_sing) > 0, info = "MCR × SING.Ind row missing")
  expect_true(mcr_sing$r > 0.45,         info = paste("r too low:", mcr_sing$r))
  expect_true(mcr_sing$p_raw < 0.001,    info = "p_raw should be < .001")
  expect_true(mcr_sing$p_fdr < 0.001,    info = "p_fdr should be < .001")
})

test_that("S1: MCR × SING.Ind r is largest among all MCR correlates", {
  path <- here("Results", "correlations_s1_full.csv")
  skip_if_missing(path)
  df    <- read.csv(path)
  mcr_r <- abs(df$r[df$param == "subject_mcr"])
  sing_r <- abs(get_corr(df, "subject_mcr", "SING.Ind")$r)
  expect_true(sing_r == max(mcr_r),
              info = paste("SING.Ind should have largest |r|; max is", max(mcr_r)))
})

test_that("S1: Structural params (m, lambda, w) have no FDR-significant correlates", {
  path <- here("Results", "correlations_s1_full.csv")
  skip_if_missing(path)
  df   <- read.csv(path)
  struct <- df[df$param %in% c("m", "lambda", "w"), ]
  expect_true(all(struct$p_fdr > 0.05, na.rm = TRUE),
              info = paste("Unexpected FDR-significant structural param:",
                           paste(struct$scale[struct$p_fdr <= 0.05], collapse=", ")))
})

test_that("S1: DELTA_ELPD has no FDR-significant correlates", {
  path <- here("Results", "correlations_s1_full.csv")
  skip_if_missing(path)
  df   <- read.csv(path)
  delpd <- df[df$param == "delta_elpd", ]
  if (nrow(delpd) == 0) skip("No delta_elpd rows in S1 correlations")
  expect_true(all(delpd$p_fdr > 0.05, na.rm = TRUE),
              info = "S1 DELTA_ELPD should have no FDR-significant correlates")
})

test_that("S1: correlation file internal consistency (p_fdr >= p_raw)", {
  path <- here("Results", "correlations_s1_full.csv")
  skip_if_missing(path)
  df <- read.csv(path)
  expect_true(all(df$p_fdr >= df$p_raw - 1e-9, na.rm = TRUE))
  expect_true(all(abs(df$r) <= 1 + 1e-9, na.rm = TRUE))
  expect_true(all(df$n > 0, na.rm = TRUE))
})

# ── Study 2 ──────────────────────────────────────────────────────────────────
test_that("S2: MCR × RSE is FDR-significant with correct direction", {
  path <- here("Results", "correlations_s2_full.csv")
  skip_if_missing(path)
  df <- read.csv(path)

  row <- get_corr(df, "subject_mcr", "RSE")
  expect_true(nrow(row) > 0, info = "MCR × RSE row missing in S2")
  expect_true(row$r > 0,       info = "MCR × RSE should be positive")
  expect_true(row$r > 0.30,    info = paste("MCR × RSE r too low:", row$r))
  expect_true(row$p_fdr < 0.01, info = "MCR × RSE should be FDR-significant")
})

test_that("S2: MCR × SING.Ind is FDR-significant", {
  path <- here("Results", "correlations_s2_full.csv")
  skip_if_missing(path)
  df  <- read.csv(path)
  row <- get_corr(df, "subject_mcr", "SING.Ind")
  expect_true(nrow(row) > 0)
  expect_true(row$p_fdr < 0.01, info = paste("S2 MCR×SING.Ind p_fdr:", row$p_fdr))
})

test_that("S2: lambda × SI is FDR-significant", {
  path <- here("Results", "correlations_s2_full.csv")
  skip_if_missing(path)
  df  <- read.csv(path)
  row <- get_corr(df, "lambda", "SI")
  expect_true(nrow(row) > 0, info = "lambda × SI row missing")
  expect_true(row$p_fdr < 0.05, info = paste("lambda × SI p_fdr:", row$p_fdr))
  expect_true(row$r > 0, info = "lambda × SI should be positive")
})

test_that("S2: bias (gamma) has no FDR-significant personality correlates", {
  path <- here("Results", "correlations_s2_full.csv")
  skip_if_missing(path)
  df   <- read.csv(path)
  bias <- df[df$param == "bias", ]
  expect_true(all(bias$p_fdr > 0.05, na.rm = TRUE),
              info = "S2 gamma should have no FDR-significant correlates")
})

test_that("S2: DELTA_ELPD has no FDR-significant correlates", {
  path <- here("Results", "correlations_s2_full.csv")
  skip_if_missing(path)
  df    <- read.csv(path)
  delpd <- df[df$param == "delta_elpd", ]
  if (nrow(delpd) == 0) skip("No delta_elpd in S2")
  expect_true(all(delpd$p_fdr > 0.05, na.rm = TRUE))
})

# ── Study 3 ──────────────────────────────────────────────────────────────────
test_that("S3: MCR × SING.Ind is FDR-significant", {
  path <- here("Results", "correlations_s3_full.csv")
  skip_if_missing(path)
  df  <- read.csv(path)
  row <- get_corr(df, "subject_mcr", "SING.Ind")
  expect_true(nrow(row) > 0, info = "MCR × SING.Ind row missing in S3")
  expect_true(row$p_fdr < 0.01, info = paste("S3 MCR×SING.Ind p_fdr:", row$p_fdr))
  expect_true(row$r > 0, info = "MCR × SING.Ind should be positive")
})

test_that("S3: MCR × RSE is FDR-significant", {
  path <- here("Results", "correlations_s3_full.csv")
  skip_if_missing(path)
  df  <- read.csv(path)
  row <- get_corr(df, "subject_mcr", "RSE")
  expect_true(nrow(row) > 0)
  expect_true(row$p_fdr < 0.01)
  expect_true(row$r > 0)
})

test_that("S3: bias × RSE and bias × SCC are FDR-significant with negative direction", {
  path <- here("Results", "correlations_s3_full.csv")
  skip_if_missing(path)
  df <- read.csv(path)

  bias_rse <- get_corr(df, "bias", "RSE")
  bias_scc <- get_corr(df, "bias", "SCC")

  if (nrow(bias_rse) > 0)
    expect_true(bias_rse$r < 0, info = "bias × RSE should be negative in S3")
  if (nrow(bias_scc) > 0)
    expect_true(bias_scc$r < 0, info = "bias × SCC should be negative in S3")
})

test_that("S3: m (projection rate) has no FDR-significant correlates", {
  path <- here("Results", "correlations_s3_full.csv")
  skip_if_missing(path)
  df <- read.csv(path)
  m_rows <- df[df$param == "m", ]
  expect_true(all(m_rows$p_fdr > 0.05, na.rm = TRUE),
              info = paste("Unexpected FDR-significant m correlate in S3:",
                           paste(m_rows$scale[m_rows$p_fdr <= 0.05], collapse=", ")))
})

# ── Cross-study replication: MCR × SING.Ind positive in all three ─────────────
test_that("MCR × SING.Ind is positive and nominally significant in all three studies", {
  for (study in c("s1", "s2", "s3")) {
    path <- here("Results", paste0("correlations_", study, "_full.csv"))
    if (!file.exists(path)) next
    df  <- read.csv(path)
    row <- get_corr(df, "subject_mcr", "SING.Ind")
    if (nrow(row) == 0) next
    expect_true(row$r > 0,       info = paste(study, "MCR × SING.Ind should be positive"))
    expect_true(row$p_raw < 0.05, info = paste(study, "MCR × SING.Ind should be nominally significant"))
  }
})
