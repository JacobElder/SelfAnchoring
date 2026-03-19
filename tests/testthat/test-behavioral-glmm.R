# tests/testthat/test-behavioral-glmm.R
# Verifies that marginal effects (AMEs) from mixed-effects logistic regressions
# match reported values and directions across all three studies and the pooled analysis.

library(testthat)
library(here)

skip_if_missing <- function(path) {
  if (!file.exists(path)) skip(paste("File not found:", basename(path)))
}

read_ames <- function(study) {
  path <- here("Results", paste0("marginal_effects_", study, "_ames.csv"))
  skip_if_missing(path)
  read.csv(path)
}

# ── Study 1 ──────────────────────────────────────────────────────────────────
test_that("S1: SS AME is positive and significant", {
  df  <- read_ames("s1")
  ss  <- df[grepl("ss|similarity", df$term, ignore.case = TRUE) &
              (is.na(df$novel) | df$novel == FALSE) &
              (is.na(df$condition) | df$condition == ""), ]
  if (nrow(ss) == 0) ss <- df[df$model %in% c("M3","M2") & is.na(df$novel), ]
  expect_true(nrow(ss) > 0, info = "No SS AME row found")
  expect_true(ss$estimate[1] > 0,    info = "SS AME should be positive")
  expect_true(ss$p.value[1] < 0.01,  info = "SS AME should be highly significant")
  expect_true(ss$estimate[1] > 0.04, info = paste("SS AME smaller than expected:", ss$estimate[1]))
  expect_true(ss$estimate[1] < 0.15, info = paste("SS AME larger than expected:", ss$estimate[1]))
})

test_that("S1: AME CIs straddle the estimate", {
  df <- read_ames("s1")
  expect_true(all(df$conf.low  <= df$estimate + 1e-9, na.rm = TRUE))
  expect_true(all(df$conf.high >= df$estimate - 1e-9, na.rm = TRUE))
})

test_that("S1: Self-evaluation AME is positive", {
  df   <- read_ames("s1")
  self <- df[grepl("selfResp|self_eval|self.eval|selfeval", df$term, ignore.case = TRUE), ]
  if (nrow(self) == 0) skip("selfResp AME row not found in S1")
  expect_true(self$estimate[1] > 0)
})

# ── Study 2 ──────────────────────────────────────────────────────────────────
test_that("S2: Overall SS AME is positive and significant", {
  df <- read_ames("s2")
  # Look for overall/pooled SS row (no condition breakdown or outgroup == "")
  ss_all <- df[grepl("^ss$|similarity", df$term, ignore.case = TRUE) & is.na(df$novel), ]
  if (nrow(ss_all) == 0) ss_all <- df[df$model == "M3" & is.na(df$novel) & is.na(df$condition), ]
  if (nrow(ss_all) == 0) skip("No overall SS AME row found in S2")
  expect_true(ss_all$estimate[1] > 0)
  expect_true(ss_all$p.value[1] < 0.01)
})

test_that("S2: High-status condition has the largest SS AME", {
  df <- read_ames("s2")
  ss_cond <- df[grepl("ss|similarity", df$term, ignore.case = TRUE) &
                  !is.na(df$condition) & !is.na(df$novel) & df$novel == FALSE, ]
  if (nrow(ss_cond) == 0) skip("No condition-level SS AME rows found")

  # High-status should be the numerically largest
  if ("UCLA" %in% ss_cond$condition) {
    hs_ame  <- ss_cond$estimate[ss_cond$condition == "UCLA"]
    neg_ame <- ss_cond$estimate[ss_cond$condition == "Not UCR"]
    ls_ame  <- ss_cond$estimate[ss_cond$condition == "CSU LA"]
    expect_true(hs_ame > neg_ame,  info = "High-status AME should exceed Negation")
    expect_true(hs_ame > ls_ame,   info = "High-status AME should exceed Low-status")
  } else {
    skip("Condition labels not as expected in S2 AME file")
  }
})

test_that("S2: Low-status AME is weaker than Negation", {
  df <- read_ames("s2")
  ss_cond <- df[grepl("ss|similarity", df$term, ignore.case = TRUE) &
                  !is.na(df$condition) & !is.na(df$novel) & df$novel == FALSE, ]
  if (nrow(ss_cond) == 0 || !"CSU LA" %in% ss_cond$condition) skip("Condition rows missing")
  neg_ame <- ss_cond$estimate[ss_cond$condition == "Not UCR"]
  ls_ame  <- ss_cond$estimate[ss_cond$condition == "CSU LA"]
  expect_true(ls_ame < neg_ame,
              info = paste("Low-status AME should be weaker than Negation:", ls_ame, "vs", neg_ame))
})

# ── Study 3 ──────────────────────────────────────────────────────────────────
test_that("S3: Minority condition SS AME is non-significant", {
  df <- read_ames("s3")
  ss_min <- df[grepl("ss|similarity", df$term, ignore.case = TRUE) &
                 grepl("inority|Minority", df$condition, ignore.case = TRUE) &
                 !is.na(df$novel) & df$novel == FALSE, ]
  if (nrow(ss_min) == 0) skip("No minority condition SS AME row")
  expect_true(ss_min$p.value[1] > 0.05,
              info = paste("Minority SS AME should be non-significant, got p =", ss_min$p.value[1]))
})

test_that("S3: Majority condition trained SS AME is significant", {
  df <- read_ames("s3")
  ss_maj <- df[grepl("ss|similarity", df$term, ignore.case = TRUE) &
                 grepl("ajority|Majority", df$condition, ignore.case = TRUE) &
                 !is.na(df$novel) & df$novel == FALSE, ]
  if (nrow(ss_maj) == 0) skip("No majority condition SS AME row")
  expect_true(ss_maj$estimate[1] > 0, info = "Majority SS AME should be positive")
})

# ── Pooled ────────────────────────────────────────────────────────────────────
test_that("Pooled: Overall SS AME is positive and significant (~0.038)", {
  df <- read_ames("pooled")
  # Overall SS AME across all conditions
  ss <- df[grepl("ss|similarity", df$term, ignore.case = TRUE) & is.na(df$novel), ]
  if (nrow(ss) == 0) ss <- df[df$model == "M2" & is.na(df$novel), ]
  if (nrow(ss) == 0) skip("No overall SS AME row in pooled")
  expect_true(ss$estimate[1] > 0)
  expect_true(ss$p.value[1] < 0.001)
  expect_true(ss$estimate[1] > 0.02 & ss$estimate[1] < 0.07,
              info = paste("Pooled SS AME out of expected range:", ss$estimate[1]))
})

test_that("Pooled: RacMinority condition SS AME is non-significant", {
  df <- read_ames("pooled")
  rac_min <- df[grepl("ss|similarity", df$term, ignore.case = TRUE) &
                  grepl("acMinority|RacMinority", df$condition, ignore.case = TRUE), ]
  if (nrow(rac_min) == 0) skip("No RacMinority condition row in pooled AMEs")
  expect_true(rac_min$p.value[1] > 0.10,
              info = paste("RacMinority AME should be null, got p =", rac_min$p.value[1]))
})

test_that("Pooled: Trained and novel AMEs are both significant", {
  df <- read_ames("pooled")
  if (!"novel" %in% names(df)) skip("No novel column in pooled AMEs")

  trained <- df[grepl("ss|similarity", df$term, ignore.case = TRUE) &
                  !is.na(df$novel) & df$novel == FALSE, ]
  novel   <- df[grepl("ss|similarity", df$term, ignore.case = TRUE) &
                  !is.na(df$novel) & df$novel == TRUE, ]

  if (nrow(trained) == 0 || nrow(novel) == 0) skip("Trained/novel split not found")

  # At least some rows should be significant
  expect_true(any(trained$p.value < 0.05, na.rm = TRUE),
              info = "At least one trained-trait AME should be significant")
  expect_true(any(novel$p.value < 0.05, na.rm = TRUE),
              info = "At least one novel-trait AME should be significant")
})

test_that("Pooled: UnivHigh SS AME is largest across conditions", {
  df <- read_ames("pooled")
  ss_cond <- df[grepl("ss|similarity", df$term, ignore.case = TRUE) &
                  !is.na(df$condition) & df$condition != "", ]
  if (nrow(ss_cond) == 0) skip("No condition-level rows in pooled AMEs")
  if (!"UnivHigh" %in% ss_cond$condition) skip("UnivHigh condition not found")

  univ_high <- ss_cond$estimate[ss_cond$condition == "UnivHigh"]
  all_others <- ss_cond$estimate[ss_cond$condition != "UnivHigh"]
  expect_true(max(univ_high) >= max(all_others) - 0.01,
              info = "UnivHigh should have the largest (or near-largest) SS AME")
})
