# tests/testthat/test-figure-helpers.R
# Unit tests for pure helper functions used across figure-making scripts.
# Functions are redefined here so tests run without sourcing the full scripts.

library(testthat)

# ── Helpers under test ────────────────────────────────────────────────────────

# From make_figures.R
bt <- function(x, scale = 1) pnorm(x) * scale

ceil2 <- function(x) ceiling(x * 20) / 20

# From make_figures.R — gradient data builder
make_grad_df <- function(cond_df, cond_col = "condition") {
  S_seq <- seq(0, 1, by = 0.05)
  do.call(rbind, lapply(seq_len(nrow(cond_df)), function(i) {
    lam    <- cond_df$M[i]
    lam_lo <- max(0.01, lam - cond_df$SD[i])
    lam_hi <- lam + cond_df$SD[i]
    data.frame(
      S         = S_seq,
      g_med     = S_seq ^ lam,
      g_lo      = S_seq ^ lam_hi,
      g_hi      = S_seq ^ lam_lo,
      condition = cond_df[[cond_col]][i]
    )
  }))
}

# From make_param_recovery_fig.R
gen_placeholder <- function(model_name, params, n = 61, seed = 42) {
  set.seed(seed)
  param_ranges <- list(
    m      = c(0.3, 9.7, 2.80),
    m_in   = c(0.3, 9.7, 3.50),
    m_out  = c(0.3, 9.7, 2.50),
    bias   = c(0.05, 0.95, 0.36),
    lambda = c(0.3, 4.9,  3.56),
    w      = c(0.10, 0.90, 0.56)
  )
  do.call(rbind, lapply(params, function(p) {
    rng       <- param_ranges[[p]]
    rng_width <- rng[2] - rng[1]
    true_pr   <- rnorm(n, qnorm((rng[3] - rng[1]) / rng_width), 0.5)
    true_vals <- pmin(rng[2], pmax(rng[1], pnorm(true_pr) * rng_width + rng[1]))
    noise     <- rnorm(n, 0, rng_width * 0.22)
    rec_vals  <- pmin(rng[2], pmax(rng[1], true_vals + noise))
    data.frame(model = model_name, param = p, subj_idx = seq_len(n),
               true = true_vals, recovered = rec_vals)
  }))
}

# ── bt() ─────────────────────────────────────────────────────────────────────
test_that("bt() maps 0 to 0.5 * scale", {
  expect_equal(bt(0),     0.5)
  expect_equal(bt(0, 10), 5.0)
  expect_equal(bt(0, 5),  2.5)
})

test_that("bt() approaches scale at +Inf and 0 at -Inf", {
  expect_equal(bt(Inf,  1),  1,  tolerance = 1e-9)
  expect_equal(bt(-Inf, 1),  0,  tolerance = 1e-9)
  expect_equal(bt(Inf,  10), 10, tolerance = 1e-9)
})

test_that("bt() round-trips through qnorm", {
  vals <- c(0.10, 0.36, 0.50, 0.75, 0.90)
  expect_equal(bt(qnorm(vals)), vals, tolerance = 1e-8)
  expect_equal(bt(qnorm(vals), 5), vals * 5, tolerance = 1e-8)
})

test_that("bt() is monotonically increasing in x", {
  x <- seq(-3, 3, by = 0.5)
  expect_true(all(diff(bt(x)) > 0))
})

# ── ceil2() ──────────────────────────────────────────────────────────────────
test_that("ceil2() rounds up to nearest 0.05", {
  expect_equal(ceil2(0.10), 0.10)
  expect_equal(ceil2(0.11), 0.15)
  expect_equal(ceil2(0.20), 0.20)
  expect_equal(ceil2(0.21), 0.25)
  expect_equal(ceil2(0.91), 0.95)
  expect_equal(ceil2(0.95), 0.95)
})

test_that("ceil2() output is always a multiple of 0.05", {
  vals <- seq(0, 1, by = 0.013)
  result <- ceil2(vals)
  expect_true(all(abs(result * 20 - round(result * 20)) < 1e-9))
})

# ── make_grad_df() ────────────────────────────────────────────────────────────
test_that("make_grad_df() returns correct columns and structure", {
  cond_df <- data.frame(condition = c("A", "B"), M = c(3.5, 2.0), SD = c(0.3, 0.2))
  out <- make_grad_df(cond_df)
  expect_true(all(c("S", "g_med", "g_lo", "g_hi", "condition") %in% names(out)))
  expect_setequal(unique(out$condition), c("A", "B"))
})

test_that("make_grad_df() g values are in [0, 1]", {
  cond_df <- data.frame(condition = "A", M = 3.5, SD = 0.5)
  out <- make_grad_df(cond_df)
  expect_true(all(out$g_med >= 0 & out$g_med <= 1))
  expect_true(all(out$g_lo  >= 0 & out$g_lo  <= 1))
  expect_true(all(out$g_hi  >= 0 & out$g_hi  <= 1))
})

test_that("make_grad_df() g_med decays with decreasing S", {
  cond_df <- data.frame(condition = "A", M = 3.5, SD = 0.1)
  out <- make_grad_df(cond_df)
  g_ordered <- out$g_med[order(out$S)]
  expect_true(all(diff(g_ordered) >= 0))  # monotonically non-decreasing with S
})

test_that("make_grad_df() CI ordering: g_lo <= g_med <= g_hi", {
  # Higher lambda → steeper decay → lower g; so g_lo uses lam_hi and g_hi uses lam_lo
  cond_df <- data.frame(condition = "A", M = 3.5, SD = 0.3)
  out <- make_grad_df(cond_df)
  # At S < 1, higher lambda means lower weight: g_lo <= g_med <= g_hi
  mid <- out[out$S > 0.1 & out$S < 0.9, ]
  expect_true(all(mid$g_lo <= mid$g_med + 1e-9))
  expect_true(all(mid$g_med <= mid$g_hi + 1e-9))
})

test_that("make_grad_df() handles SD = 0 (no ribbon spread)", {
  cond_df <- data.frame(condition = "A", M = 3.5, SD = 0)
  out <- make_grad_df(cond_df)
  expect_equal(out$g_lo, out$g_med, tolerance = 1e-9)
  expect_equal(out$g_hi, out$g_med, tolerance = 1e-9)
})

# ── gen_placeholder() ─────────────────────────────────────────────────────────
test_that("gen_placeholder() returns correct columns", {
  out <- gen_placeholder("Sym_Lambda", c("m", "bias", "lambda", "w"), n = 30)
  expect_true(all(c("model", "param", "subj_idx", "true", "recovered") %in% names(out)))
})

test_that("gen_placeholder() returns n rows per parameter", {
  n     <- 40
  params <- c("m", "bias", "lambda", "w")
  out   <- gen_placeholder("Sym_Lambda", params, n = n)
  expect_equal(nrow(out), n * length(params))
  for (p in params) {
    expect_equal(nrow(out[out$param == p, ]), n)
  }
})

test_that("gen_placeholder() values stay within parameter ranges", {
  out <- gen_placeholder("Sym_Lambda", c("m", "bias", "lambda", "w"), n = 200, seed = 1)
  m_rows     <- out[out$param == "m", ]
  bias_rows  <- out[out$param == "bias", ]
  lam_rows   <- out[out$param == "lambda", ]
  w_rows     <- out[out$param == "w", ]

  expect_true(all(m_rows$true     >= 0.3  & m_rows$true     <= 9.7))
  expect_true(all(bias_rows$true  >= 0.05 & bias_rows$true  <= 0.95))
  expect_true(all(lam_rows$true   >= 0.3  & lam_rows$true   <= 4.9))
  expect_true(all(w_rows$true     >= 0.10 & w_rows$true     <= 0.90))
})

test_that("gen_placeholder() is reproducible with same seed", {
  out1 <- gen_placeholder("Sym_Lambda", c("m", "bias"), n = 20, seed = 99)
  out2 <- gen_placeholder("Sym_Lambda", c("m", "bias"), n = 20, seed = 99)
  expect_equal(out1$true, out2$true)
})
