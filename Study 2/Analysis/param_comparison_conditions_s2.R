# param_comparison_conditions_s2.R
# Compare model-derived parameters across outgroup contrast conditions — Study 2
#
# Conditions: "Not UCR" (Negation), "UCLA" (High-Status), "CSU LA" (Low-Status)
#
# Prerequisites:
#   - loo_delta_indiff_s2.R must have run first (produces ind_diffs_s2_full.csv)
#   - model_comparison_s2_results.csv must be current
#
# Analyses:
#   1. Group-level back-transformed parameters from winning model summary
#   2. Individual-level: lm(param ~ condition) for each parameter
#   3. Pairwise contrasts (High-Status vs Negation, Low-Status vs Negation)
#   4. λ by condition (MANUSCRIPT PENDING)
#   5. r(λ, outgroup_warmth) (MANUSCRIPT PENDING)
#   6. m_in vs m_out by condition (if asym_lambda wins)
#   7. MCR by condition
#
# Output:
#   Results/param_by_condition_s2_desc.csv
#   Results/param_by_condition_s2_lm.csv
#   Results/ind_diffs_s2_full_enriched.csv

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# ── 0. Determine winning model ─────────────────────────────────────────────────
comp_file <- here("Results", "model_comparison_s2_results.csv")
if (!file.exists(comp_file)) stop("model_comparison_s2_results.csv not found. Run run_model_comparison_s2_remaining.R first.")
comp   <- read.csv(comp_file, row.names = 1)
winner <- rownames(comp)[1]
cat(sprintf("Winning model: %s\n", winner))
cat("LOO comparison:\n")
print(round(comp[, c("elpd_loo","se_elpd_loo","elpd_diff","se_diff")], 2))

# ── 1. Load individual differences data ────────────────────────────────────────
ind_file <- here("Results", "ind_diffs_s2_full.csv")
if (!file.exists(ind_file)) stop("ind_diffs_s2_full.csv not found. Run loo_delta_indiff_s2.R first.")
ind_df <- read.csv(ind_file)
cat(sprintf("\nN subjects: %d\n", nrow(ind_df)))

# ── 2. Ensure outgroup column + recode to readable labels ──────────────────────
if (!"outgroup" %in% names(ind_df)) {
  fulldf_cond <- read.csv(here("Study 2/Cleaning/output/fullTest.csv")) |>
    filter(!is.na(ingChoiceN)) |>
    distinct(subID, outgroup)
  ind_df <- left_join(ind_df, fulldf_cond, by = "subID")
}

ind_df$condition <- factor(ind_df$outgroup,
                            levels = c("Not UCR", "UCLA", "CSU LA"),
                            labels = c("Negation", "High-Status", "Low-Status"))

cat("\nN by condition:\n")
print(table(ind_df$condition))

# ── 3. Group-level parameters from winning model summary ──────────────────────
# Back-transform mu_pr (probit scale) to original parameter scale.
# Scales: S_Sym_Lambda:  m × 10, bias × 1 (Phi), lambda × 5, w × 1
#         S_Asym_Lambda: m_in × 10, m_out × 10, bias × 1, lambda × 5, w × 1
summary_win <- tryCatch(
  read.csv(here("Results", paste0("summary_s2_", winner, ".csv"))),
  error = function(e) { message("Cannot read summary CSV for winner: ", winner); NULL }
)

if (!is.null(summary_win)) {
  mu_pr_rows <- summary_win |> filter(grepl("^\"?mu_pr\\[", variable))

  if (winner == "sym_lambda") {
    scales <- c(10, 1, 5, 1)
    param_nm <- c("m", "bias", "lambda", "w")
  } else if (winner == "asym_lambda") {
    scales <- c(10, 10, 1, 5, 1)
    param_nm <- c("m_in", "m_out", "bias", "lambda", "w")
  } else if (winner == "symmetric") {
    scales <- c(10, 1, 1)
    param_nm <- c("m", "bias", "w")
  } else {
    scales <- c(1, 1)
    param_nm <- c("bias", "w")
  }

  n_p <- min(nrow(mu_pr_rows), length(scales))
  mu_pr_rows <- mu_pr_rows[seq_len(n_p), ]
  mu_pr_rows$param    <- param_nm[seq_len(n_p)]
  mu_pr_rows$med_bt   <- pnorm(mu_pr_rows$median) * scales[seq_len(n_p)]
  mu_pr_rows$q5_bt    <- pnorm(mu_pr_rows$q5)     * scales[seq_len(n_p)]
  mu_pr_rows$q95_bt   <- pnorm(mu_pr_rows$q95)    * scales[seq_len(n_p)]

  cat("\n=== Group-level parameters (back-transformed from probit scale) ===\n")
  print(mu_pr_rows[, c("param", "med_bt", "q5_bt", "q95_bt")], row.names = FALSE, digits = 3)
}

# ── 4. Individual-level parameters: descriptives by condition ─────────────────
param_cols <- intersect(
  c("m", "m_in", "m_out", "lambda", "bias", "w", "subject_mcr", "delta_elpd"),
  names(ind_df)
)
cat(sprintf("\nParameter columns: %s\n", paste(param_cols, collapse = ", ")))

desc_list <- list()
for (p in param_cols) {
  d <- ind_df |>
    group_by(condition) |>
    summarise(
      M      = round(mean(.data[[p]], na.rm = TRUE), 3),
      SD     = round(sd(.data[[p]],   na.rm = TRUE), 3),
      Mdn    = round(median(.data[[p]], na.rm = TRUE), 3),
      n      = sum(!is.na(.data[[p]])),
      .groups = "drop"
    ) |>
    mutate(param = p)
  desc_list[[p]] <- d

  cat(sprintf("\n── %s by condition ──\n", p))
  print(d, row.names = FALSE)
}
desc_df <- do.call(rbind, desc_list)

# ── 5. Linear models: param ~ condition (reference = Negation) ────────────────
lm_results <- list()
for (p in param_cols) {
  if (all(is.na(ind_df[[p]]))) next
  fit <- lm(reformulate("condition", response = p), data = ind_df)
  s   <- summary(fit)
  f   <- s$fstatistic
  p_f <- pf(f[1], f[2], f[3], lower.tail = FALSE)

  coefs <- as.data.frame(coef(s))
  coefs$term  <- rownames(coefs)
  coefs$param <- p
  coefs$F_stat <- f[1]
  coefs$F_df1  <- f[2]
  coefs$F_df2  <- f[3]
  coefs$p_F    <- p_f
  coefs$R2     <- s$r.squared
  lm_results[[p]] <- coefs

  cat(sprintf("\n── lm(%s ~ condition) ──\n", p))
  cat(sprintf("  F(%d, %d) = %.3f, p = %.4f, R² = %.3f\n",
              as.integer(f[2]), as.integer(f[3]), f[1], p_f, s$r.squared))
  print(round(coef(s)[, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")], 4))
}
lm_df <- do.call(rbind, lm_results)

# ── 6. MANUSCRIPT PENDING: λ by condition ─────────────────────────────────────
if ("lambda" %in% names(ind_df)) {
  cat("\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat("MANUSCRIPT PENDING: lambda (λ) by condition\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")

  fit_lam <- lm(lambda ~ condition, data = ind_df)
  s_lam   <- summary(fit_lam)
  f_lam   <- s_lam$fstatistic
  p_f_lam <- pf(f_lam[1], f_lam[2], f_lam[3], lower.tail = FALSE)

  lam_means <- ind_df |>
    group_by(condition) |>
    summarise(M = round(mean(lambda, na.rm = TRUE), 3),
              SD = round(sd(lambda, na.rm = TRUE), 3),
              .groups = "drop")

  cat("\nλ by condition (individual medians):\n")
  for (i in seq_len(nrow(lam_means))) {
    cat(sprintf("  %s: M = %.3f, SD = %.3f\n",
                lam_means$condition[i], lam_means$M[i], lam_means$SD[i]))
  }
  cat(sprintf("\nOverall ANOVA: F(%d, %d) = %.3f, p = %.4f, R² = %.3f\n",
              as.integer(f_lam[2]), as.integer(f_lam[3]),
              f_lam[1], p_f_lam, s_lam$r.squared))
  cat("Coefficients (reference = Negation):\n")
  print(round(coef(s_lam)[, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")], 4))

  # Save lambda means for manuscript insertion
  write.csv(lam_means, here("Results", "lambda_by_condition_s2.csv"), row.names = FALSE)
}

# ── 7. MANUSCRIPT PENDING: r(λ, outgroup_warmth) ─────────────────────────────
cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("MANUSCRIPT PENDING: r(lambda, outgroup warmth)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

fulldf_warm <- read.csv(here("Study 2/Cleaning/output/fullTest.csv")) |>
  filter(!is.na(ingChoiceN)) |>
  distinct(subID, outgroup, InUCLATherm, InCSULATherm)

fulldf_warm <- fulldf_warm |>
  mutate(outgroup_warmth = case_when(
    outgroup == "UCLA"   ~ InUCLATherm,
    outgroup == "CSU LA" ~ InCSULATherm,
    TRUE                 ~ NA_real_
  )) |>
  mutate(condition = factor(outgroup,
                            levels = c("Not UCR", "UCLA", "CSU LA"),
                            labels = c("Negation", "High-Status", "Low-Status")))

ind_df2 <- left_join(ind_df, fulldf_warm[, c("subID", "outgroup_warmth")], by = "subID")

if ("lambda" %in% names(ind_df2)) {
  # Overall (UCLA + CSU LA conditions where warmth is available)
  ct_all <- cor.test(ind_df2$lambda, ind_df2$outgroup_warmth, use = "complete.obs")
  n_ct   <- sum(complete.cases(ind_df2[, c("lambda", "outgroup_warmth")]))
  cat(sprintf("\nOverall r(λ, warmth) [High-Status + Low-Status]: r = %.3f, p = %.4f, n = %d\n",
              ct_all$estimate, ct_all$p.value, n_ct))
  cat(sprintf("  95%% CI [%.3f, %.3f]\n", ct_all$conf.int[1], ct_all$conf.int[2]))

  # By condition
  for (cond in c("High-Status", "Low-Status")) {
    sub <- filter(ind_df2, condition == cond)
    n_c <- sum(complete.cases(sub[, c("lambda", "outgroup_warmth")]))
    if (n_c > 5) {
      ct_c <- cor.test(sub$lambda, sub$outgroup_warmth, use = "complete.obs")
      cat(sprintf("  %s: r = %.3f, p = %.4f, n = %d (CI [%.3f, %.3f])\n",
                  cond, ct_c$estimate, ct_c$p.value, n_c,
                  ct_c$conf.int[1], ct_c$conf.int[2]))
    }
  }
}

# ── 8. γ (ingroup bias) × outgroup preference: lm(bias ~ condition * warmth) ──
cat("\n=== γ × outgroup preference thermometer (High-Status + Low-Status only) ===\n")
cat("outgroup_warmth = UCR thermometer minus outgroup thermometer (positive = prefer UCR over outgroup)\n\n")

if ("outgroup_warmth" %in% names(ind_df2) && "bias" %in% names(ind_df2)) {
  df_bias <- ind_df2 |>
    filter(!is.na(outgroup_warmth), !is.na(bias),
           condition %in% c("High-Status", "Low-Status")) |>
    mutate(
      condition = droplevels(condition),
      outgroup_warmth.Z = as.numeric(scale(outgroup_warmth))
    )

  cat(sprintf("N = %d (High-Status: %d, Low-Status: %d)\n",
    nrow(df_bias),
    sum(df_bias$condition == "High-Status"),
    sum(df_bias$condition == "Low-Status")))

  mod_int  <- lm(bias ~ condition * outgroup_warmth.Z, data = df_bias)
  mod_main <- lm(bias ~ condition + outgroup_warmth.Z, data = df_bias)
  av <- anova(mod_main, mod_int)

  cat(sprintf("\nInteraction: F(%.0f, %.0f) = %.3f, p = %.4f\n",
    av$Df[2], av$Res.Df[2], av$F[2], av$`Pr(>F)`[2]))

  cat("\nFull interaction model coefficients:\n")
  print(round(summary(mod_int)$coefficients, 4))

  cat("\nOverall warmth slope (main effects model):\n")
  print(round(summary(mod_main)$coefficients["outgroup_warmth.Z", ], 4))

  cat("\nWithin-condition slopes:\n")
  for (cond in c("High-Status", "Low-Status")) {
    sub <- filter(df_bias, condition == cond)
    m <- lm(bias ~ outgroup_warmth.Z, data = sub)
    co <- summary(m)$coefficients["outgroup_warmth.Z", ]
    cat(sprintf("  %s (n=%d): b=%.4f, SE=%.4f, t(%.0f)=%.3f, p=%.4f\n",
      cond, nrow(sub), co[1], co[2], df.residual(m), co[3], co[4]))
  }
}

# ── 10. Projection rate by condition: α (or α_in / α_out) ────────────────────
cat("\n=== Projection rate (α) by condition ===\n")
m_cols <- intersect(c("m", "m_in", "m_out"), names(ind_df))
for (mc in m_cols) {
  fit_m <- lm(reformulate("condition", response = mc), data = ind_df)
  s_m   <- summary(fit_m)
  f_m   <- s_m$fstatistic
  p_m   <- pf(f_m[1], f_m[2], f_m[3], lower.tail = FALSE)
  cat(sprintf("\n%s ~ condition: F(%d,%d) = %.3f, p = %.4f, R² = %.3f\n",
              mc, as.integer(f_m[2]), as.integer(f_m[3]), f_m[1], p_m, s_m$r.squared))
  print(round(coef(s_m)[, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")], 4))
}

# If asym_lambda wins: within-condition m_in vs m_out comparison
if (winner == "asym_lambda" && all(c("m_in","m_out") %in% names(ind_df))) {
  cat("\n=== m_in vs m_out (within-condition paired tests) ===\n")
  for (cond in levels(ind_df$condition)) {
    sub <- filter(ind_df, condition == cond)
    if (nrow(sub) > 5 && !all(is.na(sub$m_in)) && !all(is.na(sub$m_out))) {
      tt <- t.test(sub$m_in, sub$m_out, paired = TRUE)
      cat(sprintf("  %s: m_in = %.3f, m_out = %.3f | t(%d) = %.3f, p = %.4f\n",
                  cond, mean(sub$m_in, na.rm = TRUE), mean(sub$m_out, na.rm = TRUE),
                  tt$parameter, tt$statistic, tt$p.value))
    }
  }
}

# ── 11. MCR by condition ──────────────────────────────────────────────────────
if ("subject_mcr" %in% names(ind_df)) {
  cat("\n=== MCR by condition ===\n")
  fit_mcr <- lm(subject_mcr ~ condition, data = ind_df)
  s_mcr   <- summary(fit_mcr)
  f_mcr   <- s_mcr$fstatistic
  p_mcr   <- pf(f_mcr[1], f_mcr[2], f_mcr[3], lower.tail = FALSE)
  cat(sprintf("F(%d,%d) = %.3f, p = %.4f, R² = %.3f\n",
              as.integer(f_mcr[2]), as.integer(f_mcr[3]),
              f_mcr[1], p_mcr, s_mcr$r.squared))
  print(round(coef(s_mcr)[, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")], 4))

  mcr_means <- ind_df |>
    group_by(condition) |>
    summarise(M = round(mean(subject_mcr, na.rm=T), 3),
              SD = round(sd(subject_mcr, na.rm=T), 3),
              .groups = "drop")
  cat("\nMCR by condition:\n")
  print(mcr_means, row.names = FALSE)
}

# ── 12. Save enriched individual differences ──────────────────────────────────
write.csv(desc_df,  here("Results", "param_by_condition_s2_desc.csv"), row.names = FALSE)
write.csv(lm_df,    here("Results", "param_by_condition_s2_lm.csv"),   row.names = FALSE)
write.csv(ind_df2,  here("Results", "ind_diffs_s2_full_enriched.csv"), row.names = FALSE)

message("\nSaved:")
message("  Results/param_by_condition_s2_desc.csv")
message("  Results/param_by_condition_s2_lm.csv")
message("  Results/ind_diffs_s2_full_enriched.csv")
