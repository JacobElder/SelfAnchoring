# param_comparison_conditions_s3.R
# Compare model-derived parameters across conditions — Study 3
# Parallel to Study 2/Analysis/param_comparison_conditions_s2.R
#
# Conditions: "Minority" (Asian vs Latino), "Majority" (Ingroup vs White)
#
# Prerequisites:
#   - loo_delta_indiff_s3.R must have run first (produces ind_diffs_s3_full.csv)
#   - model_comparison_s3_results.csv must be current
#
# Analyses:
#   1. Group-level back-transformed parameters from winning model summary
#   2. Individual-level: lm(param ~ condition) for each parameter
#   3. λ by condition (MANUSCRIPT PENDING)
#   4. r(λ, outgroup_warmth) by condition (MANUSCRIPT PENDING)
#   5. m_in vs m_out within each condition (if asym_lambda wins)
#   6. MCR by condition
#
# Output:
#   Results/param_by_condition_s3_desc.csv
#   Results/param_by_condition_s3_lm.csv
#   Results/ind_diffs_s3_full_enriched.csv

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# ── 0. Determine winning model ─────────────────────────────────────────────────
comp_file <- here("Results", "model_comparison_s3_results.csv")
if (!file.exists(comp_file)) stop("model_comparison_s3_results.csv not found. Run run_model_comparison_s3.R first.")
comp   <- read.csv(comp_file, row.names = 1)
winner <- rownames(comp)[1]
cat(sprintf("Winning model: %s\n", winner))
cat("LOO comparison:\n")
print(round(comp[, c("elpd_loo","se_elpd_loo","elpd_diff","se_diff")], 2))

# ── 1. Load individual differences data ────────────────────────────────────────
ind_file <- here("Results", "ind_diffs_s3_full.csv")
if (!file.exists(ind_file)) stop("ind_diffs_s3_full.csv not found. Run loo_delta_indiff_s3.R first.")
ind_df <- read.csv(ind_file)
cat(sprintf("\nN subjects: %d\n", nrow(ind_df)))

# ── 2. Ensure condition column is present ──────────────────────────────────────
# condition should already be in ind_diffs from loo_delta_indiff_s3.R
if (!"condition" %in% names(ind_df)) {
  fulldf_cond <- read.csv(here("Study 3/Cleaning/output/fullTest_fixed.csv")) |>
    filter(!is.na(ingChoiceN)) |>
    distinct(subID, condition)
  ind_df <- left_join(ind_df, fulldf_cond, by = "subID")
}

ind_df$condition <- factor(ind_df$condition, levels = c("Minority", "Majority"))
cat("\nN by condition:\n"); print(table(ind_df$condition))

# ── 3. Group-level parameters from winning model summary ──────────────────────
# Back-transform mu_pr (probit scale) to original parameter scale.
summary_win <- tryCatch(
  read.csv(here("Results", paste0("summary_s3_", winner, ".csv"))),
  error = function(e) { message("Cannot read summary CSV for: ", winner); NULL }
)

if (!is.null(summary_win)) {
  mu_pr_rows <- summary_win |> filter(grepl("^\"?mu_pr\\[", variable))

  if (winner == "sym_lambda") {
    scales <- c(10, 1, 5, 1); param_nm <- c("m", "bias", "lambda", "w")
  } else if (winner == "asym_lambda") {
    scales <- c(10, 10, 1, 5, 1); param_nm <- c("m_in", "m_out", "bias", "lambda", "w")
  } else if (winner == "symmetric") {
    scales <- c(10, 1, 1); param_nm <- c("m", "bias", "w")
  } else {
    scales <- c(1, 1); param_nm <- c("bias", "w")
  }

  n_p <- min(nrow(mu_pr_rows), length(scales))
  mu_pr_rows <- mu_pr_rows[seq_len(n_p), ]
  mu_pr_rows$param  <- param_nm[seq_len(n_p)]
  mu_pr_rows$med_bt <- pnorm(mu_pr_rows$median) * scales[seq_len(n_p)]
  mu_pr_rows$q5_bt  <- pnorm(mu_pr_rows$q5)     * scales[seq_len(n_p)]
  mu_pr_rows$q95_bt <- pnorm(mu_pr_rows$q95)    * scales[seq_len(n_p)]

  cat("\n=== Group-level parameters (back-transformed) ===\n")
  print(mu_pr_rows[, c("param","med_bt","q5_bt","q95_bt")], row.names = FALSE, digits = 3)
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
    summarise(M = round(mean(.data[[p]], na.rm = TRUE), 3),
              SD = round(sd(.data[[p]], na.rm = TRUE), 3),
              Mdn = round(median(.data[[p]], na.rm = TRUE), 3),
              n = sum(!is.na(.data[[p]])), .groups = "drop") |>
    mutate(param = p)
  desc_list[[p]] <- d
  cat(sprintf("\n── %s by condition ──\n", p)); print(d, row.names = FALSE)
}
desc_df <- do.call(rbind, desc_list)

# ── 5. Linear models: param ~ condition (reference = Minority) ────────────────
# With only 2 conditions, this is equivalent to a t-test
lm_results <- list()
for (p in param_cols) {
  if (all(is.na(ind_df[[p]]))) next
  fit <- lm(reformulate("condition", response = p), data = ind_df)
  s   <- summary(fit)
  f   <- s$fstatistic
  p_f <- pf(f[1], f[2], f[3], lower.tail = FALSE)

  coefs <- as.data.frame(coef(s))
  coefs$term <- rownames(coefs); coefs$param <- p
  coefs$F_stat <- f[1]; coefs$F_df1 <- f[2]; coefs$F_df2 <- f[3]
  coefs$p_F <- p_f; coefs$R2 <- s$r.squared
  lm_results[[p]] <- coefs

  cat(sprintf("\n── lm(%s ~ condition) ──\n", p))
  cat(sprintf("  F(%d, %d) = %.3f, p = %.4f, R² = %.3f\n",
              as.integer(f[2]), as.integer(f[3]), f[1], p_f, s$r.squared))
  print(round(coef(s)[, c("Estimate","Std. Error","t value","Pr(>|t|)")], 4))
}
lm_df <- do.call(rbind, lm_results)

# ── 6. MANUSCRIPT PENDING: λ by condition ─────────────────────────────────────
if ("lambda" %in% names(ind_df)) {
  cat("\n"); cat(paste(rep("=", 60), collapse = ""), "\n")
  cat("MANUSCRIPT PENDING: lambda (λ) by condition\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")

  fit_lam <- lm(lambda ~ condition, data = ind_df)
  s_lam   <- summary(fit_lam)
  f_lam   <- s_lam$fstatistic
  p_f_lam <- pf(f_lam[1], f_lam[2], f_lam[3], lower.tail = FALSE)

  lam_means <- ind_df |>
    group_by(condition) |>
    summarise(M = round(mean(lambda, na.rm=T), 3),
              SD = round(sd(lambda, na.rm=T), 3), .groups = "drop")

  cat("\nλ by condition:\n")
  for (i in seq_len(nrow(lam_means))) {
    cat(sprintf("  %s: M = %.3f, SD = %.3f\n",
                lam_means$condition[i], lam_means$M[i], lam_means$SD[i]))
  }
  cat(sprintf("\nANOVA: F(%d, %d) = %.3f, p = %.4f, R² = %.3f\n",
              as.integer(f_lam[2]), as.integer(f_lam[3]),
              f_lam[1], p_f_lam, s_lam$r.squared))
  print(round(coef(s_lam)[, c("Estimate","Std. Error","t value","Pr(>|t|)")], 4))

  write.csv(lam_means, here("Results", "lambda_by_condition_s3.csv"), row.names = FALSE)
}

# ── 7. MANUSCRIPT PENDING: r(λ, outgroup_warmth) ─────────────────────────────
# Study 3 outgroup warmth variables:
#   Minority condition (Asian vs Latino): AsianLatinoTherm
#   Majority condition (Ingroup vs White): inMajorityTherm
#   (Inspect column names; adjust as needed)
cat("\n"); cat(paste(rep("=", 60), collapse = ""), "\n")
cat("MANUSCRIPT PENDING: r(lambda, outgroup warmth)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

fulldf_warm <- read.csv(here("Study 3/Cleaning/output/fullTest_fixed.csv")) |>
  filter(!is.na(ingChoiceN)) |>
  distinct(subID, condition, AsianLatinoTherm, inMajorityTherm) |>
  mutate(outgroup_warmth = case_when(
    condition == "Minority" ~ AsianLatinoTherm,
    condition == "Majority" ~ inMajorityTherm,
    TRUE                    ~ NA_real_
  ))

ind_df2 <- left_join(ind_df, fulldf_warm[, c("subID","outgroup_warmth")], by = "subID")

if ("lambda" %in% names(ind_df2)) {
  ct_all <- cor.test(ind_df2$lambda, ind_df2$outgroup_warmth, use = "complete.obs")
  n_ct   <- sum(complete.cases(ind_df2[, c("lambda","outgroup_warmth")]))
  cat(sprintf("\nOverall r(λ, warmth): r = %.3f, p = %.4f, n = %d\n",
              ct_all$estimate, ct_all$p.value, n_ct))
  cat(sprintf("  95%% CI [%.3f, %.3f]\n", ct_all$conf.int[1], ct_all$conf.int[2]))

  for (cond in c("Minority", "Majority")) {
    sub <- filter(ind_df2, condition == cond)
    n_c <- sum(complete.cases(sub[, c("lambda","outgroup_warmth")]))
    if (n_c > 5) {
      ct_c <- cor.test(sub$lambda, sub$outgroup_warmth, use = "complete.obs")
      cat(sprintf("  %s: r = %.3f, p = %.4f, n = %d (CI [%.3f, %.3f])\n",
                  cond, ct_c$estimate, ct_c$p.value, n_c,
                  ct_c$conf.int[1], ct_c$conf.int[2]))
    }
  }
} else {
  # Note: AsianLatinoTherm or inMajorityTherm column names may differ.
  # Inspect fullTest_fixed.csv and adjust the mutate() above.
  cat("  NOTE: verify warmth column names in Study 3/Cleaning/output/fullTest_fixed.csv\n")
}

# ── 8. m/α by condition ────────────────────────────────────────────────────────
cat("\n=== Projection rate (α) by condition ===\n")
m_cols <- intersect(c("m", "m_in", "m_out"), names(ind_df))
for (mc in m_cols) {
  fit_m <- lm(reformulate("condition", response = mc), data = ind_df)
  s_m   <- summary(fit_m)
  f_m   <- s_m$fstatistic; p_m <- pf(f_m[1], f_m[2], f_m[3], lower.tail = FALSE)
  cat(sprintf("\n%s ~ condition: F(%d,%d) = %.3f, p = %.4f, R² = %.3f\n",
              mc, as.integer(f_m[2]), as.integer(f_m[3]), f_m[1], p_m, s_m$r.squared))
  print(round(coef(s_m)[, c("Estimate","Std. Error","t value","Pr(>|t|)")], 4))
}

# Paired m_in vs m_out if asym_lambda wins
if (winner == "asym_lambda" && all(c("m_in","m_out") %in% names(ind_df))) {
  cat("\n=== m_in vs m_out (paired, within-condition) ===\n")
  for (cond in levels(ind_df$condition)) {
    sub <- filter(ind_df, condition == cond)
    if (nrow(sub) > 5 && !all(is.na(sub$m_in))) {
      tt <- t.test(sub$m_in, sub$m_out, paired = TRUE)
      cat(sprintf("  %s: m_in = %.3f, m_out = %.3f | t(%d) = %.3f, p = %.4f\n",
                  cond, mean(sub$m_in,na.rm=T), mean(sub$m_out,na.rm=T),
                  tt$parameter, tt$statistic, tt$p.value))
    }
  }
}

# ── 9. MCR by condition ────────────────────────────────────────────────────────
if ("subject_mcr" %in% names(ind_df)) {
  cat("\n=== MCR by condition ===\n")
  fit_mcr <- lm(subject_mcr ~ condition, data = ind_df)
  s_mcr   <- summary(fit_mcr)
  f_mcr   <- s_mcr$fstatistic; p_mcr <- pf(f_mcr[1], f_mcr[2], f_mcr[3], lower.tail=FALSE)
  cat(sprintf("F(%d,%d) = %.3f, p = %.4f, R² = %.3f\n",
              as.integer(f_mcr[2]), as.integer(f_mcr[3]), f_mcr[1], p_mcr, s_mcr$r.squared))
  print(round(coef(s_mcr)[, c("Estimate","Std. Error","t value","Pr(>|t|)")], 4))

  mcr_means <- ind_df |>
    group_by(condition) |>
    summarise(M = round(mean(subject_mcr,na.rm=T),3),
              SD = round(sd(subject_mcr,na.rm=T),3), .groups="drop")
  cat("MCR by condition:\n"); print(mcr_means, row.names=FALSE)
}

# ── 10. Save results ────────────────────────────────────────────────────────────
write.csv(desc_df, here("Results","param_by_condition_s3_desc.csv"), row.names=FALSE)
write.csv(lm_df,   here("Results","param_by_condition_s3_lm.csv"),   row.names=FALSE)
write.csv(ind_df2, here("Results","ind_diffs_s3_full_enriched.csv"), row.names=FALSE)

message("\nSaved:")
message("  Results/param_by_condition_s3_desc.csv")
message("  Results/param_by_condition_s3_lm.csv")
message("  Results/ind_diffs_s3_full_enriched.csv")
