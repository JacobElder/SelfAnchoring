# elpd_by_condition_s2.R
# Per-condition LOO-ELPD comparison across all 4 models — Study 2
#
# Key question: Does asym_lambda show a differential fit advantage in one
# condition over others (e.g., better in High-Status where the outgroup is
# positively-regarded, potentially triggering asymmetric in/out projection)?
#
# Uses stride-based indexing: idx = ((i-1)*maxTrials + 1):(i-1)*maxTrials + nTrials[i]
# Consistent with loo_delta_indiff_s2.R; avoids padding artifacts.
#
# Prerequisites: run_model_comparison_s2_remaining.R (all 4 LOO files saved)
#
# Output:
#   Results/elpd_by_condition_s2_subj.csv     — per-subject ELPD × 4 models + condition
#   Results/elpd_by_condition_s2_summary.csv  — per-condition model totals + ΔELPD
#   Results/elpd_by_condition_s2_delta.csv    — per-condition asym vs sym t-test

suppressPackageStartupMessages({
  library(tidyverse)
  library(loo)
  library(here)
})

# ── 1. Reconstruct subject order and trial counts ─────────────────────────────
fulldf  <- read.csv(here("Study 2/Cleaning/output/fullTest.csv"))  |> filter(!is.na(ingChoiceN))
traindf <- read.csv(here("Study 2/Cleaning/output/fullTrain_fixed.csv")) |> filter(!is.na(selfResp))

common_ids  <- sort(intersect(unique(fulldf$subID), unique(traindf$subID)))
fulldf      <- filter(fulldf, subID %in% common_ids)
uIds        <- sort(common_ids)
nSubjects   <- length(uIds)
maxTrials   <- max(fulldf$trialTotalT2)
nTrials_vec <- sapply(uIds, function(id) nrow(filter(fulldf, subID == id)))

cat(sprintf("Subjects: %d, maxTrials: %d\n", nSubjects, maxTrials))

# ── 2. Condition per subject ───────────────────────────────────────────────────
cond_df <- fulldf |>
  distinct(subID, outgroup) |>
  mutate(condition = factor(outgroup,
                            levels = c("Not UCR", "UCLA", "CSU LA"),
                            labels = c("Negation", "High-Status", "Low-Status")))

subj_meta <- data.frame(subj_idx = seq_along(uIds), subID = uIds,
                        n_trials = nTrials_vec) |>
  left_join(cond_df, by = "subID")

cat("\nN by condition:\n"); print(table(subj_meta$condition))

# ── 3. Extract per-subject ELPD for all 4 models (stride indexing) ────────────
loo_files <- list(
  bias        = here("Fits", "loo_s2_bias.rds"),
  symmetric   = here("Fits", "loo_s2_symmetric.rds"),
  sym_lambda  = here("Fits", "loo_s2_sym_lambda.rds"),
  asym_lambda = here("Fits", "loo_s2_asym_lambda.rds")
)

get_subj_elpd_stride <- function(loo_file, nSubj, nTrials_vec, maxTrials) {
  l  <- readRDS(loo_file)
  pw <- l$pointwise[, "elpd_loo"]
  sapply(seq_len(nSubj), function(i) {
    idx <- ((i - 1) * maxTrials + 1):((i - 1) * maxTrials + nTrials_vec[i])
    sum(pw[idx])
  })
}

elpd_mat <- vapply(names(loo_files), function(nm) {
  f <- loo_files[[nm]]
  if (!file.exists(f)) {
    warning("LOO file not found: ", f)
    return(rep(NA_real_, nSubjects))
  }
  get_subj_elpd_stride(f, nSubjects, nTrials_vec, maxTrials)
}, numeric(nSubjects))

subj_elpd <- cbind(subj_meta, as.data.frame(elpd_mat))

# ── 4. Verify overall totals match loo_compare ────────────────────────────────
cat("\n=== Overall ELPD totals (should match loo_compare) ===\n")
overall <- sort(colSums(elpd_mat, na.rm = TRUE), decreasing = TRUE)
print(round(overall, 2))

# ── 5. Per-condition model totals ────────────────────────────────────────────
model_labels <- c(bias="Bias", symmetric="Symmetric",
                  sym_lambda="Sym+\u03bb", asym_lambda="Asym+\u03bb")

cond_totals <- subj_elpd |>
  group_by(condition) |>
  summarise(n = n(),
            across(all_of(names(loo_files)), ~sum(.x, na.rm = TRUE)),
            .groups = "drop")

cond_long <- cond_totals |>
  pivot_longer(cols = all_of(names(loo_files)),
               names_to = "model", values_to = "elpd_total") |>
  group_by(condition) |>
  mutate(
    delta_elpd  = elpd_total - max(elpd_total, na.rm = TRUE),
    rank        = rank(-elpd_total),
    model_label = model_labels[model]
  ) |>
  ungroup()

cat("\n=== Model ranking within each condition ===\n")
for (cond in levels(subj_elpd$condition)) {
  n_c <- filter(subj_elpd, condition == cond) |> nrow()
  cat(sprintf("\n── %s (N = %d) ──\n", cond, n_c))
  sub <- filter(cond_long, condition == cond) |>
    arrange(rank) |>
    select(model_label, elpd_total, delta_elpd)
  print(sub, digits = 3)
}

# ── 6. Per-condition asym_lambda vs sym_lambda advantage ──────────────────────
cat("\n=== Per-subject ΔELPD (asym_lambda − sym_lambda) by condition ===\n")
subj_elpd$delta_asym_sym <- subj_elpd$asym_lambda - subj_elpd$sym_lambda

delta_by_cond <- subj_elpd |>
  group_by(condition) |>
  summarise(
    n            = n(),
    mean_delta   = round(mean(delta_asym_sym, na.rm = TRUE), 4),
    sd_delta     = round(sd(delta_asym_sym, na.rm = TRUE), 4),
    pct_favors_asym = round(100 * mean(delta_asym_sym > 0, na.rm = TRUE), 1),
    t_stat       = {tt <- t.test(delta_asym_sym, mu = 0); round(tt$statistic, 3)},
    df_t         = {tt <- t.test(delta_asym_sym, mu = 0); round(tt$parameter, 1)},
    p_val        = {tt <- t.test(delta_asym_sym, mu = 0); round(tt$p.value, 4)},
    .groups      = "drop"
  )

cat("\nMean ΔELPD (asym − sym_lambda) per condition:\n")
print(delta_by_cond, digits = 3)

cat("\nInterpretation:\n")
cat("  + = asym_lambda fits better in that condition\n")
cat("  - = sym_lambda fits better (simpler model sufficient)\n")
cat("  Interest: Is High-Status different from Negation/Low-Status?\n")

# Pairwise condition comparison of deltas (does the advantage differ by condition?)
if (length(levels(subj_elpd$condition)) == 3) {
  cat("\n── Pairwise t-tests: ΔELPD differs by condition? ──\n")
  conds <- levels(subj_elpd$condition)
  for (i in 1:(length(conds) - 1)) {
    for (j in (i+1):length(conds)) {
      d1 <- filter(subj_elpd, condition == conds[i])$delta_asym_sym
      d2 <- filter(subj_elpd, condition == conds[j])$delta_asym_sym
      tt <- t.test(d1, d2)
      cat(sprintf("  %s vs %s: ΔΔELPD = %.3f, t(%d) = %.3f, p = %.4f\n",
                  conds[i], conds[j],
                  mean(d1, na.rm=TRUE) - mean(d2, na.rm=TRUE),
                  round(tt$parameter), tt$statistic, tt$p.value))
    }
  }
}

# ── 7. Save ───────────────────────────────────────────────────────────────────
write.csv(subj_elpd,   here("Results", "elpd_by_condition_s2_subj.csv"),    row.names = FALSE)
write.csv(cond_long,   here("Results", "elpd_by_condition_s2_summary.csv"), row.names = FALSE)
write.csv(delta_by_cond, here("Results", "elpd_by_condition_s2_delta.csv"), row.names = FALSE)

message("\nSaved:")
message("  Results/elpd_by_condition_s2_subj.csv")
message("  Results/elpd_by_condition_s2_summary.csv")
message("  Results/elpd_by_condition_s2_delta.csv")
