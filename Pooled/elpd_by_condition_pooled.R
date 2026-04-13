# elpd_by_condition_pooled.R
# Per-study / per-condition LOO-ELPD comparison: sym_lambda vs asym_lambda --- Pooled model
#
# Key question: Does asym_lambda show a fit advantage in any particular study
# or condition, even when the pooled sym_lambda wins overall?
#
# The pooled LOO log-likelihood flat vector uses stride-based indexing:
#   idx_i = ((i-1)*maxTrials + 1) : ((i-1)*maxTrials + nTrials[i])
# Subjects are ordered by sorted uniqueSubID = "study_subID".
#
# Prerequisites:
#   Fits/loo_pooled_sym_lambda.rds
#   Fits/loo_pooled_asym_lambda.rds
#   Fits/pooled_stan_data.rds  (contains nTrials, maxTrials, subjStudy)
#
# Condition info reconstructed from the same source files as prepare_pooled_data.R:
#   S1: all subjects coded as "MinGrp" (single condition)
#   S2: outgroup column -> Negation / High-Status / Low-Status
#   S3: condition column -> Minority / Majority
#
# Output:
#   Results/elpd_by_condition_pooled_subj.csv    --- per-subject
#   Results/elpd_by_condition_pooled_summary.csv --- per-study totals
#   Results/elpd_by_condition_pooled_delta.csv   --- per-study x condition t-tests

suppressPackageStartupMessages({
  library(tidyverse)
  library(loo)
  library(here)
})

# -- 1. Load pooled stan_data -------------------------------------------------
stan_data <- readRDS(here("Fits", "pooled_stan_data.rds"))

nSubjects   <- stan_data$nSubjects
maxTrials   <- max(stan_data$maxTrials)
nTrials_vec <- stan_data$nTrials
subjStudy   <- stan_data$subjStudy   # integer 1/2/3 per subject

cat(sprintf("Total pooled subjects: %d\n", nSubjects))
cat(sprintf("maxTrials: %d\n", maxTrials))
cat("\nSubjects by study:\n")
print(table(subjStudy))

# -- 2. Build condition lookup per study --------------------------------------

# S1: single condition
s1_test  <- read.csv(here("Study 1/Cleaning/output/fullTest.csv"))  |> filter(!is.na(ingChoiceN))
s1_train <- read.csv(here("Study 1/Cleaning/output/fullTrain.csv")) |> filter(!is.na(selfResp))
common_s1 <- sort(intersect(unique(s1_test$subID), unique(s1_train$subID)))
cond_s1 <- data.frame(
  uniqueSubID = paste0("1_", common_s1),
  study_label = "S1",
  condition   = "MinGrp"
)

# S2: condition from outgroup column
s2_test  <- read.csv(here("Study 2/Cleaning/output/fullTest.csv"))       |> filter(!is.na(ingChoiceN))
s2_train <- read.csv(here("Study 2/Cleaning/output/fullTrain_fixed.csv")) |> filter(!is.na(selfResp))
common_s2 <- sort(intersect(unique(s2_test$subID), unique(s2_train$subID)))
cond_s2 <- s2_test |>
  filter(subID %in% common_s2) |>
  distinct(subID, outgroup) |>
  mutate(
    uniqueSubID = paste0("2_", subID),
    study_label = "S2",
    condition   = recode(outgroup,
                         "Not UCR" = "Negation",
                         "UCLA"    = "High-Status",
                         "CSU LA"  = "Low-Status")
  ) |>
  select(uniqueSubID, study_label, condition)

# S3: condition column directly
s3_test  <- read.csv(here("Study 3/Cleaning/output/fullTest_fixed.csv"))  |> filter(!is.na(ingChoiceN))
s3_train <- read.csv(here("Study 3/Cleaning/output/fullTrain_fixed.csv")) |> filter(!is.na(selfResp))
common_s3 <- sort(intersect(unique(s3_test$subID), unique(s3_train$subID)))
cond_s3 <- s3_test |>
  filter(subID %in% common_s3) |>
  distinct(subID, condition) |>
  mutate(
    uniqueSubID = paste0("3_", subID),
    study_label = "S3"
  ) |>
  select(uniqueSubID, study_label, condition)

cond_all <- bind_rows(cond_s1, cond_s2, cond_s3)

# -- 3. Build subject metadata table ------------------------------------------
# Reconstruct uIds in the same order as prepare_pooled_data.R:
#   sort(unique(c("1_subID", "2_subID", "3_subID")))
all_ids <- c(
  paste0("1_", common_s1),
  paste0("2_", common_s2),
  paste0("3_", common_s3)
)
uIds <- sort(unique(all_ids))

stopifnot(length(uIds) == nSubjects)   # sanity check

subj_meta <- data.frame(
  subj_idx    = seq_len(nSubjects),
  uniqueSubID = uIds,
  study_num   = subjStudy,
  n_trials    = nTrials_vec
) |>
  left_join(cond_all, by = "uniqueSubID")

cat("\nN by study x condition:\n")
print(table(subj_meta$study_label, subj_meta$condition, useNA = "ifany"))

# -- 4. Extract per-subject ELPD (stride indexing) ----------------------------
get_subj_elpd_stride <- function(loo_file, nSubj, nTrials_vec, maxTrials) {
  l  <- readRDS(loo_file)
  pw <- l$pointwise[, "elpd_loo"]
  cat(sprintf("  LOO pointwise length: %d  (expected: %d)\n",
              length(pw), nSubj * maxTrials))
  sapply(seq_len(nSubj), function(i) {
    idx <- ((i - 1) * maxTrials + 1):((i - 1) * maxTrials + nTrials_vec[i])
    sum(pw[idx])
  })
}

cat("\n-- Loading sym_lambda LOO --\n")
elpd_sym  <- get_subj_elpd_stride(here("Fits", "loo_pooled_sym_lambda.rds"),
                                   nSubjects, nTrials_vec, maxTrials)
cat("-- Loading asym_lambda LOO --\n")
elpd_asym <- get_subj_elpd_stride(here("Fits", "loo_pooled_asym_lambda.rds"),
                                   nSubjects, nTrials_vec, maxTrials)

# -- 5. Verify totals match loo_compare ---------------------------------------
cat("\n=== Overall ELPD totals (should match loo_compare) ===\n")
cat(sprintf("  sym_lambda:  %.2f\n", sum(elpd_sym)))
cat(sprintf("  asym_lambda: %.2f\n", sum(elpd_asym)))
cat(sprintf("  DELTA_ELPD (asym - sym): %.2f\n", sum(elpd_asym) - sum(elpd_sym)))

# -- 6. Per-subject delta -----------------------------------------------------
subj_elpd <- subj_meta |>
  mutate(
    elpd_sym_lambda  = elpd_sym,
    elpd_asym_lambda = elpd_asym,
    delta_asym_sym   = elpd_asym - elpd_sym
  )

# -- 7. Per-study x condition summary -----------------------------------------
cond_totals <- subj_elpd |>
  group_by(study_label, condition) |>
  summarise(
    n               = n(),
    elpd_sym        = sum(elpd_sym_lambda,  na.rm = TRUE),
    elpd_asym       = sum(elpd_asym_lambda, na.rm = TRUE),
    delta_total     = sum(delta_asym_sym,   na.rm = TRUE),
    mean_delta      = round(mean(delta_asym_sym, na.rm = TRUE), 4),
    sd_delta        = round(sd(delta_asym_sym,   na.rm = TRUE), 4),
    pct_favors_asym = round(100 * mean(delta_asym_sym > 0, na.rm = TRUE), 1),
    t_stat          = {tt <- t.test(delta_asym_sym, mu = 0); round(tt$statistic, 3)},
    df_t            = {tt <- t.test(delta_asym_sym, mu = 0); round(tt$parameter, 1)},
    p_val           = {tt <- t.test(delta_asym_sym, mu = 0); round(tt$p.value, 4)},
    .groups         = "drop"
  )

cat("\n=== Per-study x condition ELPD results ===\n")
print(cond_totals, digits = 4)

# -- 8. Per-study totals -------------------------------------------------------
study_totals <- subj_elpd |>
  group_by(study_label) |>
  summarise(
    n               = n(),
    elpd_sym        = sum(elpd_sym_lambda,  na.rm = TRUE),
    elpd_asym       = sum(elpd_asym_lambda, na.rm = TRUE),
    delta_total     = sum(delta_asym_sym,   na.rm = TRUE),
    mean_delta      = round(mean(delta_asym_sym, na.rm = TRUE), 4),
    sd_delta        = round(sd(delta_asym_sym,   na.rm = TRUE), 4),
    pct_favors_asym = round(100 * mean(delta_asym_sym > 0, na.rm = TRUE), 1),
    t_stat          = {tt <- t.test(delta_asym_sym, mu = 0); round(tt$statistic, 3)},
    df_t            = {tt <- t.test(delta_asym_sym, mu = 0); round(tt$parameter, 1)},
    p_val           = {tt <- t.test(delta_asym_sym, mu = 0); round(tt$p.value, 4)},
    .groups         = "drop"
  )

cat("\n=== Per-study ELPD totals ===\n")
print(study_totals, digits = 4)

# Within-S2: pairwise between-condition comparisons
s2_deltas <- filter(subj_elpd, study_label == "S2")
s2_conds  <- sort(unique(s2_deltas$condition))
cat("\n=== S2: Between-condition pairwise comparisons of DELTA_ELPD ===\n")
for (ca in s2_conds) {
  for (cb in s2_conds) {
    if (ca < cb) {
      da <- filter(s2_deltas, condition == ca)$delta_asym_sym
      db <- filter(s2_deltas, condition == cb)$delta_asym_sym
      tt <- t.test(da, db)
      cat(sprintf("  %s vs %s: DELTA_DELTA=%.3f, t(%d)=%.3f, p=%.4f\n",
                  ca, cb,
                  mean(da, na.rm=TRUE) - mean(db, na.rm=TRUE),
                  round(tt$parameter), tt$statistic, tt$p.value))
    }
  }
}

# Within-S3: between-condition comparison
s3_deltas <- filter(subj_elpd, study_label == "S3")
cat("\n=== S3: Between-condition comparison of DELTA_ELPD ===\n")
d_minor <- filter(s3_deltas, condition == "Minority")$delta_asym_sym
d_major <- filter(s3_deltas, condition == "Majority")$delta_asym_sym
tt_s3   <- t.test(d_minor, d_major)
cat(sprintf("  Minority vs Majority: DELTA_DELTA=%.3f, t(%d)=%.3f, p=%.4f\n",
            mean(d_minor, na.rm=TRUE) - mean(d_major, na.rm=TRUE),
            round(tt_s3$parameter), tt_s3$statistic, tt_s3$p.value))

# -- 9. Save outputs -----------------------------------------------------------
write.csv(subj_elpd,    here("Results", "elpd_by_condition_pooled_subj.csv"),    row.names = FALSE)
write.csv(cond_totals,  here("Results", "elpd_by_condition_pooled_delta.csv"),   row.names = FALSE)
write.csv(study_totals, here("Results", "elpd_by_condition_pooled_summary.csv"), row.names = FALSE)

message("\nSaved:")
message("  Results/elpd_by_condition_pooled_subj.csv")
message("  Results/elpd_by_condition_pooled_delta.csv")
message("  Results/elpd_by_condition_pooled_summary.csv")
