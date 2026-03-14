# race_moderation_s3.R
# Supplementary Analysis: Does participant race (Asian vs. Latino) moderate
# Study 3 effects? Tests condition × race interactions for all key outcomes.
#
# Rationale: Study 3 pools Asian and Latino participants, assuming the two groups
# are interchangeable within each contrast condition. This is not guaranteed:
#   - Minority condition: one group contrasts against the other; they play opposite
#     roles (ingroup vs. outgroup) within the same "condition" cell.
#   - Majority condition: both groups contrast against White, but societal
#     standing of Asians vs. Latinos may differ in ways that affect self-anchoring.
# Reviewers will likely probe whether findings replicate within each racial group.
#
# Race identifier: racid (Asian / Latino), from fullTest_fixed.csv
# Conditions: Minority (Asian vs Latino), Majority (Ingroup vs White)
#
# Analyses:
#   1. N by condition × race
#   2. Parameters by condition × race — lm(param ~ condition * racid)
#   3. Condition × race ANOVA for each parameter (test 3-way and 2-way effects)
#   4. Within-racid: condition effects (Minority vs Majority) for Asian and Latino separately
#   5. MCR and ELPD by condition × race
#   6. Warmth predictors × race interactions
#   7. Individual differences (MCR × personality) by race
#
# Prerequisites: param_comparison_conditions_s3.R, warmth_analysis_s3.R
#
# Output:
#   Results/race_moderation_params_s3.csv  — condition*race lm for each parameter
#   Results/race_moderation_corr_s3.csv    — MCR × personality by race

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# ── 1. Load data ───────────────────────────────────────────────────────────────
ind_df <- read.csv(here("Results", "ind_diffs_s3_full_enriched.csv")) |>
  mutate(condition = factor(condition, levels = c("Minority", "Majority")))

fulldf <- read.csv(here("Study 3/Cleaning/output/fullTest_fixed.csv")) |>
  filter(!is.na(ingChoiceN)) |>
  distinct(subID, racid, inMinorityTherm, inMajorityTherm)

df <- left_join(ind_df, fulldf, by = "subID") |>
  mutate(racid = factor(racid, levels = c("Asian", "Latino")))

cat("N by condition × race:\n")
print(table(df$condition, df$racid))

param_cols <- intersect(c("m","bias","lambda","w","subject_mcr","delta_elpd"), names(df))

# ── 2. Parameter descriptives by condition × race ──────────────────────────────
cat("\n=== Parameter descriptives by condition × race ===\n")
desc_2way <- df |>
  group_by(condition, racid) |>
  summarise(across(all_of(param_cols),
                   list(M = ~round(mean(.x, na.rm=T), 3),
                        SD = ~round(sd(.x, na.rm=T), 3)),
                   .names = "{.col}_{.fn}"),
            n = n(), .groups = "drop")
print(desc_2way, width = 120)

# ── 3. lm(param ~ condition * racid) for each parameter ───────────────────────
cat("\n=== lm(param ~ condition * racid) — interaction tests ===\n")
lm_results <- list()
for (p in param_cols) {
  sub <- df[!is.na(df[[p]]), ]
  mod_int  <- lm(reformulate(c("condition","racid","condition:racid"), p), sub)
  mod_cond <- lm(reformulate(c("condition","racid"), p), sub)
  av <- anova(mod_cond, mod_int)
  int_F <- av$F[2]; int_p <- av$`Pr(>F)`[2]; int_df <- av$Res.Df[2]

  cat(sprintf("\n── %s ──\n", p))
  cat(sprintf("  Condition × Race interaction: F(1,%d) = %.3f, p = %.4f\n",
              int_df, int_F, int_p))
  coefs <- round(summary(mod_int)$coefficients, 4)
  print(coefs)
  cat(sprintf("  Main effect condition: F = %.3f, p = %.4f\n",
              summary(mod_cond)$fstatistic[1],
              pf(summary(mod_cond)$fstatistic[1],
                 summary(mod_cond)$fstatistic[2],
                 summary(mod_cond)$fstatistic[3], lower.tail=FALSE)))

  res <- as.data.frame(coefs)
  res$term <- rownames(res); res$param <- p
  res$int_F <- int_F; res$int_p <- int_p
  lm_results[[p]] <- res
}
lm_df <- do.call(rbind, lm_results)

# ── 4. Within-race condition effects ──────────────────────────────────────────
cat("\n=== Within-race: condition effect (Minority vs Majority) ===\n")
for (race in c("Asian","Latino")) {
  cat(sprintf("\n──── %s participants ────\n", race))
  sub_r <- filter(df, racid == race)
  for (p in param_cols) {
    sub2 <- sub_r[!is.na(sub_r[[p]]), ]
    if (nrow(sub2) < 10) next
    tt <- t.test(reformulate("condition", response = p), data = sub2)
    means <- tapply(sub2[[p]], sub2$condition, mean, na.rm=TRUE)
    cat(sprintf("  %s: Minority M=%.3f, Majority M=%.3f | t(%d)=%.3f, p=%.4f\n",
                p, means["Minority"], means["Majority"],
                round(tt$parameter), tt$statistic, tt$p.value))
  }
}

# ── 5. ELPD condition × race ───────────────────────────────────────────────────
if ("delta_elpd" %in% names(df)) {
  cat("\n=== ELPD model advantage by condition × race ===\n")
  # delta_elpd = asym − sym_lambda per subject (from loo_delta_indiff)
  desc_elpd <- df |>
    group_by(condition, racid) |>
    summarise(mean_delta = round(mean(delta_elpd, na.rm=T), 4),
              sd_delta = round(sd(delta_elpd, na.rm=T), 4),
              n = n(), .groups="drop")
  print(desc_elpd)
}

# ── 6. Warmth × race interactions ─────────────────────────────────────────────
cat("\n=== Warmth × Race interactions for key parameters ===\n")
warmth_vars <- c("inMinorityTherm","inMajorityTherm")
for (p in c("bias","lambda","subject_mcr")) {
  for (wv in warmth_vars) {
    sub <- df[!is.na(df[[p]]) & !is.na(df[[wv]]), ]
    if (nrow(sub) < 20) next
    mod_int  <- lm(reformulate(c(wv,"racid",paste0(wv,":racid")), p), sub)
    mod_main <- lm(reformulate(c(wv,"racid"), p), sub)
    av <- anova(mod_main, mod_int)
    int_F <- av$F[2]; int_p <- av$`Pr(>F)`[2]; int_df <- av$Res.Df[2]
    if (int_p < .10) {
      cat(sprintf("\n  %s ~ %s * race: F(1,%d)=%.3f, p=%.4f\n",
                  p, wv, int_df, int_F, int_p))
      print(round(summary(mod_int)$coefficients, 4))
    }
  }
}
cat("  (Only showing p < .10; remaining interactions null)\n")

# ── 7. MCR × personality correlations by race ─────────────────────────────────
cat("\n=== MCR × personality scales by race ===\n")
scale_vars <- intersect(c("DS","Proto","SCC","SI","RSE","NTB","NFC","SING.Ind","SING.Inter"),
                        names(df))
cat(sprintf("Scales available: %s\n", paste(scale_vars, collapse=", ")))

cor_race <- list()
for (race in c("Asian","Latino")) {
  sub_r <- filter(df, racid == race)
  cat(sprintf("\n──── %s (N=%d) ────\n", race, nrow(sub_r)))
  for (sv in scale_vars) {
    sub2 <- sub_r[!is.na(sub_r$subject_mcr) & !is.na(sub_r[[sv]]), ]
    if (nrow(sub2) < 5) next
    ct <- cor.test(sub2$subject_mcr, sub2[[sv]])
    cor_race[[length(cor_race)+1]] <- data.frame(
      race = race, scale = sv,
      r = round(unname(ct$estimate), 3),
      p = round(ct$p.value, 4), n = nrow(sub2)
    )
    if (ct$p.value < 0.05) {
      cat(sprintf("  MCR × %s: r=%.3f, p=%.4f, n=%d\n", sv,
                  unname(ct$estimate), ct$p.value, nrow(sub2)))
    }
  }
}
cor_race_df <- do.call(rbind, cor_race)
cat("\nAll MCR × scale correlations:\n")
print(cor_race_df, row.names=FALSE)

# Fisher z-test for race moderation of MCR × scale
cat("\n=== Fisher z: does race moderate MCR correlations? ===\n")
for (sv in scale_vars) {
  r_asian  <- filter(cor_race_df, race=="Asian",  scale==sv)$r
  r_latino <- filter(cor_race_df, race=="Latino", scale==sv)$r
  n_asian  <- filter(cor_race_df, race=="Asian",  scale==sv)$n
  n_latino <- filter(cor_race_df, race=="Latino", scale==sv)$n
  if (length(r_asian)==0 || length(r_latino)==0) next
  z1 <- 0.5*log((1+r_asian)/(1-r_asian))
  z2 <- 0.5*log((1+r_latino)/(1-r_latino))
  se <- sqrt(1/(n_asian-3) + 1/(n_latino-3))
  z_diff <- (z1-z2)/se
  p_diff <- 2*pnorm(-abs(z_diff))
  if (p_diff < .10) {
    cat(sprintf("  %s: Asian r=%.3f, Latino r=%.3f | Fisher z=%.3f, p=%.4f\n",
                sv, r_asian, r_latino, z_diff, p_diff))
  }
}
cat("  (Only showing p < .10; remaining null)\n")

# ── 8. Condition main effect replication within race ──────────────────────────
cat("\n=== Summary: condition × race cell means for key params ===\n")
key_params <- intersect(c("bias","lambda","w","subject_mcr"), names(df))
cell_means <- df |>
  group_by(condition, racid) |>
  summarise(across(all_of(key_params),
                   ~round(mean(.x, na.rm=T), 3)),
            n = n(), .groups="drop")
print(cell_means, width=120)

# ── 9. Save ───────────────────────────────────────────────────────────────────
write.csv(lm_df, here("Results","race_moderation_params_s3.csv"), row.names=FALSE)
write.csv(cor_race_df, here("Results","race_moderation_corr_s3.csv"), row.names=FALSE)

message("\nSaved:")
message("  Results/race_moderation_params_s3.csv")
message("  Results/race_moderation_corr_s3.csv")
