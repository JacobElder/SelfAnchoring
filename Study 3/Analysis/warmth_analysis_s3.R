# warmth_analysis_s3.R
# Intergroup preference correlates of model parameters — Study 3
#
# Thermometer variables are all DIFFERENCE SCORES (not absolute warmth):
#   Therm_1 = feeling toward Whites
#   Therm_2 = feeling toward Asians
#   Therm_4 = feeling toward Latinos
#
#   AsianLatinoTherm  = Therm_2 − Therm_4  (Asian − Latino, group-level)
#   WhiteAsianTherm   = Therm_1 − Therm_2  (White − Asian)
#   WhiteLatinoTherm  = Therm_1 − Therm_4  (White − Latino)
#
#   inMinorityTherm:
#     Asian participants → Therm_2 − Therm_4 (Asian − Latino)
#     Latino participants → Therm_4 − Therm_2 (Latino − Asian)
#     = INGROUP-OVER-OTHER-MINORITY preference (positive = prefer own over other minority)
#
#   inMajorityTherm:
#     Asian participants → Therm_2 − Therm_1 (Asian − White)
#     Latino participants → Therm_4 − Therm_1 (Latino − White)
#     = INGROUP-OVER-MAJORITY preference (positive = prefer own minority over Whites;
#       negative = prefer Whites over own minority)
#
# Key constructs:
#   inMinorityTherm  = Minority-condition ingroup preference (own vs other-minority)
#   inMajorityTherm  = Cross-condition ingroup preference (own minority vs Whites; available both conditions)
#   AsianLatinoTherm = Group-level Asian−Latino difference (NOT person-specific — avoid for correlations)
#
# Analyses:
#   1. Descriptives by condition
#   2. r(param, inMinorityTherm) — person-specific ingroup preference over other-minority
#      Relevant primarily for Minority condition; analogous to Study 2 r(λ, warmth)
#   3. r(param, inMajorityTherm) — person-specific preference over White majority
#      Both conditions; tests whether anti-White affect predicts params
#   4. lm(param ~ inMajorityTherm * condition) — condition moderation
#   5. Within-condition focused tests for Minority condition
#
# Prerequisites: param_comparison_conditions_s3.R (produces ind_diffs_s3_full_enriched.csv)
#
# Output:
#   Results/warmth_correlations_s3.csv
#   Results/warmth_lm_s3.csv

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# ── 1. Load data ───────────────────────────────────────────────────────────────
ind_df <- read.csv(here("Results", "ind_diffs_s3_full_enriched.csv")) |>
  mutate(condition = factor(condition, levels = c("Minority", "Majority")))

fulldf_raw <- read.csv(here("Study 3/Cleaning/output/fullTest_fixed.csv")) |>
  filter(!is.na(ingChoiceN)) |>
  distinct(subID, condition, inMinorityTherm, inMajorityTherm, AsianLatinoTherm,
           WhiteAsianTherm, WhiteLatinoTherm)

df <- left_join(ind_df, fulldf_raw, by = c("subID", "condition"))

cat(sprintf("N = %d; Minority = %d, Majority = %d\n",
            nrow(df),
            sum(df$condition == "Minority"),
            sum(df$condition == "Majority")))

# ── 2. Descriptives by condition ───────────────────────────────────────────────
cat("\n=== Thermometer descriptives by condition ===\n")
cat("Reminder: all scores are person-specific DIFFERENCE scores (ingroup − outgroup or ingroup − majority)\n")
cat("  inMinorityTherm  = own-minority preference over other-minority (positive = prefer own)\n")
cat("  inMajorityTherm  = own-minority preference over White majority  (positive = prefer own over Whites)\n\n")

df |>
  group_by(condition) |>
  summarise(
    inMinority_M  = round(mean(inMinorityTherm, na.rm=T), 1),
    inMinority_SD = round(sd(inMinorityTherm,  na.rm=T), 1),
    inMinority_n  = sum(!is.na(inMinorityTherm)),
    inMajority_M  = round(mean(inMajorityTherm, na.rm=T), 1),
    inMajority_SD = round(sd(inMajorityTherm,  na.rm=T), 1),
    inMajority_n  = sum(!is.na(inMajorityTherm)),
    .groups = "drop"
  ) |> print()

# ── 3. Helper ──────────────────────────────────────────────────────────────────
do_cor <- function(x_var, y_var, data, label = NULL) {
  sub <- data[!is.na(data[[x_var]]) & !is.na(data[[y_var]]), ]
  if (nrow(sub) < 5) return(NULL)
  ct <- cor.test(sub[[x_var]], sub[[y_var]])
  lbl <- if (is.null(label)) y_var else label
  data.frame(
    param      = x_var,
    thermometer= lbl,
    condition  = "Overall",
    r   = round(unname(ct$estimate), 3),
    p   = round(ct$p.value, 4),
    n   = nrow(sub),
    ci_lo = round(ct$conf.int[1], 3),
    ci_hi = round(ct$conf.int[2], 3)
  )
}

do_cor_cond <- function(x_var, y_var, data, cond, label = NULL) {
  sub <- filter(data, condition == cond)
  sub <- sub[!is.na(sub[[x_var]]) & !is.na(sub[[y_var]]), ]
  if (nrow(sub) < 5) return(NULL)
  ct <- cor.test(sub[[x_var]], sub[[y_var]])
  lbl <- if (is.null(label)) y_var else label
  data.frame(
    param      = x_var,
    thermometer= lbl,
    condition  = cond,
    r   = round(unname(ct$estimate), 3),
    p   = round(ct$p.value, 4),
    n   = nrow(sub),
    ci_lo = round(ct$conf.int[1], 3),
    ci_hi = round(ct$conf.int[2], 3)
  )
}

fisher_z_test <- function(r1, r2, n1, n2) {
  z1 <- 0.5 * log((1+r1)/(1-r1))
  z2 <- 0.5 * log((1+r2)/(1-r2))
  se  <- sqrt(1/(n1-3) + 1/(n2-3))
  z_diff <- (z1-z2)/se
  list(z = z_diff, p = 2*pnorm(-abs(z_diff)))
}

param_cols <- intersect(c("m","bias","lambda","w","subject_mcr","delta_elpd"), names(df))

# ── 4. Correlations: param × inMinorityTherm (overall and by condition) ────────
cat("\n=== r(param, inMinorityTherm) — own-minority over other-minority preference ===\n")
rows_min <- list()
for (p in param_cols) {
  rows_min[[length(rows_min)+1]] <- do_cor(p, "inMinorityTherm", df, "inMinority")
  for (cond in c("Minority","Majority")) {
    rows_min[[length(rows_min)+1]] <- do_cor_cond(p, "inMinorityTherm", df, cond, "inMinority")
  }
}
cor_min <- do.call(rbind, Filter(Negate(is.null), rows_min))
print(cor_min, row.names=FALSE, digits=3)

# ── 5. Correlations: param × inMajorityTherm (overall and by condition) ────────
cat("\n=== r(param, inMajorityTherm) — own-minority over White majority preference ===\n")
rows_maj <- list()
for (p in param_cols) {
  rows_maj[[length(rows_maj)+1]] <- do_cor(p, "inMajorityTherm", df, "inMajority")
  for (cond in c("Minority","Majority")) {
    rows_maj[[length(rows_maj)+1]] <- do_cor_cond(p, "inMajorityTherm", df, cond, "inMajority")
  }
}
cor_maj <- do.call(rbind, Filter(Negate(is.null), rows_maj))
print(cor_maj, row.names=FALSE, digits=3)

# ── 6. Combined and nominally significant ─────────────────────────────────────
cor_all <- rbind(cor_min, cor_maj)
write.csv(cor_all, here("Results","warmth_correlations_s3.csv"), row.names=FALSE)

cat("\n=== Nominally significant (p < .05) ===\n")
print(filter(cor_all, p < 0.05) |> arrange(p), row.names=FALSE)

# ── 7. Moderation: lm(param ~ inMinorityTherm * condition) ────────────────────
cat("\n\n=== Moderation: param ~ inMinorityTherm * condition ===\n")
lm_results <- list()
for (p in param_cols) {
  sub <- df[!is.na(df[[p]]) & !is.na(df$inMinorityTherm), ]
  if (nrow(sub) < 20) next
  mod_int  <- lm(reformulate(c("inMinorityTherm","condition","inMinorityTherm:condition"), p), sub)
  mod_main <- lm(reformulate(c("inMinorityTherm","condition"), p), sub)
  av <- anova(mod_main, mod_int)
  int_p <- av$`Pr(>F)`[2]; int_F <- av$F[2]; int_df <- av$Res.Df[2]
  cat(sprintf("\n  %s ~ inMinority * condition: F(1,%d) = %.3f, p = %.4f\n",
              p, int_df, int_F, int_p))
  if (int_p < .10) {
    print(round(summary(mod_int)$coefficients, 4))
  }
  coefs <- as.data.frame(summary(mod_int)$coefficients)
  coefs$term <- rownames(coefs); coefs$param <- p
  coefs$interaction_F <- int_F; coefs$interaction_p <- int_p
  lm_results[[p]] <- coefs
}

# ── 8. Moderation: lm(param ~ inMajorityTherm * condition) ────────────────────
cat("\n\n=== Moderation: param ~ inMajorityTherm * condition ===\n")
for (p in param_cols) {
  sub <- df[!is.na(df[[p]]) & !is.na(df$inMajorityTherm), ]
  if (nrow(sub) < 20) next
  mod_int  <- lm(reformulate(c("inMajorityTherm","condition","inMajorityTherm:condition"), p), sub)
  mod_main <- lm(reformulate(c("inMajorityTherm","condition"), p), sub)
  av <- anova(mod_main, mod_int)
  int_p <- av$`Pr(>F)`[2]; int_F <- av$F[2]; int_df <- av$Res.Df[2]
  cat(sprintf("\n  %s ~ inMajority * condition: F(1,%d) = %.3f, p = %.4f\n",
              p, int_df, int_F, int_p))
  if (int_p < .10) {
    print(round(summary(mod_int)$coefficients, 4))
  }
  coefs <- as.data.frame(summary(mod_int)$coefficients)
  coefs$term <- rownames(coefs); coefs$param <- p
  coefs$interaction_F <- int_F; coefs$interaction_p <- int_p
  lm_results[[paste0(p,"_inMaj")]] <- coefs
}

lm_df <- do.call(rbind, Filter(Negate(is.null), lm_results))

# ── 9. KEY: γ and MCR focused tests ───────────────────────────────────────────
cat("\n\n=== KEY: γ (ingroup bias) correlations ===\n")

for (tvar in c("inMinorityTherm","inMajorityTherm")) {
  tlabel <- if (tvar=="inMinorityTherm") "inMinority (own vs other-minority)" else "inMajority (own vs White)"
  ct_all <- cor.test(df$bias, df[[tvar]], use = "complete.obs")
  n_all  <- sum(complete.cases(df[,c("bias",tvar)]))
  cat(sprintf("\nr(γ, %s):\n", tlabel))
  cat(sprintf("  Overall: r = %.3f, p = %.4f, n = %d [CI %.3f, %.3f]\n",
              unname(ct_all$estimate), ct_all$p.value, n_all,
              ct_all$conf.int[1], ct_all$conf.int[2]))

  r_by_cond <- sapply(c("Minority","Majority"), function(cond) {
    sub <- filter(df, condition==cond)
    sub2 <- sub[complete.cases(sub[,c("bias",tvar)]),]
    ct_c <- cor.test(sub2$bias, sub2[[tvar]])
    cat(sprintf("  %s: r = %.3f, p = %.4f, n = %d\n",
                cond, unname(ct_c$estimate), ct_c$p.value, nrow(sub2)))
    c(r=unname(ct_c$estimate), n=nrow(sub2))
  })
  fz <- fisher_z_test(r_by_cond["r","Minority"], r_by_cond["r","Majority"],
                      r_by_cond["n","Minority"], r_by_cond["n","Majority"])
  cat(sprintf("  Fisher z condition diff: z = %.3f, p = %.4f\n", fz$z, fz$p))
}

cat("\n\n=== KEY: MCR correlations ===\n")
for (tvar in c("inMinorityTherm","inMajorityTherm")) {
  tlabel <- if (tvar=="inMinorityTherm") "inMinority" else "inMajority"
  ct_all <- cor.test(df$subject_mcr, df[[tvar]], use = "complete.obs")
  n_all  <- sum(complete.cases(df[,c("subject_mcr",tvar)]))
  cat(sprintf("r(MCR, %s): r = %.3f, p = %.4f, n = %d\n",
              tlabel, unname(ct_all$estimate), ct_all$p.value, n_all))
  for (cond in c("Minority","Majority")) {
    sub <- filter(df, condition==cond)
    sub2 <- sub[complete.cases(sub[,c("subject_mcr",tvar)]),]
    ct_c <- cor.test(sub2$subject_mcr, sub2[[tvar]])
    cat(sprintf("  %s: r = %.3f, p = %.4f, n = %d\n",
                cond, unname(ct_c$estimate), ct_c$p.value, nrow(sub2)))
  }
}

# ── 10. Save ───────────────────────────────────────────────────────────────────
write.csv(lm_df, here("Results","warmth_lm_s3.csv"), row.names=FALSE)

message("\nSaved:")
message("  Results/warmth_correlations_s3.csv")
message("  Results/warmth_lm_s3.csv")
