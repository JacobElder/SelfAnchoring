# Pooled Behavioral GLMMs — Similarity-to-Self and Generalization Across Studies 1–3
#
# Pools all three studies into a single GLMM to formally test whether
# similarity-based self-anchoring and generalization effects differ across
# intergroup contexts (conditions × studies), using Study 1's minimal groups
# as the reference category.
#
# Random effects:
#   (1 | study)            — study-level intercept (3 levels; random to acknowledge
#                             between-study variation beyond fixed condition contrasts)
#   (focal | study:subID)  — subject intercepts and slopes nested within study
#   (1 | study:trait)      — trait intercepts nested within study (same 148 traits
#                             appear in each study but their ingroup-attribution
#                             baselines differ across intergroup contexts)
#
# Conditions (reference = "Minimal_Groups"):
#   S1: Minimal_Groups
#   S2: Univ_Negation, Univ_HighStatus, Univ_LowStatus
#   S3: Racial_Minority, Racial_Majority
#
# OUTPUT: Results/marginal_effects_pooled_ames.csv
#         Results/marginal_effects_pooled_predictions.csv

library(lme4)
library(marginaleffects)
library(tidyverse)
library(here)

options(marginaleffects_safe = FALSE)

# ── 1. Load and harmonize data ────────────────────────────────────────────────

s1 <- read.csv(here("Study 1/Cleaning/output/fullTest.csv")) |>
  filter(!is.na(ingChoiceN)) |>
  mutate(
    study     = "S1",
    condition = "Minimal_Groups"
  )

s2 <- read.csv(here("Study 2/Cleaning/output/fullTest.csv")) |>
  filter(!is.na(ingChoiceN)) |>
  mutate(
    study     = "S2",
    condition = recode(outgroup,
      "Not UCR" = "Univ_Negation",
      "UCLA"    = "Univ_HighStatus",
      "CSU LA"  = "Univ_LowStatus"
    )
  )

s3 <- read.csv(here("Study 3/Cleaning/output/fullTest_fixed.csv")) |>
  filter(!is.na(ingChoiceN)) |>
  mutate(
    study     = "S3",
    condition = recode(condition,
      "Minority" = "Racial_Minority",
      "Majority" = "Racial_Majority"
    )
  )

# Common columns to keep
keep_cols <- c("study", "subID", "trait", "ingChoiceN",
               "selfResp", "predicted", "desirability", "novel", "condition")

pooled <- bind_rows(
  s1 |> select(all_of(keep_cols)),
  s2 |> select(all_of(keep_cols)),
  s3 |> select(all_of(keep_cols))
)

# ── 2. Prepare variables ──────────────────────────────────────────────────────

# study:subID interaction used directly in formulas for explicit nesting notation

# Z-score predictors within study (so 1 unit = within-study SD for comparability)
pooled <- pooled |>
  group_by(study) |>
  mutate(
    selfResp.Z    = as.numeric(scale(selfResp)),
    predicted.Z   = as.numeric(scale(predicted)),
    desirability.Z = as.numeric(scale(desirability))
  ) |>
  ungroup()

pooled$ingChoiceN <- as.integer(pooled$ingChoiceN)
pooled$novel      <- factor(pooled$novel, levels = c(0, 1),
                            labels = c("Trained", "Held-Out"))

# Condition: reference = Minimal_Groups (S1)
cond_levels <- c("Minimal_Groups",
                 "Univ_Negation", "Univ_HighStatus", "Univ_LowStatus",
                 "Racial_Minority", "Racial_Majority")
pooled$condition <- factor(pooled$condition, levels = cond_levels)

message("Pooled N = ", length(unique(paste(pooled$study, pooled$subID))), " subjects, ",
        nrow(pooled), " observations")
message("Conditions: ", paste(levels(pooled$condition), collapse = ", "))

ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 300000))

# ── 3. Models ─────────────────────────────────────────────────────────────────

# M1: Desirability × Condition
message("Fitting M1: Desirability × Condition...")
m1 <- glmer(
  ingChoiceN ~ desirability.Z * condition + (1 | study) +
    (desirability.Z | study:subID) + (1 | study:trait),
  data = pooled, family = binomial, control = ctrl, nAGQ = 1
)

# M2: Self-Evaluations × Condition
message("Fitting M2: Self-Evaluations × Condition...")
m2 <- glmer(
  ingChoiceN ~ selfResp.Z * condition + desirability.Z + (1 | study) +
    (selfResp.Z | study:subID) + (1 | study:trait),
  data = pooled, family = binomial, control = ctrl, nAGQ = 1
)

# M3: Similarity-to-Self × Condition
message("Fitting M3: Similarity-to-Self × Condition...")
m3 <- glmer(
  ingChoiceN ~ predicted.Z * condition + desirability.Z + (1 | study) +
    (predicted.Z | study:subID) + (1 | study:trait),
  data = pooled, family = binomial, control = ctrl, nAGQ = 1
)

# M4: Similarity-to-Self × Novel (averaged over condition)
message("Fitting M4: Similarity-to-Self × Novel...")
m4 <- suppressMessages(glmer(
  ingChoiceN ~ predicted.Z * novel + condition + desirability.Z + (1 | study) +
    (predicted.Z + novel | study:subID) + (1 | study:trait),
  data = pooled, family = binomial, control = ctrl, nAGQ = 1
))

# M5: Similarity-to-Self × Novel × Condition (3-way — generalization by context)
message("Fitting M5: Similarity-to-Self × Novel × Condition (3-way)...")
m5 <- suppressMessages(glmer(
  ingChoiceN ~ predicted.Z * novel * condition + desirability.Z + (1 | study) +
    (predicted.Z + novel | study:subID) + (1 | study:trait),
  data = pooled, family = binomial, control = ctrl, nAGQ = 1
))

# ── 4. Average Marginal Effects ───────────────────────────────────────────────
message("\nComputing average marginal effects...")

# Overall AMEs (marginalizing over condition)
ame_desirability <- avg_comparisons(m1, variables = "desirability.Z")
ame_selfResp     <- avg_comparisons(m2, variables = "selfResp.Z")
ame_predicted    <- avg_comparisons(m3, variables = "predicted.Z")

# AMEs by condition
ame_des_cond   <- avg_comparisons(m1, variables = "desirability.Z", by = "condition")
ame_self_cond  <- avg_comparisons(m2, variables = "selfResp.Z",     by = "condition")
ame_pred_cond  <- avg_comparisons(m3, variables = "predicted.Z",    by = "condition")

# Generalization (Trained vs. Held-Out), averaged over condition
ame_novel      <- avg_comparisons(m4, variables = "predicted.Z", by = "novel")

# 3-way: Similarity-to-Self × Novel × Condition
ame_novel_cond <- avg_comparisons(m5, variables = "predicted.Z",
                                  by = c("novel", "condition"))

# Predicted P(ingroup) at -1/0/+1 SD Similarity-to-Self, by condition
avg_pred_cond <- avg_predictions(
  m3,
  variables = list(predicted.Z = c(-1, 0, 1),
                   condition   = levels(pooled$condition))
)

# ── 5. Print results ──────────────────────────────────────────────────────────

cat("\n=== Overall AMEs (pooled across all conditions) ===\n")
for (nm in c("desirability.Z", "selfResp.Z", "predicted.Z")) {
  obj <- switch(nm,
    desirability.Z = ame_desirability,
    selfResp.Z     = ame_selfResp,
    predicted.Z    = ame_predicted)
  d <- as.data.frame(obj)
  cat(sprintf("%-15s  AME=%+.4f  SE=%.4f  CI=[%.4f, %.4f]  p=%.4f\n",
              nm, d$estimate, d$std.error, d$conf.low, d$conf.high, d$p.value))
}

cat("\n=== AMEs by Condition ===\n")
cat("--- Similarity-to-Self ---\n")
print(as.data.frame(ame_pred_cond)[, c("condition", "estimate", "std.error",
                                        "conf.low", "conf.high", "p.value")],
      digits = 3, row.names = FALSE)

cat("--- Self-Evaluations ---\n")
print(as.data.frame(ame_self_cond)[, c("condition", "estimate", "std.error",
                                        "conf.low", "conf.high", "p.value")],
      digits = 3, row.names = FALSE)

cat("--- Desirability ---\n")
print(as.data.frame(ame_des_cond)[, c("condition", "estimate", "std.error",
                                       "conf.low", "conf.high", "p.value")],
      digits = 3, row.names = FALSE)

cat("\n=== Generalization (Trained vs. Held-Out, pooled over condition) ===\n")
print(as.data.frame(ame_novel)[, c("novel", "estimate", "std.error",
                                    "conf.low", "conf.high", "p.value")],
      digits = 3, row.names = FALSE)

cat("\n=== 3-way: Similarity-to-Self × Novel × Condition ===\n")
print(as.data.frame(ame_novel_cond)[, c("novel", "condition", "estimate",
                                          "std.error", "conf.low", "conf.high", "p.value")],
      digits = 3, row.names = FALSE)

cat("\n=== Predicted P(ingroup) at -1/0/+1 SD Similarity-to-Self, by Condition ===\n")
print(as.data.frame(avg_pred_cond)[, c("predicted.Z", "condition",
                                        "estimate", "conf.low", "conf.high")],
      digits = 3, row.names = FALSE)

# ── 6. Save ───────────────────────────────────────────────────────────────────
ames <- bind_rows(
  as.data.frame(ame_desirability) |> mutate(model = "M1_overall"),
  as.data.frame(ame_des_cond)     |> mutate(model = "M1_by_cond"),
  as.data.frame(ame_selfResp)     |> mutate(model = "M2_overall"),
  as.data.frame(ame_self_cond)    |> mutate(model = "M2_by_cond"),
  as.data.frame(ame_predicted)    |> mutate(model = "M3_overall"),
  as.data.frame(ame_pred_cond)    |> mutate(model = "M3_by_cond"),
  as.data.frame(ame_novel)        |> mutate(model = "M4_novel"),
  as.data.frame(ame_novel_cond)   |> mutate(model = "M5_novel_by_cond")
)

write.csv(ames,
          here("Results", "marginal_effects_pooled_ames.csv"),
          row.names = FALSE)
write.csv(as.data.frame(avg_pred_cond),
          here("Results", "marginal_effects_pooled_predictions.csv"),
          row.names = FALSE)

message("\nSaved to Results/marginal_effects_pooled_ames.csv")
message("       and Results/marginal_effects_pooled_predictions.csv")
