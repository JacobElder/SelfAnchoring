# Marginal Effects — Study 1 Behavioral GLMMs
# Uses marginaleffects::avg_comparisons() to express GLMM results in probability units.
#
# WHY: Log-odds (β) and odds ratios (OR) are hard to communicate. avg_comparisons()
# gives Average Marginal Effects (AMEs): the average change in P(ingroup choice) for
# a 1-unit (= 1-SD, since predictors are Z-scored) increase in each predictor.
#
# MODELS (matching manuscript narrative):
#   M1: Desirability → ingroup choice
#   M2: Self-evaluations (selfResp.Z) → ingroup choice
#   M3: Similarity-to-Self (predicted.Z) → ingroup choice  ← primary finding
#   M4: Similarity-to-Self × Novel (trained vs. held-out generalization test)
#
# OUTPUT: CSV + console table for user verification before manuscript insertion.

library(lme4)
library(marginaleffects)
library(tidyverse)
library(here)

# ── Load data ─────────────────────────────────────────────────────────────────
fullTest  <- read.csv(here("Study 1/Cleaning/output/fullTest.csv"))
fullTrain <- read.csv(here("Study 1/Cleaning/output/fullTrain.csv"))

fullTest <- fullTest |> filter(!is.na(ingChoiceN))

# Predictors (Z-scored so 1 unit = 1 SD)
fullTest$ingChoiceN    <- as.integer(fullTest$ingChoiceN)  # 0/1
fullTest$selfResp.Z    <- scale(fullTest$selfResp)[,1]
fullTest$predicted.Z   <- scale(fullTest$predicted)[,1]
fullTest$desirability.Z <- scale(fullTest$desirability)[,1]
fullTest$novel         <- factor(fullTest$novel, levels = c(0, 1),
                                  labels = c("Trained", "Held-Out"))

ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))

# ── M1: Desirability ──────────────────────────────────────────────────────────
message("Fitting M1: Desirability...")
m1 <- glmer(
  ingChoiceN ~ desirability.Z +
    (desirability.Z | subID) + (1 | trait),
  data = fullTest, family = binomial, control = ctrl, nAGQ = 1
)

# ── M2: Self-evaluations ──────────────────────────────────────────────────────
message("Fitting M2: Self-evaluations...")
m2 <- glmer(
  ingChoiceN ~ selfResp.Z + desirability.Z +
    (selfResp.Z | subID) + (1 | trait),
  data = fullTest, family = binomial, control = ctrl, nAGQ = 1
)

# ── M3: Similarity-to-Self (cross-validated predicted) ────────────────────────
message("Fitting M3: Similarity-to-Self...")
m3 <- glmer(
  ingChoiceN ~ predicted.Z + desirability.Z +
    (predicted.Z | subID) + (1 | trait),
  data = fullTest, family = binomial, control = ctrl, nAGQ = 1
)

# ── M4: Similarity-to-Self × Novel ────────────────────────────────────────────
message("Fitting M4: Similarity-to-Self x Novel...")
m4 <- glmer(
  ingChoiceN ~ predicted.Z * novel + desirability.Z +
    (predicted.Z + novel | subID) + (1 | trait),
  data = fullTest, family = binomial, control = ctrl, nAGQ = 1
)

# ── Average Marginal Effects ──────────────────────────────────────────────────
message("\nComputing average marginal effects...")

# For M1-M3: AME of focal predictor (1-SD increase → change in P(ingroup))
ame_desirability <- avg_comparisons(m1, variables = "desirability.Z")
ame_selfResp     <- avg_comparisons(m2, variables = "selfResp.Z")
ame_predicted    <- avg_comparisons(m3, variables = "predicted.Z")

# For M4: AMEs by novel status (separate slopes) + interaction (difference in AMEs)
ame_predicted_by_novel <- avg_comparisons(
  m4,
  variables  = "predicted.Z",
  by         = "novel"         # AME of SS separately for Trained vs. Held-Out
)

ame_novel_interaction <- avg_comparisons(
  m4,
  variables  = "predicted.Z",
  by         = "novel",
  hypothesis = "b1 = b2"       # test: AME_trained == AME_held_out
)

# ── Dense-grid predictions for smooth plotting (Figure 7) ────────────────────
# Use observed data range (not ±1 SD) so curves span the full predictor spread,
# matching the original dissertation ggpredict() style.
x_des  <- seq(min(fullTest$desirability.Z, na.rm=TRUE),
               max(fullTest$desirability.Z, na.rm=TRUE), length.out = 50)
x_self <- seq(min(fullTest$selfResp.Z,     na.rm=TRUE),
               max(fullTest$selfResp.Z,     na.rm=TRUE), length.out = 50)
x_ss   <- seq(min(fullTest$predicted.Z,    na.rm=TRUE),
               max(fullTest$predicted.Z,    na.rm=TRUE), length.out = 50)

message("Computing dense-grid predictions for M1, M2, M3, M4...")
avg_pred_desirability <- avg_predictions(m1, variables = list(desirability.Z = x_des))
avg_pred_selfResp     <- avg_predictions(m2, variables = list(selfResp.Z     = x_self))
avg_pred_predicted    <- avg_predictions(m3, variables = list(predicted.Z    = x_ss))
avg_pred_novel        <- avg_predictions(m4, variables = list(predicted.Z    = x_ss,
                                                               novel = c("Trained", "Held-Out")))

# ── Print AME summary ─────────────────────────────────────────────────────────
fmt_ame <- function(x, label) {
  cat(sprintf("\n=== %s ===\n", label))
  print(as.data.frame(x)[, c("term", "contrast", "estimate", "std.error", "conf.low", "conf.high", "p.value")],
        digits = 3, row.names = FALSE)
}

fmt_ame(ame_desirability,       "M1: AME of Desirability (1-SD → ΔP)")
fmt_ame(ame_selfResp,           "M2: AME of Self-Evaluations (1-SD → ΔP)")
fmt_ame(ame_predicted,          "M3: AME of Similarity-to-Self (1-SD → ΔP)")
fmt_ame(ame_predicted_by_novel, "M4: AME of Similarity-to-Self by Novel Status")
cat("\n=== M4: Interaction test (Trained AME == Held-Out AME?) ===\n")
print(as.data.frame(ame_novel_interaction)[,
  intersect(c("term","estimate","std.error","conf.low","conf.high","p.value"),
            names(as.data.frame(ame_novel_interaction)))],
  digits = 3, row.names = FALSE)

# ── Save results ──────────────────────────────────────────────────────────────
results <- bind_rows(
  as.data.frame(ame_desirability)       |> mutate(model = "M1_desirability"),
  as.data.frame(ame_selfResp)           |> mutate(model = "M2_selfResp"),
  as.data.frame(ame_predicted)          |> mutate(model = "M3_SimilarityToSelf"),
  as.data.frame(ame_predicted_by_novel) |> mutate(model = "M4_SimilarityToSelf_by_novel")
)
write.csv(results, here("Results", "marginal_effects_s1_ames.csv"), row.names = FALSE)

preds <- bind_rows(
  as.data.frame(avg_pred_desirability) |>
    rename(predictor_val = desirability.Z) |>
    mutate(model = "M1", predictor = "desirability"),
  as.data.frame(avg_pred_selfResp) |>
    rename(predictor_val = selfResp.Z) |>
    mutate(model = "M2", predictor = "selfResp"),
  as.data.frame(avg_pred_predicted) |>
    rename(predictor_val = predicted.Z) |>
    mutate(model = "M3", predictor = "ss"),
  as.data.frame(avg_pred_novel) |>
    rename(predictor_val = predicted.Z) |>
    mutate(model = "M4", predictor = "ss")
)
write.csv(preds, here("Results", "marginal_effects_s1_predictions.csv"), row.names = FALSE)

message("\nResults saved to Results/marginal_effects_s1_ames.csv")
message("           and Results/marginal_effects_s1_predictions.csv")
