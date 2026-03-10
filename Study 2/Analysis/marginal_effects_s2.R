# Marginal Effects — Study 2 Behavioral GLMMs
# Parallel structure to marginal_effects_s1.R.
# Condition: outgroup ("Not UCR" = negation ref, "UCLA" = high-status, "CSU LA" = low-status)
# Covariates: propCorrLOO.Z (combined ingroup+outgroup LOO accuracy, parallel to S1)

library(lme4)
library(marginaleffects)
library(tidyverse)
library(here)
options(marginaleffects_safe = FALSE)

fullTest <- read.csv(here("Study 2/Cleaning/output/fullTest.csv")) |>
  filter(!is.na(ingChoiceN))

fullTest$ingChoiceN      <- as.integer(fullTest$ingChoiceN)
fullTest$selfResp.Z      <- scale(fullTest$selfResp)[,1]
fullTest$predicted.Z     <- scale(fullTest$predicted)[,1]
fullTest$desirability.Z  <- scale(fullTest$desirability)[,1]
fullTest$novel    <- factor(fullTest$novel,   levels = c(0,1), labels = c("Trained","Held-Out"))
fullTest$outgroup <- factor(fullTest$outgroup)
fullTest$outgroup <- relevel(fullTest$outgroup, ref = "Not UCR")

ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))

# M1: Desirability
message("Fitting M1: Desirability...")
m1 <- glmer(
  ingChoiceN ~ desirability.Z * outgroup +
    (desirability.Z | subID) + (1 | trait),
  data = fullTest, family = binomial, control = ctrl, nAGQ = 1
)

# M2: Self-evaluations
message("Fitting M2: Self-evaluations...")
m2 <- glmer(
  ingChoiceN ~ selfResp.Z * outgroup + desirability.Z +
    (selfResp.Z | subID) + (1 | trait),
  data = fullTest, family = binomial, control = ctrl, nAGQ = 1
)

# M3: Similarity-to-Self
message("Fitting M3: Similarity-to-Self...")
m3 <- glmer(
  ingChoiceN ~ predicted.Z * outgroup + desirability.Z +
    (predicted.Z | subID) + (1 | trait),
  data = fullTest, family = binomial, control = ctrl, nAGQ = 1
)

# M4: Similarity-to-Self × Novel
message("Fitting M4: Similarity-to-Self x Novel...")
m4 <- suppressMessages(glmer(
  ingChoiceN ~ predicted.Z * novel + outgroup + desirability.Z +
    (predicted.Z + novel | subID) + (1 | trait),
  data = fullTest, family = binomial, control = ctrl, nAGQ = 1
))

# M5: Similarity-to-Self × Novel × Outgroup — tests whether generalization
# differs by outgroup contrast condition.
message("Fitting M5: Similarity-to-Self x Novel x Outgroup...")
m5 <- suppressMessages(glmer(
  ingChoiceN ~ predicted.Z * novel * outgroup + desirability.Z +
    (predicted.Z + novel | subID) + (1 | trait),
  data = fullTest, family = binomial, control = ctrl, nAGQ = 1
))

message("\nComputing average marginal effects...")

# Overall AMEs (marginalizing over outgroup)
ame_desirability <- avg_comparisons(m1, variables = "desirability.Z")
ame_selfResp     <- avg_comparisons(m2, variables = "selfResp.Z")
ame_predicted    <- avg_comparisons(m3, variables = "predicted.Z")

# AMEs by condition
ame_des_cond  <- avg_comparisons(m1, variables = "desirability.Z", by = "outgroup")
ame_self_cond <- avg_comparisons(m2, variables = "selfResp.Z",     by = "outgroup")
ame_pred_cond <- avg_comparisons(m3, variables = "predicted.Z",    by = "outgroup")

# Generalization test (averaged over condition)
ame_novel     <- avg_comparisons(m4, variables = "predicted.Z", by = "novel")

# Generalization test by outgroup condition (3-way)
ame_novel_cond <- avg_comparisons(m5, variables = "predicted.Z", by = c("novel", "outgroup"))

# ── Dense-grid predictions for smooth plotting (Figure 7) ────────────────────
x_des  <- seq(min(fullTest$desirability.Z, na.rm=TRUE),
               max(fullTest$desirability.Z, na.rm=TRUE), length.out = 50)
x_self <- seq(min(fullTest$selfResp.Z,     na.rm=TRUE),
               max(fullTest$selfResp.Z,     na.rm=TRUE), length.out = 50)
x_ss   <- seq(min(fullTest$predicted.Z,    na.rm=TRUE),
               max(fullTest$predicted.Z,    na.rm=TRUE), length.out = 50)

outgroups <- c("Not UCR","UCLA","CSU LA")

message("Computing dense-grid predictions by condition...")
avg_pred_des_cond  <- avg_predictions(m1, variables = list(desirability.Z = x_des,  outgroup = outgroups))
avg_pred_self_cond <- avg_predictions(m2, variables = list(selfResp.Z     = x_self, outgroup = outgroups))
avg_pred_cond      <- avg_predictions(m3, variables = list(predicted.Z    = x_ss,   outgroup = outgroups))

cat("\n=== Overall AMEs ===\n")
for (nm in c("desirability.Z","selfResp.Z","predicted.Z")) {
  obj <- switch(nm,
    desirability.Z = ame_desirability,
    selfResp.Z     = ame_selfResp,
    predicted.Z    = ame_predicted)
  d <- as.data.frame(obj)
  cat(sprintf("%-15s  AME=%+.4f  SE=%.4f  CI=[%.4f, %.4f]  p=%.4f\n",
              nm, d$estimate, d$std.error, d$conf.low, d$conf.high, d$p.value))
}

cat("\n=== AMEs by Condition ===\n")
cat("--- Desirability ---\n")
print(as.data.frame(ame_des_cond)[, c("outgroup","estimate","std.error","conf.low","conf.high","p.value")], digits=3, row.names=FALSE)
cat("--- Self-Evaluations ---\n")
print(as.data.frame(ame_self_cond)[, c("outgroup","estimate","std.error","conf.low","conf.high","p.value")], digits=3, row.names=FALSE)
cat("--- Similarity-to-Self ---\n")
print(as.data.frame(ame_pred_cond)[, c("outgroup","estimate","std.error","conf.low","conf.high","p.value")], digits=3, row.names=FALSE)

cat("\n=== Generalization (Trained vs Held-Out, averaged over condition) ===\n")
print(as.data.frame(ame_novel)[, c("novel","estimate","std.error","conf.low","conf.high","p.value")], digits=3, row.names=FALSE)

cat("\n=== Generalization by Condition (3-way: SS x Novel x Outgroup) ===\n")
print(as.data.frame(ame_novel_cond)[, c("novel","outgroup","estimate","std.error","conf.low","conf.high","p.value")], digits=3, row.names=FALSE)

# Save
ames <- bind_rows(
  as.data.frame(ame_desirability) |> mutate(model="M1_overall"),
  as.data.frame(ame_des_cond)     |> mutate(model="M1_by_cond"),
  as.data.frame(ame_selfResp)     |> mutate(model="M2_overall"),
  as.data.frame(ame_self_cond)    |> mutate(model="M2_by_cond"),
  as.data.frame(ame_predicted)    |> mutate(model="M3_overall"),
  as.data.frame(ame_pred_cond)    |> mutate(model="M3_by_cond"),
  as.data.frame(ame_novel)        |> mutate(model="M4_novel"),
  as.data.frame(ame_novel_cond)   |> mutate(model="M5_novel_by_cond")
)
write.csv(ames, here("Results","marginal_effects_s2_ames.csv"), row.names=FALSE)

preds <- bind_rows(
  as.data.frame(avg_pred_des_cond)  |>
    rename(predictor_val = desirability.Z) |>
    mutate(model = "M1", predictor = "desirability"),
  as.data.frame(avg_pred_self_cond) |>
    rename(predictor_val = selfResp.Z) |>
    mutate(model = "M2", predictor = "selfResp"),
  as.data.frame(avg_pred_cond)      |>
    rename(predictor_val = predicted.Z) |>
    mutate(model = "M3", predictor = "ss")
)
write.csv(preds, here("Results","marginal_effects_s2_predictions.csv"), row.names=FALSE)
message("Saved to Results/marginal_effects_s2_ames.csv and marginal_effects_s2_predictions.csv")
