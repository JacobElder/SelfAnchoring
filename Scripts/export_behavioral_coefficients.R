library(lme4)
library(broom.mixed)
library(dplyr)
library(readr)
library(here)

ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 300000))

rename_terms <- function(x) {
  recode(
    x,
    `(Intercept)` = "Intercept",
    `desirability.Z` = "Trait Desirability (z)",
    `selfResp.Z` = "Self-Evaluation (z)",
    `predicted.Z` = "Similarity-to-Self (z)",
    `novelHeld-Out` = "Novel Trait (vs. trained)",
    `outgroupUCLA` = "High-Status Outgroup (vs. negation)",
    `outgroupCSU LA` = "Low-Status Outgroup (vs. negation)",
    `conditionMajority` = "Majority Outgroup Condition (vs. minority)",
    `desirability.Z:outgroupUCLA` = "Trait Desirability x High-Status Outgroup",
    `desirability.Z:outgroupCSU LA` = "Trait Desirability x Low-Status Outgroup",
    `selfResp.Z:outgroupUCLA` = "Self-Evaluation x High-Status Outgroup",
    `selfResp.Z:outgroupCSU LA` = "Self-Evaluation x Low-Status Outgroup",
    `predicted.Z:outgroupUCLA` = "Similarity-to-Self x High-Status Outgroup",
    `predicted.Z:outgroupCSU LA` = "Similarity-to-Self x Low-Status Outgroup",
    `desirability.Z:conditionMajority` = "Trait Desirability x Majority Outgroup Condition",
    `selfResp.Z:conditionMajority` = "Self-Evaluation x Majority Outgroup Condition",
    `predicted.Z:conditionMajority` = "Similarity-to-Self x Majority Outgroup Condition",
    `predicted.Z:novelHeld-Out` = "Similarity-to-Self x Novel Trait",
    `novelHeld-Out:outgroupUCLA` = "Novel Trait x High-Status Outgroup",
    `novelHeld-Out:outgroupCSU LA` = "Novel Trait x Low-Status Outgroup",
    `novelHeld-Out:conditionMajority` = "Novel Trait x Majority Outgroup Condition",
    `predicted.Z:novelHeld-Out:outgroupUCLA` = "Similarity-to-Self x Novel Trait x High-Status Outgroup",
    `predicted.Z:novelHeld-Out:outgroupCSU LA` = "Similarity-to-Self x Novel Trait x Low-Status Outgroup",
    `predicted.Z:novelHeld-Out:conditionMajority` = "Similarity-to-Self x Novel Trait x Majority Outgroup Condition",
    .default = x
  )
}

tidy_model <- function(model, study, model_id, model_label) {
  broom.mixed::tidy(model, effects = "fixed") |>
    transmute(
      study = study,
      model = model_id,
      model_label = model_label,
      term = rename_terms(term),
      estimate = estimate,
      std_error = std.error,
      statistic = statistic,
      p_value = p.value
    )
}

fit_study_1 <- function() {
  dat <- read.csv(here("Study 1/Cleaning/output/fullTest.csv")) |>
    filter(!is.na(ingChoiceN)) |>
    mutate(
      ingChoiceN = as.integer(ingChoiceN),
      selfResp.Z = as.numeric(scale(selfResp)),
      predicted.Z = as.numeric(scale(predicted)),
      desirability.Z = as.numeric(scale(desirability)),
      novel = factor(novel, levels = c(0, 1), labels = c("Trained", "Held-Out"))
    )

  m1 <- glmer(
    ingChoiceN ~ desirability.Z + (desirability.Z | subID) + (1 | trait),
    data = dat, family = binomial, control = ctrl, nAGQ = 1
  )
  m2 <- glmer(
    ingChoiceN ~ selfResp.Z + desirability.Z + (selfResp.Z | subID) + (1 | trait),
    data = dat, family = binomial, control = ctrl, nAGQ = 1
  )
  m3 <- glmer(
    ingChoiceN ~ predicted.Z + desirability.Z + (predicted.Z | subID) + (1 | trait),
    data = dat, family = binomial, control = ctrl, nAGQ = 1
  )
  m4 <- glmer(
    ingChoiceN ~ predicted.Z * novel + desirability.Z + (predicted.Z + novel | subID) + (1 | trait),
    data = dat, family = binomial, control = ctrl, nAGQ = 1
  )

  bind_rows(
    tidy_model(m1, "Study 1", "M1", "Trait Desirability"),
    tidy_model(m2, "Study 1", "M2", "Self-Evaluation"),
    tidy_model(m3, "Study 1", "M3", "Similarity-to-Self"),
    tidy_model(m4, "Study 1", "M4", "Similarity-to-Self x Novel Trait")
  )
}

fit_study_2 <- function() {
  dat <- read.csv(here("Study 2/Cleaning/output/fullTest.csv")) |>
    filter(!is.na(ingChoiceN)) |>
    mutate(
      ingChoiceN = as.integer(ingChoiceN),
      selfResp.Z = as.numeric(scale(selfResp)),
      predicted.Z = as.numeric(scale(predicted)),
      desirability.Z = as.numeric(scale(desirability)),
      novel = factor(novel, levels = c(0, 1), labels = c("Trained", "Held-Out")),
      outgroup = factor(outgroup),
      outgroup = relevel(outgroup, ref = "Not UCR")
    )

  m1 <- glmer(
    ingChoiceN ~ desirability.Z * outgroup + (desirability.Z | subID) + (1 | trait),
    data = dat, family = binomial, control = ctrl, nAGQ = 1
  )
  m2 <- glmer(
    ingChoiceN ~ selfResp.Z * outgroup + desirability.Z + (selfResp.Z | subID) + (1 | trait),
    data = dat, family = binomial, control = ctrl, nAGQ = 1
  )
  m3 <- glmer(
    ingChoiceN ~ predicted.Z * outgroup + desirability.Z + (predicted.Z | subID) + (1 | trait),
    data = dat, family = binomial, control = ctrl, nAGQ = 1
  )
  m4 <- glmer(
    ingChoiceN ~ predicted.Z * novel + outgroup + desirability.Z + (predicted.Z + novel | subID) + (1 | trait),
    data = dat, family = binomial, control = ctrl, nAGQ = 1
  )
  m5 <- glmer(
    ingChoiceN ~ predicted.Z * novel * outgroup + desirability.Z + (predicted.Z + novel | subID) + (1 | trait),
    data = dat, family = binomial, control = ctrl, nAGQ = 1
  )

  bind_rows(
    tidy_model(m1, "Study 2", "M1", "Trait Desirability x Outgroup"),
    tidy_model(m2, "Study 2", "M2", "Self-Evaluation x Outgroup"),
    tidy_model(m3, "Study 2", "M3", "Similarity-to-Self x Outgroup"),
    tidy_model(m4, "Study 2", "M4", "Similarity-to-Self x Novel Trait"),
    tidy_model(m5, "Study 2", "M5", "Similarity-to-Self x Novel Trait x Outgroup")
  )
}

fit_study_3 <- function() {
  dat <- read.csv(here("Study 3/Cleaning/output/fullTest_fixed.csv")) |>
    filter(!is.na(ingChoiceN)) |>
    mutate(
      ingChoiceN = as.integer(ingChoiceN),
      selfResp.Z = as.numeric(scale(selfResp)),
      predicted.Z = as.numeric(scale(predicted)),
      desirability.Z = as.numeric(scale(desirability)),
      novel = factor(novel, levels = c(0, 1), labels = c("Trained", "Held-Out")),
      condition = factor(condition),
      condition = relevel(condition, ref = "Minority")
    )

  m1 <- glmer(
    ingChoiceN ~ desirability.Z * condition + (desirability.Z | subID) + (1 | trait),
    data = dat, family = binomial, control = ctrl, nAGQ = 1
  )
  m2 <- glmer(
    ingChoiceN ~ selfResp.Z * condition + desirability.Z + (selfResp.Z | subID) + (1 | trait),
    data = dat, family = binomial, control = ctrl, nAGQ = 1
  )
  m3 <- glmer(
    ingChoiceN ~ predicted.Z * condition + desirability.Z + (predicted.Z | subID) + (1 | trait),
    data = dat, family = binomial, control = ctrl, nAGQ = 1
  )
  m4 <- glmer(
    ingChoiceN ~ predicted.Z * novel + condition + desirability.Z + (predicted.Z + novel | subID) + (1 | trait),
    data = dat, family = binomial, control = ctrl, nAGQ = 1
  )
  m5 <- glmer(
    ingChoiceN ~ predicted.Z * novel * condition + desirability.Z + (predicted.Z + novel | subID) + (1 | trait),
    data = dat, family = binomial, control = ctrl, nAGQ = 1
  )

  bind_rows(
    tidy_model(m1, "Study 3", "M1", "Trait Desirability x Condition"),
    tidy_model(m2, "Study 3", "M2", "Self-Evaluation x Condition"),
    tidy_model(m3, "Study 3", "M3", "Similarity-to-Self x Condition"),
    tidy_model(m4, "Study 3", "M4", "Similarity-to-Self x Novel Trait"),
    tidy_model(m5, "Study 3", "M5", "Similarity-to-Self x Novel Trait x Condition")
  )
}

dir.create(here("Results"), showWarnings = FALSE)

all_coefs <- bind_rows(
  fit_study_1(),
  fit_study_2(),
  fit_study_3()
)

write_csv(all_coefs, here("Results", "behavioral_model_coefficients.csv"))

split(all_coefs, all_coefs$study) |>
  lapply(function(df) {
    write_csv(df, here("Results", paste0(gsub(" ", "_", tolower(unique(df$study))), "_behavioral_model_coefficients.csv")))
  })

message("Saved behavioral coefficient tables to Results/.")
