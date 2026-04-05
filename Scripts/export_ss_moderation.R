library(lme4)
library(broom.mixed)
library(dplyr)
library(readr)
library(here)
library(parallel)

ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))

moderators <- c(
  "DS" = "Dialectical Self-Views",
  "SCC" = "Self-Concept Clarity",
  "NFC" = "Need for Cognition",
  "NTB" = "Need to Belong",
  "RSE" = "Self-Esteem",
  "SI" = "Social Identification",
  "Proto" = "Perceived Prototypicality",
  "SING.Ind" = "Independent Self-Construal",
  "SING.Inter" = "Interdependent Self-Construal"
)

tidy_select <- function(model, study, moderator_label, terms_keep) {
  broom.mixed::tidy(model, effects = "fixed") |>
    filter(term %in% terms_keep) |>
    mutate(study = study, moderator = moderator_label) |>
    select(study, moderator, term, estimate, std.error, statistic, p.value)
}

fit_s1 <- function() {
  dat <- read.csv(here("Study 1/Cleaning/output/fullTest.csv")) |>
    filter(!is.na(ingChoiceN)) |>
    mutate(
      ingChoiceN = as.integer(ingChoiceN),
      predicted.Z = as.numeric(scale(predicted)),
      desirability.Z = as.numeric(scale(desirability))
    )

  bind_rows(mclapply(names(moderators), function(v) {
    message("Study 1: ", moderators[[v]])
    dat$mod.Z <- as.numeric(scale(dat[[v]]))
    m <- glmer(
      ingChoiceN ~ predicted.Z * mod.Z + desirability.Z + (predicted.Z | subID) + (1 | trait),
      data = dat, family = binomial, control = ctrl, nAGQ = 1
    )
    tidy_select(m, "Study 1", moderators[[v]], "predicted.Z:mod.Z")
  }, mc.cores = 4))
}

fit_s2 <- function() {
  dat <- read.csv(here("Study 2/Cleaning/output/fullTest.csv")) |>
    filter(!is.na(ingChoiceN)) |>
    mutate(
      ingChoiceN = as.integer(ingChoiceN),
      predicted.Z = as.numeric(scale(predicted)),
      desirability.Z = as.numeric(scale(desirability)),
      outgroup = factor(outgroup),
      outgroup = relevel(outgroup, ref = "Not UCR")
    )

  bind_rows(mclapply(names(moderators), function(v) {
    message("Study 2: ", moderators[[v]])
    dat$mod.Z <- as.numeric(scale(dat[[v]]))
    m <- glmer(
      ingChoiceN ~ predicted.Z * outgroup * mod.Z + desirability.Z +
        (predicted.Z | subID) + (1 | trait),
      data = dat, family = binomial, control = ctrl, nAGQ = 1
    )
    tidy_select(
      m,
      "Study 2",
      moderators[[v]],
      c("predicted.Z:mod.Z", "predicted.Z:outgroupUCLA:mod.Z", "predicted.Z:outgroupCSU LA:mod.Z")
    )
  }, mc.cores = 4))
}

fit_s3 <- function() {
  dat <- read.csv(here("Study 3/Cleaning/output/fullTest_fixed.csv")) |>
    filter(!is.na(ingChoiceN)) |>
    mutate(
      ingChoiceN = as.integer(ingChoiceN),
      predicted.Z = as.numeric(scale(predicted)),
      desirability.Z = as.numeric(scale(desirability)),
      condition = factor(condition),
      condition = relevel(condition, ref = "Minority")
    )

  bind_rows(mclapply(names(moderators), function(v) {
    message("Study 3: ", moderators[[v]])
    dat$mod.Z <- as.numeric(scale(dat[[v]]))
    m <- glmer(
      ingChoiceN ~ predicted.Z * condition * mod.Z + desirability.Z +
        (predicted.Z | subID) + (1 | trait),
      data = dat, family = binomial, control = ctrl, nAGQ = 1
    )
    tidy_select(
      m,
      "Study 3",
      moderators[[v]],
      c("predicted.Z:mod.Z", "predicted.Z:conditionMajority:mod.Z")
    )
  }, mc.cores = 4))
}

out <- bind_rows(fit_s1(), fit_s2(), fit_s3()) |>
  mutate(
    term = recode(
      term,
      `predicted.Z:mod.Z` = "Similarity-to-Self x Moderator",
      `predicted.Z:outgroupUCLA:mod.Z` = "Similarity-to-Self x High-Status Outgroup x Moderator",
      `predicted.Z:outgroupCSU LA:mod.Z` = "Similarity-to-Self x Low-Status Outgroup x Moderator",
      `predicted.Z:conditionMajority:mod.Z` = "Similarity-to-Self x Majority Outgroup Condition x Moderator"
    )
  )

write_csv(out, here("Results", "ss_moderation_results.csv"))
print(out)
