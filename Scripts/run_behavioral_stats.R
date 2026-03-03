# Statistical Analysis Script for Self-Anchoring
# This script fills in the behavioral STATS placeholders.

library(lme4)
library(lmerTest)
library(marginaleffects)
library(dplyr)
library(tidyverse)
library(broom.mixed)
library(emmeans)

results <- list()

# --- STUDY 2 ---
message("Analyzing Study 2...")
fullTest2 <- read.csv("Study 2/Cleaning/output/fullTest.csv")
indDiffs2 <- fullTest2 %>% filter(!duplicated(subID))

# Segregation (groupHomoph) by condition
indDiffs2$outgroup <- factor(indDiffs2$outgroup, levels = c("Not UCR", "UCLA", "CSU LA"))
m_gh2 <- lm(scale(groupHomoph) ~ outgroup, data = indDiffs2)
results$s2_gh_anova <- anova(m_gh2)
results$s2_gh_comparisons <- avg_comparisons(m_gh2)
results$s2_gh_means <- avg_predictions(m_gh2, variables = "outgroup")

# --- STUDY 3 ---
message("Analyzing Study 3...")
s3_path <- "Study 3/Cleaning/output/fullTest_fixed.csv"
if (file.exists(s3_path)) {
  fullTest3 <- read.csv(s3_path)
  indDiffs3 <- fullTest3 %>% filter(!duplicated(subID))
  
  # Segregation (groupHomoph) by condition
  # In Study 3, 'condition' column exists
  indDiffs3$condition <- factor(indDiffs3$condition, levels = c("Minority", "Majority"))
  if ("groupHomoph" %in% names(indDiffs3)) {
    m_gh3 <- lm(scale(groupHomoph) ~ condition, data = indDiffs3)
    results$s3_gh_anova <- anova(m_gh3)
    results$s3_gh_comparisons <- avg_comparisons(m_gh3)
    results$s3_gh_means <- avg_predictions(m_gh3, variables = "condition")
  }
  
  # Check for other behavioral interactions in Study 3
  # Similarity-to-self (WSR) * condition on ingChoiceN
  fullTest3$condition <- factor(fullTest3$condition, levels = c("Minority", "Majority"))
  m_wsr3 <- glmer(ingChoiceN ~ scale(WSR) * condition + (1 | subID) + (1 | trait), 
                  data = fullTest3, family = binomial)
  results$s3_wsr_interaction <- summary(m_wsr3)$coefficients
  results$s3_wsr_comparisons <- avg_comparisons(m_wsr3, variables = "WSR", by = "condition")
}

# Print results
print(results)
saveRDS(results, "behavioral_results.rds")
