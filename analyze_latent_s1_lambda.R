# Latent Parameter Analysis for Study 1 (Minimal Groups) - Asymmetric + Lambda Model
library(cmdstanr)
library(tidyverse)

# 1. Load Model Fit
file_path <- "fit_s1_asym_lambda.rds"
fit_s1 <- readRDS(file_path)

draws <- fit_s1$draws(format = "df")

# 2. Group-level Parameters
# mu_pr[1] = tau, [2] = m_in, [3] = m_out, [4] = bias, [5] = lambda
mu_samples <- draws %>% select(starts_with("mu_pr"))

mu_transformed <- mu_samples %>%
  summarise(
    tau = median(pnorm(`mu_pr[1]`)) * 10,
    m_in = median(pnorm(`mu_pr[2]`)) * 10,
    m_out = median(pnorm(`mu_pr[3]`)) * 10,
    bias = median(pnorm(`mu_pr[4]`)),
    lambda = median(pnorm(`mu_pr[5]`)) * 5
  )

print("Group-level Transformed Medians:")
print(mu_transformed)

# 3. Individual-level parameter extraction for correlations
sub_params <- draws %>%
  select(starts_with("tau["), starts_with("m_in["), starts_with("m_out["), 
         starts_with("bias["), starts_with("lambda[")) %>%
  summarise(across(everything(), median)) %>%
  pivot_longer(everything(), names_to = "param", values_to = "median") %>%
  extract(param, into = c("parameter", "index"), regex = "(.*)\\[(\\d+)\\]", convert = TRUE) %>%
  pivot_wider(names_from = parameter, values_from = median)

fulldf <- read.csv("Study 1/Cleaning/output/fullTest.csv")
uIds <- sort(unique(fulldf$subID))
sub_params$subID <- uIds[sub_params$index]

indDiffs <- fulldf %>%
  filter(!duplicated(subID)) %>%
  select(subID, SI, Proto, RSE, NFC, NTB, DS, SING.Ind, SING.Inter)

combined_s1 <- sub_params %>% left_join(indDiffs, by = "subID")

# 4. Correlations
cor_results <- combined_s1 %>%
  select(m_in, m_out, lambda, bias, SI, Proto, RSE) %>%
  cor(use = "pairwise.complete.obs")

print("Correlation Matrix:")
print(cor_results)

saveRDS(combined_s1, "study1_asym_lambda_params_inddiffs.rds")
