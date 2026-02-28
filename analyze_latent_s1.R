# Latent Parameter Analysis for Study 1 (Minimal Groups)
library(cmdstanr)
library(tidyverse)
library(marginaleffects)

# 1. Load Study 1 Winner
file_path <- "fit_s1_symmetric.rds"
fit_s1 <- readRDS(file_path)

draws <- fit_s1$draws(format = "df")

# 2. Extract Latent Parameters (Posterior Medians)
sub_params <- draws %>%
  select(starts_with("tau["), starts_with("m["), starts_with("bias[")) %>%
  summarise(across(everything(), median)) %>%
  pivot_longer(everything(), names_to = "param", values_to = "median") %>%
  extract(param, into = c("parameter", "index"), regex = "(.*)\\[(\\d+)\\]", convert = TRUE) %>%
  pivot_wider(names_from = parameter, values_from = median)

# 3. Load Study 1 Individual Differences
fulldf <- read.csv("Study 1/Cleaning/output/fullTest.csv")
uIds <- sort(unique(fulldf$subID))
sub_params$subID <- uIds[sub_params$index]

indDiffs <- fulldf %>%
  filter(!duplicated(subID)) %>%
  select(subID, SI, Proto, RSE, NFC, NTB, DS, SING.Ind, SING.Inter)

# Merge
combined_s1 <- sub_params %>%
  left_join(indDiffs, by = "subID")

# 4. Associations with m (Projection Rate)
message("Study 1: Association between m (Projection Rate) and Social Identification")
m_si <- lm(scale(m) ~ scale(SI), data = combined_s1)
print(summary(m_si))

message("Study 1: Association between m (Projection Rate) and Self-Prototypicality")
m_proto <- lm(scale(m) ~ scale(Proto), data = combined_s1)
print(summary(m_proto))

# 5. Latent Marginal Effects
message("Calculating latent marginal effects (AME of SI on m)...")
print(avg_slopes(m_si))

saveRDS(combined_s1, "study1_latent_params_inddiffs.rds")
message("Done.")
