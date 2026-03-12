library(lme4)
library(simr)
library(tidyverse)
library(here)

# 1. Load Study 1 Data as the base structure
s1_df <- read.csv(here("Study 1/Cleaning/output/fullTest.csv")) %>%
  filter(!is.na(ingChoiceN))

s1_df$predicted.Z <- scale(s1_df$predicted)[,1]
s1_df$desirability.Z <- scale(s1_df$desirability)[,1]

# 2. Create a 3-condition structure
s1_ext <- bind_rows(
  s1_df %>% mutate(outgroup = "Negation", subID = subID),
  s1_df %>% mutate(outgroup = "High-Status", subID = subID + 1000),
  s1_df %>% mutate(outgroup = "Low-Status", subID = subID + 2000)
)
s1_ext$outgroup <- factor(s1_ext$outgroup, levels = c("Negation", "High-Status", "Low-Status"))

message("Fitting base model...")
m_base <- glmer(
  ingChoiceN ~ predicted.Z * outgroup + desirability.Z +
    (predicted.Z | subID) + (1 | trait),
  data = s1_ext, family = binomial, 
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
)

# 3. Set Effect Sizes for the Interaction
# We use the focal interaction effect sizes (0.35) 
fixed_effects <- fixef(m_base)
fixed_effects["predicted.Z:outgroupHigh-Status"] <- 0.35
fixed_effects["predicted.Z:outgroupLow-Status"]  <- -0.35
fixef(m_base) <- fixed_effects

# 4. Use Likelihood Ratio Test for the interaction term (2 degrees of freedom)
# This avoids the error with fixed() on a factor interaction.
message("Calculating a priori power for Interaction (Similarity x Condition) at N=183 via LR test...")
p_inter <- powerSim(m_base, test = fcompare(~ predicted.Z + outgroup + desirability.Z), nsim = 100)

print(p_inter)

# 5. Generate Power Curve
message("Generating power curve for interaction...")
pc <- powerCurve(m_base, test = fcompare(~ predicted.Z + outgroup + desirability.Z), 
                 along = "subID", breaks = c(60, 90, 120, 150, 183), nsim = 100)

print(pc)

# Save curve data
curve_data <- data.frame(
  n = pc$x,
  power = pc$ps
)
write.csv(curve_data, here("Results", "power_curve_s2_a_priori_interaction.csv"), row.names = FALSE)
