# General Parameter Recovery Template
# Adapted from existing PR_S_Logistic_1mOppose_Bias.R
# This script is designed to recover parameters for Asymmetric + Lambda models.

library(tidyverse)
library(cmdstanr)
library(parallel)

simulate_behavior <- function(self, trainIdx, testIdx, simMat, params, model_type = "asym_lambda") {
  # params: vector of [tau, m_in, m_out, bias, lambda, w]
  tau <- params[1]
  m_in <- params[2]
  m_out <- params[3]
  bias <- params[4]
  lambda <- params[5]
  w <- ifelse(length(params) > 5, params[6], 1.0) # default w=1 if not provided
  
  # 1. Transform Self-Beliefs
  GPin <- 1 / (1 + exp(-m_in * (self - 4)))
  GPout <- 1 / (1 + exp(m_out * (self - 4)))
  
  choices <- integer(length(testIdx))
  
  for (t in seq_along(testIdx)) {
    # 2. Generalization with Sensitivity (Lambda)
    PS <- simMat[trainIdx, testIdx[t]] ^ lambda
    
    # 3. Evidence Integration
    simW_in <- sum(GPin * PS)
    simW_out <- sum(GPout * PS)
    
    # 4. Choice Rule with Weight (w)
    # p(In) prop to bias * (simW_in * w + (1-w))^tau
    ev_in <- (simW_in * w + (1 - w)) ^ tau
    ev_out <- (simW_out * w + (1 - w)) ^ tau
    
    prob_in <- (bias * ev_in) / (bias * ev_in + (1 - bias) * ev_out)
    
    # Ensure no NAs
    if (is.na(prob_in)) prob_in <- 0.5
    
    choices[t] <- sample(c(1, 2), size = 1, prob = c(1 - prob_in, prob_in))
  }
  return(choices)
}

# --- RECOVERY WORKFLOW ---
# 1. Load Real Data/Similarity Matrix
# 2. Generate Synthetic Parameters (based on group-level fits)
# 3. Simulate choices for each synthetic subject
# 4. Fit Stan model to synthetic data
# 5. Correlate True vs. Recovered parameters

message("Parameter recovery template created. Update simulate_behavior() once final model is confirmed.")
