data {
  int<lower=1> nSubjects; 
  int<lower=1> nStudies;
  array[nSubjects] int<lower=1, upper=nStudies> subjStudy;
  int<lower = 1> maxTrials; 
  int<lower = 1> maxTrain; 
  array[nSubjects] int<lower = 1> nTrain; 
  array[nSubjects] int<lower = 1> nTrials; 
  array[nSubjects, maxTrials] int<lower = 0, upper = 2> groupChoice; 
  array[nSubjects] matrix[maxTrials, maxTrain] prevSim; 
  array[nSubjects] vector[maxTrain] prevSelf; 
}

parameters {
  vector[2] global_mu_pr;
  vector<lower=0.01>[2] global_sigma;
  array[nStudies] vector[2] study_mu_raw;
  vector<lower=0.01>[2] study_sigma;
  vector[nSubjects] bias_pr;    
  vector[nSubjects] w_pr;
}

transformed parameters {
  vector<lower=0, upper=1>[nSubjects] bias;
  vector<lower=0, upper=1>[nSubjects] w;
  array[nStudies] vector[2] study_mu_pr;
  for (s in 1:nStudies) {
    study_mu_pr[s] = global_mu_pr + global_sigma .* study_mu_raw[s];
  }
  for (i in 1:nSubjects) {
    bias[i] = Phi_approx(study_mu_pr[subjStudy[i], 1] + study_sigma[1] * bias_pr[i]);
    w[i]    = Phi_approx(study_mu_pr[subjStudy[i], 2] + study_sigma[2] * w_pr[i]);
  }
}

model {
  global_mu_pr ~ normal(0, 1);
  global_sigma ~ cauchy(0, 1);
  for (s in 1:nStudies) study_mu_raw[s] ~ normal(0, 1);
  study_sigma ~ cauchy(0, 1);
  bias_pr ~ normal(0, 1);
  w_pr ~ normal(0, 1);

  for (s in 1:nSubjects) {
    for (t in 1:nTrials[s]) {
      if (groupChoice[s, t] > 0) {
        real p_model = bias[s];
        target += log_sum_exp(log(w[s]) + bernoulli_lpmf(groupChoice[s, t] - 1 | 0.5),
                             log1m(w[s]) + bernoulli_lpmf(groupChoice[s, t] - 1 | p_model));
      }
    }
  }
}

generated quantities {
  // ================================================================
  // GROUP-LEVEL PARAMETERS (back-transformed global means)
  //   fit$draws("mu_bias")   # ingroup response bias [0, 1]
  //   fit$draws("mu_w")      # lapse rate [0, 1]
  // ================================================================
  real<lower=0, upper=1> mu_bias = Phi_approx(global_mu_pr[1]);
  real<lower=0, upper=1> mu_w    = Phi_approx(global_mu_pr[2]);

  // ================================================================
  // TRIAL-LEVEL QUANTITIES (index = (s-1)*maxTrials + t)
  //   fit$draws("log_lik")   # log-likelihood for LOO-CV
  //   fit$draws("p_pred")    # P(ingroup) per trial (= bias[s]; no projection)
  //   Note: MCR not applicable — Bias model has no similarity-based projection
  // ================================================================
  vector[nSubjects * maxTrials] log_lik = rep_vector(0.0, nSubjects * maxTrials);
  vector[nSubjects * maxTrials] p_pred  = rep_vector(0.0, nSubjects * maxTrials);

  for (s in 1:nSubjects) {
    for (t in 1:nTrials[s]) {
      if (groupChoice[s, t] > 0) {
        real p_model = bias[s];
        log_lik[(s-1)*maxTrials + t] = log_sum_exp(log(w[s]) + bernoulli_lpmf(groupChoice[s, t] - 1 | 0.5),
                                                  log1m(w[s]) + bernoulli_lpmf(groupChoice[s, t] - 1 | p_model));
        p_pred[(s-1)*maxTrials + t]  = w[s] * 0.5 + (1.0 - w[s]) * p_model;
      }
    }
  }
}
