data {
  int<lower=1> nSubjects; 
  int<lower=1> maxTrials; 
  array[nSubjects] int<lower=1> nTrials; 
  array[nSubjects, maxTrials] int<lower=0, upper=2> groupChoice; 
}

parameters {
  vector[2] mu_pr;
  vector<lower=0>[2] sigma;
  vector[nSubjects] bias_pr;    
  vector[nSubjects] w_pr;
}

transformed parameters {
  vector<lower=0, upper=1>[nSubjects] bias;
  vector<lower=0, upper=1>[nSubjects] w;

  for (i in 1:nSubjects) {
    bias[i] = Phi_approx(mu_pr[1] + sigma[1] * bias_pr[i]);
    w[i]    = Phi_approx(mu_pr[2] + sigma[2] * w_pr[i]);
  }
}

model {
  mu_pr  ~ normal(0, 1);
  sigma  ~ normal(0, 0.3);
  bias_pr ~ normal(0, 1);
  w_pr    ~ normal(0, 1);

  for (s in 1:nSubjects) {
    // Pre-calculate subject-level mixture likelihoods
    real lp_mix_ingroup  = log_mix(w[s], bernoulli_lpmf(0 | 0.5), bernoulli_lpmf(0 | bias[s]));
    real lp_mix_outgroup = log_mix(w[s], bernoulli_lpmf(1 | 0.5), bernoulli_lpmf(1 | bias[s]));
    
    for (t in 1:nTrials[s]) {
      if (groupChoice[s, t] == 1) {
        target += lp_mix_ingroup;
      } else if (groupChoice[s, t] == 2) {
        target += lp_mix_outgroup;
      }
    }
  }
}

generated quantities {
  real<lower=0, upper=1> mu_bias = Phi_approx(mu_pr[1]);
  real<lower=0, upper=1> mu_w    = Phi_approx(mu_pr[2]);

  vector[nSubjects * maxTrials] log_lik = rep_vector(0.0, nSubjects * maxTrials);
  vector[nSubjects * maxTrials] p_pred  = rep_vector(0.0, nSubjects * maxTrials);

  for (s in 1:nSubjects) {
    real lp_mix_ingroup  = log_mix(w[s], bernoulli_lpmf(0 | 0.5), bernoulli_lpmf(0 | bias[s]));
    real lp_mix_outgroup = log_mix(w[s], bernoulli_lpmf(1 | 0.5), bernoulli_lpmf(1 | bias[s]));
    real p_obs = w[s] * 0.5 + (1.0 - w[s]) * bias[s];
    
    for (t in 1:nTrials[s]) {
      if (groupChoice[s, t] == 1) {
        log_lik[(s-1)*maxTrials + t] = lp_mix_ingroup;
        p_pred[(s-1)*maxTrials + t]  = p_obs;
      } else if (groupChoice[s, t] == 2) {
        log_lik[(s-1)*maxTrials + t] = lp_mix_outgroup;
        p_pred[(s-1)*maxTrials + t]  = p_obs;
      }
    }
  }
}
