data {
  int<lower=1> nSubjects; 
  int<lower=1> maxTrials; 
  array[nSubjects] int<lower=1> nTrials; 
  array[nSubjects, maxTrials] int<lower=0, upper=2> groupChoice; 
}

parameters {
  vector[1] mu_pr;
  vector<lower=0>[1] sigma;
  vector[nSubjects] bias_pr;    
}

transformed parameters {
  vector<lower=0, upper=1>[nSubjects] bias;

  for (i in 1:nSubjects) {
    bias[i] = Phi_approx(mu_pr[1] + sigma[1] * bias_pr[i]); 
  }
}

model {
  mu_pr  ~ normal(0, 1);
  sigma  ~ normal(0, 0.3);
  bias_pr ~ normal(0, 1);

  for (s in 1:nSubjects) {
    real logit_p = log(bias[s]) - log(1-bias[s]);
    for (t in 1:nTrials[s]) {
      if (groupChoice[s, t] > 0) {
        (groupChoice[s, t] - 1) ~ bernoulli_logit(logit_p);
      }
    }
  }    
}

generated quantities {
  vector[nSubjects] log_lik;
  for (s in 1:nSubjects) {
    real logit_p = log(bias[s]) - log(1-bias[s]);
    log_lik[s] = 0;
    for (t in 1:nTrials[s]) {
      if (groupChoice[s, t] > 0) {
        log_lik[s] += bernoulli_logit_lpmf(groupChoice[s, t] - 1 | logit_p);
      }
    }
  }
}
