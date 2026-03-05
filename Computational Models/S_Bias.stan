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
    for (t in 1:nTrials[s]) {
      if (groupChoice[s, t] > 0) {
        real p_model = bias[s];
        target += log_sum_exp(log(w[s]) + bernoulli_lpmf(groupChoice[s, t] - 1 | 0.5),
                             log1m(w[s]) + bernoulli_lpmf(groupChoice[s, t] - 1 | p_model));
      }
    }
  }
}
