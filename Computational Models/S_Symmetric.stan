data {
  int<lower=1> nSubjects; 
  int<lower=1> maxTrials; 
  int<lower=1> maxTrain; 
  array[nSubjects] int<lower=1> nTrain; 
  array[nSubjects] int<lower=1> nTrials; 
  array[nSubjects, maxTrials] int<lower=0, upper=2> groupChoice; 
  array[nSubjects] matrix[maxTrials, maxTrain] prevSim; 
  array[nSubjects] vector[maxTrain] prevSelf; 
}

parameters {
  vector[4] mu_pr;
  vector<lower=0>[4] sigma;
  vector[nSubjects] tau_pr;     
  vector[nSubjects] m_pr;    
  vector[nSubjects] bias_pr;    
  vector[nSubjects] w_pr;
}

transformed parameters {
  vector<lower=0, upper=10>[nSubjects] tau;
  vector<lower=0, upper=10>[nSubjects] m;
  vector<lower=0, upper=1>[nSubjects] bias;
  vector<lower=0, upper=1>[nSubjects] w;

  for (i in 1:nSubjects) {
    tau[i]   = Phi_approx(mu_pr[1] + sigma[1] * tau_pr[i]) * 10; 
    m[i]     = Phi_approx(mu_pr[2] + sigma[2] * m_pr[i]) * 10; 
    bias[i]  = Phi_approx(mu_pr[3] + sigma[3] * bias_pr[i]); 
    w[i]     = Phi_approx(mu_pr[4] + sigma[4] * w_pr[i]);
  }
}

model {
  mu_pr  ~ normal(0, 1);
  sigma  ~ normal(0, 0.3);
  tau_pr    ~ normal(0, 1);
  m_pr      ~ normal(0, 1);
  bias_pr   ~ normal(0, 1);
  w_pr      ~ normal(0, 1);

  for (s in 1:nSubjects) {
    vector[nTrain[s]] GP;
    matrix[nTrials[s], nTrain[s]] PS;
    vector[nTrials[s]] simW_in;
    vector[nTrials[s]] simW_out;
    vector[nTrials[s]] logit_p;

    GP[1:nTrain[s]] = inv_logit(m[s] * (prevSelf[s, 1:nTrain[s]] - 4));
    PS = prevSim[s, 1:nTrials[s], 1:nTrain[s]];
    
    simW_in  = PS * GP[1:nTrain[s]] + 1e-9;
    simW_out = PS * (1.0 - GP[1:nTrain[s]]) + 1e-9;
    
    logit_p = log(bias[s] + 1e-9) - log(1.0 - bias[s] + 1e-9) + tau[s] * (log(simW_in) - log(simW_out));

    for (t in 1:nTrials[s]) {
      if (groupChoice[s, t] > 0) {
        target += log_sum_exp(log(w[s]) + bernoulli_lpmf(groupChoice[s, t] - 1 | 0.5),
                             log1m(w[s]) + bernoulli_logit_lpmf(groupChoice[s, t] - 1 | logit_p[t]));
      }
    }
  }
}
