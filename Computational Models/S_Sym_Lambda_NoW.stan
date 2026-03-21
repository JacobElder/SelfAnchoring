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
  vector[3] mu_pr;
  vector<lower=0>[3] sigma;
  vector[nSubjects] m_pr;
  vector[nSubjects] bias_pr;
  vector[nSubjects] lambda_pr;
}

transformed parameters {
  vector<lower=0, upper=10>[nSubjects] m;
  vector<lower=0, upper=1>[nSubjects] bias;
  vector<lower=0, upper=5>[nSubjects] lambda;

  for (i in 1:nSubjects) {
    m[i]      = Phi_approx(mu_pr[1] + sigma[1] * m_pr[i]) * 10;
    bias[i]   = Phi_approx(mu_pr[2] + sigma[2] * bias_pr[i]);
    lambda[i] = Phi_approx(mu_pr[3] + sigma[3] * lambda_pr[i]) * 5;
  }
}

model {
  mu_pr  ~ normal(0, 1);
  sigma[1] ~ normal(0, 0.3);  // m
  sigma[2] ~ normal(0, 0.3);  // bias
  sigma[3] ~ normal(0, 0.1);  // lambda: tighter prior to reduce m-lambda confound
  m_pr      ~ normal(0, 1);
  bias_pr   ~ normal(0, 1);
  lambda_pr ~ normal(0, 1);

  for (s in 1:nSubjects) {
    vector[nTrain[s]] GP;
    matrix[nTrials[s], nTrain[s]] PS;
    vector[nTrials[s]] simW_in;
    vector[nTrials[s]] simW_out;
    vector[nTrials[s]] logit_p;

    GP[1:nTrain[s]] = inv_logit(m[s] * (prevSelf[s, 1:nTrain[s]] - 4.0));
    PS = pow(prevSim[s, 1:nTrials[s], 1:nTrain[s]], lambda[s]);

    simW_in  = PS * GP[1:nTrain[s]] + 1e-9;
    simW_out = PS * (1.0 - GP[1:nTrain[s]]) + 1e-9;

    logit_p = log(bias[s] + 1e-9) - log(1.0 - bias[s] + 1e-9) + (log(simW_in) - log(simW_out));

    for (t in 1:nTrials[s]) {
      if (groupChoice[s, t] > 0) {
        target += bernoulli_logit_lpmf(groupChoice[s, t] - 1 | logit_p[t]);
      }
    }
  }
}

generated quantities {
  real<lower=0, upper=10> mu_m      = Phi_approx(mu_pr[1]) * 10;
  real<lower=0, upper=1>  mu_bias   = Phi_approx(mu_pr[2]);
  real<lower=0, upper=5>  mu_lambda = Phi_approx(mu_pr[3]) * 5;

  vector[nSubjects * maxTrials] log_lik = rep_vector(0.0, nSubjects * maxTrials);
  vector[nSubjects * maxTrials] p_pred  = rep_vector(0.0, nSubjects * maxTrials);
  vector[nSubjects * maxTrials] mcr     = rep_vector(0.0, nSubjects * maxTrials);

  for (s in 1:nSubjects) {
    vector[nTrain[s]] GP;
    matrix[nTrials[s], nTrain[s]] PS;
    vector[nTrials[s]] simW_in;
    vector[nTrials[s]] simW_out;
    vector[nTrials[s]] logit_p;

    GP[1:nTrain[s]] = inv_logit(m[s] * (prevSelf[s, 1:nTrain[s]] - 4.0));
    PS = pow(prevSim[s, 1:nTrials[s], 1:nTrain[s]], lambda[s]);

    simW_in  = PS * GP[1:nTrain[s]] + 1e-9;
    simW_out = PS * (1.0 - GP[1:nTrain[s]]) + 1e-9;

    logit_p = log(bias[s] + 1e-9) - log(1.0 - bias[s] + 1e-9) + (log(simW_in) - log(simW_out));

    for (t in 1:nTrials[s]) {
      if (groupChoice[s, t] > 0) {
        log_lik[(s-1)*maxTrials + t] = bernoulli_logit_lpmf(groupChoice[s, t] - 1 | logit_p[t]);
        p_pred[(s-1)*maxTrials + t]  = inv_logit(logit_p[t]);
        mcr[(s-1)*maxTrials + t]     = simW_in[t] / simW_out[t];
      }
    }
  }
}
