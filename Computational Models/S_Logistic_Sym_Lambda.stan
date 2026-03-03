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
  vector[nSubjects] lambda_pr;  
}

transformed parameters {
  vector<lower=0, upper=10>[nSubjects] tau;
  vector<lower=0, upper=10>[nSubjects] m;
  vector<lower=0, upper=1>[nSubjects] bias;
  vector<lower=0, upper=5>[nSubjects] lambda;

  for (i in 1:nSubjects) {
    tau[i]   = Phi_approx(mu_pr[1] + sigma[1] * tau_pr[i]) * 10; 
    m[i]     = Phi_approx(mu_pr[2] + sigma[2] * m_pr[i]) * 10; 
    bias[i]  = Phi_approx(mu_pr[3] + sigma[3] * bias_pr[i]); 
    lambda[i] = Phi_approx(mu_pr[4] + sigma[4] * lambda_pr[i]) * 5;
  }
}

model {
  mu_pr  ~ normal(0, 1);
  sigma  ~ normal(0, 0.3);
  tau_pr    ~ normal(0, 1);
  m_pr      ~ normal(0, 1);
  bias_pr   ~ normal(0, 1);
  lambda_pr ~ normal(0, 1);

  for (s in 1:nSubjects) {
    vector[nTrain[s]] GPin;
    vector[nTrain[s]] GPout;
    matrix[nTrials[s], nTrain[s]] PS;
    vector[nTrials[s]] simW_in;
    vector[nTrials[s]] simW_out;
    vector[nTrials[s]] logit_p;

    GPin[1:nTrain[s]]  = inv_logit(m[s] * (prevSelf[s, 1:nTrain[s]] - 4));
    GPout[1:nTrain[s]] = inv_logit(-m[s] * (prevSelf[s, 1:nTrain[s]] - 4));
    
    for (t in 1:nTrials[s]) {
      PS[t] = pow(prevSim[s, t, 1:nTrain[s]], lambda[s]);
    }
    
    simW_in = PS * GPin[1:nTrain[s]] + 1e-9;
    simW_out = PS * GPout[1:nTrain[s]] + 1e-9;
    
    logit_p = log(bias[s]) - log(1-bias[s]) + tau[s] * (log(simW_in) - log(simW_out));

    for (t in 1:nTrials[s]) {
      if (groupChoice[s, t] > 0) {
        (groupChoice[s, t] - 1) ~ bernoulli_logit(logit_p[t]);
      }
    }
  }    
}

generated quantities {
  vector[nSubjects] log_lik;
  for (s in 1:nSubjects) {
    vector[nTrain[s]] GPin;
    vector[nTrain[s]] GPout;
    matrix[nTrials[s], nTrain[s]] PS;
    vector[nTrials[s]] simW_in;
    vector[nTrials[s]] simW_out;
    vector[nTrials[s]] logit_p;
    
    log_lik[s] = 0;
    
    GPin[1:nTrain[s]]  = inv_logit(m[s] * (prevSelf[s, 1:nTrain[s]] - 4));
    GPout[1:nTrain[s]] = inv_logit(-m[s] * (prevSelf[s, 1:nTrain[s]] - 4));
    
    for (t in 1:nTrials[s]) {
      PS[t] = pow(prevSim[s, t, 1:nTrain[s]], lambda[s]);
    }
    
    simW_in = PS * GPin[1:nTrain[s]] + 1e-9;
    simW_out = PS * GPout[1:nTrain[s]] + 1e-9;
    logit_p = log(bias[s]) - log(1-bias[s]) + tau[s] * (log(simW_in) - log(simW_out));
    
    for (t in 1:nTrials[s]) {
      if (groupChoice[s, t] > 0) {
        log_lik[s] += bernoulli_logit_lpmf(groupChoice[s, t] - 1 | logit_p[t]);
      }
    }
  }
}
