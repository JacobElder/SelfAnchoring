data {
  int<lower=1> nSubjects; // number of subjects total
  int<lower = 1> maxTrials; //number of testing/generalization trials
  int<lower = 1> maxTrain; // max number of training trials
  array[nSubjects] int<lower = 1> nTrain; // per participant number of training trials
  array[nSubjects] int<lower = 1> nTrials; //per participant number of testing/generalization trials
  array[nSubjects, maxTrials] int<lower = 0, upper = 2> groupChoice; // which group chosen
  
  array[nSubjects] matrix[maxTrials, maxTrain] prevSim; // Matrix of subject's training similarities to testing traits
  array[nSubjects] vector[maxTrain] prevSelf; // matrix of training self-evaluations
}

parameters {
  // Hyper(group)-parameters
  vector[3] mu_pr;
  vector<lower=0>[3] sigma;

  // Subject-level raw parameters (for Matt trick)
  vector[nSubjects] tau_pr;  // inverse temperature
  vector[nSubjects] m_pr;  // slope for ingroup
  vector[nSubjects] bias_pr;  // slope for bias
}

transformed parameters {
  // subject-level parameters
  vector<lower=0, upper=10>[nSubjects] tau; 
  vector<lower=0, upper=10>[nSubjects] m; 
  vector<lower=0, upper=1>[nSubjects] bias; 

  for (i in 1:nSubjects) {
    tau[i] = Phi_approx(mu_pr[1] + sigma[1] * tau_pr[i]) * 10; 
    m[i] = Phi_approx(mu_pr[2] + sigma[2] * m_pr[i]) * 10; 
    bias[i] = Phi_approx(mu_pr[3] + sigma[3] * bias_pr[i]); 
  }
}

model {
  // Hyperparameters
  mu_pr  ~ normal(0, 1);
  sigma ~ normal(0, .3);

  // individual parameters
  tau_pr ~ normal(0, 1);
  m_pr ~ normal(0, 2); 
  bias_pr ~ normal(0, 1);

  for (s in 1:nSubjects) {
    vector[nTrain[s]] GPin;
    vector[nTrain[s]] GPout;
    vector[nTrials[s]] simW_in;
    vector[nTrials[s]] simW_out;
    vector[nTrials[s]] logit_p;

    // Logistic transformation of eval ratings
    GPin[1:nTrain[s]] = inv_logit(m[s] * (prevSelf[s, 1:nTrain[s]] - 4));
    GPout[1:nTrain[s]] = inv_logit(-m[s] * (prevSelf[s, 1:nTrain[s]] - 4));
    
    // Matrix multiplication for all trials at once
    simW_in = prevSim[s, 1:nTrials[s], 1:nTrain[s]] * GPin[1:nTrain[s]] + 1e-9;
    simW_out = prevSim[s, 1:nTrials[s], 1:nTrain[s]] * GPout[1:nTrain[s]] + 1e-9;
    
    // Logit for Ingroup choice (Choice 2)
    logit_p = log(bias[s]) - log(1-bias[s]) + tau[s] * (log(simW_in) - log(simW_out));
    
    for (t in 1:nTrials[s]) {
      if (groupChoice[s, t] > 0) {
        (groupChoice[s, t] - 1) ~ bernoulli_logit(logit_p[t]);
      }
    }
  }    
}

generated quantities {
  real<lower=0, upper=10> mu_tau = Phi_approx(mu_pr[1]) * 10;
  real<lower=0, upper=10> mu_m = Phi_approx(mu_pr[2]) * 10;
  real<lower=0, upper=1> mu_bias = Phi_approx(mu_pr[3]);

  vector[nSubjects] log_lik;
  
  for (s in 1:nSubjects) {
    vector[nTrain[s]] GPin;
    vector[nTrain[s]] GPout;
    vector[nTrials[s]] simW_in;
    vector[nTrials[s]] simW_out;
    vector[nTrials[s]] logit_p;
    
    log_lik[s] = 0;
    
    GPin[1:nTrain[s]] = inv_logit(m[s] * (prevSelf[s, 1:nTrain[s]] - 4));
    GPout[1:nTrain[s]] = inv_logit(-m[s] * (prevSelf[s, 1:nTrain[s]] - 4));
    
    simW_in = prevSim[s, 1:nTrials[s], 1:nTrain[s]] * GPin[1:nTrain[s]] + 1e-9;
    simW_out = prevSim[s, 1:nTrials[s], 1:nTrain[s]] * GPout[1:nTrain[s]] + 1e-9;
    logit_p = log(bias[s]) - log(1-bias[s]) + tau[s] * (log(simW_in) - log(simW_out));
    
    for (t in 1:nTrials[s]) {
       if (groupChoice[s, t] > 0) {
          log_lik[s] += bernoulli_logit_lpmf(groupChoice[s, t] - 1 | logit_p[t]);
       }
    }
  }
}
