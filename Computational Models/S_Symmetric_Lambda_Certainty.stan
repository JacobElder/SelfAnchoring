// Symmetric + Lambda + Certainty Model for Self-Anchoring
// This model adds a 'gamma_cert' parameter to the Symmetric+Lambda architecture.
// 'gamma_cert' scales the evidence from similarity, capturing individual differences
// in how decisively a person uses their self-concept to make group judgments.
// gamma_cert > 1: Over-weights evidence (more decisive)
// gamma_cert < 1: Under-weights evidence (less decisive)

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
  // Group-level parameters (probit-scale)
  vector[5] mu_pr;      // m, bias, lambda, w, gamma_cert
  vector<lower=0>[5] sigma;

  // Subject-level parameters (raw, non-centered)
  vector[nSubjects] m_pr;
  vector[nSubjects] bias_pr;    
  vector[nSubjects] lambda_pr;
  vector[nSubjects] w_pr;
  vector[nSubjects] gamma_cert_pr;
}

transformed parameters {
  // Subject-level parameters (transformed to constrained scale)
  vector<lower=0, upper=10>[nSubjects] m;
  vector<lower=0, upper=1>[nSubjects] bias;
  vector<lower=0, upper=5>[nSubjects] lambda;
  vector<lower=0, upper=1>[nSubjects] w;
  vector<lower=0>[nSubjects] gamma_cert; // Certainty parameter (lognormal prior implies > 0)

  for (i in 1:nSubjects) {
    m[i]          = Phi_approx(mu_pr[1] + sigma[1] * m_pr[i]) * 10;
    bias[i]       = Phi_approx(mu_pr[2] + sigma[2] * bias_pr[i]);
    lambda[i]     = Phi_approx(mu_pr[3] + sigma[3] * lambda_pr[i]) * 5;
    w[i]          = Phi_approx(mu_pr[4] + sigma[4] * w_pr[i]);
    gamma_cert[i] = exp(mu_pr[5] + sigma[5] * gamma_cert_pr[i]); // Log-normal
  }
}

model {
  // Priors
  mu_pr  ~ normal(0, 1);
  // Set a slightly more informative prior for gamma_cert's mean to center it near 1.
  mu_pr[5] ~ normal(0, 0.1); 
  sigma  ~ normal(0, 0.3);
  
  m_pr          ~ normal(0, 1);
  bias_pr       ~ normal(0, 1);
  lambda_pr     ~ normal(0, 1);
  w_pr          ~ normal(0, 1);
  gamma_cert_pr ~ normal(0, 1);

  // Likelihood
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
    
    // The 'gamma_cert' parameter scales the evidence from similarity.
    logit_p = logit(bias[s]) + gamma_cert[s] * (log(simW_in) - log(simW_out));

    for (t in 1:nTrials[s]) {
      if (groupChoice[s, t] > 0) {
        target += log_mix(w[s], 
                         bernoulli_lpmf(groupChoice[s, t] - 1 | 0.5),
                         bernoulli_logit_lpmf(groupChoice[s, t] - 1 | logit_p[t]));
      }
    }
  }
}

generated quantities {
  // Group-level means on natural scale
  real<lower=0, upper=10> mu_m      = Phi_approx(mu_pr[1]) * 10;
  real<lower=0, upper=1>  mu_bias   = Phi_approx(mu_pr[2]);
  real<lower=0, upper=5>  mu_lambda = Phi_approx(mu_pr[3]);
  real<lower=0, upper=1>  mu_w      = Phi_approx(mu_pr[4]);
  real<lower=0>           mu_gamma_cert = exp(mu_pr[5]);

  // Trial-level quantities for posterior predictive checks and LOO
  vector[nSubjects * maxTrials] log_lik = rep_vector(0.0, nSubjects * maxTrials);
  
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

    logit_p = logit(bias[s]) + gamma_cert[s] * (log(simW_in) - log(simW_out));

    for (t in 1:nTrials[s]) {
      if (groupChoice[s, t] > 0) {
        log_lik[(s-1)*maxTrials + t] = log_mix(w[s], 
                                              bernoulli_lpmf(groupChoice[s, t] - 1 | 0.5),
                                              bernoulli_logit_lpmf(groupChoice[s, t] - 1 | logit_p[t]));
      }
    }
  }
}
