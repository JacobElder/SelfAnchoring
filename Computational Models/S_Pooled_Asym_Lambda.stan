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
  // Global hyperparameters (Across Studies)
  vector[5] global_mu_pr;
  vector<lower=0>[5] global_sigma;

  // Study-level raw parameters
  array[nStudies] vector[5] study_mu_pr;
  vector<lower=0>[5] study_sigma;

  // Subject-level raw parameters (Matt trick)
  vector[nSubjects] tau_pr;     
  vector[nSubjects] m_in_pr;    
  vector[nSubjects] m_out_pr;   
  vector[nSubjects] bias_pr;    
  vector[nSubjects] lambda_pr;  
}

transformed parameters {
  vector<lower=0, upper=10>[nSubjects] tau;
  vector<lower=0, upper=10>[nSubjects] m_in;
  vector<lower=0, upper=10>[nSubjects] m_out;
  vector<lower=0, upper=1>[nSubjects] bias;
  vector<lower=0, upper=5>[nSubjects] lambda;

  for (i in 1:nSubjects) {
    int s = subjStudy[i];
    tau[i]    = Phi_approx(study_mu_pr[s, 1] + study_sigma[1] * tau_pr[i]) * 10; 
    m_in[i]   = Phi_approx(study_mu_pr[s, 2] + study_sigma[2] * m_in_pr[i]) * 10; 
    m_out[i]  = Phi_approx(study_mu_pr[s, 3] + study_sigma[3] * m_out_pr[i]) * 10; 
    bias[i]   = Phi_approx(study_mu_pr[s, 4] + study_sigma[4] * bias_pr[i]); 
    lambda[i] = Phi_approx(study_mu_pr[s, 5] + study_sigma[5] * lambda_pr[i]) * 5;
  }
}

model {
  // Global Priors
  global_mu_pr ~ normal(0, 1);
  global_sigma ~ normal(0, 0.3);

  // Study-level Priors (Partial Pooling)
  for (s in 1:nStudies) {
    study_mu_pr[s] ~ normal(global_mu_pr, global_sigma);
  }
  study_sigma ~ normal(0, 0.3);

  // Subject-level priors
  tau_pr    ~ normal(0, 1);
  m_in_pr   ~ normal(0, 1);
  m_out_pr  ~ normal(0, 1);
  bias_pr   ~ normal(0, 1);
  lambda_pr ~ normal(0, 1);

  for (s in 1:nSubjects) {
    vector[nTrain[s]] GPin;
    vector[nTrain[s]] GPout;
    matrix[nTrials[s], nTrain[s]] PS;
    vector[nTrials[s]] simW_in;
    vector[nTrials[s]] simW_out;
    vector[nTrials[s]] logit_p;

    GPin[1:nTrain[s]]  = inv_logit(m_in[s] * (prevSelf[s, 1:nTrain[s]] - 4));
    GPout[1:nTrain[s]] = inv_logit(-m_out[s] * (prevSelf[s, 1:nTrain[s]] - 4));
    
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
  array[nStudies] real mu_tau;
  array[nStudies] real mu_m_in;
  array[nStudies] real mu_m_out;
  array[nStudies] real mu_bias;
  array[nStudies] real mu_lambda;

  for (s in 1:nStudies) {
    mu_tau[s]    = Phi_approx(study_mu_pr[s, 1]) * 10;
    mu_m_in[s]   = Phi_approx(study_mu_pr[s, 2]) * 10;
    mu_m_out[s]  = Phi_approx(study_mu_pr[s, 3]) * 10;
    mu_bias[s]   = Phi_approx(study_mu_pr[s, 4]);
    mu_lambda[s] = Phi_approx(study_mu_pr[s, 5]) * 5;
  }
}
