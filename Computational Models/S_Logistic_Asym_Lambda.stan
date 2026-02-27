data {
  int<lower=1> nSubjects; // number of subjects total
  int<lower = 1> maxTrials; //number of testing/generalization trials
  int<lower = 1> maxTrain; // max number of training trials
  int<lower = 1> nTrain[nSubjects]; // per participant number of training trials
  int<lower = 1> nTrials[nSubjects]; //per participant number of testing/generalization trials
  int<lower = 0, upper = 2> groupChoice[nSubjects,maxTrials]; // which group chosen
  
  vector[maxTrain] prevSim[nSubjects, maxTrials]; // similarities from training to testing traits
  vector[maxTrain] prevSelf[nSubjects]; // matrix of training self-evaluations
}

parameters {
  // Hyper(group)-parameters
  vector[5] mu_pr;
  vector<lower=0>[5] sigma;

  // Subject-level raw parameters (Matt trick)
  vector[nSubjects] tau_pr;     // temperature
  vector[nSubjects] m_in_pr;    // projection rate (ingroup)
  vector[nSubjects] m_out_pr;   // repulsion rate (outgroup)
  vector[nSubjects] bias_pr;    // classification bias
  vector[nSubjects] lambda_pr;  // generalization sensitivity
}

transformed parameters {
  vector<lower=0, upper=10>[nSubjects] tau;
  vector<lower=0, upper=10>[nSubjects] m_in;
  vector<lower=0, upper=10>[nSubjects] m_out;
  vector<lower=0, upper=1>[nSubjects] bias;
  vector<lower=0, upper=5>[nSubjects] lambda;

  for (i in 1:nSubjects) {
    tau[i]    = Phi_approx(mu_pr[1] + sigma[1] * tau_pr[i]) * 10; 
    m_in[i]   = Phi_approx(mu_pr[2] + sigma[2] * m_in_pr[i]) * 10; 
    m_out[i]  = Phi_approx(mu_pr[3] + sigma[3] * m_out_pr[i]) * 10; 
    bias[i]   = Phi_approx(mu_pr[4] + sigma[4] * bias_pr[i]); 
    lambda[i] = Phi_approx(mu_pr[5] + sigma[5] * lambda_pr[i]) * 5;
  }
}

model {
  // Hyperparameters
  mu_pr  ~ normal(0, 1);
  sigma  ~ normal(0, 0.3);

  // Individual-level priors
  tau_pr    ~ normal(0, 1);
  m_in_pr   ~ normal(0, 1);
  m_out_pr  ~ normal(0, 1);
  bias_pr   ~ normal(0, 1);
  lambda_pr ~ normal(0, 1);

  for (s in 1:nSubjects) {
    vector[2] simW;
    vector[2] prob;
    vector[nTrain[s]] GPin;
    vector[nTrain[s]] GPout;
    vector[nTrain[s]] PS;

    // Decoupled projection and repulsion rates
    GPin[1:nTrain[s]]  = rep_vector(1, nTrain[s]) ./ (1 + exp((-m_in[s])  * (prevSelf[s, 1:nTrain[s]] - 4)));
    GPout[1:nTrain[s]] = rep_vector(1, nTrain[s]) ./ (1 + exp((m_out[s]) * (prevSelf[s, 1:nTrain[s]] - 4)));
    
    for (t in 1:nTrials[s]) {
      // Generalization sensitivity applied to semantic similarity
      for (i in 1:nTrain[s]) {
        PS[i] = pow(prevSim[s, t, i], lambda[s]);
      }
      
      simW[1] = dot_product(GPout[1:nTrain[s]], PS); // Evidence for outgroup
      simW[2] = dot_product(GPin[1:nTrain[s]], PS);  // Evidence for ingroup
      
      // Categorical choice using evidence weighted by bias and temperature
      prob[1] = ((1 - bias[s]) * pow(simW[1], tau[s])) / (((1 - bias[s]) * pow(simW[1], tau[s])) + (bias[s] * pow(simW[2], tau[s])));
      prob[2] = (bias[s] * pow(simW[2], tau[s])) / (((1 - bias[s]) * pow(simW[1], tau[s])) + (bias[s] * pow(simW[2], tau[s])));
      
      groupChoice[s, t] ~ categorical(prob);
    }
  }    
}

generated quantities {
  real<lower=0, upper=10> mu_tau;
  real<lower=0, upper=10> mu_m_in;
  real<lower=0, upper=10> mu_m_out;
  real<lower=0, upper=1>  mu_bias;
  real<lower=0, upper=5>  mu_lambda;

  real log_lik[nSubjects];
  real y_pred[nSubjects, maxTrials];

  for (i in 1:nSubjects) {
    for (t in 1:maxTrials) {
      y_pred[i, t] = -1;
    }
  }

  mu_tau    = Phi_approx(mu_pr[1]) * 10;
  mu_m_in   = Phi_approx(mu_pr[2]) * 10;
  mu_m_out  = Phi_approx(mu_pr[3]) * 10;
  mu_bias   = Phi_approx(mu_pr[4]);
  mu_lambda = Phi_approx(mu_pr[5]) * 5;

  {
    for (s in 1:nSubjects) {
      vector[2] simW;
      vector[2] prob;
      vector[nTrain[s]] GPin;
      vector[nTrain[s]] GPout;
      vector[nTrain[s]] PS;
      
      log_lik[s] = 0;
      
      GPin[1:nTrain[s]]  = rep_vector(1, nTrain[s]) ./ (1 + exp((-m_in[s])  * (prevSelf[s, 1:nTrain[s]] - 4)));
      GPout[1:nTrain[s]] = rep_vector(1, nTrain[s]) ./ (1 + exp((m_out[s]) * (prevSelf[s, 1:nTrain[s]] - 4)));
      
      for (t in 1:nTrials[s]) {
        for (i in 1:nTrain[s]) {
          PS[i] = pow(prevSim[s, t, i], lambda[s]);
        }
        
        simW[1] = dot_product(GPout[1:nTrain[s]], PS);
        simW[2] = dot_product(GPin[1:nTrain[s]], PS);
        
        prob[1] = ((1 - bias[s]) * pow(simW[1], tau[s])) / (((1 - bias[s]) * pow(simW[1], tau[s])) + (bias[s] * pow(simW[2], tau[s])));
        prob[2] = (bias[s] * pow(simW[2], tau[s])) / (((1 - bias[s]) * pow(simW[1], tau[s])) + (bias[s] * pow(simW[2], tau[s])));
          
        log_lik[s] += categorical_lpmf(groupChoice[s, t] | prob);
        y_pred[s, t] = categorical_rng(prob);
      }
    }   
  }
}
