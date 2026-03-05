data {
  int<lower=1> nSubjects;
  int<lower = 1> maxTrials; //number of trials
  int<lower = 1> nTrials[nSubjects]; //number of trials per sub
  int<lower = 0, upper = 2> groupChoice[nSubjects,maxTrials]; //index of which side chosen
  real sg[nSubjects, maxTrials];
}

// transformed data {
//   vector[296] initV;  // initial values for V
//   initV = rep_vector(0.5, 296);
// }

parameters {
  // Declare all parameters as vectors for vectorizing
  // Hyper(group)-parameters
  vector[5] mu_pr;
  vector<lower=0>[5] sigma;

  // Subject-level raw parameters (for Matt trick)
  //vector[nSubjects] A_pr;    // learning rate
  vector[nSubjects] tau_pr;  // inverse temperature
  vector[nSubjects] SAi_pr;  // inverse temperature
  vector[nSubjects] bias_pr;
  vector[nSubjects] SAo_pr;
  vector[nSubjects] expP_pr;
  
}

transformed parameters {
  // subject-level parameters
  vector<lower=0, upper=10>[nSubjects] tau;
  vector<lower=0, upper=2>[nSubjects] SAi;
  vector<lower=0, upper=2>[nSubjects] bias;
  vector<lower=0, upper=2>[nSubjects] SAo;
  vector<lower=0, upper=10>[nSubjects] expP;

  for (i in 1:nSubjects) {
    tau[i] = Phi_approx(mu_pr[1] + sigma[1] * tau_pr[i]) * 10;
    SAi[i]   = Phi_approx(mu_pr[2]  + sigma[2]  * SAi_pr[i]) * 2;
    bias[i]   = Phi_approx(mu_pr[3]  + sigma[3]  * bias_pr[i]) * 2;
    SAo[i]   = Phi_approx(mu_pr[4]  + sigma[4]  * SAo_pr[i]) * 2;
    expP[i]   = Phi_approx(mu_pr[5]  + sigma[5]  * expP_pr[i]) * 10;
    
  }
}

model {
  // Hyperparameters
  mu_pr  ~ normal(0, 1);
  sigma ~ normal(0, 0.2);

  // individual parameters
  tau_pr ~ normal(0, 1);
  SAi_pr ~ normal(0, 1);
  bias_pr ~ normal(0, 1);
  SAo_pr ~ normal(0, 1);
  expP_pr ~ normal(0, 1);

  for (s in 1:nSubjects) {

    vector[2] simW;

    for (t in 1:nTrials[s]) {
      
      simW[1] = pow( (sg[s,t] * (SAo[s] - 1)), expP[s]);
      simW[2] = pow( (sg[s,t] * (SAi[s] - 1)), expP[s]);
      simW[2] = simW[2] + (bias[s] - 1);
      
      groupChoice[s,t] ~ categorical_logit( tau[s] * simW );
      
    }
  }    
}

generated quantities {
  // For group level parameters
  real<lower=0, upper=10> mu_tau;
  real<lower=0, upper=2> mu_SAi;
  real<lower=0, upper=2> mu_bias;
  real<lower=0, upper=2> mu_SAo;
  real<lower=0, upper=10> mu_expP;

  // For log likelihood calculation
  real log_lik[nSubjects];

  // For posterior predictive check
  real y_pred[nSubjects, maxTrials];

  // Set all posterior predictions to 0 (avoids NULL values)
  for (i in 1:nSubjects) {
    for (t in 1:maxTrials) {
      y_pred[i, t] = -1;
    }
  }

  //mu_A   = Phi_approx(mu_pr[1]);
  mu_tau = Phi_approx(mu_pr[1]) * 10;
  mu_SAi   = Phi_approx(mu_pr[2]) * 2;
  mu_bias = Phi_approx(mu_pr[3]) * 2;
  mu_SAo   = Phi_approx(mu_pr[4]) * 2;
  mu_expP   = Phi_approx(mu_pr[5]) * 10;

  { // local section, this saves time and space
    
    for (s in 1:nSubjects) {

      vector[2] simW;

      log_lik[s] = 0;

      for (t in 1:(nTrials[s])) {
        
      simW[1] = pow( (sg[s,t] * (SAo[s] - 1)), expP[s]);
      simW[2] = pow( (sg[s,t] * (SAi[s] - 1)), expP[s]);
      simW[2] = simW[2] + (bias[s] - 1);
        
        // compute log likelihood of current trial
        log_lik[s] += categorical_logit_lpmf( groupChoice[s, t] | tau[s] * simW );

        // generate posterior prediction for current trial
        y_pred[s, t] = categorical_rng(softmax( tau[s] * simW ));
        
      }

  }   
    
  }
}

