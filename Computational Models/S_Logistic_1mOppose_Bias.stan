data {
  int<lower=1> nSubjects; // number of subjects total
  int<lower = 1> maxTrials; //number of testing/generalization trials
  int<lower = 1> maxTrain; // max number of training trials
  int<lower = 1> nTrain[nSubjects]; // per participant number of training trials
  int<lower = 1> nTrials[nSubjects]; //per participant number of testing/generalization trials
  int<lower = 0, upper = 2> groupChoice[nSubjects,maxTrials]; // which group chosen
  
  vector[maxTrain] prevSim[nSubjects, maxTrials]; // 3-dimensional matrix of subject's training traits similarities to testing traits
  //int<lower = 0, upper = 7> prevSelf[nSubjects,maxTrain]; // matrix of training self-evaluations
  vector[maxTrain] prevSelf[nSubjects]; // matrix of training self-evaluations
  
}


parameters {
  // Declare all parameters as vectors for vectorizing
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
  vector<lower=0, upper=10>[nSubjects] tau; # temperature
  vector<lower=0, upper=10>[nSubjects] m; # growth rate
  vector<lower=0, upper=1>[nSubjects] bias; # bias

  for (i in 1:nSubjects) {
    
    tau[i] = Phi_approx(mu_pr[1] + sigma[1] * tau_pr[i]) * 10; 
    m[i] = Phi_approx(mu_pr[2] + sigma[2] * m_pr[i]) * 10; 
    bias[i] = Phi_approx(mu_pr[3] + sigma[3] * bias_pr[i]); 
    
  }
  
}

model {
  // Hyperparameters
  mu_pr  ~ normal(0, 1);
  // sigma ~ normal(0, 0.2);
  sigma ~ normal(0, .3);

  // individual parameters
  tau_pr ~ normal(0, 1);
  // m_pr ~ normal(0, 1);
  m_pr ~ normal(0, 2); // making more diffuse
  bias_pr ~ normal(0, 1);

  for (s in 1:nSubjects) {

    vector[2] simW;
    vector[2] prob;
    vector[nTrain[s]] curSelf;
    vector[nTrain[s]] GPin;
    vector[nTrain[s]] GPout;
    vector[nTrain[s]] PS;

    GPin[1:nTrain[s]] = rep_vector(1,nTrain[s])./(1 + exp((-m[s])*(prevSelf[s,1:nTrain[s]]-4))); # logistic function for transforming eval ratings to ingroup
    GPout[1:nTrain[s]] = rep_vector(1,nTrain[s])./(1 + exp((m[s])*(prevSelf[s,1:nTrain[s]]-4))); # logistic function for transforming eval ratings to outgroup
    
    for (t in 1:nTrials[s]) {
      
      PS[1:nTrain[s]] = prevSim[s,t,1:nTrain[s]]; // similarities from training to current generalization trait
      
      simW[1] = dot_product(GPout[1:nTrain[s]],PS); # dot-product of ingroup and sim
      simW[2] = dot_product(GPin[1:nTrain[s]],PS); #dot-product of ingroup and sim
      
      prob[1] = ( (1 - bias[s]) * pow( ( simW[1] ) ,tau[s] ) ) / ( ( (1 - bias[s]) *  pow( simW[1]  , tau[s] ) ) + ( (bias[s]) * pow( ( simW[2] ) ,tau[s] ) ) ); // convert to probabilities
      prob[2] = ( (bias[s]) * pow( ( simW[2] ) ,tau[s] ) ) / ( ( (1 - bias[s]) *  pow( simW[1]  , tau[s] ) ) + ( (bias[s]) * pow( ( simW[2] ) ,tau[s] ) ) ); // convert to probabilities
      
      groupChoice[s,t] ~ categorical( prob );
      
    }
  }    
}

generated quantities {
  // For group level parameters
  real<lower=0, upper=10> mu_tau;
  real<lower=0, upper=10> mu_m;
  real<lower=0, upper=1> mu_bias;

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
  mu_m   = Phi_approx(mu_pr[2]) * 10;
  mu_bias   = Phi_approx(mu_pr[3]);

  { // local section, this saves time and space
    
    for (s in 1:nSubjects) {

    vector[2] simW;
    vector[2] prob;
    vector[nTrain[s]] curSelf;
    vector[nTrain[s]] GPin;
    vector[nTrain[s]] GPout;
    vector[nTrain[s]] PS;
    
    log_lik[s] = 0;
    
    GPin[1:nTrain[s]] = rep_vector(1,nTrain[s])./(1 + exp((-m[s])*(prevSelf[s,1:nTrain[s]]-4)));
    GPout[1:nTrain[s]] = rep_vector(1,nTrain[s])./(1 + exp((m[s])*(prevSelf[s,1:nTrain[s]]-4)));
    
    for (t in 1:nTrials[s]) {

      PS[1:nTrain[s]] = prevSim[s,t,1:nTrain[s]]; // similarities from training to current generalization trait
      
      simW[1] = dot_product(GPout[1:nTrain[s]],PS);
      simW[2] = dot_product(GPin[1:nTrain[s]],PS);
      
      prob[1] = ( (1 - bias[s]) * pow( ( simW[1] ) ,tau[s] ) ) / ( ( (1 - bias[s]) *  pow( simW[1]  , tau[s] ) ) + ( (bias[s]) * pow( ( simW[2] ) ,tau[s] ) ) ); // convert to probabilities
      prob[2] = ( (bias[s]) * pow( ( simW[2] ) ,tau[s] ) ) / ( ( (1 - bias[s]) *  pow( simW[1]  , tau[s] ) ) + ( (bias[s]) * pow( ( simW[2] ) ,tau[s] ) ) ); // convert to probabilities
        
        // compute log likelihood of current trial
        log_lik[s] += categorical_lpmf( groupChoice[s, t] | prob );

        // generate posterior prediction for current trial
        y_pred[s, t] = categorical_rng(prob);
        
      }

  }   
    
  }
}

