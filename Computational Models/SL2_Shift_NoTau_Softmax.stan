data {
  int<lower=1> nSubjects; // number of subjects total
  int<lower = 1> maxTrials; //number of testing/generalization trials
  int<lower = 1> maxTrain; // max number of training trials
  int<lower = 1> nTrain[nSubjects]; // per participant number of training trials
  int<lower = 1> nTrials[nSubjects]; //per participant number of testing/generalization trials
  int<lower = 0, upper = 2> groupChoice[nSubjects,maxTrials]; // which group chosen
  real sg[nSubjects, maxTrials]; // self-weighted similarity
  
  vector[maxTrain] prevSim[nSubjects, maxTrials]; // 3-dimensional matrix of subject's training traits similarities to testing traits
  vector[maxTrain] prevSelf[nSubjects]; // matrix of training self-evaluations
  
}


parameters {
  // Declare all parameters as vectors for vectorizing
  // Hyper(group)-parameters
  vector[3] mu_pr;
  vector<lower=0>[3] sigma;

  // Subject-level raw parameters (for Matt trick)
  // vector[nSubjects] tau_pr;  // inverse temperature
  vector[nSubjects] m_in_pr;  // slope for ingroup
  vector[nSubjects] m_out_pr;  // slope for outgroup
  vector[nSubjects] shift_pr;
  
}

transformed parameters {
  // subject-level parameters
  // vector<lower=0, upper=10>[nSubjects] tau;
  vector[nSubjects] m_in;
  vector[nSubjects] m_out;
  vector<lower=0, upper=6>[nSubjects] shift;

  for (i in 1:nSubjects) {
    
    // tau[i] = Phi_approx(mu_pr[1] + sigma[1] * tau_pr[i]) * 10;
    
    shift[i] = Phi_approx(mu_pr[3] + sigma[3] * shift_pr[i]) * 6;
    
  }
  
   m_in   = mu_pr[1]  + sigma[1]  * m_in_pr;
   m_out   = mu_pr[2]  + sigma[2]  * m_out_pr;
  
}

model {
  // Hyperparameters
  mu_pr  ~ normal(0, 1);
  // sigma[1] ~ normal(0, 0.2);
  //sigma[3:4] ~ normal(0, 1.0);
  //sigma[3:4] ~ cauchy(0, 0.35);
  sigma[1:2] ~ cauchy(0, 1.0);
  sigma[3] ~ normal(0, 0.2);

  // individual parameters
  // tau_pr ~ normal(0, 1);
  m_in_pr ~ normal(0, 1);
  m_out_pr ~ normal(0, 1);
  shift_pr ~ normal(0, 1);

  for (s in 1:nSubjects) {

    vector[2] simW;
    vector[2] prob;
    vector[nTrain[s]] curSelf;
    vector[nTrain[s]] GPin;
    vector[nTrain[s]] GPout;
    vector[nTrain[s]] PS;

    GPin[1:nTrain[s]] = rep_vector(1,nTrain[s])./(1 + exp((-m_in[s])*(prevSelf[s,1:nTrain[s]] -  (shift[s]+1) )));
    GPout[1:nTrain[s]] = rep_vector(1,nTrain[s])./(1 + exp((-m_out[s])*(prevSelf[s,1:nTrain[s]] - (shift[s]+1) )));
    
    for (t in 1:nTrials[s]) {
      
      // print(GPin)
      // print(GPout)
      
      PS[1:nTrain[s]] = prevSim[s,t,1:nTrain[s]]; // similarities from training to current generalization trait
      
      // simW[1] = dot_product(GPout[1:nTrain[s]],PS); // summation and multiplication of group probabilities with similarities
      // simW[2] = dot_product(GPin[1:nTrain[s]],PS);
      simW[1] = dot_product(GPout[1:nTrain[s]],PS);
      simW[2] = dot_product(GPin[1:nTrain[s]],PS);
      
      // print("Weights")
      // print(simW);
      // 
      // print("Self")
      // print(prevSelf[s,1:nTrain[s]]);
      // print("Similarity")
      // print(PS)
      // print("Inroup Probability")
      // print(GPin[1:nTrain[s]])
      // print("Outroup Probability")
      // print(GPout[1:nTrain[s]])
      // 
      // print("Parameters")
      // print(tau[s]);
      // print(m_in[s]);
      // print(m_out[s]);
      // print(bias[s]);
      
      // prob[1] = ( pow( ( simW[1] ) ,tau[s] ) ) / ( ( pow( simW[1]  , tau[s] ) ) + ( pow( ( simW[2] ) ,tau[s] ) ) ); // convert to probabilities
      // prob[2] = ( pow( ( simW[2] ) ,tau[s] ) ) / ( ( pow( simW[1]  , tau[s] ) ) + ( pow( ( simW[2] ) ,tau[s] ) ) ); // convert to probabilities
      
      // print(simW[1])
      // print(simW[2])
      
      prob[1] = simW[1] / (simW[1] + simW[2]);
      prob[2] = 1 - prob[1];
      
      // print("Choice Probability")
      // print(prob);
      
      // groupChoice[s,t] ~ categorical( prob );
      groupChoice[s,t] ~ categorical_logit( simW );
      
    }
  }    
}

generated quantities {
  // For group level parameters
  // real<lower=0, upper=10> mu_tau;
  real mu_m_in;
  real mu_m_out;
  real<lower=0, upper=6> mu_shift;

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
  // mu_tau = Phi_approx(mu_pr[1]) * 10;
  mu_m_in   = mu_pr[1];
  mu_m_out = mu_pr[2];
  mu_shift = Phi_approx(mu_pr[3]) * 6;

  { // local section, this saves time and space
    
    for (s in 1:nSubjects) {

    vector[2] simW;
    vector[2] prob;
    vector[nTrain[s]] curSelf;
    vector[nTrain[s]] GPin;
    vector[nTrain[s]] GPout;
    vector[nTrain[s]] PS;
    
    log_lik[s] = 0;
    
    GPin[1:nTrain[s]] = rep_vector(1,nTrain[s])./(1 + exp((-m_in[s])*(prevSelf[s,1:nTrain[s]] -  (shift[s]+1) )));
    GPout[1:nTrain[s]] = rep_vector(1,nTrain[s])./(1 + exp((-m_out[s])*(prevSelf[s,1:nTrain[s]] - (shift[s]+1) )));
    
    for (t in 1:nTrials[s]) {

      PS[1:nTrain[s]] = prevSim[s,t,1:nTrain[s]]; // similarities from training to current generalization trait
      // simW[1] = dot_product(GPout[1:nTrain[s]],PS); // summation and multiplication of group probabilities with similarities
      // simW[2] = dot_product(GPin[1:nTrain[s]],PS);
      simW[1] = dot_product(GPout[1:nTrain[s]],PS);
      simW[2] = dot_product(GPin[1:nTrain[s]],PS);
      
      // prob[1] = ( pow( ( simW[1] ) ,tau[s] ) ) / ( ( pow( simW[1]  , tau[s] ) ) + ( pow( ( simW[2] ) ,tau[s] ) ) ); // convert to probabilities
      // prob[2] = ( pow( ( simW[2] ) ,tau[s] ) ) / ( ( pow( simW[1]  , tau[s] ) ) + ( pow( ( simW[2] ) ,tau[s] ) ) ); // convert to probabilities
      
      prob[1] = simW[1] / (simW[1] + simW[2]);
      prob[2] = 1 - prob[1];

      // print(prob)
      // print(sum(prob))
      // print(groupChoice[s,t])
      // print(categorical_lpmf( groupChoice[s, t] | prob ))
        
        // compute log likelihood of current trial
        log_lik[s] += categorical_lpmf( groupChoice[s, t] | prob );
        
        // print(log_lik[s])

        // generate posterior prediction for current trial
        //y_pred[s, t] = categorical_rng(prob);
        y_pred[s, t] = categorical_rng(softmax(simW));
        
      }

  }   
    
  }
}

