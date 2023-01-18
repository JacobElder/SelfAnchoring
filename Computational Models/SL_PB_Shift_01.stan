data {
  int<lower=1> nSubjects; // number of subjects total
  int<lower = 1> maxTrials; //number of testing/generalization trials
  int<lower = 1> maxTrain; // max number of training trials
  int<lower = 1> nTrain[nSubjects]; // per participant number of training trials
  int<lower = 1> nTrials[nSubjects]; //per participant number of testing/generalization trials
  int<lower = 0, upper = 2> groupChoice[nSubjects,maxTrials]; // which group chosen
  real sg[nSubjects, maxTrials]; // self-weighted similarity
  
  vector[maxTrain] prevSim[nSubjects, maxTrials]; // 3-dimensional matrix of subject's training traits similarities to testing traits
  //int<lower = 0, upper = 7> prevSelf[nSubjects,maxTrain]; // matrix of training self-evaluations
  vector[maxTrain] prevSelf[nSubjects]; // matrix of training self-evaluations
  
  // real prevSim[nSubjects, maxTrials, maxTrain];
  // int<lower = 0, upper = 7> prevSelf[nSubjects,maxTrain]; // matrix of training self-evaluations
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
  vector[nSubjects] bias_pr;  // inverse temperature
  vector[nSubjects] tau_pr;  // inverse temperature
  vector[nSubjects] m_in_pr;  // slope for ingroup
  vector[nSubjects] m_out_pr;  // slope for outgroup
  vector[nSubjects] shift_pr;  // slope for outgroup
  
}

transformed parameters {
  // subject-level parameters
  vector<lower=0, upper=1>[nSubjects] bias;
  vector<lower=0, upper=10>[nSubjects] tau;
  vector<lower=0, upper=5>[nSubjects] m_in;
  vector<lower=0, upper=5>[nSubjects] m_out;
  vector<lower=0, upper=6>[nSubjects] shift;

  for (i in 1:nSubjects) {
    
    bias[i] = Phi_approx(mu_pr[1] + sigma[1] * bias_pr[i]);
    tau[i] = Phi_approx(mu_pr[2] + sigma[2] * tau_pr[i]) * 10;
    m_in[i] = Phi_approx(mu_pr[3]  + sigma[3]  * m_in_pr[i]) * 5;
    m_out[i] = Phi_approx(mu_pr[4]  + sigma[4]  * m_out_pr[i]) * 5;
    shift[i] = Phi_approx(mu_pr[5] + sigma[5] * shift_pr[i]) * 6;
    
  }
  
}

model {
  // Hyperparameters
  mu_pr  ~ normal(0, 1);
  sigma ~ normal(0, 0.2);

  // individual parameters
  bias_pr ~ normal(0, 1);
  tau_pr ~ normal(0, 1);
  m_in_pr ~ normal(0, 1);
  m_out_pr ~ normal(0, 1);
  shift_pr ~ normal(0, 1);
  
  # More explicit priors
  // bias_pr ~ normal(0, 1);
  // tau_pr ~ normal(0, 1);
  // m_in_pr ~ normal(1, 2.5);
  // m_out_pr ~ normal(1, 2.5);
  // shift_pr ~ normal(0, 1);

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
      
      prob[1] = ( (1 - bias[s]) * pow( ( simW[1] ) ,tau[s] ) ) / ( ( (1 - bias[s]) *  pow( simW[1]  , tau[s] ) ) + ( (bias[s]) * pow( ( simW[2] ) ,tau[s] ) ) ); // convert to probabilities
      prob[2] = ( (bias[s]) * pow( ( simW[2] ) ,tau[s] ) ) / ( ( (1 - bias[s]) *  pow( simW[1]  , tau[s] ) ) + ( (bias[s]) * pow( ( simW[2] ) ,tau[s] ) ) ); // convert to probabilities
      
      // print("Choice Probability")
      // print(prob);
      
      groupChoice[s,t] ~ categorical( prob );
      
    }
  }    
}

generated quantities {
  // For group level parameters
  real<lower=0, upper=1> mu_bias;
  real<lower=0, upper=10> mu_tau;
  real<lower=0, upper=5> mu_m_in;
  real<lower=0, upper=5> mu_m_out;
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
  mu_bias = Phi_approx(mu_pr[1]);
  mu_tau = Phi_approx(mu_pr[2]) * 10;
  mu_m_in   = Phi_approx(mu_pr[3]) * 5;
  mu_m_out = Phi_approx(mu_pr[4]) * 5;
  mu_shift = Phi_approx(mu_pr[5]) * 6;

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
      
      prob[1] = ( (1 - bias[s]) * pow( ( simW[1] ) ,tau[s] ) ) / ( ( (1 - bias[s]) *  pow( simW[1]  , tau[s] ) ) + ( (bias[s]) * pow( ( simW[2] ) ,tau[s] ) ) ); // convert to probabilities
      prob[2] = ( (bias[s]) * pow( ( simW[2] ) ,tau[s] ) ) / ( ( (1 - bias[s]) *  pow( simW[1]  , tau[s] ) ) + ( (bias[s]) * pow( ( simW[2] ) ,tau[s] ) ) ); // convert to probabilities
      
      // print(prob)
      // print(sum(prob))
      // print(groupChoice[s,t])
      // print(categorical_lpmf( groupChoice[s, t] | prob ))
        
        // compute log likelihood of current trial
        log_lik[s] += categorical_lpmf( groupChoice[s, t] | prob );
        
        // print(log_lik[s])

        // generate posterior prediction for current trial
        y_pred[s, t] = categorical_rng(prob);
        
      }

  }   
    
  }
}

