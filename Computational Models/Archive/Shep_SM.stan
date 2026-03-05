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
  
}


parameters {
  // Declare all parameters as vectors for vectorizing
  // Hyper(group)-parameters
  vector[5] mu_pr;
  vector<lower=0>[5] sigma;

  // Subject-level raw parameters (for Matt trick)
  vector[nSubjects] tau_pr;  // inverse temperature
  vector[nSubjects] c_in_pr;  // slope for ingroup
  vector[nSubjects] c_out_pr;  // slope for outgroup
  vector[nSubjects] p_in_pr;  // slope for outgroup
  vector[nSubjects] p_out_pr;  // slope for outgroup
  
}

transformed parameters {
  // subject-level parameters
  vector<lower=0, upper=10>[nSubjects] tau;
  vector<lower=0, upper=2>[nSubjects] c_in;
  vector<lower=0, upper=2>[nSubjects] c_out;
  vector<lower=0, upper=2>[nSubjects] p_in;
  vector<lower=0, upper=2>[nSubjects] p_out;

  for (i in 1:nSubjects) {
    
    tau[i] = Phi_approx(mu_pr[1] + sigma[1] * tau_pr[i]) * 10;
    c_in[i] = Phi_approx(mu_pr[2] + sigma[2] * c_in_pr[i]) * 2;
    c_out[i] = Phi_approx(mu_pr[3] + sigma[3] * c_out_pr[i]) * 2;
    p_in[i] = Phi_approx(mu_pr[4] + sigma[4] * p_in_pr[i]) * 2;
    p_out[i] = Phi_approx(mu_pr[5] + sigma[5] * p_out_pr[i]) * 2;
    
  }
  
}

model {
  // Hyperparameters
  mu_pr  ~ normal(0, 1);
  sigma ~ normal(0, 0.2);
  //sigma[3:4] ~ normal(0, 1.0);
  //sigma[3:4] ~ cauchy(0, 0.35);
  //sigma[2:3] ~ cauchy(0, 1.0);

  // individual parameters
  tau_pr ~ normal(0, 1);
  c_in_pr ~ normal(0, 1);
  c_out_pr ~ normal(0, 1);
  p_in_pr ~ normal(0, 1);
  p_out_pr ~ normal(0, 1);

  for (s in 1:nSubjects) {

    vector[2] simW;
    vector[2] prob;
    vector[nTrain[s]] curSelf;
    vector[nTrain[s]] GPin;
    vector[nTrain[s]] GPout;
    vector[nTrain[s]] PS;

    #GPin[1:nTrain[s]] = rep_vector(1,nTrain[s])./(1 + exp((-m_in[s])*(prevSelf[s,1:nTrain[s]]-4)));
    #GPout[1:nTrain[s]] = rep_vector(1,nTrain[s])./(1 + exp((-m_out[s])*(prevSelf[s,1:nTrain[s]]-4)));
    
    for (t in 1:nTrials[s]) {
      
      // print(GPin)
      // print(GPout)
      
      simW[2] = exp(-c_in[s] * pow( (1/sg[s,t]), p_in[s]) );
      simW[1] = exp(-c_out[s] * pow( (1/sg[s,t]), p_out[s]) );
      // print("INGROUP WEIGHT")
      // print(simW[1])
      // print("OUTGROUP WEIGHT")
      // print(simW[2])
      
      // PS[1:nTrain[s]] = prevSim[s,t,1:nTrain[s]]; // similarities from training to current generalization trait
      
      // simW[1] = dot_product(GPout[1:nTrain[s]],PS); // summation and multiplication of group probabilities with similarities
      // simW[2] = dot_product(GPin[1:nTrain[s]],PS);
      
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
      
      groupChoice[s,t] ~ categorical_logit( tau[s] * simW );
      
    }
  }    
}

generated quantities {
  // For group level parameters
  real<lower=0, upper=10> mu_tau;
  real<lower=0, upper=2> mu_c_in;
  real<lower=0, upper=2> mu_c_out;
  real<lower=0, upper=2> mu_p_in;
  real<lower=0, upper=2> mu_p_out;

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
  mu_c_in = Phi_approx(mu_pr[2]) * 2;
  mu_c_out = Phi_approx(mu_pr[3]) * 2;
  mu_p_in = Phi_approx(mu_pr[4]) * 2;
  mu_p_out = Phi_approx(mu_pr[5]) * 2;

  { // local section, this saves time and space
    
    for (s in 1:nSubjects) {

    vector[2] simW;
    vector[2] prob;
    vector[nTrain[s]] curSelf;
    vector[nTrain[s]] GPin;
    vector[nTrain[s]] GPout;
    vector[nTrain[s]] PS;
    
    log_lik[s] = 0;
    
    // GPin[1:nTrain[s]] = rep_vector(1,nTrain[s])./(1 + exp((-m_in[s])*(prevSelf[s,1:nTrain[s]]-4)));
    // GPout[1:nTrain[s]] = rep_vector(1,nTrain[s])./(1 + exp((-m_out[s])*(prevSelf[s,1:nTrain[s]]-4)));
    
    for (t in 1:nTrials[s]) {

      // PS[1:nTrain[s]] = prevSim[s,t,1:nTrain[s]]; // similarities from training to current generalization trait
      // simW[1] = dot_product(GPout[1:nTrain[s]],PS); // summation and multiplication of group probabilities with similarities
      // simW[2] = dot_product(GPin[1:nTrain[s]],PS);
      // simW[1] = dot_product(GPout[1:nTrain[s]],PS);
      // simW[2] = dot_product(GPin[1:nTrain[s]],PS);
      
      simW[2] = exp(-c_in[s] * pow( (1/sg[s,t]), p_in[s]) );
      simW[1] = exp(-c_out[s] * pow( (1/sg[s,t]), p_out[s]) );
      
      prob[1] = ( pow( ( simW[1] ) ,tau[s] ) ) / ( ( pow( simW[1]  , tau[s] ) ) + ( pow( ( simW[2] ) ,tau[s] ) ) ); // convert to probabilities
      prob[2] = ( pow( ( simW[2] ) ,tau[s] ) ) / ( ( pow( simW[1]  , tau[s] ) ) + ( pow( ( simW[2] ) ,tau[s] ) ) ); // convert to probabilities
      
      // print(prob)
      // print(sum(prob))
      // print(groupChoice[s,t])
      // print(categorical_lpmf( groupChoice[s, t] | prob ))
        
        // compute log likelihood of current trial
        log_lik[s] += categorical_logit_lpmf( groupChoice[s, t] | tau[s] * simW );

        // generate posterior prediction for current trial
        y_pred[s, t] = categorical_rng(softmax( tau[s] * simW ));
        
      }

  }   
    
  }
}

