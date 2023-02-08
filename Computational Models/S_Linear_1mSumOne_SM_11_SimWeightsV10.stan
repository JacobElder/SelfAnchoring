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
  
  real inGsumG[nSubjects, maxTrials]; // self-weighted similarity
  real outGsumG[nSubjects, maxTrials]; // self-weighted similarity
  
}


parameters {
  // Declare all parameters as vectors for vectorizing
  // Hyper(group)-parameters
  vector[3] mu_pr;
  vector<lower=0>[3] sigma;

  // Subject-level raw parameters (for Matt trick)
  vector[nSubjects] tau_pr;  // inverse temperature
  vector[nSubjects] m_in_pr;  // slope for ingroup
  vector[nSubjects] m_out_pr;  // slope for ingroup
  
}

transformed parameters {
  // subject-level parameters
  vector<lower=0, upper=10>[nSubjects] tau;
  vector<lower=0, upper=10>[nSubjects] m_in;
  vector<lower=0, upper=10>[nSubjects] m_out;

  for (i in 1:nSubjects) {
    
    tau[i] = Phi_approx(mu_pr[1] + sigma[1] * tau_pr[i]) * 10;
    m_in[i] = Phi_approx(mu_pr[2] + sigma[2] * m_in_pr[i]) * 10;
    m_out[i] = Phi_approx(mu_pr[3] + sigma[3] * m_out_pr[i]) * 10;
    
  }
  
}

model {
  // Hyperparameters
  mu_pr  ~ normal(0, 1);
  sigma ~ normal(0, 0.2);
  //sigma[3:4] ~ cauchy(0, 0.35);

  // individual parameters
  tau_pr ~ normal(0, 1);
  m_in_pr ~ normal(0, 1);
  m_out_pr ~ normal(0, 1);

  for (s in 1:nSubjects) {

    vector[2] simW;
    vector[2] curW;
    vector[2] W;
    vector[2] prob;
    vector[nTrain[s]] curSelf;
    vector[nTrain[s]] GPin;
    vector[nTrain[s]] GPout;
    vector[nTrain[s]] PS;

    // GPin[1:nTrain[s]] = rep_vector(1,nTrain[s])./(1 + exp((-m_in[s])*(prevSelf[s,1:nTrain[s]]-4)));
    // GPout[1:nTrain[s]] = rep_vector(1,nTrain[s])./(1 + exp((-m_out[s])*(prevSelf[s,1:nTrain[s]]-4)));
    // GPin[1:nTrain[s]] = (m[s])*(prevSelf[s,1:nTrain[s]]); // Convert evals to ingroup evals
    // GPout[1:nTrain[s]] = (1-m[s])*(prevSelf[s,1:nTrain[s]]); // Convert evals to outgroup evals
    
    for (t in 1:nTrials[s]) {
      
      GPin[1:nTrain[s]] = m_in[s]*(inGsumG[s,t])*(prevSelf[s,1:nTrain[s]]); // Convert evals to ingroup evals
      GPout[1:nTrain[s]] = m_out[s]*(outGsumG[s,t])*(prevSelf[s,1:nTrain[s]]); // Convert evals to outgroup evals

      PS[1:nTrain[s]] = prevSim[s,t,1:nTrain[s]]; // similarities from training to current generalization trait

      simW[1] = dot_product(GPout[1:nTrain[s]],PS);
      simW[2] = dot_product(GPin[1:nTrain[s]],PS);
      
      prob[1] = ( pow( ( simW[1] ) ,tau[s] ) ) / ( ( pow( simW[1]  , tau[s] ) ) + ( pow( ( simW[2] ) ,tau[s] ) ) ); // convert to probabilities
      prob[2] = ( pow( ( simW[2] ) ,tau[s] ) ) / ( ( pow( simW[1]  , tau[s] ) ) + ( pow( ( simW[2] ) ,tau[s] ) ) ); // convert to probabilities
      
      groupChoice[s,t] ~ categorical( prob );
      
    }
  }    
}

generated quantities {
  // For group level parameters
  real<lower=0, upper=10> mu_tau;
  real<lower=0, upper=10> mu_m_in;
  real<lower=0, upper=10> mu_m_out;

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
  mu_m_in   = Phi_approx(mu_pr[2]) * 10;
  mu_m_out  = Phi_approx(mu_pr[3]) * 10;

  { // local section, this saves time and space
    
    for (s in 1:nSubjects) {

    vector[2] simW;
    vector[2] curW;
    vector[2] W;
    vector[2] prob;
    vector[nTrain[s]] curSelf;
    vector[nTrain[s]] GPin;
    vector[nTrain[s]] GPout;
    vector[nTrain[s]] PS;
    
    log_lik[s] = 0;
    
    // GPin[1:nTrain[s]] = (m[s])*(prevSelf[s,1:nTrain[s]]);
    // GPout[1:nTrain[s]] = (1-m[s])*(prevSelf[s,1:nTrain[s]]);
    
    for (t in 1:nTrials[s]) {

      GPin[1:nTrain[s]] = m_in[s]*(inGsumG[s,t])*(prevSelf[s,1:nTrain[s]]); // Convert evals to ingroup evals
      GPout[1:nTrain[s]] = m_out[s]*(outGsumG[s,t])*(prevSelf[s,1:nTrain[s]]); // Convert evals to outgroup evals

      PS[1:nTrain[s]] = prevSim[s,t,1:nTrain[s]]; // similarities from training to current generalization trait

      simW[1] = dot_product(GPout[1:nTrain[s]],PS);
      simW[2] = dot_product(GPin[1:nTrain[s]],PS);
      
      prob[1] = ( pow( ( simW[1] ) ,tau[s] ) ) / ( ( pow( simW[1]  , tau[s] ) ) + ( pow( ( simW[2] ) ,tau[s] ) ) ); // convert to probabilities
      prob[2] = ( pow( ( simW[2] ) ,tau[s] ) ) / ( ( pow( simW[1]  , tau[s] ) ) + ( pow( ( simW[2] ) ,tau[s] ) ) ); // convert to probabilities
        
        // compute log likelihood of current trial
        log_lik[s] += categorical_lpmf( groupChoice[s, t] | prob );

        // generate posterior prediction for current trial
        y_pred[s, t] = categorical_rng(prob);
        
      }

  }   
    
  }
}

