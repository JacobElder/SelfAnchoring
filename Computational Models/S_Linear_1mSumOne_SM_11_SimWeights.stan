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
  
  real inGsum[nSubjects, maxTrials]; // self-weighted similarity
  real outGsum[nSubjects, maxTrials]; // self-weighted similarity
  
}


parameters {
  // Declare all parameters as vectors for vectorizing
  // Hyper(group)-parameters
  vector[5] mu_pr;
  vector<lower=0>[5] sigma;

  // Subject-level raw parameters (for Matt trick)
  vector[nSubjects] tau_pr;  // inverse temperature
  vector[nSubjects] m_pr;  // slope for ingroup
  vector[nSubjects] wIn_pr;  // slope for ingroup
  vector[nSubjects] wOut_pr;  // slope for ingroup
  vector[nSubjects] mix_pr;  // slope for ingroup
  
}

transformed parameters {
  // subject-level parameters
  vector<lower=0, upper=10>[nSubjects] tau;
  vector<lower=0, upper=1>[nSubjects] m;
  vector<lower=0, upper=2>[nSubjects] wOut;
  vector<lower=0, upper=2>[nSubjects] wIn;
  vector<lower=0, upper=1>[nSubjects] mix;

  for (i in 1:nSubjects) {
    
    tau[i] = Phi_approx(mu_pr[1] + sigma[1] * tau_pr[i]) * 10;
    m[i] = Phi_approx(mu_pr[2] + sigma[2] * m_pr[i]);
    wIn[i] = Phi_approx(mu_pr[3] + sigma[3] * wIn_pr[i]) * 2;
    wOut[i] = Phi_approx(mu_pr[4] + sigma[4] * wOut_pr[i]) * 2;
    mix[i] = Phi_approx(mu_pr[5] + sigma[5] * mix_pr[i]);
    
  }
  
}

model {
  // Hyperparameters
  mu_pr  ~ normal(0, 1);
  sigma ~ normal(0, 0.2);
  //sigma[3:4] ~ cauchy(0, 0.35);

  // individual parameters
  tau_pr ~ normal(0, 1);
  m_pr ~ normal(0, 1);
  wIn_pr ~ normal(0, 1);
  wOut_pr ~ normal(0, 1);
  mix_pr ~ normal(0, 1);

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
    GPin[1:nTrain[s]] = (m[s])*(prevSelf[s,1:nTrain[s]]/7); // Convert evals to ingroup evals
    GPout[1:nTrain[s]] = (1-m[s])*(prevSelf[s,1:nTrain[s]]/7); // Convert evals to outgroup evals
    
    for (t in 1:nTrials[s]) {
      
      PS[1:nTrain[s]] = prevSim[s,t,1:nTrain[s]]; // similarities from training to current generalization trait
      
      simW[1] = dot_product(GPout[1:nTrain[s]],PS)/10; // multiply and sum outgroup evals with sim
      simW[2] = dot_product(GPin[1:nTrain[s]],PS)/10; // multiply and sum ingroup evals with sim
      
      curW[1] = (wOut[s]-1) * outGsum[s,t]; // multiply weight by summed sim-to-outgroup choices thus far
      curW[2] = (wIn[s]-1) * inGsum[s,t]; // multiply weight by summed sim-to-ingroup choices thus far
      
      W = curW * (mix[s]) + simW * (1-mix[s]); // Mixture
      
      groupChoice[s,t] ~ categorical_logit( tau[s] * W ); // decision rule
      
    }
  }    
}

generated quantities {
  // For group level parameters
  real<lower=0, upper=10> mu_tau;
  real<lower=0, upper=1> mu_m;
  real<lower=0, upper=2> mu_inW;
  real<lower=0, upper=2> mu_outW;
  real<lower=0, upper=1> mu_mix;

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
  mu_m   = Phi_approx(mu_pr[2]);
  mu_inW   = Phi_approx(mu_pr[3]) * 2;
  mu_outW   = Phi_approx(mu_pr[4]) * 2;
  mu_mix   = Phi_approx(mu_pr[5]);

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
    
    GPin[1:nTrain[s]] = (m[s])*(prevSelf[s,1:nTrain[s]]/7);
    GPout[1:nTrain[s]] = (1-m[s])*(prevSelf[s,1:nTrain[s]]/7);
    
    for (t in 1:nTrials[s]) {

      PS[1:nTrain[s]] = prevSim[s,t,1:nTrain[s]]; // similarities from training to current generalization trait
      
      simW[1] = dot_product(GPout[1:nTrain[s]],PS)/10;
      simW[2] = dot_product(GPin[1:nTrain[s]],PS)/10;
      
      curW[1] = (wOut[s]-1) * outGsum[s,t];
      curW[2] = (wIn[s]-1) * inGsum[s,t];
      
      W = curW * (mix[s]) + simW * (1-mix[s]);
        
        // compute log likelihood of current trial
        log_lik[s] += categorical_logit_lpmf( groupChoice[s, t] | tau[s] * W );

        // generate posterior prediction for current trial
        y_pred[s, t] = categorical_rng(softmax( tau[s] * W ));
        
      }

  }   
    
  }
}

