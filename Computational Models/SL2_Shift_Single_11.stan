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
  vector[3] mu_pr;
  vector<lower=0>[3] sigma;

  // Subject-level raw parameters (for Matt trick)
  vector[nSubjects] tau_pr;  // inverse temperature
  vector[nSubjects] m_pr;  // slope for ingroup
  vector[nSubjects] shift_pr;  // slope for shift
  
}

transformed parameters {
  // subject-level parameters
  vector<lower=0, upper=10>[nSubjects] tau;
  vector<lower=0, upper=20>[nSubjects] m;
  vector<lower=0, upper=4>[nSubjects] shift;

  for (i in 1:nSubjects) {
    
    tau[i] = Phi_approx(mu_pr[1] + sigma[1] * tau_pr[i]) * 10;
    m[i] = Phi_approx(mu_pr[2] + sigma[2] * m_pr[i]) * 20;
    shift[i] = Phi_approx(mu_pr[3] + sigma[3] * shift_pr[i]) * 4;
    
  }
  
}

model {
  // Hyperparameters
  mu_pr  ~ normal(0, 1);
  sigma ~ normal(0, 0.2);
  // sigma[1:2] ~ normal(0, 0.2);
  //sigma[3:4] ~ normal(0, 1.0);
  //sigma[3:4] ~ cauchy(0, 0.35);
  //sigma[3] ~ cauchy(0, 1.0);

  // individual parameters
  tau_pr ~ normal(0, 1);
  m_pr ~ normal(0, 1);
  shift_pr ~ normal(0, 1);

  for (s in 1:nSubjects) {

    vector[2] simW;
    vector[2] prob;
    vector[nTrain[s]] curSelf;
    vector[nTrain[s]] GPin;
    vector[nTrain[s]] PS;

    // GPin[1:nTrain[s]] = rep_vector(L[s],nTrain[s])./(1 + exp((-(m[s]-10) )*(prevSelf[s,1:nTrain[s]] -  (shift[s]+1) )));

    for (t in 1:nTrials[s]) {
      
      // print(GPin)
      // print(GPout)
      
      PS[1:nTrain[s]] = prevSim[s,t,1:nTrain[s]]; // similarities from training to current generalization trait
      
      // prob[2] = (1)/(1 + exp((-(m[s]-10) )*( (( dot_product(prevSelf[s,1:nTrain[s]],PS) / sum(prevSelf[s,1:nTrain[s]]) ) )  - (shift[s]) ) )) ;
      // prob[2] = (1)/(1 + exp((-(m[s]-10) )*( (( dot_product(prevSelf[s,1:nTrain[s]],PS)/sum(prevSelf[s,1:nTrain[s]]) ) )  - (shift[s]) ) )) ;
      prob[2] = (1)/(1 + exp((-(m[s]-10) )*( (( ( sg[s,t] ) - .2499395 )/ .06519279 )  - (shift[s]-2) ) )) ;
      // print("Probabilities")
      // print(prob)
      // print(sum(prob))
      prob[1] = 1 - prob[2];

      prob[2] = ( pow( ( prob[2] ) ,tau[s] ) ) / ( ( pow( prob[1]  , tau[s] ) ) + ( pow( ( prob[2] ) ,tau[s] ) ) ); // convert to probabilities
      // print("ProbabilitiesTau")
      // print(prob)
      // print(sum(prob))
      prob[1] = 1 - prob[2];
   
      //groupChoice[s,t] ~ categorical( prob );
      (groupChoice[s,t]-1) ~ bernoulli( prob[2] );
      
    }
  }    
}

generated quantities {
  // For group level parameters
  real<lower=0, upper=10> mu_tau;
  real<lower=0, upper=20> mu_m;
  real<lower=0, upper=4> mu_shift;

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
  mu_m = Phi_approx(mu_pr[2]) * 20;
  mu_shift = Phi_approx(mu_pr[3]) * 4;

  { // local section, this saves time and space
    
    for (s in 1:nSubjects) {

    vector[2] simW;
    vector[2] prob;
    vector[nTrain[s]] curSelf;
    vector[nTrain[s]] GPin;
    vector[nTrain[s]] PS;
    
    log_lik[s] = 0;
    
    for (t in 1:nTrials[s]) {

      PS[1:nTrain[s]] = prevSim[s,t,1:nTrain[s]]; // similarities from training to current generalization trait
      
      //prob[2] = (1)/(1 + exp((-(m[s]-10) )*( ( dot_product(prevSelf[s,1:nTrain[s]],PS) / sum(prevSelf[s,1:nTrain[s]]) )   - (shift[s]) ) )) ;
      // prob[2] = (1)/(1 + exp((-(m[s]-10) )*( (( dot_product(prevSelf[s,1:nTrain[s]],PS)/sum(prevSelf[s,1:nTrain[s]]) ) )  - (shift[s]) ) )) ;
      prob[2] = (1)/(1 + exp((-(m[s]-10) )*( (( ( sg[s,t] ) - .2499395 )/ .06519279 )  - (shift[s]-2) ) )) ;
      prob[1] = 1 - prob[2];

      prob[2] = ( pow( ( prob[2] ) ,tau[s] ) ) / ( ( pow( prob[1]  , tau[s] ) ) + ( pow( ( prob[2] ) ,tau[s] ) ) ); // convert to probabilities
      prob[1] = 1 - prob[2];
        
        // compute log likelihood of current trial
        // log_lik[s] += categorical_lpmf( groupChoice[s, t] | prob );
        log_lik[s] += bernoulli_lpmf( (groupChoice[s, t]-1) | prob[2] );
        

        // generate posterior prediction for current trial
        // y_pred[s, t] = categorical_rng(prob);
        y_pred[s, t] = bernoulli_rng(prob[2]);
        
      }

  }   
    
  }
}

