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
  vector[1] mu_pr;
  vector<lower=0>[1] sigma;

  // Subject-level raw parameters (for Matt trick)
  //vector[nSubjects] A_pr;    // learning rate
  vector[nSubjects] bias_pr;  // inverse temperature
  
}

transformed parameters {
  // subject-level parameters
  vector<lower=0, upper=1>[nSubjects] bias;

  for (i in 1:nSubjects) {
    
    bias[i] = Phi_approx(mu_pr[1] + sigma[1] * bias_pr[i]);
    
  }
  
}

model {
  // Hyperparameters
  mu_pr  ~ normal(0, 1);
  sigma ~ normal(0, 0.2);

  // individual parameters
  bias_pr ~ normal(0, 1);

  for (s in 1:nSubjects) {

    vector[2] prob;
    
    for (t in 1:nTrials[s]) {
      
      prob[1] = (1 - bias[s])/( (1 - bias[s]) + (bias[s]) ); // convert to probabilities
      prob[2] = (bias[s])/( (1 - bias[s]) + (bias[s]) ); // convert to probabilities
      
      groupChoice[s,t] ~ categorical( prob );
      
    }
  }    
}

generated quantities {
  // For group level parameters
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

  mu_bias = Phi_approx(mu_pr[1]);

  { // local section, this saves time and space
    
    for (s in 1:nSubjects) {

    vector[2] prob;
    
    log_lik[s] = 0;
    
    for (t in 1:nTrials[s]) {
      
      prob[1] = (1 - bias[s])/( (1 - bias[s]) + (bias[s]) ); // convert to probabilities
      prob[2] = (bias[s])/( (1 - bias[s]) + (bias[s]) ); // convert to probabilities
        
        // compute log likelihood of current trial
        log_lik[s] += categorical_lpmf( groupChoice[s, t] | prob );

        // generate posterior prediction for current trial
        y_pred[s, t] = categorical_rng(prob);
        
      }

  }   
    
  }
}

