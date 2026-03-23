data {
  int<lower=1> nSubjects;
  int<lower=1> nStudies;
  array[nSubjects] int<lower=1, upper=nStudies> subjStudy;
  int<lower=1> maxTrials;
  int<lower=1> maxTrain;
  array[nSubjects] int<lower=1> nTrain;
  array[nSubjects] int<lower=1> nTrials;
  array[nSubjects, maxTrials] int<lower=0, upper=2> groupChoice;
}

parameters {
  real global_mu_pr;
  real<lower=0> global_sigma;
  array[nStudies] real study_mu_raw;
  array[nStudies] real<lower=0> study_sigma;
  vector[nSubjects] bias_pr;
}

transformed parameters {
  vector<lower=0, upper=1>[nSubjects] bias;
  array[nStudies] real study_mu_pr;

  for (s in 1:nStudies) {
    study_mu_pr[s] = global_mu_pr + global_sigma * study_mu_raw[s];
  }
  for (i in 1:nSubjects) {
    int s = subjStudy[i];
    bias[i] = Phi_approx(study_mu_pr[s] + study_sigma[s] * bias_pr[i]);
  }
}

model {
  global_mu_pr  ~ normal(0, 1);
  global_sigma  ~ normal(0, 0.3);

  for (s in 1:nStudies) {
    study_mu_raw[s] ~ normal(0, 1);
    study_sigma[s]  ~ normal(0, 0.3);
  }

  bias_pr ~ normal(0, 1);

  for (s in 1:nSubjects) {
    for (t in 1:nTrials[s]) {
      if (groupChoice[s, t] == 1) {
        target += bernoulli_lpmf(0 | bias[s]);
      } else if (groupChoice[s, t] == 2) {
        target += bernoulli_lpmf(1 | bias[s]);
      }
    }
  }
}

generated quantities {
  real<lower=0, upper=1> mu_bias = Phi_approx(global_mu_pr);

  vector[nSubjects * maxTrials] log_lik = rep_vector(0.0, nSubjects * maxTrials);
  vector[nSubjects * maxTrials] p_pred  = rep_vector(0.0, nSubjects * maxTrials);

  for (s in 1:nSubjects) {
    for (t in 1:nTrials[s]) {
      if (groupChoice[s, t] == 1) {
        log_lik[(s-1)*maxTrials + t] = bernoulli_lpmf(0 | bias[s]);
        p_pred[(s-1)*maxTrials + t]  = bias[s];
      } else if (groupChoice[s, t] == 2) {
        log_lik[(s-1)*maxTrials + t] = bernoulli_lpmf(1 | bias[s]);
        p_pred[(s-1)*maxTrials + t]  = bias[s];
      }
    }
  }
}
