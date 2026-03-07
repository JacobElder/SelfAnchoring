data {
  int<lower=1> nSubjects;
  int<lower=1> nStudies;
  array[nSubjects] int<lower=1, upper=nStudies> subjStudy;
  int<lower=1> maxTrials;
  int<lower=1> maxTrain;
  array[nSubjects] int<lower=1> nTrain;
  array[nSubjects] int<lower=1> nTrials;
  array[nSubjects, maxTrials] int<lower=0, upper=2> groupChoice;
  array[nSubjects] matrix[maxTrials, maxTrain] prevSim;
  array[nSubjects] vector[maxTrain] prevSelf;
}

parameters {
  // tau (choice temperature) removed: unidentifiable when lapse rate (w) is high.
  // tau fixed at 1 implicitly.
  vector[5] global_mu_pr;            // [1]=m_in, [2]=m_out, [3]=bias, [4]=lambda, [5]=w
  vector<lower=0.01>[5] global_sigma;
  array[nStudies] vector[5] study_mu_raw;
  vector<lower=0.01>[5] study_sigma;
  vector[nSubjects] m_in_pr;
  vector[nSubjects] m_out_pr;
  vector[nSubjects] bias_pr;
  vector[nSubjects] lambda_pr;
  vector[nSubjects] w_pr;
}

transformed parameters {
  vector<lower=0, upper=10>[nSubjects] m_in;
  vector<lower=0, upper=10>[nSubjects] m_out;
  vector<lower=0, upper=1>[nSubjects] bias;
  vector<lower=0, upper=5>[nSubjects] lambda;
  vector<lower=0, upper=1>[nSubjects] w;

  array[nStudies] vector[5] study_mu_pr;
  for (s in 1:nStudies) {
    study_mu_pr[s] = global_mu_pr + global_sigma .* study_mu_raw[s];
  }

  for (i in 1:nSubjects) {
    int s = subjStudy[i];
    m_in[i]   = Phi_approx(study_mu_pr[s, 1] + study_sigma[1] * m_in_pr[i]) * 10;
    m_out[i]  = Phi_approx(study_mu_pr[s, 2] + study_sigma[2] * m_out_pr[i]) * 10;
    bias[i]   = Phi_approx(study_mu_pr[s, 3] + study_sigma[3] * bias_pr[i]);
    lambda[i] = Phi_approx(study_mu_pr[s, 4] + study_sigma[4] * lambda_pr[i]) * 5;
    w[i]      = Phi_approx(study_mu_pr[s, 5] + study_sigma[5] * w_pr[i]);
  }
}

model {
  // Global Priors
  global_mu_pr ~ normal(0, 1);
  global_sigma ~ cauchy(0, 1);

  // Study-level Priors (Non-centered)
  for (s in 1:nStudies) {
    study_mu_raw[s] ~ normal(0, 1);
  }
  study_sigma ~ cauchy(0, 1);

  // Subject-level priors (Non-centered)
  m_in_pr   ~ normal(0, 1);
  m_out_pr  ~ normal(0, 1);
  bias_pr   ~ normal(0, 1);
  lambda_pr ~ normal(0, 1);
  w_pr      ~ normal(0, 1);

  for (s in 1:nSubjects) {
    vector[nTrain[s]] GPin;
    vector[nTrain[s]] GPout;
    matrix[nTrials[s], nTrain[s]] PS;
    vector[nTrials[s]] simW_in;
    vector[nTrials[s]] simW_out;
    vector[nTrials[s]] logit_p;

    GPin[1:nTrain[s]]  = inv_logit(m_in[s] * (prevSelf[s, 1:nTrain[s]] - 4));
    GPout[1:nTrain[s]] = inv_logit(-m_out[s] * (prevSelf[s, 1:nTrain[s]] - 4));
    PS = pow(prevSim[s, 1:nTrials[s], 1:nTrain[s]], lambda[s]);

    simW_in  = PS[1:nTrials[s], 1:nTrain[s]] * GPin[1:nTrain[s]] + 1e-9;
    simW_out = PS[1:nTrials[s], 1:nTrain[s]] * GPout[1:nTrain[s]] + 1e-9;

    logit_p = log(bias[s] + 1e-9) - log(1.0 - bias[s] + 1e-9) + (log(simW_in) - log(simW_out));

    for (t in 1:nTrials[s]) {
      if (groupChoice[s, t] > 0) {
        target += log_sum_exp(log(w[s]) + bernoulli_lpmf(groupChoice[s, t] - 1 | 0.5),
                             log1m(w[s]) + bernoulli_logit_lpmf(groupChoice[s, t] - 1 | logit_p[t]));
      }
    }
  }
}

generated quantities {
  // ================================================================
  // GROUP-LEVEL PARAMETERS (back-transformed global means)
  //   fit$draws("mu_m_in")   # ingroup projection slope [0, 10]
  //   fit$draws("mu_m_out")  # outgroup contrast slope [0, 10]
  //   fit$draws("mu_bias")   # ingroup response bias [0, 1]
  //   fit$draws("mu_lambda") # generalization sensitivity [0, 5]
  //   fit$draws("mu_w")      # lapse rate [0, 1]
  // ================================================================
  real<lower=0, upper=10> mu_m_in   = Phi_approx(global_mu_pr[1]) * 10;
  real<lower=0, upper=10> mu_m_out  = Phi_approx(global_mu_pr[2]) * 10;
  real<lower=0, upper=1>  mu_bias   = Phi_approx(global_mu_pr[3]);
  real<lower=0, upper=5>  mu_lambda = Phi_approx(global_mu_pr[4]) * 5;
  real<lower=0, upper=1>  mu_w      = Phi_approx(global_mu_pr[5]);

  // ================================================================
  // TRIAL-LEVEL QUANTITIES (index = (s-1)*maxTrials + t)
  //   fit$draws("log_lik")   # log-likelihood for LOO-CV
  //   fit$draws("p_pred")    # P(ingroup) per trial
  //   fit$draws("mcr")       # model-derived Metacontrast Ratio per trial
  //                          # MCR_t = simW_in[t] / simW_out[t]
  //                          # (asymmetric GPin/GPout; lambda-sharpened similarity)
  // SUBJECT-LEVEL: fit$draws("subject_mcr")  # mean MCR per subject
  // ================================================================
  vector[nSubjects * maxTrials] log_lik     = rep_vector(0.0, nSubjects * maxTrials);
  vector[nSubjects * maxTrials] p_pred      = rep_vector(0.0, nSubjects * maxTrials);
  vector[nSubjects * maxTrials] mcr         = rep_vector(0.0, nSubjects * maxTrials);
  vector[nSubjects]             subject_mcr = rep_vector(0.0, nSubjects);

  for (s in 1:nSubjects) {
    vector[nTrain[s]] GPin;
    vector[nTrain[s]] GPout;
    matrix[nTrials[s], nTrain[s]] PS;
    vector[nTrials[s]] simW_in;
    vector[nTrials[s]] simW_out;
    vector[nTrials[s]] logit_p;
    int valid_t = 0;
    real mcr_sum = 0.0;

    GPin[1:nTrain[s]]  = inv_logit(m_in[s] * (prevSelf[s, 1:nTrain[s]] - 4));
    GPout[1:nTrain[s]] = inv_logit(-m_out[s] * (prevSelf[s, 1:nTrain[s]] - 4));
    PS = pow(prevSim[s, 1:nTrials[s], 1:nTrain[s]], lambda[s]);

    simW_in  = PS[1:nTrials[s], 1:nTrain[s]] * GPin[1:nTrain[s]] + 1e-9;
    simW_out = PS[1:nTrials[s], 1:nTrain[s]] * GPout[1:nTrain[s]] + 1e-9;

    logit_p = log(bias[s] + 1e-9) - log(1.0 - bias[s] + 1e-9) + (log(simW_in) - log(simW_out));

    for (t in 1:nTrials[s]) {
      if (groupChoice[s, t] > 0) {
        real mcr_t = simW_in[t] / simW_out[t];
        log_lik[(s-1)*maxTrials + t] = log_sum_exp(log(w[s]) + bernoulli_lpmf(groupChoice[s, t] - 1 | 0.5),
                                                  log1m(w[s]) + bernoulli_logit_lpmf(groupChoice[s, t] - 1 | logit_p[t]));
        p_pred[(s-1)*maxTrials + t]  = w[s] * 0.5 + (1.0 - w[s]) * inv_logit(logit_p[t]);
        mcr[(s-1)*maxTrials + t]     = mcr_t;
        mcr_sum  += mcr_t;
        valid_t  += 1;
      }
    }
    subject_mcr[s] = valid_t > 0 ? mcr_sum / valid_t : 0.0;
  }
}
