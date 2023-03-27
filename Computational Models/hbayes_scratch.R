# hbayes scratch
model <- S_Logistic_1mOppose_Biasfit
y_pred <- rstan::extract(model, pars='y_pred')$y_pred
dim(y_pred)

# y_pred --> 6000 (MCMC samples) x 58 (subjects) x 148 (trials)

y_pred=y_pred-1
y_pred_mean = apply(y_pred, c(2,3), mean)  # average of 4000 MCMC samples

dim(y_pred_mean)  # y_pred_mean --> 58 (subjects) x 148 (trials)

numSubjs = dim(y_pred)[2]  # number of subjects

subjList = uIds  # list of subject IDs
maxT = maxTrials  # maximum number of trials
true_y = array(NA, c(numSubjs, maxT)) # true data (`true_y`)

## true data for each subject
for (i in 1:numSubjs) {
  tmpID = subjList[i]
  tmpData = subset(fulldf, subID == tmpID)
  true_y[i, ] = model_data$groupChoice[i, ]-1
  #true_y[i, ] = c((tmpData$ingChoiceN+1),rep(NA,(maxT-length(tmpData$ingChoiceN))))  # only for data with a 'choice' column
}

y_pred_mean[y_pred_mean==-1] <- NA

for(i in 1:nrow(true_y)){
  print(i)
  print( cor.test(true_y[i,1:model_data$nTrials[i] ],y_pred_mean[i,1:model_data$nTrials[i] ]) )
}

## Subject #1
plot(true_y[24, 1:model_data$nTrials[24] ], type="l", xlab="Trial", ylab="Choice (0 or 1)", yaxt="n")
lines(y_pred_mean[24, 1:model_data$nTrials[24] ], col="red", lty=2)
axis(side=2, at = c(1,2) )
legend("bottomleft", legend=c("True", "PPC"), col=c("black", "red"), lty=1:2)

## Subject #1
plot(true_y[39, 1:model_data$nTrials[39] ], type="l", xlab="Trial", ylab="Choice (0 or 1)", yaxt="n")
lines(y_pred_mean[39, 1:model_data$nTrials[39] ], col="red", lty=2)
axis(side=2, at = c(1,2) )
legend("bottomleft", legend=c("True", "PPC"), col=c("black", "red"), lty=1:2)

plot(true_y[33, 1:model_data$nTrials[33] ], type="l", xlab="Trial", ylab="Choice (0 or 1)", yaxt="n")
lines(y_pred_mean[33, 1:model_data$nTrials[33] ], col="red", lty=2)
axis(side=2, at = c(1,2) )
legend("bottomleft", legend=c("True", "PPC"), col=c("black", "red"), lty=1:2)

library(bayesplot)
ppc_dens_overlay(true_y[39, 1:model_data$nTrials[39]], y_pred[,39, 1:model_data$nTrials[39]])

ppc_ecdf_overlay(true_y[1,1:148], y_pred[,1,1:148], discrete=T)

ppc_dens_overlay_grouped(true_y[1,1:100], y_pred[,1,1:100], group=as.factor(y_pred[1:58,,]))

ppc_ecdf_overlay(true_y[1,1:100], y_pred[,1,1:100], discrete=T)

ppc_stat_2d(true_y[1,1:148], y_pred[,1,1:148], stat = c("mean", "sd"))
