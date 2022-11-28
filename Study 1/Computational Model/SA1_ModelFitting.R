library("rstan") # observe startup messages
library("tidyverse")
library(doParallel)
library(igraph)
library(loo)
library(here)

fulldf <- read.csv(here("Study 1/Cleaning/output/fullTest.csv"))
fulldf <- fulldf[!is.na(fulldf$ingChoiceN),]

traindf <- read.csv(here("Study 1/Cleaning/output/fullTrain.csv"))
traindf <- traindf[!is.na(traindf$selfResp),]

allPosCents <- read.csv("/Volumes/GoogleDrive/My Drive/Volumes/Research Project/Trait Network_Behaviral/generating network/output/allPosCents.csv")
allNegCents <-read.csv("/Volumes/GoogleDrive/My Drive/Volumes/Research Project/Trait Network_Behaviral/generating network/output/allNegCents.csv")
allCombCents <- rbind(allPosCents, allNegCents)

setwd("~/Google Drive/Volumes/")
posDf <- read.csv("./Research Project/Trait Network_Behaviral/generating network/output/adjacencyMatrix_p.csv")
posMat <- as.matrix(posDf)
posGraph <- graph.adjacency(posMat)
negDf <- read.csv("./Research Project/Trait Network_Behaviral/generating network/output/adjacencyMatrix_n.csv")
negMat <- as.matrix(negDf)
negGraph <- graph.adjacency(negMat)
simMat <- similarity.dice(posGraph)

Idxmat<-cbind(fulldf$IdxLeft, fulldf$IdxRight)
fulldf$IdxChoose <- Idxmat[cbind(seq_along(fulldf$choice), fulldf$choice)]

# uIds <- unique(fulldf$subID)
# Idxmat<-matrix(ncol=2, nrow=nrow(fulldf), fulldf$IdxLeft, fulldf$IdxRight)
# fulldf$Idxchoose <-Idxmat[cbind(seq_along(fulldf$choice), fulldf$choice)]
# uIds<-unique(fulldf$subID)
# fulldf$V_S <- NA
# for(i in uIds){
#   
#   subDf <- subset(fulldf, subID==i)
#   for(n in 1:nrow(subDf)){
#     
#     V_S <- .50
#     ind <- subDf$Idxchoose[n]
#     prevs <- fulldf$Idxchoose[1:(n-1)]
#     
#     if(ind<149 & sum(prevs<149) > 0){
#       indC <- ind
#       prevsP <- prevs[which(prevs<149)]
#       curSim <- similarity.dice(posGraph)[indC,prevsP]
#       
#       prevFeed <- subDf$choiceProp[1:(n-1)]
#       prevFeed2 <- prevFeed[which(prevs<149)]
#     }else if(ind>148 & sum(prevs>148) > 0 ){
#       indC <- ind-148
#       prevsN <- prevs[which(prevs>148)] - 148
#       curSim <- similarity.dice(negGraph)[indC,prevsN]
#       
#       prevFeed <- subDf$choiceProp[1:(n-1)]
#       prevFeed2 <- prevFeed[which(prevs>148)]
#     }
#     
#     V_S = sum(curSim * prevFeed2) / sum(curSim)
#     
#     if(is.na(V_S)){
#       break
#     }
#     
#     fulldf$V_S[fulldf$subID == i & fulldf$trialTotal==subDf$trialTotal[n]] <- V_S
#   }
#   
# }

uIds<-unique(fulldf$subID)
maxSubjs=length(uIds)
maxTrials=max(fulldf$trialTotal)
maxTrain=91

groupChoice=array(0,c(maxSubjs, maxTrials))
sg=array(0,c(maxSubjs,maxTrials))
lengthArray <- matrix(nrow=length(uIds))
lenTrain <- matrix(nrow=length(uIds))
prevSim=array(0,c(maxSubjs,maxTrials,maxTrain))
prevSelf=array(0,c(maxSubjs,maxTrain))

for(i in 1:length(uIds) ){
  df <- subset(fulldf, subID==uIds[i])
  t <- nrow(df)
  lengthArray[i] <- t
  groupChoice[i,1:t]<-(df$ingChoiceN+1)
  sg[i,1:t] <- (df$SE)
  
  lenTrain[i] <- nrow(traindf[traindf$subID==uIds[i],])
  prevSelf[i,1:lenTrain[i]] <- traindf$selfResp[traindf$subID==uIds[i]]
  for(trial in 1:t){
    subDf <- fulldf[fulldf$subID==uIds[i],]
    curIdx <- subDf$Idx[trial]
    curSims <- simMat[traindf$Idx[traindf$subID==uIds[i]], subDf$Idx[trial]]
    prevSim[i,trial,1:lenTrain[i]]=curSims
  }
  
}

model_data <- list( nSubjects=maxSubjs,
                    #nArms = length(unique(df$choiceCue)),
                    maxTrials = maxTrials,
                    nTrials = as.numeric(lengthArray),
                    groupChoice = groupChoice,
                    sg = sg,
                    prevSelf=prevSelf,
                    prevSim=prevSim,
                    maxTrain=91,
                    nTrain=as.numeric(lenTrain))



#############################

## SELF-ANCHOR WITH LINEAR DEPENDENCE MODEL ##

#############################

iter=3000
warmup=floor(iter/2)
modelFile <- here("Computational Models/Anchor.stan")
cores<-detectCores()
anchorfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(biasanchorfit)
traceplot(anchorfit)
anchor_summary <- summary(anchorfit, pars = c("SA", "tau"), probs = c(0.1, 0.9))$summary
print(anchor_summary)
get_posterior_mean(anchorfit, pars=c('SA','tau'))[,5]
anchorparams <- data.frame(SA=get_posterior_mean(basefit, pars=c('SA'))[,5],
                          Temp=get_posterior_mean(basefit, pars=c('tau'))[,5],
                          LL=get_posterior_mean(basefit, pars=c('log_lik'))[,5])
k <- 2
anchorparams$BIC <- log(lengthArray) * k - 2 * (anchorparams$LL)
anchorparams$AIC <- 2 * k - 2 * (anchorparams$LL)
anchor_LL <- extract_log_lik(anchorfit)
anchor_LOO <- loo(anchor_LL)
anchor_WAIC <- waic(anchor_LL)

#############################

## BIAS AND SELF-ANCHOR WITH LINEAR DEPENDENCE MODEL ##

#############################

iter=3000
warmup=floor(iter/2)
modelFile <- here("Computational Models/BiasAnchor.stan")
cores<-detectCores()
biasanchorfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(biasanchorfit)
biasanchor_summary <- summary(biasanchorfit, pars = c("SA", "tau", "bias"), probs = c(0.1, 0.9))$summary
print(biasanchor_summary)
get_posterior_mean(biasanchorfit, pars=c('SA','tau', 'bias'))[,5]
biasanchorparams <- data.frame(SA=get_posterior_mean(biasanchorfit, pars=c('SA'))[,5],
                           Temp=get_posterior_mean(biasanchorfit, pars=c('tau'))[,5],
                           bias=get_posterior_mean(biasanchorfit, pars=c('bias'))[,5],
                           LL=get_posterior_mean(biasanchorfit, pars=c('log_lik'))[,5])
k <- 3
biasanchorparams$BIC <- log(lengthArray) * k - 2 * (biasanchorparams$LL)
biasanchorparams$AIC <- 2 * k - 2 * (biasanchorparams$LL)
biasanchor_LL <- extract_log_lik(biasanchorfit)
biasanchor_LOO <- loo(biasanchor_LL)
biasanchor_WAIC <- waic(biasanchor_LL)

#############################

## BIAS AND SELF-ANCHOR WITH INDDEPENDENCE MODEL ##

#############################

iter=3000
warmup=floor(iter/2)
modelFile <- here("Computational Models/BiasAnchor2.stan")
cores<-detectCores()
biasanchor2fit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(biasanchor2fit)
biasanchor2_summary <- summary(biasanchor2fit, pars = c("SAi","SAo", "tau", "bias"), probs = c(0.1, 0.9))$summary
print(biasanchor2_summary)
get_posterior_mean(biasanchor2fit, pars=c('SAi','SAo','tau', 'bias'))[,5]
biasanchor2params <- data.frame(SAi=get_posterior_mean(biasanchor2fit, pars=c('SAi'))[,5],
                                SAo=get_posterior_mean(biasanchor2fit, pars=c('SAo'))[,5],
                               Temp=get_posterior_mean(biasanchor2fit, pars=c('tau'))[,5],
                               bias=get_posterior_mean(biasanchor2fit, pars=c('bias'))[,5],
                               LL=get_posterior_mean(biasanchor2fit, pars=c('log_lik'))[,5])
k <- 4
biasanchor2params$BIC <- log(lengthArray) * k - 2 * (biasanchor2params$LL)
biasanchor2params$AIC <- 2 * k - 2 * (biasanchor2params$LL)
biasanchor2_LL <- extract_log_lik(biasanchor2fit)
biasanchor2_LOO <- loo(biasanchor2_LL)
biasanchor2_WAIC <- waic(biasanchor2_LL)

#############################

## BIAS AND SELF-ANCHOR WITH INDEPENDENCE MODEL ##

#############################

iter=3000
warmup=floor(iter/2)
modelFile <- here("Computational Models/BiasAnchor2Exp.stan")
cores<-detectCores()
biasanchor2expfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(biasanchor2expfit)
biasanchor2exp_summary <- summary(biasanchor2expfit, pars = c("SAi","SAo", "tau", "bias","expP"), probs = c(0.1, 0.9))$summary
print(biasanchor2exp_summary)
get_posterior_mean(biasanchor2expfit, pars=c('SAi','SAo','tau', 'bias'))[,5]
biasanchor2expparams <- data.frame(SAi=get_posterior_mean(biasanchor2expfit, pars=c('SAi'))[,5],
                                SAo=get_posterior_mean(biasanchor2expfit, pars=c('SAo'))[,5],
                                Temp=get_posterior_mean(biasanchor2expfit, pars=c('tau'))[,5],
                                bias=get_posterior_mean(biasanchor2expfit, pars=c('bias'))[,5],
                                expP=get_posterior_mean(biasanchor2expfit, pars=c('expP'))[,5],
                                LL=get_posterior_mean(biasanchor2expfit, pars=c('log_lik'))[,5])
k <- 4
biasanchor2expparams$BIC <- log(lengthArray) * k - 2 * (biasanchor2expparams$LL)
biasanchor2expparams$AIC <- 2 * k - 2 * (biasanchor2expparams$LL)
biasanchor2exp_LL <- extract_log_lik(biasanchor2expfit)
biasanchor2exp_LOO <- loo(biasanchor2exp_LL)
biasanchor2exp_WAIC <- waic(biasanchor2exp_LL)

#############################

## SIMILARITY LOGISTIC MODEL ##

#############################

iter=3000
warmup=floor(iter/2)
modelFile <- here("Computational Models/SimLogistic.stan")
cores<-detectCores()
SLfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(SLfit)
SL_summary <- summary(SLfit, pars = c("m_in", "m_out", "tau"), probs = c(0.1, 0.9))$summary
print(SL_summary)
get_posterior_mean(SLfit, pars=c('m_in','m_out','tau'))[,5]
SLparams <- data.frame(m_in=get_posterior_mean(SLfit, pars=c('m_in'))[,5],
                       m_out=get_posterior_mean(SLfit, pars=c('m_out'))[,5],
                           Temp=get_posterior_mean(SLfit, pars=c('tau'))[,5],
                           LL=get_posterior_mean(SLfit, pars=c('log_lik'))[,5])
k <- 2
SLparams$BIC <- log(lengthArray) * k - 2 * (SLparams$LL)
SLparams$AIC <- 2 * k - 2 * (SLparams$LL)
SL_LL <- extract_log_lik(SLfit)
SL_LOO <- loo(SL_LL)
SL_WAIC <- waic(SL_LL)


############################

# MODEL VALIDATION

############################

print(loo_compare(list("SA1"=biasanchor_LOO,
                       "SA2"=biasanchor2_LOO
)),simplify = F)

y_pred <- rstan::extract(biasanchor2fit, pars='y_pred')$y_pred
dim(y_pred)

# y_pred --> 6000 (MCMC samples) x 58 (subjects) x 148 (trials)

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
  true_y[i, ] = c((tmpData$ingChoiceN+1),rep(NA,(maxT-length(tmpData$ingChoiceN))))  # only for data with a 'choice' column
}

y_pred_mean[y_pred_mean==-1] <- NA

## Subject #1
plot(true_y[2, ], type="l", xlab="Trial", ylab="Choice (0 or 1)", yaxt="n")
lines(y_pred_mean[2,], col="red", lty=2)
axis(side=2, at = c(0,1) )
legend("bottomleft", legend=c("True", "PPC"), col=c("black", "red"), lty=1:2)

library(bayesplot)
ppc_dens_overlay(true_y[1,1:100], y_pred[3001:6000,1,1:100])