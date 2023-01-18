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
    curIdx <- df$Idx[trial]
    curSims <- simMat[traindf$Idx[traindf$subID==uIds[i]], df$Idx[trial]]
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

## Bias MODEL ##

#############################

iter=3000
warmup=floor(iter/2)
modelFile <- here("Computational Models/Bias.stan")
cores<-detectCores()
biasfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(biasfit)
traceplot(biasfit)
bias_summary <- summary(biasfit, pars = c("bias", "tau"), probs = c(0.1, 0.9))$summary
print(bias_summary)
get_posterior_mean(biasfit, pars=c('SA','tau'))[,5]
biasparams <- data.frame(bias=get_posterior_mean(biasfit, pars=c('bias'))[,5],
                           Temp=get_posterior_mean(biasfit, pars=c('tau'))[,5],
                           LL=get_posterior_mean(biasfit, pars=c('log_lik'))[,5])
k <- 2
biasparams$BIC <- log(lengthArray) * k - 2 * (biasparams$LL)
biasparams$AIC <- 2 * k - 2 * (biasparams$LL)
bias_LL <- extract_log_lik(biasfit)
bias_LOO <- loo(bias_LL)
bias_WAIC <- waic(bias_LL)

#############################

## PROBABILITY BIAS MODEL ##

#############################

iter=2000
warmup=floor(iter/2)
modelFile <- here("Computational Models/ProbBias.stan")
cores<-detectCores()
PBfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52)
traceplot(PBfit)
PB_summary <- summary(PBfit, pars = c("bias", "tau"), probs = c(0.1, 0.9))$summary
print(PB_summary)
get_posterior_mean(PBfit, pars=c('bias','tau'))[,5]
PBparams <- data.frame(bias=get_posterior_mean(PBfit, pars=c('bias'))[,5],
                       Temp=get_posterior_mean(PBfit, pars=c('tau'))[,5],
                       LL=get_posterior_mean(PBfit, pars=c('log_lik'))[,5])
k <- 2
PBparams$BIC <- log(lengthArray) * k - 2 * (PBparams$LL)
PBparams$AIC <- 2 * k - 2 * (PBparams$LL)
PB_LL <- extract_log_lik(PBfit)
PB_LOO <- loo(PB_LL)
PB_WAIC <- waic(PB_LL)

#############################

## PROBABILITY BIAS MODEL ##

#############################

iter=3000
warmup=floor(iter/2)
modelFile <- here("Computational Models/SL2_01.stan")
cores<-detectCores()
SL2fit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(SL2fit)
SL2_summary <- summary(SL2fit, pars = c("tau", "m_in", "m_out"), probs = c(0.1, 0.9))$summary
print(SL2_summary)
get_posterior_mean(SL2fit, pars=c('tau', 'm_in', 'm_out'))[,5]
SL2params <- data.frame(Temp=get_posterior_mean(SL2fit, pars=c('tau'))[,5],
                          m_in=get_posterior_mean(SL2fit, pars=c('m_in'))[,5],
                          m_out=get_posterior_mean(SL2fit, pars=c('m_out'))[,5],
                          LL=get_posterior_mean(SL2fit, pars=c('log_lik'))[,5])
k <- 2
SL2params$BIC <- log(lengthArray) * k - 2 * (SL2params$LL)
SL2params$AIC <- 2 * k - 2 * (SL2params$LL)
SL2_LL <- extract_log_lik(SL2fit)
SL2_LOO <- loo(SL2_LL)
SL2_WAIC <- waic(SL2_LL)

#############################

## PROBABILITY BIAS MODEL ##

#############################

iter=2000
warmup=floor(iter/2)
modelFile <- here("Computational Models/SL2_Shift_01.stan")
cores<-detectCores()
SL2_Shiftfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(SL2_Shiftfit)
SL2_Shift_summary <- summary(SL2_Shiftfit, pars = c("tau", "m_in", "m_out","shift"), probs = c(0.1, 0.9))$summary
print(SL2_Shift_summary)
get_posterior_mean(SL2_Shiftfit, pars=c('tau', 'm_in', 'm_out','shift'))[,5]
SL2_Shiftparams <- data.frame(Temp=get_posterior_mean(SL2_Shiftfit, pars=c('tau'))[,5],
                        m_in=get_posterior_mean(SL2_Shiftfit, pars=c('m_in'))[,5],
                        m_out=get_posterior_mean(SL2_Shiftfit, pars=c('m_out'))[,5],
                        shift=get_posterior_mean(SL2_Shiftfit, pars=c('shift'))[,5],
                        LL=get_posterior_mean(SL2_Shiftfit, pars=c('log_lik'))[,5])
k <- 4
SL2_Shiftparams$BIC <- log(lengthArray) * k - 2 * (SL2_Shiftparams$LL)
SL2_Shiftparams$AIC <- 2 * k - 2 * (SL2_Shiftparams$LL)
SL2_Shift_LL <- extract_log_lik(SL2_Shiftfit)
SL2_Shift_LOO <- loo(SL2_Shift_LL)
SL2_Shift_WAIC <- waic(SL2_Shift_LL)

#############################

## PROBABILITY BIAS MODEL ##

#############################

iter=2000
warmup=floor(iter/2)
modelFile <- here("Computational Models/SL2_2Shift_01.stan")
cores<-detectCores()
SL2_2Shiftfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(SL2_2Shiftfit)
SL2_2Shift_summary <- summary(SL2_2Shiftfit, pars = c("tau", "m_in", "m_out","shift_in", "shift_out"), probs = c(0.1, 0.9))$summary
print(SL2_2Shift_summary)
get_posterior_mean(SL2_2Shiftfit, pars=c('tau', 'm_in', 'm_out','shift'))[,5]
SL2_2Shiftparams <- data.frame(Temp=get_posterior_mean(SL2_2Shiftfit, pars=c('tau'))[,5],
                              m_in=get_posterior_mean(SL2_2Shiftfit, pars=c('m_in'))[,5],
                              m_out=get_posterior_mean(SL2_2Shiftfit, pars=c('m_out'))[,5],
                              shift_in=get_posterior_mean(SL2_2Shiftfit, pars=c('shift_in'))[,5],
                              shift_out=get_posterior_mean(SL2_2Shiftfit, pars=c('shift_out'))[,5],
                              LL=get_posterior_mean(SL2_2Shiftfit, pars=c('log_lik'))[,5])
k <- 5
SL2_2Shiftparams$BIC <- log(lengthArray) * k - 2 * (SL2_2Shiftparams$LL)
SL2_2Shiftparams$AIC <- 2 * k - 2 * (SL2_2Shiftparams$LL)
SL2_2Shift_LL <- extract_log_lik(SL2_2Shiftfit)
SL2_2Shift_LOO <- loo(SL2_2Shift_LL)
SL2_2Shift_WAIC <- waic(SL2_2Shift_LL)

#############################

## SIMILARITY LOGISTIC WITH SHIFT ##

#############################

iter=3000
warmup=floor(iter/2)
modelFile <- here("Computational Models/SL2_Shift_Softmax.stan")
cores<-detectCores()
SL2_Shiftfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(SL2_Shiftfit)
SL2_Shift_summary <- summary(SL2_Shiftfit, pars = c("tau", "m_in", "m_out", "shift"), probs = c(0.1, 0.9))$summary
print(SL2_Shift_summary)
get_posterior_mean(SL2_Shiftfit, pars=c('tau', 'm_in', 'm_out'))[,5]
SL2_Shiftparams <- data.frame(Temp=get_posterior_mean(SL2_Shiftfit, pars=c('tau'))[,5],
                        m_in=get_posterior_mean(SL2_Shiftfit, pars=c('m_in'))[,5],
                        m_out=get_posterior_mean(SL2_Shiftfit, pars=c('m_out'))[,5],
                        shift=get_posterior_mean(SL2_Shiftfit, pars=c('shift'))[,5],
                        LL=get_posterior_mean(SL2_Shiftfit, pars=c('log_lik'))[,5])
k <- 4
SL2_Shiftparams$BIC <- log(lengthArray) * k - 2 * (SL2_Shiftparams$LL)
SL2_Shiftparams$AIC <- 2 * k - 2 * (SL2_Shiftparams$LL)
SL2_Shift_LL <- extract_log_lik(SL2_Shiftfit)
SL2_Shift_LOO <- loo(SL2_Shift_LL)
SL2_Shift_WAIC <- waic(SL2_Shift_LL)

#############################

## SIMILARITY LOGISTIC WITH SHIFT AND BIAS ##

#############################

iter=3000
warmup=floor(iter/2)
modelFile <- here("Computational Models/SL_PB_Shift_01.stan")
cores<-detectCores()
SL_PB_Shiftfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52, control = list(max_treedepth = 15, adapt_delta = 0.99))
traceplot(SL_PB_Shiftfit)
SL_PB_Shift_summary <- summary(SL_PB_Shiftfit, pars = c("tau", "m_in", "m_out", "shift","bias"), probs = c(0.1, 0.9))$summary
print(SL_PB_Shift_summary)
get_posterior_mean(SL_PB_Shiftfit, pars=c('tau', 'm_in', 'm_out','shift','bias'))[,5]
SL_PB_Shiftparams <- data.frame(Temp=get_posterior_mean(SL_PB_Shiftfit, pars=c('tau'))[,5],
                              m_in=get_posterior_mean(SL_PB_Shiftfit, pars=c('m_in'))[,5],
                              m_out=get_posterior_mean(SL_PB_Shiftfit, pars=c('m_out'))[,5],
                              shift=get_posterior_mean(SL_PB_Shiftfit, pars=c('shift'))[,5],
                              bias=get_posterior_mean(SL_PB_Shiftfit, pars=c('bias'))[,5],
                              LL=get_posterior_mean(SL_PB_Shiftfit, pars=c('log_lik'))[,5])
k <- 5
SL_PB_Shiftparams$BIC <- log(lengthArray) * k - 2 * (SL_PB_Shiftparams$LL)
SL_PB_Shiftparams$AIC <- 2 * k - 2 * (SL_PB_Shiftparams$LL)
SL_PB_Shift_LL <- extract_log_lik(SL_PB_Shiftfit)
SL_PB_Shift_LOO <- loo(SL_PB_Shift_LL)
SL_PB_Shift_WAIC <- waic(SL_PB_Shift_LL)

#############################

## SIMILARITY LOGISTIC WITH SHIFT ##

#############################

iter=3000
warmup=floor(iter/2)
modelFile <- here("Computational Models/SL2_Shift_NoTau_Softmax.stan")
cores<-detectCores()
SL2_Shift_NoTaufit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(SL2_Shift_NoTaufit)
SL2_Shift_NoTau_summary <- summary(SL2_Shift_NoTaufit, pars = c("m_in", "m_out", "shift"), probs = c(0.1, 0.9))$summary
print(SL2_Shift_NoTau_summary)
get_posterior_mean(SL2_Shift_NoTaufit, pars=c( 'm_in', 'm_out', 'shift'))[,5]
SL2_Shift_NoTauparams <- data.frame(
                              m_in=get_posterior_mean(SL2_Shift_NoTaufit, pars=c('m_in'))[,5],
                              m_out=get_posterior_mean(SL2_Shift_NoTaufit, pars=c('m_out'))[,5],
                              shift=get_posterior_mean(SL2_Shift_NoTaufit, pars=c('shift'))[,5],
                              LL=get_posterior_mean(SL2_Shift_NoTaufit, pars=c('log_lik'))[,5])
k <- 4
SL2_Shift_NoTauparams$BIC <- log(lengthArray) * k - 2 * (SL2_Shift_NoTauparams$LL)
SL2_Shift_NoTauparams$AIC <- 2 * k - 2 * (SL2_Shift_NoTauparams$LL)
SL2_Shift_NoTau_LL <- extract_log_lik(SL2_Shift_NoTaufit)
SL2_Shift_NoTau_LOO <- loo(SL2_Shift_NoTau_LL)
SL2_Shift_NoTau_WAIC <- waic(SL2_Shift_NoTau_LL)

#############################

## PROBABILITY BIAS MODEL ##

#############################

iter=3000
warmup=floor(iter/2)
modelFile <- here("Computational Models/SL_PB.stan")
cores<-detectCores()
SL_PBfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(SL_PBfit)
SL_PB_summary <- summary(SL_PBfit, pars = c("bias", "tau", "m_in", "m_out"), probs = c(0.1, 0.9))$summary
print(SL_PB_summary)
get_posterior_mean(SL_PBfit, pars=c('bias','tau', 'm_in', 'm_out'))[,5]
SL_PBparams <- data.frame(bias=get_posterior_mean(SL_PBfit, pars=c('bias'))[,5],
                       Temp=get_posterior_mean(SL_PBfit, pars=c('tau'))[,5],
                       m_in=get_posterior_mean(SL_PBfit, pars=c('m_in'))[,5],
                       m_out=get_posterior_mean(SL_PBfit, pars=c('m_out'))[,5],
                       LL=get_posterior_mean(SL_PBfit, pars=c('log_lik'))[,5])
k <- 2
SL_PBparams$BIC <- log(lengthArray) * k - 2 * (SL_PBparams$LL)
SL_PBparams$AIC <- 2 * k - 2 * (SL_PBparams$LL)
SL_PB_LL <- extract_log_lik(SL_PBfit)
SL_PB_LOO <- loo(SL_PB_LL)
SL_PB_WAIC <- waic(SL_PB_LL)

#############################

## TESTING ##

#############################

iter=3000
warmup=floor(iter/2)
modelFile <- here("Computational Models/test.stan")
cores<-detectCores()
testfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(SL_PBfit)
SL_PB_summary <- summary(SL_PBfit, pars = c("bias", "tau", "m_in", "m_out"), probs = c(0.1, 0.9))$summary
print(SL_PB_summary)
get_posterior_mean(SL_PBfit, pars=c('bias','tau', 'm_in', 'm_out'))[,5]
SL_PBparams <- data.frame(bias=get_posterior_mean(SL_PBfit, pars=c('bias'))[,5],
                          Temp=get_posterior_mean(SL_PBfit, pars=c('tau'))[,5],
                          m_in=get_posterior_mean(SL_PBfit, pars=c('m_in'))[,5],
                          m_out=get_posterior_mean(SL_PBfit, pars=c('m_out'))[,5],
                          LL=get_posterior_mean(SL_PBfit, pars=c('log_lik'))[,5])
k <- 2
SL_PBparams$BIC <- log(lengthArray) * k - 2 * (SL_PBparams$LL)
SL_PBparams$AIC <- 2 * k - 2 * (SL_PBparams$LL)
SL_PB_LL <- extract_log_lik(SL_PBfit)
SL_PB_LOO <- loo(SL_PB_LL)
SL_PB_WAIC <- waic(SL_PB_LL)

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
SLfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52, control = list(max_treedepth = 12, adapt_delta = 0.975))
traceplot(SLfit)
SL_summary <- summary(SLfit, pars = c("m_in", "m_out", "tau"), probs = c(0.1, 0.9))$summary
print(SL_summary)
get_posterior_mean(SLfit, pars=c('m_in','m_out','tau'))[,5]
SLparams <- data.frame(m_in=get_posterior_mean(SLfit, pars=c('m_in'))[,5],
                       m_out=get_posterior_mean(SLfit, pars=c('m_out'))[,5],
                           Temp=get_posterior_mean(SLfit, pars=c('tau'))[,5],
                           LL=get_posterior_mean(SLfit, pars=c('log_lik'))[,5])
k <- 3
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

model <- SL_PBfit
y_pred <- rstan::extract(model, pars='y_pred')$y_pred
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
