library("rstan") # observe startup messages
library("tidyverse")
library(doParallel)
library(igraph)
library(loo)
library(here)

#here::i_am("./Study 1/Computational Model/SA1_ModelFitting.R")

fulldf <- read.csv(here("Study 1/Cleaning/output/fullTest.csv"))
fulldf <- fulldf[!is.na(fulldf$ingChoiceN),]

traindf <- read.csv(here("Study 1/Cleaning/output/fullTrain.csv"))
traindf <- traindf[!is.na(traindf$selfResp),]

allPosCents <- read.csv("~/Google Drive/Volumes/Research Project/Trait Network_Behaviral/generating network/output/allPosCents.csv")
allNegCents <-read.csv("~/Google Drive/Volumes/Research Project/Trait Network_Behaviral/generating network/output/allNegCents.csv")
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

inGssv=array(0,c(maxSubjs,maxTrials))
outGssv=array(0,c(maxSubjs,maxTrials))
inGmsv=array(0,c(maxSubjs,maxTrials))
outGmsv=array(0,c(maxSubjs,maxTrials))

inGssvG=array(0,c(maxSubjs,maxTrials))
outGssvG=array(0,c(maxSubjs,maxTrials))
inGmsvG=array(0,c(maxSubjs,maxTrials))
outGmsvG=array(0,c(maxSubjs,maxTrials))

fulldf <- fulldf[order(fulldf$subID, fulldf$trialTotalT2),]

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
    
    # Similarity to Ingroup and Outgroup Classifications thus far
    inGchoices <-which(df$ingChoiceN[df$trialTotalT2 < df$trialTotalT2[trial]]==1)
    outGchoices <-which(df$ingChoiceN[df$trialTotalT2 < df$trialTotalT2[trial]]==0)
    testSim <- simMat[df$Idx[df$trialTotalT2 < df$trialTotalT2[trial]], df$Idx[trial]]
    inGSim <- testSim[inGchoices]
    outGSim <- testSim[outGchoices]
    inGSim <- inGSim[!is.na(inGSim)]
    inGsumsim <- sum(inGSim)
    outGsumsim <- sum(outGSim)
    inGmeansim <- mean(inGSim)
    outGmeansim <- mean(outGSim)
    inGmsv[i,trial] <- inGmeansim
    outGmsv[i,trial] <- outGmeansim
    inGssv[i,trial] <- inGsumsim
    outGssv[i,trial] <- outGsumsim
    inGmsv[is.na(inGmsv)] <- 0
    outGmsv[is.na(outGmsv)] <- 0
    inGssv[is.na(inGssv)] <- 0
    outGssv[is.na(outGssv)] <- 0
    
    # Similarity to Ingroup and Outgroup Classifications disregarding trial
    inGchoicesGEN <-which(df$ingChoiceN==1)
    outGchoicesGEN <-which(df$ingChoiceN==0)
    testSimGEN <- simMat[df$Idx, df$Idx[trial]]
    inGSimGEN <- testSimGEN[inGchoicesGEN]
    outGSimGEN <- testSimGEN[outGchoicesGEN]
    inGSimGEN <- inGSimGEN[!is.na(inGSimGEN)]
    inGsumsimGEN <- sum(inGSimGEN)
    outGsumsimGEN <- sum(outGSimGEN)
    inGmeansimGEN <- mean(inGSimGEN)
    outGmeansimGEN <- mean(outGSimGEN)
    inGmsvG[i,trial] <- inGmeansimGEN
    outGmsvG[i,trial] <- outGmeansimGEN
    inGssvG[i,trial] <- inGsumsimGEN
    outGssvG[i,trial] <- outGsumsimGEN
    
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

model_data2 <- list( nSubjects=maxSubjs,
                    #nArms = length(unique(df$choiceCue)),
                    maxTrials = maxTrials,
                    nTrials = as.numeric(lengthArray),
                    groupChoice = groupChoice,
                    sg = sg,
                    prevSelf=prevSelf,
                    prevSim=prevSim,
                    maxTrain=91,
                    nTrain=as.numeric(lenTrain),
                    inGsum = inGssv,
                    outGsum = outGssv,
                    inGmean = inGmsv,
                    outGmean = outGmsv,
                    inGsumG = inGssvG,
                    outGsumG = outGssvG,
                    inGmeanG = inGmsvG,
                    outGmeanG = outGmsvG)

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
SL2_2Shiftfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52)#, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(SL2_2Shiftfit)
SL2_2Shift_summary <- summary(SL2_2Shiftfit, pars = c("tau", "m_in", "m_out","shift_in", "shift_out"), probs = c(0.1, 0.9))$summary
print(SL2_2Shift_summary)
get_posterior_mean(SL2_2Shiftfit, pars=c('tau', 'm_in', 'm_out','shift_in',"shift_out"))[,5]
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

## PROBABILITY BIAS MODEL ##

#############################

iter=2000
warmup=floor(iter/2)
modelFile <- here("Computational Models/SL2_2Shift_Bias.stan")
cores<-detectCores()
SL2_2Shift_Biasfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52, control = list(max_treedepth = 15, adapt_delta = 0.99))
traceplot(SL2_2Shift_Biasfit)
SL2_2Shift_Bias_summary <- summary(SL2_2Shift_Biasfit, pars = c("tau", "m_in", "m_out","shift_in", "shift_out","bias"), probs = c(0.1, 0.9))$summary
print(SL2_2Shift_Bias_summary)
get_posterior_mean(SL2_2Shift_Biasfit, pars=c('tau', 'm_in', 'm_out','shift_in',"shift_out","bias"))[,5]
SL2_2Shift_Biasparams <- data.frame(Temp=get_posterior_mean(SL2_2Shift_Biasfit, pars=c('tau'))[,5],
                                      m_in=get_posterior_mean(SL2_2Shift_Biasfit, pars=c('m_in'))[,5],
                                      m_out=get_posterior_mean(SL2_2Shift_Biasfit, pars=c('m_out'))[,5],
                                      shift_in=get_posterior_mean(SL2_2Shift_Biasfit, pars=c('shift_in'))[,5],
                                      shift_out=get_posterior_mean(SL2_2Shift_Biasfit, pars=c('shift_out'))[,5],
                                      bias=get_posterior_mean(SL2_2Shift_Biasfit, pars=c('bias'))[,5],
                                      LL=get_posterior_mean(SL2_2Shift_Biasfit, pars=c('log_lik'))[,5])
k <- 6
SL2_2Shift_Biasparams$BIC <- log(lengthArray) * k - 2 * (SL2_2Shift_Biasparams$LL)
SL2_2Shift_Biasparams$AIC <- 2 * k - 2 * (SL2_2Shift_Biasparams$LL)
SL2_2Shift_Bias_LL <- extract_log_lik(SL2_2Shift_Biasfit)
SL2_2Shift_Bias_LOO <- loo(SL2_2Shift_Bias_LL)
SL2_2Shift_Bias_WAIC <- waic(SL2_2Shift_Bias_LL)

#############################

## PROBABILITY BIAS MODEL ##

#############################

iter=2000
warmup=floor(iter/2)
modelFile <- here("Computational Models/SL2_2Shift_Bias_L_01.stan")
cores<-detectCores()
SL2_2Shift_Bias_Lfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52)#, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(SL2_2Shift_Bias_Lfit)
SL2_2Shift_Bias_L_summary <- summary(SL2_2Shift_Bias_Lfit, pars = c("tau", "m_in", "m_out","shift_in", "shift_out","bias","L"), probs = c(0.1, 0.9))$summary
print(SL2_2Shift_Bias_L_summary)
get_posterior_mean(SL2_2Shift_Bias_Lfit, pars=c('tau', 'm_in', 'm_out','shift_in',"shift_out"))[,5]
SL2_2Shift_Bias_Lparams <- data.frame(Temp=get_posterior_mean(SL2_2Shift_Bias_Lfit, pars=c('tau'))[,5],
                               m_in=get_posterior_mean(SL2_2Shift_Bias_Lfit, pars=c('m_in'))[,5],
                               m_out=get_posterior_mean(SL2_2Shift_Bias_Lfit, pars=c('m_out'))[,5],
                               shift_in=get_posterior_mean(SL2_2Shift_Bias_Lfit, pars=c('shift_in'))[,5],
                               shift_out=get_posterior_mean(SL2_2Shift_Bias_Lfit, pars=c('shift_out'))[,5],
                               bias=get_posterior_mean(SL2_2Shift_Bias_Lfit, pars=c('bias'))[,5],
                               L=get_posterior_mean(SL2_2Shift_Bias_Lfit, pars=c('L'))[,5],
                               LL=get_posterior_mean(SL2_2Shift_Bias_Lfit, pars=c('log_lik'))[,5])
k <- 7
SL2_2Shift_Bias_Lparams$BIC <- log(lengthArray) * k - 2 * (SL2_2Shift_Bias_Lparams$LL)
SL2_2Shift_Bias_Lparams$AIC <- 2 * k - 2 * (SL2_2Shift_Bias_Lparams$LL)
SL2_2Shift_Bias_L_LL <- extract_log_lik(SL2_2Shift_Bias_Lfit)
SL2_2Shift_Bias_L_LOO <- loo(SL2_2Shift_Bias_L_LL)
SL2_2Shift_Bias_L_WAIC <- waic(SL2_2Shift_Bias_L_LL)

#############################

## PROBABILITY BIAS MODEL ##

#############################

iter=2000
warmup=floor(iter/2)
modelFile <- here("Computational Models/SL2_Shift_Single_NoTau_11.stan")
cores<-detectCores()
SL2_2Shiftfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52)#, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(SL2_2Shiftfit)
SL2_2Shift_summary <- summary(SL2_2Shiftfit, pars = c("tau", "m","shift"), probs = c(0.1, 0.9))$summary
print(SL2_2Shift_summary)
get_posterior_mean(SL2_2Shiftfit, pars=c('tau', 'm','shift'))[,5]
SL2_2Shiftparams <- data.frame(Temp=get_posterior_mean(SL2_2Shiftfit, pars=c('tau'))[,5],
                                    m=get_posterior_mean(SL2_2Shiftfit, pars=c('m'))[,5],
                                    shift=get_posterior_mean(SL2_2Shiftfit, pars=c('shift'))[,5],
                                    LL=get_posterior_mean(SL2_2Shiftfit, pars=c('logik'))[,5])
k <- 3
SL2_2Shiftparams$BIC <- log(lengthArray) * k - 2 * (SL2_2Shiftparams$LL)
SL2_2Shiftparams$AIC <- 2 * k - 2 * (SL2_2Shiftparams$LL)
SL2_2Shift <- extract_log_lik(SL2_2Shiftfit)
SL2_2ShiftOO <- loo(SL2_2Shift)
SL2_2Shift_WAIC <- waic(SL2_2Shift)

#############################

## PROBABILITY BIAS MODEL ##

#############################

iter=2000
warmup=floor(iter/2)
modelFile <- here("Computational Models/SL2_Shift_Single_NoTau_11_SM.stan")
cores<-detectCores()
SL2_Shift_NoTaufit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52)#, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(SL2_Shift_NoTaufit)
SL2_Shift_NoTau_summary <- summary(SL2_Shift_NoTaufit, pars = c( "m","shift"), probs = c(0.1, 0.9))$summary
print(SL2_Shift_NoTau_summary)
get_posterior_mean(SL2_Shift_NoTaufit, pars=c('m','shift'))[,5]
SL2_Shift_NoTauparams <- data.frame(
                               m=get_posterior_mean(SL2_Shift_NoTaufit, pars=c('m'))[,5],
                               shift=get_posterior_mean(SL2_Shift_NoTaufit, pars=c('shift'))[,5],
                               LL=get_posterior_mean(SL2_Shift_NoTaufit, pars=c('logik'))[,5])
k <- 3
SL2_Shift_NoTauparams$BIC <- log(lengthArray) * k - 2 * (SL2_Shift_NoTauparams$LL)
SL2_Shift_NoTauparams$AIC <- 2 * k - 2 * (SL2_Shift_NoTauparams$LL)
SL2_Shift_NoTau <- extract_log_lik(SL2_Shift_NoTaufit)
SL2_Shift_NoTauOO <- loo(SL2_Shift_NoTau)
SL2_Shift_NoTau_WAIC <- waic(SL2_Shift_NoTau)

#############################

## PROBABILITY BIAS MODEL ##

#############################

iter=2000
warmup=floor(iter/2)
modelFile <- here("Computational Models/SL2_NoShift_NoTau_Single_11.stan")
cores<-detectCores()
SL2_NoShift_NoTaufit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52)#, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(SL2_NoShift_NoTaufit)
SL2_NoShift_NoTau_summary <- summary(SL2_NoShift_NoTaufit, pars = c( "m"), probs = c(0.1, 0.9))$summary
print(SL2_NoShift_NoTau_summary)
get_posterior_mean(SL2_NoShift_NoTaufit, pars=c( 'm'))[,5]
SL2_NoShift_NoTauparams <- data.frame(m=get_posterior_mean(SL2_NoShift_NoTaufit, pars=c('m'))[,5],
                               LL=get_posterior_mean(SL2_NoShift_NoTaufit, pars=c('log_lik'))[,5])
k <- 3
SL2_NoShift_NoTauparams$BIC <- log(lengthArray) * k - 2 * (SL2_NoShift_NoTauparams$LL)
SL2_NoShift_NoTauparams$AIC <- 2 * k - 2 * (SL2_NoShift_NoTauparams$LL)
SL2_NoShift_NoTau <- extract_log_lik(SL2_NoShift_NoTaufit)
SL2_NoShift_NoTauOO <- loo(SL2_NoShift_NoTau)
SL2_NoShift_NoTau_WAIC <- waic(SL2_NoShift_NoTau)

#############################

## PROBABILITY BIAS MODEL ##

#############################

iter=2000
warmup=floor(iter/2)
modelFile <- here("Computational Models/SL2_2Shift_Bias_Single_11.stan")
cores<-detectCores()
SL2_2Shift_Biasfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52)#, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(SL2_2Shift_Biasfit)
SL2_2Shift_Bias_summary <- summary(SL2_2Shift_Biasfit, pars = c("tau", "m","shift","bias"), probs = c(0.1, 0.9))$summary
print(SL2_2Shift_Bias_summary)
get_posterior_mean(SL2_2Shift_Biasfit, pars=c('tau', 'm','shift', 'bias'))[,5]
SL2_2Shift_Biasparams <- data.frame(Temp=get_posterior_mean(SL2_2Shift_Biasfit, pars=c('tau'))[,5],
                                      m=get_posterior_mean(SL2_2Shift_Biasfit, pars=c('m'))[,5],
                                      shift=get_posterior_mean(SL2_2Shift_Biasfit, pars=c('shift'))[,5],
                                      bias=get_posterior_mean(SL2_2Shift_Biasfit, pars=c('bias'))[,5],
                                      LL=get_posterior_mean(SL2_2Shift_Biasfit, pars=c('logik'))[,5])
k <- 4
SL2_2Shift_Biasparams$BIC <- log(lengthArray) * k - 2 * (SL2_2Shift_Biasparams$LL)
SL2_2Shift_Biasparams$AIC <- 2 * k - 2 * (SL2_2Shift_Biasparams$LL)
SL2_2Shift_BiasL <- extract_log_lik(SL2_2Shift_Biasfit)
SL2_2Shift_BiasOO <- loo(SL2_2Shift_BiasL)
SL2_2Shift_Bias_WAIC <- waic(SL2_2Shift_BiasL)

#############################

## PROBABILITY BIAS MODEL ##

#############################

iter=2000
warmup=floor(iter/2)
modelFile <- here("Computational Models/SL2_2Shift_Bias_L_Single_11.stan")
cores<-detectCores()
SL2_2Shift_Bias_Lfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52)#, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(SL2_2Shift_Bias_Lfit)
SL2_2Shift_Bias_L_summary <- summary(SL2_2Shift_Bias_Lfit, pars = c("tau", "m","shift","bias","L"), probs = c(0.1, 0.9))$summary
print(SL2_2Shift_Bias_L_summary)
get_posterior_mean(SL2_2Shift_Bias_Lfit, pars=c('tau', 'm','shift', 'bias', 'L'))[,5]
SL2_2Shift_Bias_Lparams <- data.frame(Temp=get_posterior_mean(SL2_2Shift_Bias_Lfit, pars=c('tau'))[,5],
                                      m=get_posterior_mean(SL2_2Shift_Bias_Lfit, pars=c('m'))[,5],
                                      shift=get_posterior_mean(SL2_2Shift_Bias_Lfit, pars=c('shift'))[,5],
                                      bias=get_posterior_mean(SL2_2Shift_Bias_Lfit, pars=c('bias'))[,5],
                                      L=get_posterior_mean(SL2_2Shift_Bias_Lfit, pars=c('L'))[,5],
                                      LL=get_posterior_mean(SL2_2Shift_Bias_Lfit, pars=c('log_lik'))[,5])
k <- 7
SL2_2Shift_Bias_Lparams$BIC <- log(lengthArray) * k - 2 * (SL2_2Shift_Bias_Lparams$LL)
SL2_2Shift_Bias_Lparams$AIC <- 2 * k - 2 * (SL2_2Shift_Bias_Lparams$LL)
SL2_2Shift_Bias_L_LL <- extract_log_lik(SL2_2Shift_Bias_Lfit)
SL2_2Shift_Bias_L_LOO <- loo(SL2_2Shift_Bias_L_LL)
SL2_2Shift_Bias_L_WAIC <- waic(SL2_2Shift_Bias_L_LL)

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

iter=2000
warmup=floor(iter/2)
modelFile <- here("Computational Models/SL_PB_Shift.stan") # Growth/m with no boundaries kind of works, a bit volatile and some unreliability tho
cores<-detectCores()
SL_PB_Shiftfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52, control = list(max_treedepth = 12, adapt_delta = 0.95))
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

## SHEPARD MODEL ##

#############################

iter=2000
warmup=floor(iter/2)
modelFile <- here("Computational Models/S_Linear_SM.stan")
cores<-detectCores()
S_Linearfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(S_Linearfit)
S_Linear_summary <- summary(S_Linearfit, pars = c("tau", "m_in", "m_out"), probs = c(0.1, 0.9))$summary
print(S_Linear_summary)
get_posterior_mean(S_Linearfit, pars=c('tau', 'm_in', 'm_out'))[,5]
S_Linearparams <- data.frame(Temp=get_posterior_mean(S_Linearfit, pars=c('tau'))[,5],
                        m_in=get_posterior_mean(S_Linearfit, pars=c('m_in'))[,5],
                        m_out=get_posterior_mean(S_Linearfit, pars=c('m_out'))[,5],
                        LL=get_posterior_mean(S_Linearfit, pars=c('log_lik'))[,5])
k <- 2
S_Linearparams$BIC <- log(lengthArray) * k - 2 * (S_Linearparams$LL)
S_Linearparams$AIC <- 2 * k - 2 * (S_Linearparams$LL)
S_Linear_LL <- extract_log_lik(S_Linearfit)
S_Linear_LOO <- loo(S_Linear_LL)
S_Linear_WAIC <- waic(S_Linear_LL)

#############################

## SHEPARD MODEL ##

#############################

iter=2000
warmup=floor(iter/2)
modelFile <- here("Computational Models/S_Linear_1mOppose_SM_11.stan")
cores<-detectCores()
S_Linear_1mOpposefit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(S_Linear_1mOpposefit)
S_Linear_1mOppose_summary <- summary(S_Linear_1mOpposefit, pars = c("tau", "m"), probs = c(0.1, 0.9))$summary
print(S_Linear_1mOppose_summary)
get_posterior_mean(S_Linear_1mOpposefit, pars=c('tau', 'm'))[,5]
S_Linear_1mOpposeparams <- data.frame(Temp=get_posterior_mean(S_Linear_1mOpposefit, pars=c('tau'))[,5],
                                      m=get_posterior_mean(S_Linear_1mOpposefit, pars=c('m'))[,5],
                                      LL=get_posterior_mean(S_Linear_1mOpposefit, pars=c('log_lik'))[,5])
k <- 2
S_Linear_1mOpposeparams$BIC <- log(lengthArray) * k - 2 * (S_Linear_1mOpposeparams$LL)
S_Linear_1mOpposeparams$AIC <- 2 * k - 2 * (S_Linear_1mOpposeparams$LL)
S_Linear_1mOppose_LL <- extract_log_lik(S_Linear_1mOpposefit)
S_Linear_1mOppose_LOO <- loo(S_Linear_1mOppose_LL)
S_Linear_1mOppose_WAIC <- waic(S_Linear_1mOppose_LL)

#############################

## SHEPARD MODEL ##

#############################

iter=2000
warmup=floor(iter/2)
modelFile <- here("Computational Models/S_Linear_1mOpposeV2_SM_11.stan")
cores<-detectCores()
S_Linear_1mOpposeV2fit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(S_Linear_1mOpposeV2fit)
S_Linear_1mOpposeV2_summary <- summary(S_Linear_1mOpposeV2fit, pars = c("tau", "m"), probs = c(0.1, 0.9))$summary
print(S_Linear_1mOpposeV2_summary)
get_posterior_mean(S_Linear_1mOpposeV2fit, pars=c('tau', 'm'))[,5]
S_Linear_1mOpposeV2params <- data.frame(Temp=get_posterior_mean(S_Linear_1mOpposeV2fit, pars=c('tau'))[,5],
                             m=get_posterior_mean(S_Linear_1mOpposeV2fit, pars=c('m'))[,5],
                             LL=get_posterior_mean(S_Linear_1mOpposeV2fit, pars=c('log_lik'))[,5])
k <- 2
S_Linear_1mOpposeV2params$BIC <- log(lengthArray) * k - 2 * (S_Linear_1mOpposeV2params$LL)
S_Linear_1mOpposeV2params$AIC <- 2 * k - 2 * (S_Linear_1mOpposeV2params$LL)
S_Linear_1mOpposeV2_LL <- extract_log_lik(S_Linear_1mOpposeV2fit)
S_Linear_1mOpposeV2_LOO <- loo(S_Linear_1mOpposeV2_LL)
S_Linear_1mOpposeV2_WAIC <- waic(S_Linear_1mOpposeV2_LL)

#############################

## SHEPARD MODEL ##

#############################

iter=2000
warmup=floor(iter/2)
modelFile <- here("Computational Models/S_Linear_1mOpposeMeta_SM_11.stan")
cores<-detectCores()
S_Linear_1mOpposeMetafit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(S_Linear_1mOpposeMetafit)
S_Linear_1mOpposeMeta_summary <- summary(S_Linear_1mOpposeMetafit, pars = c("tau", "m","meta"), probs = c(0.1, 0.9))$summary
print(S_Linear_1mOpposeMeta_summary)
get_posterior_mean(S_Linear_1mOpposeMetafit, pars=c('tau', 'm','meta'))[,5]
S_Linear_1mOpposeMetaparams <- data.frame(Temp=get_posterior_mean(S_Linear_1mOpposeMetafit, pars=c('tau'))[,5],
                                        m=get_posterior_mean(S_Linear_1mOpposeMetafit, pars=c('m'))[,5],
                                        meta=get_posterior_mean(S_Linear_1mOpposeMetafit, pars=c('meta'))[,5],
                                        LL=get_posterior_mean(S_Linear_1mOpposeMetafit, pars=c('log_lik'))[,5])
k <- 3
S_Linear_1mOpposeMetaparams$BIC <- log(lengthArray) * k - 2 * (S_Linear_1mOpposeMetaparams$LL)
S_Linear_1mOpposeMetaparams$AIC <- 2 * k - 2 * (S_Linear_1mOpposeMetaparams$LL)
S_Linear_1mOpposeMeta_LL <- extract_log_lik(S_Linear_1mOpposeMetafit)
S_Linear_1mOpposeMeta_LOO <- loo(S_Linear_1mOpposeMeta_LL)
S_Linear_1mOpposeMeta_WAIC <- waic(S_Linear_1mOpposeMeta_LL)

#############################

## SHEPARD MODEL ##

#############################

iter=2000
warmup=floor(iter/2)
modelFile <- here("Computational Models/S_Linear_1mSumOne_SM_11.stan")
cores<-detectCores()
S_Linear_1mSumOnefit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(S_Linear_1mSumOnefit)
S_Linear_1mSumOne_summary <- summary(S_Linear_1mSumOnefit, pars = c("tau", "m"), probs = c(0.1, 0.9))$summary
print(S_Linear_1mSumOne_summary)
get_posterior_mean(S_Linear_1mSumOnefit, pars=c('tau', 'm'))[,5]
S_Linear_1mSumOneparams <- data.frame(Temp=get_posterior_mean(S_Linear_1mSumOnefit, pars=c('tau'))[,5],
                                      m=get_posterior_mean(S_Linear_1mSumOnefit, pars=c('m'))[,5],
                                      LL=get_posterior_mean(S_Linear_1mSumOnefit, pars=c('log_lik'))[,5])
k <- 2
S_Linear_1mSumOneparams$BIC <- log(lengthArray) * k - 2 * (S_Linear_1mSumOneparams$LL)
S_Linear_1mSumOneparams$AIC <- 2 * k - 2 * (S_Linear_1mSumOneparams$LL)
S_Linear_1mSumOne_LL <- extract_log_lik(S_Linear_1mSumOnefit)
S_Linear_1mSumOne_LOO <- loo(S_Linear_1mSumOne_LL)
S_Linear_1mSumOne_WAIC <- waic(S_Linear_1mSumOne_LL)

#############################

## SHEPARD MODEL ##

#############################

iter=2000
warmup=floor(iter/2)
modelFile <- here("Computational Models/S_Linear_1mSumOne_SM_11_SimWeights.stan")
cores<-detectCores()
S_Linear_1mSumOne_SimWeightsfit <- stan(modelFile, data = model_data2, iter = iter, warmup = warmup, cores = cores-1, seed = 52)#, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(S_Linear_1mSumOne_SimWeightsfit)
S_Linear_1mSumOne_SimWeights_summary <- summary(S_Linear_1mSumOne_SimWeightsfit, pars = c("tau", "m", "wOut", "wIn", "mix"), probs = c(0.1, 0.9))$summary
print(S_Linear_1mSumOne_SimWeights_summary)
get_posterior_mean(S_Linear_1mSumOne_SimWeightsfit, pars=c('tau', 'm', 'wOut', 'wIn', 'mix'))[,5]
S_Linear_1mSumOne_SimWeightsparams <- data.frame(Temp=get_posterior_mean(S_Linear_1mSumOne_SimWeightsfit, pars=c('tau'))[,5],
                                      m=get_posterior_mean(S_Linear_1mSumOne_SimWeightsfit, pars=c('m'))[,5],
                                      wOut=get_posterior_mean(S_Linear_1mSumOne_SimWeightsfit, pars=c('wOut'))[,5],
                                      wIn=get_posterior_mean(S_Linear_1mSumOne_SimWeightsfit, pars=c('wIn'))[,5],
                                      mix=get_posterior_mean(S_Linear_1mSumOne_SimWeightsfit, pars=c('mix'))[,5],
                                      LL=get_posterior_mean(S_Linear_1mSumOne_SimWeightsfit, pars=c('log_lik'))[,5])
k <- 2
S_Linear_1mSumOne_SimWeightsparams$BIC <- log(lengthArray) * k - 2 * (S_Linear_1mSumOne_SimWeightsparams$LL)
S_Linear_1mSumOne_SimWeightsparams$AIC <- 2 * k - 2 * (S_Linear_1mSumOne_SimWeightsparams$LL)
S_Linear_1mSumOne_SimWeights_LL <- extract_log_lik(S_Linear_1mSumOne_SimWeightsfit)
S_Linear_1mSumOne_SimWeights_LOO <- loo(S_Linear_1mSumOne_SimWeights_LL)
S_Linear_1mSumOne_SimWeights_WAIC <- waic(S_Linear_1mSumOne_SimWeights_LL)

#############################

## SHEPARD MODEL ##

#############################

iter=2000
warmup=floor(iter/2)
modelFile <- here("Computational Models/S_Linear_2b_bias_SM_11.stan")
cores<-detectCores()
S_Linear_2b_biasfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(S_Linear_2b_biasfit)
S_Linear_2b_bias_summary <- summary(S_Linear_2b_biasfit, pars = c("tau", "m_in", "m_out","b_in","b_out","bias"), probs = c(0.1, 0.9))$summary
print(S_Linear_2b_bias_summary)
get_posterior_mean(S_Linear_2b_biasfit, pars=c('tau', 'm_in', 'm_out',"b_in","b_out","bias"))[,5]
S_Linear_2b_biasparams <- data.frame(Temp=get_posterior_mean(S_Linear_2b_biasfit, pars=c('tau'))[,5],
                             m_in=get_posterior_mean(S_Linear_2b_biasfit, pars=c('m_in'))[,5],
                             m_out=get_posterior_mean(S_Linear_2b_biasfit, pars=c('m_out'))[,5],
                             b_in=get_posterior_mean(S_Linear_2b_biasfit, pars=c('b_in'))[,5],
                             b_out=get_posterior_mean(S_Linear_2b_biasfit, pars=c('b_out'))[,5],
                             bias=get_posterior_mean(S_Linear_2b_biasfit, pars=c('bias'))[,5],
                             LL=get_posterior_mean(S_Linear_2b_biasfit, pars=c('log_lik'))[,5])
k <- 6
S_Linear_2b_biasparams$BIC <- log(lengthArray) * k - 2 * (S_Linear_2b_biasparams$LL)
S_Linear_2b_biasparams$AIC <- 2 * k - 2 * (S_Linear_2b_biasparams$LL)
S_Linear_2b_bias_LL <- extract_log_lik(S_Linear_2b_biasfit)
S_Linear_2b_bias_LOO <- loo(S_Linear_2b_bias_LL)
S_Linear_2b_bias_WAIC <- waic(S_Linear_2b_bias_LL)

#############################

## SHEPARD MODEL ##

#############################

iter=2000
warmup=floor(iter/2)
modelFile <- here("Computational Models/Shep_SM.stan")
cores<-detectCores()
shepfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(shepfit)
traceplot(shepfit)
shep_summary <- summary(shepfit, pars = c("c_in", "c_out", "p_in", "p_out", "tau"), probs = c(0.1, 0.9))$summary
print(shep_summary)
get_posterior_mean(shepfit, pars=c("c_in", "c_out", "p_in", "p_out", "tau"))[,5]
shepparams <- data.frame(c_in=get_posterior_mean(shepfit, pars=c('c_in'))[,5],
                         c_out=get_posterior_mean(shepfit, pars=c('c_out'))[,5],
                         p_in=get_posterior_mean(shepfit, pars=c('p_in'))[,5],
                         p_out=get_posterior_mean(shepfit, pars=c('p_out'))[,5],
                           Temp=get_posterior_mean(shepfit, pars=c('tau'))[,5],
                           LL=get_posterior_mean(shepfit, pars=c('log_lik'))[,5])
k <- 5
shepparams$BIC <- log(lengthArray) * k - 2 * (shepparams$LL)
shepparams$AIC <- 2 * k - 2 * (shepparams$LL)
shep_LL <- extract_log_lik(shepfit)
shep_LOO <- loo(shep_LL)
shep_WAIC <- waic(shep_LL)

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

print(loo_compare(list("Opp"=S_Linear_1mOppose_LOO,
                       "Sum"=S_Linear_1mSumOne_LOO,
                       "Sim"=S_Linear_1mSumOne_SimWeights_LOO
)),simplify = F
)

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
