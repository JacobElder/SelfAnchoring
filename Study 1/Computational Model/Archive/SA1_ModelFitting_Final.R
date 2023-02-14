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

des=array(0,c(maxSubjs,maxTrials))

fulldf <- fulldf[order(fulldf$subID, fulldf$trialTotalT2),]

for(i in 1:length(uIds) ){
  df <- subset(fulldf, subID==uIds[i])
  t <- nrow(df)
  lengthArray[i] <- t
  groupChoice[i,1:t]<-(df$ingChoiceN+1)
  sg[i,1:t] <- (df$SE)
  des[i,1:t] <- (df$desirability)*100
  
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
    inGchoicesGEN <-which(df$ingChoiceN[-trial]==1)
    outGchoicesGEN <-which(df$ingChoiceN[-trial]==0)
    testSimGEN <- simMat[df$Idx[-trial], df$Idx[trial]]
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

iter=2000
warmup=floor(iter/2)
modelFile <- here("Computational Models/Bias.stan")
cores<-detectCores()
biasfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52)#, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(biasfit)
traceplot(biasfit)
bias_summary <- summary(biasfit, pars = c("bias"), probs = c(0.1, 0.9))$summary
print(bias_summary)
get_posterior_mean(biasfit, pars=c('bias'))[,5]
biasparams <- data.frame(bias=get_posterior_mean(biasfit, pars=c('bias'))[,5],
                           LL=get_posterior_mean(biasfit, pars=c('log_lik'))[,5])
k <- 1
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

## SELF-TO-GROUP LINEAR MAPPING MODEL ##

#############################

iter=2000
warmup=floor(iter/2)
modelFile <- here("Computational Models/S_Linear_1mSumOne.stan")
cores<-detectCores()
S_Linear_1mSumOnefit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52)#, control = list(max_treedepth = 12, adapt_delta = 0.95))
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

## SELF-TO-GROUP LINEAR MAPPING WITH PROBABILITY BIAS MODEL ##

#############################

iter=3000
warmup=floor(iter/2)
modelFile <- here("Computational Models/S_Linear_1mSumOne_BiasV2.stan")
cores<-detectCores()
S_Linear_1mSumOne_Biasfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(S_Linear_1mSumOne_Biasfit)
S_Linear_1mSumOne_Bias_summary <- summary(S_Linear_1mSumOne_Biasfit, pars = c("tau", "m", "bias"), probs = c(0.1, 0.9))$summary
print(S_Linear_1mSumOne_Bias_summary)
get_posterior_mean(S_Linear_1mSumOne_Biasfit, pars=c('tau', 'm', 'bias'))[,5]
S_Linear_1mSumOne_Biasparams <- data.frame(Temp=get_posterior_mean(S_Linear_1mSumOne_Biasfit, pars=c('tau'))[,5],
                                      m=get_posterior_mean(S_Linear_1mSumOne_Biasfit, pars=c('m'))[,5],
                                      bias=get_posterior_mean(S_Linear_1mSumOne_Biasfit, pars=c('bias'))[,5],
                                      LL=get_posterior_mean(S_Linear_1mSumOne_Biasfit, pars=c('log_lik'))[,5])
k <- 3
S_Linear_1mSumOne_Biasparams$BIC <- log(lengthArray) * k - 2 * (S_Linear_1mSumOne_Biasparams$LL)
S_Linear_1mSumOne_Biasparams$AIC <- 2 * k - 2 * (S_Linear_1mSumOne_Biasparams$LL)
S_Linear_1mSumOne_Bias_LL <- extract_log_lik(S_Linear_1mSumOne_Biasfit)
S_Linear_1mSumOne_Bias_LOO <- loo(S_Linear_1mSumOne_Bias_LL)
S_Linear_1mSumOne_Bias_WAIC <- waic(S_Linear_1mSumOne_Bias_LL)

#############################

## SELF-TO-GROUP LOGISTIC MAPPING WITH PROBABILITY BIAS MODEL ##

#############################

iter=2000
warmup=floor(iter/2)
modelFile <- here("Computational Models/S_Logistic_1mOppose_Bias.stan")
cores<-detectCores()
S_Logistic_1mOppose_Biasfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52)#, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(S_Logistic_1mOppose_Biasfit)
S_Logistic_1mOppose_Bias_summary <- summary(S_Logistic_1mOppose_Biasfit, pars = c("tau", "m", "bias"), probs = c(0.1, 0.9))$summary
print(S_Logistic_1mOppose_Bias_summary)
get_posterior_mean(S_Logistic_1mOppose_Biasfit, pars=c('tau', 'm', 'bias'))[,5]
S_Logistic_1mOppose_Biasparams <- data.frame(Temp=get_posterior_mean(S_Logistic_1mOppose_Biasfit, pars=c('tau'))[,5],
                                           m=get_posterior_mean(S_Logistic_1mOppose_Biasfit, pars=c('m'))[,5],
                                           bias=get_posterior_mean(S_Logistic_1mOppose_Biasfit, pars=c('bias'))[,5],
                                           LL=get_posterior_mean(S_Logistic_1mOppose_Biasfit, pars=c('log_lik'))[,5])
k <- 3
S_Logistic_1mOppose_Biasparams$BIC <- log(lengthArray) * k - 2 * (S_Logistic_1mOppose_Biasparams$LL)
S_Logistic_1mOppose_Biasparams$AIC <- 2 * k - 2 * (S_Logistic_1mOppose_Biasparams$LL)
S_Logistic_1mOppose_Bias_LL <- extract_log_lik(S_Logistic_1mOppose_Biasfit)
S_Logistic_1mOppose_Bias_LOO <- loo(S_Logistic_1mOppose_Bias_LL)
S_Logistic_1mOppose_Bias_WAIC <- waic(S_Logistic_1mOppose_Bias_LL)

#############################

## SELF-TO-GROUP LOGISTIC MAPPING WITH PROBABILITY BIAS AND CONSTRAINT MODEL ##

#############################

iter=2000
warmup=floor(iter/2)
modelFile <- here("Computational Models/S_Logistic_1mOppose_Bias_L.stan")
cores<-detectCores()
S_Logistic_1mOppose_Bias_Lfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52)#, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(S_Logistic_1mOppose_Bias_Lfit)
S_Logistic_1mOppose_Bias_L_summary <- summary(S_Logistic_1mOppose_Bias_Lfit, pars = c("tau", "m", "bias", "L"), probs = c(0.1, 0.9))$summary
print(S_Logistic_1mOppose_Bias_L_summary)
get_posterior_mean(S_Logistic_1mOppose_Bias_Lfit, pars=c('tau', 'm', 'bias', 'L'))[,5]
S_Logistic_1mOppose_Bias_Lparams <- data.frame(Temp=get_posterior_mean(S_Logistic_1mOppose_Bias_Lfit, pars=c('tau'))[,5],
                                             m=get_posterior_mean(S_Logistic_1mOppose_Bias_Lfit, pars=c('m'))[,5],
                                             bias=get_posterior_mean(S_Logistic_1mOppose_Bias_Lfit, pars=c('bias'))[,5],
                                             L=get_posterior_mean(S_Logistic_1mOppose_Bias_Lfit, pars=c('L'))[,5],
                                             LL=get_posterior_mean(S_Logistic_1mOppose_Bias_Lfit, pars=c('log_lik'))[,5])
k <- 4
S_Logistic_1mOppose_Bias_Lparams$BIC <- log(lengthArray) * k - 2 * (S_Logistic_1mOppose_Bias_Lparams$LL)
S_Logistic_1mOppose_Bias_Lparams$AIC <- 2 * k - 2 * (S_Logistic_1mOppose_Bias_Lparams$LL)
S_Logistic_1mOppose_Bias_L_LL <- extract_log_lik(S_Logistic_1mOppose_Bias_Lfit)
S_Logistic_1mOppose_Bias_L_LOO <- loo(S_Logistic_1mOppose_Bias_L_LL)
S_Logistic_1mOppose_Bias_L_WAIC <- waic(S_Logistic_1mOppose_Bias_L_LL)

#############################

## SELF-TO-GROUP LOGISTIC MAPPING WITH PROBABILITY BIAS AND CONSTRAINT AND OPPOSING SHIFT MODEL ##

#############################

iter=2000
warmup=floor(iter/2)
modelFile <- here("Computational Models/S_Logistic_1mOppose_Bias_L_Shift2.stan")
cores<-detectCores()
S_Logistic_1mOppose_Bias_L_Shiftfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52)#, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(S_Logistic_1mOppose_Bias_L_Shiftfit)
S_Logistic_1mOppose_Bias_L_Shift_summary <- summary(S_Logistic_1mOppose_Bias_L_Shiftfit, pars = c("tau", "m", "bias", "shift"), probs = c(0.1, 0.9))$summary
print(S_Logistic_1mOppose_Bias_L_Shift_summary)
get_posterior_mean(S_Logistic_1mOppose_Bias_L_Shiftfit, pars=c('tau', 'm', 'bias', 'L', 'shift'))[,5]
S_Logistic_1mOppose_Bias_L_Shiftparams <- data.frame(Temp=get_posterior_mean(S_Logistic_1mOppose_Bias_L_Shiftfit, pars=c('tau'))[,5],
                                               m=get_posterior_mean(S_Logistic_1mOppose_Bias_L_Shiftfit, pars=c('m'))[,5],
                                               bias=get_posterior_mean(S_Logistic_1mOppose_Bias_L_Shiftfit, pars=c('bias'))[,5],
                                               L=get_posterior_mean(S_Logistic_1mOppose_Bias_L_Shiftfit, pars=c('L'))[,5],
                                               shift=get_posterior_mean(S_Logistic_1mOppose_Bias_L_Shiftfit, pars=c('shift'))[,5],
                                               LL=get_posterior_mean(S_Logistic_1mOppose_Bias_L_Shiftfit, pars=c('log_lik'))[,5])
k <- 5
S_Logistic_1mOppose_Bias_L_Shiftparams$BIC <- log(lengthArray) * k - 2 * (S_Logistic_1mOppose_Bias_L_Shiftparams$LL)
S_Logistic_1mOppose_Bias_L_Shiftparams$AIC <- 2 * k - 2 * (S_Logistic_1mOppose_Bias_L_Shiftparams$LL)
S_Logistic_1mOppose_Bias_L_Shift_LL <- extract_log_lik(S_Logistic_1mOppose_Bias_L_Shiftfit)
S_Logistic_1mOppose_Bias_L_Shift_LOO <- loo(S_Logistic_1mOppose_Bias_L_Shift_LL)
S_Logistic_1mOppose_Bias_L_Shift_WAIC <- waic(S_Logistic_1mOppose_Bias_L_Shift_LL)

#############################

## SELF-TO-GROUP LOGISTIC MAPPING WITH PROBABILITY BIAS AND SHIFT MODEL ##

#############################

iter=2000
warmup=floor(iter/2)
modelFile <- here("Computational Models/S_Logistic_1mOppose_Bias_Shift.stan")
cores<-detectCores()
S_Logistic_1mOppose_Bias_Shiftfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52)#, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(S_Logistic_1mOppose_Bias_Shiftfit)
S_Logistic_1mOppose_Bias_Shift_summary <- summary(S_Logistic_1mOppose_Bias_Shiftfit, pars = c("tau", "m", "bias", "shift"), probs = c(0.1, 0.9))$summary
print(S_Logistic_1mOppose_Bias_Shift_summary)
get_posterior_mean(S_Logistic_1mOppose_Bias_Shiftfit, pars=c('tau', 'm', 'bias', 'shift'))[,5]
S_Logistic_1mOppose_Bias_Shiftparams <- data.frame(Temp=get_posterior_mean(S_Logistic_1mOppose_Bias_Shiftfit, pars=c('tau'))[,5],
                                             m=get_posterior_mean(S_Logistic_1mOppose_Bias_Shiftfit, pars=c('m'))[,5],
                                             bias=get_posterior_mean(S_Logistic_1mOppose_Bias_Shiftfit, pars=c('bias'))[,5],
                                             shift=get_posterior_mean(S_Logistic_1mOppose_Bias_Shiftfit, pars=c('shift'))[,5],
                                             LL=get_posterior_mean(S_Logistic_1mOppose_Bias_Shiftfit, pars=c('log_lik'))[,5])
k <- 3
S_Logistic_1mOppose_Bias_Shiftparams$BIC <- log(lengthArray) * k - 2 * (S_Logistic_1mOppose_Bias_Shiftparams$LL)
S_Logistic_1mOppose_Bias_Shiftparams$AIC <- 2 * k - 2 * (S_Logistic_1mOppose_Bias_Shiftparams$LL)
S_Logistic_1mOppose_Bias_Shift_LL <- extract_log_lik(S_Logistic_1mOppose_Bias_Shiftfit)
S_Logistic_1mOppose_Bias_Shift_LOO <- loo(S_Logistic_1mOppose_Bias_Shift_LL)
S_Logistic_1mOppose_Bias_Shift_WAIC <- waic(S_Logistic_1mOppose_Bias_Shift_LL)

############################

# MODEL VALIDATION

############################

print(loo_compare(list("Bias"=PB_LOO,
                       "SLinS"=S_Linear_1mSumOne_LOO,
                       "SLinS.Bias"=S_Linear_1mSumOne_Bias_LOO,
                       "SLogO.Bias"=S_Logistic_1mOppose_Bias_LOO,
                       "SLogO.Bias.L"=S_Logistic_1mOppose_Bias_L_LOO
)),simplify = F
)

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


######

y_correct <- y_pred

for(i in 1:numSubjs){
  y_correct[,i,]<-apply(y_pred[,i,], 1, function(x) x==true_y[i,])
}

y_correct_mean = apply(y_correct, c(2,3), mean)  # average of 4000 MCMC samples

dim(y_correct_mean)  # y_pred_mean --> 58 (subjects) x 148 (trials)

numSubjs = dim(y_correct_mean)[2]  # number of subjects

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

min(model_data$nTrials)

qplot(y_correct_mean[,1:min(model_data$nTrials)],bins=100)
mean(y_correct>.50)

qplot(y_correct[,,1:min(model_data$nTrials)],bins=100)
      