library("rstan") # observe startup messages
library("tidyverse")
library(doParallel)
library(igraph)
library(loo)
library(here)

here::i_am("./Study 1/Computational Model/SA1_Modeling_Final.R")

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

## PROBABILITY BIAS MODEL ##

#############################

iter=2000
warmup=floor(iter/2)
modelFile <- here("Computational Models/ProbBias.stan")
cores<-detectCores()
PBfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52)
traceplot(PBfit)
PB_summary <- summary(PBfit, pars = c("bias"), probs = c(0.1, 0.9))$summary
print(PB_summary)
get_posterior_mean(PBfit, pars=c('bias'))[,5]
PBparams <- data.frame(bias=get_posterior_mean(PBfit, pars=c('bias'))[,5],
                       LL=get_posterior_mean(PBfit, pars=c('log_lik'))[,5])
k <- 1
PBparams$BIC <- log(lengthArray) * k - 2 * (PBparams$LL)
PBparams$AIC <- 2 * k - 2 * (PBparams$LL)
PB_LL <- extract_log_lik(PBfit)
PB_LOO <- loo(PB_LL)
PB_WAIC <- waic(PB_LL)

#############################

## SELF-TO-GROUP LOGISTIC MAPPING WITH PROBABILITY BIAS MODEL ##

#############################

iter=3000
warmup=floor(iter/2)
modelFile <- here("Computational Models/S_Logistic_1mOppose_Bias.stan")
cores<-detectCores()
S_Logistic_1mOppose_Biasfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52, control = list(max_treedepth = 12, adapt_delta = 0.95))
traceplot(S_Logistic_1mOppose_Biasfit)
S_Logistic_1mOppose_Bias_summary <- summary(S_Logistic_1mOppose_Biasfit, pars = c("tau", "m", "bias"), probs = c(0.1, 0.9))$summary
print(S_Logistic_1mOppose_Bias_summary)
get_posterior_mean(S_Logistic_1mOppose_Biasfit, pars=c('tau', 'm', 'bias'))[,5]
S_Logistic_1mOppose_Biasparams <- data.frame(Temp=get_posterior_mean(S_Logistic_1mOppose_Biasfit, pars=c('tau'))[,5],
                                             m=get_posterior_mean(S_Logistic_1mOppose_Biasfit, pars=c('m'))[,5],
                                             bias=get_posterior_mean(S_Logistic_1mOppose_Biasfit, pars=c('bias'))[,5],
                                             LL=get_posterior_mean(S_Logistic_1mOppose_Biasfit, pars=c('log_lik'))[,5])
k <- 3
S_Logistic_1mOppose_Biasparams$BIC <- log(model_data$nTrials) * k - 2 * (S_Logistic_1mOppose_Biasparams$LL)
S_Logistic_1mOppose_Biasparams$AIC <- 2 * k - 2 * (S_Logistic_1mOppose_Biasparams$LL)
S_Logistic_1mOppose_Bias_LL <- extract_log_lik(S_Logistic_1mOppose_Biasfit)
S_Logistic_1mOppose_Bias_LOO <- loo(S_Logistic_1mOppose_Bias_LL)
S_Logistic_1mOppose_Bias_WAIC <- waic(S_Logistic_1mOppose_Bias_LL)

############################

# MODEL VALIDATION

############################

print(loo_compare(list("Bias"=PB_LOO,
                       "SLogO.Bias"=S_Logistic_1mOppose_Bias_LOO
)),simplify = F
)

save.image("/Volumes/Research Project/Self-Anchoring/Study 1/Analysis/ModelFittingSA1_Environment.RData")

# Strict threshold

any(S_Logistic_1mOppose_Bias_summary[,7]>1.01)

S_Logistic_1mOppose_Bias_summary[which(S_Logistic_1mOppose_Bias_summary[,7]>1.01),]

# Liberal threshold

any(S_Logistic_1mOppose_Bias_summary[,7]>1.05)

S_Logistic_1mOppose_Biasparams$subID <- uIds
write.csv(S_Logistic_1mOppose_Biasparams, here::here("Study 1/Cleaning/output/S_Logistic_1m_Oppose_Biasparams.csv"), row.names=F)
arrow::write_parquet(S_Logistic_1mOppose_Biasparams, here::here("Study 1/Cleaning/output/S_Logistic_1m_Oppose_Biasparams.parquet"))
      
