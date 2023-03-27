library("rstan") # observe startup messages
library("tidyverse")
library(doParallel)
library(igraph)
library(loo)
library(here)

here::i_am("./Study 3/Computational Model/SA3_Modeling_Final.R")

fulldf <- as.data.frame( arrow::read_parquet(here("Study 3/Cleaning/output/fullTest.parquet")) )
fulldf <- fulldf[!is.na(fulldf$ingChoiceN),]

traindf <- as.data.frame( arrow::read_parquet(here("Study 3/Cleaning/output/fullTrain.parquet")) )
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

uConds <- unique(fulldf$condition)
for(c in uConds){
  conddf <- subset(fulldf, condition==c)
  uIds<-unique(conddf$subID)
  maxSubjs=length(uIds)
  maxTrials=max(conddf$trialTotal)
  
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
    df <- subset(conddf, subID==uIds[i])
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
  
  assign(paste0(gsub(" ","",c),"_model_data"),model_data)
  assign(paste0(gsub(" ","",c),"_model_data2"),model_data2)
  
}

for(i in uConds){
  model_data <- get(paste0(i,"_model_data"))

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
PBparams$BIC <- log(model_data$nTrials) * k - 2 * (PBparams$LL)
PBparams$AIC <- 2 * k - 2 * (PBparams$LL)
PB_LL <- extract_log_lik(PBfit)
PB_LOO <- loo(PB_LL)
PB_WAIC <- waic(PB_LL)

assign(paste0("PBfit.",i), PBfit)
assign(paste0("PB_summary.",i), PB_summary)
assign(paste0("PBparams.",i), PBparams)
assign(paste0("PB_LL.",i), PB_LL)
assign(paste0("PB_LOO.",i), PB_LOO)
assign(paste0("PB_WAIC.",i), PB_WAIC)

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

assign(paste0("S_Logistic_1m_Oppose_Biasfit.",i), S_Logistic_1mOppose_Biasfit)
assign(paste0("S_Logistic_1m_Oppose_Bias_summary.",i), S_Logistic_1mOppose_Bias_summary)
assign(paste0("S_Logistic_1m_Oppose_Biasparams.",i), S_Logistic_1mOppose_Biasparams)
assign(paste0("S_Logistic_1m_Oppose_Bias_LL.",i), S_Logistic_1mOppose_Bias_LL)
assign(paste0("S_Logistic_1m_Oppose_Bias_LOO.",i), S_Logistic_1mOppose_Bias_LOO)
assign(paste0("S_Logistic_1m_Oppose_Bias_WAIC.",i), S_Logistic_1mOppose_Bias_WAIC)

#############################

## SELF-TO-GROUP LOGISTIC MAPPING WITH PROBABILITY BIAS AND CONSTRAINT MODEL ##

#############################

# iter=3000
# warmup=floor(iter/2)
# modelFile <- here("Computational Models/S_Logistic_1mOppose_Bias_L.stan")
# cores<-detectCores()
# S_Logistic_1mOppose_Bias_Lfit <- stan(modelFile, data = model_data, iter = iter, warmup = warmup, cores = cores-1, seed = 52, control = list(max_treedepth = 12, adapt_delta = 0.95))
# traceplot(S_Logistic_1mOppose_Bias_Lfit)
# S_Logistic_1mOppose_Bias_L_summary <- summary(S_Logistic_1mOppose_Bias_Lfit, pars = c("tau", "m", "bias", "L"), probs = c(0.1, 0.9))$summary
# print(S_Logistic_1mOppose_Bias_L_summary)
# get_posterior_mean(S_Logistic_1mOppose_Bias_Lfit, pars=c('tau', 'm', 'bias', 'L'))[,5]
# S_Logistic_1mOppose_Bias_Lparams <- data.frame(Temp=get_posterior_mean(S_Logistic_1mOppose_Bias_Lfit, pars=c('tau'))[,5],
#                                              m=get_posterior_mean(S_Logistic_1mOppose_Bias_Lfit, pars=c('m'))[,5],
#                                              bias=get_posterior_mean(S_Logistic_1mOppose_Bias_Lfit, pars=c('bias'))[,5],
#                                              L=get_posterior_mean(S_Logistic_1mOppose_Bias_Lfit, pars=c('L'))[,5],
#                                              LL=get_posterior_mean(S_Logistic_1mOppose_Bias_Lfit, pars=c('log_lik'))[,5])
# k <- 4
# S_Logistic_1mOppose_Bias_Lparams$BIC <- log(model_data$nTrials) * k - 2 * (S_Logistic_1mOppose_Bias_Lparams$LL)
# S_Logistic_1mOppose_Bias_Lparams$AIC <- 2 * k - 2 * (S_Logistic_1mOppose_Bias_Lparams$LL)
# S_Logistic_1mOppose_Bias_L_LL <- extract_log_lik(S_Logistic_1mOppose_Bias_Lfit)
# S_Logistic_1mOppose_Bias_L_LOO <- loo(S_Logistic_1mOppose_Bias_L_LL)
# S_Logistic_1mOppose_Bias_L_WAIC <- waic(S_Logistic_1mOppose_Bias_L_LL)
# 
# assign(paste0("S_Logistic_1m_Oppose_Bias_Lfit.",i), S_Logistic_1mOppose_Bias_Lfit)
# assign(paste0("S_Logistic_1m_Oppose_Bias_L_summary.",i), S_Logistic_1mOppose_Bias_L_summary)
# assign(paste0("S_Logistic_1m_Oppose_Bias_Lparams.",i), S_Logistic_1mOppose_Bias_Lparams)
# assign(paste0("S_Logistic_1m_Oppose_Bias_L_LL.",i), S_Logistic_1mOppose_Bias_L_LL)
# assign(paste0("S_Logistic_1m_Oppose_Bias_L_LOO.",i), S_Logistic_1mOppose_Bias_L_LOO)
# assign(paste0("S_Logistic_1m_Oppose_Bias_L_WAIC.",i), S_Logistic_1mOppose_Bias_L_WAIC)

############################

# MODEL VALIDATION

############################

comparison <- print(loo_compare(list("Bias"=PB_LOO,
                       "SLogO.Bias"=S_Logistic_1mOppose_Bias_LOO
#                       "SLogO.Bias.L"=S_Logistic_1mOppose_Bias_L_LOO
)),simplify = F
)

assign(paste0("comparison",i),comparison)

}

save.image("/Volumes/Research Project/Self-Anchoring/Study 3/Analysis/ModelFittingSA3_Environment.RData")

any(S_Logistic_1m_Oppose_Bias_summary.Majority[,7]>1.01)
any(S_Logistic_1m_Oppose_Bias_summary.Minority[,7]>1.01)
S_Logistic_1m_Oppose_Bias_summary.Minority[which(S_Logistic_1m_Oppose_Bias_summary.Minority[,7]>1.01),]

write.csv(S_Logistic_1m_Oppose_Biasparams.Majority, here("Study 3/Cleaning/output/S_Logistic_1m_Oppose_Biasparams.Majority.csv"), row.names=F)
arrow::write_parquet(S_Logistic_1m_Oppose_Biasparams.Majority, here("Study 3/Cleaning/output/S_Logistic_1m_Oppose_Biasparams.Majority.parquet"))
write.csv(S_Logistic_1m_Oppose_Biasparams.Minority, here("Study 3/Cleaning/output/S_Logistic_1m_Oppose_Biasparams.Minority.csv"), row.names=F)
arrow::write_parquet(S_Logistic_1m_Oppose_Biasparams.Minority, here("Study 3/Cleaning/output/S_Logistic_1m_Oppose_Biasparams.Minority.parquet"))

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

S_Logistic_1m_Oppose_Bias_LL.CSULA


#########
# GROUP DIFFERENCES
#########

Majority_group_m_posterior <- (pnorm(rstan::extract(S_Logistic_1m_Oppose_Biasfit.Majority, pars = "mu_pr")$mu_pr[,2])*10)
Minority_group_m_posterior <- (pnorm(rstan::extract(S_Logistic_1m_Oppose_Biasfit.Minority, pars = "mu_pr")$mu_pr[,2])*10)

hBayesDM::plotHDI(Majority_group_m_posterior - Minority_group_m_posterior)

Majority_group_tau_posterior <- (pnorm(rstan::extract(S_Logistic_1m_Oppose_Biasfit.Majority, pars = "mu_pr")$mu_pr[,1])*10)
Minority_group_tau_posterior <- (pnorm(rstan::extract(S_Logistic_1m_Oppose_Biasfit.Minority, pars = "mu_pr")$mu_pr[,1])*10)

hBayesDM::plotHDI(Majority_group_tau_posterior - Minority_group_tau_posterior)

Majority_group_bias_posterior <- (pnorm(rstan::extract(S_Logistic_1m_Oppose_Biasfit.Majority, pars = "mu_pr")$mu_pr[,3]))
Minority_group_bias_posterior <- (pnorm(rstan::extract(S_Logistic_1m_Oppose_Biasfit.Minority, pars = "mu_pr")$mu_pr[,3]))

hBayesDM::plotHDI(Majority_group_bias_posterior - Minority_group_bias_posterior)


# ROPE Range
ROPErange <- function(posterior){
  rang <- c(-.1 * sd(posterior), .1 * sd(posterior) )
  return(rang)
}

# for(i in 1:length(paramnames)){
#   diff <- demPost[,i] - repPost[,i]
#   assign(paste0(paramnames[i],"Diff"),diff)
#   hdiofdiff <- hBayesDM::plotHDI(diff)
#   assign(paste0(paramnames[i],"_HDI"),diff)
#   ro <- bayestestR::rope(diff,range=ROPErange(diff))
#   assign(paste0(paramnames[i],"_ROPE"),ro)
#   pd <- bayestestR::p_direction(diff)
#   assign(paste0( paramnames[i],"_pd"), pd )
#   
# }

      