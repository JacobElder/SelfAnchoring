# To test code:
# self <- round(runif(91,min=1,max=7))
# trainIdx <- sample(1:148, 91, replace=F)
# testIdx <- sample(1:148, 148, replace=F)
# param <- c(6, .8, .6)

PR_S_Logistic_1mOppose_Bias <- function(self, trainIdx, testIdx, simMat, param){
   
   tau <- param[1]
   m <- param[2]
   bias <- param[3]
   
   GPin = 1/(1 + exp((-m)*(self-4))); # logistic function for transforming eval ratings to ingroup
   GPout = 1/(1 + exp((m)*(self-4))); # logistic function for transforming eval ratings to outgroup
   
   simW <- c(NA,NA)
   prob <- c(NA,NA)
   
   process <- function(t){
     
     PS = simMat[trainIdx, testIdx[t]]
     
     simW[1] = sum(GPout*PS); # dot-product of ingroup and sim
     simW[2] = sum(GPin*PS); #dot-product of outgroup and sim
     
     prob[1] = ( (1 - bias) * ( simW[1]^tau ) ) / ( ( (1 - bias) *  ( simW[1]^tau ) ) + ( (bias) * ( simW[2]^tau ) ) ); # convert to probabilities
     prob[2] = ( (bias) * ( simW[2]^tau ) ) / ( ( (1 - bias) *  ( simW[1]^tau ) ) + ( (bias) * ( simW[2]^tau ) ) ); # convert to probabilities
     
     groupChoice = sample(c(1,2),size=1,prob=prob);
     
   }
   
   sapply(1:148, function(x) 
          process(x)
          )
   
}
