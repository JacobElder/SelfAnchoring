S_Logistic_1mOppose_Bias_L.Simulate <- function(simMat, trainTraits, testTraits, params){
  m<-params[1]
  bias<-params[2]
  tau<-param[3]
  L<-param[4]
  p = c()
  SEs <- round( runif(length(trainTraits), min=1, max=7) )
  GPin = L/(1 + exp((-m)*(SEs-4))) + ((1-L)/2)
  GPout = L/(1 + exp((m)*(SEs-4))) + ((1-L)/2)
  
  dataframe = matrix(nrow=testTraits,ncol=1)
  
  for (t in 1:length(testTraits)) {
    
    PS = simMat[trainTraits,t]
    
    simW[1] = dot_product(GPout[1:nTrain[s]],PS)
    simW[2] = dot_product(GPin[1:nTrain[s]],PS)
    
    p[1] = ( (1 - bias[s]) * ( simW[1]^tau[s] ) ) / ( ( (1 - bias[s]) *  ( simW[1]^ tau[s] ) ) + ( (bias[s]) * ( simW[2]^tau[s] ) ) )
    p[2] = ( (bias[s]) * ( simW[2]^tau[s] ) ) / ( ( (1 - bias[s]) *  ( simW[1]^ tau[s] ) ) + ( (bias[s]) * ( simW[2]^tau[s] ) ) )
    
    groupChoice[t,] = sample(c(1,2),1,p)
    
  }
  
}