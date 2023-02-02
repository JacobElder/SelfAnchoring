shep <- function(self, sim, c=1, p=1){
  ss <- sum(self*sim)/sum(self)
  ess <- exp(-c * (1/ss)^p)
  return(ess)
}

shep <- function(self, sim, c=1, p=1){
  simshep <- exp(-c * (1/ss)^p)
  output <- sum(self*simshep)/sum(self)
  return(output)
}

set.seed(25)
shep(rep(2,10),runif(10, .6, .9),c=.7,p=1)
set.seed(25)
shep(rep(2,10),runif(10, .1, 4),c=.5,p=1)

ssins <- seq(from=.06,to=.40, by=.02)
shep <- function(simin, c=1, p=1){
  simshep <- exp(-c * (1/simin)^p)
  return(simshep)
}

Cinputs = seq(from=.1,to=2,by=.1)
Pinputs = seq(from=.25,to=2,by=.25)
ssins <- seq(from=.06,to=.40, by=.02)
df <- expand_grid(Ss = ssins, Cs = Cinputs, Ps = Pinputs)
df <- df[order(df$Ss, df$Cs, df$Ps),]

df$trans <- shep(df$Ss, df$Cs, df$Ps)

FSplot<-ggplot(data=df, aes(x=Ss, y=trans, group=as.factor(Cs), color= as.factor(Cs) )) +
  geom_line() + xlab("Similarity-to-Self") + ylab("Group Generalization")  + theme_classic() +
  theme(legend.position="top") + labs(color="Sensitivity Parameter, C")  + 
  facet_wrap(~ Ps)
FSplot

# #Minputs=seq(from=.1,to=5,by=.25)
# sim <- seq(from=0,to=1,by=.1)
# #df<-data.frame(Minputs=sort(rep(Minputs, 7)), eval=rep(1:7) )
# df<-data.frame(Minputs=sort(rep(sim, 7)), eval=rep(1:7) )
# df$trans<-shep(self = df$eval, sim = df$sim, c=.25, p = 1)

library(ggplot2)
# Softmax choice function
logsumexp <- function (x) {
  y <- max(x)
  y + log(sum(exp(x - y)))
}
softmax <- function (x) {
  exp(x - logsumexp(x))
}

ssins <-seq(from=.05,to=.40,by=.001)
Pin <- 1.269011
Cin <- 0.7088053
Temp <- 7.91
indf <- data.frame(Prob=apply(cbind(shep(ssins,c=Cin,p=Pin), shep(ssins,c=Cout,p=Pout)), 1, function(x) softmax(x*Temp))[1,],trans=shep(ssins,c=Cin,p=Pin), Ss=ssins, C, P, group="ingroup")
Pout <- 1.426963
Cout <- 0.8405059
outdf <- data.frame(Prob=apply(cbind(shep(ssins,c=Cin,p=Pin), shep(ssins,c=Cout,p=Pout)), 1, function(x) softmax(x*Temp))[2,],trans=shep(ssins,c=C,p=P), Ss=ssins, C, P, group="outgroup")
df<-rbind(indf,outdf)
colnames(df) <- c("Prob","trans","Ss","Cs","Ps","group")


FSplot<-ggplot(data=df, aes(x=Ss, y=Prob, group=as.factor(group), color= as.factor(group) )) +
  geom_line() + xlab("Similarity-to-Self") + ylab("Group Generalization")  + theme_classic() +
  theme(legend.position="top") + labs(color="Group")
FSplot
