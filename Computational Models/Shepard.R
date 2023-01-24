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
FSplot<-ggplot(data=df, aes(x=eval, y=trans, group=as.factor(Minputs), color= as.factor(Minputs) )) +
  geom_line() + xlab("Original Self-Evaluation") + ylab("Probability of Group Normativity")  + theme_classic() +
  theme(legend.position="top") + labs(color="Growth Parameter")
FSplot
