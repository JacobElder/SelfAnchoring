


logFunc <- function(input, slope=1, shift=0){
  input = input - 4
  output = 1 / (1 + exp(-slope * (input - shift ) ))
  return(output)
}

Minputs=c(.5,.6,.7,.8,.9,1,1.25,1.5,2,2.5,3,3.5,4,4.5,5,5.5)
df<-data.frame(Minputs=sort(rep(Minputs, 7)), eval=rep(1:7) )
df$trans<-logFunc(df$eval, df$Minputs, shift = 0)

library(ggplot2)
FSplot<-ggplot(data=df, aes(x=eval, y=trans, group=as.factor(Minputs), color= as.factor(Minputs) )) +
  geom_line() + xlab("Original Self-Evaluation") + ylab("Probability of Group Normativity")  + theme_classic() +
  theme(legend.position="top") + labs(color="Growth Parameter")
FSplot
#ggsave("~/Documents/UC Riverside/Studies/Feedback Study/Manuscript/Figures/FS_demonstration.jpg", units = "in", dpi = 600, width = 4, height = 4)

Minputs <- c(4.407416, 2.48609)
df<-data.frame(Minputs=sort(rep(Minputs, 7)), eval=rep(1:7) )
df$trans<-logFunc(df$eval, df$Minputs, shift = 0)
FSplot<-ggplot(data=df, aes(x=eval, y=trans, group=as.factor(Minputs), color= as.factor(Minputs) )) +
  geom_line() + geom_point() + xlab("Original Self-Evaluation") + ylab("Probability of Group Normativity")  + theme_classic() +
  theme(legend.position="top") + labs(color="Growth Parameter")  + 
  scale_color_manual(labels=c('Outgroup','Ingroup'), values = wes_palette("Darjeeling1"))
FSplot





logFunc <- function(input, slope=1, shift=0, L=1){
  input = input - 4
  output = (L / (1 + exp(-slope * (input - shift ) )) ) + (1-L)/2
  return(output)
}

Minputs=c(.5,.6,.7,.8,.9,1,1.25,1.5,2,2.5,3,3.5,4,4.5,5,5.5)
df<-data.frame(Minputs=sort(rep(Minputs, 7)), eval=rep(1:7) )
df$trans<-logFunc(df$eval, df$Minputs, shift = 0, L = .5)

library(ggplot2)
FSplot<-ggplot(data=df, aes(x=eval, y=trans, group=as.factor(Minputs), color= as.factor(Minputs) )) +
  geom_line() + xlab("Original Self-Evaluation") + ylab("Probability of Group Normativity")  + theme_classic() +
  theme(legend.position="top") + labs(color="Growth Parameter")
FSplot

###

logFunc <- function(input, slope=1, shift=0, L=1, B=0){
  input = input - 4
  output = (L / (1 + exp(-slope * (input - shift ) )) ) + (B)
  return(output)
}

Minputs=c(.5,.6,.7,.8,.9,1,1.25,1.5,2,2.5,3,3.5,4,4.5,5,5.5)
df<-data.frame(Minputs=sort(rep(Minputs, 7)), eval=rep(1:7) )
df$trans<-logFunc(df$eval, df$Minputs, shift = 0, L = .75, B = .1)

library(ggplot2)
FSplot<-ggplot(data=df, aes(x=eval, y=trans, group=as.factor(Minputs), color= as.factor(Minputs) )) +
  geom_line() + xlab("Original Self-Evaluation") + ylab("Probability of Group Normativity")  + theme_classic() +
  theme(legend.position="top") + labs(color="Growth Parameter")
FSplot


DP <- sum(model_data$prevSelf[1,] * model_data$prevSim[1,1,])
logFunc <- function(input, slope=1, shift=0){
  #input = input - 4
  output = 1 / (1 + exp(-slope * (input - shift ) ))
  return(output)
}
logFunc(DP,shift=140:145)

logFunc <- function(sim, self, slope=1, shift=0){
  output = 1 / (1 + exp(-slope * (sim * (self - shift) ) ))
  return(output)
}
logFunc(rep(seq(from=0,to=.7,by=.05),7),1:7)
