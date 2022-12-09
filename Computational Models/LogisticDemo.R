


logFunc <- function(input, slope=1, shift=0){
  input = input - 4
  output = 1 / (1 + exp(-slope * (input - shift ) ))
  return(output)
}

Minputs=c(.5,.6,.7,.8,.9,1,1.25,1.5,2,2.5,3,3.5,4,4.5,5,5.5)
df<-data.frame(Minputs=sort(rep(Minputs, 7)), eval=rep(1:7) )
df$trans<-logFunc(df$eval, df$Minputs, shift = 1)

library(ggplot2)
FSplot<-ggplot(data=df, aes(x=eval, y=trans, group=as.factor(Minputs), color= as.factor(Minputs) )) +
  geom_line() + xlab("Original Self-Evaluation") + ylab("Probability of Group Normativity")  + theme_classic() +
  theme(legend.position="top") + labs(color=Growth) 
FSplot
ggsave("~/Documents/UC Riverside/Studies/Feedback Study/Manuscript/Figures/FS_demonstration.jpg", units = "in", dpi = 600, width = 4, height = 4)
