library(ggplot2)
require(zoo)
library(plyr)
library(dplyr)
library(tidyr)
library(RColorBrewer)
test blablalbal
options(scipen=999)
list<-list()
for(i in c(1:998, 1001)){
  data<-read.csv(paste(i, "/density.csv", sep=""))
  list[[length(list)+1]] = data
}

all<-as.data.frame(do.call(rbind.data.frame, list))

all<-aggregate(all$n, by=list(all$density, all$Stress), FUN=sum)
all$percentage<-NULL
colnames(all)<-c("density", "Stress", "n")
windows<-aggregate(all$n, by=list(all$Stress), FUN=sum)
colnames(windows)<-c("Stress", "total")
percent<-merge(all, windows, by="Stress")
percent$percentage<-percent$n/percent$total*100
percent<-percent[(percent$density > 0),]
plot<-ggplot(percent, aes(x=as.factor(density), y=percentage)) + geom_point(size=2) + 
  facet_wrap(~Stress, ncol=2) + theme_bw() + scale_y_log10() + geom_hline(yintercept=0.01) +
  xlab("Density") + ylab ("Percentage chance") +
  theme(text = element_text(size=15)) 
ggsave(plot, file="density.pdf")
write.csv(percent, "Density_percentage.csv")

list<-list()
for(i in 1:100){
  data<-read.csv(paste(i, "/consec.csv", sep=""))
  list[[length(list)+1]] = data
}

all<-as.data.frame(do.call(rbind.data.frame, list))
all<-aggregate(all$n, by=list(all$Consecutive, all$Stress), FUN=sum)
all$percentage<-NULL
colnames(all)<-c("consecutive", "Stress", "n")
windows<-aggregate(all$n, by=list(all$Stress), FUN=sum)
colnames(windows)<-c("Stress", "total")
percent<-merge(all, windows, by="Stress")
percent$percentage<-percent$n/percent$total*100
percent<-percent[(percent$consecutive > 0),]
plot<-ggplot(percent, aes(x=as.factor(consecutive), y=percentage)) + geom_point(size=2) + 
  facet_wrap(~Stress, ncol=2) + theme_bw() + scale_y_log10() + geom_hline(yintercept=0.01) +
  xlab("Consecutive") + ylab ("Percentage chance") +
  theme(text = element_text(size=15)) 
ggsave(plot, file="consecutive.pdf")
write.csv(percent, "consecutive_percent.csv")

