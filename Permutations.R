library(ggplot2)
require(zoo)
library(plyr)
library(dplyr)
library(tidyr)
options(scipen = 999)

# Read in gene coordinate file, containing Chromosome name, gene start position, gene end
# position, gene id, and strand. Do not need strand.

all <- read.csv("/home/hbenbow/Broad_spectrum_hotspots/wheat_all.csv", header=T)
colnames(all)<-c( "Chromosome", "start", "end", "GeneID", "Score", "strand")
# Read in expression scores. This is a .csv file in which each row is a gene
# an d order must correspond to the order in the gene coordinate file. Each column after 
# the gene column is 
expression_scores<-read.csv("/home/hbenbow/Broad_spectrum_hotspots/expression_scores.csv", header=T)


diseases<-colnames(expression_scores)[-1]
permutations=10
dir.create("Permutations")
# this section reads in all files and does the window analysis
# test
correlations_consec<-list() 
correlations_density<-list()
stress_together<-list()
conseq<-list()
conseq2<-list()
positives2<-list()
Density<-list()
predictionsc<-list()
predictionsd<-list()

for (d in diseases){
  perm = list() 
  List<-list()
  selection<-as.data.frame(expression_scores[,1])
  selection[[d]]<-expression_scores[[d]]
  colnames(selection)<-c("GeneID", "Disease")
  bed <- join(all, selection, by="GeneID")
  bed[is.na(bed)]<-0
  shuff<-transform(bed, Disease=sample(Disease))
  shuff$density<-rollapply(shuff$Disease, width=10, FUN=mean, by.column=FALSE, fill=0)
  shuff$Consecutive<-sequence(rle(as.character(shuff$Disease))$lengths)
  shuff[is.na(shuff)]<-0
  perm[[length(perm)+1]] = shuff
  
  for(i in seq_along(perm)){
    shuffled<-perm[[i]]
    positives<-shuffled
    positives<-as.data.frame(positives$density)
    colnames(positives)<-"density"
    counts_density<-count(positives, density)
    positives2[[length(positives2)+1]] = counts_density
    consec<-shuffled[!(shuffled$Disease == 0.0 ),]
    counts_consec<-count(consec, Consecutive)
    conseq[[length(conseq)+1]] = counts_consec
  }
  
  positives<-as.data.frame(do.call(rbind.data.frame, positives2))
  counts_density<-aggregate(positives$n, by=list(positives$density), FUN=sum)
  colnames(counts_density)<-c("density", "n")
  consecs<-as.data.frame(do.call(rbind.data.frame, conseq))
  counts_consec<-aggregate(consecs$n, by=list(consecs$Consecutive), FUN=sum)
  colnames(counts_consec)<-c("Consecutive", "n")
  counts_consec$Stress<-paste(d)
  counts_density$Stress<-paste(d)
  conseq2[[length(conseq2)+1]] = counts_consec
  Density[[length(Density)+1]] = counts_density
  
}

densitys<-do.call(rbind.data.frame, Density)
all<-aggregate(densitys$n, by=list(densitys$density, densitys$Stress), FUN=sum)
all$percentage<-NULL
colnames(all)<-c("density", "Stress", "n")
windows<-aggregate(all$n, by=list(all$Stress), FUN=sum)
colnames(windows)<-c("Stress", "total")
percent<-merge(all, windows, by="Stress")
percent$percentage<-percent$n/percent$total*100
percent<-percent[(percent$density > 0),]
ggplot(percent, aes(x=as.factor(density), y=percentage)) + geom_point(size=2) + 
  facet_wrap(~Stress) + theme_bw() + scale_y_log10() + geom_hline(yintercept=0.01) +
  xlab("Density") + ylab ("Percentage chance") +
  theme(text = element_text(size=15)) 

conseqs<-do.call(rbind.data.frame, conseq2)
all<-aggregate(conseqs$n, by=list(conseqs$Consecutive, conseqs$Stress), FUN=sum)
all$percentage<-NULL
colnames(all)<-c("consecutive", "Stress", "n")
windows<-aggregate(all$n, by=list(all$Stress), FUN=sum)
colnames(windows)<-c("Stress", "total")
percent<-merge(all, windows, by="Stress")
percent$percentage<-percent$n/percent$total*100
percent<-percent[(percent$consecutive > 0),]
ggplot(percent, aes(x=as.factor(consecutive), y=percentage)) + geom_point(size=2) + 
  facet_wrap(~Stress) + theme_bw() + scale_y_log10() + geom_hline(yintercept=0.01) +
  xlab("Density") + ylab ("Percentage chance") +
  theme(text = element_text(size=15)) 


write.csv(conseqs, file="Permutations/consec.csv", row.names=F)
write.csv(densitys, file="Permutations/density.csv", row.names=F)

