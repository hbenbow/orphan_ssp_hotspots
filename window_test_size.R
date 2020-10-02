library(ggplot2)
require(zoo)
library(plyr)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# This is the threshold for a hotspot that we decide based on the permutation testing
threshold=0.3
# ==========================================================================================
# This section reads in the files
# which directory for the merged data
setwd("~/Documents/STB_hotspot/")
all <- read.csv("~/Documents/STB_hotspot/wheat_all.csv", header=T)
colnames(all)<-c( "Chromosome", "start", "end", "GeneID","Score", "strand")
expression_scores<-read.csv("~/Documents/STB_hotspot/expression_scores.csv", header=T)

dir.create("Gene_cluster_analysis")

clusters<-function(all, expression_scores, threshold, window){
  # ==========================================================================================
  # this section reads in all files and does the window analysis
  dir.create("Gene_cluster_analysis")
  diseases<-colnames(expression_scores)[-1]
  hotspots_together<-list()
  stress_together<-list()
  for (d in diseases){
    
    # =========================================================================================
    List<-list()
    for(chromosome in (paste(unique(all$Chromosome)))){
      bed<-all[(all$Chromosome == chromosome),]
      selection<-as.data.frame(expression_scores[,1])
      selection[[d]]<-expression_scores[[d]]
      colnames(selection)<-c("GeneID", "Disease")
      bed <- join(bed, selection, by="GeneID")
      bed[is.na(bed)]<-0
      bed$density<-rollapply(bed$Disease, width=window, FUN=mean, by.column=FALSE, fill=0)
      bed$Consecutive<-sequence(rle(as.character(bed$Disease))$lengths)
      bed$Stress<-paste(d)
      List[[length(List)+1]] = bed
    }
    df<-do.call(rbind.data.frame, List)
    # df is a data frame that is the same as the ben input file but it has the window analysis density colum
    # included and it has the name of the stress appended to a column
    
    # pick hotspots out of the large window analysis file
    
    d2<-list()
    Density<-list()
    for(i in 1:nrow(df)){
      data<-df[i,]
      if(data$density >= threshold){
        peak=paste(data$GeneID)
        min=i-4
        max=i+4
        hotspot<-df[(min:max),]
        Density[[length(Density)+1]]= hotspot
      }}
    density_together<-do.call(rbind.data.frame, Density)
    density_together<-density_together[!duplicated(density_together$GeneID),]
    
    if(length(density_together) >=1){
      density_together$Hotspot <- "Hotspot"
      hotspots_together[[length(hotspots_together)+1]]<-density_together
      write.csv(density_together, paste("Gene_cluster_analysis/", "density", d, ".csv", sep=""))
      df<-merge(df, density_together, all.x=T)
      df[is.na(df)]<-"Not hotspot"
      stress_together[[length(stress_together)+1]] = df
      df$Hotspot<-factor(df$Hotspot, levels=c("Not hotspot" , "Hotspot"))
    } else{ 
      df$Hotspot<-"Not hotspot"  
      stress_together[[length(stress_together)+1]] = df
    }
    
    plot<-
      ggplot(df[!(df$Chromosome=="chrUn"),], aes(x=end, y=density)) +geom_jitter(size=1.6, aes(colour=Hotspot), alpha=0.6) + facet_wrap(~Chromosome, ncol=3) + xlab("Position (bp)") + 
      ylab("Gene Density") + ggtitle(paste(d)) + theme_bw() + 
      geom_hline(yintercept=threshold, alpha=0.7) +
      theme(text = element_text(size=16, colour="black")) + 
      coord_cartesian(ylim=c(0,1))+
      scale_color_manual( values=c("grey60", "orangered2"))
    write.csv(df, file=paste("Gene_cluster_analysis/", d, "analysed.csv", sep="_"), row.names=F)
    ggsave(plot, file=paste("Gene_cluster_analysis/", d, ".pdf"), width=300, height=300, unit="mm", dpi=400)
    
  }
  
  all_stress_hotspots<-as.data.frame(do.call(rbind.data.frame,hotspots_together))
  all_stress_together<-as.data.frame(do.call(rbind.data.frame, stress_together))
  write.csv(all_stress_together, file="Gene_cluster_analysis/all_stress_together.csv", row.names=F)
  write.csv(all_stress_hotspots, file="Gene_cluster_analysis/all_stress_hotspots.csv", row.names=F)
  all_stress_together$Hotspot<-factor(all_stress_together$Hotspot, levels=c("Not hotspot", "Hotspot"))
  
  
  
  # ====================================================================================================
  for(i in unique(all_stress_together$Chromosome)){
    data<-all_stress_together[(all_stress_together$Chromosome==i),]
    write.csv(data, file=paste("Gene_cluster_analysis/", i, ".csv", sep=""), row.names=F)
    plot<-
      ggplot(data, aes(x=end, y=density)) +geom_jitter(size=1.6, aes(colour=Hotspot), alpha=0.6) +facet_grid(Stress~.) + xlab("Position (bp)") + 
      ylab("Gene Density") + ggtitle(paste(i)) + theme_bw() +
      geom_hline(yintercept=threshold, alpha=0.7) +
      theme(text = element_text(size=16, colour="black")) +
      scale_color_manual( values=c("grey60", "orangered2")) +
      coord_cartesian(ylim=c(0,1))
    ggsave(plot, file=paste("Gene_cluster_analysis/", i, "all_stress.pdf"), width=300, height=300, unit="mm", dpi=400)
  }
}

for(i in 10:20){
  dir.create(paste("window_", i, sep=""))
  setwd(paste("window_", i, sep=""))
  clusters(all, expression_scores, threshold, i)
}
