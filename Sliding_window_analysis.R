library(ggplot2)
require(zoo)
library(plyr)
library(dplyr)
library(tidyr)
library(RColorBrewer)

#  ID thresholds for each stress
density_thresholds<-read.csv("~/Documents/Orphan_ssp_hotspots/Data/Density_percentage.csv")
consec_thresholds<-read.csv("~/Documents/Orphan_ssp_hotspots/Data/consecutive_percent.csv")
dt<-list()
ct<-list()
for(i in unique(density_thresholds$Stress)){
  data=density_thresholds[(density_thresholds$Stress==i),]
  min = data[(data$percentage <=0.05),]
  if(nrow(min) > 0){
    min=as.data.frame(min(min$density))
    min$stress<-paste(i)
    colnames(min)<-c("min", "stress")
    dt[[length(dt)+1]] = min
  }else{
    min=as.data.frame(0)
    min$stress<-paste(i)
    colnames(min)<-c("min", "stress")
    ct[[length(ct)+1]] = min
  }
}

for(i in unique(density_thresholds$Stress)){
  data=consec_thresholds[(consec_thresholds$Stress==i),]
  min = data[(data$percentage <=0.05),]
  if(nrow(min) > 1){
    min=as.data.frame(min(min$consecutive))
    min$stress<-paste(i)
    colnames(min)<-c("min", "stress")
    ct[[length(ct)+1]] = min
  }else{
    min=as.data.frame(max(data$consecutive))
    min$stress<-paste(i)
    colnames(min)<-c("min", "stress")
    ct[[length(ct)+1]] = min
  }
}
den_thr<-do.call(rbind.data.frame, dt)
consec_thr<-do.call(rbind.data.frame, ct)

# ==========================================================================================
# This section reads in the files
# which directory for the merged data
setwd("~/Documents/Orphan_ssp_hotspots/")
all <- read.csv("~/Documents/Hotspots/Paper_version_4/wheat_all.csv", header=T)
colnames(all)<-c( "Chromosome", "start", "end", "GeneID","Score", "strand")
expression_scores<-read.csv("~/Documents/Orphan_ssp_hotspots/Data/Wheat_SSPs.csv", header=T)
dir.create("Gene_cluster_analysis")

clusters<-function(all, expression_scores){
  # ==========================================================================================
  # this section reads in all files and does the window analysis
  diseases<-colnames(expression_scores)[-1]
  hotspots_together<-list()
  stress_together<-list()
  no_hots<-list()
  for (d in diseases){
    thres_d<-den_thr[(den_thr$stress==d),1]
    thres_c<-as.numeric(consec_thr[(consec_thr$stress==d),1])
    # =========================================================================================
    List<-list()
    for(chromosome in (paste(unique(all$Chromosome)))){
      bed<-all[(all$Chromosome == chromosome),]
      selection<-as.data.frame(expression_scores[,1])
      selection[[d]]<-expression_scores[[d]]
      colnames(selection)<-c("GeneID", "Disease")
      bed <- join(bed, selection, by="GeneID")
      bed[is.na(bed)]<-0
      bed$density<-rollapply(bed$Disease, width=10, FUN=mean, by.column=FALSE, fill=0)
      bed$Consecutive<-sequence(rle(as.character(bed$Disease))$lengths)
      bed$Stress<-paste(d)
      List[[length(List)+1]] = bed
    }
    df<-do.call(rbind.data.frame, List)
    # df is a data frame that is the same as the bed input file but it has the window analysis density colum
    # included and it has the name of the stress appended to a column
    
    # pick hotspots out of the large window analysis file
    
    d2<-list()
    Density<-list()
    for(i in 1:nrow(df)){
      data<-df[i,]
      if(data$density >= thres_d){
        peak=paste(data$GeneID)
        min=i-5
        max=i+5
        hotspot<-df[(min:max),]
        Density[[length(Density)+1]]= hotspot
      }}
    density_together<-do.call(rbind.data.frame, Density)
    density_together<-density_together[!duplicated(density_together$GeneID),]
    
    df$position<-as.numeric(row.names(df))
    nod<-df[(df$Disease==0),]
    yesd<-df[(df$Disease==1),]
    nod$Consecutive=0
    dd<-rbind(nod, yesd)
    df<-dd[order(dd$position),]
    
    c2<-list()
    consec<-list()
    for(i in 1:nrow(df)){
      data<-df[i,]
      if(data$Consecutive >= thres_c){
        peak=paste(data$GeneID)
        min=i-data$Consecutive
        max=i
        hotspot<-df[(min:max),]
        consec[[length(consec)+1]]= hotspot
      }}
    consec_together<-do.call(rbind.data.frame, consec)
    consec_together<-consec_together[!duplicated(consec_together$GeneID),]
    consec_together<-consec_together[(consec_together$Disease==1),]
    
    
    
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
      no_hots[[length(no_hots)+1]] = df
    }
    
    plot<-
      ggplot(df[!(df$Chromosome=="chrUn"),], aes(x=end, y=density)) +geom_jitter(size=1.6, aes(colour=Hotspot), alpha=0.6) + facet_wrap(~Chromosome, ncol=3) + xlab("Position (bp)") + 
      ylab("Gene Density") + ggtitle(paste(d)) + theme_bw() + 
      geom_hline(yintercept=thres_d, alpha=0.7) +
      theme(text = element_text(size=16, colour="black")) + 
      coord_cartesian(ylim=c(0,1.1))+
      scale_color_manual( values=c("grey60", "orangered2"))
    write.csv(df, file=paste("Gene_cluster_analysis/", d, "analysed.csv", sep="_"), row.names=F)
    ggsave(plot, file=paste("Gene_cluster_analysis/", d, ".pdf"), width=400, height=300, unit="mm", dpi=400)
    
    
    
    
  }
  
  all_stress_hotspots<-as.data.frame(do.call(rbind.data.frame,hotspots_together))
  all_stress_together<-as.data.frame(do.call(rbind.data.frame, stress_together))
  write.csv(all_stress_together, file="Gene_cluster_analysis/all_stress_together.csv", row.names=F)
  write.csv(all_stress_hotspots, file="Gene_cluster_analysis/all_stress_hotspots.csv", row.names=F)
  all_stress_together$Hotspot<-factor(all_stress_together$Hotspot, levels=c("Not hotspot", "Hotspot"))


  # ====================================================================================================
  all_stress_together_graph<-all_stress_together
  for(i in unique(all_stress_together_graph$Chromosome)){
    data<-all_stress_together_graph[(all_stress_together_graph$Chromosome==i),]
    write.csv(data, file=paste("Gene_cluster_analysis/", i, ".csv", sep=""), row.names=F)
    if(length(unique(data$Hotspot)) >1){
      plot<-
        ggplot(data, aes(x=end, y=density)) +geom_jitter(size=1.6, aes(colour=Hotspot), alpha=0.6) +facet_grid(Stress~.) + xlab("Position (bp)") + 
        ylab("Gene Density") + ggtitle(paste(i)) + theme_bw() +
        theme(text = element_text(size=16, colour="black")) +
        scale_color_manual( values=c("grey60", "orangered2")) +
        coord_cartesian(ylim=c(0,1.2))
      ggsave(plot, file=paste("Gene_cluster_analysis/", i, "all_stress.pdf"), width=300, height=300, unit="mm", dpi=400)
    }else{
      plot<-
        ggplot(data, aes(x=end, y=density)) +geom_jitter(size=1.6, aes(colour=Hotspot), alpha=0.6) +facet_grid(Stress~.) + xlab("Position (bp)") + 
        ylab("Gene Density") + ggtitle(paste(i)) + theme_bw() +
        theme(text = element_text(size=16, colour="black")) +
        scale_color_manual( values=c("grey60")) +
        coord_cartesian(ylim=c(0,1.2))
      ggsave(plot, file=paste("Gene_cluster_analysis/", i, "all_stress.pdf"), width=300, height=300, unit="mm", dpi=400)
    }
  }
}
