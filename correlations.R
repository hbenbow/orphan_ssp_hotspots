library(gplots)
library(tidyr)
library(Hmisc)
library(corrplot)

my_palette <- colorRampPalette(c("red", "grey", "green"))(n = 50)
final_hotspots <- read.csv("~/Documents/Hotspots/Paper_version_4/final_hotspots.csv")
Fg_final <- read.csv("~/Documents/Hotspots/Paper_version_4/Fg_final.csv")
ERP013829_tpm <- read.delim("~/Documents/Hotspots/Paper_version_4/ERP013829_tpm.tsv")
dir.create("Correlations_FR_only")
setwd("Correlations_FR_only")

for(i in unique(final_hotspots$Hotspot)){
  data<-final_hotspots[final_hotspots$Hotspot== i,]
  peak<-round(nrow(data)/2)
  peak<-as.character(data[peak,5])
  peak2<-as.numeric(row.names(Fg_final[(Fg_final$GeneID==peak),]))
  min=peak2-10
  max=peak2+10
  hotspot<-Fg_final[(min:max),]
  hotspot_genes<-as.character(hotspot$GeneID)
  exp_data<-subset(ERP013829_tpm, ERP013829_tpm$gene %in% hotspot_genes)
  exp_data$colors<-as.character(hotspot$Colour)
  row.names(exp_data)<-exp_data$gene
  exp_data$gene<-NULL
  exp_data$sum<-rowSums(exp_data[,1:72])
  exp_data<-exp_data[!(exp_data$sum==0),]
  exp_data$row<-seq(1:nrow(exp_data))
  redmin<-min(exp_data[(exp_data$colors=="red"),75])
  redmax<-max(exp_data[(exp_data$colors=="red"),75])
  rows<-seq((redmin-1):(redmax+1))
  exp_data<-exp_data[(redmin-2):(redmax+2),]
  colors<-exp_data$colors
  exp_data$colors<-NULL
  exp_data$sum<-NULL
  exp_data$row<-NULL
  write.csv(exp_data, file=paste(i, "expression.csv"))
  exp_data<-as.matrix(t(exp_data))
  cor<-rcorr(exp_data)
  write.csv(cor$r, file=paste(i, "Rcor.csv"))
  write.csv(cor$P, file=paste(i, "Pcor.csv"))
  pdf(paste(i, ".pdf"))
  heatmap.2(cor$r, margins=c(15,15), Rowv = NA, Colv = NA, symm=T, revC=F, RowSideColors=colors, ColSideColors = colors, dendrogram = "none", trace = "none", scale="none", col=my_palette)
  dev.off()
  pdf(paste(i, "corr_plot.pdf"))
  corrplot(cor$r, type="lower", order="original",p.mat = cor$P, 
          sig.level = 0.05, insig = "blank", tl.col="black", tl.cex = 0.7, tl.srt = 45)
   dev.off()   
}
  
# do we want to change the size of the flanking region?
# maybe to equal the size of the hotspot?
kruskall<-list()
shapiro<-list()
stats<-list()
for(i in unique(final_hotspots$Hotspot)){
  # calculate correlation of hotspot only
  data<-final_hotspots[final_hotspots$Hotspot== i,]
  data<-data[(data$Fg.response==1),]
  hotspot_genes<-as.character(data$GeneID)
  exp_data<-subset(ERP013829_tpm, ERP013829_tpm$gene %in% hotspot_genes)
  row.names(exp_data)<-exp_data$gene
  exp_data$gene<-NULL
  exp_data$sum<-rowSums(exp_data)
  exp_data<-exp_data[!(exp_data$sum==0),]
  exp_data$sum<-NULL
  exp_data<-as.matrix(t(exp_data))
  cor<-rcorr(exp_data)
  cor_hotspot<-cor
  # calculate correlation of flank regions
  f1<-as.character(data$GeneID[1])
  f2<-as.character(data$GeneID[nrow(data)])
  flank1<-as.numeric(row.names(Fg_final[(Fg_final$GeneID==f1),]))
  flank2<-as.numeric(row.names(Fg_final[(Fg_final$GeneID==f2),]))
  min=flank1-30
  max=flank2+30
  flankregion<-Fg_final[(min:max),]
  flankregion<-flankregion[!(flankregion$GeneID %in% hotspot_genes),]
  flank_genes<-as.character(flankregion$GeneID)
  flank_exp_data<-subset(ERP013829_tpm, ERP013829_tpm$gene %in% flank_genes)
  row.names(flank_exp_data)<-flank_exp_data$gene
  flank_exp_data$gene<-NULL
  flank_exp_data$sum<-rowSums(flank_exp_data)
  flank_exp_data<-flank_exp_data[!(flank_exp_data$sum==0),]
  flank_exp_data$sum<-NULL
  flank_exp_data<-as.matrix(t(flank_exp_data))
  cor<-rcorr(flank_exp_data)
  cor_flank<-cor
  # convert correlations into long format
  cor_hotspot_long<-gather(as.data.frame(cor_hotspot$r), key="gene", value ="Cc")
  cor_hotspot_long$Data<-"Hotspot"
  cor_flank_long<-gather(as.data.frame(cor_flank$r), key="gene", value ="Cc")
  cor_flank_long$Data<-"Flank"
  all_cor<-rbind(cor_flank_long, cor_hotspot_long)
  all_cor$absolute<-abs(all_cor$Cc)
  all_cor<-all_cor[!(all_cor$absolute==1),]
  all_cor$hotspot<-paste(i)
  agg<-aggregate(all_cor$absolute, FUN=mean, by=list(all_cor$Data))
  agg$median<-aggregate(all_cor$absolute, FUN=median, by=list(all_cor$Data))[,2]
  agg$sd<-aggregate(all_cor$absolute, FUN=sd, by=list(all_cor$Data))[,2]
  agg$n<-aggregate(all_cor$absolute, FUN=length, by=list(all_cor$Data))[,2]
  agg$SE<-agg$sd/sqrt(agg$n)
  agg$hotspot<-paste(i)
  stats[[length(stats) + 1]] <-  agg
  write.csv(all_cor, file=paste(i, "hotspot_cc.csv"))
  # calculate statistics
  # check for normality
  all_cor$Data<-as.factor(all_cor$Data)
  shap<-as.data.frame(unlist(shapiro.test(all_cor$absolute)))
  shap<-as.data.frame(t(shap))
  shap$hotspot<-paste(i)
  shap<-shap[,c(2,5)]
  shapiro[[length(shapiro) + 1]] <-  shap
  krusk<-as.data.frame(unlist(kruskal.test(all_cor$absolute ~ all_cor$Data)))
  krus2<-as.data.frame(t(krusk))
  krus2$hotspot<-paste(i)
  krus2<-krus2[,c(3, 6)]
  kruskall[[length(kruskall) + 1]] <-  krus2
  plot<-
    ggplot(all_cor, aes(x=Data, y=absolute)) + 
    geom_boxplot(fill="grey60") + 
    theme_classic() + 
    xlab("Data set") + 
    ylab("Absolute correlation coefficient") + 
    theme(text =  element_text(size=15))
  ggsave(plot, filename = paste(i, "boxplot.pdf"))
}

kruskall_test <- do.call("rbind", kruskall)
write.csv(kruskall_test, file="kruskall-wallis.csv")
shapiro_test <- do.call("rbind", shapiro)
write.csv(shapiro_test, "shapiro.csv")
stats<-do.call("rbind", stats)
write.csv(stats, "Stats.csv")

