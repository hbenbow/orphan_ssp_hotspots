setwd("~/Documents/Broad_spectrum_clusters/Raw_data/")
files=dir(pattern="RES")
all <- read.csv("~/Documents/Hotspots/Paper_version_4/wheat_all.csv", header=T)
colnames(all)<-c( "Chromosome", "start", "end", "GeneID","Score", "strand")
expression_scores<-as.data.frame(all$GeneID)
colnames(expression_scores)<-"Gene"
for(i in files){
  study<-paste(i)
  study<-substr(study, 5, nchar(study)-4)
  data<-read.csv(i)
  data=data[(data$padj<0.05),]
  data=data[(abs(data$log2FoldChange)>1),]
  genes<-data.frame(do.call('rbind', strsplit(data$Gene,'.',fixed=TRUE)))[,1]
  genes<-as.data.frame(unique(genes))
  genes$DE<-1
  colnames(genes)<-c("Gene", paste(study))
  expression_scores<-left_join(expression_scores, genes, by="Gene", all.x=T)
  expression_scores[is.na(expression_scores)]<-0
}


write.csv(expression_scores, file="~/Documents/Broad_spectrum_clusters/Raw_data/expression_scores.csv", row.names=F)
