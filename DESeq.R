setwd("~/Documents/STB_hotspot/")
# DESeq2 differentially expression testing
library(DESeq2)
i="ERP009837"

data<-read.table(paste(i, "_count.tsv.gz", sep=""), header=T)
names<-data$transcript
data[1]<-NULL
data<-apply( data, 2, as.integer)
rownames(data)<-names
order<-as.data.frame(colnames(data))
colnames(order)<-"run_accession"
meta<-metadata[(metadata$secondary_study_accession == i),]
meta<-merge(order, meta, by="run_accession", sort=F)

genotypes<-paste(unique(meta$Variety))
Treatments<-paste(unique(meta$Stress.disease))
paste(unique(genotypes))
paste(unique(Treatments))
paste(unique(meta$Factor))

riband_meta <- read.csv("~/Documents/STB_hotspot/riband_meta.csv", row.names=1)
# if(length(genotypes) == 1){
dds <- DESeqDataSetFromMatrix(countData = data, colData = riband_meta, design = ~ Stress.disease)
dds<-DESeq(dds)


RES1<-as.data.frame(results(dds, contrast=c("Stress.disease", "Zymoseptoria.tritici.inoculation.1.day", "mock.inoculation.1.day"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES2<-as.data.frame(results(dds, contrast=c("Stress.disease", "Zymoseptoria.tritici.inoculation.4.days", "mock.inoculation.4.days"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES3<-as.data.frame(results(dds, contrast=c("Stress.disease", "Zymoseptoria.tritici.inoculation.9.days", "mock.inoculation.9.days"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES4<-as.data.frame(results(dds, contrast=c("Stress.disease", "Zymoseptoria.tritici.inoculation.14.days", "mock.inoculation.14.days"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES5<-as.data.frame(results(dds, contrast=c("Stress.disease", "Zymoseptoria.tritici.inoculation.21.days", "mock.inoculation.21.days"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))

RES1$Timepoint<-1
RES2$Timepoint<-4
RES3$Timepoint<-9
RES4$Timepoint<-14
RES5$Timepoint<-21

all<-rbind(RES1, RES2, RES3, RES4, RES5)
all_sig<-all[(all$padj<0.05),]
all_sig<-all_sig[abs(all_sig$log2FoldChange)>.5,]
all_sig<-na.omit(all_sig)

write.csv(all_sig, file="~/Documents/STB_hotspot/riband_STB_expression_scores.csv")


