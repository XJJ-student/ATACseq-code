rm(list=ls())
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
setwd("~/xjj/drug")
ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,!(colnames(ic501)%in%c("PC.100","PC.34"))]
colnames(ic50)<-sample_id[match(colnames(ic50),sample_id[,2]),1]
drug_info<-read.csv("compiled_drug_info.csv",header = TRUE,sep = ",")
drug_result<-drug_info[match(rownames(ic50),drug_info$drug_id),]

HDAC1<-drug_result[drug_result$target=="HDAC",]
ccc<-c("S2759","S1194","S1047")
HDAC2<-drug_result[drug_result$drug_id%in%ccc,]
HDAC3<-rbind(HDAC1,HDAC2)
HDAC<-HDAC3[order(HDAC3[,1],decreasing = T),]
HDAC_ic50<-ic50[match(as.character(HDAC[,1]),rownames(ic50)),]
library(pheatmap)
a=cor(t(HDAC_ic50))
result=pheatmap(a,scale = "none",main = "The correlation coefficients between HDACi and chemotherapy' AUC",show_rownames=T,show_colnames=T,
                clustering_distance_rows = "correlation",
                clustering_distance_cols = "correlation",clustering_method = "complete",cutree_rows=2,cutree_cols=2)

bb<-t(apply(HDAC_ic50,1,function(x) scale(x)))
colnames(bb)<-colnames(HDAC_ic50)
norma_result<-pheatmap(bb,scale = "none",main = "The AUC of HDACi and chemotherapy are normalized",show_rownames=T,show_colnames=T,
                       clustering_distance_rows = "correlation",
                       clustering_distance_cols = "euclidean",clustering_method = "ward.D2")
clusters_filter<-cutree(norma_result$tree_col,k=2)
sensitivity=colnames(HDAC_ic50[,clusters_filter==1])
resistance=colnames(HDAC_ic50[,clusters_filter==2])

sensitivity_AUC<-HDAC_ic50[,match(sensitivity,colnames(HDAC_ic50))]
resistance_AUC<-HDAC_ic50[,match(resistance,colnames(HDAC_ic50))]
### T.test
geneid<-rownames(HDAC_ic50)
tresult=matrix(0,length(geneid),4)
for (i in 1:length(geneid)){
  ttest=t.test(sensitivity_AUC[i,],resistance_AUC[i,],alternative = "two.sided")
  tresult[i,1:3]=c(geneid[i],ttest$statistic,ttest$p.value)
}
tresult[,4]=p.adjust(tresult[,3],method="BH")
colnames(tresult)<-c("geneid","statistic","p.value","FDR")
which(format(as.numeric(tresult[,4]), scientific=F)<0.05)
sum(format(as.numeric(tresult[,4]), scientific=F)<0.05)


### wilcox.test
tresult=matrix(0,length(geneid),4)
for (i in 1:length(geneid)){
  ttest=wilcox.test(as.numeric(sensitivity_AUC[i,]),as.numeric(resistance_AUC[i,]),alternative ="two.sided")
  tresult[i,1:3]=c(geneid[i],ttest$statistic,ttest$p.value)
}
tresult[,4]=p.adjust(tresult[,3],method="BH")
colnames(tresult)<-c("geneid","statistic","p.value","FDR")
which(format(as.numeric(tresult[,4]), scientific=F)<0.05)
sum(format(as.numeric(tresult[,4]), scientific=F)<0.05)
