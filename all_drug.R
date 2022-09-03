#######  根据药敏IC50数据将样本分为两类
rm(list=ls())
setwd("~/xjj/drug")
ic501<-read.csv("changhai_ic50.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,-(which(colnames(ic501)==c("PC.100","PC.34")))]
na_count<-apply(ic50,1,function(x) sum(is.na(x)))
d_ic50<-ic50[-which(na_count>30),]##去掉NA值大于30的
for(i in 1:nrow(d_ic50)){
  d_ic50[i,which(is.na(d_ic50[i,]))]<-(max(d_ic50[i,which(!is.na(d_ic50[i,]))]))*2
}

#没有做标准化
ic2<-log2(d_ic50+1)
tresult=matrix(0,length(rownames(ic2)),6)
for (i in 1:length(rownames(ic2))){
  median_value<-median(as.numeric(ic2[i,]))
  res<-ic2[i,which(ic2[i,]>median_value)]
  sen<-ic2[i,which(ic2[i,]<=median_value)]
  ttest=t.test(res,sen)#wilcox.test
  res_names<-colnames(ic2)[which(ic2[i,]>median_value)]
  sen_names<-colnames(ic2)[which(ic2[i,]<=median_value)]
  res_names1<-paste0(res_names, collapse = ",")
  sen_names1<-paste0(sen_names, collapse = ",")
  tresult[i,1:5]=c(rownames(ic2)[i],ttest$statistic,ttest$p.value,res_names1,sen_names1)
}
tresult[,6]=p.adjust(tresult[,3],method="BH")
colnames(tresult)<-c("drugid","statistic","p.value","res_names","sen_names","FDR")
sum(as.numeric(tresult[,6])<0.05)
#tresult_want<-tresult[as.numeric(tresult[,6])<0.05,]
tresult_want<-tresult[as.numeric(tresult[,2])>7,] #FC较大的留下
ic22<-ic2[match(tresult_want[,1],rownames(ic2)),]
setwd("~/xjj/drug/drug_result/durg62")
write.table(ic22,"drug42_IC50.txt",sep="\t",col.names=T,row.names=T)


which(tresult_want[,4]==tresult_want[-match(unique(tresult_want[,4]),tresult_want[,4]),4])
common<-tresult_want[which(tresult_want[,4]==tresult_want[-match(unique(tresult_want[,4]),tresult_want[,4]),4]),]
setwd("~/xjj/drug/drug_result/durg62")
write.table(tresult_want,"drug62_sample.txt",sep="\t",col.names=T,row.names=F)

setwd("~/xjj/drug/drug_result/durg62")
write.table(tresult_want,"drug42_sample.txt",sep="\t",col.names=T,row.names=F)


i3<-tresult_want[,1]
ic3<-ic2[match(i3,rownames(ic2)),]
library(pheatmap)
heatmap=pheatmap(ic3,scale = "none",main = "log2 IC50",show_rownames=T,show_colnames=T,
                 cluster_rows = TRUE,cluster_cols = TRUE, clustering_distance_rows = "correlation",
                 clustering_distance_cols = "correlation",clustering_method = "ward.D2")

common_sample<-tresult_want[match(common[,1],tresult_want[,1]),]
i3<-c("5-FU","GEM","IRI","OXA","PAC")
i3<-common[,1]
ic3<-ic2[match(i3,rownames(ic2)),]
heatmap=pheatmap(ic3,scale = "none",main = "log2 IC50",show_rownames=T,show_colnames=T,
                 cluster_rows = TRUE,cluster_cols = TRUE, clustering_distance_rows = "correlation",
                 clustering_distance_cols = "correlation",clustering_method = "ward.D2")
################
wresult=matrix(0,length(rownames(ic2)),6)
for (i in 1:length(rownames(ic2))){
  median_value<-median(as.numeric(ic2[i,]))
  res<-as.numeric(ic2[i,which(ic2[i,]>median_value)])
  sen<-as.numeric(ic2[i,which(ic2[i,]<=median_value)])
  wtest=wilcox.test(res,sen)#wilcox.test
  res_names<-colnames(ic2)[which(ic2[i,]>median_value)]
  sen_names<-colnames(ic2)[which(ic2[i,]<=median_value)]
  res_names1<-paste0(res_names, collapse = ",")
  sen_names1<-paste0(sen_names, collapse = ",")
  wresult[i,1:5]=c(rownames(ic2)[i],wtest$statistic,wtest$p.value,res_names1,sen_names1)
}
wresult[,6]=p.adjust(wresult[,3],method="BH")
colnames(wresult)<-c("drugid","statistic","p.value","res_names","sen_names","FDR")
sum(as.numeric(wresult[,6])<0.05)
wresult_want<-wresult[as.numeric(wresult[,6])<0.05,]
############## HDAC
drug_info<-read.csv("compiled_drug_info.csv",header = TRUE,sep = ",")
drug_result<-drug_info[match(rownames(d_ic50),drug_info$drug_id),]
HDAC1<-drug_result[drug_result$target=="HDAC",]
aaa<-c("S1848","S2759","S1194","S1047")
HDAC2<-drug_result[drug_result$drug_id%in%aaa,]
HDAC3<-rbind(HDAC1,HDAC2)
HDAC<-HDAC3[order(HDAC3[,1],decreasing = T),]

int_HDAC<-intersect(as.character(HDAC[,1]),tresult_want[,1])
HDAC_ic<-tresult_want[match(int_HDAC,tresult_want[,1]),]
HDAC_ic50<-d_ic50[match(HDAC_ic[,1],rownames(d_ic50)),]
ic2<-log2(HDAC_ic50+1)
library(pheatmap)
heatmap=pheatmap(ic2,scale = "none",main = "log2 IC50",show_rownames=T,show_colnames=T,
                 cluster_rows = TRUE,cluster_cols = TRUE, clustering_distance_rows = "correlation",
                 clustering_distance_cols = "euclidean",clustering_method = "ward.D2")
HDAC_ic501<-d_ic50[match("S8043",rownames(d_ic50)),]
HDAC_ic502<-d_ic50[match("S1194",rownames(d_ic50)),]
corre<-cor.test(as.numeric(HDAC_ic501),as.numeric(HDAC_ic502))
library(broom)
tidy(corre)$p.value[1]
tidy(corre)$estimate[1]

common_sample<-tresult_want[match("S8043",tresult_want[,1]),]
library(dplyr)
sensitivity<-unlist(strsplit(common_sample[5], ","))
resistance<-unlist(strsplit(common_sample[4], ","))

common_sample<-tresult_want[match("S1194",tresult_want[,1]),]
library(dplyr)
sensitivity1<-unlist(strsplit(common_sample[5], ","))
resistance1<-unlist(strsplit(common_sample[4], ","))

intersect(sensitivity,sensitivity1)
length(sensitivity)
length(sensitivity1)
length(intersect(sensitivity,sensitivity1))

intersect(resistance,resistance1)
length(resistance)
length(resistance1)
length(intersect(resistance,resistance1))

i3<-c("S1047","S8043","S1194")
ic3<-ic2[match(i3,rownames(ic2)),]
heatmap=pheatmap(ic3,scale = "none",main = "log2 IC50",show_rownames=T,show_colnames=T,
                 cluster_rows = TRUE,cluster_cols = TRUE, clustering_distance_rows = "correlation",
                 clustering_distance_cols = "correlation",clustering_method = "complete")


##########################  用0-1表示敏感和耐药 ########################
rm(list=ls())
setwd("~/xjj/drug")
ic501<-read.csv("changhai_ic50.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,-(which(colnames(ic501)==c("PC.100","PC.34")))]
na_count<-apply(ic50,1,function(x) sum(is.na(x)))
d_ic50<-ic50[-which(na_count>30),]##去掉NA值大于30的
for(i in 1:nrow(d_ic50)){
  d_ic50[i,which(is.na(d_ic50[i,]))]<-(max(d_ic50[i,which(!is.na(d_ic50[i,]))]))*10
}

ic2<-log2(d_ic50+1)

result<-matrix(2,nrow(ic2),ncol(ic2))
for (i in 1:length(rownames(ic2))){
  median_value<-median(as.numeric(ic2[i,]))
  result[i,which(ic2[i,]>median_value)]<-1
  result[i,which(ic2[i,]<=median_value)]<-0
}
heatmap=pheatmap(result,scale = "none",main = "log2 IC50",show_rownames=T,show_colnames=T,
                 cluster_rows = TRUE,cluster_cols = TRUE, clustering_distance_rows = "correlation",
                 clustering_distance_cols = "correlation",clustering_method = "ward.D2")

#没有做标准化
ic2<-log2(d_ic50+1)
tresult=matrix(0,length(rownames(ic2)),6)
for (i in 1:length(rownames(ic2))){
  median_value<-median(as.numeric(ic2[i,]))
  res<-ic2[i,which(ic2[i,]>median_value)]
  sen<-ic2[i,which(ic2[i,]<=median_value)]
  ttest=t.test(res,sen)#wilcox.test
  res_names<-colnames(ic2)[which(ic2[i,]>median_value)]
  sen_names<-colnames(ic2)[which(ic2[i,]<=median_value)]
  res_names1<-paste0(res_names, collapse = ",")
  sen_names1<-paste0(sen_names, collapse = ",")
  tresult[i,1:5]=c(rownames(ic2)[i],ttest$statistic,ttest$p.value,res_names1,sen_names1)
}
tresult[,6]=p.adjust(tresult[,3],method="BH")
colnames(tresult)<-c("drugid","statistic","p.value","res_names","sen_names","FDR")
sum(as.numeric(tresult[,6])<0.05)
tresult_want<-tresult[as.numeric(tresult[,6])<0.05,]



#标准化,最大最小值标准化
data_min_max<-apply(d_ic50,1,function(x) (x-min(x)) / (max(x)-min(x)))
ic2<-data_min_max
library(broom)
tresult=matrix(0,length(colnames(ic2)),6)
for (i in 1:length(colnames(ic2))){
  median_value<-median(as.numeric(ic2[,i]))
  res<-ic2[which(ic2[,i]>median_value),i]
  sen<-ic2[which(ic2[,i]<=median_value),i]
  ttest=t.test(res,sen,alternative = "two.sided")#wilcox.test
  res_names<-rownames(ic2)[which(ic2[,i]>median_value)]
  sen_names<-rownames(ic2)[which(ic2[,i]<=median_value)]
  res_names1<-paste0(res_names, collapse = ",")
  sen_names1<-paste0(sen_names, collapse = ",")
  tresult[i,1:5]=c(colnames(ic2)[i],tidy(ttest)$statistic[1],tidy(ttest)$p.value[1],res_names1,sen_names1)
}
tresult[,6]=p.adjust(tresult[,3],method="BH")
colnames(tresult)<-c("drugid","statistic","p.value","res_names","sen_names","FDR")
sum(as.numeric(tresult[,6])<0.05)
tresult_want<-tresult[as.numeric(tresult[,6])<0.05,]

which(tresult_want[,4]==tresult_want[-match(unique(tresult_want[,4]),tresult_want[,4]),4])
common<-tresult_want[which(tresult_want[,4]==tresult_want[-match(unique(tresult_want[,4]),tresult_want[,4]),4]),]

i3<-tresult_want[,1]
ic21<-t(ic2)
ic3<-ic21[match(i3,rownames(ic21)),]
heatmap=pheatmap(ic3,scale = "none",main = "log2 IC50",show_rownames=T,show_colnames=T,
                 cluster_rows = TRUE,cluster_cols = TRUE, clustering_distance_rows = "euclidean",
                 clustering_distance_cols = "euclidean",clustering_method = "ward.D")


i3<-c("S1047","S8043")
ic21<-t(ic2)
ic3<-ic21[match(i3,rownames(ic21)),]
heatmap=pheatmap(ic3,scale = "none",main = "log2 IC50",show_rownames=T,show_colnames=T,
                 cluster_rows = TRUE,cluster_cols = TRUE, clustering_distance_rows = "correlation",
                 clustering_distance_cols = "euclidean",clustering_method = "ward.D2")





###############################################################
###########   clinical information
rm(list=ls())
common_sample<-tresult_want[match(common[,1],tresult_want[,1]),]
library(dplyr)
sensitivity<-unlist(strsplit(common_sample[1,5], ","))
resistance<-unlist(strsplit(common_sample[1,4], ","))
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]
library("data.table")
clinical<-fread("HDAC_IC50_information.csv",header=T,data.table=F)
sen<-clinical[match(sensitivity_new,clinical$sample_new),]
res<-clinical[match(resistance_new,clinical$sample_new),]
median(sen$age)
median(res$age)
sum(res$sex==1)
sum(res$location==1)
sum(sen$`Degree of differentiation`==1)
sum(sen$`Neoadjuvant therapy`==0)
unique(sen$`Tumor stage`)


###############################################################
#############   mutation landscape
rm(list=ls())
setwd("~/xjj/WGS_CNV/mutation")
library("data.table")
mut<-fread("Organoid_mut_binary.txt",header=T,data.table=F)
mut1<-mut[,-1]

common_sample<-tresult_want[match(common[,1],tresult_want[,1]),]
library(dplyr)
sensitivity<-unlist(strsplit(common_sample[1,5], ","))
resistance<-unlist(strsplit(common_sample[1,4], ","))
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
resistance_new<-as.character(sample_id[match(resistance,sample_id[,2]),1])
sensitivity_new<-as.character(sample_id[match(sensitivity,sample_id[,2]),1])
res_int<-intersect(resistance_new,colnames(mut))
resistance_mut<-mut[,match(res_int,colnames(mut))]
sen_int<-intersect(sensitivity_new,colnames(mut))
sensitivity_mut<-mut[,match(sen_int,colnames(mut))]

muty<-cbind(resistance_mut,sensitivity_mut)
delete<-apply(muty,1,function(x) sum(x==1))
length(delete[which(delete>=10)])

geneid<-mut[,1]
geneid[which(delete>=10)]
data<-data.frame()
for(i in 1:length(geneid)){
  a<-sum(sensitivity_mut[i,]==1)
  c<-sum(sensitivity_mut[i,]==0)
  b<-sum(resistance_mut[i,]==1)
  d<-sum(resistance_mut[i,]==0)
  tmp<-matrix(c(a, c,b, d),
              nrow = 2,
              dimnames = list(Truth = c("mutation", "wild"),
                              Guess = c("sensitivity", "resistance")))
  P<-fisher.test(tmp,alternative = "two.sided",conf.level = 0.95)$p.value
  data[i,1]<-geneid[i]
  data[i,2]<-P
}
data[which(data[,2]<0.05),]  ###只有一个基因HYDIN
#############  突变谱中耐药相关基因
mutp<-cbind(sensitivity_mut,resistance_mut)
label<-c(rep(1,ncol(sensitivity_mut)),rep(0,ncol(resistance_mut)))
cor_mut_gene<-data.frame()
for(i in 1:nrow(mutp)){
  cot<-cor.test(as.numeric(mutp[i,]),label,method="spearman")
  cor_mut_gene[i,1]<-geneid[i]
  cor_mut_gene[i,2]<-as.numeric(cot[4])
  cor_mut_gene[i,3]<-as.numeric(cot[3])
}
cor_mut_gene[,4]<-p.adjust(cor_mut_gene[,3],method="BH")
colnames(cor_mut_gene)<-c("geneid","cor","Pvalue","FDR")
which(cor_mut_gene[,3]<0.05)



###############################################################
#############   CNV landscape
rm(list=ls())
library("data.table")
library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(stringr)
G="5-FU"
file_name="5-FU"

CNV<-function(G,file_name){
setwd("~/xjj/WGS_CNV/CNV_FACETS")
mut<-fread("CNV_FACETS_median_CNA_matrix.txt",header=T,data.table=F)
setwd("~/xjj/drug/drug_result/durg62")
tresult_want<-read.table("drug62_sample.txt",sep="\t",header=T)
common_sample<-as.data.frame(tresult_want[match(G,tresult_want[,1]),])
sensitivity<-unlist(strsplit(as.character(common_sample[1,5]), ","))
resistance<-unlist(strsplit(as.character(common_sample[1,4]), ","))
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
resistance_new<-as.character(sample_id[match(resistance,sample_id[,2]),1])
sensitivity_new<-as.character(sample_id[match(sensitivity,sample_id[,2]),1])
res_int<-intersect(resistance_new,colnames(mut))
resistance_mut<-mut[,match(res_int,colnames(mut))]
sen_int<-intersect(sensitivity_new,colnames(mut))
sensitivity_mut<-mut[,match(sen_int,colnames(mut))]
geneid<-mut[,1]
data_amp<-data.frame()
for(i in 1:length(geneid)){
  a<-sum(sensitivity_mut[i,]==1)
  c<-sum(sensitivity_mut[i,]!=1)
  b<-sum(resistance_mut[i,]==1)
  d<-sum(resistance_mut[i,]!=1)
  tmp<-matrix(c(a, c,b, d),
              nrow = 2,
              dimnames = list(Truth = c("amplification", "wild"),
                              Guess = c("sensitivity", "resistance")))
  P<-fisher.test(tmp,alternative = "two.sided",conf.level = 0.95)$p.value
  data_amp[i,1]<-geneid[i]
  data_amp[i,2]<-P
}
data_amp[,3]<-p.adjust(data_amp[,2],method="BH")
amp_gene<-data_amp[which(data_amp[,2]<0.05),]  
colnames(amp_gene)<-c("geneid","Pvalue","FDR")
cat(paste("The number of amp gene is ",nrow(amp_gene),sep=""),"\n")
cat(paste(amp_gene[,1],sep=""),"\n")
sensitivity1<-apply(sensitivity_mut,1,function(x) mean(x))
resistance1<-apply(resistance_mut,1,function(x) mean(x))
s<-match(amp_gene[,1],geneid)
cat(paste("sen-res>0 is ",sum(sensitivity1[s]-resistance1[s]>0),sep=""),"\n")

data_censored<-data.frame()
for(i in 1:length(geneid)){
  a<-sum(sensitivity_mut[i,]==(-1))
  c<-sum(sensitivity_mut[i,]!=(-1))
  b<-sum(resistance_mut[i,]==(-1))
  d<-sum(resistance_mut[i,]!=(-1))
  tmp<-matrix(c(a, c,b, d),
              nrow = 2,
              dimnames = list(Truth = c("censored", "wild"),
                              Guess = c("sensitivity", "resistance")))
  P<-fisher.test(tmp,alternative = "two.sided",conf.level = 0.95)$p.value
  data_censored[i,1]<-geneid[i]
  data_censored[i,2]<-P
}
data_censored[,3]<-p.adjust(data_censored[,2],method="BH")
censored_gene<-data_censored[which(data_censored[,2]<0.05),]  ###控制P值
colnames(censored_gene)<-c("geneid","Pvalue","FDR")
cat(paste("The number of censored gene is ",nrow(censored_gene),sep=""),"\n")
cat(paste(censored_gene[,1],sep=""),"\n")

sensitivity1<-apply(sensitivity_mut,1,function(x) mean(x))
resistance1<-apply(resistance_mut,1,function(x) mean(x))
s<-match(censored_gene[,1],geneid)
sum(sensitivity1[s]-resistance1[s]<0)
cat(paste("res-sen>0 is ",sum(resistance1[s]-sensitivity1[s]>0),sep=""),"\n")

int_CNV<-intersect(censored_gene[,1],amp_gene[,1])
censored_gene1<-censored_gene[-match(int_CNV,censored_gene[,1]),]#去掉既是扩增也是删失的基因
amp_gene1<-amp_gene[-match(int_CNV,amp_gene[,1]),]
cat(paste("The number of amp int censored is ",length(int_CNV),sep=""),"\n")
cat(paste("The number of final amp is ",nrow(amp_gene1),sep=""),"\n")
cat(paste("The number of final censored is ",nrow(censored_gene1),sep=""),"\n")
if (nrow(amp_gene1)>0){
  setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/7_CNV/fisher",sep=""))
  write.table(amp_gene1,paste(file_name,"sensitivity-resistance-amp-CNV-gene.txt",sep="-"),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
}
if (nrow(amp_gene1)>0){
  setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/7_CNV/fisher",sep=""))
  write.table(censored_gene1,paste(file_name,"sensitivity-resistance-censored-CNV-gene.txt",sep="-"),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
}
######### CNV KEGG
genea=bitr(amp_gene1[,1],fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
kegga<-enrichKEGG(genea[,2],organism = 'hsa',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'none',
                 minGSSize = 10,maxGSSize = 500,qvalueCutoff = 1,use_internal_data = FALSE)
barplot(kegga,showCategory=10,drop=T,x = "GeneRatio",color = "pvalue")

eea<-kegga@result
eea1<-eea[which(eea$pvalue<0.05),]
setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/7_CNV/fisher",sep=""))
write.table(eea1,"kegg_CNV_amp_pathway_P0.05.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

genec=bitr(censored_gene1[,1],fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
keggc<-enrichKEGG(genec[,2],organism = 'hsa',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'none',
                  minGSSize = 10,maxGSSize = 500,qvalueCutoff = 1,use_internal_data = FALSE)
barplot(keggc,showCategory=10,drop=T,x = "GeneRatio",color = "pvalue")

eec<-keggc@result
eec1<-eec[which(eec$pvalue<0.05),]
setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/7_CNV/fisher",sep=""))
write.table(eec1,"kegg_CNV_censored_pathway_P0.05.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
}
G="S1181"
file_name="S1181"
CNV(G,file_name)

#########################  一次性都做了，结果用表汇总
rm(list=ls())
library("data.table")
library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(stringr)
setwd("~/xjj/WGS_CNV/CNV_FACETS")
mut<-fread("CNV_FACETS_median_CNA_matrix.txt",header=T,data.table=F)
setwd("~/xjj/drug/drug_result/durg62")
tresult_want<-read.table("drug42_sample.txt",sep="\t",header=T)
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")

CNV_result<-data.frame(matrix(0,nrow(tresult_want),7))
for(k in 1:nrow(tresult_want)){
  G=as.character(tresult_want[k,1])
  cat(G,"\n")
  file_name=as.character(tresult_want[k,1])
  common_sample<-as.data.frame(tresult_want[match(G,tresult_want[,1]),])
  sensitivity<-unlist(strsplit(as.character(common_sample[1,5]), ","))
  resistance<-unlist(strsplit(as.character(common_sample[1,4]), ","))
  resistance_new<-as.character(sample_id[match(resistance,sample_id[,2]),1])
  sensitivity_new<-as.character(sample_id[match(sensitivity,sample_id[,2]),1])
  res_int<-intersect(resistance_new,colnames(mut))
  resistance_mut<-mut[,match(res_int,colnames(mut))]
  sen_int<-intersect(sensitivity_new,colnames(mut))
  sensitivity_mut<-mut[,match(sen_int,colnames(mut))]
  geneid<-mut[,1]
  data_amp<-data.frame()
  for(i in 1:length(geneid)){
    a<-sum(sensitivity_mut[i,]==1)
    c<-sum(sensitivity_mut[i,]!=1)
    b<-sum(resistance_mut[i,]==1)
    d<-sum(resistance_mut[i,]!=1)
    tmp<-matrix(c(a, c,b, d),
                nrow = 2,
                dimnames = list(Truth = c("amplification", "wild"),
                                Guess = c("sensitivity", "resistance")))
    P<-fisher.test(tmp,alternative = "two.sided",conf.level = 0.95)$p.value
    data_amp[i,1]<-geneid[i]
    data_amp[i,2]<-P
  }
  data_amp[,3]<-p.adjust(data_amp[,2],method="BH")
  amp_gene<-data_amp[which(data_amp[,2]<0.05),]  
  colnames(amp_gene)<-c("geneid","Pvalue","FDR")
  cat(paste("The number of amp gene is ",nrow(amp_gene),sep=""),"\n")
  cat(paste(amp_gene[,1],sep=""),"\n")
  CNV_result[k,1]<-G
  CNV_result[k,2]<-nrow(amp_gene)
  if(nrow(amp_gene)>0){
    CNV_result[k,3]<-paste0(amp_gene[,1],collapse =",")
  }
  if(nrow(amp_gene)==0){
    CNV_result[k,3]<-0
  }
  
  data_censored<-data.frame()
  for(j in 1:length(geneid)){
    a<-sum(sensitivity_mut[j,]==(-1))
    c<-sum(sensitivity_mut[j,]!=(-1))
    b<-sum(resistance_mut[j,]==(-1))
    d<-sum(resistance_mut[j,]!=(-1))
    tmp<-matrix(c(a, c,b, d),
                nrow = 2,
                dimnames = list(Truth = c("censored", "wild"),
                                Guess = c("sensitivity", "resistance")))
    P<-fisher.test(tmp,alternative = "two.sided",conf.level = 0.95)$p.value
    data_censored[j,1]<-geneid[i]
    data_censored[j,2]<-P
  }
  data_censored[,3]<-p.adjust(data_censored[,2],method="BH")
  censored_gene<-data_censored[which(data_censored[,2]<0.05),]  ###控制P值
  colnames(censored_gene)<-c("geneid","Pvalue","FDR")
  cat(paste("The number of censored gene is ",nrow(censored_gene),sep=""),"\n")
  cat(paste(censored_gene[,1],sep=""),"\n")
  CNV_result[k,4]<-nrow(censored_gene)
  if(nrow(censored_gene)>0){
    CNV_result[k,5]<-paste0(censored_gene[,1],collapse =",")
  }
  if(nrow(censored_gene)==0){
    CNV_result[k,5]<-0
  }
  
  int_CNV<-intersect(censored_gene[,1],amp_gene[,1])
  censored_gene1<-censored_gene[-match(int_CNV,censored_gene[,1]),]#去掉既是扩增也是删失的基因
  amp_gene1<-amp_gene[-match(int_CNV,amp_gene[,1]),]
  cat(paste("The number of amp int censored is ",length(int_CNV),sep=""),"\n")
  cat(paste("The number of final amp is ",nrow(amp_gene1),sep=""),"\n")
  cat(paste("The number of final censored is ",nrow(censored_gene1),sep=""),"\n")
  CNV_result[k,6]<-length(int_CNV)
  if(length(int_CNV)>0){
    CNV_result[k,7]<-paste0(int_CNV,collapse =",")
  }
  if(length(int_CNV)==0){
    CNV_result[k,7]<-0
  }
  
  if (nrow(amp_gene1)>0){
    setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/7_CNV/fisher",sep=""))
    write.table(amp_gene1,paste(file_name,"sensitivity-resistance-amp-CNV-gene.txt",sep="-"),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  }
  if (nrow(amp_gene1)>0){
    setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/7_CNV/fisher",sep=""))
    write.table(censored_gene1,paste(file_name,"sensitivity-resistance-censored-CNV-gene.txt",sep="-"),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  }
  
}
colnames(CNV_result)<-c("drug","Number of amp","amp_name","number of delete","delete_name","Number of int gene","int_gene_name")
setwd("~/xjj/drug/drug_result/durg62")
write.table(CNV_result,"CNV_result_42drug.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出


#######################################################################
############    RNAseq DEGs  一键式运行，结果整理成一张表
rm(list=ls())
freq<-0.01

setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
library("data.table")
exp1<-fread("star_rsem.GeneSymbol.readscount.xls",header=T,data.table=F)
exp2<-floor(exp1[-1])
delete<-apply(exp2,1,function(x) mean(x==0))
cou<-which(delete>0.5)####在超过一半样本以上的基因删去
exp<-exp2[(-cou),]
setwd("~/xjj/drug/drug_result/durg62")
tresult_want<-read.table("drug42_sample.txt",sep="\t",header=T)
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
library(DESeq2)
library(dplyr)

RNAseq_result<-data.frame(matrix(0,nrow(tresult_want),6))
for(k in 1:nrow(tresult_want)){
  G=as.character(tresult_want[k,1])
  cat(G,"\n")
  file_name=as.character(tresult_want[k,1])
  up_name<-paste("P0.01_",G,"_up",sep="")
  down_name<-paste("P0.01_",G,"_down",sep="")
  common_sample<-as.data.frame(tresult_want[match(G,tresult_want[,1]),])
  sensitivity<-unlist(strsplit(as.character(common_sample[1,5]), ","))
  resistance<-unlist(strsplit(as.character(common_sample[1,4]), ","))
  RNAseq_result[k,1]<-G
  RNAseq_result[k,2]<-as.character(common_sample[1,5])
  RNAseq_result[k,3]<-as.character(common_sample[1,4])
  resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
  sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]
  resistance_exp<-exp[,match(resistance_new,colnames(exp))]
  sensitivity_exp<-exp[,match(sensitivity_new,colnames(exp))]
  geneid<-exp1[(-cou),1]
  colDate<-data.frame(row.names = c(as.vector(resistance_new),as.vector(sensitivity_new)),
                        condition=factor(c(rep("resistance",length(resistance_new)),rep("sensitivity",length(sensitivity_new))))
  )
  datexpr<-cbind(resistance_exp,sensitivity_exp)
  counts <- apply(datexpr,2,as.numeric)###矩阵中必须是数值
  dds<-DESeqDataSetFromMatrix(countData = counts,colData = colDate,design = ~condition)
  dds<-DESeq(dds)##进行标准化分析
  res<-results(dds)##将结果输出
  res<-as.data.frame(res)
  res<-cbind(geneid,res)  ##对数据增加一列
  colnames(res)<- c('gene_id',"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj" )##对列名进行修改
  setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/1_RNAseq_DEGs",sep=""))
  write.table(res,paste("sensitivity-resistance-all-DESeq2_",G,"_RNAseq_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  resSig<-res[which(res$pvalue<freq),]
  #pvalue padj
  resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
  resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
  cat(paste("up_DEGs number is",sum(resSig$up_down=='up'),sep=" "),"\n")
  cat(paste("down_DEGs number is",sum(resSig$up_down=='down'),sep=" "),"\n")
  cat(paste("all_DEGs number is",nrow(resSig),sep=" "),"\n")
  RNAseq_result[k,4]<-sum(resSig$up_down=='up')
  RNAseq_result[k,5]<-sum(resSig$up_down=='down')
  RNAseq_result[k,6]<-nrow(resSig)
  up_gene<-resSig[which(resSig$log2FoldChange>0),]
  down_gene<-resSig[which(resSig$log2FoldChange<0),]
  setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/1_RNAseq_DEGs",sep=""))
  write.table(up_gene,paste("sensitivity-resistance-all-DESeq2_RNAseq_",up_name,"_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  write.table(down_gene,paste("sensitivity-resistance-all-DESeq2_RNAseq_",down_name,"_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
}
colnames(RNAseq_result)<-c("drug","sensitivity","resistance","up-DEGs","down-DEGs","all-DEGs")
setwd("~/xjj/drug/drug_result/durg62")
write.table(RNAseq_result,"RNAseq_result_result_42drug_P0.05.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出


#############################  每个药物单独运行 ########################
rm(list=ls())
freq<-0.05
G<-"S1030"
name<-"S1030"
up_name<-"P0.05_S1030_up"
down_name<-"P0.05_S1030_down"
file_name<-"S1030"

RNA_DEGs<-function(G,name,freq,up_name,down_name,file_name){
setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
library("data.table")
exp1<-fread("star_rsem.GeneSymbol.readscount.xls",header=T,data.table=F)
exp2<-floor(exp1[-1])
delete<-apply(exp2,1,function(x) mean(x==0))
cou<-which(delete>0.5)####在超过一半样本以上的基因删去
exp<-exp2[(-cou),]
setwd("~/xjj/drug/drug_result/durg62")
tresult_want<-read.table("drug62_sample.txt",sep="\t",header=T)
common_sample<-as.data.frame(tresult_want[match(G,tresult_want[,1]),])
library(dplyr)
sensitivity<-unlist(strsplit(as.character(common_sample[1,5]), ","))
resistance<-unlist(strsplit(as.character(common_sample[1,4]), ","))
cat(sensitivity,"\n")
cat(resistance,"\n")
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]
resistance_exp<-exp[,match(resistance_new,colnames(exp))]
sensitivity_exp<-exp[,match(sensitivity_new,colnames(exp))]
geneid<-exp1[(-cou),1]
library(DESeq2)
colDate<-data.frame(row.names = c(as.vector(resistance_new),as.vector(sensitivity_new)),
                    condition=factor(c(rep("resistance",length(resistance_new)),rep("sensitivity",length(sensitivity_new))))
)
datexpr<-cbind(resistance_exp,sensitivity_exp)
counts <- apply(datexpr,2,as.numeric)###矩阵中必须是数值
dds<-DESeqDataSetFromMatrix(countData = counts,colData = colDate,design = ~condition)
dds<-DESeq(dds)##进行标准化分析
res<-results(dds)##将结果输出
res<-as.data.frame(res)
res<-cbind(geneid,res)  ##对数据增加一列
colnames(res)<- c('gene_id',"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj" )##对列名进行修改
setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/1_RNAseq_DEGs",sep=""))
write.table(res,paste("sensitivity-resistance-all-DESeq2_",name,"_RNAseq_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
resSig<-res[which(res$pvalue<freq),]
#pvalue padj
resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
cat(paste("up_DEGs number is",sum(resSig$up_down=='up'),sep=" "),"\n")
cat(paste("down_DEGs number is",sum(resSig$up_down=='down'),sep=" "),"\n")
cat(paste("all_DEGs number is",nrow(resSig),sep=" "),"\n")
up_gene<-resSig[which(resSig$log2FoldChange>0),]
down_gene<-resSig[which(resSig$log2FoldChange<0),]
resistance_exp1 <- apply(resistance_exp,2,as.numeric)
sensitivity_exp1 <- apply(sensitivity_exp,2,as.numeric)
resistance_exp11<-apply(resistance_exp1,1,mean)
sensitivity_exp11<-apply(sensitivity_exp1,1,mean)
up1<-match(up_gene[,1],geneid)
down1<-match(down_gene[,1],geneid)
cat(paste("up-down value is",sum((sensitivity_exp11[up1]-resistance_exp11[up1])>0),sep=" "),"\n")
cat(paste("down-up value is",sum((sensitivity_exp11[down1]-resistance_exp11[down1])<0),sep=" "),"\n")
setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/1_RNAseq_DEGs",sep=""))
write.table(up_gene,paste("sensitivity-resistance-all-DESeq2_RNAseq_",up_name,"_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(down_gene,paste("sensitivity-resistance-all-DESeq2_RNAseq_",down_name,"_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
}
RNA_DEGs(G,name,freq,up_name,down_name,file_name)

#########  RNAseq DEGs Volcano Plot
rm(list=ls())
library(ggrepel) 
library(ggplot2)
file_name<-"S1030"
name<-"S1030"
freq<-0.05 # pvalue
FC<-0
text_freq<-0.01  #FDR
text_FC<-2
RNA_Volcano<-function(file_name,name,freq,FC,text_freq,text_FC){
setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/1_RNAseq_DEGs",sep=""))
df<-read.table(paste("sensitivity-resistance-all-DESeq2_",name,"_RNAseq_DEGs.txt",sep=""),sep = '\t',header= T,row.names = 1)
df$threshold = factor(ifelse(df$pvalue < freq & abs(df$log2FoldChange) >= FC, ifelse(df$log2FoldChange >= FC ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
df$gene <- row.names(df) #添加一列基因名，以便备注
p<-ggplot(df,aes(x=log2FoldChange,y= -log10(pvalue),color=threshold))+
  geom_point(data = df[df$pvalue<freq & abs(df$log2FoldChange)>FC,],size = 1)+ 
  geom_point(data = df[df$pvalue>freq | abs(df$log2FoldChange)<FC,],size = 1)+
  scale_color_manual(values=c('blue','grey','red'))+#确定点的颜色
  geom_text_repel(
    data = df[df$padj<text_freq & abs(df$log2FoldChange)>text_FC,],
    aes(label = gene),
    size = 3,
    color = "black",
    segment.color = "black", show.legend = FALSE )+#添加关注的点的基因名
  ylab('-log10 (pvalue)')+#修改y轴名称
  xlab('log2 (FoldChange)')+#修改x轴名称
  geom_hline(yintercept=-log10(freq),linetype=(text_FC*2))#添加横线|logFoldChange|>0.25
  geom_vline(xintercept=c(-(text_FC),text_FC),linetype=4)#添加竖线padj<0.05
plot(p)
DEGs<-df[which(df$padj<text_freq & abs(df$log2FoldChange)>text_FC),]
setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/1_RNAseq_DEGs/Volcano",sep=""))
write.table(DEGs,paste("FDR_",text_freq,"_FC_",text_FC,"_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
}
RNA_Volcano(file_name,name,freq,FC,text_freq,text_FC)
#################################################################################
################  DEGs与组蛋白修饰酶的交叠
rm(list=ls())
library(gplots)
library(VennDiagram)
file_name<-"S1030"
name<-"S1030"
freq<-0.05 #差异基因的阈值
FC<-0
freq_name<-"P0.05"
DEGs_Acetylase<-function(file_name,name,freq,FC,freq_name){
setwd("~/xjj/drug/drug_result/HDAC_frontiers/5_histone-modifying enzymes")
Acetylase<-read.table("histone-modifying-enzymes.txt",header=T,sep="\t")
setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/1_RNAseq_DEGs",sep=""))
df<-read.table(paste("sensitivity-resistance-all-DESeq2_",name,"_RNAseq_DEGs.txt",sep=""),sep = '\t',header= T,row.names = 1)
DEGs<-df[df$pvalue<freq,]
DEGs[which(DEGs$log2FoldChange>FC),'up_down']<-'up'
DEGs[which(DEGs$log2FoldChange<FC),'up_down']<-'down'
#pvalue padj
Acetylase_DEGs<-intersect(toupper(Acetylase[,1]),toupper(rownames(DEGs)))
Acetylase_info1<-Acetylase[match(Acetylase_DEGs,toupper(Acetylase[,1])),]
DEGs$gene<-toupper(rownames(DEGs))
Acetylase_info2<-DEGs[match(Acetylase_DEGs,toupper(rownames(DEGs))),7:8]
Acetylase_info<-cbind(Acetylase_info1,Acetylase_info2)
Acetylase_name<- unique(toupper(Acetylase[,1]))
DEGs_name <- unique(toupper(rownames(DEGs)))
input  <-list(Acetylase_name,DEGs_name)
ven<-venn(input,showSetLogicLabel=TRUE)
plot(ven)
tmp <- venn(input)
int<-attr(tmp, "intersections")
setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/1_RNAseq_DEGs/Acetylase",sep=""))
write.table(Acetylase_info,paste("Acetylase_info_",freq_name,".txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
venn.diagram(x=list(Acetylase=Acetylase_name,DEGs=DEGs_name),cex = 1,margin = 0.1, paste("Acetylase_info_",freq_name,".png",sep=""),fill=c("red","blue"))
}
DEGs_Acetylase(file_name,name,freq,FC,freq_name)

#################################################################################
################  clusterProfiler
rm(list=ls())
library(org.Hs.eg.db)
library(clusterProfiler)
library(stringr)

file_name<-"S1030"
name<-"P0.05_S1030_down"
kegg_freq_name<-"FDR0.05_down"
go_freq_name<-"FDR0.01_down"
keggfreq<-0.05 #kegg p<0.05
gofreq<-0.01    #go FDR<0.05
clusterProfiler<-function(file_name,name,keggfreq,kegg_freq_name,gofreq,go_freq_name){
  setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/1_RNAseq_DEGs",sep=""))
  up_gene<-read.table(paste("sensitivity-resistance-all-DESeq2_RNAseq_",name,"_DEGs.txt",sep=""),sep = '\t',header= T)
  gene=bitr(up_gene[,1],fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
  kegg<-enrichKEGG(gene[,2],organism = 'hsa',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'BH',
                 minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.05,use_internal_data = FALSE)
  ke<-barplot(kegg,showCategory=10,drop=T,x = "GeneRatio",color = "p.adjust")
  plot(ke)
  #em<-emapplot(kegg,showCategory = 30)
  #plot(em)
  ee<-kegg@result
  ee1<-ee[which(ee$p.adjust<keggfreq),]
  cat(paste("The pathway of KEGG is",nrow(ee1),sep=" "),"\n")
  #p.adjust
  enrich_gene<-ee1$geneID
  pathway_gene<-unique(unlist(strsplit(enrich_gene,split="/")))
  pathway_gene2=bitr(pathway_gene,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
  setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/1_RNAseq_DEGs/pathway",sep=""))
  write.table(ee1,paste("kegg_DEG_pathway_",kegg_freq_name,".txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  write.table(pathway_gene2,paste("kegg_DEG_pathway_",kegg_freq_name,"_gene.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  go<-enrichGO(gene[,2],OrgDb = org.Hs.eg.db,ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.01,keyType = 'ENTREZID')
  godot<-dotplot(go,showCategory=10)
  plot(godot)
  a<-go@result
  go_BP<-a[which(a$p.adjust<gofreq),]
  cat(paste("The pathway of GO is",nrow(go_BP),sep=" "),"\n")
  enrich_genego<-go_BP$geneID
  pathway_genego<-unique(unlist(strsplit(enrich_genego,split="/")))
  pathway_genego2=bitr(pathway_genego,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
  write.table(pathway_genego2,paste("GO_DEG_pathway_",go_freq_name,"_gene.txt",sep=""),sep = '\t',col.names = F,row.names = F,quote = FALSE)##数据输出
  write.table(go_BP,paste("GO_DEG_pathway_",go_freq_name,".txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
}
clusterProfiler(file_name,name,keggfreq,kegg_freq_name,gofreq,go_freq_name)
#######################    ATACseq   ################################################
#############  DApeaks
rm(list=ls())
G<-"S7575"
file_name<-"S7575"
freq<-0.05
out_file<-"sen-res-S7575-0.05-P"
DApeaks<-function(G,file_name,freq,out_file){
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_count_name<-read.table("peak_count.txt",sep="\t",header=T)
setwd("~/xjj/drug/drug_result/durg62")
tresult_want<-read.table("drug62_sample.txt",sep="\t",header=T)
common_sample<-as.data.frame(tresult_want[match(G,tresult_want[,1]),])
library(dplyr)
sensitivity<-unlist(strsplit(as.character(common_sample[1,5]), ","))
resistance<-unlist(strsplit(as.character(common_sample[1,4]), ","))
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]
gsub("\\.","-",colnames(peak_count_name))
resistance_peak<-peak_count_name[,match(resistance_new,gsub("\\.","-",colnames(peak_count_name)))]
sensitivity_peak<-peak_count_name[,match(sensitivity_new,gsub("\\.","-",colnames(peak_count_name)))]
peak_name<-peak_count_name[,1]
resistance_peak1 <- apply(resistance_peak,2,as.numeric)
sensitivity_peak1 <- apply(sensitivity_peak,2,as.numeric)
library(DESeq2)
colDate<-data.frame(row.names = c(as.vector(resistance),as.vector(sensitivity)),
                    condition=factor(c(rep("resistance",length(resistance)),rep("sensitivity",length(sensitivity))))
)
datexpr<-cbind(resistance_peak,sensitivity_peak)##前面的相对于后面的上下调
counts <- apply(datexpr,2,as.numeric)   ###矩阵中必须是数值
dds<-DESeqDataSetFromMatrix(countData = counts,colData = colDate,design = ~condition)
dds<-DESeq(dds)##进行标准化分析
sizeFactors(dds)##查看每个主成分的标准化值
res<-results(dds)##将结果输出
res<-as.data.frame(res)
res<-cbind(peak_count_name[,1],res)  ##对数据增加一列
colnames(res)<- c('peak_id',"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj" )##对列名进行修改
setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/2_ATACseq_DApeaks",sep=""))
write.table(res,paste("sensitivity-resistance-all-DESeq2_",file_name,"_ATAC.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
resSig<-res[which(res$pvalue<freq),]
# pvalue padj
resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
cat(paste("up DApeaks is",sum(resSig$up_down=='up'),sep=" "),"\n")
cat(paste("down DApeaks is",sum(resSig$up_down=='down'),sep=" "),"\n")
cat(paste("all DApeaks is",nrow(resSig),sep=" "),"\n")
up_gene<-as.data.frame(resSig[which(resSig$log2FoldChange>0),1])
down_gene<-as.data.frame(resSig[which(resSig$log2FoldChange<0),1])
up1<-match(up_gene[,1],peak_count_name[,1])
down1<-match(down_gene[,1],peak_count_name[,1])
library(dplyr)
a_up<-apply(up_gene,1,function(x) unlist(strsplit(as.character(x), "_")))
d_up<-t(a_up)
e_up<-data.frame(rep("+",nrow(d_up)))
peaks_up<-cbind(up_gene,d_up,e_up)
colnames(peaks_up)<-c("peak_id","chr","start","end","strand")
setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/2_ATACseq_DApeaks",sep=""))
write.table(peaks_up,file=paste(out_file,"-logFC0-DESeq2-up.txt",sep=""),sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')
write.table(d_up,file=paste(out_file,"-logFC0-DESeq2-up.gff",sep=""),sep = '\t',
            col.names = F,row.names = F,quote = FALSE)
a_down<-apply(down_gene,1,function(x) unlist(strsplit(as.character(x), "_")))
d_down<-t(a_down)
e_down<-data.frame(rep("+",nrow(d_down)))
peaks_down<-cbind(down_gene,d_down,e_down)
colnames(peaks_down)<-c("peak_id","chr","start","end","strand")
write.table(peaks_down,file=paste(out_file,"-logFC0-DESeq2-down.txt",sep=""),sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')
write.table(d_down,file=paste(out_file,"-logFC0-DESeq2-down.gff",sep=""),sep = '\t',
            col.names = F,row.names = F,quote = FALSE)
#####  bed
f_up<-data.frame(rep(1,nrow(d_up)))
peaks_bed_up<-cbind(d_up,up_gene,f_up,e_up)
colnames(peaks_bed_up)<-c("chr","start","end","peak_id","score","strand")
write.table(peaks_bed_up,file=paste(out_file,"-logFC0-DESeq2-up.bed",sep=""),sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')
f_down<-data.frame(rep(1,nrow(d_down)))
peaks_bed_down<-cbind(d_down,down_gene,f_down,e_down)
colnames(peaks_bed_down)<-c("chr","start","end","peak_id","score","strand")
write.table(peaks_bed_down,file=paste(out_file,"-logFC0-DESeq2-down.bed",sep=""),sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')
}
DApeaks(G,file_name,freq,out_file)


#########################  ATACseq 一体式   ###############################
#############  DApeaks
rm(list=ls())
freq<-0.05
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_count_name<-read.table("peak_count.txt",sep="\t",header=T)
setwd("~/xjj/drug/drug_result/durg62")
tresult_want<-read.table("drug42_sample.txt",sep="\t",header=T)
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
library(DESeq2)
library(dplyr)


ATACseq_result<-data.frame(matrix(0,nrow(tresult_want),6))
for(k in 33:nrow(tresult_want)){
  G=as.character(tresult_want[k,1])
  cat(G,"\n")
  file_name=as.character(tresult_want[k,1])
  out_file<-paste("sen-res-",G,"-0.05-P",sep="")
  
  common_sample<-as.data.frame(tresult_want[match(G,tresult_want[,1]),])
  sensitivity<-unlist(strsplit(as.character(common_sample[1,5]), ","))
  resistance<-unlist(strsplit(as.character(common_sample[1,4]), ","))
  ATACseq_result[k,1]<-G
  ATACseq_result[k,2]<-as.character(common_sample[1,5])
  ATACseq_result[k,3]<-as.character(common_sample[1,4])
  resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
  sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]
  gsub("\\.","-",colnames(peak_count_name))
  resistance_peak<-peak_count_name[,match(resistance_new,gsub("\\.","-",colnames(peak_count_name)))]
  sensitivity_peak<-peak_count_name[,match(sensitivity_new,gsub("\\.","-",colnames(peak_count_name)))]
  peak_name<-peak_count_name[,1]
  colDate<-data.frame(row.names = c(as.vector(resistance),as.vector(sensitivity)),
                      condition=factor(c(rep("resistance",length(resistance)),rep("sensitivity",length(sensitivity))))
  )
  datexpr<-cbind(resistance_peak,sensitivity_peak)##前面的相对于后面的上下调
  counts <- apply(datexpr,2,as.numeric)   ###矩阵中必须是数值
  dds<-DESeqDataSetFromMatrix(countData = counts,colData = colDate,design = ~condition)
  dds<-DESeq(dds)##进行标准化分析
  sizeFactors(dds)##查看每个主成分的标准化值
  res<-results(dds)##将结果输出
  res<-as.data.frame(res)
  res<-cbind(peak_count_name[,1],res)  ##对数据增加一列
  colnames(res)<- c('peak_id',"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj" )##对列名进行修改
  setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/2_ATACseq_DApeaks",sep=""))
  write.table(res,paste("sensitivity-resistance-all-DESeq2_",file_name,"_ATAC.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  resSig<-res[which(res$pvalue<freq),]
  # pvalue padj
  resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
  resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
  cat(paste("up DApeaks is",sum(resSig$up_down=='up'),sep=" "),"\n")
  cat(paste("down DApeaks is",sum(resSig$up_down=='down'),sep=" "),"\n")
  cat(paste("all DApeaks is",nrow(resSig),sep=" "),"\n")
  ATACseq_result[k,4]<-sum(resSig$up_down=='up')
  ATACseq_result[k,5]<-sum(resSig$up_down=='down')
  ATACseq_result[k,6]<-nrow(resSig)
  up_gene<-as.data.frame(resSig[which(resSig$log2FoldChange>0),1])
  down_gene<-as.data.frame(resSig[which(resSig$log2FoldChange<0),1])
  up1<-match(up_gene[,1],peak_count_name[,1])
  down1<-match(down_gene[,1],peak_count_name[,1])
  a_up<-apply(up_gene,1,function(x) unlist(strsplit(as.character(x), "_")))
  d_up<-t(a_up)
  e_up<-data.frame(rep("+",nrow(d_up)))
  peaks_up<-cbind(up_gene,d_up,e_up)
  colnames(peaks_up)<-c("peak_id","chr","start","end","strand")
  setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/2_ATACseq_DApeaks",sep=""))
  write.table(peaks_up,file=paste(out_file,"-logFC0-DESeq2-up.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  write.table(d_up,file=paste(out_file,"-logFC0-DESeq2-up.gff",sep=""),sep = '\t',
              col.names = F,row.names = F,quote = FALSE)
  a_down<-apply(down_gene,1,function(x) unlist(strsplit(as.character(x), "_")))
  d_down<-t(a_down)
  e_down<-data.frame(rep("+",nrow(d_down)))
  peaks_down<-cbind(down_gene,d_down,e_down)
  colnames(peaks_down)<-c("peak_id","chr","start","end","strand")
  write.table(peaks_down,file=paste(out_file,"-logFC0-DESeq2-down.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  write.table(d_down,file=paste(out_file,"-logFC0-DESeq2-down.gff",sep=""),sep = '\t',
              col.names = F,row.names = F,quote = FALSE)
  #####  bed
  f_up<-data.frame(rep(1,nrow(d_up)))
  peaks_bed_up<-cbind(d_up,up_gene,f_up,e_up)
  colnames(peaks_bed_up)<-c("chr","start","end","peak_id","score","strand")
  write.table(peaks_bed_up,file=paste(out_file,"-logFC0-DESeq2-up.bed",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  f_down<-data.frame(rep(1,nrow(d_down)))
  peaks_bed_down<-cbind(d_down,down_gene,f_down,e_down)
  colnames(peaks_bed_down)<-c("chr","start","end","peak_id","score","strand")
  write.table(peaks_bed_down,file=paste(out_file,"-logFC0-DESeq2-down.bed",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
}
colnames(ATACseq_result)<-c("drug","sensitivity","resistance","up-peaks","down-peaks","all-peaks")
setwd("~/xjj/drug/drug_result/durg62")
write.table(ATACseq_result,"ATACseq_result_result_42drug",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出


#######################################################################################
####  homer注释完的信息,画饼图
rm(list=ls())
anno_name<-"PAC-0.05-P-up"
file_name<-"PAC"

ATAC_pie<-function(file_name,anno_name){
setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/2_ATACseq_DApeaks",sep=""))
library(data.table)
homer_anno<-fread(paste("sen-res-",anno_name,"-annotation.txt",sep=""),header=T,data.table=F)
cat(length(unique(homer_anno$`Gene Name`)))
colnames(homer_anno)<-c("peakid",colnames(homer_anno)[-1])
library(stringr)
Annotation<-apply(as.data.frame(homer_anno$Annotation),1,function(x) str_replace_all(as.character(x),"\\((.*?)\\)",""))##取出括号前的字符
Annotation1<-as.data.frame(Annotation)
Annotation2<-data.frame(lapply(strsplit(as.character(Annotation1[,1]),'\\.'), function(x) x[1])%>%unlist())
position<-as.character(unique(Annotation2[,1]))
result<-data.frame()
for(i in 1:length(position)){
  result[i,1]<-position[i]
  result[i,2]<-sum(Annotation %in% position[i])
  result[i,3]<-(sum(Annotation %in% position[i])/(length(Annotation)))*100
}
library(ggplot2) # 加载包
dt = data.frame(A = result[,3],B = result[,1])#建立数据框
dt = dt[order(dt$A, decreasing = TRUE),]   ## 用 order() 让数据框的数据按 A 列数据从大到小排序
myLabel = as.vector(dt$B)   ## 转成向量，否则图例的标签可能与实际顺序不一致
myLabel = paste(myLabel, "(", round(dt$A / 1, 4), "%)", sep = "")   ## 用 round() 对结果保留两位小数
p = ggplot(dt, aes(x = "", y = A, fill = B)) + #创建坐标轴
  geom_bar(stat = "identity") + 
  geom_bar(stat = "identity", width = 1) +   #当width >= 1 时中心的杂点将消失
  coord_polar(theta = "y") +  # 把柱状图折叠成饼图（极坐标）
  labs(x = "", y = "", title = "") +  # 将横纵坐标的标签设为空
  theme(axis.ticks = element_blank()) +  # 将左上角边框的刻度去掉
  theme(legend.title = element_blank(), legend.position = "left")+   ## 将图例标题设为空，并把图例方放在左边位置
  scale_fill_discrete(breaks = dt$B, labels = myLabel)+   # 将原来的图例标签换成现在的myLabel
  theme(axis.text.x = element_blank())   ## 去掉饼图的外框上的数值，即去除原柱状图的X轴，把X轴的刻度文字去掉
  #geom_text(aes(y = A/2 + c(0, cumsum(A)[-length(A)]), x = sum(A)/20, label = myLabel), size = 5)   # 在图中加上百分比：x 调节标签到圆心的距离, y 调节标签的左右位置
print(p) #显示饼图
}
ATAC_pie(file_name,anno_name)

#####################  符合在TSS100KB之内条件的注释基因
rm(list=ls())
library(gplots)
library(VennDiagram)
file_name<-"PAC"
anno_name_down<-"PAC-0.05-P-down" #注释文件
anno_name_up<-"PAC-0.05-P-up"
up_name<-"P0.05_PAC_up"  #差异基因文件
down_name<-"P0.05_PAC_down"
freq_name<-"P0.05" #保存文件名中的阈值

DEGs_int_DApeaks<-function(file_name,anno_name_down,anno_name_up,up_name,down_name,freq_name){
setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/2_ATACseq_DApeaks",sep=""))
library(data.table)
homer_anno_down<-fread(paste("sen-res-",anno_name_down,"-annotation.txt",sep=""),header=T,data.table=F)
colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]
homer_anno_up<-fread(paste("sen-res-",anno_name_up,"-annotation.txt",sep=""),header=T,data.table=F)
colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/1_RNAseq_DEGs",sep=""))
up_gene<-read.table(paste("sensitivity-resistance-all-DESeq2_RNAseq_",up_name,"_DEGs.txt",sep=""),header=T,sep="\t")
down_gene<-read.table(paste("sensitivity-resistance-all-DESeq2_RNAseq_",down_name,"_DEGs.txt",sep=""),header=T,sep="\t")
all<-rbind(up_gene,down_gene)
DEGs<-as.character(all[,1])
down_DEGs_DApeaks<-intersect(Anno_gene_100Kb_down$`Gene Name`,DEGs)
cat(paste("Down peak int DEGs is",length(down_DEGs_DApeaks),sep=" "),"\n")
up_DEGs_DApeaks<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)
cat(paste("Up peak int DEGs is",length(up_DEGs_DApeaks),sep=" "),"\n")
down_peaks_gene <- Anno_gene_100Kb_down$`Gene Name`
up_peaks_gene <- Anno_gene_100Kb_up$`Gene Name`
input  <-list(unique(down_peaks_gene),unique(up_peaks_gene),unique(DEGs))
plot(venn(input,showSetLogicLabel=TRUE))
tmp <- venn(input)
int<-attr(tmp, "intersections")
setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/3_intersect/venn",sep=""))
venn.diagram(x=list(downpeaks=unique(down_peaks_gene),uppeaks=unique(up_peaks_gene),RNAseq=unique(DEGs)), paste("DEGs_int_DApeaks_",file_name,"_",freq_name,".png",sep=""),fill=c("red","green","blue"),margin = 0.1)
#write.table(a,paste("DApeaks_DEGs_",file_name,"_",freq_name,".txt",sep=""),col.names=F,row.names=F)
up1<-int[["B:C"]]
up2<-int[["A:B:C"]]
up3<-c(up1,up2)
int_up<-intersect(up3,up_gene$gene_id)
cat(paste("int up DEGs is",length(int_up),sep=" "),"\n")
down1<-int[["A:C"]]
down2<-int[["A:B:C"]]
down3<-c(down1,down2)
int_down<-intersect(down3,down_gene$gene_id)
cat(paste("int down DEGs is",length(int_down),sep=" "),"\n")
}
DEGs_int_DApeaks(file_name,anno_name_down,anno_name_up,up_name,down_name,freq_name)
  


######## 差异上调的peaks与DEGs对应的基因在两类样本中的mRNA水平
###真的是上调的基因
up1<-int[["B:C"]]
up2<-int[["A:B:C"]]
up3<-c(up1,up2)

down1<-int[["A:C"]]
down2<-int[["A:B:C"]]
down3<-c(down1,down2)
up3<-down3
setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
library("data.table")
exp4<-fread("star_rsem.GeneSymbol.readscount.xls",header=T,data.table=F)
exp1<-exp4[match(up3,exp4[,1]),]
exp2<-floor(exp1[-1])
delete<-apply(exp2,1,function(x) mean(x==0))
cou<-which(delete>0.5)####在超过一半样本以上的基因删去
exp<-exp2

sensitivity<-c("PC.64","PC.81","PC.116","PC.L","PC.2","PC.115","PC.16","PC.97","PC.101","PC.27","PC.8","PC.109","PC.136","PC.13","PC.78","PC.98","PC.130","PC.117","PC.139","PC.52")  ##20
resistance<-c("PC.G","PC.111","PC.134","PC.14","PC.40","PC.135","PC.22","PC.18","PC.105","PC.104","PC.56","PC.119","PC.112","PC.102","PC.121","PC.5","PC.I")  ##17
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]
resistance_exp1<-exp[,match(resistance_new,colnames(exp))]
sensitivity_exp1<-exp[,match(sensitivity_new,colnames(exp))]
resistance_exp <- apply(resistance_exp1,2,as.numeric)
sensitivity_exp <- apply(sensitivity_exp1,2,as.numeric)
geneid<-exp1[,1]
### T.test
tresult=matrix(0,nrow(resistance_exp),4)
for (i in 1:nrow(resistance_exp)){
  ttest=t.test(sensitivity_exp[i,],resistance_exp[i,],alternative = "greater")
  tresult[i,1:3]=c(geneid[i],ttest$statistic,ttest$p.value)
}
tresult[,4]=p.adjust(tresult[,3],method="BH")
colnames(tresult)<-c("geneid","statistic","p.value","FDR")
want_result_up<-tresult[which(tresult[,3]<0.05 & tresult[,2]>0),]

tresult=matrix(0,nrow(resistance_exp),4)
for (i in 1:nrow(resistance_exp)){
  ttest=t.test(sensitivity_exp[i,],resistance_exp[i,],alternative = "less")
  tresult[i,1:3]=c(geneid[i],ttest$statistic,ttest$p.value)
}
tresult[,4]=p.adjust(tresult[,3],method="BH")
colnames(tresult)<-c("geneid","statistic","p.value","FDR")
want_result_down<-tresult[which(tresult[,3]<0.05 & tresult[,2]<0),]

##### wilcox.test
wresult=matrix(0,nrow(resistance_exp),4)
for (i in 1:nrow(resistance_exp)){
  wtest=wilcox.test(sensitivity_exp[i,],resistance_exp[i,],alternative = "greater")
  wresult[i,1:3]=c(geneid[i],wtest$statistic,wtest$p.value)
}
wresult[,4]=p.adjust(wresult[,3],method="BH")
colnames(wresult)<-c("geneid","statistic","p.value","FDR")
want_result_upw<-wresult[which(wresult[,3]<0.05 & wresult[,2]>0),]

wresult=matrix(0,nrow(resistance_exp),4)
for (i in 1:nrow(resistance_exp)){
  wtest=wilcox.test(sensitivity_exp[i,],resistance_exp[i,],alternative = "two.sided")
  wresult[i,1:3]=c(geneid[i],wtest$statistic,wtest$p.value)
}
wresult[,4]=p.adjust(wresult[,3],method="BH")
colnames(wresult)<-c("geneid","statistic","p.value","FDR")
want_result_doww<-wresult[which(wresult[,3]<0.05 & wresult[,2]<0),]

####  DESeq2
up1<-int[["B:C"]]
up2<-int[["A:B:C"]]
up3<-c(up1,up2)
int_up<-intersect(up3,up_gene$gene_id)

down1<-int[["A:C"]]
down2<-int[["A:B:C"]]
down3<-c(down1,down2)
int_down<-intersect(down3,down_gene$gene_id)

##########################################################
################  ATACseq-CNV
setwd("~/xjj/drug/drug_result/HDAC_frontiers/7_CNV")
amp_CNV<-read.table("sensitivity-resistance-amp-CNV-gene.txt",sep = '\t',header=T)##数据输出
censored_CNV<-read.table("sensitivity-resistance-censored-CNV-gene.txt",sep = '\t',header=T)##数据输出
setwd("~/xjj/WGS_CNV/CNV_FACETS")
gene_position<-read.table("gene_position.txt",sep = '\t',header=T)##数据输出
amp_position<-gene_position[match(amp_CNV[,1],gene_position[,5]),]
censored_position<-gene_position[match(censored_CNV[,1],gene_position[,5]),]
setwd("~/xjj/drug/drug_result/HDAC_frontiers/2_ATACseq_DApeaks")
library(data.table)
homer_anno_down<-fread("sen-res-HDAC13-heatmap-0.05-P-down-annotation.txt",header=T,data.table=F)
colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]

homer_anno_up<-fread("sen-res-HDAC13-heatmap-0.05-P-up-annotation.txt",header=T,data.table=F)
colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
amp_CNV_DApeak<-intersect(amp_position[,5],Anno_gene_100Kb_up$`Gene Name`)
amp_CNV_DApeak1<-intersect(amp_position[,5],Anno_gene_100Kb_down$`Gene Name`)

censored_CNV_DApeak1<-intersect(censored_position[,5],Anno_gene_100Kb_up$`Gene Name`)

#cn<-as.character(amp_position[match(amp_CNV_DApeak,amp_position[,5]),4])
#library(dplyr)
#a_up<-data.frame(lapply(strsplit(as.character(cn),'='), function(x) x[2])%>%unlist())


##################################################################################
##############  DEpeaks在不同样本变化的倍数与DEGs在不同样本变化的倍数之间的相关性
rm(list=ls())
library(data.table)
file_name<-"PAC"
anno_name_down<-"PAC-0.05-P-down" #注释文件
anno_name_up<-"PAC-0.05-P-up"
up_name<-"P0.05_PAC_up"  #差异基因文件
down_name<-"P0.05_PAC_down"
freq_name<-"P0.05"  #保存文件名时的阈值
method<-"spearman"
DEGs_corre_DApeaks<-function(file_name,anno_name_down,anno_name_up,up_name,down_name,method){
setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/2_ATACseq_DApeaks",sep=""))
homer_anno_down<-fread(paste("sen-res-",anno_name_down,"-annotation.txt",sep=""),header=T,data.table=F)
colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]
homer_anno_up<-fread(paste("sen-res-",anno_name_up,"-annotation.txt",sep=""),header=T,data.table=F)
colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/1_RNAseq_DEGs",sep=""))
up_gene<-read.table(paste("sensitivity-resistance-all-DESeq2_RNAseq_",up_name,"_DEGs.txt",sep=""),header=T,sep="\t")
down_gene<-read.table(paste("sensitivity-resistance-all-DESeq2_RNAseq_",down_name,"_DEGs.txt",sep=""),header=T,sep="\t")
all<-rbind(down_gene,up_gene)
DEGs<-as.character(all[,1])
down_DEGs_DApeaks<-intersect(Anno_gene_100Kb_down$`Gene Name`,DEGs)
up_DEGs_DApeaks<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)
int_all_DEGs<-unique(c(down_DEGs_DApeaks,up_DEGs_DApeaks))
all_peaks_anno<-rbind(homer_anno_down,homer_anno_up)
#DEGs_peaks<-all_peaks_anno[all_peaks_anno$`Gene Name` %in% int_all_DEGs,]
DEGs_peaks<-all_peaks_anno[match(int_all_DEGs,all_peaks_anno$`Gene Name`),1]
setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/2_ATACseq_DApeaks",sep=""))
DApeak<-fread(paste("sensitivity-resistance-all-DESeq2_",file_name,"_ATAC.txt",sep=""),header=T,data.table=F)
DEGs_peaks_FC<-DApeak[match(DEGs_peaks,DApeak$peak_id),3]
DEGs_FC<-all[match(int_all_DEGs,DEGs),3]
gene<-as.character(all[match(int_all_DEGs,DEGs),1])
#pearson", "kendall", "spearman
cor_result<-cor.test(DEGs_peaks_FC, DEGs_FC,alternative = "two.sided",method = method)
cat(paste("correlation is",cor_result$estimate,sep=" "),"\n")
cat(paste("pvalue is",cor_result$p.value,sep=" "),"\n")
################################   画相关性分析图
a1<-matrix(DEGs_FC,ncol=1)
a2<-matrix(DEGs_peaks_FC,ncol=1)
dat<-cbind(a1,a2)
dat1<-as.data.frame(dat)
rownames(dat1)<-gene
colnames(dat1)<-c("DEGslog2FC","DApeakslog2FC")
library(ggplot2)
library(ggpubr)
library(ggrepel)
a11<-c(1:10)
a21<-c((nrow(dat1)-9):nrow(dat1))
a3<-c(a11,a21)
want_dat<-dat1[order(dat1[,1])[a3],]
p<-ggplot(data=dat1, aes(x=DEGslog2FC, y=DApeakslog2FC))+geom_point(color="red",size=0.1)+stat_smooth(method="lm",se=FALSE)+stat_cor(data=dat1, method = method)+
       ggtitle(paste(method,"-P0.05",sep="")) +
       theme(plot.title = element_text(hjust = 0.5))+
       geom_text_repel(
       data = want_dat[,c(1:2)],
       aes(label = rownames(want_dat)),
       size = 3,
       color = "black",
       segment.color = "black", show.legend = FALSE )
plot(p)
}
#pearson", "kendall", "spearman
DEGs_corre_DApeaks(file_name,anno_name_down,anno_name_up,up_name,down_name,method)
  

#############################  对高低可及peaks-DEGs进行KEGG通路富集分析
rm(list=ls())
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(stringr)

file_name<-"PAC"
anno_name_down<-"PAC-0.05-P-down" #注释文件
down_name<-"P0.05_PAC_down"
freq_name<-"FDR0.05_down"

peaks_DEGs_KEGG<-function(file_name,anno_name_down,down_name,freq_name){
setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/2_ATACseq_DApeaks",sep=""))
homer_anno_down<-fread(paste("sen-res-",anno_name_down,"-annotation.txt",sep=""),header=T,data.table=F)
colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]
setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/1_RNAseq_DEGs",sep=""))
down_gene<-read.table(paste("sensitivity-resistance-all-DESeq2_RNAseq_",down_name,"_DEGs.txt",sep=""),header=T,sep="\t")
DEGs<-as.character(down_gene[,1])
down_DEGs_DApeaks<-intersect(Anno_gene_100Kb_down$`Gene Name`,DEGs)
cat(paste("Int gene is ",length(down_DEGs_DApeaks),sep=""),"\n")
gene=bitr(down_DEGs_DApeaks,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
kegg<-enrichKEGG(gene[,2],organism = 'hsa',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'BH',
                 minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.05,use_internal_data = FALSE)
ee<-kegg@result
ee1<-ee[which(ee$p.adjust<0.05),]
cat(paste("KEGG pathway is",nrow(ee1),sep=" "),"\n")
#ee1<-ee[which(ee$p.adjust<0.05),]
ke<-barplot(kegg,showCategory=10,drop=T,x = "GeneRatio",color = "pvalue")
plot(ke)
ee2<-ee1[1:10,]
enrich_gene<-ee1$geneID
pathway_gene<-unique(unlist(strsplit(enrich_gene,split="/")))
cat(paste("Top pathway gene is ",length(pathway_gene),sep=""),"\n")
pathway_gene2=bitr(pathway_gene,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/3_intersect/pathway",sep=""))
write.table(ee1,paste("kegg_pathway_DEGs_DApeaks_",freq_name,".txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(pathway_gene,paste("kegg_pathway_DEGs_DApeaks_",freq_name,"_gene.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
}
peaks_DEGs_KEGG(file_name,anno_name_down,down_name,freq_name)

#############################################################################
######################  使用homer进行motif富集，并对应到TFs
rm(list=ls())
TF_names<-"sen-res-PAC-0.05-P-logFC0-DESeq2-up"
file_name<-"PAC"
freq<-0.05 #TF的选择，FDR<0.05 
freq_name<-"up_FDR0.05"  
a<-10  #最显著的前a个TFs 

TFs<-function(file_name,TF_names,freq,freq_name,a){
setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/4_TFs/homer_result/",TF_names,sep=""))
knownResults<-read.delim("knownResults.txt",header=T,sep='\t')
FDR05<-knownResults[knownResults$q.value..Benjamini.<freq,]
#FDR05<-knownResults[knownResults$P.value<0.01,]
library(dplyr)
homer_TF1<-data.frame(lapply(strsplit(as.character(FDR05[,1]),'/'), function(x) x[3])%>%unlist())
FDR051<-FDR05[which(homer_TF1=="Homer"),]
homer_TF2<-data.frame(lapply(strsplit(as.character(FDR051[,1]),'/'), function(x) x[1])%>%unlist())
library(stringr)
homer_TF3<-apply(homer_TF2,1,function(x) str_replace_all(as.character(x),"\\((.*?)\\)",""))##取出括号前的字符
cat(paste("The number of motif is",length(homer_TF3),sep=" "),"\n")
homer_TF4<-unique(homer_TF3)
cat(paste("The number of TFs is",length(homer_TF4),sep=" "),"\n")
new_homer_result<-cbind(homer_TF3,FDR051[,-1])
colnames(new_homer_result)<-c("TF_name",colnames(new_homer_result)[-1])
write.table(new_homer_result,paste("DApeaks_P0.05-logFC0-DESeq2-",freq_name,"_TFs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
#####  画柱状图
data1<-cbind(as.data.frame(new_homer_result[,1]),new_homer_result$P.value)
data1[,2]<-(-log2(data1[,2]))
colnames(data1)<-c("TFname","Pvalue")
data11<-data1[match(unique(data1[,1]),data1[,1]),]##重叠的TFqu取第一个
data2<-data11[1:a,]
data2$TFname=factor(data2$TFname,levels = data2[,1])
new_homer_result1<-new_homer_result[match(unique(data1[,1]),data1[,1]),]
val<-new_homer_result1$X..of.Target.Sequences.with.Motif[1:a]
library(ggplot2)
p<-ggplot(data=data2,mapping=aes(x=TFname,y=as.numeric(Pvalue)))+
  geom_bar(stat="identity",fill="blue4")+
  #geom_text(aes(label = val, vjust = -0.8, hjust = 0.5), show.legend = TRUE)+
  ylab("-log2(Pvalue)")+
  ggtitle(paste("DApeaks_P0.05-logFC0-DESeq2-",freq_name,"_TFs",sep="")) +
  theme(plot.title = element_text(hjust = 0.5))
plot(p)
}
TFs(file_name,TF_names,freq,freq_name,a)
#############################################################################
######################  使用FIMO进行motif富集，并对应到TFs
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDAC_frontiers/4_TFs/FIMO_result")
library(data.table)
fimo<-fread("up05fimo.tsv",header=T,data.table=F)
colnames(fimo)<-c(colnames(fimo)[1:8],"FDR","matched_sequence")
library(dplyr)
TF_fimo<-fimo %>% distinct(motif_id,motif_alt_id,FDR, .keep_all = TRUE)
df<-TF_fimo[which(TF_fimo$FDR<0.05),]
length(unique(df$motif_alt_id))

TF_motif<-data.frame()
for(i in 1:length(unique(df[,2]))){
  TF_motif[i,1]<-unique(df[,2])[i]
  TF_motif[i,2]<-sum(df[,2] %in% unique(df[,2])[i])
}
TF_motif_new<-TF_motif[order(TF_motif[,2],decreasing = T),]
colnames(TF_motif_new)<-c("TF","number")

TF_motif_new$TF=factor(TF_motif_new$TF,levels = TF_motif_new[,1])
ggplot(data = TF_motif_new[1:20,],mapping = aes(x = TF[1:20], y =number[1:20]))+ 
  geom_bar(stat = 'identity',position="dodge")+
  geom_text(aes(label = number, vjust = -0.8, hjust = 0.5), show.legend = TRUE)+
  ylab("Top20-TFs")+
  xlab("TFs-names")+
  ggtitle("up05fimo_TFs") +
  theme(plot.title = element_text(hjust = 0.5))





######################################################################################################
##################  差异表达转录因子(DETF)及其与DApeaks相关的目标DEGs  上下调的转录因子是分开看的
rm(list=ls())
library(VennDiagram)
library(dplyr)
library(stringr)
file_name<-"S8495"
DE_gene<-"P0.05_S8495_up" #or down
TF_names<-"sen-res-S8495-0.05-P-logFC0-DESeq2-up"
freq<-0.05  #FDR<0.05
venn_name<-"up_DETFs_P005"
TF_target_name<-"TRANSFAC" #JASPAR, MotifMap, ENCODE
up_name<-"P0.05_S8495_up"
down_name<-"P0.05_S8495_down"
DE_name<-"P0.05_S8495_up"
regulation<-"up"

DETFs<-function(file_name,DE_name,DE_gene,TF_names,freq,venn_name,TF_target_name,up_name,down_name,regulation){
setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/1_RNAseq_DEGs",sep=""))
DE_gene<-read.table(paste("sensitivity-resistance-all-DESeq2_RNAseq_",DE_name,"_DEGs.txt",sep=""),header=T,sep="\t")
cat(paste("The DEGs number is",nrow(DE_gene),sep=" "),"\n")
setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/4_TFs/homer_result/",TF_names,sep=""))
knownResults<-read.delim("knownResults.txt",header=T,sep='\t')
FDR05<-knownResults[knownResults$P.value<freq,]
#FDR05<-knownResults[knownResults$q.value..Benjamini.<freq,]
homer_TF1<-data.frame(lapply(strsplit(as.character(FDR05[,1]),'/'), function(x) x[3])%>%unlist())
FDR051<-FDR05[which(homer_TF1=="Homer"),]
homer_TF2<-data.frame(lapply(strsplit(as.character(FDR051[,1]),'/'), function(x) x[1])%>%unlist())
homer_TF3<-apply(homer_TF2,1,function(x) str_replace_all(as.character(x),"\\((.*?)\\)",""))##取出括号前的字符
cat(paste("The number of motif is",length(homer_TF3),sep=" "),"\n")
homer_TF4<-unique(homer_TF3)
cat(paste("The number of TFs is",length(homer_TF4),sep=" "),"\n")
new_homer_result<-cbind(homer_TF3,FDR051[,-1])
colnames(new_homer_result)<-c("TF_name",colnames(new_homer_result)[-1])
DETFs<-intersect(toupper(new_homer_result$TF_name),toupper(DE_gene$gene_id))
cat(paste("The DETFs number is",length(DETFs),sep=" "),"\n")
cat(DETFs,"\n")
setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/3_intersect/DETFs_venn",sep=""))
venn.diagram(x=list(TF_enrich_hypo_DApeaks=toupper(new_homer_result$TF_name),down_DEGs=toupper(DE_gene$gene_id)), paste(venn_name,".png",sep=""),fill=c("red","blue"),margin = 0.3)
###############   TF找到对应的靶基因
setwd("~/xjj/drug/drug_result/HDAC_frontiers/4_TFs/TRANSFAC_TF")
target_TF<-read.delim(paste(TF_target_name,"-gene_attribute_edges.txt",sep=""),header=T,sep='\t')
target_TF_gene<-target_TF[,c(1,4)]
target_TF_gene<-target_TF_gene[-1,]
colnames(target_TF_gene)<-c("Target","TFs")
cat(intersect(DETFs,unique(target_TF_gene[,2])),"\n")
DETF_want<-intersect(DETFs,unique(target_TF_gene[,2]))
result<-NULL
for(i in 1:length(DETF_want)){
  result1<-target_TF_gene[target_TF_gene[,2] %in% DETF_want[i],]
  result<-rbind(result,result1)
}
### 查看靶标基因中有多少是DEGs
setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/1_RNAseq_DEGs",sep=""))
up_gene<-read.table(paste("sensitivity-resistance-all-DESeq2_RNAseq_",up_name,"_DEGs.txt",sep=""),header=T,sep="\t")
down_gene<-read.table(paste("sensitivity-resistance-all-DESeq2_RNAseq_",down_name,"_DEGs.txt",sep=""),header=T,sep="\t")
data<-NULL
for(i in 1:unique(result[,2])){
  target<-result[result[,2] %in% unique(result[,2])[i],1]
  DETFs_regulation_UP<-intersect(toupper(target),toupper(up_gene[,1]))
  cat(length(DETFs_regulation_UP),"\t")
  DETFs_regulation_DOWN<-intersect(toupper(target),toupper(down_gene[,1]))
  cat(length(DETFs_regulation_DOWN),"\t")
  cat(length(DETFs_regulation_UP)+length(DETFs_regulation_DOWN),"\t")
  one<-matrix(rep(unique(result[,2])[i],length(DETFs_regulation_UP)+length(DETFs_regulation_DOWN)),ncol=1)
  two<-matrix(c(DETFs_regulation_UP,DETFs_regulation_DOWN),ncol=1)
  three<-matrix(c(rep("up",length(DETFs_regulation_UP)),rep("down",length(DETFs_regulation_DOWN))),ncol=1)
  data1<-cbind(one,two,three)
  data<-rbind(data,data1)
}
colnames(data)<-c("DETFs","target-DEGs","target-DEGs-dysregulation")
cat(paste("The DETFs contain target is ",unique(data[,1]),sep=""),"\n")
setwd(paste("~/xjj/drug/drug_result/durg62/",file_name,"/4_TFs/DETFs-target-DEGs",sep=""))
write.table(data,paste("P0.05_",regulation,"_DEGs_DApeaks_FDR",freq,"_DETFs_",TF_target_name,"-gene.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
}
TF_target_name<-"TRANSFAC" #TRANSFAC,JASPAR, MotifMap, ENCODE
DETFs(file_name,DE_name,DE_gene,TF_names,freq,venn_name,TF_target_name,up_name,down_name,regulation)


#######################################################################
#####################    TCGA  生存分析
rm(list=ls())
library(org.Hs.eg.db)
library(stringr)
library(clusterProfiler)
library(survival)
library(survminer)
g<-"HOXA11"

setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/6_TCGA")
chemotherapy_drug_sample1<-read.table("drug_paad1.txt",header=T,sep="\t") #化疗药物信息
#chemotherapy_drug_sample<-chemotherapy_drug_sample1[which(chemotherapy_drug_sample1$pharmaceutical_therapy_drug_name=="Oxaliplatin"),]
chemotherapy_drug_sample<-chemotherapy_drug_sample1[which(chemotherapy_drug_sample1$pharmaceutical_therapy_drug_name %in% c("5-FU","5 FU","5FU","5-Fluorouracil","5-fu","5-fluorouracil","5-Fluorouracil?")),]
#chemotherapy_drug_sample<-read.table("drug_paad1.txt",header=T,sep="\t") #化疗药物信息
TCGA_exp<-read.table("gene_exp.txt",header=T,sep="\t",row.names=1) #TCGA的FPKM
TCGA_exp_colname1<-substring(colnames(TCGA_exp),1,12)
TCGA_exp_colname<-gsub("\\.","-",TCGA_exp_colname1) #提取符合条件的样本名
drug_sample<-intersect(chemotherapy_drug_sample[,1],TCGA_exp_colname) #有用药的样本
clinical_follow_up1<-read.table("clinical_follow_up.txt",header=T,sep="\t") #化疗药物信息
clinical_follow_up2<-clinical_follow_up1[-c(1:2),c(2,11,12,13)]
clinical_follow_up21<-clinical_follow_up2[which(clinical_follow_up2$vital_status!="[Not Available]"),]
clinical_follow_up3<-clinical_follow_up21[-61,]
a<-which(clinical_follow_up3[,3]=="[Not Available]")
b<-which(clinical_follow_up3[,3]!="[Not Available]")
death_time<-clinical_follow_up3[a,]
colnames(death_time)<-c("sample","state","nowant","time")
follow_time1<-clinical_follow_up3[b,]
follow_time<-follow_time1[,c(1,2,4,3)]
colnames(follow_time)<-c("sample","state","nowant","time")
state_time<-rbind(death_time,follow_time)##生存时间和生存状态
drug_sample_follow<-intersect(drug_sample,state_time[,1])##有生存有用药的样本
drug_sample_follow_exp<-TCGA_exp[,match(drug_sample_follow,TCGA_exp_colname)]
gene=bitr(rownames(drug_sample_follow_exp),fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") 
int_gene<-intersect(gene[,1],rownames(drug_sample_follow_exp))
exp1<-drug_sample_follow_exp[match(int_gene,rownames(drug_sample_follow_exp)),]
rownames(exp1)<-gene[match(int_gene,gene[,1]),2]
#median<-mean(as.numeric(exp1[match(g,rownames(exp1)),]))###均值
median<-median(as.numeric(exp1[match(g,rownames(exp1)),]))###中位数，基因需要换
up<-which(as.numeric(exp1[match(g,rownames(exp1)),])>=median)
down<-which(as.numeric(exp1[match(g,rownames(exp1)),])<median)
class_ind<-matrix(0,ncol(exp1),1)
class_ind[up,1]<-"high"   #基因表达更高，越敏感，理论上生存越好
class_ind[down,1]<-"low" #基因表达更低，越耐药，理论上生存越差
exp1_colname1<-substring(colnames(exp1),1,12)
exp1_colname<-gsub("\\.","-",exp1_colname1) #提取符合条件的样本名
surv_info<-state_time[match(exp1_colname,state_time[,1]),c(2,4)]#生存时间和生存状态
label<-as.matrix(class_ind)
TTP=as.numeric(surv_info[,2])  ##生存时间
status_TTP=as.matrix(surv_info[,1])  ##生存状态
status_TTP[which(status_TTP=="Alive")]=0
status_TTP[which(status_TTP=="Dead")]=1
status_TTP<-as.numeric(status_TTP)
surv_info1<-as.data.frame(cbind(TTP,status_TTP))
summary(coxph(Surv(TTP, status_TTP) ~ label,data=surv_info1)) 
surv_TTP<-survfit(Surv(TTP, status_TTP) ~ label,data=surv_info1)
ggsurvplot(surv_TTP,
           pval = TRUE, conf.int = TRUE, #是否显示p值/显示风险区域
           risk.table = TRUE # 将风险表显示在生存曲线下面
)


##########################################################################
#################     生存相关基因
rm(list=ls())
setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/6_TCGA")
chemotherapy_drug_sample1<-read.table("drug_paad1.txt",header=T,sep="\t") #化疗药物信息
chemotherapy_drug_sample<-chemotherapy_drug_sample1[which(chemotherapy_drug_sample1$pharmaceutical_therapy_drug_name %in% c("5-FU","5 FU","5FU","5-Fluorouracil","5-fu","5-fluorouracil","5-Fluorouracil?")),]

TCGA_exp<-read.table("gene_exp.txt",header=T,sep="\t",row.names=1) #TCGA的FPKM
TCGA_exp_colname1<-substring(colnames(TCGA_exp),1,12)
TCGA_exp_colname<-gsub("\\.","-",TCGA_exp_colname1) #提取符合条件的样本名
drug_sample<-intersect(chemotherapy_drug_sample[,1],TCGA_exp_colname) #有用药的样本
clinical_follow_up1<-read.table("clinical_follow_up.txt",header=T,sep="\t") #化疗药物信息
clinical_follow_up2<-clinical_follow_up1[-c(1:2),c(2,11,12,13)]
clinical_follow_up21<-clinical_follow_up2[which(clinical_follow_up2$vital_status!="[Not Available]"),]
clinical_follow_up3<-clinical_follow_up21[-61,]
a<-which(clinical_follow_up3[,3]=="[Not Available]")
b<-which(clinical_follow_up3[,3]!="[Not Available]")
death_time<-clinical_follow_up3[a,]
colnames(death_time)<-c("sample","state","nowant","time")
follow_time1<-clinical_follow_up3[b,]
follow_time<-follow_time1[,c(1,2,4,3)]
colnames(follow_time)<-c("sample","state","nowant","time")
state_time<-rbind(death_time,follow_time)##生存时间和生存状态

drug_sample_follow<-intersect(drug_sample,state_time[,1])##有生存有用药的样本
drug_sample_follow_exp<-TCGA_exp[,match(drug_sample_follow,TCGA_exp_colname)]

library(org.Hs.eg.db)
library(stringr)
library(clusterProfiler)
gene=bitr(rownames(drug_sample_follow_exp),fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") 
int_gene<-intersect(gene[,1],rownames(drug_sample_follow_exp))
exp1<-drug_sample_follow_exp[match(int_gene,rownames(drug_sample_follow_exp)),]
rownames(exp1)<-gene[match(int_gene,gene[,1]),2]
delete<-apply(exp1,1,function(x) mean(x==0))
cou<-which(delete>0.5)####在超过一半样本以上的基因删去
exp<-exp1[(-cou),]  #####0值小于一半样本的基因，且有药物信息，有生存信息
#exp<-exp1
exp1_colname1<-substring(colnames(exp),1,12)
exp1_colname<-gsub("\\.","-",exp1_colname1) #提取符合条件的样本名
surv_info<-state_time[match(exp1_colname,state_time[,1]),c(2,4)]#生存时间和生存状态
TTP=as.numeric(surv_info[,2])  ##生存时间
status_TTP=as.matrix(surv_info[,1])  ##生存状态
status_TTP[which(status_TTP=="Alive")]=0
status_TTP[which(status_TTP=="Dead")]=1
status_TTP<-as.numeric(status_TTP)
surv_info1<-as.data.frame(cbind(TTP,status_TTP))
library(survival)
library(survminer)
library(broom)
result1<-data.frame()
for(i in 1:nrow(exp)){
  median<-mean(as.numeric(exp[i,]))###中位数，基因需要换result1是均值，result是中位数
  up<-which(as.numeric(exp[i,])>=median)
  down<-which(as.numeric(exp[i,])<median)
  class_ind<-matrix(0,ncol(exp),1)
  class_ind[up,1]<-"high"   #基因表达更高，越敏感，理论上生存越好
  class_ind[down,1]<-"low" #基因表达更低，越耐药，理论上生存越差
  label<-as.matrix(class_ind)
  coxp<-coxph(Surv(TTP, status_TTP) ~ label,data=surv_info1)
  result1[i,1]<-rownames(exp)[i]
  result1[i,2]<-tidy(coxp)$estimate
  result1[i,3]<-tidy(coxp)$p.value
}
result1[,4]=p.adjust(result1[,3],method="BH")
colnames(result1)<-c("geneid","regression coefficient","pvalue","FDR")
surv_coff_gene1<-result1[which(result1[,3]<0.05),]#卡FDR时没有，所以卡P

setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/1_RNAseq_DEGs")
up_gene<-read.table("sensitivity-resistance-all-DESeq2_RNAseq_P0.05_up_DEGs.txt",header=T,sep="\t")
down_gene<-read.table("sensitivity-resistance-all-DESeq2_RNAseq_P0.05_down_DEGs.txt",header=T,sep="\t")
all_DEGs<-rbind(up_gene,down_gene)
DEGs<-all_DEGs$gene_id

setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/2_ATACseq_DApeaks")
library(data.table)
homer_anno_down<-fread("sen-res-chemotherapy-0.05-P-down-annotation.txt",header=T,data.table=F)
colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]

homer_anno_up<-fread("sen-res-chemotherapy-0.05-P-up-annotation.txt",header=T,data.table=F)
colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
all_DApeaks_gene<-rbind(Anno_gene_100Kb_up,Anno_gene_100Kb_down)
DApeaks_gene<-unique(all_DApeaks_gene$`Gene Name`)

#input  <-list(unique(surv_coff_gene1[,1]),unique(DEGs),unique(DApeaks_gene))
input  <-list(as.character(unique(surv_coff_gene1[,1])),as.character(unique(DEGs)))
venn(input,showSetLogicLabel=TRUE)
tmp <- venn(input)
int<-attr(tmp, "intersections")

setwd("~/xjj/drug/drug_result/durg62/5-FU/5_TCGA")
library(VennDiagram)
venn.diagram(x=list(Survival_related_genes=unique(surv_coff_gene1[,1]),RNA_DEGs=unique(DEGs)), "surv_coff_DEGs0.05.png",fill=c("red","blue"))

##对交叠的基因使用survfit()函数进行KM生存分析
intt<-int[["A:B"]]
expkm<-exp[match(intt,rownames(exp)),]
library(survival)
library(survminer)
library(broom)
result2<-data.frame()
for(i in 1:nrow(expkm)){
  median<-median(as.numeric(expkm[i,]))###中位数，基因需要换result1是均值，result是中位数
  up<-which(as.numeric(expkm[i,])>=median)
  down<-which(as.numeric(expkm[i,])<median)
  class_ind<-matrix(0,ncol(expkm),1)
  class_ind[up,1]<-"high"   #基因表达更高，越敏感，理论上生存越好
  class_ind[down,1]<-"low" #基因表达更低，越耐药，理论上生存越差
  label<-as.matrix(class_ind)
  surv_TTP<-survdiff(Surv(TTP, status_TTP) ~ label,data=surv_info1)
  p.value <- 1 - pchisq(surv_TTP$chisq, length(surv_TTP$n) -1)
  result2[i,1]<-rownames(expkm)[i]
  result2[i,2]<-p.value
}
colnames(result2)<-c("geneid","pvalue")
surv_coff_gene2<-result2[which(result2[,2]<0.05),]
setwd("~/xjj/drug/drug_result/durg62/5-FU/5_TCGA")
write.table(surv_coff_gene2,"surv_coff_DEGs.txt",sep="\t",col.names=T,row.names=F)



###################################################################################
################   分类器 randomForest
rm(list=ls())
setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
library("data.table")
exp1<-fread("star_rsem.GeneSymbol.TPM.xls",header=T,data.table=F)
exp2<-log2(exp1[-1]+1)
delete<-apply(exp2,1,function(x) mean(x==0))
cou<-which(delete>0.5)####在超过一半样本以上的基因删去
exp<-exp2
sensitivity<-c("PC.64","PC.81","PC.116","PC.L","PC.2","PC.115","PC.16","PC.97","PC.101","PC.27","PC.8","PC.109","PC.136","PC.13","PC.78","PC.98","PC.130","PC.117","PC.139","PC.52")  ##20
resistance<-c("PC.G","PC.111","PC.134","PC.14","PC.40","PC.135","PC.22","PC.18","PC.105","PC.104","PC.56","PC.119","PC.112","PC.102","PC.121","PC.5","PC.I")  ##17
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]
resistance_exp<-exp[,match(resistance_new,colnames(exp))]
sensitivity_exp<-exp[,match(sensitivity_new,colnames(exp))]
geneid<-exp1[,1]

data<-cbind(geneid,resistance_exp,sensitivity_exp) ####37个样本的表达谱
symbol<-as.matrix(data[,1])####基因名
glist<-symbol
exp<-t(data[,-1])  ###纯表达谱
#colnames(exp)<-paste("g",1:length(symbol),sep="")
colnames(exp)<-symbol
exp<-as.data.frame(exp)
gt<-rep(c(0,1),c(17,20)) #ground_truth 耐药17 敏感20

setwd("~/xjj/drug/drug_result/HDAC_frontiers/1_RNAseq_DEGs/Acetylase")
Acetylase_info<-read.table("Acetylase_info_P0.05.txt",sep = '\t',header = T)
k=100
loc<-matrix(1,37,k)  ###构建一张37行100列的矩阵，用于做测试集的样本位置
rownames(loc)<-rownames(exp) ##37个样本名
set.seed(123)
for(i in 1:k){
  rr<-sample(1:17,4)  ##耐药中随机选取2个作为作为验证集
  sr<-sample(18:37,4) ##敏感中随机选取2个样本作为验证集
  loc[c(rr,sr),i]<-2 #测试集
}
IMP<-matrix(0,length(symbol),2*k)  ###基因名长度，有200列
colnames(IMP)[seq(2,2*k,2)]<-"MeanDecreaseGini"
colnames(IMP)[seq(1,2*k,2)]<-"DRG"
rownames(IMP)<-symbol
ability<-matrix(NA,6,k)
rownames(ability)<-c("sen","spe","acc","val_sen","val_sep","val_acc")

####
for(j in 1:k)
{
  print(j)
  train.data<-exp[which(loc[,j]==1),]
  train.gt<-gt[which(loc[,j]==1)]
  test.data<-exp[which(loc[,j]==2),]
  test.gt<-gt[which(loc[,j]==2)]
  
  dat<-cbind(train.gt,train.data)
  Acetylase_DEGs<-Acetylase_info[,1]
  IMP[match(Acetylase_DEGs,colnames(exp)),2*j-1]<-1

  dat1<-dat[,c(1,match(Acetylase_DEGs,colnames(dat)))]
  dat1$train.gt[dat1$train.gt==0]="resistant"
  dat1$train.gt[dat1$train.gt==1]="sensitive"
  dat1$train.gt<-as.factor(dat1$train.gt)
  
  library(randomForest)
  set.seed(123)
  dat1.rf <- randomForest(train.gt ~ ., data=dat1, importance=T,proximity=T)
  print(dat1.rf) #展示随机森林模型简要信息
  ability[1,j]<-dat1.rf$ confusion[1,1]/13
  ability[2,j]<-dat1.rf$ confusion[2,2]/16
  ability[3,j]<-(dat1.rf$ confusion[1,1]+dat1.rf$ confusion[2,2])/33
  imp<-dat1.rf$importance
  IMP[match(rownames(imp),colnames(exp)),2*j]<-imp[,4]
  #提取随机森林模型中以准确率递减方法得到维度重要性值。type=2为基尼系数方法
  varImpPlot(dat1.rf,sort=TRUE,n.var=nrow(dat1.rf$importance), main = "variable importance")
  hist(treesize(dat1.rf))   #展示随机森林模型中每棵决策树的节点数
  max(treesize(dat1.rf));min(treesize(dat1.rf))
  MDSplot(dat1.rf, dat1$train.gt, palette=rep(1, 2), pch=as.numeric(dat1$train.gt)) #展示数据集在二维情况下各类别的具体分布情况
  pre<- predict(dat1.rf,newdata=test.data)
  print(pre)
  pred_out_1<-predict(object=dat1.rf,newdata=test.data,type="prob")
  print(pred_out_1)
  ability[4,j]<-length(grep("res",pre[1:4]))/4
  ability[5,j]<-length(grep("sen",pre[5:8]))/4
  ability[6,j]<-(length(grep("res",pre[1:4]))+length(grep("sen",pre[5:8])))/8
}
IMP1<-IMP[match(Acetylase_info[,1],colnames(exp)),]
setwd("~/xjj/drug/drug_result/HDAC_frontiers/6_classifier/randomForest/Acetylase_DEGs")
write.table(loc,"Acetylase_DEGs_train_test4.txt",sep="\t",col.names=F)
write.table(IMP1,"Acetylase_DEGs_MeanDecreaseGini4.txt",sep="\t",col.names=F)
write.table(ability,"Acetylase_DEGs_train_test_ability4.txt",sep="\t",col.names=F)

###############    分类器  SVM
rm(list=ls())
setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
library("data.table")
exp1<-fread("star_rsem.GeneSymbol.TPM.xls",header=T,data.table=F)
exp2<-log2(exp1[-1]+1)
delete<-apply(exp2,1,function(x) mean(x==0))
cou<-which(delete>0.5)####在超过一半样本以上的基因删去
exp<-exp2
sensitivity<-c("PC.64","PC.81","PC.116","PC.L","PC.2","PC.115","PC.16","PC.97","PC.101","PC.27","PC.8","PC.109","PC.136","PC.13","PC.78","PC.98","PC.130","PC.117","PC.139","PC.52")  ##20
resistance<-c("PC.G","PC.111","PC.134","PC.14","PC.40","PC.135","PC.22","PC.18","PC.105","PC.104","PC.56","PC.119","PC.112","PC.102","PC.121","PC.5","PC.I")  ##17
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]
resistance_exp<-exp[,match(resistance_new,colnames(exp))]
sensitivity_exp<-exp[,match(sensitivity_new,colnames(exp))]
geneid<-exp1[,1]

data<-cbind(geneid,resistance_exp,sensitivity_exp) ####37个样本的表达谱
symbol<-as.matrix(data[,1])####基因名
glist<-symbol
exp<-t(data[,-1])  ###纯表达谱
#colnames(exp)<-paste("g",1:length(symbol),sep="")
colnames(exp)<-symbol
exp<-as.data.frame(exp)
gt<-rep(c(0,1),c(17,20)) #ground_truth 耐药17 敏感20

setwd("~/xjj/drug/drug_result/HDAC_frontiers/1_RNAseq_DEGs/Acetylase")
Acetylase_info<-read.table("Acetylase_info_P0.05.txt",sep = '\t',header = T)

#随机100次
k=100
loc<-matrix(1,37,k)
rownames(loc)<-rownames(exp)
set.seed(123)
for(i in 1:k){
  rr<-sample(1:17,4)
  sr<-sample(18:37,4)
  loc[c(rr,sr),i]<-2 #测试集
}

IMP<-matrix(0,length(symbol),k)
rownames(IMP)<-symbol
ability<-matrix(NA,6,k)
rownames(ability)<-c("sen","spe","acc","val_sen","val_sep","val_acc")

library(e1071)
for(j in 1:k){
  print(j)
  train.data<-exp[which(loc[,j]==1),]
  train.gt<-gt[which(loc[,j]==1)]
  test.data<-exp[which(loc[,j]==2),]
  test.gt<-gt[which(loc[,j]==2)]
  
  dat<-cbind(train.gt,train.data)
  Acetylase_DEGs<-Acetylase_info[,1]
  IMP[match(Acetylase_DEGs,colnames(exp)),j]<-1
  dat1<-dat[,c(1,match(Acetylase_DEGs,colnames(dat)))]
  dat1$train.gt<-as.factor(dat1$train.gt)
  
  set.seed(100) # for reproducing results
  tuned<-tune.svm(train.gt~.,data = dat1,gamma = 5^(-6:-1),cost = 5^(0:4))
  summary(tuned)
  model.tuned<-svm(train.gt~.,data = dat1,gamma=tuned$best.parameters$gamma,cost=tuned$best.parameters$cost)
  print(summary(model.tuned))
  pred<-predict(model.tuned,dat1)
  
  ability[1,j]<-length(grep(0,pred[1:13]))/13
  ability[2,j]<-length(grep(1,pred[14:29]))/16
  ability[3,j]<-(length(grep(0,pred[1:13]))+length(grep(1,pred[14:29])))/29
  
  pre<- predict(model.tuned,newdata=test.data)
  print(pre)
  ability[4,j]<-length(grep(0,pre[1:4]))/4
  ability[5,j]<-length(grep(1,pre[5:8]))/4
  ability[6,j]<-(length(grep(0,pre[1:4]))+length(grep(1,pre[5:8])))/8
}
IMP1<-IMP[match(Acetylase_info[,1],colnames(exp)),]
setwd("~/xjj/drug/drug_result/HDAC_frontiers/6_classifier/SVM/Acetylase_DEGs")
write.table(loc,"Acetylase_DEGs_train_test4.txt",sep="\t",col.names=T)
write.table(IMP1,"Acetylase_DEGs_MeanDecreaseGini4.txt",sep="\t",col.names=T)
write.table(ability,"Acetylase_DEGs_train_test_ability4.txt",sep="\t",col.names=T)

###########  分类器   Logistic regression model











