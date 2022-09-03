#### HDAC药物聚类
rm(list=ls())
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
setwd("~/xjj/drug")
ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,!(colnames(ic501)%in%c("PC.100","PC.34"))]## all of PDAC
colnames(ic50)<-sample_id[match(colnames(ic50),sample_id[,2]),1]

drug_info<-read.csv("compiled_drug_info.csv",header = TRUE,sep = ",")
drug_result<-drug_info[match(rownames(ic50),drug_info$drug_id),]
HDAC1<-drug_result[drug_result$target=="HDAC",]
aaa<-c("S1848","S2759","S1194","S1047")
HDAC2<-drug_result[drug_result$drug_id%in%aaa,]
HDAC3<-rbind(HDAC1,HDAC2)
HDAC4<-HDAC3[order(HDAC3[,1],decreasing = T),]
HDAC<-HDAC4[-which(HDAC4$drug_id=="S1848"|HDAC4$drug_id=="S8495"),]

HDAC_ic50<-ic50[match(as.character(HDAC[,1]),rownames(ic50)),]
library(pheatmap)
a=cor(t(HDAC_ic50))
result=pheatmap(a,scale = "none",main = "The correlation coefficients between HDACi and chemotherapy' AUC",show_rownames=T,show_colnames=T,
                clustering_distance_rows = "correlation",
                clustering_distance_cols = "correlation",clustering_method = "complete",cutree_rows=2,cutree_cols=2)

##########   20种HDAC+5种化疗药
bb<-t(apply(HDAC_ic50,1,function(x) scale(x)))
colnames(bb)<-colnames(HDAC_ic50)
norma_result<-pheatmap(HDAC_ic50,scale = "none",main = "The AUC of GEM and S8495 are no-normalized",show_rownames=T,show_colnames=T,
                       clustering_distance_rows = "correlation",cluster_rows=FALSE,
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
which(tresult[,4]<0.05)
length(which(tresult[,4]<0.05))
### wilcox.test
tresult=matrix(0,length(geneid),4)
for (i in 1:length(geneid)){
  ttest=wilcox.test(as.numeric(sensitivity_AUC[i,]),as.numeric(resistance_AUC[i,]),alternative ="two.sided")
  tresult[i,1:3]=c(geneid[i],ttest$statistic,ttest$p.value)
}
tresult[,4]=p.adjust(tresult[,3],method="BH")
colnames(tresult)<-c("geneid","statistic","p.value","FDR")
which(tresult[,4]<0.05)
length(which(tresult[,4]<0.05))

################################## 化疗药物在一个药物以上敏感划分为敏感组
rm(list=ls())
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
setwd("~/xjj/drug")
ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,!(colnames(ic501)%in%c("PC.100","PC.34"))]## all of PDAC
colnames(ic50)<-sample_id[match(colnames(ic50),sample_id[,2]),1]

setwd("~/xjj/drug/drug_result/drug64_AUC/newresult20220515")
tresult_want<-read.table("drug64_sample220515.txt",sep="\t",header=T)
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")

setwd("~/xjj/drug")
drug_info<-read.csv("compiled_drug_info.csv",header = TRUE,sep = ",")
ccc<-c("5-FU","GEM","IRI","OXA","PAC")
result<-matrix(0,length(ccc),length(colnames(ic50)))
for(j in 1:length(ccc)){
  G<-ccc[j]
  common_sample<-as.data.frame(tresult_want[match(G,tresult_want[,1]),])
  sensitivity<-unlist(strsplit(as.character(common_sample[1,5]), ","))
  resistance<-unlist(strsplit(as.character(common_sample[1,4]), ","))
  result[j,match(sensitivity,colnames(ic50))]<-1
}
sen_sum<-apply(result,2,sum)
colnames(result)<-colnames(ic50)

sensitivity<-colnames(result)[which((sen_sum>=1))] #### 19 samples
resistance<-colnames(result)[-match(sensitivity,colnames(result))]


################################## HDACI药物在一个药物以上敏感划分为敏感组
rm(list=ls())
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
setwd("~/xjj/drug")
ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,!(colnames(ic501)%in%c("PC.100","PC.34"))]## all of PDAC
colnames(ic50)<-sample_id[match(colnames(ic50),sample_id[,2]),1]

setwd("~/xjj/drug/drug_result/drug64_AUC/newresult20220515")
tresult_want<-read.table("drug64_sample220515.txt",sep="\t",header=T)
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")

setwd("~/xjj/drug")
drug_info<-read.csv("compiled_drug_info.csv",header = TRUE,sep = ",")
HDAC1<-drug_info[drug_info$target=="HDAC",]
ccc1<-c(as.character(HDAC1[,1]),"S1848","S2759","S1194","S1047")##去掉S1848，因为它在两类样本中AUC差异不显著
ccc<-ccc1[-match(c("S8495","S1096"),ccc1)]
ccc<-ccc1[-match(c("S8495","S1848"),ccc1)]

result<-matrix(0,length(ccc),length(colnames(ic50)))
for(j in 1:length(ccc)){
  G<-ccc[j]
  common_sample<-as.data.frame(tresult_want[match(G,tresult_want[,1]),])
  sensitivity<-unlist(strsplit(as.character(common_sample[1,5]), ","))
  resistance<-unlist(strsplit(as.character(common_sample[1,4]), ","))
  result[j,match(sensitivity,colnames(ic50))]<-1
}
sen_sum<-apply(result,2,sum)
colnames(result)<-colnames(ic50)

sensitivity<-colnames(result)[which((sen_sum>=1))] #### 30 samples
resistance<-colnames(result)[-match(sensitivity,colnames(result))]
#resistance<-colnames(result)[which((sen_sum==0))]
sensitivity_HDAC<-sensitivity
resistance_HDAC<-resistance
HDAC_ic50<-ic50[match(ccc,rownames(ic50)),]
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
which(tresult[,4]<0.05) # all of 18 drugs are different
length(which(tresult[,4]<0.05))
### wilcox.test
tresult=matrix(0,length(geneid),4)
for (i in 1:length(geneid)){
  ttest=wilcox.test(as.numeric(sensitivity_AUC[i,]),as.numeric(resistance_AUC[i,]),alternative ="two.sided")
  tresult[i,1:3]=c(geneid[i],ttest$statistic,ttest$p.value)
}
tresult[,4]=p.adjust(tresult[,3],method="BH")
colnames(tresult)<-c("geneid","statistic","p.value","FDR")
which(tresult[,4]<0.05) # all of 18 drugs are different
length(which(tresult[,4]<0.05))


###################  设置一个值作为判定标准  ####################
rm(list=ls())
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
setwd("~/xjj/drug")
ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,!(colnames(ic501)%in%c("PC.100","PC.34"))]## all of PDAC
colnames(ic50)<-sample_id[match(colnames(ic50),sample_id[,2]),1]

setwd("~/xjj/drug/drug_result/drug64_AUC/newresult20220515")
tresult_want<-read.table("drug64_sample220515.txt",sep="\t",header=T)
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
library(DESeq2)
library(dplyr)
library(foreach)

setwd("~/xjj/drug")
drug_info<-read.csv("compiled_drug_info.csv",header = TRUE,sep = ",")
HDAC1<-drug_info[drug_info$target=="HDAC",]
ccc<-c(as.character(HDAC1[,1]),"S1848","S2759","S1194","S1047")##去掉S1848，因为它在两类样本中AUC差异不显著
bbb<-c("5-FU","GEM","IRI","OXA","PAC")  ##5种化疗药
ddd<-c(ccc,bbb)
want_ic50<-ic50[match(ddd,rownames(ic50)),]

result<-NULL
for(i in 1:nrow(want_ic50)){
  a<-t(want_ic50[i,])
  b<-matrix(rep(rownames(want_ic50)[i],length(a)),ncol=1)
  result1<-cbind(b,a)
  result<-rbind(result,result1)
}
colnames(result)<-c("group","value")
library(foreign)
library(ggplot2)
dat<-as.data.frame(result)
ggplot(dat,aes(x=group,y=as.numeric(as.character(value))))+
  geom_boxplot()+
  theme_bw()+
  geom_rect(aes(xmin=0.5,xmax=5.5,ymin=0,ymax=Inf),
            fill='grey80',color='grey80')+
  geom_boxplot()+
  geom_point(size=0.8)

#假定AUC的值为0.5

want_ic50<-want_ic50[order(rownames(want_ic50)),]
result_sen<-matrix(0,25,2)
for(k in 1:25){
  result_sen[k,1]<-rownames(want_ic50)[k]
  result_sen[k,2]<-paste0(colnames(want_ic50)[which(as.numeric(want_ic50[k,])>0.55)],collapse =",")
  
}
result_sen1<-result_sen[which(result_sen[,2]!=""),2]
result_sen_sample<-unique(unlist(strsplit(as.character(result_sen1), ",")))
length(result_sen_sample)##22
sensitivity=result_sen_sample
resistance=colnames(want_ic50)[-match(sensitivity,colnames(want_ic50))]


#####  单独在5种化疗药物中看，大于特定的AUC则样本为敏感组，否则为耐药组
want_ic50<-want_ic50[order(rownames(want_ic50)),]
want_ic501<-want_ic50[rownames(want_ic50)%in%bbb,]
result_chemo_sen<-matrix(0,nrow(want_ic501),2)
for(k in 1:nrow(want_ic501)){
  result_chemo_sen[k,1]<-rownames(want_ic501)[k]
  result_chemo_sen[k,2]<-paste0(colnames(want_ic501)[which(as.numeric(want_ic501[k,])>0.5)],collapse =",")
}
result_chemo_sen1<-result_chemo_sen[which(result_chemo_sen[,2]!=""),2]
result_chemo_sen_sample<-unique(unlist(strsplit(as.character(result_chemo_sen1), ",")))
length(result_chemo_sen_sample)##24

#####  单独在20种HDAC药物中看，大于特定的AUC则样本为敏感组，否则为耐药组
want_ic50<-want_ic50[order(rownames(want_ic50)),]
want_ic503<-want_ic50[rownames(want_ic50)%in%ccc,]
want_ic502<-want_ic503[-which(rownames(want_ic503)=="S1848"|rownames(want_ic503)=="S8495"),]

result_HDAC_sen<-matrix(0,nrow(want_ic502),2)
for(k in 1:nrow(want_ic502)){
  result_HDAC_sen[k,1]<-rownames(want_ic502)[k]
  result_HDAC_sen[k,2]<-paste0(colnames(want_ic502)[which(as.numeric(want_ic502[k,])>0.5)],collapse =",")
}
result_HDAC_sen1<-result_HDAC_sen[which(result_HDAC_sen[,2]!=""),2]
result_HDAC_sen_sample<-unique(unlist(strsplit(as.character(result_HDAC_sen1), ",")))
length(result_HDAC_sen_sample)##25
sen<-intersect(result_HDAC_sen_sample,result_chemo_sen_sample)
res<-colnames(want_ic50)[-match(sen,colnames(want_ic50))]

sen<-union(result_HDAC_sen_sample,result_chemo_sen_sample)
res<-colnames(want_ic50)[-match(sen,colnames(want_ic50))]


sensitivity=result_HDAC_sen_sample
resistance=colnames(want_ic50)[-match(sensitivity,colnames(want_ic50))]

sensitivity_AUC<-want_ic502[,match(sensitivity,colnames(want_ic502))]
resistance_AUC<-want_ic502[,match(resistance,colnames(want_ic502))]
### T.test
geneid<-rownames(want_ic502)
tresult=matrix(0,length(geneid),4)
for (i in 1:length(geneid)){
  ttest=t.test(sensitivity_AUC[i,],resistance_AUC[i,],alternative = "two.sided")
  tresult[i,1:3]=c(geneid[i],ttest$statistic,ttest$p.value)
}
tresult[,4]=p.adjust(tresult[,3],method="BH")
colnames(tresult)<-c("geneid","statistic","p.value","FDR")
which(tresult[,4]<0.05)
### wilcox.test
tresult=matrix(0,length(geneid),4)
for (i in 1:length(geneid)){
  ttest=wilcox.test(as.numeric(sensitivity_AUC[i,]),as.numeric(resistance_AUC[i,]),alternative ="two.sided")
  tresult[i,1:3]=c(geneid[i],ttest$statistic,ttest$p.value)
}
tresult[,4]=p.adjust(tresult[,3],method="BH")
colnames(tresult)<-c("geneid","statistic","p.value","FDR")
which(tresult[,4]<0.05) #S1848

##############  将样本分成四个区
#散点图
rm(list=ls())
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
setwd("~/xjj/drug")
ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,!(colnames(ic501)%in%c("PC.100","PC.34"))]## all of PDAC
colnames(ic50)<-sample_id[match(colnames(ic50),sample_id[,2]),1]

setwd("~/xjj/drug/drug_result/drug64_AUC/newresult20220515")
tresult_want<-read.table("drug64_sample220515.txt",sep="\t",header=T)
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")

setwd("~/xjj/drug")
drug_info<-read.csv("compiled_drug_info.csv",header = TRUE,sep = ",")
HDAC1<-drug_info[drug_info$target=="HDAC",]
ccc<-c(as.character(HDAC1[,1]),"S1848","S2759","S1194","S1047")##去掉S1848，因为它在两类样本中AUC差异不显著
bbb<-c("5-FU","GEM","IRI","OXA","PAC")  ##5种化疗药
want_AUC_chemo<-ic50[rownames(ic50)%in%bbb,]
want_AUC_HDAC<-ic50[rownames(ic50)%in%ccc,]

#####每个样本在一类药物中最敏感的AUC值
AUC_chemo<-apply(want_AUC_chemo,2,max)
AUC_HDAC<-apply(want_AUC_HDAC,2,max)
data <- data.frame(x=AUC_chemo,y=AUC_HDAC)
#作散点图

freq<-0.55
want<-rbind(data[data$y>freq&data$x<freq,],data[data$y<freq&data$x>freq,])
ggplot(data, aes(x=x, y=y)) + 
  geom_point()+
  geom_hline(yintercept= freq)+
  geom_vline(xintercept=freq)+
  geom_point(data = data[data$y>freq&data$x<freq,],size = 4)+
  geom_point(data = data[data$y<freq&data$x>freq,],size = 4)+
  xlab("chemotherapeutic") +
  ylab("HDACIs")+
  geom_text_repel(
    data = want,
    aes(label = rownames(want)),
    size = 3,
    color = "black",
    segment.color = "black", show.legend = FALSE )

HDAC_sen<-data[data$y>freq&data$x<freq,]
chemo_sen<-data[data$y<freq&data$x>freq,]
all_sen<-data[data$y>freq&data$x>freq,]
all_res<-data[data$y<freq&data$x<freq,]
##画热图
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
setwd("~/xjj/drug")
ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,!(colnames(ic501)%in%c("PC.100","PC.34"))]## all of PDAC
colnames(ic50)<-sample_id[match(colnames(ic50),sample_id[,2]),1]

drug_info<-read.csv("compiled_drug_info.csv",header = TRUE,sep = ",")
drug_result<-drug_info[match(rownames(ic50),drug_info$drug_id),]
HDAC1<-drug_result[drug_result$target=="HDAC",]
aaa<-c("S1848","S2759","S1194","S1047","5-FU","GEM","IRI","OXA","PAC")
HDAC2<-drug_result[drug_result$drug_id%in%aaa,]
HDAC3<-rbind(HDAC1,HDAC2)
HDAC<-HDAC3[order(HDAC3[,1],decreasing = T),]

HDAC_ic50<-ic50[match(as.character(HDAC[,1]),rownames(ic50)),]
##########   20种HDAC+5种化疗药
HDAC_sen<-data[data$y>freq&data$x<freq,]
chemo_sen<-data[data$y<freq&data$x>freq,]
all_sen<-data[data$y>freq&data$x>freq,]
all_res<-data[data$y<freq&data$x<freq,]

chemo_sen_AUC<-HDAC_ic50[,match(rownames(chemo_sen),colnames(HDAC_ic50))]
chemo_result<-pheatmap(chemo_sen_AUC,scale = "none",main = "Sensitive to chemotherapy and resistant to HDACIs",show_rownames=T,show_colnames=T,
                       clustering_distance_rows = "correlation",cluster_rows=FALSE,cluster_cols=FALSE,
                       clustering_distance_cols = "euclidean",clustering_method = "ward.D2")
HDAC_sen_AUC<-HDAC_ic50[,match(rownames(HDAC_sen),colnames(HDAC_ic50))]
HDAC_result<-pheatmap(HDAC_sen_AUC,scale = "none",main = "Sensitive to HDACIs and resistant to chemotherapy",show_rownames=T,show_colnames=T,
                      clustering_distance_rows = "correlation",cluster_rows=FALSE,cluster_cols=FALSE,
                      clustering_distance_cols = "euclidean",clustering_method = "ward.D2")


all_sen_AUC<-HDAC_ic50[,match(rownames(all_sen),colnames(HDAC_ic50))]
all_sen_result<-pheatmap(all_sen_AUC,scale = "none",main = "Sensitive to HDACIs and chemotherapy",show_rownames=T,show_colnames=T,
                      clustering_distance_rows = "correlation",cluster_rows=FALSE,cluster_cols=FALSE,
                      clustering_distance_cols = "euclidean",clustering_method = "ward.D2")

all_res_AUC<-HDAC_ic50[,match(rownames(all_res),colnames(HDAC_ic50))]
all_res_result<-pheatmap(all_res_AUC,scale = "none",main = "Resistant to HDACIs and chemotherapy",show_rownames=T,show_colnames=T,
                         clustering_distance_rows = "correlation",cluster_rows=FALSE,cluster_cols=FALSE,
                         clustering_distance_cols = "euclidean",clustering_method = "ward.D2",
                         color = colorRampPalette(c("#4169E1", "white"))(100))


################## 筛选HDAC's biomarker ###########################
#######################    ATACseq   ################################################
#############  DApeaks
rm(list=ls())
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_count_name<-read.table("peak_count.txt",sep="\t",header=T) #没有全0的

setwd("~/xjj/drug/drug_result/HDAC20220524/0_sample")
result<-read.table("HDAC_19drug_37sample.txt",sep = '\t',header = T)
colnames(result)<-gsub("\\.","-",colnames(result))
sen_sum<-apply(result,2,sum)
sensitivity<-colnames(result)[which((sen_sum>=10))] #### 19 samples
resistance<-colnames(result)[-match(sensitivity,colnames(result))] #### 18 samples

resistance_peak<-peak_count_name[,match(resistance,gsub("\\.","-",colnames(peak_count_name)))]
sensitivity_peak<-peak_count_name[,match(sensitivity,gsub("\\.","-",colnames(peak_count_name)))]
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
dds##这是一个关于各种内容的矩阵
dds<-DESeq(dds)##进行标准化分析
#rld <- rlogTransformation(dds)  ## 得到经过DESeq2软件normlization的表达矩阵！
#exprSet_new=assay(rld)
#rownames(exprSet_new)<-peak_name
sizeFactors(dds)##查看每个主成分的标准化值
res<-results(dds)##将结果输出
head(res)
class(res)##可以看出它是DESeq的属性，要转化为表格
res<-as.data.frame(res)
head(res)
res<-cbind(peak_count_name[,1],res)  ##对数据增加一列
head(res)
colnames(res)<- c('peak_id',"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj" )##对列名进行修改
setwd("~/xjj/drug/drug_result/HDAC20220524/2_ATACseq")
write.table(res,"sensitivity-resistance-all-DESeq2_HDAC_AUC_ATAC.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
#write.table(exprSet_new,"sensitivity-resistance-HDAC-DESeq2_normlization.txt",sep = '\t',col.names = T,row.names = T,quote = FALSE)##数据输出

#setwd("~/xjj/drug/drug_result/HDAC_drug/heatmap")
res<-read.table("sensitivity-resistance-all-DESeq2_HDAC_chemotherapy_AUC_ATAC.txt",header = TRUE,sep = "\t")
#resSig<-res[which(res$padj<0.05 & abs(res$log2FoldChange>1)),]
resSig<-res[which(res$pvalue<0.05 & abs(res$log2FoldChange)>1),]
resSig<-res[which(res$pvalue<0.05),]

# pvalue padj
##新增一列，将log2FoldChange>0标注为up，<0标准为down
resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
sum(resSig$up_down=='up')
sum(resSig$up_down=='down')
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
out_file<-"sen-res-HDAC-AUC-0.05-P-logFC0"
write.table(peaks_up,file=paste(out_file,"-DESeq2-up.txt",sep=""),sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')
write.table(d_up,file=paste(out_file,"-DESeq2-up.gff",sep=""),sep = '\t',
            col.names = F,row.names = F,quote = FALSE)
a_down<-apply(down_gene,1,function(x) unlist(strsplit(as.character(x), "_")))
d_down<-t(a_down)
e_down<-data.frame(rep("+",nrow(d_down)))
peaks_down<-cbind(down_gene,d_down,e_down)
colnames(peaks_down)<-c("peak_id","chr","start","end","strand")
write.table(peaks_down,file=paste(out_file,"-DESeq2-down.txt",sep=""),sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')
write.table(d_down,file=paste(out_file,"-DESeq2-down.gff",sep=""),sep = '\t',
            col.names = F,row.names = F,quote = FALSE)
#####  bed
f_up<-data.frame(rep(1,nrow(d_up)))
peaks_bed_up<-cbind(d_up,up_gene,f_up,e_up)
colnames(peaks_bed_up)<-c("chr","start","end","peak_id","score","strand")
write.table(peaks_bed_up,file=paste(out_file,"-DESeq2-up.bed",sep=""),sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')
f_down<-data.frame(rep(1,nrow(d_down)))
peaks_bed_down<-cbind(d_down,down_gene,f_down,e_down)
colnames(peaks_bed_down)<-c("chr","start","end","peak_id","score","strand")
write.table(peaks_bed_down,file=paste(out_file,"-DESeq2-down.bed",sep=""),sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')

####################################################################
rm(list=ls())
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_count_name<-read.table("peak_count.txt",sep="\t",header=T)
a<-as.matrix(peak_count_name[,1])
a_up<-apply(a,1,function(x) unlist(strsplit(as.character(x), "_")))
d_up<-t(a_up)
e_up<-data.frame(rep("+",nrow(d_up)))
f_up<-data.frame(rep(1,nrow(d_up)))
peaks_bed_up<-cbind(d_up,a,f_up,e_up)
colnames(peaks_bed_up)<-c("chr","start","end","peak_id","score","strand")
dim(peaks_bed_up)
setwd("~/xjj/drug/drug_result/HDAC20220524/2_ATACseq")
write.table(peaks_bed_up,file="all_peaks.bed",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')


#######################################################################################
####  homer注释完的信息,画饼图
rm(list=ls())
setwd("~/xjj/drug/drug_result/chemotherapy/2_ATACseq/annotation")
library(data.table)
homer_anno<-fread("sen-res-chemotherapy_AUC-0.01-P-down-annotation.txt",header=T,data.table=F)
length(unique(homer_anno$`Gene Name`))
colnames(homer_anno)<-c("peakid",colnames(homer_anno)[-1])

library(stringr)
Annotation<-apply(as.data.frame(homer_anno$Annotation),1,function(x) str_replace_all(as.character(x),"\\((.*?)\\)",""))##取出括号前的字符
Annotation1<-as.data.frame(Annotation)
Annotation2<-data.frame(lapply(strsplit(as.character(Annotation1[,1]),'\\.'), function(x) x[1])%>%unlist())
position<-as.character(unique(Annotation2[,1]))
result<-data.frame()
for(i in 1:length(position)){
  result[i,1]<-position[i]
  result[i,2]<-sum(Annotation2[,1] %in% position[i])
  result[i,3]<-(sum(Annotation2[,1] %in% position[i])/(length(Annotation2[,1])))*100
}
colnames(result)<-c("Type","Count","Ratio")
library(ggplot2) # 加载包
dt = data.frame(A = result[,3],B = result[,1])#建立数据框
dt = dt[order(dt$A, decreasing = TRUE),]   ## 用 order() 让数据框的数据按 A 列数据从大到小排序
myLabel = as.vector(dt$B)   ## 转成向量，否则图例的标签可能与实际顺序不一致
myLabel = paste(myLabel, "(", round(dt$A / 1, 2), "%)", sep = "")   ## 用 round() 对结果保留两位小数
p = ggplot(dt, aes(x = "", y = A,fill = B)) + #创建坐标轴，fill = B
  geom_bar(stat = "identity") + 
  geom_bar(stat = "identity", width = 1) +   #当width >= 1 时中心的杂点将消失
  coord_polar(theta = "y") +  # 把柱状图折叠成饼图（极坐标）
  labs(x = "", y = "", title = "") +  # 将横纵坐标的标签设为空
  theme(axis.ticks = element_blank()) +  # 将左上角边框的刻度去掉
  theme(legend.title = element_blank(), legend.position = "left")+   ## 将图例标题设为空，并把图例方放在左边位置
  scale_fill_discrete(breaks = dt$B, labels = myLabel)+   # 将原来的图例标签换成现在的myLabel
  theme(axis.text.x = element_blank())+   ## 去掉饼图的外框上的数值，即去除原柱状图的X轴，把X轴的刻度文字去掉
  labs(title = "all-peaks")+
  theme(plot.title = element_text(hjust = 0.5))
#geom_text(aes(y = A/2 + c(0, cumsum(A)[-length(A)]), x = sum(A)/20, label = myLabel), size = 5)   # 在图中加上百分比：x 调节标签到圆心的距离, y 调节标签的左右位置
print(p) #显示饼图
setwd("~/xjj/drug/drug_result/chemotherapy/2_ATACseq/annotation")
ggsave("all_annotation.pdf",p,width = 6, height = 5)
write.table(result,file="down_P0.01_FC0_annotation_ratio.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)

########   堆积图   ###################
rm(list=ls())
library(reshape2)
setwd("~/xjj/drug/drug_result/chemotherapy/2_ATACseq/annotation")
all<-read.delim("all_annotation_ratio.txt",sep = '\t',header = T)
down<-read.delim("down_P0.01_FC0_annotation_ratio.txt",sep = '\t',header = T)
up<-read.delim("up_P0.01_FC0_annotation_ratio.txt",sep = '\t',header = T)
all1<-all[order(all[,1]),]
down1<-down[order(down[,1]),]
up1<-up[order(up[,1]),]
dat_m1<-rbind(all1[,c(1,3)],down1[,c(1,3)],up1[,c(1,3)])
m1<-matrix(c(rep("all",8),rep("down",8),rep("up",8)),ncol=1)
dat_m<-cbind(m1,dat_m1)
colnames(dat_m) = c('Group','Type','value')
#定义`Group`列的出图顺序
dat_m$Group = factor(dat_m$Group, levels = c("all","down","up"))
library(RColorBrewer)
display.brewer.all()

library(ggplot2)
library(scales)
p<-ggplot(data = dat_m, aes(x = Group, y = value, fill = Type)) + 
  theme_bw()+
  geom_bar(stat= 'identity', width = 0.7,position="fill")+ #堆叠图，position = fill 表示堆叠图
  labs(x = 'Group',y = 'Ratio',title = "Genomic location of DARs annotation",vjust = 0.5)+ #定义坐标轴以及图例标题
  scale_fill_brewer(palette = 'RdBu') +#自定义颜色，可通过`library(RColorBrewer);display.brewer.all()`来展示所有备选项
  #scale_y_continuous(labels = percent) +  ## 百分比坐标轴（需加载scales包）
  guides(fill = guide_legend(ncol = 1, bycol = TRUE, override.aes = list(size = 5))) +#定义图例的布局，1列，排序，图例中色块的大小增大5倍
  theme(axis.title.y = element_text(color = 'black',size = 12),
        axis.title.x = element_text(color = 'black',size = 12,vjust = -1.2),
        axis.text.y = element_text(color = 'black',size = 10),
        axis.text.x = element_text(color = 'black',size = 12), #x轴标签偏转45°，并下降0.5
        panel.grid = element_blank(),
        legend.position = 'right',
        legend.key.height = unit(0.6,'cm'),#定义图例中色块的高度
        legend.text = element_text(color = 'black',size = 10))
setwd("~/xjj/drug/drug_result/chemotherapy/2_ATACseq/annotation")
ggsave(p, filename = 'annotation-ratio-all-down-up.pdf', width = 8, height = 6, dpi = 600)
###############  peaks都处于什么位置  ######################    
rm(list=ls())
setwd("~/xjj/drug/drug_result/chemotherapy/2_ATACseq/annotation")
library(data.table)
homer_anno1<-fread("all_peaks-annotation.txt",header=T,data.table=F)

setwd("~/xjj/drug/drug_result/chemotherapy/2_ATACseq")
res<-read.table("sensitivity-resistance-all-DESeq2_chemotherapy_AUC_ATAC.txt",header = TRUE,sep = "\t")
resSig<-res[which(res$pvalue<0.01),]
int_peak<-intersect(resSig[,1],homer_anno1[,1])

resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
up_gene<-as.data.frame(resSig[which(resSig$log2FoldChange>0),1])
down_gene<-as.data.frame(resSig[which(resSig$log2FoldChange<0),1])
int_peak<-intersect(up_gene[,1],homer_anno1[,1])


homer_anno<-homer_anno1[match(int_peak,homer_anno1[,1]),]
one<-sum(abs(homer_anno$`Distance to TSS`)<=100)
two<-length(which(abs(homer_anno$`Distance to TSS`)>=100 & abs(homer_anno$`Distance to TSS`)<=1000))
three<-length(which(abs(homer_anno$`Distance to TSS`)>=1000 & abs(homer_anno$`Distance to TSS`)<=10000))
four<-length(which(abs(homer_anno$`Distance to TSS`)>=10000 & abs(homer_anno$`Distance to TSS`)<=100000))
five<-sum(abs(homer_anno$`Distance to TSS`)>=100000)
ma<-matrix(rep("allpeaks-peaks",5),ncol=1)
value<-matrix(c(one,two,three,four,five),ncol=1)
type<-matrix(c("0-100","100-1000","1000-10000","10000-100000",">100000"),ncol=1)
data<-cbind(ma,type,value)
colnames(data)<-c("sample","type","value")

library(ggplot2) # 加载包
dt = data.frame(A = (as.numeric(data[,3])/sum(as.numeric(data[,3])))*100,B = data[,2])#建立数据框
myLabel = as.vector(dt$B)   ## 转成向量，否则图例的标签可能与实际顺序不一致
myLabel = paste(myLabel, "(", round(dt$A / 1, 2), "%)", sep = "")   ## 用 round() 对结果保留两位小数

p = ggplot(dt, aes(x = "", y = A, fill = B)) + #创建坐标轴
  geom_bar(stat = "identity") + 
  geom_bar(stat = "identity", width = 1) +   #当width >= 1 时中心的杂点将消失
  coord_polar(theta = "y") +  # 把柱状图折叠成饼图（极坐标）
  labs(x = "", y = "", title = "downpeaks-P0.01-position") +  # 将横纵坐标的标签设为空
  theme(axis.ticks = element_blank()) +  # 将左上角边框的刻度去掉
  theme(legend.title = element_blank(), legend.position = "left")+   ## 将图例标题设为空，并把图例方放在左边位置
  scale_fill_discrete(breaks = dt$B, labels = myLabel)+   # 将原来的图例标签换成现在的myLabel
  theme(axis.text.x = element_blank())+   ## 去掉饼图的外框上的数值，即去除原柱状图的X轴，把X轴的刻度文字去掉
  theme(plot.title = element_text(hjust = 0.5))
#geom_text(aes(y = A/2 + c(0, cumsum(A)[-length(A)]), x = sum(A)/20, label = myLabel), size = 5)   # 在图中加上百分比：x 调节标签到圆心的距离, y 调节标签的左右位置
print(p) #显示饼图
setwd("~/xjj/drug/drug_result/chemotherapy/2_ATACseq/annotation")
ggsave("uppeaks-P0.01-position.pdf",p,width = 6, height = 5)
write.table(cbind(data,dt),file="allpeaks_position.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')

###### 堆积图
rm(list=ls())
library(reshape2)
setwd("~/xjj/drug/drug_result/chemotherapy/2_ATACseq/annotation")
all<-read.delim("allpeaks_position.txt",sep = '\t',header = T)
#DA<-read.delim("DApeaks_P0.01_FC0_position.txt",sep = '\t',header = T)
down<-read.delim("downpeaks_P0.01_FC0_position.txt",sep = '\t',header = T)
up<-read.delim("uppeaks_P0.01_FC0_position.txt",sep = '\t',header = T)
dat_m<-rbind(all[,c(1,2,4)],down[,c(1,2,4)],up[,c(1,2,4)])
colnames(dat_m) = c('Group','Type','value')
#定义`Group`列的出图顺序
dat_m$Group = factor(dat_m$Group, levels = c("allpeaks-peaks","down-P0.01-peaks","up-P0.01-peaks"))
dat_m$Type = factor(dat_m$Type, levels = unique(dat_m$Type))

library(RColorBrewer)
display.brewer.all()

library(ggplot2)
library(scales)
p<-ggplot(data = dat_m, aes(x = Group, y = value, fill = Type)) + 
  theme_bw()+
  geom_bar(stat= 'identity', width = 0.5,position="fill")+ #堆叠图，position = fill 表示堆叠图
  labs(x = '',y = 'Percentage',title = "Distance to the closest transcription start site (TSS)",vjust = 0.5)+ #定义坐标轴以及图例标题
  scale_fill_brewer(palette = 'RdBu') +#自定义颜色，可通过`library(RColorBrewer);display.brewer.all()`来展示所有备选项
  #scale_y_continuous(labels = percent) +  ## 百分比坐标轴（需加载scales包）
  guides(fill = guide_legend(ncol = 1, bycol = TRUE, override.aes = list(size = 5))) +#定义图例的布局，1列，排序，图例中色块的大小增大5倍
  theme(axis.title.y = element_text(color = 'black',size = 12),
        axis.title.x = element_text(color = 'black',size = 12,vjust = -1.2),
        axis.text.y = element_text(color = 'black',size = 10),
        axis.text.x = element_text(color = 'black',size = 12,angle = 30,vjust = 0.5), #x轴标签偏转45°，并下降0.5
        panel.grid = element_blank(),
        legend.position = 'right',
        legend.key.height = unit(0.6,'cm'),#定义图例中色块的高度
        legend.text = element_text(color = 'black',size = 10))
setwd("~/xjj/drug/drug_result/chemotherapy/2_ATACseq/annotation")
ggsave(p, filename = 'Distance to the closest transcription start site.pdf', width = 6, height = 5, dpi = 600)

#####################  符合在TSS100KB之内条件的注释基因
rm(list=ls())
setwd("~/xjj/drug/drug_result/chemotherapy/2_ATACseq/annotation")
library(data.table)
homer_anno_down<-fread("sen-res-chemotherapy_AUC-0.01-P-logFC0-down-annotation.txt",header=T,data.table=F)
colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]
length(unique(Anno_gene_100Kb_down$`Gene Name`))

homer_anno_up<-fread("sen-res-chemotherapy_AUC-0.01-P-logFC0-up-annotation.txt",header=T,data.table=F)
colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
length(unique(Anno_gene_100Kb_up$`Gene Name`))

setwd("~/xjj/drug/drug_result/chemotherapy/1_RNAseq")
up_gene<-read.table("sensitivity-resistance-all-DESeq2_AUC_RNAseq_P0.05_FC0_up_DEGs.txt",header=T,sep="\t")
down_gene<-read.table("sensitivity-resistance-all-DESeq2_AUC_RNAseq_P0.05_FC0_down_DEGs.txt",header=T,sep="\t")
all<-rbind(up_gene,down_gene)
DEGs<-as.character(all[,1])
down_DEGs_DApeaks<-intersect(Anno_gene_100Kb_down$`Gene Name`,DEGs)
up_DEGs_DApeaks<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)
cat(length(down_DEGs_DApeaks),"\n")
cat(length(up_DEGs_DApeaks),"\n")

library(gplots)
library(VennDiagram)
down_peaks_gene <- Anno_gene_100Kb_down$`Gene Name`
up_peaks_gene <- Anno_gene_100Kb_up$`Gene Name`
input  <-list(unique(down_peaks_gene),unique(up_peaks_gene),unique(DEGs))
venn(input,showSetLogicLabel=TRUE)
tmp <- venn(input)
int<-attr(tmp, "intersections")

setwd("~/xjj/drug/drug_result/chemotherapy/1_2_intersect/venn")
library(VennDiagram)
venn.diagram(x=list(downpeaks=unique(down_peaks_gene),uppeaks=unique(up_peaks_gene),RNAseq=unique(DEGs)), "DEGP005_DApeakP001.png",fill=c("red","green","blue"),margin = 0.1)

down_peaks_gene1<-as.data.frame(unique(down_peaks_gene))
colnames(down_peaks_gene1)<-"down_peaks_gene"
up_peaks_gene1<-as.data.frame(unique(up_peaks_gene))
colnames(up_peaks_gene1)<-"up_peaks_gene"
write.table(down_peaks_gene1,"down_peaks0.01_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(up_peaks_gene1,"up_peaks0.01_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

######## 差异上调的peaks与DEGs对应的基因在两类样本中的mRNA水平
###真的是上调的基因
####  DESeq2
up1<-int[["B:C"]]
up2<-int[["A:B:C"]]
up3<-c(up1,up2)
int_up<-as.data.frame(intersect(up3,up_gene$gene_id))
colnames(int_up)<-"int_up"
int_up_peak<-as.data.frame(homer_anno_up[homer_anno_up$`Gene Name`%in%int_up[,1],1])
colnames(int_up_peak)<-"int_up_peak"

down1<-int[["A:C"]]
down2<-int[["A:B:C"]]
down3<-c(down1,down2)
int_down<-as.data.frame(intersect(down3,down_gene$gene_id))
colnames(int_down)<-"int_down"
int_down_peak<-as.data.frame(homer_anno_down[homer_anno_down$`Gene Name`%in%int_down[,1],1])
colnames(int_down_peak)<-"int_down_peak"

down_DEGs_down_DApeaks<-as.data.frame(intersect(Anno_gene_100Kb_down$`Gene Name`,down_gene[,1]))
colnames(down_DEGs_down_DApeaks)<-"down_DEGs_down_DApeaks"
up_DEGs_up_DApeaks<-as.data.frame(intersect(Anno_gene_100Kb_up$`Gene Name`,up_gene[,1]))
colnames(up_DEGs_up_DApeaks)<-"up_DEGs_up_DApeaks"

setwd("~/xjj/drug/drug_result/chemotherapy/1_2_intersect/venn")
write.table(int_up,"int_up_peak0.01_all_DEG0.05_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(int_down,"int_up_peak0.01_all_DEG0.05_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(int_up_peak,"int_up_peak0.01_all_DEG0.05_peak.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(int_down_peak,"int_down_peak0.01_all_DEG0.05_peak.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(up_DEGs_up_DApeaks,"int_up_peak0.01_up_DEG0.05_peak.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(down_DEGs_down_DApeaks,"int_down_peak0.01_down_DEG0.05_peak.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出


##################################################################################
##############  DEpeaks在不同样本变化的倍数与DEGs在不同样本变化的倍数之间的相关性
rm(list=ls())
setwd("~/xjj/drug/drug_result/chemotherapy/2_ATACseq/annotation")
library(data.table)
homer_anno_down<-fread("sen-res-chemotherapy_AUC-0.01-P-logFC0-down-annotation.txt",header=T,data.table=F)
colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]
length(unique(Anno_gene_100Kb_down$`Gene Name`))

homer_anno_up<-fread("sen-res-chemotherapy_AUC-0.01-P-logFC0-up-annotation.txt",header=T,data.table=F)
colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
length(unique(Anno_gene_100Kb_up$`Gene Name`))

setwd("~/xjj/drug/drug_result/chemotherapy/1_RNAseq")
up_gene<-read.table("sensitivity-resistance-all-DESeq2_AUC_RNAseq_P0.05_FC0_up_DEGs.txt",header=T,sep="\t")
down_gene<-read.table("sensitivity-resistance-all-DESeq2_AUC_RNAseq_P0.05_FC0_down_DEGs.txt",header=T,sep="\t")
all<-rbind(up_gene,down_gene)
DEGs<-as.character(all[,1])
down_DEGs_DApeaks<-intersect(Anno_gene_100Kb_down$`Gene Name`,DEGs)
up_DEGs_DApeaks<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)

int_all_DEGs<-unique(c(down_DEGs_DApeaks,up_DEGs_DApeaks))
all_peaks_anno<-rbind(homer_anno_down,homer_anno_up)
#DEGs_peaks<-all_peaks_anno[all_peaks_anno$`Gene Name` %in% int_all_DEGs,]
DEGs_peaks<-all_peaks_anno[match(int_all_DEGs,all_peaks_anno$`Gene Name`),1] #match只取出对应上的第一个peak

DEGs_peaks1<-all_peaks_anno[all_peaks_anno$`Gene Name` %in% int_all_DEGs,]
setwd("~/xjj/drug/drug_result/chemotherapy/2_ATACseq")
DApeak<-fread("sensitivity-resistance-all-DESeq2_chemotherapy_AUC_ATAC.txt",header=T,data.table=F)
DEGs_peaks_FC<-DApeak[match(DEGs_peaks,DApeak$peak_id),3]
DEGs_FC<-all[match(int_all_DEGs,DEGs),3]

#pearson", "kendall", "spearman
cor_result<-cor.test(DEGs_peaks_FC, DEGs_FC,alternative = "two.sided",method = "pearson")
cor_result$estimate
cor_result$p.value
################################   画相关性分析图
gene<-as.character(all[match(int_all_DEGs,DEGs),1])
a1<-matrix(DEGs_FC,ncol=1)
a2<-matrix(DEGs_peaks_FC,ncol=1)
dat<-cbind(a1,a2)
dat1<-as.data.frame(dat)
rownames(dat1)<-gene

colnames(dat1)<-c("DEGslog2FC","DApeakslog2FC")
setwd("~/xjj/drug/drug_result/chemotherapy/1_2_intersect/correlation")
write.table(DEG_DApeak_FC,"DEGs_FC_DEGs_peaks_FC_correlation_comprehensive.txt",sep = '\t',col.names = T,row.names = T,quote = FALSE)##数据输出

library(ggplot2)
library(ggpubr)
a11<-c(1:10)
a21<-c((nrow(dat1)-9):nrow(dat1))
a3<-c(a11,a21)
want_dat<-dat1[order(dat1[,1])[a3],]
ps<-ggplot(data=dat1, aes(x=DEGslog2FC, y=DApeakslog2FC))+geom_point(color="red",size=0.1)+stat_smooth(method="lm",se=FALSE)+stat_cor(data=dat1, method = "pearson")+
  ggtitle("pearson-DApeaksP0.01-DEGsP0.05") +
  theme(plot.title = element_text(hjust = 0.5))+
  geom_text_repel(
    data = want_dat[,c(1:2)],
    aes(label = different_FC[order(dat1[,1])[a3],1]),
    size = 3,
    color = "black",
    segment.color = "black", show.legend = FALSE )
setwd("~/xjj/drug/drug_result/chemotherapy/1_2_intersect/correlation")
ggsave("correlation-DApeaksP0.01-DEGsP0.05-pearson-comprehensive.pdf",ps,width = 10, height = 10)

#####################################################################
rm(list=ls())
setwd("~/xjj/drug/drug_result/chemotherapy/2_ATACseq/annotation")
library(data.table)
homer_anno_down<-fread("sen-res-chemotherapy_AUC-0.01-P-logFC0-down-annotation.txt",header=T,data.table=F)
colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]
length(unique(Anno_gene_100Kb_down$`Gene Name`))

homer_anno_up<-fread("sen-res-chemotherapy_AUC-0.01-P-logFC0-up-annotation.txt",header=T,data.table=F)
colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
length(unique(Anno_gene_100Kb_up$`Gene Name`))

setwd("~/xjj/drug/drug_result/chemotherapy/1_RNAseq")
up_gene<-read.table("sensitivity-resistance-all-DESeq2_AUC_RNAseq_P0.05_FC0_up_DEGs.txt",header=T,sep="\t")
down_gene<-read.table("sensitivity-resistance-all-DESeq2_AUC_RNAseq_P0.05_FC0_down_DEGs.txt",header=T,sep="\t")
all<-rbind(up_gene,down_gene)
DEGs<-as.character(all[,1])
down_DEGs_DApeaks<-intersect(Anno_gene_100Kb_down$`Gene Name`,DEGs)
up_DEGs_DApeaks<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)

int_all_DEGs<-unique(c(down_DEGs_DApeaks,up_DEGs_DApeaks))
all_peaks_anno<-rbind(homer_anno_down,homer_anno_up)
setwd("~/xjj/drug/drug_result/chemotherapy/2_ATACseq")
DApeak<-fread("sensitivity-resistance-all-DESeq2_chemotherapy_AUC_ATAC.txt",header=T,data.table=F)


different_FC<-NULL
for(i in 1:length(int_all_DEGs)){
  DEGs_peaks1<-all_peaks_anno[all_peaks_anno$`Gene Name` %in% int_all_DEGs[i],1]
  DEGs_peaks_FC<-as.data.frame(DApeak[match(DEGs_peaks1,DApeak$peak_id),3])
  gene_name<-matrix(rep(int_all_DEGs[i],length(DEGs_peaks1)),ncol=1)
  DEGs_FC<-matrix(rep(all[match(int_all_DEGs[i],DEGs),3],ncol=1))
  different_FC1<-cbind(gene_name,DEGs_FC,DEGs_peaks_FC)
  different_FC<-rbind(different_FC,different_FC1)
}
DEG_DApeak_FC<-different_FC[,2:3]
colnames(DEG_DApeak_FC)<-c("DEGslog2FC","DApeakslog2FC")
rownames(DEG_DApeak_FC)<-different_FC[,1]
setwd("~/xjj/drug/drug_result/chemotherapy/1_2_intersect/correlation")
write.table(DEG_DApeak_FC,"DEGs_FC_DEGs_peaks_FC_correlation_comprehensive.txt",sep = '\t',col.names = T,row.names = T,quote = FALSE)##数据输出

#pearson", "kendall", "spearman
cor_result<-cor.test(DEG_DApeak_FC[,1],DEG_DApeak_FC[,2],alternative = "two.sided",method = "spearman")
cor_result$estimate
cor_result$p.value

################################   画相关性分析图
dat1<-as.data.frame(DEG_DApeak_FC)
library(ggplot2)
library(ggpubr)
a11<-c(1:10)
a21<-c((nrow(dat1)-9):nrow(dat1))
a3<-c(a11,a21)
want_dat<-dat1[order(dat1[,1])[a3],]
ps<-ggplot(data=dat1, aes(x=DEGslog2FC, y=DApeakslog2FC))+geom_point(color="red",size=0.1)+stat_smooth(method="lm",se=FALSE)+stat_cor(data=dat1, method = "spearman")+
  ggtitle("spearman-DApeaksP0.01-DEGsP0.05") +
  theme(plot.title = element_text(hjust = 0.5))+
  geom_text_repel(
    data = want_dat[,c(1:2)],
    aes(label = different_FC[order(dat1[,1])[a3],1]),
    size = 3,
    color = "black",
    segment.color = "black", show.legend = FALSE )
setwd("~/xjj/drug/drug_result/chemotherapy/1_2_intersect/correlation")
ggsave("correlation-DApeaksP0.01-DEGsP0.05-spearman-comprehensive.pdf",ps,width = 10, height = 10)

#pearson", "kendall", "spearman
##画好看的图ggstatsplot

#### 使用超几何分布看看交叠是否随机出现
#N： 总样本数
#m： 总样本中“特定类别”的数量
#n: 从总样本中随机抽取的数量
#k: 抽取样本中属于“特定类别”的数量
# P值计算公式
1-phyper(k-1,m, N-m, n)
1-phyper(212, 746, 25421, 317) #up P=0
1-phyper(346, 1456, 24711, 394) #down P=0
# 用phyper(k-1,M, N-M, n, lower.tail=F)代替 1-phyper(k-1,m, N-m, n)
phyper(k-1,M, N-M, n, lower.tail=F)
# up k=155(交叠的DEGs),m=3718 (上调peaks中注释到的基因，100KB之内)，n= 601 (upDEGs的个数)
# N 在RNAseq中所有的基因,nrow(DEGs)=26167
setwd("~/xjj/drug/drug_result/chemotherapy/1_RNAseq")
DEGs<-read.table("sensitivity-resistance-all-DESeq2_chemotherapy_AUC_RNAseq_DEGs.txt",header=T,sep="\t")

phyper(161,3718, 22449, 653, lower.tail=F) #up P=2.680927e-13
phyper(227,3504, 22663, 916, lower.tail=F) #down P=1.277207e-21
