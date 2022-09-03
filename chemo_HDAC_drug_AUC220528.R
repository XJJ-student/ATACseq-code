#散点图
rm(list=ls())
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
setwd("~/xjj/drug")
ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,!(colnames(ic501)%in%c("PC.100","PC.34"))]## all of PDAC
colnames(ic50)<-sample_id[match(colnames(ic50),sample_id[,2]),1]
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
drug_info<-read.csv("compiled_drug_info.csv",header = TRUE,sep = ",")
HDAC1<-drug_info[drug_info$target=="HDAC",]
aaa<-c(as.character(HDAC1[,1]),"S1848","S2759","S1194","S1047","5-FU","GEM","IRI","OXA","PAC")
want_AUC<-ic50[rownames(ic50)%in%aaa,]
label<-c(rep("chemo",5),rep("HDAC",20))
#####每个样本在一类药物中最敏感的AUC值
AUC_chemo_and_HDAC<-apply(want_AUC,2,max)
result<-matrix(0,length(AUC_chemo_and_HDAC),3)
for(i in 1:length(AUC_chemo_and_HDAC)){
  result[i,1]<-colnames(want_AUC)[i]
  result[i,2]<-rownames(want_AUC)[want_AUC[,i]%in%AUC_chemo_and_HDAC[i]]
  result[i,3]<-label[want_AUC[,i]%in%AUC_chemo_and_HDAC[i]]
}

ccc<-c(as.character(HDAC1[,1]),"S1848","S2759","S1194","S1047")##去掉S1848，因为它在两类样本中AUC差异不显著
bbb<-c("5-FU","GEM","IRI","OXA","PAC")  ##5种化疗药
want_AUC_chemo<-ic50[rownames(ic50)%in%bbb,]
want_AUC_HDAC<-ic50[rownames(ic50)%in%ccc,]
AUC_chemo<-apply(want_AUC_chemo,2,max)
AUC_HDAC<-apply(want_AUC_HDAC,2,max)
abs(AUC_chemo-AUC_HDAC)


library(pheatmap)
mark_result1<-matrix(0,5,37)
  for(i in 1:length(colnames(want_AUC_chemo))){
    a<-which(want_AUC_chemo[,i]==AUC_chemo[i])
    mark_result1[a,i]<-5
  }
mark_result2<-matrix(0,20,37)
for(j in 1:length(colnames(want_AUC_HDAC))){
  b<-which(want_AUC_HDAC[,j]==AUC_HDAC[j])
  mark_result2[b,j]<-5
}
mark_result<-rbind(mark_result1,mark_result2)
  
for(i in 1:25){
  for(j in 1:37){
    mark_result[i,j]=ifelse(mark_result[i,j]>0, "*", "")# 用"*"代替>0的相对丰度，用""代替<=0.1的相对丰度
  }
}
pheatmap(want_AUC,scale="none",cluster_cols=FALSE, cluster_rows=FALSE,gaps_row = 5,display_numbers=mark_result)

library(pheatmap)
mark_result<-matrix(0,25,37)
for(i in 1:length(colnames(want_AUC))){
  a<-which(rownames(want_AUC)==result[which(result[,1]==colnames(want_AUC)[i]),2])
  mark_result[a,i]<-5
}

for(i in 1:25){
  for(j in 1:37){
    mark_result[i,j]=ifelse(mark_result[i,j]>0, "*", "")# 用"*"代替>0的相对丰度，用""代替<=0.1的相对丰度
  }
}
pheatmap(want_AUC,scale="none",cluster_cols=FALSE, cluster_rows=FALSE,gaps_row = 5,display_numbers=mark_result)



geneid<-colnames(want_AUC_chemo)
tresult=matrix(0,length(geneid),4)
for (i in 1:length(geneid)){
  ttest=t.test(want_AUC_HDAC[,i],want_AUC_chemo[,i],alternative = "two.sided")
  tresult[i,1:3]=c(geneid[i],ttest$statistic,ttest$p.value)
}
tresult[,4]=p.adjust(tresult[,3],method="BH")
colnames(tresult)<-c("geneid","statistic","p.value","FDR")
which(tresult[,4]<0.05)
length(which(tresult[,3]<0.05))

####  数量严重不均衡，没法用
geneid<-colnames(want_AUC_chemo)
tresult=matrix(0,length(geneid),4)
for (i in 1:length(geneid)){
  ttest=wilcox.test(as.numeric(want_AUC_HDAC[i,]),as.numeric(want_AUC_chemo[i,]),alternative ="two.sided")
  tresult[i,1:3]=c(geneid[i],ttest$statistic,ttest$p.value)
}
tresult[,4]=p.adjust(tresult[,3],method="BH")
colnames(tresult)<-c("geneid","statistic","p.value","FDR")
which(tresult[,4]<0.05)


#######################################################################################
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
want_AUC_chemo<-ic50[rownames(ic50)%in%"GEM",]
want_AUC_HDAC<-ic50[rownames(ic50)%in%"S8495",]

#####每个样本在一类药物中最敏感的AUC值
AUC_chemo<-apply(want_AUC_chemo,2,max)
AUC_HDAC<-apply(want_AUC_HDAC,2,max)
data <- data.frame(x=AUC_chemo,y=AUC_HDAC)
data <- data.frame(x=t(want_AUC_chemo),y=t(want_AUC_HDAC))
#作散点图
freq<-0.5
ggplot(data, aes(x=GEM, y=S8495)) + 
  geom_point()+
  geom_hline(yintercept= freq)+
  geom_vline(xintercept=freq)+
  xlab("chemotherapeutic") +
  ylab("HDACIs")


want<-rbind(data[data$y>freq&data$x<freq,],data[data$y<freq&data$x>freq,])
ggplot(data, aes(x=x, y=y)) + 
  geom_point()+
  geom_hline(yintercept= freq)+
  geom_vline(xintercept=freq)+
  #geom_point(data = data[data$y>freq&data$x<freq,],size = 4)+
  #geom_point(data = data[data$y<freq&data$x>freq,],size = 4)+
  xlab("chemotherapeutic") +
  ylab("HDACIs")+
  geom_text_repel(
    data = want,
    aes(label = rownames(want)),
    size = 3,
    color = "black",
    segment.color = "black", show.legend = FALSE )

HDAC_sen<-data[data$y>freq&data$x<freq,]
colnames(HDAC_sen)<-c("AUC_chemo","AUC_HDAC")

chemo_sen<-data[data$y<freq&data$x>freq,]
colnames(chemo_sen)<-c("AUC_chemo","AUC_HDAC")

all_sen<-data[data$y>freq&data$x>freq,]
colnames(all_sen)<-c("AUC_chemo","AUC_HDAC")

all_res<-data[data$y<freq&data$x<freq,]
colnames(all_res)<-c("AUC_chemo","AUC_HDAC")


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

####################################配对
# I:HDAC_sen II:chemo_sen III:all_sen IV:all_res
#### I II+III+IV
HDAChigh_chemolow<-rownames(HDAC_sen)
other1<-c(rownames(chemo_sen),rownames(all_sen),rownames(all_res))
#### II I+III+IV
chemohigh_HDAClow<-rownames(chemo_sen)
other2<-c(rownames(HDAC_sen),rownames(all_sen),rownames(all_res))

#### III I+II+IV
allhigh<-rownames(all_sen)
other3<-c(rownames(HDAC_sen),rownames(chemo_sen),rownames(all_res))

#### IV I+II+III
alllow<-rownames(all_res)
other4<-c(rownames(HDAC_sen),rownames(chemo_sen),rownames(all_sen))

HDAChigh_chemolow<-c("DAC-23","DAC-26","DAC-31","DAC-35","DAC-7","DAC-8","DAC-2","DAC-14")
chemohigh_HDAClow<-c("DAC-33","DAC-34","DAC-5","DAC-11","DAC-18","DAC-37","DAC-39")
allhigh<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-36","DAC-6","DAC-1","DAC-9","DAC-12","DAC-15","DAC-3","DAC-16","DAC-19","DAC-38")
alllow<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-13")


want<-rbind(want_AUC_chemo,want_AUC_HDAC)
sensitivity_AUC<-want[,match(chemohigh_HDAClow,colnames(want))]
resistance_AUC<-want[,match(other2,colnames(want))]
geneid<-rownames(want)
tresult=matrix(0,length(geneid),4)
for (i in 1:length(geneid)){
  ttest=t.test(sensitivity_AUC[i,],resistance_AUC[i,],alternative = "two.sided")
  tresult[i,1:3]=c(geneid[i],ttest$statistic,ttest$p.value)
}
tresult[,4]=p.adjust(tresult[,3],method="BH")
colnames(tresult)<-c("geneid","statistic","p.value","FDR")
which(tresult[,4]<0.05)
length(which(tresult[,3]<0.05))

############################################################
############    RNAseq DEGs
rm(list=ls())
setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
library("data.table")
exp1<-fread("star_rsem.GeneSymbol.readscount.xls",header=T,data.table=F)
exp2<-floor(exp1[,-1])
delete<-apply(exp2,1,function(x) mean(x==0))
cou<-which(delete>0.75)####在超过75%样本以上的都是0的基因删去
exp<-exp2[(-cou),]

resistance<-other4
sensitivity<-alllow

resistance_exp<-exp[,match(resistance,colnames(exp))]
sensitivity_exp<-exp[,match(sensitivity,colnames(exp))]
geneid<-exp1[(-cou),1]
library(DESeq2)
colDate<-data.frame(row.names = c(as.vector(resistance),as.vector(sensitivity)),
                    condition=factor(c(rep("resistance",length(resistance)),rep("sensitivity",length(sensitivity))))
)
datexpr<-cbind(resistance_exp,sensitivity_exp)
counts <- apply(datexpr,2,as.numeric)###矩阵中必须是数值
dds<-DESeqDataSetFromMatrix(countData = counts,colData = colDate,design = ~condition)
dds##这是一个关于各种内容的矩阵
dds<-DESeq(dds)##进行标准化分析
sizeFactors(dds)##查看每个主成分的标准化值
res<-results(dds)##将结果输出
head(res)
class(res)##可以看出它是DESeq的属性，要转化为表格
res<-as.data.frame(res)
head(res)
res<-cbind(geneid,res)  ##对数据增加一列
head(res)
colnames(res)<- c('gene_id',"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj" )##对列名进行修改
setwd("~/xjj/drug/drug_result/HDACi_chemo/alllow/1_RNAseq")
write.table(res,"alllow-other4-all-DESeq2_AUC_RNAseq_DEGs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
#resSig<-res[which(res$pvalue<0.05 & abs(res$log2FoldChange>1)),]
#resSig<-res[which(res$padj<0.05),]
resSig<-res[which(res$pvalue<0.05),]

#pvalue padj pvalue<0.05 1441个DEGs
##新增一列，将log2FoldChange>0标注为up，<0标准为down
resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
sum(resSig$up_down=='up')
sum(resSig$up_down=='down')
up_gene<-resSig[which(resSig$log2FoldChange>0),]
down_gene<-resSig[which(resSig$log2FoldChange<0),]

resistance_exp1 <- apply(resistance_exp,2,as.numeric)
sensitivity_exp1 <- apply(sensitivity_exp,2,as.numeric)
resistance_exp11<-apply(resistance_exp1,1,mean)
sensitivity_exp11<-apply(sensitivity_exp1,1,mean)
up1<-match(up_gene[,1],geneid)
down1<-match(down_gene[,1],geneid)
sum((sensitivity_exp11[up1]-resistance_exp11[up1])>0)
sum((sensitivity_exp11[down1]-resistance_exp11[down1])<0)
setwd("~/xjj/drug/drug_result/HDACi_chemo/alllow/1_RNAseq")
write.table(up_gene,"alllow_other4_RNAseq_P0.05_FC0_up_DEGs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(down_gene,"alllow_other4_RNAseq_P0.05_FC0_down_DEGs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

##################################### 转录组区分亚组  转录组的特征不是很完美的将四个亚组分开#################
setwd("~/xjj/drug/drug_result/HDACi_chemo/alllow/1_RNAseq")
alllow<-read.table("alllow_other4_RNAseq_P0.05_FC0_up_DEGs.txt",sep="\t",header=T)
setwd("~/xjj/drug/drug_result/HDACi_chemo/allhigh/1_RNAseq")
allhigh<-read.table("allhigh_other3_RNAseq_P0.05_FC0_up_DEGs.txt",sep="\t",header=T)
setwd("~/xjj/drug/drug_result/HDACi_chemo/HDAChigh_chemolow/1_RNAseq")
HDAChigh_chemolow<-read.table("HDAChigh_chemolow_other1_RNAseq_P0.05_FC0_up_DEGs.txt",sep="\t",header=T)
setwd("~/xjj/drug/drug_result/HDACi_chemo/chemohigh_HDAClow/1_RNAseq")
chemohigh_HDAClow<-read.table("chemohigh_HDAClow_other2_RNAseq_P0.05_FC0_up_DEGs.txt",sep="\t",header=T)

one1<-intersect(as.character(alllow[,1]),as.character(allhigh[,1]))#0
one2<-intersect(as.character(alllow[,1]),as.character(HDAChigh_chemolow[,1]))#3
one3<-intersect(as.character(alllow[,1]),as.character(chemohigh_HDAClow[,1]))#2
two1<-intersect(as.character(allhigh[,1]),as.character(HDAChigh_chemolow[,1]))#0
two2<-intersect(as.character(allhigh[,1]),as.character(chemohigh_HDAClow[,1]))#2
three1<-intersect(as.character(HDAChigh_chemolow[,1]),as.character(chemohigh_HDAClow[,1]))#2

##用这些特征画热图
HDAChigh_chemolow1<-rownames(HDAC_sen)
chemohigh_HDAClow1<-rownames(chemo_sen)
allhigh1<-rownames(all_sen)
alllow1<-rownames(all_res)
sample_order<-c(HDAChigh_chemolow1,chemohigh_HDAClow1,allhigh1,alllow1)

setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
library("data.table")
exp<-fread("star_rsem.GeneSymbol.TPM.xls",header=T,data.table=F)
peaks<-c(as.character(HDAChigh_chemolow[,1]),as.character(chemohigh_HDAClow[,1]),as.character(allhigh[,1]),as.character(alllow[,1]))
peak_RPKM1<-exp[match(peaks,exp[,1]),]
peaks1<-matrix(c(rep("Class1",nrow(HDAChigh_chemolow)),rep("Class2",nrow(chemohigh_HDAClow)),rep("Class3",nrow(allhigh)),rep("Class4",nrow(alllow))),ncol=1)
peak_RPKM_order<-peak_RPKM1[,match(sample_order,colnames(peak_RPKM1))]

col_cut=c(length(HDAChigh_chemolow1),length(HDAChigh_chemolow1)+length(chemohigh_HDAClow1),length(HDAChigh_chemolow1)+length(chemohigh_HDAClow1)+length(allhigh1))
row_cut=c(nrow(HDAChigh_chemolow),nrow(HDAChigh_chemolow)+nrow(chemohigh_HDAClow),nrow(HDAChigh_chemolow)+nrow(chemohigh_HDAClow)+nrow(allhigh))
library(pheatmap)
data0=cbind(peaks1,peak_RPKM_order)

anno1=c(rep("Class1",length(HDAChigh_chemolow1)),rep("Class2",length(chemohigh_HDAClow1)),rep("Class3",length(allhigh1)),rep("Class4",length(alllow1)))
anno1=data.frame(anno1)

rownames(anno1)=as.character(t(colnames(data0))[2:38])
ann_colors = list(
  anno1 = c(Class1 = "#4dbbd5", Class2 = "#00a087",Class3 = "#f39b7f",Class4="#8491b4"), anno=c("Class1"= "#4dbbd5","Class2"= "#00a087","Class3"= "#f39b7f","Class4"="#8491b4"))
H=6
W=5
library(RColorBrewer)
color1=colorRampPalette(c("#3363aa","#1384d5","#23b2ae"))(4)
color10=colorRampPalette(c(color1[1],color1[2]))(35)
color11=colorRampPalette(c(color1[2],color1[4]))(15)
color2=colorRampPalette(c( "#23b2ae","#c5bc5e","#faf513"))(4)
color20=colorRampPalette(c(color2[1],color2[3]))(15)
color21=colorRampPalette(c(color2[3],color2[4]))(35)
mycolor=c(color10,color11,color20,color21)
anno=data.frame(data0$peaks1)
rownames(data0)=as.character(1:dim(data0)[1])
rownames(anno)=rownames(data0)
colnames(anno)="anno"
setwd("~/xjj/drug/drug_result/HDACi_chemo/allpeak")
pdf(paste0("expheatmap.pdf"),height=H,width=W)
pheatmap(data0[2:38],scale="row",color = mycolor,fontsize=7,fontsize_row = 40*H/dim(data0)[1], fontsize_col = 7,cluster_cols=FALSE, cluster_rows=FALSE,annotation_col=anno1,gaps_col = col_cut,gaps_row = row_cut, annotation_row=anno,annotation_colors = ann_colors)
dev.off()


#############################################################################
#########  RNAseq DEGs Volcano Plot
rm(list=ls())
library(ggrepel) 
library(ggplot2)
setwd("~/xjj/drug/drug_result/HDACi_chemo/chemohigh_HDAClow/1_RNAseq")
df<-read.table("chemohigh_HDAClow-other2-all-DESeq2_AUC_RNAseq_DEGs.txt",sep = '\t',header= T,row.names = 1)
#确定是上调还是下调，用于给图中点上色
df$threshold = factor(ifelse(df$pvalue  < 0.05 & abs(df$log2FoldChange) >= 1, ifelse(df$log2FoldChange >= 0 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
df$gene <- row.names(df) #添加一列基因名，以便备注
p<-ggplot(df,aes(x=log2FoldChange,y= -log10(pvalue),color=threshold))+
  geom_point(data = df[df$pvalue<0.05&abs(df$log2FoldChange)>1,],size = 1)+ 
  geom_point(data = df[df$pvalue>0.05|abs(df$log2FoldChange)<1,],size = 1)+
  scale_color_manual(values=c('blue','grey','red'))+#确定点的颜色
  geom_text_repel(
    data = df[df$padj<0.05&abs(df$log2FoldChange)>2,],
    aes(label = gene),
    size = 3,
    color = "black",
    segment.color = "black", show.legend = FALSE )+#添加关注的点的基因名
  ylab('-log10 (pvalue)')+#修改y轴名称
  xlab('log2 (FoldChange)')+#修改x轴名称
  geom_hline(yintercept=-log10(0.05),linetype=4)#添加横线|logFoldChange|>0.25
#geom_vline(xintercept=c(-2,2),linetype=4)#添加竖线padj<0.05
plot(p)

setwd("~/xjj/drug/drug_result/HDACi_chemo/chemohigh_HDAClow/1_RNAseq/Volcano")
DEGs<-df[which(df$padj<0.01 & abs(df$log2FoldChange)>2),]
write.table(DEGs,"chemohigh_HDAClow_FDR01_FC2_DEGs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
ggsave("Volcano_chemohigh_HDAClow_FDR0.01_FC2.pdf",p,width = 8, height = 6)
DEGs<-df[which(df$padj<0.05 & abs(df$log2FoldChange)>2),]
write.table(DEGs,"chemohigh_HDAClow_FDR05_FC2_DEGs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
ggsave("Volcano_chemohigh_HDAClow_FDR0.05_FC2.pdf",p,width = 8, height = 6)

p1<-ggplot(df,aes(x=log2FoldChange,y= -log10(pvalue),color=threshold))+
  geom_point(data = df[df$pvalue<0.05&abs(df$log2FoldChange)>1,],size = 1)+ 
  geom_point(data = df[df$pvalue>0.05|abs(df$log2FoldChange)<1,],size = 1)+
  scale_color_manual(values=c('blue','grey','red'))+#确定点的颜色
  ylab('-log10 (pvalue)')+#修改y轴名称
  xlab('log2 (FoldChange)')+#修改x轴名称
  geom_hline(yintercept=-log10(0.05),linetype=4)#添加横线|logFoldChange|>0.25
#geom_vline(xintercept=c(-2,2),linetype=4)#添加竖线padj<0.05
plot(p1)
ggsave("Volcano_chemohigh_HDAClow_P0.05_FC1.pdf",p1,width = 8, height = 6)



#############  TCGA_PAAD   ############################# 有符合条件的基因
rm(list=ls())
library(org.Hs.eg.db)
library(stringr)
library(clusterProfiler)
library(survival)
library(survminer)
setwd("~/xjj/drug/drug_result/HDACi_chemo/alllow/1_RNAseq")
df<-read.table("alllow-other4-all-DESeq2_AUC_RNAseq_DEGs.txt",sep = '\t',header= T)
#g2<-as.character(df[which(df$pvalue<0.05),1])
g2 = as.character(df[which(df$pvalue<0.05),1])
#g2 =as.character(df[which(df$padj<0.05 & abs(df$log2FoldChange)>2),1])

setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/6_TCGA") #老位置不变
chemotherapy_drug_sample<-read.table("drug_paad1.txt",header=T,sep="\t") #就用这个文件，化疗药物信息
#chemotherapy_drug_sample<-chemotherapy_drug_sample1[which(chemotherapy_drug_sample1$pharmaceutical_therapy_drug_name=="Oxaliplatin"),]
#chemotherapy_drug_sample<-chemotherapy_drug_sample1[which(chemotherapy_drug_sample1$pharmaceutical_therapy_drug_name %in% c("5-FU","5 FU","5FU","5-Fluorouracil","5-fu","5-fluorouracil","5-Fluorouracil?")),]
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
##三等分
g1<-intersect(g2,rownames(exp1))
result<-matrix(0,length(g1),3)
for(i in 1:length(g1)){
  g<-g1[i]
  result[i,1]<-g1[i]
  gene_expp<-as.numeric(exp1[match(g,rownames(exp1)),])
  exp1_colname1<-substring(colnames(exp1),1,12)
  exp1_colname2<-gsub("\\.","-",exp1_colname1)
  up<-exp1_colname2[order(gene_expp)[(ceiling(length(gene_expp)/3)*2):length(gene_expp)]]
  down<-exp1_colname2[order(gene_expp)[1:ceiling(length(gene_expp)/3)]]
  exp1_colname<-c(up,down)
  class_ind<-matrix(0,length(exp1_colname),1)
  class_ind[1:length(up),1]<-"high"   #基因表达更高，越敏感，理论上生存越好
  class_ind[(length(up)+1):length(exp1_colname),1]<-"low" #基因表达更低，越耐药，理论上生存越差
  surv_info<-state_time[match(exp1_colname,state_time[,1]),c(2,4)]#生存时间和生存状态
  label<-as.matrix(class_ind)
  TTP=as.numeric(surv_info[,2])  ##生存时间
  status_TTP=as.matrix(surv_info[,1])  ##生存状态
  status_TTP[which(status_TTP=="Alive")]=0
  status_TTP[which(status_TTP=="Dead")]=1
  status_TTP<-as.numeric(status_TTP)
  surv_info1<-as.data.frame(cbind(TTP,status_TTP))
  
  surv_TTP<-survfit(Surv(TTP, status_TTP) ~ label,data=surv_info1)
  ggsurvplot(surv_TTP,
             pval = TRUE, conf.int = TRUE, #是否显示p值/显示风险区域
             risk.table = TRUE,# 将风险表显示在生存曲线下面
             ggtheme = theme_bw(),
             tables.theme = theme_cleantable(),
             title=paste("    ",g,"up-regulation in sensitivity",sep=" "))
  
  coxp<-coxph(Surv(TTP, status_TTP) ~ label,data=surv_info1)
  result[i,1]<-g1[i]
  result[i,2]<-tidy(coxp)$estimate
  result[i,3]<-tidy(coxp)$p.value
}
colnames(result)<-c("geneid","estimate","p.value")
#fdrresult=p.adjust(result[,3],method="BH")

a<-result[which(result[,3]<0.05),]
sum(result[,3]<0.05)
setwd("~/xjj/drug/drug_result/HDACi_chemo/alllow/1_RNAseq/TCGA_PAAD")
library(ggplot2)
##手动保存
#ggsave(paste(g,"up-regulation in sensitive.pdf",sep=" "),p,width = 5, height = 5)
write.table(a,"alllow_DEGs_P0.05_logrank_P0.05_TCGA_survived.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

##################################################################################
################  DEGs与组蛋白修饰酶的交叠
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDAC_frontiers/5_histone-modifying enzymes")
Acetylase<-read.table("histone-modifying-enzymes.txt",header=T,sep="\t")
setwd("~/xjj/drug/drug_result/HDACi_chemo/chemohigh_HDAClow/1_RNAseq")
df<-read.table("chemohigh_HDAClow-other2-all-DESeq2_AUC_RNAseq_DEGs.txt",sep = '\t',header= T,row.names = 1)
DEGs<-df[df$pvalue<0.05,]
DEGs[which(DEGs$log2FoldChange>0),'up_down']<-'up'
DEGs[which(DEGs$log2FoldChange<0),'up_down']<-'down'
#pvalue padj
Acetylase_DEGs<-intersect(toupper(Acetylase[,1]),toupper(rownames(DEGs)))
Acetylase_info1<-Acetylase[match(Acetylase_DEGs,toupper(Acetylase[,1])),]
DEGs$gene<-toupper(rownames(DEGs))
Acetylase_info2<-DEGs[match(Acetylase_DEGs,toupper(rownames(DEGs))),7:8]
Acetylase_info<-cbind(Acetylase_info1,Acetylase_info2)
library(gplots)
library(VennDiagram)
Acetylase_name<- unique(toupper(Acetylase[,1]))
DEGs_name <- unique(toupper(rownames(DEGs)))
input  <-list(Acetylase_name,DEGs_name)
venn(input,showSetLogicLabel=TRUE)
tmp <- venn(input)
int<-attr(tmp, "intersections")
setwd("~/xjj/drug/drug_result/HDACi_chemo/chemohigh_HDAClow/1_RNAseq/Acetylase")
write.table(Acetylase_info,"Acetylase_info_P0.05.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
library(VennDiagram)
venn.diagram(x=list(Acetylase=Acetylase_name,DEGs=DEGs_name),cex = 1,margin = 0.1, "Acetylase_info_P005.png",fill=c("red","blue"))

setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
library("data.table")
exp<-fread("star_rsem.GeneSymbol.TPM.xls",header=T,data.table=F)
sensitivity<-chemohigh_HDAClow
resistance<-other2

resistance_exp<-exp[,match(resistance,colnames(exp))]
sensitivity_exp<-exp[,match(sensitivity,colnames(exp))]
geneid<-exp[,1]
resistance_exp_acetylase<-resistance_exp[match(Acetylase_info[,1],geneid),]
sensitivity_exp_acetylase<-sensitivity_exp[match(Acetylase_info[,1],geneid),]
resistance_exp11<-apply(resistance_exp_acetylase,1,mean)
sensitivity_exp11<-apply(sensitivity_exp_acetylase,1,mean)
resistance_sd<-apply(resistance_exp_acetylase,1,sd)
sensitivity_sd<-apply(sensitivity_exp_acetylase,1,sd)

pvalue<-df[match(Acetylase_info[,1],rownames(df)),5]
Acetylase_info_value<-cbind(Acetylase_info,sensitivity_exp11,resistance_exp11,pvalue)
colnames(Acetylase_info_value)<-c(colnames(Acetylase_info),"sensitive","resistant","pvalue")
setwd("~/xjj/drug/drug_result/HDACi_chemo/chemohigh_HDAClow/1_RNAseq/Acetylase")
write.table(Acetylase_info_value,"chemohigh_HDAClow_Acetylase_info_P0.05_value.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

library(ggplot2)
df2<-cbind(rbind(as.matrix(Acetylase_info_value[,1]),as.matrix(Acetylase_info_value[,1])),
           rbind(as.matrix(sensitivity_exp11),as.matrix(resistance_exp11)),
           rbind(as.matrix(rep("sensitive",length(Acetylase_info_value[,1]))),as.matrix(rep("resistant",length(Acetylase_info_value[,1])))),
           rbind(as.matrix(Acetylase_info_value$sensitive+sensitivity_sd),as.matrix(Acetylase_info_value$resistant+resistance_sd)),
           rbind(as.matrix(Acetylase_info_value$sensitive-sensitivity_sd),as.matrix(Acetylase_info_value$resistant-resistance_sd)),
           rbind(as.matrix(Acetylase_info_value[,8]),as.matrix(Acetylase_info_value[,8])))
colnames(df2)<-c("gene","value","type","valuesd1","valuesd2","pvalue")
setwd("~/xjj/drug/drug_result/HDACi_chemo/chemohigh_HDAClow/1_RNAseq/Acetylase")
write.table(df2,"chemohigh_HDAClow_Acetylase_info_P0.05_value_plot.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

df3<-read.table("chemohigh_HDAClow_Acetylase_info_P0.05_value_plot.txt",sep = '\t',header= T)
temp=df3[df3[,3]=="sensitive",1]
df3[,1]=factor(df3[,1],levels=temp)
df3[which(df3$pvalue>0.01),'pva']<-"*"
df3[which(df3$pvalue<0.01 & df3$pvalue>0.001),'pva']<-"**"
#df3[(length(Acetylase_info_value[,1])+1):34,7]=NA
df3$valuesd1<-df3$value+1
p3 <- ggplot(df3, aes(x = gene, y = value, fill = type)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge(),width=0.4)+ 
  geom_errorbar(aes(ymin = valuesd1, ymax = value), width = 0.2, position = position_dodge(0.4))+
  labs(x = "Differently expression of histone-modifying enzymes", y = "mRNA average expression levels (TPM)")+
  geom_text(aes(y=60,label = pva),size = 3)
plot(p3)

setwd("~/xjj/drug/drug_result/HDACi_chemo/chemohigh_HDAClow/1_RNAseq/Acetylase")
ggsave("chemohigh_HDAClow_Acetylase_info_P0.05_value_plot.pdf",p3,width = 10, height = 6)



#################################################################################
################  clusterProfiler
rm(list=ls())
library(org.Hs.eg.db)
library(clusterProfiler)
setwd("~/xjj/drug/drug_result/HDACi_chemo/HDAChigh_chemolow/1_RNAseq")
up_gene<-read.table("HDAChigh_chemolow_other1_RNAseq_P0.05_FC0_down_DEGs.txt",header=T,sep="\t")
library(stringr)
gene=bitr(up_gene[,1],fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
kegg<-enrichKEGG(gene[,2],organism = 'hsa',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'BH',
                 minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.05,use_internal_data = FALSE)
#keu<-barplot(kegg,showCategory=15,drop=T,x = "GeneRatio",color = "pvalue")
keu<-barplot(kegg,showCategory=15,drop=T,x = "GeneRatio",color = "p.adjust")

ee<-kegg@result
#ee1<-ee[which(ee$pvalue<0.05),]
ee1<-ee[which(ee$p.adjust<0.05),]

#p.adjust
#dotplot(kegg,showCategory=8,x = "GeneRatio",color = "p.adjust",title = "DEGs-P05-KEGG")
#enrichplot::gseaplot2(ee1,1,pvalue_table=T,color="#086538")
enrich_gene<-ee1$geneID
pathway_gene<-unique(unlist(strsplit(enrich_gene,split="/")))
pathway_gene2=bitr(pathway_gene,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
for(y in 1:nrow(ee1)){
  ee1[y,8]
  b1<-matrix(unlist(strsplit(ee1[y,8],split="/")),ncol=1)
  gene6=bitr(b1,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
  ee1[y,10]<-paste0(gene6[,2],collapse ="/")
}
setwd("~/xjj/drug/drug_result/HDACi_chemo/HDAChigh_chemolow/1_RNAseq/pathway")
write.table(ee1,"HDAChigh_chemolow_kegg_DEGP0.05_pathway_FDR0.05_up_drug_target.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(pathway_gene2,"HDAChigh_chemolow_kegg_DEGP0.05_pathway_FDR0.05_up_drug_target_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
ggsave("HDAChigh_chemolow_kegg_DEGP0.05_pathway_FDR0.05_up_drug_target_gene.pdf",keu,width = 8, height = 6)


setwd("~/xjj/drug/drug_result/HDACi_chemo/HDAChigh_chemolow/1_RNAseq/pathway")
write.table(ee1,"HDAChigh_chemolow_kegg_DEGP0.05_pathway_FDR0.05_down_drug_target.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(pathway_gene2,"HDAChigh_chemolow_kegg_DEGP0.05_pathway_FDR0.05_down_drug_target_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
ggsave("HDAChigh_chemolow_kegg_DEGP0.05_pathway_FDR0.05_down_drug_target_gene.pdf",keu,width = 8, height = 6)

go<-enrichGO(gene[,2],OrgDb = org.Hs.eg.db,ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.05,keyType = 'ENTREZID')
ked<-dotplot(go,showCategory=10)
a<-go@result
go_BP<-a[which(a$p.adjust<0.05),]

enrich_genego<-go_BP$geneID
pathway_genego<-unique(unlist(strsplit(enrich_genego,split="/")))
pathway_genego2=bitr(pathway_genego,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
setwd("~/xjj/drug/drug_result/HDACi_chemo/HDAChigh_chemolow/1_RNAseq/pathway")
write.table(pathway_genego2,"HDAChigh_chemolow_GO_DEG_pathway_FDR0.05_up_drug_target_gene.txt",sep = '\t',col.names = F,row.names = F,quote = FALSE)##数据输出
write.table(go_BP,"HDAChigh_chemolow_GO_DEG_pathway_FDR0.05_up_drug_target.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
ggsave("HDAChigh_chemolow_GO_DEG_pathway_FDR0.05_up_drug_target_gene.pdf",ked,width = 8, height = 6)

setwd("~/xjj/drug/drug_result/HDACi_chemo/HDAChigh_chemolow/1_RNAseq/pathway")
write.table(pathway_genego2,"HDAChigh_chemolow_GO_DEG_pathway_FDR0.05_down_drug_target_gene.txt",sep = '\t',col.names = F,row.names = F,quote = FALSE)##数据输出
write.table(go_BP,"HDAChigh_chemolow_GO_DEG_pathway_FDR0.05_down_drug_target.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
ggsave("HDAChigh_chemolow_GO_DEG_pathway_FDR0.05_down_drug_target_gene.pdf",ked,width = 8, height = 6)

#######################    ATACseq   ################################################
#############  DApeaks
# I:HDAC_sen II:chemo_sen III:all_sen IV:all_res
#### I II+III+IV
HDAChigh_chemolow<-rownames(HDAC_sen)
other1<-c(rownames(chemo_sen),rownames(all_sen),rownames(all_res))
#### II I+III+IV
chemohigh_HDAClow<-rownames(chemo_sen)
other2<-c(rownames(HDAC_sen),rownames(all_sen),rownames(all_res))

#### III I+II+IV
allhigh<-rownames(all_sen)
other3<-c(rownames(HDAC_sen),rownames(chemo_sen),rownames(all_res))

#### IV I+II+III
alllow<-rownames(all_res)
other4<-c(rownames(HDAC_sen),rownames(chemo_sen),rownames(all_sen))


rm(list=ls())
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_count_name<-read.table("peak_count.txt",sep="\t",header=T) #没有全0的

resistance<-other4
sensitivity<-alllow

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
setwd("~/xjj/drug/drug_result/HDACi_chemo/alllow/2_ATACseq")
write.table(res,"alllow-other4-all-DESeq2_chemotherapy_AUC_ATAC.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

#setwd("~/xjj/drug/drug_result/HDAC_drug/heatmap")
#res<-read.table("sensitivity-resistance-all-DESeq2_HDAC_chemotherapy_AUC_ATAC.txt",header = TRUE,sep = "\t")
#resSig<-res[which(res$padj<0.05 & abs(res$log2FoldChange>1)),]
resSig<-res[which(res$pvalue<0.05 & abs(res$log2FoldChange)>1),]
resSig<-res[which(res$pvalue<0.01),]

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
resistance_peak11<-apply(resistance_peak1,1,mean)
sensitivity_peak11<-apply(sensitivity_peak1,1,mean)
sum((sensitivity_peak11[up1]-resistance_peak11[up1])>0)
sum((sensitivity_peak11[down1]-resistance_peak11[down1])<0)


library(dplyr)
a_up<-apply(up_gene,1,function(x) unlist(strsplit(as.character(x), "_")))
d_up<-t(a_up)
e_up<-data.frame(rep("+",nrow(d_up)))
peaks_up<-cbind(up_gene,d_up,e_up)
colnames(peaks_up)<-c("peak_id","chr","start","end","strand")
out_file<-"alllow-other4-AUC-P0.05-logFC1"
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
setwd("~/xjj/drug/drug_result/HDACi_chemo/allpeak")
write.table(peaks_bed_up,file="all_peaks.bed",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')


###############################四类中共有的peaks和特有的peaks
folder<-c("HDAChigh_chemolow","chemohigh_HDAClow","allhigh","alllow")
i=1
setwd(paste("~/xjj/drug/drug_result/HDACi_chemo/",folder[i],"/2_ATACseq",sep=""))
oneup<-read.table(paste(folder[i],"-other",i,"-AUC-P0.01-logFC0-DESeq2-up.txt",sep=""),sep="\t",header=T)
for(i in 2:4){
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo/",folder[i],"/2_ATACseq",sep=""))
  twoup<-read.table(paste(folder[i],"-other",i,"-AUC-P0.01-logFC0-DESeq2-up.txt",sep=""),sep="\t",header=T)
  oneup<-as.data.frame(intersect(as.character(oneup[,1]),as.character(twoup[,1])))
}



setwd("~/xjj/drug/drug_result/HDACi_chemo/alllow/2_ATACseq")
alllow<-read.table("alllow-other4-AUC-P0.01-logFC0-DESeq2-up.txt",sep="\t",header=T)
setwd("~/xjj/drug/drug_result/HDACi_chemo/allhigh/2_ATACseq")
allhigh<-read.table("allhigh-other3-AUC-P0.01-logFC0-DESeq2-up.txt",sep="\t",header=T)
setwd("~/xjj/drug/drug_result/HDACi_chemo/HDAChigh_chemolow/2_ATACseq")
HDAChigh_chemolow<-read.table("HDAChigh_chemolow-other1-AUC-P0.01-logFC0-DESeq2-up.txt",sep="\t",header=T)
setwd("~/xjj/drug/drug_result/HDACi_chemo/chemohigh_HDAClow/2_ATACseq")
chemohigh_HDAClow<-read.table("chemohigh_HDAClow-other2-AUC-P0.01-logFC0-DESeq2-up.txt",sep="\t",header=T)

one1<-intersect(as.character(alllow[,1]),as.character(allhigh[,1]))#0
one2<-intersect(as.character(alllow[,1]),as.character(HDAChigh_chemolow[,1]))#0
one3<-intersect(as.character(alllow[,1]),as.character(chemohigh_HDAClow[,1]))#1

two1<-intersect(as.character(allhigh[,1]),as.character(HDAChigh_chemolow[,1]))#0
two2<-intersect(as.character(allhigh[,1]),as.character(chemohigh_HDAClow[,1]))#0

three1<-intersect(as.character(HDAChigh_chemolow[,1]),as.character(chemohigh_HDAClow[,1]))#18

##用这些特征画热图
HDAChigh_chemolow1<-rownames(HDAC_sen)
chemohigh_HDAClow1<-rownames(chemo_sen)
allhigh1<-rownames(all_sen)
alllow1<-rownames(all_res)
sample_order<-c(HDAChigh_chemolow1,chemohigh_HDAClow1,allhigh1,alllow1)
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_RPKM<-read.table("peak_RPKM.txt",sep="\t",header=T)
colnames(peak_RPKM)<-gsub("\\.","-",colnames(peak_RPKM))
peaks<-c(as.character(HDAChigh_chemolow[,1]),as.character(chemohigh_HDAClow[,1]),as.character(allhigh[,1]),as.character(alllow[,1]))
peak_RPKM1<-peak_RPKM[match(peaks,peak_RPKM[,1]),]
peaks1<-matrix(c(rep("Class1",nrow(HDAChigh_chemolow)),rep("Class2",nrow(chemohigh_HDAClow)),rep("Class3",nrow(allhigh)),rep("Class4",nrow(alllow))),ncol=1)
peak_RPKM_order<-peak_RPKM1[,match(sample_order,colnames(peak_RPKM1))]

col_cut=c(length(HDAChigh_chemolow1),length(HDAChigh_chemolow1)+length(chemohigh_HDAClow1),length(HDAChigh_chemolow1)+length(chemohigh_HDAClow1)+length(allhigh1))
row_cut=c(nrow(HDAChigh_chemolow),nrow(HDAChigh_chemolow)+nrow(chemohigh_HDAClow),nrow(HDAChigh_chemolow)+nrow(chemohigh_HDAClow)+nrow(allhigh))
library(pheatmap)
data0=cbind(peaks1,peak_RPKM_order)

anno1=c(rep("Class1",length(HDAChigh_chemolow1)),rep("Class2",length(chemohigh_HDAClow1)),rep("Class3",length(allhigh1)),rep("Class4",length(alllow1)))
anno1=data.frame(anno1)

rownames(anno1)=as.character(t(colnames(data0))[2:38])
ann_colors = list(
  anno1 = c(Class1 = "#4dbbd5", Class2 = "#00a087",Class3 = "#f39b7f",Class4="#8491b4"), anno=c("Class1"= "#4dbbd5","Class2"= "#00a087","Class3"= "#f39b7f","Class4"="#8491b4"))
H=6
W=5
library(RColorBrewer)
color1=colorRampPalette(c("#3363aa","#1384d5","#23b2ae"))(4)
color10=colorRampPalette(c(color1[1],color1[2]))(35)
color11=colorRampPalette(c(color1[2],color1[4]))(15)
color2=colorRampPalette(c( "#23b2ae","#c5bc5e","#faf513"))(4)
color20=colorRampPalette(c(color2[1],color2[3]))(15)
color21=colorRampPalette(c(color2[3],color2[4]))(35)
mycolor=c(color10,color11,color20,color21)
anno=data.frame(data0$peaks1)
rownames(data0)=as.character(1:dim(data0)[1])
rownames(anno)=rownames(data0)
colnames(anno)="anno"
setwd("~/xjj/drug/drug_result/HDACi_chemo/allpeak")
pdf(paste0("finalheatmap.pdf"),height=H,width=W)
p<-pheatmap(data0[2:38],scale="row",color = mycolor,fontsize=7,fontsize_row = 70*H/dim(data0)[1], fontsize_col = 7,cluster_cols=FALSE, cluster_rows=FALSE,annotation_col=anno1,gaps_col = col_cut,gaps_row = row_cut, annotation_row=anno,annotation_colors = ann_colors)
ggsave("ATACseqheatmap.pdf",p,width = 10, height = 8)
dev.off()
#### 注意要对每一行的数据进行标准化

#######################################################################################
####  homer注释完的信息,画饼图
rm(list=ls())
folder<-c("HDAChigh_chemolow","chemohigh_HDAClow","allhigh","alllow")
library(data.table)
library(foreach)
library(stringr)
library(ggplot2)

ATAC_pie_result<-data.frame(matrix(0,length(folder),5))
ratio<-NULL
foreach(k=1:length(folder))%do%{
  file_name=as.character(folder[k])
  cat(file_name,"\n")
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo/",file_name,"/2_ATACseq/annotation",sep=""))
  homer_anno_up<-fread(paste(file_name,"-other",k,"-AUC-P0.01-logFC0-DESeq2-up-annotation.txt",sep=""),header=T,data.table=F)
  cat(length(unique(homer_anno_up$`Gene Name`)),"\n")
  ATAC_pie_result[k,1]<-file_name
  ATAC_pie_result[k,2]<-length(unique(homer_anno_up$`Gene Name`))
  colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
  Annotation<-apply(as.data.frame(homer_anno_up$Annotation),1,function(x) str_replace_all(as.character(x),"\\((.*?)\\)",""))##取出括号前的字符
  Annotation1<-as.data.frame(Annotation)
  Annotation2<-data.frame(lapply(strsplit(as.character(Annotation1[,1]),'\\.'), function(x) x[1])%>%unlist())
  position<-as.character(unique(Annotation2[,1]))
  result<-data.frame()
  for(i in 1:length(position)){
    result[i,1]<-position[i]
    result[i,2]<-sum(Annotation %in% position[i])
    result[i,3]<-(sum(Annotation %in% position[i])/(length(Annotation)))*100
  }
  z1<-c()
  for(z in 1:nrow(result)){
    z11<-result[z,1]
    z12<-round(as.numeric(result[z,3]), digits = 2)
    z13<-c(z11,z12)
    z1<-c(z1,z13)
  }
  ATAC_pie_result[k,3]<-paste0(z1,collapse =",")
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
  #print(p) #显示饼图
  ggsave(paste("ATAC_pie_",file_name,"-0.01-P-up.pdf",sep=""),p,width = 10, height = 4)
  
  homer_anno_down<-fread(paste(file_name,"-other",k,"-AUC-P0.01-logFC0-DESeq2-down-annotation.txt",sep=""),header=T,data.table=F)
  cat(length(unique(homer_anno_down$`Gene Name`)),"\n")
  ATAC_pie_result[k,4]<-length(unique(homer_anno_down$`Gene Name`))
  colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
  Annotationd<-apply(as.data.frame(homer_anno_down$Annotation),1,function(x) str_replace_all(as.character(x),"\\((.*?)\\)",""))##取出括号前的字符
  Annotation1d<-as.data.frame(Annotationd)
  Annotation2d<-data.frame(lapply(strsplit(as.character(Annotation1d[,1]),'\\.'), function(x) x[1])%>%unlist())
  positiond<-as.character(unique(Annotation2d[,1]))
  resultd<-data.frame()
  for(i in 1:length(positiond)){
    resultd[i,1]<-positiond[i]
    resultd[i,2]<-sum(Annotationd %in% positiond[i])
    resultd[i,3]<-(sum(Annotationd %in% positiond[i])/(length(Annotationd)))*100
  }
  z1d<-c()
  for(z in 1:nrow(resultd)){
    z11<-resultd[z,1]
    z12<-round(as.numeric(resultd[z,3]), digits = 2)
    z13<-c(z11,z12)
    z1d<-c(z1d,z13)
  }
  ATAC_pie_result[k,5]<-paste0(z1d,collapse =",")
  dtd = data.frame(A = resultd[,3],B = resultd[,1])#建立数据框
  dtd = dtd[order(dtd$A, decreasing = TRUE),]   ## 用 order() 让数据框的数据按 A 列数据从大到小排序
  myLabel = as.vector(dtd$B)   ## 转成向量，否则图例的标签可能与实际顺序不一致
  myLabel = paste(myLabel, "(", round(dtd$A / 1, 4), "%)", sep = "")   ## 用 round() 对结果保留两位小数
  p = ggplot(dtd, aes(x = "", y = A, fill = B)) + #创建坐标轴
    geom_bar(stat = "identity") + 
    geom_bar(stat = "identity", width = 1) +   #当width >= 1 时中心的杂点将消失
    coord_polar(theta = "y") +  # 把柱状图折叠成饼图（极坐标）
    labs(x = "", y = "", title = "") +  # 将横纵坐标的标签设为空
    theme(axis.ticks = element_blank()) +  # 将左上角边框的刻度去掉
    theme(legend.title = element_blank(), legend.position = "left")+   ## 将图例标题设为空，并把图例方放在左边位置
    scale_fill_discrete(breaks = dtd$B, labels = myLabel)+   # 将原来的图例标签换成现在的myLabel
    theme(axis.text.x = element_blank())   ## 去掉饼图的外框上的数值，即去除原柱状图的X轴，把X轴的刻度文字去掉
  #geom_text(aes(y = A/2 + c(0, cumsum(A)[-length(A)]), x = sum(A)/20, label = myLabel), size = 5)   # 在图中加上百分比：x 调节标签到圆心的距离, y 调节标签的左右位置
  ggsave(paste("ATAC_pie_",file_name,"-0.01-P-down.pdf",sep=""),p,width = 10, height = 4)
  
  down1<-resultd[order(resultd[,1]),]
  up1<-result[order(result[,1]),]
  dat_m1<-rbind(down1[,c(1,3)],up1[,c(1,3)])
  m1<-matrix(c(rep(paste(file_name,"down",sep="_"),8),rep(paste(file_name,"up",sep="_"),8)),ncol=1)
  dat_m<-cbind(m1,dat_m1)
  ratio<-rbind(ratio,dat_m)
  
}
colnames(ratio) = c('Group','Type','value')

colnames(ATAC_pie_result)<-c("drug","number of anno_up_gene","anno_up_gene","number of anno_down_gene","anno_down_gene")
setwd("~/xjj/drug/drug_result/HDACi_chemo/allpeak")
write.table(ATAC_pie_result,"ATAC_pie_result_P0.01.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

########   堆积图   ###################
rm(list=ls())
library(reshape2)
dat_m = ratio
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
        axis.text.x = element_text(color = 'black',size = 12,angle = 45,vjust = 0.6), #x轴标签偏转45°，并下降0.5
        panel.grid = element_blank(),
        legend.position = 'right',
        legend.key.height = unit(0.6,'cm'),#定义图例中色块的高度
        legend.text = element_text(color = 'black',size = 10))
plot(p)
setwd("~/xjj/drug/drug_result/HDACi_chemo/allpeak/ATACseq")
ggsave(p, filename = 'annotation-ratio-all4group-down-up.pdf', width = 8, height = 6, dpi = 600)


#########################################################################
#####################  符合在TSS100KB之内条件的注释基因
rm(list=ls())
library(gplots)
library(VennDiagram)
library(data.table)
library(foreach)
folder<-c("HDAChigh_chemolow","chemohigh_HDAClow","allhigh","alllow")

TSS100KB_anno_result<-data.frame(matrix(0,length(folder),7))
foreach(k=1:length(folder))%do%{
  file_name=as.character(folder[k])
  cat(file_name,"\n")
  TSS100KB_anno_result[k,1]<-file_name
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo/",file_name,"/2_ATACseq/annotation",sep=""))
  homer_anno_down<-fread(paste(file_name,"-other",k,"-AUC-P0.01-logFC0-DESeq2-down-annotation.txt",sep=""),header=T,data.table=F)
  colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
  Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]
  
  homer_anno_up<-fread(paste(file_name,"-other",k,"-AUC-P0.01-logFC0-DESeq2-up-annotation.txt",sep=""),header=T,data.table=F)
  colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
  Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
  
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo/",file_name,"/1_RNAseq",sep=""))
  
  up_gene<-read.table(paste(file_name,"_other",k,"_RNAseq_P0.05_FC0_up_DEGs.txt",sep=""),header=T,sep="\t")
  down_gene<-read.table(paste(file_name,"_other",k,"_RNAseq_P0.05_FC0_down_DEGs.txt",sep=""),header=T,sep="\t")
  all<-rbind(up_gene,down_gene)
  DEGs<-as.character(all[,1])
  down_DEGs_DApeaks<-intersect(Anno_gene_100Kb_down$`Gene Name`,DEGs)
  cat(paste("Down peak int DEGs is",length(down_DEGs_DApeaks),sep=" "),"\n")
  up_DEGs_DApeaks<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)
  cat(paste("Up peak int DEGs is",length(up_DEGs_DApeaks),sep=" "),"\n")
  TSS100KB_anno_result[k,2]<-length(up_DEGs_DApeaks)
  TSS100KB_anno_result[k,5]<-length(down_DEGs_DApeaks)
  down_peaks_gene <- Anno_gene_100Kb_down$`Gene Name`
  up_peaks_gene <- Anno_gene_100Kb_up$`Gene Name`
  input  <-list(unique(down_peaks_gene),unique(up_peaks_gene),unique(DEGs))
  plot(venn(input,showSetLogicLabel=TRUE))
  tmp <- venn(input)
  int<-attr(tmp, "intersections")
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo/",file_name,"/3_intersect/venn",sep=""))
  venn.diagram(x=list(downpeaks=unique(down_peaks_gene),uppeaks=unique(up_peaks_gene),RNAseq=unique(DEGs)), paste("DEGs0.05_int_DApeaks0.01_",file_name,".png",sep=""),fill=c("red","green","blue"),margin = 0.1)
  
  up1<-int[["B:C"]]
  up2<-int[["A:B:C"]]
  up3<-c(up1,up2)
  int_up<-as.data.frame(intersect(up3,as.character(up_gene$gene_id)))
  colnames(int_up)<-"int_up"
  int_up_peak<-as.data.frame(homer_anno_up[homer_anno_up$`Gene Name`%in%int_up[,1],1])
  colnames(int_up_peak)<-"int_up_peak"
  
  down1<-int[["A:C"]]
  down2<-int[["A:B:C"]]
  down3<-c(down1,down2)
  int_down<-as.data.frame(intersect(down3,as.character(down_gene$gene_id)))
  colnames(int_down)<-"int_down"
  int_down_peak<-as.data.frame(homer_anno_down[homer_anno_down$`Gene Name`%in%int_down[,1],1])
  colnames(int_down_peak)<-"int_down_peak"
  
  down_DEGs_down_DApeaks<-as.data.frame(intersect(Anno_gene_100Kb_down$`Gene Name`,down_gene[,1]))
  colnames(down_DEGs_down_DApeaks)<-"down_DEGs_down_DApeaks"
  up_DEGs_up_DApeaks<-as.data.frame(intersect(Anno_gene_100Kb_up$`Gene Name`,up_gene[,1]))
  colnames(up_DEGs_up_DApeaks)<-"up_DEGs_up_DApeaks"
  
  TSS100KB_anno_result[k,3]<-nrow(int_up)
  TSS100KB_anno_result[k,6]<-nrow(int_down)
  
  setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
  library("data.table")
  exp1<-fread("star_rsem.GeneSymbol.readscount.xls",header=T,data.table=F)
  exp2<-floor(exp1[,-1])
  delete<-apply(exp2,1,function(x) mean(x==0))
  cou<-which(delete>0.75)####在超过75%样本以上的都是0的基因删去
  exp<-exp2[(-cou),]
  Ku=nrow(int_up)
  Nu=nrow(exp) #all gene
  Mu=length(unique(Anno_gene_100Kb_up$`Gene Name`)) #all DARs annotation gene
  nu=nrow(up_gene) #all DEGs
  pu<-phyper(Ku-1,Mu, Nu-Mu, nu, lower.tail=F)
  TSS100KB_anno_result[k,4]<-pu
  
  Kd=nrow(int_down)
  Nd=nrow(exp) #all gene
  Md=length(unique(Anno_gene_100Kb_down$`Gene Name`)) #all DARs annotation gene
  nd=nrow(down_gene) #all DEGs
  pd<-phyper(Kd-1,Md, Nd-Md, nd, lower.tail=F)
  TSS100KB_anno_result[k,7]<-pd
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo/",file_name,"/3_intersect/venn",sep=""))
  write.table(int_up,"int_up_peak0.01_all_DEG0.05_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  write.table(int_down,"int_down_peak0.01_all_DEG0.05_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  write.table(int_up_peak,"int_up_peak0.01_all_DEG0.05_peak.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  write.table(int_down_peak,"int_down_peak0.01_all_DEG0.05_peak.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  write.table(up_DEGs_up_DApeaks,"int_up_peak0.01_up_DEG0.05_peak.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  write.table(down_DEGs_down_DApeaks,"int_down_peak0.01_down_DEG0.05_peak.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  
}
colnames(TSS100KB_anno_result)<-c("group","Up peak int DEGs","Up peak int up-DEGs","Hypergeometric-up","Down peak int DEGs","Down peak int down-DEGs","Hypergeometric-down")
setwd("~/xjj/drug/drug_result/HDACi_chemo/allpeak/ATACseq")
write.table(TSS100KB_anno_result,"TSS100KB_anno_result_result_4group.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

##################################################################################
##############  DEpeaks在不同样本变化的倍数与DEGs在不同样本变化的倍数之间的相关性(一个基因对应多个peaks)
rm(list=ls())
folder<-c("HDAChigh_chemolow","chemohigh_HDAClow","allhigh","alllow")
library(data.table)
library(foreach)
library(ggplot2)
library(ggpubr)
library(ggrepel)

DEGs_corre_DApeaks<-data.frame(matrix(0,length(folder),5))
foreach(k=1:length(folder))%do%{
  file_name=as.character(folder[k])
  cat(file_name,"\n")
  DEGs_corre_DApeaks[k,1]<-file_name
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo/",file_name,"/2_ATACseq/annotation",sep=""))
  homer_anno_down<-fread(paste(file_name,"-other",k,"-AUC-P0.01-logFC0-DESeq2-down-annotation.txt",sep=""),header=T,data.table=F)
  colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
  Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]
  
  homer_anno_up<-fread(paste(file_name,"-other",k,"-AUC-P0.01-logFC0-DESeq2-up-annotation.txt",sep=""),header=T,data.table=F)
  colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
  Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
  
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo/",file_name,"/1_RNAseq",sep=""))
  up_gene<-read.table(paste(file_name,"_other",k,"_RNAseq_P0.05_FC0_up_DEGs.txt",sep=""),header=T,sep="\t")
  down_gene<-read.table(paste(file_name,"_other",k,"_RNAseq_P0.05_FC0_down_DEGs.txt",sep=""),header=T,sep="\t")
  all<-rbind(up_gene,down_gene)
  DEGs<-as.character(all[,1])
  down_DEGs_DApeaks<-intersect(Anno_gene_100Kb_down$`Gene Name`,DEGs)
  up_DEGs_DApeaks<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)
  
  int_all_DEGs<-unique(c(down_DEGs_DApeaks,up_DEGs_DApeaks))
  all_peaks_anno<-rbind(homer_anno_down,homer_anno_up)
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo/",file_name,"/2_ATACseq",sep=""))
  DApeak<-fread(paste(file_name,"-other",k,"-all-DESeq2_chemotherapy_AUC_ATAC.txt",sep=""),header=T,data.table=F)
  
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
  #rownames(DEG_DApeak_FC)<-different_FC[,1]
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo/",file_name,"/3_intersect/correlation",sep=""))
  write.table(DEG_DApeak_FC,"DEGs_FC_DApeaks_FC_correlation_comprehensive.txt",sep = '\t',col.names = T,row.names = T,quote = FALSE)##数据输出
  
  #pearson", "kendall", "spearman
  cor_result<-cor.test(DEG_DApeak_FC[,1],DEG_DApeak_FC[,2],alternative = "two.sided",method = "spearman")
  DEGs_corre_DApeaks[k,2]<-cor_result$estimate
  DEGs_corre_DApeaks[k,3]<-cor_result$p.value
  
  cor_result1<-cor.test(DEG_DApeak_FC[,1],DEG_DApeak_FC[,2],alternative = "two.sided",method = "pearson")
  DEGs_corre_DApeaks[k,4]<-cor_result1$estimate
  DEGs_corre_DApeaks[k,5]<-cor_result1$p.value
  
  ################################   画相关性分析图
  dat1<-as.data.frame(DEG_DApeak_FC)
  library(ggplot2)
  library(ggpubr)
  #a11<-c(1:10)
  #a21<-c((nrow(dat1)-9):nrow(dat1))
  #a3<-c(a11,a21)
  #want_dat<-dat1[order(dat1[,1])[a3],]
  ps<-ggplot(data=dat1, aes(x=DEGslog2FC, y=DApeakslog2FC))+geom_point(color="red",size=0.3)+stat_smooth(method="lm",se=FALSE)+stat_cor(data=dat1, method = "spearman")+
    ggtitle("spearman-DApeaksP0.01-DEGsP0.05") +
    theme(plot.title = element_text(hjust = 0.5))
    #geom_text_repel(
    #  data = want_dat[,c(1:2)],
    #  aes(label = different_FC[order(dat1[,1])[a3],1]),
    #  size = 3,
    #  color = "black",
    #  segment.color = "black", show.legend = FALSE )
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo/",file_name,"/3_intersect/correlation",sep=""))
  ggsave(paste(file_name,"correlation-DApeaksP0.01-DEGsP0.05-spearman-comprehensive.pdf",sep=""),ps,width = 10, height = 10)
  
  pp<-ggplot(data=dat1, aes(x=DEGslog2FC, y=DApeakslog2FC))+geom_point(color="red",size=0.1)+stat_smooth(method="lm",se=FALSE)+stat_cor(data=dat1, method = "pearson")+
    ggtitle("pearson-DApeaksP0.01-DEGsP0.05") +
    theme(plot.title = element_text(hjust = 0.5))
    #geom_text_repel(
    #  data = want_dat[,c(1:2)],
    #  aes(label = different_FC[order(dat1[,1])[a3],1]),
    #  size = 3,
    #  color = "black",
    #  segment.color = "black", show.legend = FALSE )
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo/",file_name,"/3_intersect/correlation",sep=""))
  ggsave(paste(file_name,"correlation-DApeaksP0.01-DEGsP0.05-pearson-comprehensive.pdf",sep=""),pp,width = 10, height = 10)
  #plot(p)
}
#pearson", "kendall", "spearman
colnames(DEGs_corre_DApeaks)<-c("group","spearman-estimate","spearman-pvalue","pearson-estimate","pearson-pvalue")
setwd("~/xjj/drug/drug_result/HDACi_chemo/allpeak/ATACseq")
write.table(DEGs_corre_DApeaks,"DEGs_FC_corre_DApeaks_FC_result_4group_comprehensive.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

###################################################################################
##############  DEpeaks_FC与DEGs_FC之间的相关性(一个基因对应一个peaks)
rm(list=ls())
folder<-c("HDAChigh_chemolow","chemohigh_HDAClow","allhigh","alllow")
library(data.table)
library(foreach)
library(ggplot2)
library(ggpubr)
library(ggrepel)

DEGs_corre_DApeaks<-data.frame(matrix(0,length(folder),5))
foreach(k=1:length(folder))%do%{
  file_name=as.character(folder[k])
  cat(file_name,"\n")
  DEGs_corre_DApeaks[k,1]<-file_name
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo/",file_name,"/2_ATACseq/annotation",sep=""))
  homer_anno_down<-fread(paste(file_name,"-other",k,"-AUC-P0.01-logFC0-DESeq2-down-annotation.txt",sep=""),header=T,data.table=F)
  colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
  Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]
  
  homer_anno_up<-fread(paste(file_name,"-other",k,"-AUC-P0.01-logFC0-DESeq2-up-annotation.txt",sep=""),header=T,data.table=F)
  colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
  Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
  
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo/",file_name,"/1_RNAseq",sep=""))
  up_gene<-read.table(paste(file_name,"_other",k,"_RNAseq_P0.05_FC0_up_DEGs.txt",sep=""),header=T,sep="\t")
  down_gene<-read.table(paste(file_name,"_other",k,"_RNAseq_P0.05_FC0_down_DEGs.txt",sep=""),header=T,sep="\t")
  all<-rbind(up_gene,down_gene)
  DEGs<-as.character(all[,1])
  down_DEGs_DApeaks<-intersect(Anno_gene_100Kb_down$`Gene Name`,DEGs)
  up_DEGs_DApeaks<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)
  
  int_all_DEGs<-unique(c(down_DEGs_DApeaks,up_DEGs_DApeaks))
  all_peaks_anno<-rbind(homer_anno_down,homer_anno_up)
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo/",file_name,"/2_ATACseq",sep=""))
  DApeak<-fread(paste(file_name,"-other",k,"-all-DESeq2_chemotherapy_AUC_ATAC.txt",sep=""),header=T,data.table=F)
  down_DEGs_DApeaks<-intersect(Anno_gene_100Kb_down$`Gene Name`,DEGs)
  up_DEGs_DApeaks<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)
  int_all_DEGs<-unique(c(down_DEGs_DApeaks,up_DEGs_DApeaks))
  all_peaks_anno<-rbind(homer_anno_down,homer_anno_up)
  DEGs_peaks<-all_peaks_anno[match(int_all_DEGs,all_peaks_anno$`Gene Name`),1]
  DEGs_peaks_FC<-DApeak[match(DEGs_peaks,DApeak$peak_id),3]
  DEGs_FC<-all[match(int_all_DEGs,DEGs),3]
  gene<-as.character(all[match(int_all_DEGs,DEGs),1])
  #pearson", "kendall", "spearman
  cor_resultp<-cor.test(DEGs_peaks_FC, DEGs_FC,alternative = "two.sided",method = "pearson")
  cat(paste("pearson correlation is",cor_resultp$estimate,sep=" "),"\n")
  cat(paste("pearson pvalue is",cor_resultp$p.value,sep=" "),"\n")
  cor_results<-cor.test(DEGs_peaks_FC, DEGs_FC,alternative = "two.sided",method = "spearman")
  cat(paste("spearman correlation is",cor_results$estimate,sep=" "),"\n")
  cat(paste("spearman pvalue is",cor_results$p.value,sep=" "),"\n")
  DEGs_corre_DApeaks[k,1]<-file_name
  DEGs_corre_DApeaks[k,2]<-cor_resultp$estimate
  DEGs_corre_DApeaks[k,3]<-cor_resultp$p.value
  DEGs_corre_DApeaks[k,4]<-cor_results$estimate
  DEGs_corre_DApeaks[k,5]<-cor_results$p.value
  ################################   画相关性分析图
  a1<-matrix(DEGs_FC,ncol=1)
  a2<-matrix(DEGs_peaks_FC,ncol=1)
  dat<-cbind(a1,a2)
  dat1<-as.data.frame(dat)
  rownames(dat1)<-gene
  colnames(dat1)<-c("DEGslog2FC","DApeakslog2FC")
  a11<-c(1:10)
  a21<-c((nrow(dat1)-9):nrow(dat1))
  a3<-c(a11,a21)
  want_dat<-dat1[order(dat1[,1])[a3],]
  pp<-ggplot(data=dat1, aes(x=DEGslog2FC, y=DApeakslog2FC))+geom_point(color="red",size=0.1)+stat_smooth(method="lm",se=FALSE)+stat_cor(data=dat1, method = "pearson")+
    ggtitle(paste("pearson","-P0.05",sep="")) +
    theme(plot.title = element_text(hjust = 0.5))+
    geom_text_repel(
      data = want_dat[,c(1:2)],
      aes(label = rownames(want_dat)),
      size = 3,
      color = "black",
      segment.color = "black", show.legend = FALSE )
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo/",file_name,"/3_intersect/correlation",sep=""))
  ggsave(paste(file_name,"correlation-DApeaksP0.01-DEGsP0.05-pearson-onlyone.pdf",sep=""),pp,width = 10, height = 10)
  ps<-ggplot(data=dat1, aes(x=DEGslog2FC, y=DApeakslog2FC))+geom_point(color="red",size=0.1)+stat_smooth(method="lm",se=FALSE)+stat_cor(data=dat1, method = "spearman")+
    ggtitle(paste("spearman","-P0.05",sep="")) +
    theme(plot.title = element_text(hjust = 0.5))+
    geom_text_repel(
      data = want_dat[,c(1:2)],
      aes(label = rownames(want_dat)),
      size = 3,
      color = "black",
      segment.color = "black", show.legend = FALSE )
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo/",file_name,"/3_intersect/correlation",sep=""))
  ggsave(paste(file_name,"correlation-DApeaksP0.01-DEGsP0.05-spearman-onlyone.pdf",sep=""),ps,width = 10, height = 10)
}
#pearson", "kendall", "spearman
colnames(DEGs_corre_DApeaks)<-c("group","pearson-estimate","pearson-pvalue","spearman-estimate","spearman-pvalue")
setwd("~/xjj/drug/drug_result/HDACi_chemo/allpeak/ATACseq")
write.table(DEGs_corre_DApeaks,"DEGs_FC_corre_DApeaks_FC_result_4group_onlyone.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

#pearson", "kendall", "spearman
##画好看的图ggstatsplot

#### 使用超几何分布看看交叠是否随机出现，超几何分析在TSS100KB中完成。
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

#############################  对高低可及peaks-DEGs进行KEGG通路富集分析
rm(list=ls())
folder<-c("HDAChigh_chemolow","chemohigh_HDAClow","allhigh","alllow")
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(stringr)
freq_name<-"P0.05"

peaks_DEGs_KEGG<-data.frame(matrix(0,length(folder),5))
foreach(k=1:length(folder))%do%{
  file_name=as.character(folder[k])
  cat(file_name,"\n")
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo/",file_name,"/2_ATACseq/annotation",sep=""))
  homer_anno_down<-fread(paste(file_name,"-other",k,"-AUC-P0.01-logFC0-DESeq2-down-annotation.txt",sep=""),header=T,data.table=F)
  colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
  Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo/",file_name,"/1_RNAseq",sep=""))
  down_gene<-read.table(paste(file_name,"_other",k,"_RNAseq_P0.05_FC0_down_DEGs.txt",sep=""),header=T,sep="\t")
  
  DEGs<-as.character(down_gene[,1])
  down_DEGs_DApeaks<-intersect(Anno_gene_100Kb_down$`Gene Name`,DEGs)
  cat(paste("Down int gene is ",length(down_DEGs_DApeaks),sep=""),"\n")
  gened=bitr(down_DEGs_DApeaks,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
  keggd<-enrichKEGG(gened[,2],organism = 'hsa',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'none',
                    minGSSize = 10,maxGSSize = 500,qvalueCutoff = 1,use_internal_data = FALSE)
  eed<-keggd@result
  ee1d<-eed[which(eed$pvalue<0.05),]
  
  for(y in 1:nrow(ee1d)){
    ee1d[y,8]
    b1<-matrix(unlist(strsplit(ee1d[y,8],split="/")),ncol=1)
    gene6=bitr(b1,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
    ee1d[y,10]<-paste0(gene6[,2],collapse ="/")
  }
  
  cat(paste("KEGG pathway is",nrow(ee1d),sep=" "),"\n")
  ked<-barplot(keggd,showCategory=10,drop=T,x = "GeneRatio",color = "pvalue")
  #plot(ke)
  ee2d<-ee1d[1:10,]
  enrich_gened<-ee1d$geneID
  pathway_gened<-unique(unlist(strsplit(enrich_gened,split="/")))
  cat(paste("Top pathway gene is ",length(pathway_gened),sep=""),"\n")
  peaks_DEGs_KEGG[k,1]<-file_name
  peaks_DEGs_KEGG[k,2]<-nrow(ee1d)
  peaks_DEGs_KEGG[k,3]<-length(pathway_gened)
  pathway_gene2d=bitr(pathway_gened,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo/",file_name,"/3_intersect/pathway",sep=""))
  write.table(ee1d,paste("kegg_pathway_DEGs_DApeaks_",file_name,freq_name,"_down.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  write.table(pathway_gened,paste("kegg_pathway_DEGs_DApeaks_",file_name,freq_name,"_down_gene.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  ggsave(paste("KEGG pathway of",file_name,"-0.05-P-down-P0.05.pdf",sep=""),ked,width = 10, height = 5)
  
  
  
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo/",file_name,"/2_ATACseq/annotation",sep=""))
  homer_anno_up<-fread(paste(file_name,"-other",k,"-AUC-P0.01-logFC0-DESeq2-up-annotation.txt",sep=""),header=T,data.table=F)
  colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
  Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo/",file_name,"/1_RNAseq",sep=""))
  up_gene<-read.table(paste(file_name,"_other",k,"_RNAseq_P0.05_FC0_up_DEGs.txt",sep=""),header=T,sep="\t")
  DEGs<-as.character(up_gene[,1])
  up_DEGs_DApeaks<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)
  cat(paste("Up int gene is ",length(up_DEGs_DApeaks),sep=""),"\n")
  geneu=bitr(up_DEGs_DApeaks,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
  keggu<-enrichKEGG(geneu[,2],organism = 'hsa',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'none',
                    minGSSize = 10,maxGSSize = 500,qvalueCutoff = 1,use_internal_data = FALSE)
  eeu<-keggu@result
  ee1u<-eeu[which(eeu$pvalue<0.05),]
  for(y in 1:nrow(ee1u)){
    ee1u[y,8]
    b7<-matrix(unlist(strsplit(ee1u[y,8],split="/")),ncol=1)
    gene7=bitr(b7,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
    ee1u[y,10]<-paste0(gene7[,2],collapse ="/")
  }
  
  cat(paste("KEGG pathway is",nrow(ee1u),sep=" "),"\n")
  keu<-barplot(keggu,showCategory=10,drop=T,x = "GeneRatio",color = "pvalue")
  #plot(ke)
  ee2u<-ee1u[1:10,]
  enrich_geneu<-ee1u$geneID
  pathway_geneu<-unique(unlist(strsplit(enrich_geneu,split="/")))
  cat(paste("Top pathway gene is ",length(pathway_geneu),sep=""),"\n")
  peaks_DEGs_KEGG[k,4]<-nrow(ee1u)
  peaks_DEGs_KEGG[k,5]<-length(pathway_geneu)
  pathway_gene2u=bitr(pathway_geneu,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo/",file_name,"/3_intersect/pathway",sep=""))
  write.table(ee1d,paste("kegg_pathway_DEGs_DApeaks_",file_name,freq_name,"_up.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  write.table(pathway_gened,paste("kegg_pathway_DEGs_DApeaks_",file_name,freq_name,"_up_gene.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  ggsave(paste("KEGG pathway of",file_name,"-0.05-P-up-P0.05.pdf",sep=""),ked,width = 10, height = 5)
}

colnames(peaks_DEGs_KEGG)<-c("drug","KEGG-down-P0.05-pathway","KEGG-down-P0.05-pathway-gene","KEGG-up-P0.05-pathway","KEGG-up-P0.05-pathway-gene")
setwd("~/xjj/drug/drug_result/durg62")
write.table(peaks_DEGs_KEGG,"peaks_DEGs_KEGG_P0.05_result_42drug.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出




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
#all<-rbind(up_gene,down_gene)
#DEGs<-as.character(all[,1])
down_DEGs_DApeaks<-intersect(Anno_gene_100Kb_down$`Gene Name`,down_gene[,1])
up_DEGs_DApeaks<-intersect(Anno_gene_100Kb_up$`Gene Name`,up_gene[,1])

library(org.Hs.eg.db)
library(clusterProfiler)
up_gene<-down_DEGs_DApeaks  ## 修改
gene=bitr(up_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
## 去重
library(stringr)
kegg<-enrichKEGG(gene[,2],organism = 'hsa',keyType = 'kegg',pvalueCutoff = 0.01,pAdjustMethod = 'none',
                 minGSSize = 10,maxGSSize = 500,qvalueCutoff = 1,use_internal_data = FALSE)
ee<-kegg@result
ee1<-ee[which(ee$pvalue<0.01),]
#ee1<-ee[which(ee$p.adjust<0.05),]

p<-barplot(kegg,showCategory=15,drop=T,x = "GeneRatio",color = "pvalue")
#ee2<-ee1[1:10,]
enrich_gene<-ee1$geneID
pathway_gene<-unique(unlist(strsplit(enrich_gene,split="/")))
a<-unlist(strsplit(enrich_gene,split="/"))
result<-data.frame()
for(j in 1:length(pathway_gene)){
  count<-sum(a==pathway_gene[j])
  result[j,1]<-pathway_gene[j]
  result[j,2]<-count
}
pathway_gene2=bitr(pathway_gene,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
symbol<-pathway_gene2[match(result[,1],pathway_gene2[,1]),]
result1<-cbind(result,symbol)
colnames(result1)<-c("gene","count",colnames(symbol))
result2<-result1[order(result1[,2],decreasing = TRUE),]

a4<-NULL
for(l in 1:nrow(ee1)){
  ee1[l,8]
  a1<-matrix(unlist(strsplit(ee1[l,8],split="/")),ncol=1)
  a2<-matrix(rep(ee1[l,2],length(a1)),ncol=1)
  gene6=bitr(a1,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
  a5<-gene6[match(a1,gene6[,1]),]
  a3<-cbind(a2,a1,a5)
  a4<-rbind(a4,a3)
}

for(y in 1:nrow(ee1)){
  ee1[y,8]
  b1<-matrix(unlist(strsplit(ee1[y,8],split="/")),ncol=1)
  gene6=bitr(b1,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
  ee1[y,10]<-paste0(gene6[,2],collapse ="/")
}

setwd("~/xjj/drug/drug_result/chemotherapy/1_2_intersect/pathway")
kegg_dysregulation<-"kegg_pathway_P0.01_down"
write.table(ee1,paste("DEGP0.05_DApeaksP0.01_",kegg_dysregulation,".txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(pathway_gene,paste("DEGP0.05_DApeaksP0.01_",kegg_dysregulation,"_gene.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
ggsave(paste("DEGP0.05_DApeaksP0.01_",kegg_dysregulation,".pdf",sep=""),p,width = 5, height = 3)
write.table(result2,paste("DEGP0.05_DApeaksP0.01_",kegg_dysregulation,"_gene_symbol.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(a4,paste("DEGP0.05_DApeaksP0.01_",kegg_dysregulation,"_pathway_symbol.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出



