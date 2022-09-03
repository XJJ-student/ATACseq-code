#######  根据药敏数据将样本分为两类
rm(list=ls())
setwd("F:\\Organoid\\药物处理\\杨老师处理过的数据")
ic501<-read.csv("changhai_ic50.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,-(which(colnames(ic501)==c("PC.100","PC.34")))]
x<-NULL
for(i in 1:nrow(ic50)){
  a<-t(ic50[i,])
  b<-a[which(!is.na(a))]
  f<-log2(b+1)
  d<-rep(rownames(ic50)[i],length(f))
  e<-cbind(f,d)
  x<-rbind(x,e)
}  
data<- data.frame(x) 
colnames(data)<-c("value","drug")
range_result<-data.frame()
for(j in 1:nrow(ic50)){
  a<-ic50[j,]
  b<-as.numeric(log2(a[which(!is.na(a))]+1))
  sd<-sd(b)
  range_result[j,1]<-sd*sd
}
range_result[,2]<-rownames(ic50)
a <- range_result[,1]
range_variance<-range_result[order(a,decreasing = T),]
drug_info<-read.csv("compiled_drug_info.csv",header = TRUE,sep = ",")
drug_result2<-drug_info[match(range_variance[,2],drug_info$drug_id),]
drug_result<-cbind(range_variance[,1],drug_result2)
HDAC1<-drug_result[drug_result$target=="HDAC",]
aaa<-c("S1848","S2759","S1194","S1047")
HDAC2<-drug_result[drug_result$drug_id%in%aaa,]
HDAC3<-rbind(HDAC1,HDAC2)
HDAC<-HDAC3[order(HDAC3[,1],decreasing = T),]
HDAC_ic50<-ic50[match(as.character(HDAC[,2]),rownames(ic50)),]
ic2<-log2(HDAC_ic50+1)
library(pheatmap)
heatmap=pheatmap(ic2,scale = "none",main = "log2 IC50",show_rownames=T,show_colnames=T,
                 cluster_rows = TRUE,cluster_cols = TRUE, clustering_distance_rows = "euclidean",
                 clustering_distance_cols = "euclidean",clustering_method = "ward.D2")

######   每种药物在两类样本之间是否有差异
sample1<-c("PC.64","PC.81","PC.116","PC.L","PC.101","PC.27","PC.78","PC.98","PC.130","PC.117","PC.8","PC.109","PC.136","PC.16","PC.97","PC.2","PC.115","PC.13")
sample2<-c("PC.G","PC.111","PC.134","PC.14","PC.40","PC.135","PC.52","PC.22","PC.139","PC.18","PC.105","PC.104","PC.56","PC.119","PC.112","PC.102","PC.121","PC.5","PC.I")
sample1_ic50<-ic2[,match(sample1,colnames(ic2))]
sample2_ic50<-ic2[,match(sample2,colnames(ic2))]

tresult=matrix(0,length(rownames(ic2)),4)
for (i in 1:length(rownames(ic2))){
  a<-sample1_ic50[i,which(!is.na(sample1_ic50[i,]))]
  b<-sample2_ic50[i,which(!is.na(sample2_ic50[i,]))]
  ttest=t.test(a,b)
  tresult[i,1:3]=c(rownames(ic2)[i],ttest$statistic,ttest$p.value)
}
tresult[,4]=p.adjust(tresult[,3],method="BH")
colnames(tresult)<-c("drugid","statistic","p.value","FDR")
sum(as.numeric(tresult[,4])<0.05)
tresult_want<-tresult[as.numeric(tresult[,4])>0.05,]

#####  将没有差异的药物取出来，看看两类样本是否会发生变化
i3<-c("S1848","S7555","S8495")
ic3<-ic2[-(match(i3,rownames(ic2))),]
dn<-matrix(ifelse(is.na(ic3) == T, "NA", ""), nrow(ic3))
heatmap=pheatmap(ic3,scale = "none",main = "log2 IC50",show_rownames=T,show_colnames=T,
                 cluster_rows = TRUE,cluster_cols = TRUE, clustering_distance_rows = "euclidean",
                 clustering_distance_cols = "euclidean",clustering_method = "ward.D2",display_numbers = dn)

Drug1<-c("S1095","S2759","S1096","S1030")
Drug2<-c("S1090","S1194","S1515","S2170","S1085","S2244")
Drug3<-c("S7596","S2693","S1047","S8043","S7569","S2779","S1122")

Drug1_ic50<-ic3[match(Drug1,rownames(ic3)),]
Drug2_ic50<-ic3[match(Drug2,rownames(ic3)),]
Drug3_ic50<-ic3[match(Drug3,rownames(ic3)),]

tresult=matrix(0,length(colnames(ic3)),4)
for (i in 1:length(colnames(ic3))){
  a<-Drug2_ic50[,which(!is.na(Drug2_ic50[,i]))]
  b<-Drug3_ic50[,which(!is.na(Drug3_ic50[,i]))]
  ttest=t.test(a,b)
  tresult[i,1:3]=c(colnames(ic3)[i],ttest$statistic,ttest$p.value)
}
tresult[,4]=p.adjust(tresult[,3],method="BH")
colnames(tresult)<-c("drugid","statistic","p.value","FDR")
sum(as.numeric(tresult[,4])<0.05)

sample11<-c("PC.64","PC.81","PC.116","PC.L","PC.2","PC.115","PC.16","PC.97","PC.101","PC.27","PC.8","PC.109","PC.136","PC.13","PC.78","PC.98","PC.130","PC.117","PC.139","PC.52")  ##20
sample21<-c("PC.G","PC.111","PC.134","PC.14","PC.40","PC.135","PC.22","PC.18","PC.105","PC.104","PC.56","PC.119","PC.112","PC.102","PC.121","PC.5","PC.I")  ##17
sample11_ic50<-ic3[,match(sample11,colnames(ic3))]
sample21_ic50<-ic3[,match(sample21,colnames(ic3))]

tresult=matrix(0,length(rownames(ic3)),4)
for (i in 1:length(rownames(ic3))){
  a<-sample11_ic50[i,which(!is.na(sample11_ic50[i,]))]
  b<-sample21_ic50[i,which(!is.na(sample21_ic50[i,]))]
  ttest=t.test(a,b)
  tresult[i,1:3]=c(rownames(ic3)[i],ttest$statistic,ttest$p.value)
}
tresult[,4]=p.adjust(tresult[,3],method="BH")
colnames(tresult)<-c("drugid","statistic","p.value","FDR")
sum(as.numeric(tresult[,4])<0.05)
tresult_want<-tresult[as.numeric(tresult[,4])>0.05,]  ###剩余的17个药物在两类样本中均差异

###############################################################
###########   clinical information
rm(list=ls())
sensitivity<-c("PC.64","PC.81","PC.116","PC.L","PC.2","PC.115","PC.16","PC.97","PC.101","PC.27","PC.8","PC.109","PC.136","PC.13","PC.78","PC.98","PC.130","PC.117","PC.139","PC.52")  ##20
resistance<-c("PC.G","PC.111","PC.134","PC.14","PC.40","PC.135","PC.22","PC.18","PC.105","PC.104","PC.56","PC.119","PC.112","PC.102","PC.121","PC.5","PC.I")  ##17
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
sensitivity<-c("PC.64","PC.81","PC.116","PC.L","PC.2","PC.115","PC.16","PC.97","PC.101","PC.27","PC.8","PC.109","PC.136","PC.13","PC.78","PC.98","PC.130","PC.117","PC.139","PC.52")  ##20
resistance<-c("PC.G","PC.111","PC.134","PC.14","PC.40","PC.135","PC.22","PC.18","PC.105","PC.104","PC.56","PC.119","PC.112","PC.102","PC.121","PC.5","PC.I")  ##17
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
setwd("~/xjj/WGS_CNV/CNV_FACETS")
library("data.table")
mut<-fread("CNV_FACETS_ANNOVAR_hg19_all_censored_amplification.txt",header=T,data.table=F)
sensitivity<-c("PC.64","PC.81","PC.116","PC.L","PC.2","PC.115","PC.16","PC.97","PC.101","PC.27","PC.8","PC.109","PC.136","PC.13","PC.78","PC.98","PC.130","PC.117","PC.139","PC.52")  ##20
resistance<-c("PC.G","PC.111","PC.134","PC.14","PC.40","PC.135","PC.22","PC.18","PC.105","PC.104","PC.56","PC.119","PC.112","PC.102","PC.121","PC.5","PC.I")  ##17
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
amp_gene<-data_amp[which(data_amp[,2]<0.05),]  ###只有一个基因HYDIN
colnames(amp_gene)<-c("geneid","Pvalue","FDR")

sensitivity1<-apply(sensitivity_mut,1,function(x) mean(x))
resistance1<-apply(resistance_mut,1,function(x) mean(x))
s<-match(amp_gene[,1],geneid)
sum(sensitivity1[s]-resistance1[s]>0)

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
censored_gene<-data_censored[which(data_censored[,2]<0.05),]  ###只有一个基因HYDIN
colnames(censored_gene)<-c("geneid","Pvalue","FDR")
sensitivity1<-apply(sensitivity_mut,1,function(x) mean(x))
resistance1<-apply(resistance_mut,1,function(x) mean(x))
s<-match(censored_gene[,1],geneid)
sum(sensitivity1[s]-resistance1[s]<0)

int_CNV<-intersect(censored_gene[,1],amp_gene[,1])
setwd("~/xjj/drug/drug_result/HDAC_frontiers/7_CNV")
write.table(amp_gene,"sensitivity-resistance-amp-CNV-gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(censored_gene,"sensitivity-resistance-censored-CNV-gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
######### CNV KEGG
library(org.Hs.eg.db)
library(clusterProfiler)
library(stringr)
gene=bitr(amp_gene[,1],fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
kegg<-enrichKEGG(gene[,2],organism = 'hsa',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'none',
                 minGSSize = 10,maxGSSize = 500,qvalueCutoff = 1,use_internal_data = FALSE)
barplot(kegg,showCategory=10,drop=T,x = "GeneRatio",color = "pvalue")

ee<-kegg@result
ee1<-ee[which(ee$pvalue<0.05),]
#p.adjust
enrich_gene<-ee1$geneID
pathway_gene<-unique(unlist(strsplit(enrich_gene,split="/")))
pathway_gene2=bitr(pathway_gene,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
setwd("~/xjj/drug/drug_result/HDAC_frontiers/7_CNV")
write.table(ee1,"kegg_CNV_censored_pathway_P0.05.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(pathway_gene2,"kegg_CNV_censored_pathway_P0.05_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

############    RNAseq DEGs
rm(list=ls())
setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
library("data.table")
exp1<-fread("star_rsem.GeneSymbol.readscount.xls",header=T,data.table=F)
exp2<-floor(exp1[-1])
delete<-apply(exp2,1,function(x) mean(x==0))
cou<-which(delete>0.5)####在超过一半样本以上的基因删去
exp<-exp2[(-cou),]
sensitivity<-c("PC.64","PC.81","PC.116","PC.L","PC.2","PC.115","PC.16","PC.97","PC.101","PC.27","PC.8","PC.109","PC.136","PC.13","PC.78","PC.98","PC.130","PC.117","PC.139","PC.52")  ##20
resistance<-c("PC.G","PC.111","PC.134","PC.14","PC.40","PC.135","PC.22","PC.18","PC.105","PC.104","PC.56","PC.119","PC.112","PC.102","PC.121","PC.5","PC.I")  ##17
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
setwd("~/xjj/drug/drug_result/HDAC_frontiers/1_RNAseq_DEGs")
write.table(res,"sensitivity-resistance-all-DESeq2_HDAC_RNAseq_DEGs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
#resSig<-res[which(res$pvalue<0.01 & abs(res$log2FoldChange>1)),]
resSig<-res[which(res$pvalue<0.05),]
#pvalue padj
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

setwd("~/xjj/drug/drug_result/HDAC_frontiers/1_RNAseq_DEGs")
write.table(up_gene,"sensitivity-resistance-all-DESeq2_RNAseq_FDR0.05_up_DEGs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(down_gene,"sensitivity-resistance-all-DESeq2_RNAseq_FDR0.05_down_DEGs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

#########  RNAseq DEGs Volcano Plot
rm(list=ls())
library(ggrepel) 
library(ggplot2)
setwd("~/xjj/drug/drug_result/HDAC_frontiers/1_RNAseq_DEGs")
df<-read.table("sensitivity-resistance-all-DESeq2_HDAC_RNAseq_DEGs.txt",sep = '\t',header= T,row.names = 1)
#确定是上调还是下调，用于给图中点上色
df$threshold = factor(ifelse(df$pvalue  < 0.05 & abs(df$log2FoldChange) >= 0, ifelse(df$log2FoldChange >= 0 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
df$gene <- row.names(df) #添加一列基因名，以便备注
ggplot(df,aes(x=log2FoldChange,y= -log10(pvalue),color=threshold))+
  geom_point(data = df[df$pvalue<0.05&abs(df$log2FoldChange)>0,],size = 3)+ 
  geom_point(data = df[df$pvalue>0.05|abs(df$log2FoldChange)<0,],size = 3)+
  scale_color_manual(values=c('blue','grey','red'))+#确定点的颜色
  geom_text_repel(
    data = df[df$padj<0.01&abs(df$log2FoldChange)>2,],
    aes(label = gene),
    size = 3,
    color = "black",
    segment.color = "black", show.legend = FALSE )+#添加关注的点的基因名
  ylab('-log10 (pvalue)')+#修改y轴名称
  xlab('log2 (FoldChange)')+#修改x轴名称
  geom_hline(yintercept=-log10(0.05),linetype=4)#添加横线|logFoldChange|>0.25
  #geom_vline(xintercept=c(-2,2),linetype=4)#添加竖线padj<0.05

DEGs<-df[which(df$padj<0.01 & abs(df$log2FoldChange)>2),]
setwd("~/xjj/drug/drug_result/HDAC_frontiers/1_RNAseq_DEGs/Volcano")
write.table(DEGs,"FDR01_FC2_DEGs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出


#################################################################################
################  DEGs与组蛋白修饰酶的交叠
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDAC_frontiers/5_histone-modifying enzymes")
Acetylase<-read.table("histone-modifying-enzymes.txt",header=T,sep="\t")
setwd("~/xjj/drug/drug_result/HDAC_frontiers/1_RNAseq_DEGs")
df<-read.table("sensitivity-resistance-all-DESeq2_HDAC_RNAseq_DEGs.txt",sep = '\t',header= T,row.names = 1)
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
setwd("~/xjj/drug/drug_result/HDAC_frontiers/1_RNAseq_DEGs/Acetylase")
write.table(Acetylase_info,"Acetylase_info_P0.01.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

setwd("~/xjj/drug/drug_result/HDAC_frontiers/1_RNAseq_DEGs/Acetylase")
library(VennDiagram)
venn.diagram(x=list(Acetylase=Acetylase_name,DEGs=DEGs_name),cex = 1,margin = 0.1, "Acetylase_info_P005.png",fill=c("red","blue"))

#################################################################################
################  clusterProfiler
rm(list=ls())
library(org.Hs.eg.db)
library(clusterProfiler)
setwd("~/xjj/drug/drug_result/HDAC_frontiers/1_RNAseq_DEGs")
up_gene<-read.table("sensitivity-resistance-all-DESeq2_RNAseq_P0.05_down_DEGs.txt",header=T,sep="\t")
library(stringr)
gene=bitr(up_gene[,1],fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
kegg<-enrichKEGG(gene[,2],organism = 'hsa',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'none',
                 minGSSize = 10,maxGSSize = 500,qvalueCutoff = 1,use_internal_data = FALSE)
barplot(kegg,showCategory=10,drop=T,x = "GeneRatio",color = "pvalue")

ee<-kegg@result
ee1<-ee[which(ee$p.adjust<0.05),]
#p.adjust
#dotplot(kegg,showCategory=8,x = "GeneRatio",color = "p.adjust",title = "DEGs-P05-KEGG")
#enrichplot::gseaplot2(ee1,1,pvalue_table=T,color="#086538")
enrich_gene<-ee1$geneID
pathway_gene<-unique(unlist(strsplit(enrich_gene,split="/")))
pathway_gene2=bitr(pathway_gene,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
setwd("~/xjj/drug/drug_result/HDAC_frontiers/1_RNAseq_DEGs/pathway")
write.table(ee1,"kegg_DEG_pathway_FDR0.05_down_drug_target.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(pathway_gene2,"kegg_DEG_pathway_FDR0.05_down_drug_target_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
#dotplot(kegg,showCategory=20)
go<-enrichGO(gene[,2],OrgDb = org.Hs.eg.db,ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.01,keyType = 'ENTREZID')
dotplot(go,showCategory=10)
a<-go@result
go_BP<-a[which(a$p.adjust<0.01),]

enrich_genego<-go_BP$geneID
pathway_genego<-unique(unlist(strsplit(enrich_genego,split="/")))
pathway_genego2=bitr(pathway_genego,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
setwd("~/xjj/drug/drug_result/HDAC_frontiers/1_RNAseq_DEGs/pathway")
write.table(pathway_genego2,"GO_DEG_pathway_FDR0.01_down_drug_target_gene.txt",sep = '\t',col.names = F,row.names = F,quote = FALSE)##数据输出
write.table(go_BP,"GO_DEG_pathway_FDR0.01_down_drug_target.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

#######################    ATACseq   ################################################
#############  DApeaks
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_count_name<-read.table("peak_count.txt",sep="\t",header=T)
sensitivity<-c("PC.64","PC.81","PC.116","PC.L","PC.2","PC.115","PC.16","PC.97","PC.101","PC.27","PC.8","PC.109","PC.136","PC.13","PC.78","PC.98","PC.130","PC.117","PC.139","PC.52")  ##20
resistance<-c("PC.G","PC.111","PC.134","PC.14","PC.40","PC.135","PC.22","PC.18","PC.105","PC.104","PC.56","PC.119","PC.112","PC.102","PC.121","PC.5","PC.I")  ##17
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
dds##这是一个关于各种内容的矩阵
dds<-DESeq(dds)##进行标准化分析
sizeFactors(dds)##查看每个主成分的标准化值
res<-results(dds)##将结果输出
head(res)
class(res)##可以看出它是DESeq的属性，要转化为表格
res<-as.data.frame(res)
head(res)
res<-cbind(peak_count_name[,1],res)  ##对数据增加一列
head(res)
colnames(res)<- c('peak_id',"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj" )##对列名进行修改
setwd("~/xjj/drug/drug_result/HDAC_frontiers/1_RNAseq_DEGs")
write.table(res,"sensitivity-resistance-all-DESeq2_heatmap13HDAC_ATAC.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
#setwd("~/xjj/drug/drug_result/HDAC_drug/heatmap")
#res<-read.table("sensitivity-resistance-all-DESeq2_HDAC_RNAseq_DEGs.txt",header = TRUE,sep = "\t")
#resSig<-res[which(res$padj<0.05 & abs(res$log2FoldChange>1)),]
resSig<-res[which(res$padj<0.05),]
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
out_file<-"sen-res-HDAC13-heatmap-0.01-P"
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
#######################################################################################
####  homer注释完的信息,画饼图
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDAC_frontiers/2_ATACseq_DApeaks")
library(data.table)
homer_anno<-fread("sen-res-HDAC13-heatmap-0.05-P-up-annotation.txt",header=T,data.table=F)
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

#####################  符合在TSS100KB之内条件的注释基因
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDAC_frontiers/2_ATACseq_DApeaks")
library(data.table)
homer_anno_down<-fread("sen-res-HDAC13-heatmap-0.05-P-down-annotation.txt",header=T,data.table=F)
colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]

homer_anno_up<-fread("sen-res-HDAC13-heatmap-0.05-P-up-annotation.txt",header=T,data.table=F)
colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]

setwd("~/xjj/drug/drug_result/HDAC_frontiers/1_RNAseq_DEGs")
up_gene<-read.table("sensitivity-resistance-all-DESeq2_RNAseq_P0.05_up_DEGs.txt",header=T,sep="\t")
down_gene<-read.table("sensitivity-resistance-all-DESeq2_RNAseq_P0.05_down_DEGs.txt",header=T,sep="\t")
all<-rbind(up_gene,down_gene)
DEGs<-as.character(all[,1])
down_DEGs_DApeaks<-intersect(Anno_gene_100Kb_down$`Gene Name`,DEGs)
up_DEGs_DApeaks<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)

library(gplots)
library(VennDiagram)
down_peaks_gene <- Anno_gene_100Kb_down$`Gene Name`
up_peaks_gene <- Anno_gene_100Kb_up$`Gene Name`
input  <-list(unique(down_peaks_gene),unique(up_peaks_gene),unique(DEGs))
venn(input,showSetLogicLabel=TRUE)
tmp <- venn(input)
int<-attr(tmp, "intersections")

setwd("~/xjj/drug/drug_result/HDAC_frontiers/3_intersect/venn")
library(VennDiagram)
venn.diagram(x=list(downpeaks=unique(down_peaks_gene),uppeaks=unique(up_peaks_gene),RNAseq=unique(DEGs)), "P005.png",fill=c("red","green","blue"),margin = 0.1)
write.table(a,"DApeaks_DEGs_p05.txt",col.names=F,row.names=F)
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
setwd("~/xjj/drug/drug_result/HDAC_frontiers/2_ATACseq_DApeaks")
library(data.table)
homer_anno_down<-fread("sen-res-HDAC13-heatmap-0.05-P-down-annotation.txt",header=T,data.table=F)
colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]

homer_anno_up<-fread("sen-res-HDAC13-heatmap-0.05-P-up-annotation.txt",header=T,data.table=F)
colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]

setwd("~/xjj/drug/drug_result/HDAC_frontiers/1_RNAseq_DEGs")
up_gene<-read.table("sensitivity-resistance-all-DESeq2_RNAseq_P0.05_up_DEGs.txt",header=T,sep="\t")
down_gene<-read.table("sensitivity-resistance-all-DESeq2_RNAseq_P0.05_down_DEGs.txt",header=T,sep="\t")
all<-rbind(down_gene,up_gene)
DEGs<-as.character(all[,1])
down_DEGs_DApeaks<-intersect(Anno_gene_100Kb_down$`Gene Name`,DEGs)
up_DEGs_DApeaks<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)
int_all_DEGs<-unique(c(down_DEGs_DApeaks,up_DEGs_DApeaks))
all_peaks_anno<-rbind(homer_anno_down,homer_anno_up)
#DEGs_peaks<-all_peaks_anno[all_peaks_anno$`Gene Name` %in% int_all_DEGs,]
DEGs_peaks<-all_peaks_anno[match(int_all_DEGs,all_peaks_anno$`Gene Name`),1]

setwd("~/xjj/drug/drug_result/HDAC_frontiers/2_ATACseq_DApeaks")
DApeak<-fread("sensitivity-resistance-all-DESeq2_heatmap13HDAC_ATAC.txt",header=T,data.table=F)
DEGs_peaks_FC<-DApeak[match(DEGs_peaks,DApeak$peak_id),3]
DEGs_FC<-all[match(int_all_DEGs,DEGs),3]
gene<-as.character(all[match(int_all_DEGs,DEGs),1])
#pearson", "kendall", "spearman
cor_result<-cor.test(DEGs_peaks_FC, DEGs_FC,alternative = "two.sided",method = "pearson")
cor_result$estimate
cor_result$p.value

################################   画相关性分析图
a1<-matrix(DEGs_FC,ncol=1)
a2<-matrix(DEGs_peaks_FC,ncol=1)
dat<-cbind(a1,a2)
dat1<-as.data.frame(dat)
rownames(dat1)<-gene
colnames(dat1)<-c("DEGslog2FC","DApeakslog2FC")
library(ggplot2)
library(ggpubr)
a11<-c(1:10)
a21<-c(1204:1213)
a3<-c(a11,a21)
want_dat<-dat1[order(dat1[,1])[a3],]
ggplot(data=dat1, aes(x=DEGslog2FC, y=DApeakslog2FC))+geom_point(color="red",size=0.1)+stat_smooth(method="lm",se=FALSE)+stat_cor(data=dat1, method = "spearman")+
       ggtitle("spearman-P0.05") +
       theme(plot.title = element_text(hjust = 0.5))+
       geom_text_repel(
       data = want_dat[,c(1:2)],
       aes(label = rownames(want_dat)),
       size = 3,
       color = "black",
       segment.color = "black", show.legend = FALSE )

#pearson", "kendall", "spearman
#############################  对高低可及peaks-DEGs进行KEGG通路富集分析
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDAC_frontiers/2_ATACseq_DApeaks")
library(data.table)
homer_anno_down<-fread("sen-res-HDAC13-heatmap-0.05-P-down-annotation.txt",header=T,data.table=F)
colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]

homer_anno_up<-fread("sen-res-HDAC13-heatmap-0.05-P-up-annotation.txt",header=T,data.table=F)
colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]

setwd("~/xjj/drug/drug_result/HDAC_frontiers/1_RNAseq_DEGs")
up_gene<-read.table("sensitivity-resistance-all-DESeq2_RNAseq_P0.05_up_DEGs.txt",header=T,sep="\t")
down_gene<-read.table("sensitivity-resistance-all-DESeq2_RNAseq_P0.05_down_DEGs.txt",header=T,sep="\t")
all<-rbind(down_gene,up_gene)
DEGs<-as.character(all[,1])
down_DEGs_DApeaks<-intersect(Anno_gene_100Kb_down$`Gene Name`,DEGs)
up_DEGs_DApeaks<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)

library(org.Hs.eg.db)
library(clusterProfiler)
up_gene<-down_DEGs_DApeaks  ## 修改
gene=bitr(up_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
## 去重
library(stringr)
kegg<-enrichKEGG(gene[,2],organism = 'hsa',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'none',
                 minGSSize = 10,maxGSSize = 500,qvalueCutoff = 1,use_internal_data = FALSE)
ee<-kegg@result
ee1<-ee[which(ee$pvalue<0.05),]
ee1<-ee[which(ee$p.adjust<0.05),]

barplot(kegg,showCategory=20,drop=T,x = "GeneRatio",color = "pvalue")
ee2<-ee1[1:20,]
enrich_gene<-ee1$geneID
pathway_gene<-unique(unlist(strsplit(enrich_gene,split="/")))
pathway_gene2=bitr(pathway_gene,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
setwd("~/xjj/drug/drug_result/HDAC_frontiers/3_intersect/pathway")
write.table(ee1,"kegg_pathway_DEGs_DApeaks_P0.05_down.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(pathway_gene,"kegg_pathway_DEGs_DApeaks_P0.05_down_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

#############################################################################
######################  使用homer进行motif富集，并对应到TFs
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDAC_frontiers/4_TFs/homer_result/sen-res-HDAC13-heatmap-0.05-P-logFC0-DESeq2-up")
knownResults<-read.delim("knownResults.txt",header=T,sep='\t')
FDR05<-knownResults[knownResults$q.value..Benjamini.<0.05,]
#FDR05<-knownResults[knownResults$P.value<0.01,]
library(dplyr)
homer_TF1<-data.frame(lapply(strsplit(as.character(FDR05[,1]),'/'), function(x) x[3])%>%unlist())
FDR051<-FDR05[which(homer_TF1=="Homer"),]
homer_TF2<-data.frame(lapply(strsplit(as.character(FDR051[,1]),'/'), function(x) x[1])%>%unlist())
library(stringr)
homer_TF3<-apply(homer_TF2,1,function(x) str_replace_all(as.character(x),"\\((.*?)\\)",""))##取出括号前的字符
homer_TF4<-unique(homer_TF3)
length(homer_TF4)
new_homer_result<-cbind(homer_TF3,FDR051[,-1])
colnames(new_homer_result)<-c("TF_name",colnames(new_homer_result)[-1])
write.table(new_homer_result,"DApeaks_P0.05-logFC0-DESeq2-up_FDR0.05_TFs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

#####  画柱状图
data1<-cbind(as.data.frame(new_homer_result[,1]),new_homer_result$P.value)
data1[,2]<-(-log2(data1[,2]))
colnames(data1)<-c("TFname","Pvalue")
data11<-data1[match(unique(data1[,1]),data1[,1]),]##重叠的TFqu取第一个
data2<-data11[1:10,]
data2$TFname=factor(data2$TFname,levels = data2[,1])
new_homer_result1<-new_homer_result[match(unique(data1[,1]),data1[,1]),]
val<-new_homer_result1$X..of.Target.Sequences.with.Motif[1:10]
library(ggplot2)
ggplot(data=data2,mapping=aes(x=TFname,y=as.numeric(Pvalue)))+
  geom_bar(stat="identity",fill="blue4")+
  #geom_text(aes(label = val, vjust = -0.8, hjust = 0.5), show.legend = TRUE)+
  ylab("-log2(Pvalue)")+
  ggtitle("DApeaks_P0.05-logFC0-DESeq2-up_FDR0.05_TFs") +
  theme(plot.title = element_text(hjust = 0.5))


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
setwd("~/xjj/drug/drug_result/HDAC_frontiers/1_RNAseq_DEGs")
up_gene<-read.table("sensitivity-resistance-all-DESeq2_RNAseq_P0.05_down_DEGs.txt",header=T,sep="\t")

setwd("~/xjj/drug/drug_result/HDAC_frontiers/4_TFs/homer_result/sen-res-HDAC13-heatmap-0.05-P-logFC0-DESeq2-down")
knownResults<-read.delim("knownResults.txt",header=T,sep='\t')
#FDR05<-knownResults[knownResults$P.value<0.05,]
FDR05<-knownResults[knownResults$q.value..Benjamini.<0.05,]

library(dplyr)
homer_TF1<-data.frame(lapply(strsplit(as.character(FDR05[,1]),'/'), function(x) x[3])%>%unlist())
FDR051<-FDR05[which(homer_TF1=="Homer"),]
homer_TF2<-data.frame(lapply(strsplit(as.character(FDR051[,1]),'/'), function(x) x[1])%>%unlist())
library(stringr)
homer_TF3<-apply(homer_TF2,1,function(x) str_replace_all(as.character(x),"\\((.*?)\\)",""))##取出括号前的字符
homer_TF4<-unique(homer_TF3)
length(homer_TF4)
new_homer_result<-cbind(homer_TF3,FDR051[,-1])
colnames(new_homer_result)<-c("TF_name",colnames(new_homer_result)[-1])
length(intersect(toupper(new_homer_result$TF_name),toupper(up_gene$gene_id)))
DETFs<-intersect(toupper(new_homer_result$TF_name),toupper(up_gene$gene_id))

setwd("~/xjj/drug/drug_result/HDAC_frontiers/3_intersect/DETFs_venn")
library(VennDiagram)
venn.diagram(x=list(TF_enrich_hypo_DApeaks=toupper(new_homer_result$TF_name),down_DEGs=toupper(up_gene$gene_id)), "down_DETFs_P005.png",fill=c("red","blue"),margin = 0.3)

###############   TF找到对应的靶基因
setwd("~/xjj/drug/drug_result/HDAC_frontiers/4_TFs/TRANSFAC_TF")
target_TF<-read.delim("TRANSFAC-gene_attribute_edges.txt",header=T,sep='\t')
target_TF<-read.delim("JASPAR-gene_attribute_edges.txt",header=T,sep='\t')
target_TF<-read.delim("MotifMap-gene_attribute_edges.txt",header=T,sep='\t')
target_TF<-read.delim("ENCODE-gene_attribute_edges.txt",header=T,sep='\t')

target_TF_gene<-target_TF[,c(1,4)]
target_TF_gene<-target_TF_gene[-1,]
colnames(target_TF_gene)<-c("Target","TFs")
intersect(DETFs,unique(target_TF_gene[,2]))
DETF_want<-intersect(DETFs,unique(target_TF_gene[,2]))
result<-NULL
for(i in 1:length(DETF_want)){
  result1<-target_TF_gene[target_TF_gene[,2] %in% DETF_want[i],]
  result<-rbind(result,result1)
}
### 查看靶标基因中有多少是DEGs
setwd("~/xjj/drug/drug_result/HDAC_frontiers/1_RNAseq_DEGs")
up_gene<-read.table("sensitivity-resistance-all-DESeq2_RNAseq_P0.05_up_DEGs.txt",header=T,sep="\t")
down_gene<-read.table("sensitivity-resistance-all-DESeq2_RNAseq_P0.05_down_DEGs.txt",header=T,sep="\t")

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
setwd("~/xjj/drug/drug_result/HDAC_frontiers/4_TFs/TRANSFAC_TF/DETFs-target-DEGs")
write.table(data,"P0.01_down_DEGs_DApeaks_FDR0.05_DETFs_TRANSFAC-gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(data,"P0.01_down_DEGs_DApeaks_FDR0.05_DETFs_JASPAR-gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(data,"P0.01_down_DEGs_DApeaks_FDR0.05_DETFs_MotifMap-gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(data,"P0.01_down_DEGs_DApeaks_FDR0.05_DETFs_ENCODE-gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

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











