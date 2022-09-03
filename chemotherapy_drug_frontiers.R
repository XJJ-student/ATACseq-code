#######  根据药敏数据将样本分为两类
rm(list=ls())
setwd("~/xjj/drug")
ic501<-read.csv("changhai_ic50.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,-(which(colnames(ic501)==c("PC.100","PC.34")))]
HDAC_ic50<-ic50[c(1:5),]

ic2<-log2(HDAC_ic50+1)
library(pheatmap)
heatmap=pheatmap(ic2,scale = "none",main = "chemotherapy log2 IC50",show_rownames=T,show_colnames=T,
                 cluster_rows = TRUE,cluster_cols = TRUE, clustering_distance_rows = "euclidean",
                 clustering_distance_cols = "euclidean",clustering_method = "ward.D")

######   每种药物在两类样本之间是否有差异
sample1<-c("PC.18","PC.81","PC.98","PC.115","PC.78","PC.135","PC.G","PC.13","PC.139","PC.117","PC.14","PC.101","PC.40","PC.I","PC.27","PC.16","PC.97","PC.52","PC.2","PC.109","PC.8")
sample2<-c("PC.56","PC.105","PC.112","PC.L","PC.104","PC.136","PC.130","PC.64","PC.111","PC.119","PC.116","PC.121","PC.134","PC.22","PC.102","PC.5")
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
ic3<-ic2[1:3,]
dn<-matrix(ifelse(is.na(ic3) == T, "NA", ""), nrow(ic3))
heatmap=pheatmap(ic3,scale ="row",main = "log2 IC50",show_rownames=T,show_colnames=T,
                 cluster_rows = TRUE,cluster_cols = TRUE, clustering_distance_rows = "euclidean",
                 clustering_distance_cols = "euclidean",clustering_method = "centroid")

###############################################################
#############   mutation landscape
rm(list=ls())
setwd("~/xjj/WGS_CNV/mutation")
library("data.table")
mut<-fread("Organoid_mut_binary.txt",header=T,data.table=F)
mut1<-mut[,-1]
sensitivity<-c("PC.18","PC.81","PC.98","PC.115","PC.78","PC.135","PC.G","PC.13","PC.139","PC.117","PC.14","PC.101","PC.40","PC.I","PC.27","PC.16","PC.97","PC.52","PC.2","PC.109","PC.8")
resistance<-c("PC.56","PC.105","PC.112","PC.L","PC.104","PC.136","PC.130","PC.64","PC.111","PC.119","PC.116","PC.121","PC.134","PC.22","PC.102","PC.5")
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
resistance_new<-as.character(sample_id[match(resistance,sample_id[,2]),1])
sensitivity_new<-as.character(sample_id[match(sensitivity,sample_id[,2]),1])
res_int<-intersect(resistance_new,colnames(mut))
resistance_mut<-mut[,match(res_int,colnames(mut))]
sen_int<-intersect(sensitivity_new,colnames(mut))
sensitivity_mut<-mut[,match(sen_int,colnames(mut))]
geneid<-mut[,1]
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
which(data[,2]<0.05)  #none



############    RNAseq DEGs
rm(list=ls())
setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
library("data.table")
exp1<-fread("star_rsem.GeneSymbol.readscount.xls",header=T,data.table=F)
exp2<-floor(exp1[-1])
delete<-apply(exp2,1,function(x) mean(x==0))
cou<-which(delete>0.5)####在超过一半样本以上的基因删去
exp<-exp2[(-cou),]
sensitivity<-c("PC.18","PC.81","PC.98","PC.115","PC.78","PC.135","PC.G","PC.13","PC.139","PC.117","PC.14","PC.101","PC.40","PC.I","PC.27","PC.16","PC.97","PC.52","PC.2","PC.109","PC.8")
resistance<-c("PC.56","PC.105","PC.112","PC.L","PC.104","PC.136","PC.130","PC.64","PC.111","PC.119","PC.116","PC.121","PC.134","PC.22","PC.102","PC.5")
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
setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/1_RNAseq_DEGs")
write.table(res,"sensitivity-resistance-all-DESeq2_chemotherapy_RNAseq_DEGs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
#resSig<-res[which(res$padj<0.05 & abs(res$log2FoldChange>1)),]
resSig<-res[which(res$pvalue<0.01),]
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

setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/1_RNAseq_DEGs")
write.table(up_gene,"sensitivity-resistance-all-DESeq2_RNAseq_P0.01_up_DEGs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(down_gene,"sensitivity-resistance-all-DESeq2_RNAseq_P0.01_down_DEGs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

#########  RNAseq DEGs Volcano Plot
rm(list=ls())
library(ggrepel) 
library(ggplot2)
setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/1_RNAseq_DEGs")
df<-read.table("sensitivity-resistance-all-DESeq2_chemotherapy_RNAseq_DEGs.txt",sep = '\t',header= T,row.names = 1)
#确定是上调还是下调，用于给图中点上色
df$threshold = factor(ifelse(df$padj  < 0.05 & abs(df$log2FoldChange) >= 2, ifelse(df$log2FoldChange >= 2 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
df$gene <- row.names(df) #添加一列基因名，以便备注
ggplot(df,aes(x=log2FoldChange,y= -log10(padj),color=threshold))+
  geom_point(data = df[df$padj<0.05&abs(df$log2FoldChange)>2,],size = 3)+ 
  geom_point(data = df[df$padj>0.05|abs(df$log2FoldChange)<2,],size = 3)+
  scale_color_manual(values=c('blue','grey','red'))+#确定点的颜色
  geom_text_repel(
    data = df[df$padj<0.05&abs(df$log2FoldChange)>2,],
    aes(label = gene),
    size = 3,
    color = "black",
    segment.color = "black", show.legend = FALSE )+#添加关注的点的基因名
  ylab('-log10 (padj)')+#修改y轴名称
  xlab('log2 (FoldChange)')+#修改x轴名称
  geom_hline(yintercept=-log10(0.05),linetype=4)+#添加横线|logFoldChange|>0.25
  geom_vline(xintercept=c(-2,2),linetype=4)#添加竖线padj<0.05

DEGs<-df[which(df$padj<0.05 & abs(df$log2FoldChange)>2),]
setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/1_RNAseq_DEGs/Volcano")
write.table(DEGs,"FDR05_FC2_DEGs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出


#################################################################################
################  DEGs与组蛋白修饰酶的交叠  
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDAC_frontiers/5_histone-modifying enzymes")
Acetylase<-read.table("histone-modifying-enzymes.txt",header=T,sep="\t")
setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/1_RNAseq_DEGs")
df<-read.table("sensitivity-resistance-all-DESeq2_chemotherapy_RNAseq_DEGs.txt",sep = '\t',header= T,row.names = 1)
DEGs<-df[df$pvalue<0.01,]
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
setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/1_RNAseq_DEGs/Acetylase")
write.table(Acetylase_info,"Acetylase_info_P0.01.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出


#################################################################################
################  clusterProfiler
rm(list=ls())
library(org.Hs.eg.db)
library(clusterProfiler)
setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/1_RNAseq_DEGs")
up_gene<-read.table("sensitivity-resistance-all-DESeq2_RNAseq_P0.05_down_DEGs.txt",header=T,sep="\t")
library(stringr)
gene=bitr(up_gene[,1],fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
kegg<-enrichKEGG(gene[,2],organism = 'hsa',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'none',
                 minGSSize = 10,maxGSSize = 500,qvalueCutoff = 1,use_internal_data = FALSE)
barplot(kegg,showCategory=10,drop=T,x = "GeneRatio",color = "qvalue")

ee<-kegg@result
ee1<-ee[which(ee$pvalue<0.05),]
#p.adjust
#dotplot(kegg,showCategory=8,x = "GeneRatio",color = "p.adjust",title = "DEGs-P05-KEGG")
#enrichplot::gseaplot2(ee1,1,pvalue_table=T,color="#086538")
enrich_gene<-ee1$geneID
pathway_gene<-unique(unlist(strsplit(enrich_gene,split="/")))
pathway_gene2=bitr(pathway_gene,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/1_RNAseq_DEGs/pathway")
write.table(ee1,"kegg_DEG_pathway_P0.05_down_drug_target.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(pathway_gene2,"kegg_DEG_pathway_P0.05_down_drug_target_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
#dotplot(kegg,showCategory=20)
go<-enrichGO(gene[,2],OrgDb = org.Hs.eg.db,ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.05,keyType = 'ENTREZID')
dotplot(go,showCategory=10)
a<-go@result
go_BP<-a[which(a$p.adjust<0.05),]

enrich_genego<-go_BP$geneID
pathway_genego<-unique(unlist(strsplit(enrich_genego,split="/")))
pathway_genego2=bitr(pathway_genego,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/1_RNAseq_DEGs/pathway")
write.table(pathway_genego2,"GO_DEG_pathway_P0.05_down_drug_target_gene.txt",sep = '\t',col.names = F,row.names = F,quote = FALSE)##数据输出
write.table(go_BP,"GO_DEG_pathway_P0.05_down_drug_target.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

#######################    ATACseq   ################################################
#############  DApeaks
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_count_name<-read.table("peak_count.txt",sep="\t",header=T)
sensitivity<-c("PC.18","PC.81","PC.98","PC.115","PC.78","PC.135","PC.G","PC.13","PC.139","PC.117","PC.14","PC.101","PC.40","PC.I","PC.27","PC.16","PC.97","PC.52","PC.2","PC.109","PC.8")
resistance<-c("PC.56","PC.105","PC.112","PC.L","PC.104","PC.136","PC.130","PC.64","PC.111","PC.119","PC.116","PC.121","PC.134","PC.22","PC.102","PC.5")
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
setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/2_ATACseq_DApeaks")
write.table(res,"sensitivity-resistance-all-DESeq2_chemotherapy_ATAC.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
#setwd("~/xjj/drug/drug_result/HDAC_drug/heatmap")
#res<-read.table("resistance-sensitivity-all-DESeq2_heatmap.txt",header = TRUE,sep = "\t")
#resSig<-res[which(res$padj<0.05 & abs(res$log2FoldChange>1)),]
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
library(dplyr)
a_up<-apply(up_gene,1,function(x) unlist(strsplit(as.character(x), "_")))
d_up<-t(a_up)
e_up<-data.frame(rep("+",nrow(d_up)))
peaks_up<-cbind(up_gene,d_up,e_up)
colnames(peaks_up)<-c("peak_id","chr","start","end","strand")
out_file<-"sen-res-chemotherapy-heatmap-0.01-P"
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
setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/2_ATACseq_DApeaks")
library(data.table)
homer_anno<-fread("sen-res-chemotherapy-0.01-P-down-annotation.txt",header=T,data.table=F)
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
setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/2_ATACseq_DApeaks")
library(data.table)
homer_anno_down<-fread("sen-res-chemotherapy-0.01-P-down-annotation.txt",header=T,data.table=F)
colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]

homer_anno_up<-fread("sen-res-chemotherapy-0.01-P-up-annotation.txt",header=T,data.table=F)
colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]

setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/1_RNAseq_DEGs")
up_gene<-read.table("sensitivity-resistance-all-DESeq2_RNAseq_P0.01_up_DEGs.txt",header=T,sep="\t")
down_gene<-read.table("sensitivity-resistance-all-DESeq2_RNAseq_P0.01_down_DEGs.txt",header=T,sep="\t")
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

setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/3_intersect/venn")
library(VennDiagram)
venn.diagram(x=list(downpeaks=unique(down_peaks_gene),uppeaks=unique(up_peaks_gene),RNAseq=unique(DEGs)), "P001.png",fill=c("red","green","blue"))

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

sensitivity<-c("PC.18","PC.81","PC.98","PC.115","PC.78","PC.135","PC.G","PC.13","PC.139","PC.117","PC.14","PC.101","PC.40","PC.I","PC.27","PC.16","PC.97","PC.52","PC.2","PC.109","PC.8")
resistance<-c("PC.56","PC.105","PC.112","PC.L","PC.104","PC.136","PC.130","PC.64","PC.111","PC.119","PC.116","PC.121","PC.134","PC.22","PC.102","PC.5")
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

##################################################################################
##############  DEpeaks在不同样本变化的倍数与DEGs在不同样本变化的倍数之间的相关性
rm(list=ls())
setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/2_ATACseq_DApeaks")
library(data.table)
homer_anno_down<-fread("sen-res-chemotherapy-0.01-P-down-annotation.txt",header=T,data.table=F)
colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]

homer_anno_up<-fread("sen-res-chemotherapy-0.01-P-up-annotation.txt",header=T,data.table=F)
colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]

setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/1_RNAseq_DEGs")
up_gene<-read.table("sensitivity-resistance-all-DESeq2_RNAseq_P0.01_up_DEGs.txt",header=T,sep="\t")
down_gene<-read.table("sensitivity-resistance-all-DESeq2_RNAseq_P0.01_down_DEGs.txt",header=T,sep="\t")
all<-rbind(down_gene,up_gene)
DEGs<-as.character(all[,1])
down_DEGs_DApeaks<-intersect(Anno_gene_100Kb_down$`Gene Name`,DEGs)
up_DEGs_DApeaks<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)
int_all_DEGs<-unique(c(down_DEGs_DApeaks,up_DEGs_DApeaks))##共交叠的基因个数
all_peaks_anno<-rbind(homer_anno_down,homer_anno_up)
#DEGs_peaks<-all_peaks_anno[all_peaks_anno$`Gene Name` %in% int_all_DEGs,]
DEGs_peaks<-all_peaks_anno[match(int_all_DEGs,all_peaks_anno$`Gene Name`),1]

setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/2_ATACseq_DApeaks")
DApeak<-fread("sensitivity-resistance-all-DESeq2_chemotherapy_ATAC.txt",header=T,data.table=F)
DEGs_peaks_FC<-DApeak[match(DEGs_peaks,DApeak$peak_id),3]
DEGs_FC<-all[match(int_all_DEGs,DEGs),3]

#pearson", "kendall", "spearman
cor_result<-cor.test(DEGs_peaks_FC, DEGs_FC,alternative = "two.sided",method = "spearman")
cor_result$estimate
cor_result$p.value

################################   画相关性分析图
a1<-matrix(DEGs_FC,ncol=1)
a2<-matrix(DEGs_peaks_FC,ncol=1)
dat<-cbind(a1,a2)
dat1<-as.data.frame(dat)
colnames(dat1)<-c("DEGslog2FC","DApeakslog2FC")
library(ggplot2)
library(ggpubr)
ggplot(data=dat1, aes(x=DEGslog2FC, y=DApeakslog2FC))+geom_point(color="red")+stat_smooth(method="lm",se=FALSE)+stat_cor(data=dat1, method = "spearman")+
       ggtitle("spearman-P0.01") +
       theme(plot.title = element_text(hjust = 0.5))

#pearson", "kendall", "spearman
#############################  对高低可及peaks-DEGs进行KEGG通路富集分析
rm(list=ls())
setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/2_ATACseq_DApeaks")
library(data.table)
homer_anno_down<-fread("sen-res-chemotherapy-0.01-P-down-annotation.txt",header=T,data.table=F)
colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]

homer_anno_up<-fread("sen-res-chemotherapy-0.01-P-up-annotation.txt",header=T,data.table=F)
colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]

setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/1_RNAseq_DEGs")
up_gene<-read.table("sensitivity-resistance-all-DESeq2_RNAseq_P0.01_up_DEGs.txt",header=T,sep="\t")
down_gene<-read.table("sensitivity-resistance-all-DESeq2_RNAseq_P0.01_down_DEGs.txt",header=T,sep="\t")
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
setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/3_intersect/pathway")
write.table(ee1,"kegg_pathway_DEGs_DApeaks_P0.01_down.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(pathway_gene,"kegg_pathway_DEGs_DApeaks_P0.01_down_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

#############################################################################
######################  使用homer进行motif富集，并对应到TFs
rm(list=ls())
setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/2_ATACseq_DApeaks/homer_result/sen-res-chemotherapy-heatmap-0.01-P-logFC0-DESeq2-down")
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
setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/4_TFs/homer_result")
write.table(new_homer_result,"DApeaks_P0.01-logFC0-DESeq2-down_FDR0.05_TFs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

#####  画柱状图
data1<-cbind(as.data.frame(new_homer_result[,1]),new_homer_result$P.value)
data1[,2]<-(-log2(data1[,2]))
colnames(data1)<-c("TFname","Pvalue")
data11<-data1[match(unique(data1[,1]),data1[,1]),]##重叠的TFqu取第一个
data2<-data11[1:20,]
data2$TFname=factor(data2$TFname,levels = data2[,1])
new_homer_result1<-new_homer_result[match(unique(data1[,1]),data1[,1]),]
val<-new_homer_result1$X..of.Target.Sequences.with.Motif[1:20]
library(ggplot2)
ggplot(data=data2,mapping=aes(x=TFname,y=as.numeric(Pvalue)))+
  geom_bar(stat="identity")+
  geom_text(aes(label = val, vjust = -0.8, hjust = 0.5), show.legend = TRUE)+
  ylab("-log2(Pvalue)")+
  ggtitle("DApeaks_P0.01-logFC0-DESeq2-down_FDR0.05_TFs") +
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
setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/1_RNAseq_DEGs")
up_gene<-read.table("sensitivity-resistance-all-DESeq2_RNAseq_P0.05_down_DEGs.txt",header=T,sep="\t")

setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/2_ATACseq_DApeaks/homer_result/sen-res-chemotherapy-heatmap-0.05-P-logFC0-DESeq2-down")
knownResults<-read.delim("knownResults.txt",header=T,sep='\t')
FDR05<-knownResults[knownResults$P.value<0.05,]
#FDR05<-knownResults[knownResults$q.value..Benjamini.<0.05,]

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
setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/1_RNAseq_DEGs")
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
setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/4_TFs/TRANSFAC_TF/DETFs-target-DEGs")
write.table(data,"P0.05_down_DEGs_DApeaks_FDR0.05_DETFs_TRANSFAC-gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(data,"P0.05_down_DEGs_DApeaks_FDR0.05_DETFs_JASPAR-gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(data,"P0.05_down_DEGs_DApeaks_FDR0.05_DETFs_MotifMap-gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(data,"P0.05_down_DEGs_DApeaks_FDR0.05_DETFs_ENCODE-gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出


#### DETFs veen
library(gplots)
library(VennDiagram)
TF<- unique(toupper(new_homer_result$TF_name))
DEGs_name <- unique(toupper(up_gene$gene_id))
input  <-list(TF,DEGs_name)
venn(input,showSetLogicLabel=TRUE)
venn(input)

###################################################################################
################   分类器 randomForest
rm(list=ls())
setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
library("data.table")
exp1<-fread("star_rsem.GeneSymbol.TPM.xls",header=T,data.table=F)
exp2<-log2(exp1[,-1]+1)
delete<-apply(exp2,1,function(x) mean(x==0))
cou<-which(delete>0.5)####在超过一半样本以上的基因删去
exp<-exp2[-cou,]
sensitivity<-c("PC.18","PC.81","PC.98","PC.115","PC.78","PC.135","PC.G","PC.13","PC.139","PC.117","PC.14","PC.101","PC.40","PC.I","PC.27","PC.16","PC.97","PC.52","PC.2","PC.109","PC.8")
resistance<-c("PC.56","PC.105","PC.112","PC.L","PC.104","PC.136","PC.130","PC.64","PC.111","PC.119","PC.116","PC.121","PC.134","PC.22","PC.102","PC.5")
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]
resistance_exp<-exp[,match(resistance_new,colnames(exp))]
sensitivity_exp<-exp[,match(sensitivity_new,colnames(exp))]
geneid<-exp1[-cou,1]

data<-cbind(geneid,resistance_exp,sensitivity_exp) ####37个样本的表达谱
symbol<-as.matrix(data[,1])####基因名
glist<-symbol
exp<-t(data[,-1])  ###纯表达谱
#colnames(exp)<-paste("g",1:length(symbol),sep="")
colnames(exp)<-symbol
exp<-as.data.frame(exp)
gt<-rep(c(0,1),c(length(resistance),length(sensitivity))) #ground_truth 耐药17 敏感20

setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/1_RNAseq_DEGs/Acetylase")
Acetylase_info<-read.table("Acetylase_info_P0.05.txt",sep = '\t',header = T)


k=100
loc<-matrix(1,37,k)  ###构建一张37行100列的矩阵，用于做测试集的样本位置
rownames(loc)<-rownames(exp) ##37个样本名
set.seed(123)
for(i in 1:k){
  rr<-sample(1:length(resistance),4)  ##耐药中随机选取2个作为作为验证集
  sr<-sample((length(resistance)+1):37,4) ##敏感中随机选取2个样本作为验证集
  loc[c(rr,sr),i]<-2 #测试集
}
IMP<-matrix(0,length(symbol),2*k)  ###基因名长度，有200列
colnames(IMP)[seq(2,2*k,2)]<-"MeanDecreaseGini"
colnames(IMP)[seq(1,2*k,2)]<-"DRG"
rownames(IMP)<-symbol
ability<-matrix(NA,6,k)
rownames(ability)<-c("sen","spe","acc","val_sen","val_sep","val_acc")

####
for(j in 1:k){
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
  ability[1,j]<-dat1.rf$ confusion[1,1]/(length(resistance)-4)
  ability[2,j]<-dat1.rf$ confusion[2,2]/(length(sensitivity)-4)
  ability[3,j]<-(dat1.rf$ confusion[1,1]+dat1.rf$ confusion[2,2])/29
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
setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/5_classifier/randomForest/Acetylase_DEGs")
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
exp<-exp2[-cou,]
sensitivity<-c("PC.18","PC.81","PC.98","PC.115","PC.78","PC.135","PC.G","PC.13","PC.139","PC.117","PC.14","PC.101","PC.40","PC.I","PC.27","PC.16","PC.97","PC.52","PC.2","PC.109","PC.8")
resistance<-c("PC.56","PC.105","PC.112","PC.L","PC.104","PC.136","PC.130","PC.64","PC.111","PC.119","PC.116","PC.121","PC.134","PC.22","PC.102","PC.5")
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]
resistance_exp<-exp[,match(resistance_new,colnames(exp))]
sensitivity_exp<-exp[,match(sensitivity_new,colnames(exp))]
geneid<-exp1[-cou,1]

data<-cbind(geneid,resistance_exp,sensitivity_exp) ####37个样本的表达谱
symbol<-as.matrix(data[,1])####基因名
glist<-symbol
exp<-t(data[,-1])  ###纯表达谱
#colnames(exp)<-paste("g",1:length(symbol),sep="")
colnames(exp)<-symbol
exp<-as.data.frame(exp)
gt<-rep(c(0,1),c(length(resistance),length(sensitivity))) #ground_truth 耐药17 敏感20

setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/1_RNAseq_DEGs/Acetylase")
Acetylase_info<-read.table("Acetylase_info_P0.05.txt",sep = '\t',header = T)

#随机100次
k=100
loc<-matrix(1,37,k)
rownames(loc)<-rownames(exp)
set.seed(123)
for(i in 1:k){
  rr<-sample(1:length(resistance),4)
  sr<-sample((length(resistance)+1):37,4)
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
  tuned<-tune.svm(train.gt~.,data = dat1,gamma = 10^(-6:-1),cost = 10^(0:4))
  summary(tuned)
  model.tuned<-svm(train.gt~.,data = dat1,gamma=tuned$best.parameters$gamma,cost=tuned$best.parameters$cost)
  print(summary(model.tuned))
  pred<-predict(model.tuned,dat1)
  
  ability[1,j]<-length(grep(0,pred[1:(length(resistance)-4)]))/(length(resistance)-4)
  ability[2,j]<-length(grep(1,pred[(length(resistance)-4+1):29]))/(29-length(resistance)+4-1)
  ability[3,j]<-(length(grep(0,pred[1:(length(resistance)-4)]))+length(grep(1,pred[(length(resistance)-4+1):29])))/29
  
  pre<- predict(model.tuned,newdata=test.data)
  print(pre)
  ability[4,j]<-length(grep(0,pre[1:4]))/4
  ability[5,j]<-length(grep(1,pre[5:8]))/4
  ability[6,j]<-(length(grep(0,pre[1:4]))+length(grep(1,pre[5:8])))/8
}
IMP1<-IMP[match(Acetylase_info[,1],colnames(exp)),]
setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/5_classifier/SVM/Acetylase_DEGs")
write.table(loc,"Acetylase_DEGs_train_test4_SVM.txt",sep="\t",col.names=T)
write.table(IMP1,"Acetylase_DEGs_MeanDecreaseGini4_SVM.txt",sep="\t",col.names=T)
write.table(ability,"Acetylase_DEGs_train_test_ability4_SVM.txt",sep="\t",col.names=T)


################   分类器 randomForest  ATAC
rm(list=ls())
setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
library("data.table")
exp1<-fread("star_rsem.GeneSymbol.TPM.xls",header=T,data.table=F)
exp2<-log2(exp1[,-1]+1)
delete<-apply(exp2,1,function(x) mean(x==0))
cou<-which(delete>0.5)####在超过一半样本以上的基因删去
exp<-exp2[-cou,]
sensitivity<-c("PC.18","PC.81","PC.98","PC.115","PC.78","PC.135","PC.G","PC.13","PC.139","PC.117","PC.14","PC.101","PC.40","PC.I","PC.27","PC.16","PC.97","PC.52","PC.2","PC.109","PC.8")
resistance<-c("PC.56","PC.105","PC.112","PC.L","PC.104","PC.136","PC.130","PC.64","PC.111","PC.119","PC.116","PC.121","PC.134","PC.22","PC.102","PC.5")
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]
resistance_exp<-exp[,match(resistance_new,colnames(exp))]
sensitivity_exp<-exp[,match(sensitivity_new,colnames(exp))]
geneid<-exp1[-cou,1]

data<-cbind(geneid,resistance_exp,sensitivity_exp) ####37个样本的表达谱
symbol<-as.matrix(data[,1])####基因名
glist<-symbol
exp<-t(data[,-1])  ###纯表达谱
#colnames(exp)<-paste("g",1:length(symbol),sep="")
colnames(exp)<-symbol
exp<-as.data.frame(exp)
gt<-rep(c(0,1),c(length(resistance),length(sensitivity))) #ground_truth 耐药17 敏感20

#setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/1_RNAseq_DEGs/Acetylase")
#Acetylase_info<-read.table("Acetylase_info_P0.05.txt",sep = '\t',header = T)
Acetylase_info<-c("ATF3","SPDEF","GLI3","CREB5","BATF","FOSL2","FOXK2","SOX2","SOX10","GATA3","MAFA","HOXA9","HOXA2","BHLHE40" )

k=100
loc<-matrix(1,37,k)  ###构建一张37行100列的矩阵，用于做测试集的样本位置
rownames(loc)<-rownames(exp) ##37个样本名
set.seed(123)
for(i in 1:k){
  rr<-sample(1:length(resistance),4)  ##耐药中随机选取2个作为作为验证集
  sr<-sample((length(resistance)+1):37,4) ##敏感中随机选取2个样本作为验证集
  loc[c(rr,sr),i]<-2 #测试集
}
IMP<-matrix(0,length(symbol),2*k)  ###基因名长度，有200列
colnames(IMP)[seq(2,2*k,2)]<-"MeanDecreaseGini"
colnames(IMP)[seq(1,2*k,2)]<-"DRG"
rownames(IMP)<-symbol
ability<-matrix(NA,6,k)
rownames(ability)<-c("sen","spe","acc","val_sen","val_sep","val_acc")

####
for(j in 1:k){
  print(j)
  train.data<-exp[which(loc[,j]==1),]
  train.gt<-gt[which(loc[,j]==1)]
  test.data<-exp[which(loc[,j]==2),]
  test.gt<-gt[which(loc[,j]==2)]
  
  dat<-cbind(train.gt,train.data)
  Acetylase_DEGs<-Acetylase_info
  IMP[match(Acetylase_DEGs,colnames(exp)),2*j-1]<-1
  
  dat1<-dat[,c(1,match(Acetylase_DEGs,colnames(dat)))]
  dat1$train.gt[dat1$train.gt==0]="resistant"
  dat1$train.gt[dat1$train.gt==1]="sensitive"
  dat1$train.gt<-as.factor(dat1$train.gt)
  
  library(randomForest)
  set.seed(123)
  dat1.rf <- randomForest(train.gt ~ ., data=dat1, importance=T,proximity=T)
  print(dat1.rf) #展示随机森林模型简要信息
  ability[1,j]<-dat1.rf$ confusion[1,1]/(length(resistance)-4)
  ability[2,j]<-dat1.rf$ confusion[2,2]/(length(sensitivity)-4)
  ability[3,j]<-(dat1.rf$ confusion[1,1]+dat1.rf$ confusion[2,2])/29
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
IMP1<-IMP[match(Acetylase_info,colnames(exp)),]
setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/5_classifier/randomForest/DETFs")
write.table(loc,"DETFs_train_test4.txt",sep="\t",col.names=F)
write.table(IMP1,"DETFs_MeanDecreaseGini4.txt",sep="\t",col.names=F)
write.table(ability,"DETFs_train_test_ability4.txt",sep="\t",col.names=F)

###############    分类器  SVM  ATAC
rm(list=ls())
setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
library("data.table")
exp1<-fread("star_rsem.GeneSymbol.TPM.xls",header=T,data.table=F)
exp2<-log2(exp1[-1]+1)
delete<-apply(exp2,1,function(x) mean(x==0))
cou<-which(delete>0.5)####在超过一半样本以上的基因删去
exp<-exp2[-cou,]
sensitivity<-c("PC.18","PC.81","PC.98","PC.115","PC.78","PC.135","PC.G","PC.13","PC.139","PC.117","PC.14","PC.101","PC.40","PC.I","PC.27","PC.16","PC.97","PC.52","PC.2","PC.109","PC.8")
resistance<-c("PC.56","PC.105","PC.112","PC.L","PC.104","PC.136","PC.130","PC.64","PC.111","PC.119","PC.116","PC.121","PC.134","PC.22","PC.102","PC.5")
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]
resistance_exp<-exp[,match(resistance_new,colnames(exp))]
sensitivity_exp<-exp[,match(sensitivity_new,colnames(exp))]
geneid<-exp1[-cou,1]

data<-cbind(geneid,resistance_exp,sensitivity_exp) ####37个样本的表达谱
symbol<-as.matrix(data[,1])####基因名
glist<-symbol
exp<-t(data[,-1])  ###纯表达谱
#colnames(exp)<-paste("g",1:length(symbol),sep="")
colnames(exp)<-symbol
exp<-as.data.frame(exp)
gt<-rep(c(0,1),c(length(resistance),length(sensitivity))) #ground_truth 耐药17 敏感20

#setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/1_RNAseq_DEGs/Acetylase")
#Acetylase_info<-read.table("Acetylase_info_P0.05.txt",sep = '\t',header = T)
Acetylase_info<-c("ATF3","SPDEF","GLI3","CREB5","BATF","FOSL2","FOXK2","SOX2","SOX10","GATA3","MAFA","HOXA9","HOXA2","BHLHE40" )

#随机100次
k=100
loc<-matrix(1,37,k)
rownames(loc)<-rownames(exp)
set.seed(123)
for(i in 1:k){
  rr<-sample(1:length(resistance),4)
  sr<-sample((length(resistance)+1):37,4)
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
  Acetylase_DEGs<-Acetylase_info
  IMP[match(Acetylase_DEGs,colnames(exp)),j]<-1
  dat1<-dat[,c(1,match(Acetylase_DEGs,colnames(dat)))]
  dat1$train.gt<-as.factor(dat1$train.gt)
  
  set.seed(100) # for reproducing results
  tuned<-tune.svm(train.gt~.,data = dat1,gamma = 10^(-6:-1),cost = 10^(0:4))
  summary(tuned)
  model.tuned<-svm(train.gt~.,data = dat1,gamma=tuned$best.parameters$gamma,cost=tuned$best.parameters$cost)
  print(summary(model.tuned))
  pred<-predict(model.tuned,dat1)
  
  ability[1,j]<-length(grep(0,pred[1:(length(resistance)-4)]))/(length(resistance)-4)
  ability[2,j]<-length(grep(1,pred[(length(resistance)-4+1):29]))/(29-length(resistance)+4-1)
  ability[3,j]<-(length(grep(0,pred[1:(length(resistance)-4)]))+length(grep(1,pred[(length(resistance)-4+1):29])))/29
  
  pre<- predict(model.tuned,newdata=test.data)
  print(pre)
  ability[4,j]<-length(grep(0,pre[1:4]))/4
  ability[5,j]<-length(grep(1,pre[5:8]))/4
  ability[6,j]<-(length(grep(0,pre[1:4]))+length(grep(1,pre[5:8])))/8
}
IMP1<-IMP[match(Acetylase_info,colnames(exp)),]
setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/5_classifier/SVM/DETFs")
write.table(loc,"DETFs_train_test4_SVM.txt",sep="\t",col.names=T)
write.table(IMP1,"DETFs_MeanDecreaseGini4_SVM.txt",sep="\t",col.names=T)
write.table(ability,"DETFs_train_test_ability4_SVM.txt",sep="\t",col.names=T)



#######################################################################
#####################    TCGA  生存分析
rm(list=ls())
setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/6_TCGA")
chemotherapy_drug_sample1<-read.table("drug_paad1.txt",header=T,sep="\t") #化疗药物信息
chemotherapy_drug_sample<-chemotherapy_drug_sample1[which(chemotherapy_drug_sample1$pharmaceutical_therapy_drug_name=="Oxaliplatin"),]
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

library(org.Hs.eg.db)
library(stringr)
library(clusterProfiler)
gene=bitr(rownames(drug_sample_follow_exp),fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") 
int_gene<-intersect(gene[,1],rownames(drug_sample_follow_exp))
exp1<-drug_sample_follow_exp[match(int_gene,rownames(drug_sample_follow_exp)),]
rownames(exp1)<-gene[match(int_gene,gene[,1]),2]
g<-"SOX10"
median<-mean(as.numeric(exp1[match(g,rownames(exp1)),]))###均值
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
library(survival)
library(survminer)
summary(coxph(Surv(TTP, status_TTP) ~ label,data=surv_info1)) 
surv_TTP<-survfit(Surv(TTP, status_TTP) ~ label,data=surv_info1)
surv_TTP
ggsurvplot(surv_TTP,
           pval = TRUE, conf.int = TRUE, #是否显示p值/显示风险区域
           risk.table = TRUE # 将风险表显示在生存曲线下面
)
##########################################################################
#################     生存相关基因
rm(list=ls())
setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/6_TCGA")
#chemotherapy_drug_sample1<-read.table("drug_paad1.txt",header=T,sep="\t") #化疗药物信息
#chemotherapy_drug_sample<-chemotherapy_drug_sample1[which(chemotherapy_drug_sample1$pharmaceutical_therapy_drug_name=="Oxaliplatin"),]
chemotherapy_drug_sample<-read.table("drug_paad1.txt",header=T,sep="\t") #化疗药物信息

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
colnames(result1)<-c("geneid","regression coefficient","pvalue")
surv_coff_gene1<-result1[which(result1[,3]<0.05),]

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

input  <-list(unique(surv_coff_gene[,1]),unique(DEGs),unique(DApeaks_gene))
venn(input,showSetLogicLabel=TRUE)
tmp <- venn(input)
int<-attr(tmp, "intersections")

setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/6_TCGA")
library(VennDiagram)
venn.diagram(x=list(Survival_related_genes=unique(surv_coff_gene[,1]),RNA_DEGs=unique(DEGs),ATAC_DEGs=unique(DApeaks_gene)), "surv_coff_DEGs.png",fill=c("red","green","blue"))

##对交叠的基因使用survfit()函数进行KM生存分析
intt<-int[["A:B:C"]]
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
