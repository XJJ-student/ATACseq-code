rm(list=ls())
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge")
first_category_name = list.files("txt")  

peak_count<-NULL
for(i in 1:length(first_category_name)){
  setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt")
  peak<-read.table(first_category_name[i],sep="\t",header=F)
  count<-peak[,4]
  peak_count<-cbind(peak_count,count)
}
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt")
peak<-read.table(first_category_name[1],sep="\t",header=F)
name<-paste(peak[,1],peak[,2],peak[,3],sep='_')
peak_count_name<-cbind(name,peak_count)

library(dplyr)
first_category_name1<-lapply(strsplit(first_category_name,'_'), function(x) x[1])%>%unlist()
colnames(peak_count_name)<-c("peak_name",first_category_name1)

setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_count_name<-read.table("peak_count.txt",sep="\t",header=T)
###################      药物敏感和抑制
#setwd("~/xjj/drug")
#ic50<-read.csv("changhai_ic50.csv",header = TRUE,sep = ",",row.names=1)
#resistance<-as.character(ic50[ic50[,2]==1,1])
#sensitivity<-as.character(ic50[ic50[,2]==2,1]) #AUC越大对药物越敏感

sensitivity<-c("PC.64","PC.81","PC.116","PC.L","PC.101","PC.27","PC.78","PC.98","PC.130","PC.117","PC.8","PC.109","PC.136","PC.16","PC.97","PC.2","PC.115","PC.13")
resistance<-c("PC.G","PC.111","PC.134","PC.14","PC.40","PC.135","PC.52","PC.22","PC.139","PC.18","PC.105","PC.104","PC.56","PC.119","PC.112","PC.102","PC.121","PC.5","PC.I")

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
setwd("~/xjj/drug/drug_result/HDAC_drug/HDAC13_heatmap")
write.table(res,"-sensitivity-resistance-all-DESeq2_heatmap13HDAC_ATAC.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出



setwd("~/xjj/drug/drug_result/HDAC_drug/heatmap")
res<-read.table("resistance-sensitivity-all-DESeq2_heatmap.txt",header = TRUE,sep = "\t")
#对差异基因的结果进行差异筛选，本例采用的是p值小于0.05,log2foldchange绝对值大于一
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

##############   单个HDAC药物数据进行实验
rm(list=ls())
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge")
first_category_name = list.files("txt")  
peak_count<-NULL
for(i in 1:length(first_category_name)){
  setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt")
  peak<-read.table(first_category_name[i],sep="\t",header=F)
  count<-peak[,4]
  peak_count<-cbind(peak_count,count)
}
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt")
peak<-read.table(first_category_name[1],sep="\t",header=F)
name<-paste(peak[,1],peak[,2],peak[,3],sep='_')
peak_count_name<-cbind(name,peak_count)

library(dplyr)
first_category_name1<-lapply(strsplit(first_category_name,'_'), function(x) x[1])%>%unlist()
colnames(peak_count_name)<-c("peak_name",first_category_name1)
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
write.table(peak_count_name,"peak_count.txt",sep="\t",row.names=TRUE,col.names=TRUE)
peak_count_name<-read.table("peak_count.txt",sep="\t",header=T)
###################      药物敏感和抑制
rm(list=ls())
setwd("~/xjj/drug")
ic50<-read.csv("changhai_ic50.csv",header = TRUE,sep = ",",row.names=1)

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

# 取前后各30%
want_row_na<-ic50[which(rownames(ic50)=="S7555"),]
colname<-colnames(ic50)[which(!is.na(want_row_na))]
value<-want_row_na[which(!is.na(want_row_na))]
sensitivity<-colname[order(value)[1:ceiling(length(value)*0.3)]]
resistance<-colname[order(value)[ceiling(length(value)*0.7):length(value)]]

setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]

resistance_peak<-peak_count_name[,match(resistance_new,colnames(peak_count_name))]
sensitivity_peak<-peak_count_name[,match(sensitivity_new,colnames(peak_count_name))]
peak_name<-peak_count_name[,1]
resistance_peak1 <- apply(resistance_peak,2,as.numeric)
sensitivity_peak1 <- apply(sensitivity_peak,2,as.numeric)

resistance_peak11<-apply(resistance_peak1,1,mean)
sensitivity_peak11<-apply(sensitivity_peak1,1,mean)

library(DESeq2)
colDate<-data.frame(row.names = c(as.vector(resistance_new),as.vector(sensitivity_new)),
                    condition=factor(c(rep("resistance",length(resistance_new)),rep("sensitivity",length(sensitivity_new))))
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
#setwd("~/xjj/drug/drug_result")
setwd("~/xjj/drug/drug_result/HDAC_drug/S1030")
write.table(res,"sensitivity-resistance-all-DESeq2_S1030_30.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

setwd("~/xjj/drug/drug_result/HDAC_drug/S2759")
res<-read.table("sensitivity-resistance-all-DESeq2_S2759_30.txt",header = TRUE,sep = "\t")
#对差异基因的结果进行差异筛选，本例采用的是p值小于0.05,log2foldchange绝对值大于一
#resSig<-res[which(res$padj<0.05 & abs(res$log2FoldChange>1)),]
resSig<-res[which(res$pvalue<0.05),]
# pvalue padj

##新增一列，将log2FoldChange>0标注为up，<0标准为down
resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
sum(resSig$up_down=='up')
sum(resSig$up_down=='down')
length(which(res$log2FoldChange<0))
length(which(res$log2FoldChange<(-1)))
##保存数据
#write.table(resSig,"resistance-vs-sensitivity-fdr-0.05-FC-0-S1648.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')

up_gene<-as.data.frame(resSig[which(resSig$log2FoldChange>0),1])
down_gene<-as.data.frame(resSig[which(resSig$log2FoldChange<0),1])

up1<-match(up_gene[,1],peak_count_name[,1])
down1<-match(down_gene[,1],peak_count_name[,1])
sum((resistance_peak11[up1]-sensitivity_peak11[up1])>0)
sum((resistance_peak11[down1]-sensitivity_peak11[down1])>0)

library(dplyr)
a_up<-apply(up_gene,1,function(x) unlist(strsplit(as.character(x), "_")))
d_up<-t(a_up)
e_up<-data.frame(rep("+",nrow(d_up)))
peaks_up<-cbind(up_gene,d_up,e_up)
colnames(peaks_up)<-c("peak_id","chr","start","end","strand")
setwd("~/xjj/drug/drug_result/HDAC_drug/S2759")
out_file<-"S2759-30-p-0.05-logFC0"
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


########  DE peaks homer gene
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDAC_drug/HDAC13_heatmap")
library(data.table)
homer_anno<-fread("sen-res-HDAC13-heatmap-0.05-P-down-annotation.txt",header=T,data.table=F)
dim(homer_anno)
aiya<-apply(as.data.frame(homer_anno$`Gene Name`),1,function(x) substring(x, 1, 4))
df56<-which(aiya=="HDAC")
homer_anno$`Gene Name`[df56]
a<-homer_anno[order(homer_anno$`Gene Name`),]

###################  HDAC DEGs ∩ Diff peaks anno genes 
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDAC_drug/HDAC13_heatmap")
library(data.table)
homer_anno1<-fread("sen-res-HDAC13-heatmap-0.05-P-down-annotation.txt",header=T,data.table=F)
homer_anno2<-fread("sen-res-HDAC13-heatmap-0.05-P-up-annotation.txt",header=T,data.table=F)
colnames(homer_anno1)<-c("peakid",colnames(homer_anno1)[-1])
colnames(homer_anno2)<-c("peakid",colnames(homer_anno2)[-1])
homer_anno<-rbind(homer_anno1,homer_anno2)
dim(homer_anno)
setwd("~/xjj/drug/drug_result/HDAC_drug/HDAC13_heatmap/DEGs")
DEGs<-fread("sensitivity-resistance-heatmap13HDAC_P0.01_down_DEGs.txt",header=T,data.table=F)
dim(DEGs)
a<-intersect(DEGs[,1],homer_anno[,16])
DEGs_gene<-DEGs[match(a,DEGs[,1]),c(1,3)]
length(a)
setwd("~/xjj/drug/drug_result/HDAC_drug/HDAC13_heatmap/DEGs_intersect_Diffpeakanno_gene")
write.table(DEGs_gene,"heatmap13HDAC_P0.01_down_anno_DEGs.txt",sep = '\t',col.names = F,row.names = F,quote = FALSE)##数据输出


###############  GSEA  cancer module
rm(list=ls())
library(clusterProfiler)
library(GSEABase)
library(org.Hs.eg.db)
library(stringr)
####  DEGs
setwd("~/xjj/drug/drug_result/HDAC_drug/HDAC13_heatmap/DEGs_intersect_Diffpeakanno_gene")
DEGs<-read.delim("heatmap13HDAC_P0.05_up_anno_DEGs.txt",header=F,sep='\t')
colnames(DEGs)<-c("symbol","logFC")
geneList<-DEGs$logFC #第二列可以是folodchange，也可以是logFC
names(geneList)=DEGs$symbol #使用转换好的ID
geneList=sort(geneList,decreasing =T) #从高到低排序
gene1<-names(geneList)
setwd("~/xjj/GSEA")
c5<-read.delim("h.all.v7.4.symbols.txt",header=F,sep='\t')
result<-NULL
for(i in 1:nrow(c5)){
  a<-c5[i,which(!is.na(c5[i,-1]))]
  a1<-a[-1]
  a2<-t(a1)
  TF_name<-t(data.frame(rep(a[1],length(a1))))
  result1<-cbind(TF_name,a2)
  result<-rbind(result,result1)
}
egmt3 <- GSEA(geneList, TERM2GENE=result, verbose=FALSE,pvalueCutoff = 1)##不卡p值
head(egmt3)
egmt4=data.frame(egmt3)
position<-which(egmt4$pvalue<0.05)
#position<-which(egmt4$pvalue<0.05 & abs(egmt4$NES)>1 & egmt4$qvalues<0.25)
egmt5<-egmt4[position,]
library(data.table)
library(dplyr)
d<-data.frame(egmt5[,11])
egmt6<-apply(d,1,function(x) unlist(strsplit(as.character(x), "/")))%>%unlist()
int_DEGs<-intersect(gene1,egmt6)
result11<-data.frame()
for(i in 1:length(unique(int_DEGs))){
  result11[i,1]<-int_DEGs[i]
  result11[i,2]<-sum(egmt6%in%int_DEGs[i])
}
result12=result11[order(result11[,2],decreasing =T),]###筛选出在cancer module中最富集的基因
colnames(result12)<-c("Gene","Number")
setwd("~/xjj/drug/drug_result/HDAC_drug/HDAC13_heatmap/GSEA_result")
write.table(result12,"onlyP0.05_heatmap13HDAC_P0.05_up_anno_DEGs_GSEA_hallmark_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(egmt5,"onlyP0.05_heatmap13HDAC_P0.05_up_anno_DEGs_GSEA_hallmark_gene_sets.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
#write.table(result12,"heatmap13HDAC_P0.05_down_anno_DEGs_GSEA_cancer_modules.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

################ clusterProfiler
rm(list=ls())
library(org.Hs.eg.db)
library(clusterProfiler)
setwd("~/xjj/drug/drug_result/HDAC_drug/HDAC13_heatmap/cancer_associated_gene_drug_target")
up_gene<-read.table("heatmap13HDAC_P0.05_up_anno_DEGs_GSEA_cancer_modules_drug_target_summary.txt",header=T,sep="\t")

gene=bitr(up_gene[,3],fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
## 去重
library(stringr)
kegg<-enrichKEGG(gene[,2],organism = 'hsa',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'BH',
                 minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.05,use_internal_data = FALSE)
ee<-kegg@result
ee1<-ee[which(ee$pvalue<0.05),]

HDACi<-c("HDAC1","HDAC2","HDAC3","HDAC4","HDAC5","HDAC6","HDAC7","HDAC8","HDAC9","HDAC10","HDAC11","HDAC12","HDAC13","HDAC14","HDAC15")
HDAC_gene=bitr(HDACi,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
kegg1<-enrichKEGG(HDAC_gene[,2],organism = 'hsa',keyType = 'kegg',pvalueCutoff = 1,pAdjustMethod = 'BH',
                  minGSSize = 10,maxGSSize = 500,qvalueCutoff = 1,use_internal_data = FALSE)
kegg_BP_HDAC_gene<-kegg1@result
int_pathway<-intersect(ee1[,2],kegg_BP_HDAC_gene[,2])
HDAC_DEG_kegg_pathway<-ee1[match(int_pathway,ee1[,2]),]
setwd("G:\\类器官\\mRNA_star_rsem\\DESeq2_conut")
write.table(HDAC_DEG_kegg_pathway,"kegg_DEG_pathway_P0.01_down_DEGs_intersect_Diffpeakanno_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

enrich_gene<-HDAC_DEG_kegg_pathway$geneID
pathway_gene<-unlist(strsplit(enrich_gene,split="/"))
pathway_gene2=bitr(pathway_gene,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
write.table(pathway_gene,"kegg_pathway_P0.01_down_gene.txt",sep = '\t',col.names = F,row.names = F,quote = FALSE)##数据输出

#dotplot(kegg,showCategory=30)
go<-enrichGO(gene[,2],OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.05,keyType = 'ENTREZID')
a<-go@result
go_BP1<-a[a[,1]=="BP",]
go_BP<-go_BP1[which(go_BP1$pvalue<0.05),]
#dotplot(go_BP,showCategory=30)

HDACi<-c("HDAC1","HDAC2","HDAC3","HDAC4","HDAC5","HDAC6","HDAC7","HDAC8","HDAC9","HDAC10","HDAC11","HDAC12","HDAC13","HDAC14","HDAC15")
HDAC_gene=bitr(HDACi,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
go1<-enrichGO(HDAC_gene[,2],OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff =1,keyType = 'ENTREZID')
a1<-go1@result
go_BP_HDAC_gene<-a1[a1[,1]=="BP",]
int_pathway<-intersect(go_BP[,3],go_BP_HDAC_gene[,3])
HDAC_DEG_pathway<-go_BP[match(int_pathway,go_BP[,3]),]
setwd("G:\\类器官\\mRNA_star_rsem\\DESeq2_conut")
write.table(HDAC_DEG_pathway,"GO_DEG_pathway_P0.01_down_DEGs_intersect_Diffpeakanno_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
enrich_genego<-HDAC_DEG_pathway$geneID
pathway_genego<-unlist(strsplit(enrich_genego,split="/"))
write.table(pathway_genego,"GO_pathway_P0.01_down_gene.txt",sep = '\t',col.names = F,row.names = F,quote = FALSE)##数据输出


#######    homer-TFs-HDACi
### homer-TFs
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDAC_drug/heatmap/homer_result/sen-res-HDAC-heatmap-0.001-p-logFC0-DESeq2-up")
knownResults<-read.delim("knownResults.txt",header=T,sep='\t')
FDR05<-knownResults[knownResults$q.value..Benjamini.<0.05,]
library(dplyr)
homer_TF1<-data.frame(lapply(strsplit(as.character(FDR05[,1]),'/'), function(x) x[3])%>%unlist())
FDR051<-FDR05[which(homer_TF1=="Homer"),]
homer_TF2<-data.frame(lapply(strsplit(as.character(FDR051[,1]),'/'), function(x) x[1])%>%unlist())
library(stringr)
homer_TF3<-apply(homer_TF2,1,function(x) str_replace_all(as.character(x),"\\((.*?)\\)",""))##取出括号前的字符
homer_TF<-toupper(homer_TF3)
length(homer_TF)
homer_TF3<-as.matrix(homer_TF3)
colnames(homer_TF3)<-"TFs"
write.table(homer_TF3,"P0.001-up-TFs.txt",sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')

### anno
setwd("~/xjj/drug/drug_result/HDAC_drug/S2759")
library(data.table)
homer_anno<-fread("S2759-30-p-0.05-down-annotation.txt",header=T,data.table=F)
dim(homer_anno)
aiya<-apply(as.data.frame(homer_anno$`Gene Name`),1,function(x) substring(x, 1, 4))
df56<-which(aiya=="HDAC")
a<-unique(homer_anno$`Gene Name`[df56])
length(homer_TF)
length(a)
a[order(a)]
aa<-a[order(a)]
### TF-target
setwd("~/xjj/drug/drug_result/HDAC_drug")
library(data.table)
TFs_targets<-fread("gene_attribute_edges.txt",header=T,data.table=F)
HDAC_TFs<-TFs_targets[TFs_targets$source %in% aa,c(1,4)]
HDAC_TFs<-TFs_targets[TFs_targets$source %in% "HDAC4",c(1,4)]
HDAC_TFs<-TFs_targets[TFs_targets$source %in% "HDAC11",c(1,4)]
HDAC_TFs<-TFs_targets[TFs_targets$source %in% "HDAC2",c(1,4)]


nrow(HDAC_TFs)
int_TFs<-intersect(homer_TF,toupper(HDAC_TFs[,2]))
length(int_TFs)
want_HDAC_TFs<-HDAC_TFs[match(int_TFs,toupper(HDAC_TFs[,2])),]
nrow(want_HDAC_TFs)

library(dplyr)
unique_want_HDAC_TFs<-want_HDAC_TFs %>% distinct(source,target, .keep_all = TRUE)
nrow(unique_want_HDAC_TFs)
unique_want_HDAC_TFs


###############   做聚类
rm(list=ls())
setwd("~/xjj/drug")
ic50<-read.csv("changhai_ic50.csv",header = TRUE,sep = ",",row.names=1)
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
heatmap=pheatmap(ic2,scale = "none",main = "log2 IC50",show_rownames=T,show_colnames=T,display_numbers = dn)






############   分类器的构建
rm(list=ls())
setwd("~/xjj/drug/drug_result") 
first_category_name1 = list.files("HDAC_drug") 
first_category_name = first_category_name1[-c(1:3)]
dir = paste("~/xjj/drug/drug_result/HDAC_drug/",first_category_name,sep="")  
n = length(dir) 
library("data.table")
result<-NULL
for(i in 1:2){      
  setwd(dir[i])
  DApeaksup<-fread(paste(first_category_name[i],"30-p-0.01-logFC0-DESeq2-up.txt",sep='-'),header=T,data.table=F)
  a<-as.data.frame(DApeaksup[,1])
  DApeaksdown<-fread(paste(first_category_name[i],"30-p-0.01-logFC0-DESeq2-down.txt",sep='-'),header=T,data.table=F)
  b<-as.data.frame(DApeaksdown[,1])
  d<-rbind(a,b)
  d1<-matrix(d,ncol=1)
  #e<-data.frame(d)
  result<-cbind(result,d1[,1])
}
resulti<-intersect(result[,1],result[,2])
resulti<-intersect(resulti,result[,20])

for(i in 3:n){
  resulti<-intersect(resulti,result[,i])
}



i=1
setwd(dir[i])
DApeaksup<-fread(paste(first_category_name[i],"30-p-0.05-logFC0-DESeq2-up.txt",sep='-'),header=T,data.table=F)
a<-DApeaksup[,1]
DApeaksdown<-fread(paste(first_category_name[i],"30-p-0.05-logFC0-DESeq2-down.txt",sep='-'),header=T,data.table=F)
b<-DApeaksdown[,1]
d<-c(a,b)
i=2
setwd(dir[i])
DApeaksup<-fread(paste(first_category_name[i],"30-p-0.05-logFC0-DESeq2-up.txt",sep='-'),header=T,data.table=F)
a<-DApeaksup[,1]
DApeaksdown<-fread(paste(first_category_name[i],"30-p-0.05-logFC0-DESeq2-down.txt",sep='-'),header=T,data.table=F)
b<-DApeaksdown[,1]
d2<-c(a,b)
resulti<-intersect(d,d2)
for(i in 3:n){      
  setwd(dir[i])
  DApeaksup<-fread(paste(first_category_name[i],"30-p-0.05-logFC0-DESeq2-up.txt",sep='-'),header=T,data.table=F)
  a<-DApeaksup[,1]
  DApeaksdown<-fread(paste(first_category_name[i],"30-p-0.05-logFC0-DESeq2-down.txt",sep='-'),header=T,data.table=F)
  b<-DApeaksdown[,1]
  d3<-c(a,b)
  resulti<-intersect(resulti,d3)
}
setwd("~/xjj/drug/drug_result")
write.table(resulti,"linshi.bed",col.names=F)
#################
i3<-c("S1848","S7555","S8495")
first_category_name<-first_category_name[-which(first_category_name=="S1848")]
first_category_name<-first_category_name[-which(first_category_name=="S7555")]
first_category_name<-first_category_name[-which(first_category_name=="S8495")]
n = length(first_category_name) 
i=1
setwd(dir[i])
DApeaksup<-fread(paste(first_category_name[i],"30-p-0.05-logFC0-DESeq2-down.txt",sep='-'),header=T,data.table=F)
d<-DApeaksup[,1]
i=2
setwd(dir[i])
DApeaksup<-fread(paste(first_category_name[i],"30-p-0.05-logFC0-DESeq2-down.txt",sep='-'),header=T,data.table=F)
d2<-DApeaksup[,1]
resulti<-intersect(d,d2)
for(i in 3:n){      
  setwd(dir[i])
  DApeaksup<-fread(paste(first_category_name[i],"30-p-0.05-logFC0-DESeq2-down.txt",sep='-'),header=T,data.table=F)
  d3<-DApeaksup[,1]
  resulti<-intersect(resulti,d3)
}


