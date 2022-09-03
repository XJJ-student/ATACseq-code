rm(list=ls())
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge")
first_category_name = list.files("txt")  

peak_count<-NULL
for(i in 1:length(first_category_name)){
  setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt")
  peak<-read.table(first_category_name[i],sep="\t",header=F)
  count<-peak[,8] #peak[,8] RPKM
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
write.table(peak_count_name,"peak_RPKM.txt",sep="\t",row.names=F,col.names=TRUE)
peak_count_name<-read.table("peak_count.txt",sep="\t",header=T)
###################      药物敏感和抑制
rm(list=ls())
setwd("~/xjj/drug")
#ic50<-read.csv("changhai_ic50.csv",header = TRUE,sep = ",",row.names=1)
ic50<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)

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
colnames(range_variance)<-c("variance","Drug")
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/6_application")
write.table(range_variance,"drug_variance.txt",sep="\t",row.names=F,col.names=TRUE)


library(ggplot2)
data$drug=factor(data$drug,levels = range_variance[,2])
ggplot(data[1:2496,])+ 
  geom_boxplot(aes(x=drug, y=as.numeric(value)))


drug_info<-read.csv("compiled_drug_info.csv",header = TRUE,sep = ",")
drug_result2<-drug_info[match(range_variance[,2],drug_info$drug_id),]
drug_result<-cbind(range_variance[,1],drug_result2)

setwd("~/xjj/drug/drug_result")
write.table(drug_result,"drug_variance.txt",sep="\t",row.names=F,col.names=TRUE)
drug_result2<-read.csv("drug_variance.txt",header = TRUE,sep = "\t")


###  得到5-FU 和 OXA 是差异最大的
#setwd("~/xjj/drug/th")
#ic50_th<-read.csv("changhai_ic50_otsu.csv",header = TRUE,sep = ",",row.names=1)
#want_row_th<-ic50_th[which(rownames(ic50_th)=="S7575"),2]

# 取前后各30%
want_row_na<-ic50[which(rownames(ic50)=="5-FU"),]
colname<-colnames(ic50)[which(!is.na(want_row_na))]
value<-want_row_na[which(!is.na(want_row_na))]
sensitivity<-colname[order(value)[1:ceiling(length(value)*0.3)]]
resistance<-colname[order(value)[floor(length(value)*0.7):length(value)]]

##中位数为区分点
#want_row_na<-ic50[which(rownames(ic50)=="S1648"),]
#colname<-colnames(ic50)[which(!is.na(want_row_na))]
#value1<-want_row_na[which(!is.na(want_row_na))]
#value2<-as.numeric(value1)
#value<-log2(value2+1)
#value<-want_row_na[which(!is.na(want_row_na))]
#resistance<-colname[which(value>median(as.numeric(value)))]   
#sensitivity<-colname[which(value<=median(as.numeric(value)))] 

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
setwd("~/xjj/drug/drug_result/other_drug")
write.table(res,"resistance-sensitivity-all-DESeq2_S7958_30.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出


rm(list=ls())
setwd("~/xjj/drug/drug_result/chemotherapy")
res<-read.table("resistance-sensitivity-all-DESeq2_PAC_30.txt",header = TRUE,sep = "\t")
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
setwd("~/xjj/drug/drug_result/other_drug/Diff_peaks_P")
out_file<-"PAC"
#write.table(peaks_up,file=paste(out_file,"-30-fdr-0.05-logFC0-DESeq2-up.txt",sep=""),sep = '\t',
#            col.names = T,row.names = F,quote = FALSE,na='')
#write.table(d_up,file=paste(out_file,"-30-fdr-0.05-logFC0-DESeq2-up.gff",sep=""),sep = '\t',
#            col.names = F,row.names = F,quote = FALSE)

a_down<-apply(down_gene,1,function(x) unlist(strsplit(as.character(x), "_")))
d_down<-t(a_down)
e_down<-data.frame(rep("+",nrow(d_down)))
peaks_down<-cbind(down_gene,d_down,e_down)
colnames(peaks_down)<-c("peak_id","chr","start","end","strand")
write.table(peaks_down,file=paste(out_file,"-30-fdr-0.05-logFC0-DESeq2-down.txt",sep=""),sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')
write.table(d_down,file=paste(out_file,"-30-fdr-0.05-logFC0-DESeq2-down.gff",sep=""),sep = '\t',
            col.names = F,row.names = F,quote = FALSE)
#####  bed
f_up<-data.frame(rep(1,nrow(d_up)))
peaks_bed_up<-cbind(d_up,up_gene,f_up,e_up)
colnames(peaks_bed_up)<-c("chr","start","end","peak_id","score","strand")
write.table(peaks_bed_up,file=paste(out_file,"-30-p-0.05-logFC0-DESeq2-up.bed",sep=""),sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')

f_down<-data.frame(rep(1,nrow(d_down)))
peaks_bed_down<-cbind(d_down,down_gene,f_down,e_down)
colnames(peaks_bed_down)<-c("chr","start","end","peak_id","score","strand")
write.table(peaks_bed_down,file=paste(out_file,"-30-p-0.05-logFC0-DESeq2-down.bed",sep=""),sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')

#####  peak.txt
b<-data.frame(resSig[,1])
a<-apply(b,1,function(x) unlist(strsplit(as.character(x), "_")))
d<-t(a)
e<-data.frame(rep("+",nrow(b)))
peaks<-cbind(b,d,e)
colnames(peaks)<-c("peak_id","chr","start","end","strand")
#setwd("~/xjj/drug/drug_result")
setwd("~/xjj/drug/drug_result/chemotherapy")
write.table(peaks,"PAC-30-fdr-0.05-FC-1.txt",sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')
write.table(d,"PAC-30-fdr-0.05-FC-1.gff",sep = '\t',
            col.names = F,row.names = F,quote = FALSE)

#####  bed
b<-data.frame(resSig[,1])
a<-apply(b,1,function(x) unlist(strsplit(as.character(x), "_")))
d<-t(a)
e<-data.frame(rep("+",nrow(b)))
f<-data.frame(rep(1,nrow(b)))
peaks_bed<-cbind(d,b,f,e)
colnames(peaks_bed)<-c("chr","start","end","peak_id","score","strand")
setwd("~/xjj/drug/drug_result/chemotherapy")
write.table(peaks_bed,"PAC-30-fdr-0.05-FC-1.bed",sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')






setwd("~/xjj/try")
library("data.table")
hg19<-fread("hg19.refseq",header=T,data.table=F)

############  DEGs
# 取前后各30%
rm(list=ls())
setwd("~/xjj/drug")
ic50<-read.csv("changhai_ic50.csv",header = TRUE,sep = ",",row.names=1)
setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
library("data.table")
exp1<-fread("star_rsem.GeneSymbol.readscount.xls",header=T,data.table=F)
#FPKM<-fread("star_rsem.GeneSymbol.FPKM.xls",header=T,data.table=F)
exp2<-floor(exp1[-1])
delete<-apply(exp2,1,function(x) mean(x==0))
cou<-which(delete>0.5)
exp<-exp2[(-cou),]
  
want_row_na<-ic50[which(rownames(ic50)=="PAC"),]
colname<-colnames(ic50)[which(!is.na(want_row_na))]
value<-want_row_na[which(!is.na(want_row_na))]
sensitivity<-colname[order(as.numeric(value))[1:ceiling(length(value)*0.3)]]
resistance<-colname[order(as.numeric(value))[floor(length(value)*0.7):length(value)]]

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

setwd("~/xjj/DESeq_DEGs")
write.table(res,"sensitivity-resistance-all-DESeq2_PAC_30_DEGs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
#res<-read.table("sensitivity-resistance-all-DESeq2_OXA_30_DEGs.txt",header = TRUE,sep = "\t")
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
sum((resistance_exp11[up1]-sensitivity_exp11[up1])>0)
sum((resistance_exp11[down1]-sensitivity_exp11[down1])<0)

setwd("~/xjj/DESeq_DEGs/DESeq2_P05")
write.table(up_gene,"sensitivity-resistance-all-DESeq2_PAC_P_up_DEGs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(down_gene,"sensitivity-resistance-all-DESeq2_PAC_P_down_DEGs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
##################
NM<-hg19[match(res[,1],hg19[,6]),1]
BSF1<-cbind(NM,res[,c(3,7)])
BSF<-BSF1[which(!is.na(NM)),]
colnames(BSF)<-c("# ID","Status","Value")
write.table(BSF,"resistance-sensitivity-OXA_30_DEGs_BSF.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

############


resSig<-res[which(res$padj<0.05),]
# pvalue padj
resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
sum(resSig$up_down=='up')
sum(resSig$up_down=='down')
##保存数据
write.table(resSig,"resistance-vs-sensitivity-fdr-0.05-FC-0-S1648.txt",sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')




