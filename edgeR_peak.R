###################      药物敏感和耐药
rm(list=ls())
setwd("~/xjj/drug")
ic50<-read.csv("changhai_ic50.csv",header = TRUE,sep = ",",row.names=1)
library(edgeR)


##中位数为区分点
want_row_na<-ic50[which(rownames(ic50)=="S1648"),]
colname<-colnames(ic50)[which(!is.na(want_row_na))]
value1<-want_row_na[which(!is.na(want_row_na))]
value2<-as.numeric(value1)
value<-log2(value2+1)
#value<-want_row_na[which(!is.na(want_row_na))]
resistance<-colname[which(value>median(as.numeric(value)))]   
sensitivity<-colname[which(value<=median(as.numeric(value)))] 

setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]

resistance_peak<-peak_count_name[,match(resistance_new,colnames(peak_count_name))]
sensitivity_peak<-peak_count_name[,match(sensitivity_new,colnames(peak_count_name))]
peak_name<-peak_count_name[,1]

library(DESeq2)
colDate<-data.frame(row.names = c(as.vector(sensitivity_new),as.vector(resistance_new)),
                    condition=factor(c(rep("sensitivity",length(sensitivity_new)),rep("resistance",length(resistance_new))))
)

datexpr<-cbind(sensitivity_peak,resistance_peak)
counts <- apply(datexpr,2,as.numeric)   ###矩阵中必须是数值
dds<-DESeqDataSetFromMatrix(countData = counts,colData = colDate,design = ~condition)
dds##这是一个关于各种内容的矩阵
dds<-DESeq(dds)##进行标准化分析
sizeFactors(dds)##查看每个主成分的标准化值
res<-results(dds)##将结果输出
res
class(res)##可以看出它是DESeq的属性，要转化为表格
res<-as.data.frame(res)
head(res)
res<-cbind(peak_count_name[,1],res)  ##对数据增加一列
head(res)
colnames(res)<- c('peak_id',"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj" )##对列名进行修改
#setwd("~/xjj/drug/drug_result")
setwd("~/xjj/drug/drug_result/log_result")
write.table(res,"resistance-sensitivity-all-DESeq2_S1648_log2median.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

setwd("~/xjj/drug/drug_result")
res<-read.table("resistance-sensitivity-all-DESeq2_S1648_log2mean.txt",header = TRUE,sep = "\t")
#对差异基因的结果进行差异筛选，本例采用的是p值小于0.05,log2foldchange绝对值大于一
resSig<-res[which(res$padj<0.05 & abs(res$log2FoldChange>1)),]
# pvalue padj
##新增一列，将log2FoldChange>0标注为up，<0标准为down
resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
head(resSig)
##保存数据
write.table(resSig,"resistance-vs-sensitivity-fdr-0.05-FC-1.median_IRI.txt",sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')


#####  peak.txt
b<-data.frame(resSig[,1])
a<-apply(b,1,function(x) unlist(strsplit(as.character(x), "_")))
d<-t(a)
e<-data.frame(rep("+",nrow(b)))
peaks<-cbind(b,d,e)
colnames(peaks)<-c("peak_id","chr","start","end","strand")
#setwd("~/xjj/drug/drug_result")
setwd("~/xjj/drug/drug_result/no_log_result")
write.table(peaks,"S1648-mean-fdr-0.05-FC-1.txt",sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')
write.table(d,"S1648-mean-fdr-0.05-FC-1.gff",sep = '\t',
            col.names = F,row.names = F,quote = FALSE)

#####  bed
b<-data.frame(resSig[,1])
a<-apply(b,1,function(x) unlist(strsplit(as.character(x), "_")))
d<-t(a)
e<-data.frame(rep("+",nrow(b)))
f<-data.frame(rep(1,nrow(b)))
peaks_bed<-cbind(d,b,f,e)
colnames(peaks_bed)<-c("chr","start","end","peak_id","score","strand")
setwd("~/xjj/drug/drug_result/log_result")
write.table(peaks_bed,"S1648-mean-fdr-0.05-FC-1.bed",sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')



