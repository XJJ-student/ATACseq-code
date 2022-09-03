############### 在GEM耐药的的样本中进一步分HDAC敏感和耐药样本 #############
#1 #####设置GEM的AUC<0.5 and S8495的阈值设置为0.5
rm(list=ls())
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
setwd("~/xjj/drug")
ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,!(colnames(ic501)%in%c("PC.100","PC.34"))]## all of PDAC
colnames(ic50)<-sample_id[match(colnames(ic50),sample_id[,2]),1]
want_AUC_chemo<-ic50[rownames(ic50)%in%"GEM",]
want_AUC_HDAC<-ic50[rownames(ic50)%in%"S8495",]

############  中位数标准
GEM_sensitivity<-colnames(want_AUC_chemo)[which(as.numeric(want_AUC_chemo)>median(as.numeric(want_AUC_chemo)))]
GEM_resistance<-colnames(want_AUC_chemo)[which(as.numeric(want_AUC_chemo)<=median(as.numeric(want_AUC_chemo)))]
S8495_sensitivity<-colnames(want_AUC_HDAC)[which(as.numeric(want_AUC_HDAC)>median(as.numeric(want_AUC_HDAC)))]
S8495_resistance<-colnames(want_AUC_HDAC)[which(as.numeric(want_AUC_HDAC)<=median(as.numeric(want_AUC_HDAC)))]

GEM_res_S8495_sen<-intersect(GEM_resistance,S8495_sensitivity)
GEM_res_S8495_res<-intersect(GEM_resistance,S8495_resistance)
GEM_res_S8495_sen<-c("DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-12","DAC-14")
GEM_res_S8495_res<-c("DAC-20","DAC-21","DAC-22","DAC-23","DAC-25","DAC-29","DAC-30","DAC-32","DAC-35","DAC-5","DAC-13","DAC-39")


############################################################
#2 ###########    RNAseq DEGs
rm(list=ls())
setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
library("data.table")
exp11<-fread("star_rsem.GeneSymbol.readscount.xls",header=T,data.table=F)
GEM_res_S8495_sen<-c("DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-12","DAC-14")
GEM_res_S8495_res<-c("DAC-20","DAC-21","DAC-22","DAC-23","DAC-25","DAC-29","DAC-30","DAC-32","DAC-35","DAC-5","DAC-13","DAC-39")
all<-c(GEM_res_S8495_sen,GEM_res_S8495_res)
exp1<-exp11[,match(all,colnames(exp11))]
exp2<-floor(exp1)
delete<-apply(exp2,1,function(x) mean(x==0))
cou<-which(delete>0.75)####在超过75%样本以上的都是0的基因删去
exp<-exp2[(-cou),]
sensitivity<-GEM_res_S8495_sen
resistance<-GEM_res_S8495_res
resistance_exp<-exp[,match(resistance,colnames(exp))]
sensitivity_exp<-exp[,match(sensitivity,colnames(exp))]
geneid<-exp11[(-cou),1]


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
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/1_RNAseq")
write.table(res,"sensitivity-resistance-all-DESeq2_S8495_median_AUC_RNAseq_DEGs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
#resSig<-res[which(res$pvalue<0.05 & abs(res$log2FoldChange>1)),]
resSig<-res[which(res$pvalue<0.01),]

#resSig<-res[which(res$padj<0.05),]

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

setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/1_RNAseq")
write.table(up_gene,"sensitivity-resistance-S8495-DESeq2_median_AUC_RNAseq_P0.01_FC0_up_DEGs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(down_gene,"sensitivity-resistance-S8495-DESeq2_median_AUC_RNAseq_P0.01_FC0_down_DEGs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

#########  RNAseq DEGs Volcano Plot
rm(list=ls())
library(ggrepel) 
library(ggplot2)
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/1_RNAseq")
df<-read.table("sensitivity-resistance-all-DESeq2_S8495_median_AUC_RNAseq_DEGs.txt",sep = '\t',header= T,row.names = 1)
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
data = df[df$padj<0.05 & abs(df$log2FoldChange)>2,]
DEGs<-na.omit(as.matrix(data))

DEGs<-df[which(df$padj<0.01 & abs(df$log2FoldChange)>2),]
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/1_RNAseq/Volcano")
write.table(DEGs,"FDR05_FC2_DEGs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
ggsave("Volcano_FDR0.05_FC2.pdf",p,width = 8, height = 6)

###################################################################################
#3 ###############  DEGs与组蛋白修饰酶的交叠
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDAC_frontiers/5_histone-modifying enzymes")
Acetylase<-read.table("histone-modifying-enzymes.txt",header=T,sep="\t")
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/1_RNAseq")
df<-read.table("sensitivity-resistance-all-DESeq2_S8495_median_AUC_RNAseq_DEGs.txt",sep = '\t',header= T,row.names = 1)
DEGs1<-df[df$pvalue<0.05,]
#DEGs<-DEGs1[-which(is.na(DEGs1$pvalue)),]
DEGs<-DEGs1
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
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/1_RNAseq/Acetylase")
write.table(Acetylase_info,"Acetylase_info_P0.05.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
library(VennDiagram)
venn.diagram(x=list(Acetylase=Acetylase_name,DEGs=DEGs_name),cex = 1,margin = 0.1, "Acetylase_info_P005.png",fill=c("red","blue"))

setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
library("data.table")
exp<-fread("star_rsem.GeneSymbol.TPM.xls",header=T,data.table=F)
sensitivity<-c("DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-12","DAC-14")
resistance<-c("DAC-20","DAC-21","DAC-22","DAC-23","DAC-25","DAC-29","DAC-30","DAC-32","DAC-35","DAC-5","DAC-13","DAC-39")


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
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/1_RNAseq/Acetylase")
write.table(Acetylase_info_value,"sen_or_res_to_S8495_Acetylase_info_P0.05_value.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

library(ggplot2)
df2<-cbind(rbind(as.matrix(Acetylase_info_value[,1]),as.matrix(Acetylase_info_value[,1])),
           rbind(as.matrix(sensitivity_exp11),as.matrix(resistance_exp11)),
           rbind(as.matrix(rep("sensitive",length(Acetylase_info_value[,1]))),as.matrix(rep("resistant",length(Acetylase_info_value[,1])))),
           rbind(as.matrix(Acetylase_info_value$sensitive+sensitivity_sd),as.matrix(Acetylase_info_value$resistant+resistance_sd)),
           rbind(as.matrix(Acetylase_info_value$sensitive-sensitivity_sd),as.matrix(Acetylase_info_value$resistant-resistance_sd)),
           rbind(as.matrix(Acetylase_info_value[,8]),as.matrix(Acetylase_info_value[,8])))
colnames(df2)<-c("gene","value","type","valuesd1","valuesd2","pvalue")
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/1_RNAseq/Acetylase")
write.table(df2,"sen_or_res_to_S8495_Acetylase_info_P0.05_value_plot.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

df3<-read.table("sen_or_res_to_S8495_Acetylase_info_P0.05_value_plot.txt",sep = '\t',header= T)
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
  geom_text(aes(y=30,label = pva),size = 3)
plot(p3)

setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/1_RNAseq/Acetylase")
ggsave("sen_or_res_to_S8495_Acetylase_info_P0.05_value_plot.pdf",p3,width = 10, height = 6)


#4 ###########  clusterProfiler  ####################################
rm(list=ls())
library(org.Hs.eg.db)
library(clusterProfiler)
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/1_RNAseq")
up_gene<-read.table("sensitivity-resistance-S8495-DESeq2_median_AUC_RNAseq_P0.05_FC0_down_DEGs.txt",header=T,sep="\t")
library(stringr)
gene=bitr(up_gene[,1],fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
kegg<-enrichKEGG(gene[,2],organism = 'hsa',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'BH',
                 minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.05,use_internal_data = FALSE)
keu<-barplot(kegg,showCategory=15,drop=T,x = "GeneRatio",color = "pvalue")

ee<-kegg@result
ee1<-ee[which(ee$pvalue<0.05),]
#ee1<-ee[which(ee$p.adjust<0.05),]
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

setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/1_RNAseq/pathway")
write.table(ee1,"kegg_DEGP0.05_pathwayP0.05_down_drug_target.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(pathway_gene2,"kegg_DEGP0.05_pathwayP0.05_down_drug_target_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
ggsave("kegg_DEGP0.05_pathwayP0.05_down_drug_target_gene.pdf",keu,width = 8, height = 6)

go<-enrichGO(gene[,2],OrgDb = org.Hs.eg.db,ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.05,keyType = 'ENTREZID')
ked<-dotplot(go,showCategory=10)
a<-go@result
go_BP<-a[which(a$p.adjust<0.05),]

enrich_genego<-go_BP$geneID
pathway_genego<-unique(unlist(strsplit(enrich_genego,split="/")))
pathway_genego2=bitr(pathway_genego,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495/1_RNAseq/pathway")
write.table(pathway_genego2,"GO_DEGP0.05_pathway_FDR0.05_down_drug_target_gene.txt",sep = '\t',col.names = F,row.names = F,quote = FALSE)##数据输出
write.table(go_BP,"GO_DEGP0.05_pathway_FDR0.05_down_drug_target.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
ggsave("GO_DEGP0.05_pathway_FDR0.05_down_drug_target_gene.pdf",ked,width = 8, height = 6)

#9 ######################    ATACseq   ################################################
#############  DApeaks
rm(list=ls())
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_count_name<-read.table("peak_count.txt",sep="\t",header=T) #没有全0的
GEM_res_S8495_sen<-c("DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-12","DAC-14")
GEM_res_S8495_res<-c("DAC-20","DAC-21","DAC-22","DAC-23","DAC-25","DAC-29","DAC-30","DAC-32","DAC-35","DAC-5","DAC-13","DAC-39")

sensitivity<-GEM_res_S8495_sen
resistance<-GEM_res_S8495_res

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
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq")
write.table(res,"sensitivity-resistance-all-DESeq2_S8495_median_AUC_ATAC.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

resSig<-res[which(res$pvalue<0.05 & abs(res$log2FoldChange)>2),]
resSig<-res[which(res$pvalue<0.01),]
resSig<-res[which(res$pvalue<0.005),]

resSig<-res[which(res$padj<0.05),] #没有
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
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq")
out_file<-"sensitivity-resistance-S8495-median-AUC-P0.05-logFC2"
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
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495/2_ATACseq")
write.table(peaks_bed_up,file="all_peaks.bed",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')


#10 ######################################################################################
####  homer注释完的信息,画饼图
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq/annotation")
library(data.table)
homer_anno<-fread("sensitivity-resistance-S8495-median-AUC-P0.01-logFC0-DESeq2-down-annotation.txt",header=T,data.table=F)
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
  labs(title = "down-peaks-p0.01")+
  theme(plot.title = element_text(hjust = 0.5))
#geom_text(aes(y = A/2 + c(0, cumsum(A)[-length(A)]), x = sum(A)/20, label = myLabel), size = 5)   # 在图中加上百分比：x 调节标签到圆心的距离, y 调节标签的左右位置
print(p) #显示饼图
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq/annotation")
ggsave("down_P0.01_FC0_annotation.pdf",p,width = 6, height = 5)
write.table(result,file="down_P0.01_FC0_annotation_ratio.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)
write.table(result,file="all_annotation_ratio.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)

########   堆积图   ###################
rm(list=ls())
library(reshape2)
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495/2_ATACseq/annotation")
all<-read.delim("all_annotation_ratio.txt",sep = '\t',header = T)
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq/annotation")
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
plot(p)
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq/annotation")
ggsave(p, filename = 'annotation-ratio-all-down-up.pdf', width = 8, height = 6, dpi = 600)

#11 ##############  peaks都处于什么位置  ######################    
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq/annotation")
library(data.table)
#homer_anno1<-fread("all_peaks-annotation.txt",header=T,data.table=F)

homer_anno1<-fread("sensitivity-resistance-S8495-median-AUC-P0.01-logFC0-DESeq2-up-annotation.txt",header=T,data.table=F)
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq")
res<-read.table("sensitivity-resistance-all-DESeq2_S8495_median_AUC_ATAC.txt",header = TRUE,sep = "\t")
resSig<-res[which(res$pvalue<0.01),]
int_peak<-intersect(resSig[,1],homer_anno1[,1])

resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
up_gene<-as.data.frame(resSig[which(resSig$log2FoldChange>0),1])
down_gene<-as.data.frame(resSig[which(resSig$log2FoldChange<0),1])
int_peak<-intersect(as.character(up_gene[,1]),homer_anno1[,1])#up
int_peak<-intersect(as.character(down_gene[,1]),homer_anno1[,1])#down
int_peak<-homer_anno1[,1]#all

homer_anno<-homer_anno1[match(int_peak,homer_anno1[,1]),]
one<-sum(abs(homer_anno$`Distance to TSS`)<=100)
two<-length(which(abs(homer_anno$`Distance to TSS`)>=100 & abs(homer_anno$`Distance to TSS`)<=1000))
three<-length(which(abs(homer_anno$`Distance to TSS`)>=1000 & abs(homer_anno$`Distance to TSS`)<=10000))
four<-length(which(abs(homer_anno$`Distance to TSS`)>=10000 & abs(homer_anno$`Distance to TSS`)<=100000))
five<-sum(abs(homer_anno$`Distance to TSS`)>=100000)
ma<-matrix(rep("uppeaks-peaks",5),ncol=1)
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
  labs(x = "", y = "", title = "uppeaks-position") +  # 将横纵坐标的标签设为空
  theme(axis.ticks = element_blank()) +  # 将左上角边框的刻度去掉
  theme(legend.title = element_blank(), legend.position = "left")+   ## 将图例标题设为空，并把图例方放在左边位置
  scale_fill_discrete(breaks = dt$B, labels = myLabel)+   # 将原来的图例标签换成现在的myLabel
  theme(axis.text.x = element_blank())+   ## 去掉饼图的外框上的数值，即去除原柱状图的X轴，把X轴的刻度文字去掉
  theme(plot.title = element_text(hjust = 0.5))
#geom_text(aes(y = A/2 + c(0, cumsum(A)[-length(A)]), x = sum(A)/20, label = myLabel), size = 5)   # 在图中加上百分比：x 调节标签到圆心的距离, y 调节标签的左右位置
print(p) #显示饼图
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq/annotation")
ggsave("uppeaks-P0.01-position.pdf",p,width = 6, height = 5)
write.table(cbind(data,dt),file="uppeaks-P0.01-position.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')
ggsave("allpeaks-position.pdf",p,width = 6, height = 5)
write.table(cbind(data,dt),file="allpeaks-position.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')

###### 堆积图
rm(list=ls())
library(reshape2)
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq/annotation")
all<-read.delim("allpeaks-position.txt",sep = '\t',header = T)
#DA<-read.delim("DApeaks_P0.01_FC0_position.txt",sep = '\t',header = T)
down<-read.delim("downpeaks-P0.01-position.txt",sep = '\t',header = T)
up<-read.delim("uppeaks-P0.01-position.txt",sep = '\t',header = T)
dat_m<-rbind(all[,c(1,2,4)],down[,c(1,2,4)],up[,c(1,2,4)])
colnames(dat_m) = c('Group','Type','value')
#定义`Group`列的出图顺序
dat_m$Group = factor(dat_m$Group, levels = c("allpeaks","downpeaks","uppeaks"))
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
plot(p)
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq/annotation")
ggsave(p, filename = 'Distance to the closest transcription start site.pdf', width = 6, height = 5, dpi = 600)

#12 ####################  符合在TSS100KB之内条件的注释基因
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq/annotation")
library(data.table)
homer_anno_down<-fread("sensitivity-resistance-S8495-median-AUC-P0.01-logFC0-DESeq2-down-annotation.txt",header=T,data.table=F)
colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]
length(unique(Anno_gene_100Kb_down$`Gene Name`))

homer_anno_up<-fread("sensitivity-resistance-S8495-median-AUC-P0.01-logFC0-DESeq2-up-annotation.txt",header=T,data.table=F)
colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
length(unique(Anno_gene_100Kb_up$`Gene Name`))

setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/1_RNAseq")
up_gene<-read.table("sensitivity-resistance-S8495-DESeq2_median_AUC_RNAseq_P0.05_FC0_up_DEGs.txt",header=T,sep="\t")
down_gene<-read.table("sensitivity-resistance-S8495-DESeq2_median_AUC_RNAseq_P0.05_FC0_down_DEGs.txt",header=T,sep="\t")
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

setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/1_2_intersect/venn")
library(VennDiagram)
venn.diagram(x=list(downpeaks=unique(down_peaks_gene),uppeaks=unique(up_peaks_gene),RNAseq=unique(DEGs)), "DEGP005_DApeakP001.png",fill=c("red","green","blue"),margin = 0.1)

down_peaks_gene1<-as.data.frame(unique(down_peaks_gene))
colnames(down_peaks_gene1)<-"down_peaks_gene"
up_peaks_gene1<-as.data.frame(unique(up_peaks_gene))
colnames(up_peaks_gene1)<-"up_peaks_gene"
write.table(down_peaks_gene1,"down_peaks0.01_annotation_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(up_peaks_gene1,"up_peaks0.01_annotation_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

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

setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
library("data.table")
exp11<-fread("star_rsem.GeneSymbol.readscount.xls",header=T,data.table=F)
GEM_res_S8495_sen<-c("DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-12","DAC-14")
GEM_res_S8495_res<-c("DAC-20","DAC-21","DAC-22","DAC-23","DAC-25","DAC-29","DAC-30","DAC-32","DAC-35","DAC-5","DAC-13","DAC-39")
all<-c(GEM_res_S8495_sen,GEM_res_S8495_res)
exp1<-exp11[,match(all,colnames(exp11))]
exp2<-floor(exp1)
delete<-apply(exp2,1,function(x) mean(x==0))
cou<-which(delete>0.75)####在超过75%样本以上的都是0的基因删去
exp<-exp2[(-cou),]

Ku=nrow(int_up)
Nu=nrow(exp) #all gene
Mu=length(unique(Anno_gene_100Kb_up$`Gene Name`)) #all DARs annotation gene
nu=nrow(up_gene) #all DEGs
pu<-phyper(Ku-1,Mu, Nu-Mu, nu, lower.tail=F)
pu
Kd=nrow(int_down)
Nd=nrow(exp) #all gene
Md=length(unique(Anno_gene_100Kb_down$`Gene Name`)) #all DARs annotation gene
nd=nrow(down_gene) #all DEGs
pd<-phyper(Kd-1,Md, Nd-Md, nd, lower.tail=F)
pd
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/1_2_intersect/venn")
write.table(int_up,"int_up_peak0.01_all_DEG0.05_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(int_down,"int_down_peak0.01_all_DEG0.05_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(int_up_peak,"int_up_peak0.01_all_DEG0.05_peak.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(int_down_peak,"int_down_peak0.01_all_DEG0.05_peak.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(up_DEGs_up_DApeaks,"int_up_peak0.01_up_DEG0.05_peak.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(down_DEGs_down_DApeaks,"int_down_peak0.01_down_DEG0.05_peak.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出


#13 #################################################################################
##############  DEpeaks在不同样本变化的倍数与DEGs在不同样本变化的倍数之间的相关性
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq/annotation")
library(data.table)
homer_anno_down<-fread("sensitivity-resistance-S8495-median-AUC-P0.01-logFC0-DESeq2-down-annotation.txt",header=T,data.table=F)
colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]
length(unique(Anno_gene_100Kb_down$`Gene Name`))

homer_anno_up<-fread("sensitivity-resistance-S8495-median-AUC-P0.01-logFC0-DESeq2-up-annotation.txt",header=T,data.table=F)
colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
length(unique(Anno_gene_100Kb_up$`Gene Name`))

setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/1_RNAseq")
up_gene<-read.table("sensitivity-resistance-S8495-DESeq2_median_AUC_RNAseq_P0.05_FC0_up_DEGs.txt",header=T,sep="\t")
down_gene<-read.table("sensitivity-resistance-S8495-DESeq2_median_AUC_RNAseq_P0.05_FC0_down_DEGs.txt",header=T,sep="\t")
all<-rbind(up_gene,down_gene)
DEGs<-as.character(all[,1])
down_DEGs_DApeaks<-intersect(Anno_gene_100Kb_down$`Gene Name`,DEGs)
up_DEGs_DApeaks<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)

int_all_DEGs<-unique(c(down_DEGs_DApeaks,up_DEGs_DApeaks))
all_peaks_anno<-rbind(homer_anno_down,homer_anno_up)
#DEGs_peaks<-all_peaks_anno[all_peaks_anno$`Gene Name` %in% int_all_DEGs,]
DEGs_peaks<-all_peaks_anno[match(int_all_DEGs,all_peaks_anno$`Gene Name`),1] #match只取出对应上的第一个peak

DEGs_peaks1<-all_peaks_anno[all_peaks_anno$`Gene Name` %in% int_all_DEGs,]
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq")
DApeak<-fread("sensitivity-resistance-all-DESeq2_S8495_median_AUC_ATAC.txt",header=T,data.table=F)
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
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/1_2_intersect/correlation")
write.table(dat1,"DEGs_FC_DARs_FC_correlation.txt",sep = '\t',col.names = T,row.names = T,quote = FALSE)##数据输出

library(ggplot2)
library(ggpubr)
a11<-c(1:10)
a21<-c((nrow(dat1)-9):nrow(dat1))
a3<-c(a11,a21)
want_dat<-dat1[order(dat1[,1])[a3],]
library(ggrepel)
ps<-ggplot(data=dat1, aes(x=DEGslog2FC, y=DApeakslog2FC))+geom_point(color="red",size=0.1)+stat_smooth(method="lm",se=FALSE)+stat_cor(data=dat1, method = "spearman")+
  ggtitle("spearman-DApeaksP0.01-DEGsP0.05") +
  theme(plot.title = element_text(hjust = 0.5))+
  geom_text_repel(
    data = want_dat[,c(1:2)],
    aes(label = rownames(want_dat)),
    size = 3,
    color = "black",
    segment.color = "black", show.legend = FALSE )
plot(ps)
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/1_2_intersect/correlation")
ggsave("correlation-DApeaksP0.01-DEGsP0.05-spearman-S8495-onlyone.pdf",ps,width = 10, height = 10)

#####################################################################
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq/annotation")
library(data.table)
homer_anno_down<-fread("sensitivity-resistance-S8495-median-AUC-P0.01-logFC0-DESeq2-down-annotation.txt",header=T,data.table=F)
colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]
length(unique(Anno_gene_100Kb_down$`Gene Name`))

homer_anno_up<-fread("sensitivity-resistance-S8495-median-AUC-P0.01-logFC0-DESeq2-up-annotation.txt",header=T,data.table=F)
colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
length(unique(Anno_gene_100Kb_up$`Gene Name`))

setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/1_RNAseq")
up_gene<-read.table("sensitivity-resistance-S8495-DESeq2_median_AUC_RNAseq_P0.05_FC0_up_DEGs.txt",header=T,sep="\t")
down_gene<-read.table("sensitivity-resistance-S8495-DESeq2_median_AUC_RNAseq_P0.05_FC0_down_DEGs.txt",header=T,sep="\t")
all<-rbind(up_gene,down_gene)
DEGs<-as.character(all[,1])
down_DEGs_DApeaks<-intersect(Anno_gene_100Kb_down$`Gene Name`,DEGs)
up_DEGs_DApeaks<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)

int_all_DEGs<-unique(c(down_DEGs_DApeaks,up_DEGs_DApeaks))
all_peaks_anno<-rbind(homer_anno_down,homer_anno_up)
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq")
DApeak<-fread("sensitivity-resistance-all-DESeq2_S8495_median_AUC_ATAC.txt",header=T,data.table=F)


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
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/1_2_intersect/correlation")
write.table(DEG_DApeak_FC,"DEGs_FC_DARs_FC_correlation_comprehensive.txt",sep = '\t',col.names = T,row.names = T,quote = FALSE)##数据输出

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
plot(ps)
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/1_2_intersect/correlation")
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

#14 ##################  对高低可及peaks-DEGs进行KEGG通路富集分析
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq/annotation")
library(data.table)
homer_anno_down<-fread("sensitivity-resistance-S8495-median-AUC-P0.01-logFC0-DESeq2-down-annotation.txt",header=T,data.table=F)
colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]
length(unique(Anno_gene_100Kb_down$`Gene Name`))

homer_anno_up<-fread("sensitivity-resistance-S8495-median-AUC-P0.01-logFC0-DESeq2-up-annotation.txt",header=T,data.table=F)
colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
length(unique(Anno_gene_100Kb_up$`Gene Name`))

setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/1_RNAseq")
up_gene<-read.table("sensitivity-resistance-S8495-DESeq2_median_AUC_RNAseq_P0.05_FC0_up_DEGs.txt",header=T,sep="\t")
down_gene<-read.table("sensitivity-resistance-S8495-DESeq2_median_AUC_RNAseq_P0.05_FC0_down_DEGs.txt",header=T,sep="\t")
#all<-rbind(up_gene,down_gene)
#DEGs<-as.character(all[,1])
down_DEGs_DApeaks<-intersect(Anno_gene_100Kb_down$`Gene Name`,down_gene[,1])
up_DEGs_DApeaks<-intersect(Anno_gene_100Kb_up$`Gene Name`,up_gene[,1])

library(org.Hs.eg.db)
library(clusterProfiler)
up_gene<-up_DEGs_DApeaks  ## 修改
gene=bitr(up_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
## 去重
library(stringr)
kegg<-enrichKEGG(gene[,2],organism = 'hsa',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'none',
                 minGSSize = 10,maxGSSize = 500,qvalueCutoff = 1,use_internal_data = FALSE)
ee<-kegg@result
ee1<-ee[which(ee$pvalue<0.05),]
#ee1<-ee[which(ee$p.adjust<0.05),]

p<-barplot(kegg,showCategory=15,drop=T,x = "GeneRatio",color = "pvalue")
plot(p)
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
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/1_2_intersect/pathway")
kegg_dysregulation<-"kegg_pathway_P0.05_up"
write.table(ee1,paste("DEGP0.05_DApeaksP0.01_",kegg_dysregulation,".txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(pathway_gene,paste("DEGP0.05_DApeaksP0.01_",kegg_dysregulation,"_gene.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
ggsave(paste("DEGP0.05_DApeaksP0.01_",kegg_dysregulation,".pdf",sep=""),p,width = 5, height = 3)
write.table(result2,paste("DEGP0.05_DApeaksP0.01_",kegg_dysregulation,"_gene_symbol.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(a4,paste("DEGP0.05_DApeaksP0.01_",kegg_dysregulation,"_pathway_symbol.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

#15 #############  ridge regression  #####因为FIMO能知道peaks与TFs的关系
rm(list=ls())
library(dplyr)
library(data.table)
library(glmnet)
library(foreach)
library(ggplot2)

setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq/FIMO/P0.01_FC0_down")
fimo_down<-fread("fimo.tsv",header=T,data.table=F,fill=TRUE)
colnames(fimo_down)<-c(colnames(fimo_down)[1:8],"FDR","matched_sequence")
fimo_d<-fimo_down[which(fimo_down$FDR<0.05),]
peak_d<-unique(fimo_d[,3])
TF_namesd<-toupper(apply(as.data.frame(fimo_d[,2]),1,function(x) str_replace_all(as.character(x),"\\((.*?)\\)","")))
TF_d<-unique(TF_namesd)
#TF_d<-unique(fimo_d[,2])
X_matrix_d<-matrix(0,length(peak_d),length(TF_d))
colnames(X_matrix_d)<-TF_d
rownames(X_matrix_d)<-peak_d
for(i in 1:length(TF_d)){
  TF_peak_d<-fimo_d[fimo_d[,2] %in% TF_d[i],3]
  X_matrix_d[match(unique(TF_peak_d),peak_d),i]<-1
}
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq")
DApeak<-fread("sensitivity-resistance-all-DESeq2_S8495_median_AUC_ATAC.txt",header=T,data.table=F)
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/1_RNAseq")
DE_gene<-fread("sensitivity-resistance-all-DESeq2_S8495_median_AUC_RNAseq_DEGs.txt",header=T,data.table=F)

a_d<-gsub(":","_",peak_d)
b_d<-gsub("-","_",a_d)
Y_d<-DApeak[match(b_d,DApeak[,1]),3]

##### ridge regression
y<-Y_d
x<-X_matrix_d
#不断尝试和调整λ，找到最优解
lambdas <- 10^seq(2, -2, by = -.1)
##使用cv.glmnet()函数来进行参数lambda的影响分析。通过指定参数alpha = 0来建立Ridge回归，如果参数alpha = 1，则建立的是Lasso回归模型，nfolds =3表示使用3折交叉验证。
fit = glmnet(x,y,alpha = 0,family = "gaussian",standardize=TRUE)
cv.fit <- cv.glmnet(x,y,alpha = 0,family = 'gaussian',standardize=TRUE,grouped=FALSE,nfolds = 5,lambda = lambdas)
pd<-plot(cv.fit)
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq/FIMO/P0.01_FC0_down")
ggsave("ridge_regression_lambda_down.pdf",pd,width = 10, height = 4)
y_predicted <- predict(cv.fit,s=cv.fit$lambda.min,newx=x)

# Sum of Squares Total and Error
sst <- sum((y - mean(y))^2)
sse <- sum((y_predicted - y)^2)
# R squared
rsq_d <- round(1 - sse / sst,3)
### 检查summary，看输出的数据是否正确
tmp_coeffs <- coef(cv.fit, s = "lambda.min")
output_coef=data.frame(features = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
output_coef<-output_coef[order(output_coef[,"coefficient"],decreasing=TRUE),]
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq/FIMO/P0.01_FC0_down")
write.table(output_coef,file="ridge_regression_down.txt",sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')
res_positive1<-output_coef[which(output_coef$coefficient>0),]

resSig<-DE_gene[which(DE_gene$pvalue<0.05),1]
down_int<-intersect(as.character(res_positive1[,1]),resSig)
#######
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq/FIMO/P0.01_FC0_up")
fimo_up<-fread("fimo.tsv",header=T,data.table=F,fill=TRUE)
colnames(fimo_up)<-c(colnames(fimo_up)[1:8],"FDR","matched_sequence")
fimo_u<-fimo_up[which(fimo_up$FDR<0.05),]
peak_u<-unique(fimo_u[,3])
TF_namesu<-toupper(apply(as.data.frame(fimo_u[,2]),1,function(x) str_replace_all(as.character(x),"\\((.*?)\\)","")))
TF_u<-unique(TF_namesu)
X_matrix_u<-matrix(0,length(peak_u),length(TF_u))
colnames(X_matrix_u)<-TF_u
rownames(X_matrix_u)<-peak_u
for(j in 1:length(TF_u)){
  TF_peak_u<-fimo_u[fimo_u[,2] %in% TF_u[j],3]
  X_matrix_u[match(unique(TF_peak_u),peak_u),j]<-1
}
a_u<-gsub(":","_",peak_u)
b_u<-gsub("-","_",a_u)
Y_u<-DApeak[match(b_u,DApeak[,1]),3]
##### ridge regression
y<-Y_u
x<-X_matrix_u
lambdas <- 10^seq(2, -2, by = -.1)
fit_u = glmnet(x,y,alpha = 0,family = "gaussian",standardize=TRUE)
cv.fit_u <- cv.glmnet(x,y,alpha = 0,family = 'gaussian',standardize=TRUE,grouped=FALSE,nfolds = 5,lambda = lambdas)
pu<-plot(cv.fit_u)
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq/FIMO/P0.01_FC0_up")
ggsave("ridge_regression_lambda_up.pdf",pu,width = 10, height = 4)
y_predicted_u <- predict(cv.fit_u,s=cv.fit_u$lambda.min,newx=x)
# Sum of Squares Total and Error
sst <- sum((y - mean(y))^2)
sse <- sum((y_predicted_u - y)^2)
# R squared
rsq_u <- round(1 - sse / sst,3)
### 检查summary，看输出的数据是否正确
tmp_coeffs_u <- coef(cv.fit_u, s = "lambda.min")
output_coef_u=data.frame(features = tmp_coeffs_u@Dimnames[[1]][tmp_coeffs_u@i + 1], coefficient = tmp_coeffs_u@x)
output_coef_u<-output_coef_u[order(output_coef_u[,"coefficient"],decreasing=TRUE),]
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq/FIMO/P0.01_FC0_up")
write.table(output_coef_u,file="ridge_regression_up.txt",sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')
res_positive1_u<-output_coef_u[which(output_coef_u$coefficient>0),]
up_int<-intersect(res_positive1_u,resSig)

##################   画回归系数柱状图  #####################
res_positive1_u[,1]<-apply(as.data.frame(res_positive1_u[,1]),1,function(x) str_replace_all(as.character(x),"\\((.*?)\\)",""))
res_positive1[,1]<-apply(as.data.frame(res_positive1[,1]),1,function(x) str_replace_all(as.character(x),"\\((.*?)\\)",""))

int_TF<-c(intersect(res_positive1_u[,1],res_positive1[,1]),"(Intercept)","")
up_TFs1<-as.data.frame(res_positive1_u[-(match(int_TF,res_positive1_u[,1])[which(!is.na(match(int_TF,res_positive1_u[,1])))]),])
down_TFs1<-as.data.frame(res_positive1[-(match(int_TF,res_positive1[,1])[which(!is.na(match(int_TF,res_positive1[,1])))]),])
## 指定绘图顺序
up_TFs1$features <- factor(up_TFs1$features,
                           levels = rev(unique(up_TFs1$features)),
                           ordered = T)
down_TFs1$features <- factor(down_TFs1$features,
                             levels = rev(unique(down_TFs1$features)),
                             ordered = T)
## 绘制右侧的条形图
left_1 <- ggplot(up_TFs1,aes(x=features,y=coefficient*(-1)))+ #*-1 change direction
  coord_flip()+
  labs(x="",y="UP (sen vs res)",title="")+
  geom_bar(fill="red",colour="white", #orange
           size=0,width = 0.8,
           stat = "identity")
#left_1
## 绘制左侧的条形图
right_1 <- ggplot(down_TFs1,aes(x=features,y=coefficient))+
  coord_flip()+labs(x="",y="DOWN (sen vs res)",title="")+
  geom_bar(fill="blue",colour="white", #purple
           size=0,width = 0.8,stat = "identity")
#right_1
#设置纵轴的位置改为右侧；
right_2<- right_1+scale_x_discrete(position = "top")
#right_2
## 图表的组合,使用patchwork包将两个图合并成一个；
library(patchwork)
pu<-left_1+right_2+labs(title="Regression coefficients showing enrichmen of TFs")+
  theme(plot.title = element_text(hjust = 0.5))
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq/FIMO")
ggsave(paste("ridge_regression_","DApeakP0.01","TF_coefficient.pdf",sep=""),pu,width = 8, height = 6)
write.table(up_TFs1,"FIMO_result_up_TFs_motifFDR0.05_DApeakP0.01.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(down_TFs1,"FIMO_result_down_TFs_motifFDR0.05_DApeakP0.01.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

##################  ridge regresion 找到的TF及其靶标DEGs  上下调的转录因子是分开看的
###############   TF找到对应的靶基因
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq/FIMO")
up_TFs1<-read.table("FIMO_result_up_TFs_motifFDR0.05_DApeakP0.01.txt",sep = '\t',header=TRUE)
down_TFs1<-read.table("FIMO_result_down_TFs_motifFDR0.05_DApeakP0.01.txt",sep = '\t',header=TRUE)

DETFs<-as.character(up_TFs1[,1])
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
# TF-target
result<-NULL
for(i in 1:length(DETF_want)){
  result1<-target_TF_gene[target_TF_gene[,2] %in% DETF_want[i],]
  result<-rbind(result,result1)
}
### 查看靶标基因中有多少是DEGs
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/1_RNAseq")
up_gene<-read.table("sensitivity-resistance-S8495-DESeq2_median_AUC_RNAseq_P0.05_FC0_up_DEGs.txt",header=T,sep="\t")
down_gene<-read.table("sensitivity-resistance-S8495-DESeq2_median_AUC_RNAseq_P0.05_FC0_down_DEGs.txt",header=T,sep="\t")

data<-NULL
for(i in 1:length(unique(result[,2]))){
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
colnames(data)<-c("TFs","target-DEGs","target-DEGs-dysregulation")

aa<-unique(data[,1])
aa
re<-data.frame()
for(m in 1:length(unique(data[,1]))){
  dat<-data[data[,1]%in%aa[m],]
  if(class(dat)=="character"){
    up<-sum(dat[3]=="up")
    down<-sum(dat[3]=="down")
    re[m,1]<-aa[m]
    re[m,2]<-up
    re[m,3]<-down
  }
  if(class(dat)!="character"){
    up<-sum(dat[,3]=="up")
    down<-sum(dat[,3]=="down")
    re1<-matrix(c(aa[m],up,down),nrow=1)
    re<-rbind(re,re1)
    re[m,1]<-aa[m]
    re[m,2]<-up
    re[m,3]<-down
  }
}
colnames(re)<-c("geneid","up","down")
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq/TF_target")
write.table(data,"sen_res_down_TFs_TRANSFAC-gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(data,"sen_res_down_TFs_JASPAR-gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(data,"sen_res_down_TFs_MotifMap-gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(data,"sen_res_down_TFs_ENCODE-gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

write.table(re,"number_of_down_TRANSFAC-gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(re,"number_of_down_JASPAR-gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(re,"number_of_down_MotifMap-gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(re,"number_of_down_ENCODE-gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

#########################################################################
####图形展示TFs-靶基因情况
rm(list=ls())
library(ggplot2)
#1
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq/TF_target")
up_ENCODE <- read.delim('number_of_up_ENCODE-gene.txt',header = T,sep="\t")
down_ENCODE <- read.delim('number_of_down_ENCODE-gene.txt',header = T,sep="\t")
data1<-rbind(cbind(as.character(up_ENCODE[,1]),matrix(rep("up",nrow(up_ENCODE)),ncol=1),up_ENCODE[,2]),cbind(as.character(up_ENCODE[,1]),matrix(rep("down",nrow(up_ENCODE)),ncol=1),up_ENCODE[,3]))
data2<-rbind(cbind(as.character(down_ENCODE[,1]),matrix(rep("up",nrow(down_ENCODE)),ncol=1),down_ENCODE[,2]),cbind(as.character(down_ENCODE[,1]),matrix(rep("down",nrow(down_ENCODE)),ncol=1),down_ENCODE[,3]))
data3<-matrix(c(rep("up_ENCODE",(nrow(up_ENCODE)*2)),rep("down_ENCODE",(nrow(down_ENCODE)*2))),ncol=1)
data<-as.data.frame(cbind(data3,rbind(data1,data2)))
colnames(data)<-c("group","gene","type","value")

data$gene=factor(data$gene,levels=unique(data$gene))
p1<-ggplot(data,aes(x = gene, y = as.numeric(as.character(value)), fill = type))+
  stat_summary(geom = 'bar',fun = 'mean',cex=1.3,width=.6,position = position_dodge())+
  theme_bw()+
  geom_rect(aes(xmin=0.5,xmax=8.5,ymin=0,ymax=Inf),
            fill='grey80',color='grey80')+
  stat_summary(geom = 'bar',fun = 'mean',cex=1.3,width=.6,position = position_dodge())+
  #labs(x = "TFs", y = "Counts of TFs-target DEGs",title="ENCODE-TFs-Target-DEGs")+
  labs(x = "", y = "Counts")+
  scale_color_manual(values = c('#4a8a53','#941319'))+
  scale_fill_manual(values = c('#4a8a53','#941319'))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
#theme(plot.title=element_text(hjust=0.5))
plot(p1)
ggsave("ENCODE-TFs-Target-DEGs.pdf",p1,width =10,height =6,dpi = 600)
#2
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq/TF_target")
up_TRANSFAC <- read.delim('number_of_up_TRANSFAC-gene.txt',header = T,sep="\t")
down_TRANSFAC <- read.delim('number_of_down_TRANSFAC-gene.txt',header = T,sep="\t")
data1<-rbind(cbind(as.character(up_TRANSFAC[,1]),matrix(rep("up",nrow(up_TRANSFAC)),ncol=1),up_TRANSFAC[,2]),cbind(as.character(up_TRANSFAC[,1]),matrix(rep("down",nrow(up_TRANSFAC)),ncol=1),up_TRANSFAC[,3]))
data2<-rbind(cbind(as.character(down_TRANSFAC[,1]),matrix(rep("up",nrow(down_TRANSFAC)),ncol=1),down_TRANSFAC[,2]),cbind(as.character(down_TRANSFAC[,1]),matrix(rep("down",nrow(down_TRANSFAC)),ncol=1),down_TRANSFAC[,3]))
data3<-matrix(c(rep("up_TRANSFAC",(nrow(up_TRANSFAC)*2)),rep("down_TRANSFAC",(nrow(down_TRANSFAC)*2))),ncol=1)
data<-as.data.frame(cbind(data3,rbind(data1,data2)))
colnames(data)<-c("group","gene","type","value")

data$gene=factor(data$gene,levels=unique(data$gene))
p2<-ggplot(data,aes(x = gene, y = as.numeric(as.character(value)), fill = type))+
  stat_summary(geom = 'bar',fun = 'mean',cex=1.3,width=.6,position = position_dodge())+
  theme_bw()+
  geom_rect(aes(xmin=0.5,xmax=3.5,ymin=0,ymax=Inf),
            fill='grey80',color='grey80')+
  stat_summary(geom = 'bar',fun = 'mean',cex=1.3,width=.6,position = position_dodge())+
  #labs(x = "TFs", y = "Counts of TFs-target DEGs",title="TRANSFAC-TFs-Target-DEGs")+
  labs(x = "", y = "Counts")+
  scale_color_manual(values = c('#4a8a53','#941319'))+
  scale_fill_manual(values = c('#4a8a53','#941319'))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
#theme(plot.title=element_text(hjust=0.5))
plot(p2)
ggsave("TRANSFAC-TFs-Target-DEGs.pdf",p2,width =10,height =6,dpi = 600)

#3  
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq/TF_target")
up_JASPAR <- read.delim('number_of_up_JASPAR-gene.txt',header = T,sep="\t")
down_JASPAR <- read.delim('number_of_down_JASPAR-gene.txt',header = T,sep="\t")
data1<-rbind(cbind(as.character(up_JASPAR[,1]),matrix(rep("up",nrow(up_JASPAR)),ncol=1),up_JASPAR[,2]),cbind(as.character(up_JASPAR[,1]),matrix(rep("down",nrow(up_JASPAR)),ncol=1),up_JASPAR[,3]))
data2<-rbind(cbind(as.character(down_JASPAR[,1]),matrix(rep("up",nrow(down_JASPAR)),ncol=1),down_JASPAR[,2]),cbind(as.character(down_JASPAR[,1]),matrix(rep("down",nrow(down_JASPAR)),ncol=1),down_JASPAR[,3]))
data3<-matrix(c(rep("up_JASPAR",(nrow(up_JASPAR)*2)),rep("down_JASPAR",(nrow(down_JASPAR)*2))),ncol=1)
data<-as.data.frame(cbind(data3,rbind(data1,data2)))
colnames(data)<-c("group","gene","type","value")

data$gene=factor(data$gene,levels=unique(data$gene))
p3<-ggplot(data,aes(x = gene, y = as.numeric(as.character(value)), fill = type))+
  stat_summary(geom = 'bar',fun = 'mean',cex=1.3,width=.6,position = position_dodge())+
  theme_bw()+
  geom_rect(aes(xmin=0.5,xmax=5.5,ymin=0,ymax=Inf),
            fill='grey80',color='grey80')+
  stat_summary(geom = 'bar',fun = 'mean',cex=1.3,width=.6,position = position_dodge())+
  #labs(x = "TFs", y = "Counts of TFs-target DEGs",title="JASPAR-TFs-Target-DEGs")+
  labs(x = "", y = "Counts")+
  scale_color_manual(values = c('#4a8a53','#941319'))+
  scale_fill_manual(values = c('#4a8a53','#941319'))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
#theme(plot.title=element_text(hjust=0.5))
plot(p3)
ggsave("JASPAR-TFs-Target-DEGs.pdf",p3,width =10,height =6,dpi = 600)

#4  
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq/TF_target")
up_MotifMap <- read.delim('number_of_up_MotifMap-gene.txt',header = T,sep="\t")
down_MotifMap <- read.delim('number_of_down_MotifMap-gene.txt',header = T,sep="\t")
data1<-rbind(cbind(as.character(up_MotifMap[,1]),matrix(rep("up",nrow(up_MotifMap)),ncol=1),up_MotifMap[,2]),cbind(as.character(up_MotifMap[,1]),matrix(rep("down",nrow(up_MotifMap)),ncol=1),up_MotifMap[,3]))
data2<-rbind(cbind(as.character(down_MotifMap[,1]),matrix(rep("up",nrow(down_MotifMap)),ncol=1),down_MotifMap[,2]),cbind(as.character(down_MotifMap[,1]),matrix(rep("down",nrow(down_MotifMap)),ncol=1),down_MotifMap[,3]))
data3<-matrix(c(rep("up_MotifMap",(nrow(up_MotifMap)*2)),rep("down_MotifMap",(nrow(down_MotifMap)*2))),ncol=1)
data<-as.data.frame(cbind(data3,rbind(data1,data2)))
colnames(data)<-c("group","gene","type","value")

data$gene=factor(data$gene,levels=unique(data$gene))
p4<-ggplot(data,aes(x = gene, y = as.numeric(as.character(value)), fill = type))+
  stat_summary(geom = 'bar',fun = 'mean',cex=1.3,width=.6,position = position_dodge())+
  theme_bw()+
  geom_rect(aes(xmin=0.5,xmax=2.5,ymin=0,ymax=Inf),
            fill='grey80',color='grey80')+
  stat_summary(geom = 'bar',fun = 'mean',cex=1.3,width=.6,position = position_dodge())+
  #labs(x = "TFs", y = "Counts of TFs-target DEGs",title="MotifMap-TFs-Target-DEGs")+
  labs(x = "", y = "Counts")+
  scale_color_manual(values = c('#4a8a53','#941319'))+
  scale_fill_manual(values = c('#4a8a53','#941319'))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
#theme(plot.title=element_text(hjust=0.5))
plot(p4)
ggsave("MotifMap-TFs-Target-DEGs.pdf",p4,width =10,height =6,dpi = 600)


library(ggfortify)
library(ggpubr)
p<-ggarrange(p1, p2,p3,p4, labels = c("ENCODE", "TRANSFAC","JASPAR","MotifMap"),
             ncol = 2, nrow = 2,hjust = -0.5,
             vjust = 7,label.x = 0,
             label.y = 1.135,widths = 1,
             heights = 1,common.legend=TRUE,font.label = list(size = 10))
plot(p)
ggsave("TFs-Target-DEGs-four-database-horizontal.pdf",p,width =15,height =10,dpi = 600)


###################################### WGS #####################################
#########  mutation  统计学上的差异(fisher's exect test)和看频率(plot)，也看了TCGA的突变
rm(list=ls())
library("data.table")
library(dplyr)
library(stringr)
library(foreach)
setwd("~/xjj/WGS_CNV/mutation")
library("data.table")
mut<-fread("Organoid_mut_binary.txt",header=T,data.table=F)
mut1<-mut[,-1]
sensitivity<-c("DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-12","DAC-14")
resistance<-c("DAC-20","DAC-21","DAC-22","DAC-23","DAC-25","DAC-29","DAC-30","DAC-32","DAC-35","DAC-5","DAC-13","DAC-39")

mutation_result<-data.frame(matrix(0,1,3))
res_int<-intersect(resistance,colnames(mut))
resistance_mut<-mut[,match(res_int,colnames(mut))]
sen_int<-intersect(sensitivity,colnames(mut))
sensitivity_mut<-mut[,match(sen_int,colnames(mut))]
muty<-cbind(resistance_mut,sensitivity_mut)
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
data_result<-data[which(data[,2]<0.05),]
FDRvalue=p.adjust(data[,2],method="BH")
which(FDRvalue<0.05)###多重检验校正下，没有基因是差异的

if(nrow(data_result)>0){
  mutation_result[1,2]<-nrow(data_result)
  mutation_result[1,3]<-paste0(data_result[,1],collapse =",")
}
if(nrow(data_result)==0){
  mutation_result[1,2]<-0
  mutation_result[1,3]<-0
}
mutation_result[1,1]<-"GEM"
colnames(mutation_result)<-c("drug","The number of gene mutation","mutation gene")
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_GEM/3_WGS")
write.table(mutation_result,"mutation_fisher_resultP0.1.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(data,"mutation_fisher_result_all.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出


###############################################################
#############   CNV landscape T-test 
rm(list=ls())
library("data.table")
library(dplyr)
library(stringr)
library(foreach)
setwd("~/xjj/WGS_CNV/CNV_FACETS")
library("data.table")
mut<-fread("CNV_FACETS_median_CNA_matrix.txt",header=T,data.table=F)# -2,2 在
mut1<-mut[,-1]

sensitivity<-c("DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-12","DAC-14")
resistance<-c("DAC-20","DAC-21","DAC-22","DAC-23","DAC-25","DAC-29","DAC-30","DAC-32","DAC-35","DAC-5","DAC-13","DAC-39")
res_int<-intersect(resistance,colnames(mut))
resistance_mut<-mut[,match(res_int,colnames(mut))]
sen_int<-intersect(sensitivity,colnames(mut))
sensitivity_mut<-mut[,match(sen_int,colnames(mut))]
for(k in 1:ncol(sensitivity_mut)){
  sensitivity_mut[which(sensitivity_mut[,k]!=1),k]<-0
}
for(l in 1:ncol(resistance_mut)){
  resistance_mut[which(resistance_mut[,l]!=1),l]<-0
}
geneid<-mut[,1]
tresult=matrix(0,length(geneid),4)
for (i in 1:length(geneid)){
  ttest=t.test(sensitivity_mut[i,],resistance_mut[i,])
  tresult[i,1:3]=c(geneid[i],ttest$statistic,ttest$p.value)
}
tresult[,4]=p.adjust(tresult[,3],method="BH")
colnames(tresult)<-c("geneid","statistic","p.value","FDR")
a1<-tresult[which(tresult[,3]<0.05),] #pvalue
a1<-tresult[which(tresult[,4]<0.05),] #fdrvalue 没有差异基因
up<-a1[which(a1[,2]>0),1] #sensitivity
down<-a1[which(a1[,2]<0),1] #resistance
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_GEM/3_WGS/CNV")
write.table(a1,"cnv_result_T_test.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出


#############   CNV landscape Fisher's exect test 
rm(list=ls())
library("data.table")
library(dplyr)
library(stringr)
library(foreach)
setwd("~/xjj/WGS_CNV/CNV_FACETS")
library("data.table")
mut<-fread("CNV_FACETS_median_CNA_matrix.txt",header=T,data.table=F)# -2,2 在
mut1<-mut[,-1]
sensitivity<-c("DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-12","DAC-14")
resistance<-c("DAC-20","DAC-21","DAC-22","DAC-23","DAC-25","DAC-29","DAC-30","DAC-32","DAC-35","DAC-5","DAC-13","DAC-39")

res_int<-intersect(resistance,colnames(mut))
resistance_mut<-mut[,match(res_int,colnames(mut))]
sen_int<-intersect(sensitivity,colnames(mut))
sensitivity_mut<-mut[,match(sen_int,colnames(mut))]

muty<-cbind(resistance_mut,sensitivity_mut)
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
data_result<-data[which(data[,2]<0.05),]  
colnames(data_result)<-c("gene","pvalue")
FDRvalue=p.adjust(data[,2],method="BH")
which(FDRvalue<0.05)###多重检验校正下，没有基因是差异的


setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_GEM/3_WGS/CNV")
write.table(data_result,"CNV_fisher_resultP0.05.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(data,"CNV_fisher_result_all.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出


##################  DEG DApeaks circos
####0.加载工具包####
#画基因组圈图
#BiocManager::install("circlize")
####0.加载工具包####
#画基因组圈图
library(circlize)
#读取gtf文件来注释差异基因
library("rtracklayer")

####1.加载依赖数据####
#这里使用import导入gtf文件， 生成一个GRangs对象
gtf_data = import('~/xjj/CUP/gtf/Homo_sapiens.GRCh37.75.gtf') #gtf的路径
gtf_data = as.data.frame(gtf_data)

### ATAC
library("data.table")
df1<-read.table("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq/sensitivity-resistance-all-DESeq2_S8495_median_AUC_ATAC.txt",sep = '\t',header= T)
resSig1<-df1[df1$pvalue<0.01,]
up_peak<-resSig1[which(resSig1$log2FoldChange>0),]
down_peak<-resSig1[which(resSig1$log2FoldChange<0),]
library(dplyr)
peak_up<-as.data.frame(t(apply(as.data.frame(up_peak$peak_id),1,function(x) unlist(strsplit(as.character(x), "_")))))
peak_down<-as.data.frame(t(apply(as.data.frame(down_peak$peak_id),1,function(x) unlist(strsplit(as.character(x), "_")))))

up_logFC = up_peak$log2FoldChange
down_logFC =down_peak$log2FoldChange
bed11=as.data.frame(cbind(peak_up,up_logFC))
colnames(bed11)<-c("chr","start","end","value1")
bed11[,2]<-as.integer(as.character(bed11[,2]))
bed11[,3]<-as.integer(as.character(bed11[,3]))

bed22=as.data.frame(cbind(peak_down,down_logFC))
colnames(bed22)<-c("chr","start","end","value1")
bed22[,2]<-as.integer(as.character(bed22[,2]))
bed22[,3]<-as.integer(as.character(bed22[,3]))
#黑色下调，红上调
set.seed(60)
circos.initializeWithIdeogram(species = "hg19",plotType = c("axis"))
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(1))
  circos.text(mean(xlim), mean(ylim), chr, cex = 0.7, col = "white",
              facing = "inside", niceFacing = TRUE)
}, track.height = 0.15, bg.border = NA)

bed_list12 = list(bed22,bed11)
### 默认情况下点的颜色
circos.genomicTrack(bed_list12, 
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicPoints(region, value, pch = 16, cex = 0.5, col = i, ...)
                    })
### 这个图会不断的叠加，所以如果重新画记得归零
###自己设置颜色，大于阈值（例如0）即上调为红色散点，下调为绿色散点
text(0, 0, "DARs-P0.01", cex = 1)
text(0, 0.3, "up-9679", cex = 1)
text(0, 0.2, "down-15080", cex = 1)
circos.genomicTrack(bed_list12, 
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicPoints(region, value, pch = 16, cex = 0.5, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), 
                                           circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040"),...)
                    })

#清除
circos.clear()

########################################################################
####  整合
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq/annotation")
library(data.table)
homer_anno_down<-fread("sensitivity-resistance-S8495-median-AUC-P0.01-logFC0-DESeq2-down-annotation.txt",header=T,data.table=F)
length(unique(homer_anno_down$`Gene Name`))
colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]
length(unique(Anno_gene_100Kb_down$`Gene Name`))

homer_anno_up<-fread("sensitivity-resistance-S8495-median-AUC-P0.01-logFC0-DESeq2-up-annotation.txt",header=T,data.table=F)
length(unique(homer_anno_up$`Gene Name`))
colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
length(unique(Anno_gene_100Kb_up$`Gene Name`))

setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/1_RNAseq")
up_gene<-read.table("sensitivity-resistance-S8495-DESeq2_median_AUC_RNAseq_P0.05_FC0_up_DEGs.txt",header=T,sep="\t")
down_gene<-read.table("sensitivity-resistance-S8495-DESeq2_median_AUC_RNAseq_P0.05_FC0_down_DEGs.txt",header=T,sep="\t")
all<-rbind(up_gene,down_gene)
DEGs<-as.character(all[,1])
down_DEGs_DApeaks<-intersect(Anno_gene_100Kb_down$`Gene Name`,DEGs)
up_DEGs_DApeaks<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)
df<-all[match(c(down_DEGs_DApeaks,up_DEGs_DApeaks),DEGs),]

downpeaks<-Anno_gene_100Kb_down[(Anno_gene_100Kb_down$`Gene Name`%in%down_DEGs_DApeaks),1]
uppeaks<-Anno_gene_100Kb_up[(Anno_gene_100Kb_up$`Gene Name`%in%up_DEGs_DApeaks),1]
allpeaks11<-c(downpeaks,uppeaks)
df11<-read.table("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495_median/2_ATACseq/sensitivity-resistance-all-DESeq2_S8495_median_AUC_ATAC.txt",sep = '\t',header= T)
df1<-df11[match(allpeaks11,df11[,1]),]

library(circlize)
library(rtracklayer)
gtf_data = import('~/xjj/CUP/gtf/Homo_sapiens.GRCh37.75.gtf') #gtf的路径
gtf_data = as.data.frame(gtf_data)

### ATAC
#library("data.table")
#setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_GEM/2_ATACseq/circos")
#df1<-read.table("DApeakall.txt",sep = '\t',header= T)
resSig1<-df1[df1$pvalue<0.01,]
up_peak<-resSig1[which(resSig1$log2FoldChange>0),]
down_peak<-resSig1[which(resSig1$log2FoldChange<0),]
library(dplyr)
peak_up<-as.data.frame(t(apply(as.data.frame(up_peak$peak_id),1,function(x) unlist(strsplit(as.character(x), "_")))))
peak_down<-as.data.frame(t(apply(as.data.frame(down_peak$peak_id),1,function(x) unlist(strsplit(as.character(x), "_")))))

up_logFC = up_peak$log2FoldChange
down_logFC =down_peak$log2FoldChange
bed11=as.data.frame(cbind(peak_up,up_logFC))
colnames(bed11)<-c("chr","start","end","value1")
bed11[,2]<-as.integer(as.character(bed11[,2]))
bed11[,3]<-as.integer(as.character(bed11[,3]))

bed22=as.data.frame(cbind(peak_down,down_logFC))
colnames(bed22)<-c("chr","start","end","value1")
bed22[,2]<-as.integer(as.character(bed22[,2]))
bed22[,3]<-as.integer(as.character(bed22[,3]))
#黑色下调，红上调
set.seed(60)
circos.initializeWithIdeogram(species = "hg19",plotType = c("axis"))
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(1))
  circos.text(mean(xlim), mean(ylim), chr, cex = 0.7, col = "white",
              facing = "inside", niceFacing = TRUE)
}, track.height = 0.15, bg.border = NA)

bed_list12 = list(bed22,bed11)
#text(0, 0, "DEG-P0.05-DARs-P0.01", cex = 1)
#text(0, 0.3, "up-9679", cex = 1)
#text(0, 0.2, "down-15080", cex = 1)
circos.genomicTrack(bed_list12, 
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicPoints(region, value, pch = 16, cex = 0.5, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), 
                                           circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040"),...)
                    })

#setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_GEM/2_ATACseq/circos")
#df<-read.table("all_DEGs.txt",sep = '\t',header= T)
resSig<-df[df$pvalue<0.05,]
up_gene<-resSig[which(resSig$log2FoldChange>0),]
down_gene<-resSig[which(resSig$log2FoldChange<0),]

gene_up=up_gene$gene_id
gene_down=down_gene$gene_id
up_logFC = up_gene$log2FoldChange
down_logFC =down_gene$log2FoldChange

gene_up_mete=gtf_data[na.omit(match(gene_up,gtf_data$gene_name)),]
bed1=cbind(gene_up_mete[,1:3],up_logFC[match(gene_up_mete$gene_name,gene_up)])
bed1[,1]<-paste(rep("chr",nrow(bed1)),bed1[,1],sep="")
colnames(bed1) = c("chr","start","end","value1")

gene_down_mete=gtf_data[na.omit(match(gene_down,gtf_data$gene_name)),]
bed2=cbind(gene_down_mete[,1:3],down_logFC[match(gene_down_mete$gene_name,gene_down)])
bed2[,1]<-paste(rep("chr",nrow(bed2)),bed2[,1],sep="")
colnames(bed2) = c("chr","start","end","value1")

####3.初始化基因组圈图####
bed_list = list(bed2,bed1)
circos.genomicTrack(bed_list, 
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicPoints(region, value, pch = 16, cex = 0.5, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), 
                                           circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040"),...)
                    })

text(0, 1, "Overlap between DEGs and the nearest genes of DApeaks", cex = 1)


