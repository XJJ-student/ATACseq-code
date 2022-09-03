####1 样本信息################################配对
GEM_sensitive<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
                 "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
                 "DAC-19","DAC-37","DAC-39","DAC-38")
GEM_res_S8495_sen<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
GEM_res_S8495_res<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")


# I:GEM_sen II:S8495_sen III:S8495_res
#### I II+III
GEM_sensitive<-GEM_sensitive
GEM_resistance<-c(GEM_res_S8495_sen,GEM_res_S8495_res)

#### II I+III
GEM_res_S8495_sen<-GEM_res_S8495_sen
other2<-c(GEM_sensitive,GEM_res_S8495_res)

#### III I+II
GEM_res_S8495_res<-GEM_res_S8495_res
other3<-c(GEM_sensitive,GEM_res_S8495_sen)

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

resistance<-other3
sensitivity<-GEM_res_S8495_res

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
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/1_RNAseq/GEM_res_S8495_res")
write.table(res,"GEM_res_S8495_res-other3-all-DESeq2_AUC_RNAseq_DEGs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
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
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/1_RNAseq/GEM_res_S8495_res")
write.table(up_gene,"GEM_res_S8495_res_other3_RNAseq_P0.05_FC0_up_DEGs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(down_gene,"GEM_res_S8495_res_other3_RNAseq_P0.05_FC0_down_DEGs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

##################################### 转录组区分亚组  转录组的特征不是很完美的将四个亚组分开#################
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_GEM/1_RNAseq")
GEM_sen<-read.table("sensitivity-resistance-GEM-DESeq2_AUC_RNAseq_P0.05_FC0_up_DEGs.txt",sep="\t",header=T)
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/1_RNAseq/GEM_res_S8495_res")
GEM_res_S8495_res<-read.table("GEM_res_S8495_res_other3_RNAseq_P0.05_FC0_up_DEGs.txt",sep="\t",header=T)
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/1_RNAseq/GEM_res_S8495_sen")
GEM_res_S8495_sen<-read.table("GEM_res_S8495_sen_other2_RNAseq_P0.05_FC0_up_DEGs.txt",sep="\t",header=T)

one1<-intersect(as.character(GEM_sen[,1]),as.character(GEM_res_S8495_res[,1]))#0
one2<-intersect(as.character(GEM_sen[,1]),as.character(GEM_res_S8495_sen[,1]))#1 "CRISP3"
one3<-intersect(as.character(GEM_res_S8495_res[,1]),as.character(GEM_res_S8495_sen[,1]))#1 "UNC93B1"

##用这些特征画热图
GEM_sensitive1<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
                 "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
                 "DAC-19","DAC-37","DAC-39","DAC-38")
GEM_res_S8495_sen1<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
GEM_res_S8495_res1<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")

sample_order<-c(GEM_sensitive1,GEM_res_S8495_sen1,GEM_res_S8495_res1)
setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
library("data.table")
exp<-fread("star_rsem.GeneSymbol.TPM.xls",header=T,data.table=F)
peaks<-c(as.character(GEM_sen[,1]),as.character(GEM_res_S8495_sen[,1]),as.character(GEM_res_S8495_res[,1]))
peak_RPKM1<-exp[match(peaks,exp[,1]),]
peaks1<-matrix(c(rep("Class1",nrow(GEM_sen)),rep("Class2",nrow(GEM_res_S8495_sen)),rep("Class3",nrow(GEM_res_S8495_res))),ncol=1)
peak_RPKM_order<-peak_RPKM1[,match(sample_order,colnames(peak_RPKM1))]

row_cut=c(nrow(GEM_sen),nrow(GEM_sen)+nrow(GEM_res_S8495_sen))
col_cut=c(length(GEM_sensitive1),length(GEM_sensitive1)+length(GEM_res_S8495_sen1))
library(pheatmap)
data0=cbind(peaks1,peak_RPKM_order)
anno1=c(rep("Class1",length(GEM_sensitive1)),rep("Class2",length(GEM_res_S8495_sen1)),rep("Class3",length(GEM_res_S8495_res1)))
anno1=data.frame(anno1)
rownames(anno1)=as.character(t(colnames(data0))[2:38])
ann_colors = list(
  anno1 = c(Class1 = "#4dbbd5", Class2 = "#00a087",Class3 = "#f39b7f"), anno=c("Class1"= "#4dbbd5","Class2"= "#00a087","Class3"= "#f39b7f"))
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
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/1_RNAseq")
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

##################################################################################
################  DEGs与组蛋白修饰酶的交叠
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDAC_frontiers/5_histone-modifying enzymes")
Acetylase<-read.table("histone-modifying-enzymes.txt",header=T,sep="\t")
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/1_RNAseq/GEM_res_S8495_res")
df<-read.table("GEM_res_S8495_res-other3-all-DESeq2_AUC_RNAseq_DEGs.txt",sep = '\t',header= T,row.names = 1)

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
intersect(Acetylase_name,DEGs_name)
input  <-list(Acetylase_name,DEGs_name)
venn(input,showSetLogicLabel=TRUE)
tmp <- venn(input)
int<-attr(tmp, "intersections")
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/1_RNAseq/GEM_res_S8495_res/Acetylase")
write.table(Acetylase_info,"Acetylase_info_P0.05.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
library(VennDiagram)
venn.diagram(x=list(Acetylase=Acetylase_name,DEGs=DEGs_name),cex = 1,margin = 0.1, "Acetylase_info_P005.png",fill=c("red","blue"))

setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
library("data.table")
exp<-fread("star_rsem.GeneSymbol.TPM.xls",header=T,data.table=F)
sensitivity<-GEM_res_S8495_res
resistance<-other3

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
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/1_RNAseq/GEM_res_S8495_res/Acetylase")
write.table(Acetylase_info_value,"GEM_res_S8495_res_Acetylase_info_P0.05_value.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

library(ggplot2)
df2<-cbind(rbind(as.matrix(Acetylase_info_value[,1]),as.matrix(Acetylase_info_value[,1])),
           rbind(as.matrix(sensitivity_exp11),as.matrix(resistance_exp11)),
           rbind(as.matrix(rep("sensitive",length(Acetylase_info_value[,1]))),as.matrix(rep("resistant",length(Acetylase_info_value[,1])))),
           rbind(as.matrix(Acetylase_info_value$sensitive+sensitivity_sd),as.matrix(Acetylase_info_value$resistant+resistance_sd)),
           rbind(as.matrix(Acetylase_info_value$sensitive-sensitivity_sd),as.matrix(Acetylase_info_value$resistant-resistance_sd)),
           rbind(as.matrix(Acetylase_info_value[,8]),as.matrix(Acetylase_info_value[,8])))
colnames(df2)<-c("gene","value","type","valuesd1","valuesd2","pvalue")
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/1_RNAseq/GEM_res_S8495_res/Acetylase")
write.table(df2,"GEM_res_S8495_res_Acetylase_info_P0.05_value_plot.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

df3<-read.table("GEM_res_S8495_res_Acetylase_info_P0.05_value_plot.txt",sep = '\t',header= T)
temp=df3[df3[,3]=="sensitive",1]
df3[,1]=factor(df3[,1],levels=temp)
df3[which(df3$pvalue>0.01),'pva']<-"*"
df3[which(df3$pvalue<0.01 & df3$pvalue>0.001),'pva']<-"**"
df3[which(df3$pvalue<0.001),'pva']<-"***"

#df3[(length(Acetylase_info_value[,1])+1):34,7]=NA
df3$valuesd1<-df3$value+1
p3 <- ggplot(df3, aes(x = gene, y = value, fill = type)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge(),width=0.4)+ 
  geom_errorbar(aes(ymin = valuesd1, ymax = value), width = 0.2, position = position_dodge(0.4))+
  labs(x = "Differently expression of histone-modifying enzymes", y = "mRNA average expression levels (TPM)")+
  geom_text(aes(y=30,label = pva),size = 3)+
  theme(axis.title.y = element_text(color = 'black',size = 12),
        axis.title.x = element_text(color = 'black',size = 12,vjust = -1.2),
        axis.text.y = element_text(color = 'black',size = 10),
        axis.text.x = element_text(color = 'black',size = 12,angle = 30,vjust = 0.5))
plot(p3)
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/1_RNAseq/GEM_res_S8495_res/Acetylase")
ggsave("GEM_res_S8495_res_Acetylase_info_P0.05_value_plot.pdf",p3,width = 10, height = 6)


#################################################################################
################  clusterProfiler
rm(list=ls())
library(org.Hs.eg.db)
library(clusterProfiler)
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/1_RNAseq/GEM_res_S8495_sen")
up_gene<-read.table("GEM_res_S8495_sen_other2_RNAseq_P0.05_FC0_down_DEGs.txt",header=T,sep="\t")
library(stringr)
gene=bitr(up_gene[,1],fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
kegg<-enrichKEGG(gene[,2],organism = 'hsa',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'BH',
                 minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.05,use_internal_data = FALSE)
#keu<-barplot(kegg,showCategory=15,drop=T,x = "GeneRatio",color = "pvalue")
keu<-barplot(kegg,showCategory=15,drop=T,x = "GeneRatio",color = "p.adjust")

ee<-kegg@result
ee1<-ee[which(ee$pvalue<0.05),]
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
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/1_RNAseq/GEM_res_S8495_sen/pathway")
write.table(ee1,"GEM_res_S8495_sen_kegg_DEGP0.05_pathway_P0.05_up_drug_target.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(pathway_gene2,"GEM_res_S8495_sen_kegg_DEGP0.05_pathway_P0.05_up_drug_target_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
ggsave("GEM_res_S8495_sen_kegg_DEGP0.05_pathway_P0.05_up_drug_target_gene.pdf",keu,width = 8, height = 6)


setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/1_RNAseq/GEM_res_S8495_sen/pathway")
write.table(ee1,"GEM_res_S8495_sen_kegg_DEGP0.05_pathway_P0.05_down_drug_target.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(pathway_gene2,"GEM_res_S8495_sen_kegg_DEGP0.05_pathway_P0.05_down_drug_target_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
ggsave("GEM_res_S8495_sen_kegg_DEGP0.05_pathway_P0.05_down_drug_target_gene.pdf",keu,width = 8, height = 6)

go<-enrichGO(gene[,2],OrgDb = org.Hs.eg.db,ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.05,keyType = 'ENTREZID')
ked<-dotplot(go,showCategory=10)
a<-go@result
go_BP<-a[which(a$p.adjust<0.05),]

enrich_genego<-go_BP$geneID
pathway_genego<-unique(unlist(strsplit(enrich_genego,split="/")))
pathway_genego2=bitr(pathway_genego,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/1_RNAseq/GEM_res_S8495_sen/pathway")
write.table(pathway_genego2,"GEM_res_S8495_sen_GO_DEG_pathway_FDR0.05_up_drug_target_gene.txt",sep = '\t',col.names = F,row.names = F,quote = FALSE)##数据输出
write.table(go_BP,"GEM_res_S8495_sen_GO_DEG_pathway_FDR0.05_up_drug_target.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
ggsave("GEM_res_S8495_sen_GO_DEG_pathway_FDR0.05_up_drug_target_gene.pdf",ked,width = 8, height = 6)

setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/1_RNAseq/GEM_res_S8495_sen/pathway")
write.table(pathway_genego2,"GEM_res_S8495_sen_GO_DEG_pathway_FDR0.05_down_drug_target_gene.txt",sep = '\t',col.names = F,row.names = F,quote = FALSE)##数据输出
write.table(go_BP,"GEM_res_S8495_sen_GO_DEG_pathway_FDR0.05_down_drug_target.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
ggsave("GEM_res_S8495_sen_GO_DEG_pathway_FDR0.05_down_drug_target_gene.pdf",ked,width = 8, height = 6)

#######################    ATACseq   ################################################
#############  DApeaks
GEM_sensitive<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
                 "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
                 "DAC-19","DAC-37","DAC-39","DAC-38")
GEM_res_S8495_sen<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
GEM_res_S8495_res<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")


# I:GEM_sen II:S8495_sen III:S8495_res
#### I II+III
GEM_sensitive<-GEM_sensitive
GEM_resistance<-c(GEM_res_S8495_sen,GEM_res_S8495_res)

#### II I+III
GEM_res_S8495_sen<-GEM_res_S8495_sen
other2<-c(GEM_sensitive,GEM_res_S8495_res)

#### III I+II
GEM_res_S8495_res<-GEM_res_S8495_res
other3<-c(GEM_sensitive,GEM_res_S8495_sen)


rm(list=ls())
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_count_name<-read.table("peak_count.txt",sep="\t",header=T) #没有全0的

resistance<-other2
sensitivity<-GEM_res_S8495_sen

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
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq/GEM_res_S8495_res")
write.table(res,"GEM_res_S8495_res-other3-all-DESeq2_chemotherapy_AUC_ATAC.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

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
out_file<-"GEM_res_S8495_sen-other2-AUC-P0.05-logFC1"
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
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq")
write.table(peaks_bed_up,file="all_peaks.bed",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')


##################ATACseq数据 四类中共有的peaks和特有的peaks
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_GEM/2_ATACseq")
GEM_sen<-read.table("sensitivity-resistance-GEM-AUC-P0.01-logFC0-DESeq2-up.txt",sep="\t",header=T)
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq/GEM_res_S8495_res")
GEM_res_S8495_res<-read.table("GEM_res_S8495_res-other3-AUC-P0.01-logFC0-DESeq2-up.txt",sep="\t",header=T)
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq/GEM_res_S8495_sen")
GEM_res_S8495_sen<-read.table("GEM_res_S8495_sen-other2-AUC-P0.01-logFC0-DESeq2-up.txt",sep="\t",header=T)

one1<-intersect(as.character(GEM_sen[,1]),as.character(GEM_res_S8495_res[,1]))#0
one2<-intersect(as.character(GEM_sen[,1]),as.character(GEM_res_S8495_sen[,1]))#0
one3<-intersect(as.character(GEM_res_S8495_res[,1]),as.character(GEM_res_S8495_sen[,1]))#1 "chr10_23105182_23105682"

##用这些特征画热图
GEM_sensitive1<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
                  "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
                  "DAC-19","DAC-37","DAC-39","DAC-38")
GEM_res_S8495_sen1<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
GEM_res_S8495_res1<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")

sample_order<-c(GEM_sensitive1,GEM_res_S8495_sen1,GEM_res_S8495_res1)

setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_RPKM<-read.table("peak_RPKM.txt",sep="\t",header=T)
colnames(peak_RPKM)<-gsub("\\.","-",colnames(peak_RPKM))
peaks<-c(as.character(GEM_sen[,1]),as.character(GEM_res_S8495_sen[,1]),as.character(GEM_res_S8495_res[,1]))
peak_RPKM1<-peak_RPKM[match(peaks,peak_RPKM[,1]),]
peaks1<-matrix(c(rep("Class1",nrow(GEM_sen)),rep("Class2",nrow(GEM_res_S8495_sen)),rep("Class3",nrow(GEM_res_S8495_res))),ncol=1)
peak_RPKM_order<-peak_RPKM1[,match(sample_order,colnames(peak_RPKM1))]


row_cut=c(nrow(GEM_sen),nrow(GEM_sen)+nrow(GEM_res_S8495_sen))
col_cut=c(length(GEM_sensitive1),length(GEM_sensitive1)+length(GEM_res_S8495_sen1))
library(pheatmap)
data0=cbind(peaks1,peak_RPKM_order)

anno1=c(rep("Class1",length(GEM_sensitive1)),rep("Class2",length(GEM_res_S8495_sen1)),rep("Class3",length(GEM_res_S8495_res1)))
anno1=data.frame(anno1)

rownames(anno1)=as.character(t(colnames(data0))[2:38])
ann_colors = list(
  anno1 = c(Class1 = "#4dbbd5", Class2 = "#00a087",Class3 = "#f39b7f"), anno=c("Class1"= "#4dbbd5","Class2"= "#00a087","Class3"= "#f39b7f"))
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
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq")
pheatmap(data0[2:38],scale="row",color = mycolor,fontsize=7,fontsize_row = 70*H/dim(data0)[1], fontsize_col = 7,cluster_cols=FALSE, cluster_rows=FALSE,annotation_col=anno1,gaps_col = col_cut,gaps_row = row_cut, annotation_row=anno,annotation_colors = ann_colors)
dev.off()

three_categories_ATACseq_pheatmap


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
folder<-c("GEM_res_S8495_sen","GEM_res_S8495_res")

TSS100KB_anno_result<-data.frame(matrix(0,length(folder),7))
foreach(k=1:length(folder))%do%{
  file_name=as.character(folder[k])
  cat(file_name,"\n")
  TSS100KB_anno_result[k,1]<-file_name
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq/",file_name,"/annotation",sep=""))
  homer_anno_down<-fread(paste(file_name,"-other",k+1,"-AUC-P0.01-logFC0-DESeq2-down-annotation.txt",sep=""),header=T,data.table=F)
  colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
  Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]
  cat(length(unique(Anno_gene_100Kb_down$`Gene Name`)),"\n")
  
  homer_anno_up<-fread(paste(file_name,"-other",k+1,"-AUC-P0.01-logFC0-DESeq2-up-annotation.txt",sep=""),header=T,data.table=F)
  colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
  Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
  cat(length(unique(Anno_gene_100Kb_up$`Gene Name`)),"\n")
  
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/1_RNAseq/",file_name,sep=""))
  up_gene<-read.table(paste(file_name,"_other",k+1,"_RNAseq_P0.05_FC0_up_DEGs.txt",sep=""),header=T,sep="\t")
  down_gene<-read.table(paste(file_name,"_other",k+1,"_RNAseq_P0.05_FC0_down_DEGs.txt",sep=""),header=T,sep="\t")
  cat(nrow(up_gene),"\n")
  cat(nrow(down_gene),"\n")
  
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
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/3_intersect/venn/",file_name,sep=""))
  venn.diagram(x=list(downpeaks=unique(down_peaks_gene),uppeaks=unique(up_peaks_gene),RNAseq=unique(DEGs)), paste("DEGs0.05_int_DApeaks0.01_",file_name,".png",sep=""),fill=c("red","green","blue"),margin = 0.1)
  
  up1<-int[["B:C"]]
  up2<-int[["A:B:C"]]
  up3<-c(up1,up2)
  int_up<-as.data.frame(intersect(up3,as.character(up_gene$gene_id)))
  colnames(int_up)<-"int_up"
  cat(paste("Up peak int up DEGs is",nrow(int_up),sep=" "),"\n")
  int_up_peak<-as.data.frame(homer_anno_up[homer_anno_up$`Gene Name`%in%int_up[,1],1])
  colnames(int_up_peak)<-"int_up_peak"
  
  down1<-int[["A:C"]]
  down2<-int[["A:B:C"]]
  down3<-c(down1,down2)
  int_down<-as.data.frame(intersect(down3,as.character(down_gene$gene_id)))
  colnames(int_down)<-"int_down"
  cat(paste("Down peak int down DEGs is",nrow(int_down),sep=" "),"\n")
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
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/3_intersect/venn/",file_name,sep=""))
  write.table(int_up,"int_up_peak0.01_all_DEG0.05_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  write.table(int_down,"int_down_peak0.01_all_DEG0.05_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  write.table(int_up_peak,"int_up_peak0.01_all_DEG0.05_peak.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  write.table(int_down_peak,"int_down_peak0.01_all_DEG0.05_peak.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  write.table(up_DEGs_up_DApeaks,"int_up_peak0.01_up_DEG0.05_peak.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  write.table(down_DEGs_down_DApeaks,"int_down_peak0.01_down_DEG0.05_peak.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  
}
colnames(TSS100KB_anno_result)<-c("group","Up peak int DEGs","Up peak int up-DEGs","Hypergeometric-up","Down peak int DEGs","Down peak int down-DEGs","Hypergeometric-down")
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/3_intersect/venn")
write.table(TSS100KB_anno_result,"TSS100KB_anno_result_result_2group.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

##################################################################################
##############  DEpeaks在不同样本变化的倍数与DEGs在不同样本变化的倍数之间的相关性(一个基因对应多个peaks)
rm(list=ls())
folder<-c("GEM_res_S8495_sen","GEM_res_S8495_res")
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
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq/",file_name,"/annotation",sep=""))
  homer_anno_down<-fread(paste(file_name,"-other",k+1,"-AUC-P0.01-logFC0-DESeq2-down-annotation.txt",sep=""),header=T,data.table=F)
  colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
  Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]
  
  homer_anno_up<-fread(paste(file_name,"-other",k+1,"-AUC-P0.01-logFC0-DESeq2-up-annotation.txt",sep=""),header=T,data.table=F)
  colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
  Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
  
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/1_RNAseq/",file_name,sep=""))
  up_gene<-read.table(paste(file_name,"_other",k+1,"_RNAseq_P0.05_FC0_up_DEGs.txt",sep=""),header=T,sep="\t")
  down_gene<-read.table(paste(file_name,"_other",k+1,"_RNAseq_P0.05_FC0_down_DEGs.txt",sep=""),header=T,sep="\t")
  all<-rbind(up_gene,down_gene)
  DEGs<-as.character(all[,1])
  down_DEGs_DApeaks<-intersect(Anno_gene_100Kb_down$`Gene Name`,DEGs)
  up_DEGs_DApeaks<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)
  
  int_all_DEGs<-unique(c(down_DEGs_DApeaks,up_DEGs_DApeaks))
  all_peaks_anno<-rbind(homer_anno_down,homer_anno_up)
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq/",file_name,sep=""))
  DApeak<-fread(paste(file_name,"-other",k+1,"-all-DESeq2_chemotherapy_AUC_ATAC.txt",sep=""),header=T,data.table=F)
  
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
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/3_intersect/correlation/",file_name,sep=""))
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
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/3_intersect/correlation/",file_name,sep=""))
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
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/3_intersect/correlation/",file_name,sep=""))
  ggsave(paste(file_name,"correlation-DApeaksP0.01-DEGsP0.05-pearson-comprehensive.pdf",sep=""),pp,width = 10, height = 10)
  #plot(p)
}
#pearson", "kendall", "spearman
colnames(DEGs_corre_DApeaks)<-c("group","spearman-estimate","spearman-pvalue","pearson-estimate","pearson-pvalue")
setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/3_intersect/correlation/",sep=""))
write.table(DEGs_corre_DApeaks,"DEGs_FC_corre_DApeaks_FC_result_2group_comprehensive.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

###################################################################################
##############  DEpeaks_FC与DEGs_FC之间的相关性(一个基因对应一个peaks)
rm(list=ls())
folder<-c("GEM_res_S8495_sen","GEM_res_S8495_res")
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
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq/",file_name,"/annotation",sep=""))
  homer_anno_down<-fread(paste(file_name,"-other",k+1,"-AUC-P0.01-logFC0-DESeq2-down-annotation.txt",sep=""),header=T,data.table=F)
  colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
  Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]
  
  homer_anno_up<-fread(paste(file_name,"-other",k+1,"-AUC-P0.01-logFC0-DESeq2-up-annotation.txt",sep=""),header=T,data.table=F)
  colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
  Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
  
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/1_RNAseq/",file_name,sep=""))
  up_gene<-read.table(paste(file_name,"_other",k+1,"_RNAseq_P0.05_FC0_up_DEGs.txt",sep=""),header=T,sep="\t")
  down_gene<-read.table(paste(file_name,"_other",k+1,"_RNAseq_P0.05_FC0_down_DEGs.txt",sep=""),header=T,sep="\t")
  all<-rbind(up_gene,down_gene)
  DEGs<-as.character(all[,1])
  down_DEGs_DApeaks<-intersect(Anno_gene_100Kb_down$`Gene Name`,DEGs)
  up_DEGs_DApeaks<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)
  
  int_all_DEGs<-unique(c(down_DEGs_DApeaks,up_DEGs_DApeaks))
  all_peaks_anno<-rbind(homer_anno_down,homer_anno_up)
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq/",file_name,sep=""))
  DApeak<-fread(paste(file_name,"-other",k+1,"-all-DESeq2_chemotherapy_AUC_ATAC.txt",sep=""),header=T,data.table=F)
  
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
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/3_intersect/correlation/",file_name,sep=""))
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
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/3_intersect/correlation/",file_name,sep=""))
  ggsave(paste(file_name,"correlation-DApeaksP0.01-DEGsP0.05-spearman-onlyone.pdf",sep=""),ps,width = 10, height = 10)
}
#pearson", "kendall", "spearman
colnames(DEGs_corre_DApeaks)<-c("group","pearson-estimate","pearson-pvalue","spearman-estimate","spearman-pvalue")
setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/3_intersect/correlation",sep=""))
write.table(DEGs_corre_DApeaks,"DEGs_FC_corre_DApeaks_FC_result_2group_onlyone.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

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
folder<-c("GEM_res_S8495_sen","GEM_res_S8495_res")
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(stringr)
library(ggplot2)
freq_name<-"FDR0.05"

peaks_DEGs_KEGG<-data.frame(matrix(0,length(folder),5))
foreach(k=1:length(folder))%do%{
  file_name=as.character(folder[k])
  cat(file_name,"\n")
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq/",file_name,"/annotation",sep=""))
  homer_anno_down<-fread(paste(file_name,"-other",k+1,"-AUC-P0.01-logFC0-DESeq2-down-annotation.txt",sep=""),header=T,data.table=F)
  colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
  Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/1_RNAseq/",file_name,sep=""))
  down_gene<-read.table(paste(file_name,"_other",k+1,"_RNAseq_P0.05_FC0_down_DEGs.txt",sep=""),header=T,sep="\t")
  
  DEGs<-as.character(down_gene[,1])
  down_DEGs_DApeaks<-intersect(Anno_gene_100Kb_down$`Gene Name`,DEGs)
  cat(paste("Down int gene is ",length(down_DEGs_DApeaks),sep=""),"\n")
  gened=bitr(down_DEGs_DApeaks,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
  keggd<-enrichKEGG(gened[,2],organism = 'hsa',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'BH',
                    minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.05,use_internal_data = FALSE)
  eed<-keggd@result
  ee1d<-eed[which(eed$p.adjust<0.05),]
  if(nrow(ee1d)>0){
    for(y in 1:nrow(ee1d)){
      ee1d[y,8]
      b1<-matrix(unlist(strsplit(ee1d[y,8],split="/")),ncol=1)
      gene6=bitr(b1,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
      ee1d[y,10]<-paste0(gene6[,2],collapse ="/")
    }
  }
  if(nrow(ee1d)==0){
    ee1d=ee1d
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
  
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/3_intersect/pathway/",file_name,sep=""))
  write.table(ee1d,paste("kegg_pathway_DEGs_DApeaks_",file_name,freq_name,"_down.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  write.table(pathway_gened,paste("kegg_pathway_DEGs_DApeaks_",file_name,freq_name,"_down_gene.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  ggsave(paste("KEGG pathway of",file_name,freq_name,"-down-P0.05.pdf",sep=""),ked,width = 10, height = 5)
  
  
  
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq/",file_name,"/annotation",sep=""))
  homer_anno_up<-fread(paste(file_name,"-other",k+1,"-AUC-P0.01-logFC0-DESeq2-up-annotation.txt",sep=""),header=T,data.table=F)
  colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
  Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/1_RNAseq/",file_name,sep=""))
  up_gene<-read.table(paste(file_name,"_other",k+1,"_RNAseq_P0.05_FC0_up_DEGs.txt",sep=""),header=T,sep="\t")
  DEGs<-as.character(up_gene[,1])
  up_DEGs_DApeaks<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)
  cat(paste("Up int gene is ",length(up_DEGs_DApeaks),sep=""),"\n")
  geneu=bitr(up_DEGs_DApeaks,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
  keggu<-enrichKEGG(geneu[,2],organism = 'hsa',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'BH',
                    minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.05,use_internal_data = FALSE)
  eeu<-keggu@result
  ee1u<-eeu[which(eeu$p.adjust<0.05),]
  if(nrow(ee1u)>0){
    for(y in 1:nrow(ee1u)){
      ee1u[y,8]
      b7<-matrix(unlist(strsplit(ee1u[y,8],split="/")),ncol=1)
      gene7=bitr(b7,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
      ee1u[y,10]<-paste0(gene7[,2],collapse ="/")
    }
  }
  if(nrow(ee1u)==0){
    ee1u=ee1u
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
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/3_intersect/pathway/",file_name,sep=""))
  write.table(ee1d,paste("kegg_pathway_DEGs_DApeaks_",file_name,freq_name,"_up.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  write.table(pathway_gened,paste("kegg_pathway_DEGs_DApeaks_",file_name,freq_name,"_up_gene.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  ggsave(paste("KEGG pathway of",file_name,freq_name,"-up-P0.05.pdf",sep=""),ked,width = 10, height = 5)
}

colnames(peaks_DEGs_KEGG)<-c("drug","KEGG-down-P0.05-pathway","KEGG-down-P0.05-pathway-gene","KEGG-up-P0.05-pathway","KEGG-up-P0.05-pathway-gene")
setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/3_intersect/pathway",sep=""))
write.table(peaks_DEGs_KEGG,"peaks_DEGs_KEGG_P0.01_result_2group.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出


#15 #############  ridge regression  #####因为FIMO能知道peaks与TFs的关系
rm(list=ls())
library(dplyr)
library(data.table)
library(glmnet)
library(foreach)
library(ggplot2)
folder<-c("GEM_res_S8495_sen","GEM_res_S8495_res")

result<-matrix(0,2,6)
for(k in 1:2){
  file_name=as.character(folder[k])
  cat(file_name,"\n")
  result[k,1]<-file_name
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq/",file_name,"/FIMO/P0.01_FC0_down",sep=""))
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
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq/",file_name,sep=""))
  DApeak<-fread(paste(file_name,"-other",k+1,"-all-DESeq2_chemotherapy_AUC_ATAC.txt",sep=""),header=T,data.table=F)
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/1_RNAseq/",file_name,sep=""))
  DE_gene<-fread(paste(file_name,"-other",k+1,"-all-DESeq2_AUC_RNAseq_DEGs.txt",sep=""),header=T,data.table=F)
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
setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq/",file_name,"/FIMO/P0.01_FC0_down",sep=""))
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
write.table(output_coef,file="ridge_regression_down.txt",sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')
res_positive1<-output_coef[which(output_coef$coefficient>0),]

resSig<-DE_gene[which(DE_gene$pvalue<0.05),1]
down_int<-intersect(as.character(res_positive1[,1]),resSig)
result[k,2]<-length(down_int)
result[k,3]<-paste0(down_int,collapse =",")
#######
setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq/",file_name,"/FIMO/P0.01_FC0_up",sep=""))
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
setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq/",file_name,"/FIMO/P0.01_FC0_up",sep=""))
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
write.table(output_coef_u,file="ridge_regression_up.txt",sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')
res_positive1_u<-output_coef_u[which(output_coef_u$coefficient>0),]
up_int<-intersect(res_positive1_u,resSig)
result[k,4]<-length(up_int)
result[k,5]<-paste0(up_int,collapse =",")
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
setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq/",file_name,sep=""))
ggsave(paste("ridge_regression_","DApeakP0.01","TF_coefficient.pdf",sep=""),pu,width = 8, height = 6)
write.table(up_TFs1,"FIMO_result_up_TFs_motifFDR0.05_DApeakP0.01.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(down_TFs1,"FIMO_result_down_TFs_motifFDR0.05_DApeakP0.01.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
result[k,6]<-nrow(up_TFs1)
result[k,7]<-nrow(down_TFs1)

}
colnames(result)<-c("Type","DEGs int downTFs","DEGs int downTFs","DEGs int upTFs","DEGs int upTFs","upTFs","downTFs")


##################  ridge regresion 找到的TF之间的相关性
##PPI网络直接在网站上面进行 https://www.string-db.org/
##PPI网络里面可以直接出来TFs富集的各个区域、位置和通路
###############  GSEA  cancer module
rm(list=ls())
library(clusterProfiler)
library(GSEABase)
library(org.Hs.eg.db)
library(stringr)
####  DEGs
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_GEM/2_ATACseq/FIMO")
up_TFs1<-read.table("FIMO_result_up_TFs_motifFDR0.05_DApeakP0.01.txt",sep = '\t',header=TRUE)
down_TFs1<-read.table("FIMO_result_down_TFs_motifFDR0.05_DApeakP0.01.txt",sep = '\t',header=TRUE)

DEGs<-down_TFs1
colnames(DEGs)<-c("symbol","coefficient")
geneList<-DEGs$coefficient #第二列可以是folodchange，也可以是logFC
names(geneList)=DEGs$symbol #使用转换好的ID
geneList=sort(geneList,decreasing =T) #从高到低排序
gene1<-names(geneList)
setwd("~/xjj/GSEA")
c5<-read.delim("c1.all.v7.4.symbols.txt",header=F,sep='\t')
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
egmt3 <- GSEA(geneList, TERM2GENE=result)
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





rm(list=ls())
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_GEM/2_ATACseq/FIMO")
up_TFs1<-read.table("FIMO_result_up_TFs_motifFDR0.05_DApeakP0.01.txt",sep = '\t',header=TRUE)
down_TFs1<-read.table("FIMO_result_down_TFs_motifFDR0.05_DApeakP0.01.txt",sep = '\t',header=TRUE)
DETFs<-c(as.character(up_TFs1[,1]),as.character(down_TFs1[,1]))

library(org.Hs.eg.db)
library(clusterProfiler)
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/1_RNAseq/GEM_res_S8495_res")
up_gene<-read.table("GEM_res_S8495_res_other3_RNAseq_P0.05_FC0_up_DEGs.txt",header=T,sep="\t")
library(stringr)
gene=bitr(up_gene[,1],fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
kegg<-enrichKEGG(gene[,2],organism = 'hsa',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'none',
                 minGSSize = 10,maxGSSize = 500,qvalueCutoff = 1,use_internal_data = FALSE)
#keu<-barplot(kegg,showCategory=15,drop=T,x = "GeneRatio",color = "pvalue")
keu<-barplot(kegg,showCategory=15,drop=T,x = "GeneRatio",color = "p.adjust")

ee<-kegg@result
ee1<-ee[which(ee$pvalue<0.05),]
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
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/1_RNAseq/GEM_res_S8495_res/pathway")
write.table(ee1,"GEM_res_S8495_res_kegg_DEGP0.05_pathway_P0.05_up_drug_target.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(pathway_gene2,"GEM_res_S8495_res_kegg_DEGP0.05_pathway_P0.05_up_drug_target_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
ggsave("GEM_res_S8495_res_kegg_DEGP0.05_pathway_P0.05_up_drug_target_gene.pdf",keu,width = 8, height = 6)


setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/1_RNAseq/GEM_res_S8495_res/pathway")
write.table(ee1,"GEM_res_S8495_res_kegg_DEGP0.05_pathway_P0.05_down_drug_target.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(pathway_gene2,"GEM_res_S8495_res_kegg_DEGP0.05_pathway_P0.05_down_drug_target_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
ggsave("GEM_res_S8495_res_kegg_DEGP0.05_pathway_P0.05_down_drug_target_gene.pdf",keu,width = 8, height = 6)

go<-enrichGO(gene[,2],OrgDb = org.Hs.eg.db,ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.05,keyType = 'ENTREZID')
ked<-dotplot(go,showCategory=10)
a<-go@result
go_BP<-a[which(a$p.adjust<0.05),]

enrich_genego<-go_BP$geneID
pathway_genego<-unique(unlist(strsplit(enrich_genego,split="/")))
pathway_genego2=bitr(pathway_genego,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/1_RNAseq/GEM_res_S8495_res/pathway")
write.table(pathway_genego2,"GEM_res_S8495_res_GO_DEG_pathway_FDR0.05_up_drug_target_gene.txt",sep = '\t',col.names = F,row.names = F,quote = FALSE)##数据输出
write.table(go_BP,"GEM_res_S8495_res_GO_DEG_pathway_FDR0.05_up_drug_target.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
ggsave("GEM_res_S8495_res_GO_DEG_pathway_FDR0.05_up_drug_target_gene.pdf",ked,width = 8, height = 6)

setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/1_RNAseq/GEM_res_S8495_res/pathway")
write.table(pathway_genego2,"GEM_res_S8495_res_GO_DEG_pathway_FDR0.05_down_drug_target_gene.txt",sep = '\t',col.names = F,row.names = F,quote = FALSE)##数据输出
write.table(go_BP,"GEM_res_S8495_res_GO_DEG_pathway_FDR0.05_down_drug_target.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
ggsave("GEM_res_S8495_res_GO_DEG_pathway_FDR0.05_down_drug_target_gene.pdf",ked,width = 8, height = 6)

##################  ridge regresion 找到的TF及其靶标DEGs  上下调的转录因子是分开看的
###############   TF找到对应的靶基因
rm(list=ls())
folder<-c("GEM_res_S8495_sen","GEM_res_S8495_res")
two_result<-NULL
for(k in 1:2){
  file_name=as.character(folder[k])
  cat(file_name,"\n")
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq/",file_name,"/FIMO",sep=""))
  up_TFs1<-read.table("FIMO_result_up_TFs_motifFDR0.05_DApeakP0.01.txt",sep = '\t',header=TRUE)
  DETFs<-as.character(up_TFs1[,1])
  T<-c("TRANSFAC","JASPAR","MotifMap","ENCODE")
  up_result<-matrix(0,4,5)
for(j in 1:4){
  file=T[j]
  up_result[j,1]<-file_name
  up_result[j,2]<-"up"
  up_result[j,3]<-file
  setwd("~/xjj/drug/drug_result/HDAC_frontiers/4_TFs/TRANSFAC_TF")
  target_TF<-read.delim(paste(file,"-gene_attribute_edges.txt",sep=""),header=TRUE,sep='\t')
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
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/1_RNAseq/",file_name,sep=""))
  up_gene<-read.table(paste(file_name,"_other",k+1,"_RNAseq_P0.05_FC0_up_DEGs.txt",sep=""),header=TRUE,sep="\t")
  down_gene<-read.table(paste(file_name,"_other",k+1,"_RNAseq_P0.05_FC0_down_DEGs.txt",sep=""),header=TRUE,sep="\t")

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
  colnames(data)<-c("TFs","target-DEGs","target-DEGs-dysregulation")
  aa<-unique(data[,1])
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
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq/",file_name,"/TF_target",sep=""))
  write.table(re,paste("number_of_up_",file,"-gene.txt",sep=""),sep = '\t',col.names = TRUE,row.names = F,quote = FALSE)##数据输出
  write.table(data,paste(file_name,"_up_TFs_",file,"-gene.txt",sep=""),sep = '\t',col.names = TRUE,row.names = F,quote = FALSE)##数据输出
  up_result[j,4]<-length(aa)
  up_result[j,5]<-paste0(aa,collapse =",")
}

setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq/",file_name,"/FIMO",sep=""))
down_TFs1<-read.table("FIMO_result_down_TFs_motifFDR0.05_DApeakP0.01.txt",sep = '\t',header=TRUE)
DETFs<-as.character(down_TFs1[,1])
T<-c("TRANSFAC","JASPAR","MotifMap","ENCODE")
down_result<-matrix(0,4,5)
for(j in 1:4){
  file=T[j]
  down_result[j,1]<-file_name
  down_result[j,2]<-"down"
  down_result[j,3]<-file
  setwd("~/xjj/drug/drug_result/HDAC_frontiers/4_TFs/TRANSFAC_TF")
  target_TF<-read.delim(paste(file,"-gene_attribute_edges.txt",sep=""),header=TRUE,sep='\t')
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
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/1_RNAseq/",file_name,sep=""))
  up_gene<-read.table(paste(file_name,"_other",k+1,"_RNAseq_P0.05_FC0_up_DEGs.txt",sep=""),header=TRUE,sep="\t")
  down_gene<-read.table(paste(file_name,"_other",k+1,"_RNAseq_P0.05_FC0_down_DEGs.txt",sep=""),header=TRUE,sep="\t")
  
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
  colnames(data)<-c("TFs","target-DEGs","target-DEGs-dysregulation")
  aa<-unique(data[,1])
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
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq/",file_name,"/TF_target",sep=""))
  write.table(re,paste("number_of_down_",file,"-gene.txt",sep=""),sep = '\t',col.names = TRUE,row.names = F,quote = FALSE)##数据输出
  write.table(data,paste(file_name,"_down_TFs_",file,"-gene.txt",sep=""),sep = '\t',col.names = TRUE,row.names = F,quote = FALSE)##数据输出
  down_result[j,4]<-length(aa)
  down_result[j,5]<-paste0(aa,collapse =",")
}
one_result<-rbind(up_result,down_result)
two_result<-rbind(two_result,one_result)
}
colnames(two_result)<-c("type","dysregulator","database of TFs","number of TF have target","list of target to TFs")
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq/GEM_res_S8495_res/TF_target")
write.table(two_result,"TFs_Targer_all",sep = '\t',col.names = TRUE,row.names = F,quote = FALSE)##数据输出

#########################################################################
####图形展示TFs-靶基因情况
rm(list=ls())
library(ggplot2)
#1
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq/GEM_res_S8495_sen/TF_target")
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
#setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq/GEM_res_S8495_sen/TF_target")
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
  geom_rect(aes(xmin=0.5,xmax=2.5,ymin=0,ymax=Inf),
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
#setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq/GEM_res_S8495_res/TF_target")
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
  geom_rect(aes(xmin=0.5,xmax=6.5,ymin=0,ymax=Inf),
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
#setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq/GEM_res_S8495_res/TF_target")
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

#######################################################################
###############  探索DApeaks与药物之间的关系  ####################################
################## ATACseq-drug correlation ######
rm(list=ls())
library(data.table)
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_RPKM<-read.table("peak_RPKM.txt",sep="\t",header=T) #没有全0的
colnames(peak_RPKM)<-gsub("\\.","-",colnames(peak_RPKM))
peak_name<-peak_RPKM[,1]

setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,!(colnames(ic501)%in%c("PC.100","PC.34"))]## all of PDAC
colnames(ic50)<-sample_id[match(colnames(ic50),sample_id[,2]),1]
want_AUC_chemo<-ic50[rownames(ic50)%in%"GEM",]
want_AUC_HDAC<-ic50[rownames(ic50)%in%"S8495",]

GEM_sensitive<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
                 "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
                 "DAC-19","DAC-37","DAC-39","DAC-38")
GEM_res_S8495_sen<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
GEM_res_S8495_res<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")
all<-c(GEM_sensitive,GEM_res_S8495_sen,GEM_res_S8495_res)
peak_RPKM37<-peak_RPKM[,match(all,colnames(peak_RPKM))]

folder<-c("GEM_res_S8495_sen","GEM_res_S8495_res")
two_result<-matrix(0,2,6)
for(m in 1:2){
  file_name=as.character(folder[m])
  cat(file_name,"\n")
  setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq/",file_name,sep=""))
  DApeak<-fread(paste(file_name,"-other",m+1,"-all-DESeq2_chemotherapy_AUC_ATAC.txt",sep=""),header=T,data.table=F)
  DAR<-unique(DApeak[DApeak$pvalue<0.01,1])
  DApeaks<-peak_name[match(DAR,peak_name)] #DApeaks and position
  DApeak_clinical<-peak_RPKM37[match(DAR,peak_name),]
  DApeak_clinical<-DApeak_clinical[-269,]
corralation_all_s<-NULL
corralation_all_p<-NULL
corralation_all1<-matrix(0,nrow(DApeak_clinical),7) #总得有26541个peaks，最终与GEM药物有显著相关的peaks只有2189个。
for(k in 1:nrow(DApeak_clinical)){
  corralation_all1[k,1]<-as.character(DApeaks[k])
  corralation_all1[k,2]<-"spearman"
  cor_result_spearman<-cor.test(as.numeric(want_AUC_HDAC),as.numeric(DApeak_clinical[k,]),alternative = "two.sided",method = "spearman")
  corralation_all1[k,3]<-cor_result_spearman$estimate
  corralation_all1[k,4]<-cor_result_spearman$p.value
  corralation_all1[k,5]<-"pearson"
  cor_result_pearson<-cor.test(as.numeric(want_AUC_HDAC),as.numeric(DApeak_clinical[k,]),alternative = "two.sided",method = "pearson")
  corralation_all1[k,6]<-cor_result_pearson$estimate
  corralation_all1[k,7]<-cor_result_pearson$p.value
}
corralation_all_s<-corralation_all1[which(as.numeric(corralation_all1[,4])<0.05),]
corralation_all_p<-corralation_all1[which(as.numeric(corralation_all1[,7])<0.05),]

colnames(corralation_all_s)<-c("Peak_name","spearman","spearman_estimate","spearman_pvalue","pearson","pearson_estimate","pearson_pvalue")
colnames(corralation_all_p)<-c("Peak_name","spearman","spearman_estimate","spearman_pvalue","pearson","pearson_estimate","pearson_pvalue")
setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/3_intersect/peaks_drug/",file_name,sep=""))
write.table(corralation_all_s,"GEM_DARs_drug_corralation_all_spearman_p0.05.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(corralation_all_p,"GEM_DARs_drug_corralation_all_pearson_p0.05.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(corralation_all1,"GEM_DARs_drug_corralation_all.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

two_result[m,1]<-file_name
two_result[m,2]<-length(DAR)
two_result[m,3]<-nrow(corralation_all_s)
two_result[m,4]<-nrow(corralation_all_s)/length(DAR)
two_result[m,5]<-nrow(corralation_all_p)
two_result[m,6]<-nrow(corralation_all_p)/length(DAR)
}
colnames(two_result)<-c("type","Total peaks","corralation_all_s","ratio of corralation_all_s","corralation_all_p","ratio of corralation_all_p")
setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/3_intersect/peaks_drug",sep=""))
write.table(two_result,"DARs_drug_corralation_all.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出


################## RNAseq-ATACseq-drug correlation ######
rm(list=ls())
library(data.table)
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_RPKM<-read.table("peak_RPKM.txt",sep="\t",header=T) #没有全0的
colnames(peak_RPKM)<-gsub("\\.","-",colnames(peak_RPKM))
peak_name<-peak_RPKM[,1]

setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,!(colnames(ic501)%in%c("PC.100","PC.34"))]## all of PDAC
colnames(ic50)<-sample_id[match(colnames(ic50),sample_id[,2]),1]
want_AUC_chemo<-ic50[rownames(ic50)%in%"GEM",]
want_AUC_HDAC<-ic50[rownames(ic50)%in%"S8495",]

#GEM_sensitive<-c("DAC-24","DAC-27","DAC-28","DAC-4","DAC-33","DAC-34","DAC-36","DAC-6","DAC-1","DAC-9","DAC-11","DAC-15","DAC-3","DAC-16","DAC-18","DAC-19","DAC-37","DAC-38")
#GEM_res_S8495_sen<-c("DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-12","DAC-14")
#GEM_res_S8495_res<-c("DAC-20","DAC-21","DAC-22","DAC-23","DAC-25","DAC-29","DAC-30","DAC-32","DAC-35","DAC-5","DAC-13","DAC-39")

GEM_sensitive<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
                 "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
                 "DAC-19","DAC-37","DAC-39","DAC-38")
GEM_res_S8495_sen<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
GEM_res_S8495_res<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")
all<-c(GEM_sensitive,GEM_res_S8495_sen,GEM_res_S8495_res)
peak_RPKM37<-peak_RPKM[,match(all,colnames(peak_RPKM))]

folder<-c("GEM_res_S8495_sen","GEM_res_S8495_res")
k=2
file_name<-folder[k]
setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq/",file_name,"/annotation",sep=""))
homer_anno_down<-fread(paste(file_name,"-other",k+1,"-AUC-P0.01-logFC0-DESeq2-down-annotation.txt",sep=""),header=T,data.table=F)
colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]

homer_anno_up<-fread(paste(file_name,"-other",k+1,"-AUC-P0.01-logFC0-DESeq2-up-annotation.txt",sep=""),header=T,data.table=F)
colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]

setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/1_RNAseq/",file_name,sep=""))
up_gene<-read.table(paste(file_name,"_other",k+1,"_RNAseq_P0.05_FC0_up_DEGs.txt",sep=""),header=T,sep="\t")
down_gene<-read.table(paste(file_name,"_other",k+1,"_RNAseq_P0.05_FC0_down_DEGs.txt",sep=""),header=T,sep="\t")

downDEGs_peaks<-c()
for(i in 1:nrow(down_gene)){
  DEGs_peaks1<-Anno_gene_100Kb_down[Anno_gene_100Kb_down$`Gene Name` %in% down_gene[i,1],1]
  downDEGs_peaks<-c(downDEGs_peaks,DEGs_peaks1)
}

upDEGs_peaks<-c()
for(i in 1:nrow(up_gene)){
  DEGs_peaks2<-Anno_gene_100Kb_up[Anno_gene_100Kb_up$`Gene Name` %in% up_gene[i,1],1]
  upDEGs_peaks<-c(upDEGs_peaks,DEGs_peaks2)
}
DAR<-unique(c(downDEGs_peaks,upDEGs_peaks))
DApeaks<-peak_name[match(DAR,peak_name)] #DApeaks and position
DApeak_clinical<-peak_RPKM37[match(DAR,peak_name),]

corralation_all_s<-NULL
corralation_all_p<-NULL
corralation_all1<-matrix(0,nrow(DApeak_clinical),7) #总得有2346个peaks，最终与GEM药物有显著相关的peaks只有201个。
for(k in 1:nrow(DApeak_clinical)){
  corralation_all1[k,1]<-as.character(DApeaks[k])
  corralation_all1[k,2]<-"spearman"
  cor_result_spearman<-cor.test(as.numeric(want_AUC_chemo),as.numeric(DApeak_clinical[k,]),alternative = "two.sided",method = "spearman")
  corralation_all1[k,3]<-cor_result_spearman$estimate
  corralation_all1[k,4]<-cor_result_spearman$p.value
  corralation_all1[k,5]<-"pearson"
  cor_result_pearson<-cor.test(as.numeric(want_AUC_chemo),as.numeric(DApeak_clinical[k,]),alternative = "two.sided",method = "pearson")
  corralation_all1[k,6]<-cor_result_pearson$estimate
  corralation_all1[k,7]<-cor_result_pearson$p.value
}
corralation_all_s<-corralation_all1[which(as.numeric(corralation_all1[,4])<0.05),]
corralation_all_p<-corralation_all1[which(as.numeric(corralation_all1[,7])<0.05),]

colnames(corralation_all_s)<-c("Peak_name","spearman","spearman_estimate","spearman_pvalue","pearson","pearson_estimate","pearson_pvalue")
colnames(corralation_all_p)<-c("Peak_name","spearman","spearman_estimate","spearman_pvalue","pearson","pearson_estimate","pearson_pvalue")
setwd(paste("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/3_intersect/peaks_drug/",file_name,sep=""))
write.table(corralation_all_s,"S8495_res_DARs_DEGs_drug_corralation_all_spearman_p0.05.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(corralation_all_p,"S8495_res_DEGs_drug_corralation_all_pearson_p0.05.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(corralation_all1,"S8495_res_DEGs_drug_corralation_all.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

#############################################################################
#########  mutation 画图
rm(list=ls())
library(maftools)
GEM_sensitive<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
                 "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
                 "DAC-19","DAC-37","DAC-39","DAC-38")
GEM_res_S8495_sen<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
GEM_res_S8495_res<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")

sensitivity_new<-GEM_res_S8495_res
resistance_new<-c(GEM_sensitive,GEM_res_S8495_sen)
setwd("~/176/Changhai_WGS/results/2_Variants/SomaticVcfs")
laml = read.maf(maf = "WGS.maf.gz")
var_maf<-subsetMaf(maf = laml, tsb = c(resistance_new,sensitivity_new)) ##取子集
write.table(laml@data,"WGS.Mut.xls",row.names = F,col.names = T,sep="\t")
####Performs Pair-wise Fisher's Exact test to detect mutually exclusive or co-occuring events.
##BrBG PiYG PRGn PuOr RdBu RdGy RdYlBu RdYlGn Spectral
output <- somaticInteractions(maf=var_maf, top=20, pvalue=c(0.05, 0.01))
luad.sig <- oncodrive(maf=var_maf, minMut=5,pvalMethod="zscore")
#显示特定基因
col = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
oncoplot(maf = var_maf, writeMatrix=T,colors ="PuOr")
oncostrip(maf = var_maf, writeMatrix=T,colors = col,removeNonMutated = F,showTumorSampleBarcodes=T)
#titv函数将SNP分类为Transitions_vs_Transversions，并以各种方式返回汇总表的列表。 汇总数据也可以显示为一个箱线图，显示六种不同转换的总体分布，并作为堆积条形图显示每个样本中的转换比例。
titv(var_maf, useSyn = FALSE, plot = TRUE, file = NULL)


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
GEM_sensitive<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
                 "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
                 "DAC-19","DAC-37","DAC-39","DAC-38")
GEM_res_S8495_sen<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
GEM_res_S8495_res<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")

sensitivity<-GEM_res_S8495_res
resistance<-c(GEM_sensitive,GEM_res_S8495_sen)

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
data_result<-data[which(data[,2]<0.1),]  
if(nrow(data_result)>0){
  mutation_result[1,2]<-nrow(data_result)
  mutation_result[1,3]<-paste0(data_result[,1],collapse =",")
}
if(nrow(data_result)==0){
  mutation_result[1,2]<-0
  mutation_result[1,3]<-0
}
mutation_result[1,1]<-"GEM_res_S8495_res"
colnames(mutation_result)<-c("drug","The number of gene mutation","mutation gene")
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/4_WGS/GEM_res_S8495_res")
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

GEM_sensitive<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
                 "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
                 "DAC-19","DAC-37","DAC-39","DAC-38")
GEM_res_S8495_sen<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
GEM_res_S8495_res<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")

sensitivity<-GEM_res_S8495_sen
resistance<-c(GEM_sensitive,GEM_res_S8495_res)
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
a1<-tresult[which(tresult[,4]<0.05),]
up<-a1[which(a1[,2]>0),1] #sensitivity
down<-a1[which(a1[,2]<0),1] #resistance
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/4_WGS/GEM_res_S8495_sen")
write.table(a1,"cnv_result_T_test.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/4_WGS/GEM_res_S8495_sen")
amp_gene1<-read.table("cnv_result_T_test.txt",sep = '\t',header = T)
library("data.table")
library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(stringr)
library(foreach)
genea=bitr(amp_gene1[,1],fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
kegga<-enrichKEGG(genea[,2],organism = 'hsa',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'none',
                  minGSSize = 10,maxGSSize = 500,qvalueCutoff = 1,use_internal_data = FALSE)
ka<-barplot(kegga,showCategory=10,drop=T,x = "GeneRatio",color = "pvalue")
eea<-kegga@result
eea1<-eea[which(eea$pvalue<0.05),]
for(y in 1:nrow(eea1)){
  eea1[y,8]
  b1<-matrix(unlist(strsplit(eea1[y,8],split="/")),ncol=1)
  gene6=bitr(b1,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
  eea1[y,10]<-paste0(gene6[,2],collapse ="/")
}
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/4_WGS/GEM_res_S8495_sen")
write.table(eea1,"kegg_CNV_amp_p0.05_pathway_P0.05.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
library(ggplot2)
ggsave(paste("kegg_CNV_amp_pathway_P0.05_Ttest_DEGs_P0.05.pdf",sep=""),ka,width = 8, height = 6)



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
GEM_sensitive<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
                 "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
                 "DAC-19","DAC-37","DAC-39","DAC-38")
GEM_res_S8495_sen<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
GEM_res_S8495_res<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")

sensitivity<-GEM_res_S8495_sen
resistance<-c(GEM_sensitive,GEM_res_S8495_res)
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
data_result<-data[which(data[,2]<0.1),]  
colnames(data_result)<-c("gene","pvalue")
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/4_WGS/GEM_res_S8495_sen")
write.table(data_result,"CNV_fisher_resultP0.1.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(data,"CNV_fisher_result_all.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/4_WGS/GEM_res_S8495_sen")
amp_gene1<-read.table("CNV_fisher_resultP0.1.txt",sep = '\t',header = T)
library("data.table")
library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(stringr)
library(foreach)
genea=bitr(amp_gene1[,1],fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
kegga<-enrichKEGG(genea[,2],organism = 'hsa',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'none',
                  minGSSize = 10,maxGSSize = 500,qvalueCutoff = 1,use_internal_data = FALSE)
ka<-barplot(kegga,showCategory=10,drop=T,x = "GeneRatio",color = "pvalue")
eea<-kegga@result
eea1<-eea[which(eea$pvalue<0.05),]
for(y in 1:nrow(eea1)){
  eea1[y,8]
  b1<-matrix(unlist(strsplit(eea1[y,8],split="/")),ncol=1)
  gene6=bitr(b1,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
  eea1[y,10]<-paste0(gene6[,2],collapse ="/")
}

setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/4_WGS/GEM_res_S8495_sen")
write.table(eea1,"kegg_CNV_fisher_amp_p0.1_pathway_P0.05.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
library(ggplot2)
ggsave(paste("kegg_CNV_fisher_amp_pathway_P0.1_DEGs_P0.05.pdf",sep=""),ka,width = 8, height = 6)


###########################################################################
############  WGS与ATACseq的联合分析 在peaks中的突变位点很少，那看看非编码突变的数据怎么样
BiocManager::install("maftools")
rm(list=ls())
library(maftools)
setwd("~/176/Changhai_WGS/results/2_Variants/SomaticVcfs")
laml = read.maf(maf = "WGS.maf.gz")
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_GEM/3_WGS/mutation")
write.table(laml@data,"WGS.Mut.xls",row.names = F,col.names = T,sep="\t")

data1<-laml@data
GEM_sensitive<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
                 "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
                 "DAC-19","DAC-37","DAC-39","DAC-38")
GEM_res_S8495_sen<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
GEM_res_S8495_res<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")
sample<-c(GEM_sensitive,GEM_res_S8495_res,GEM_res_S8495_sen)
sample1<-intersect(data1$Tumor_Sample_Barcode,sample)
data<-data1[data1$Tumor_Sample_Barcode%in%sample1,]



df<-data[,c(2,3,8)]
df <- df %>% distinct(Chromosome,Start_Position, .keep_all = T)
unique(df[,1])

setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq/GEM_res_S8495_sen")
library("data.table")
res<-fread("GEM_res_S8495_sen-other2-all-DESeq2_chemotherapy_AUC_ATAC.txt",header=T,data.table=F)
resSig<-res[which(res$pvalue<0.01),]
# pvalue padj
##新增一列，将log2FoldChange>0标注为up，<0标准为down
resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
sum(resSig$up_down=='up')
sum(resSig$up_down=='down')
up_gene<-as.data.frame(resSig[which(resSig$log2FoldChange>0),1])
down_gene<-as.data.frame(resSig[which(resSig$log2FoldChange<0),1])

library(dplyr)
a_up<-apply(up_gene,1,function(x) unlist(strsplit(as.character(x), "_")))
peaks_up<-t(a_up)
colnames(peaks_up)<-c("chr","start","end")

result<-NULL
for(m in 1:nrow(df)){
  mu<-df[m,]
  want_peak<-peaks_up[peaks_up[,1]%in%as.character(mu[1,1]),]
  pos<-sum(as.numeric(want_peak[,2])-as.numeric(mu[1,2])<0 & as.numeric(want_peak[,3])-as.numeric(mu[1,2])>0)
  result1<-matrix(c(as.character(mu[1,1]),as.numeric(mu[1,2]),as.character(as.data.frame(mu[1,3])[1,1]),pos),nrow=1)
  result<-rbind(result,result1)
}
unique(result[,4])
sum(result[,4]=="1")###有24个突变位点位于上调的peaks中
up_result<-result[which(result[,4]=="1"),] 
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/3_intersect/WGS_ATACseq/GEM_res_S8495_sen")
write.table(up_result,"GEM_res_S8495_sen_mutation_ATACseq_up_result.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

library(dplyr)
a_down<-apply(down_gene,1,function(x) unlist(strsplit(as.character(x), "_")))
peaks_down<-t(a_down)
colnames(peaks_down)<-c("chr","start","end")

result_down<-NULL
for(k in 1:nrow(df)){
  mu<-df[k,]
  want_peak<-peaks_down[peaks_down[,1]%in%as.character(mu[1,1]),]
  pos<-sum(as.numeric(want_peak[,2])-as.numeric(mu[1,2])<0 & as.numeric(want_peak[,3])-as.numeric(mu[1,2])>0)
  result1<-matrix(c(as.character(mu[1,1]),as.numeric(mu[1,2]),as.character(as.data.frame(mu[1,3])[1,1]),pos),nrow=1)
  result_down<-rbind(result_down,result1)
}
unique(result_down[,4])
sum(result_down[,4]=="1") ###有24个突变位点位于上调的peaks中
down_result<-result_down[which(result_down[,4]=="1"),] 
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/3_intersect/WGS_ATACseq/GEM_res_S8495_sen")
write.table(down_result,"GEM_res_S8495_sen_mutation_ATACseq_down_result.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

#所占的比例，上下调中分别是0.3%，总的是占0.6%


####   去吃饭的时候跑
############### 非编码基因突变与ATACseq之间的联系，数据从网上下载
rm(list=ls())
library(data.table)
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_GEM/3_WGS/mutation")
data1<-fread("nocoding_mutation.txt",header=T,data.table=F)
GEM_sensitive<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
                 "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
                 "DAC-19","DAC-37","DAC-39","DAC-38")
GEM_res_S8495_sen<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
GEM_res_S8495_res<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")
sample<-c(GEM_sensitive,GEM_res_S8495_res,GEM_res_S8495_sen)
sample1<-intersect(gsub("CAS-","",data1$Sample_ID),sample)
data<-data1[gsub("CAS-","",data1$Sample_ID)%in%sample1,]

df<-data[,c(1,2,6)]
df <- df %>% distinct(Chr,Start, .keep_all = T)
unique(df[,1])
df<-df[-which(df[,1]=="chrM"),]

setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/2_ATACseq/GEM_res_S8495_sen")
library("data.table")
res<-fread("GEM_res_S8495_sen-other2-all-DESeq2_chemotherapy_AUC_ATAC.txt",header=T,data.table=F)
resSig<-res[which(res$pvalue<0.01),]
# pvalue padj
##新增一列，将log2FoldChange>0标注为up，<0标准为down
resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
sum(resSig$up_down=='up')
sum(resSig$up_down=='down')
up_gene<-as.data.frame(resSig[which(resSig$log2FoldChange>0),1])
down_gene<-as.data.frame(resSig[which(resSig$log2FoldChange<0),1])

library(dplyr)
a_up<-apply(up_gene,1,function(x) unlist(strsplit(as.character(x), "_")))
peaks_up<-t(a_up)
colnames(peaks_up)<-c("chr","start","end")

result<-matrix(0,nrow(df),4)
for(k in 1:nrow(df)){
  mu<-df[k,]
  want_peak<-peaks_up[peaks_up[,1]%in%as.character(mu[1,1]),]
  pos<-sum(as.numeric(want_peak[,2])-as.numeric(mu[1,2])<0 & as.numeric(want_peak[,3])-as.numeric(mu[1,2])>0)
  result[k,1]<-as.character(mu[1,1])
  result[k,2]<-as.numeric(mu[1,2])
  result[k,3]<-mu[1,3]
  result[k,4]<-pos
}
sum(result[,4]=="1")###有24个突变位点位于上调的peaks中

non_up<-result[which(result[,4]=="1"),]
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/3_intersect/WGS_ATACseq/GEM_res_S8495_sen")
write.table(non_up,"GEM_res_S8495_sen_non_mutation_ATACseq_up_result.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出


library(dplyr)
a_down<-apply(down_gene,1,function(x) unlist(strsplit(as.character(x), "_")))
peaks_down<-t(a_down)
colnames(peaks_down)<-c("chr","start","end")

result_down<-matrix(0,nrow(df),4)
for(k in 1:nrow(df)){
  mu<-df[k,]
  want_peak<-peaks_down[peaks_down[,1]%in%as.character(mu[1,1]),]
  pos<-sum(as.numeric(want_peak[,2])-as.numeric(mu[1,2])<0 & as.numeric(want_peak[,3])-as.numeric(mu[1,2])>0)
  result_down[k,1]<-as.character(mu[1,1])
  result_down[k,2]<-as.numeric(mu[1,2])
  result_down[k,3]<-as.character(mu[1,3])
  result_down[k,4]<-pos
}
sum(result_down[,4]=="1") ###有24个突变位点位于上调的peaks中
#所占的比例，上下调中分别是0.3%，总的是占0.6%
non_down<-result_down[which(result_down[,4]=="1"),]
setwd("~/xjj/drug/drug_result/HDACi_chemo618/GEM_S8495_three_categories/3_intersect/WGS_ATACseq/GEM_res_S8495_sen")
write.table(result,"GEM_res_S8495_sen_non_mutation_ATACseq_up_result.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出







