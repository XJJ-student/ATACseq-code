############## ATACseq 热图
##用这些特征画热图，还有一点点问题，就是每一类不是单独那一类peak值比较大
rm(list=ls())
library(data.table)
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_RPKM<-read.table("peak_RPKM.txt",sep="\t",header=T) #没有全0的
colnames(peak_RPKM)<-gsub("\\.","-",colnames(peak_RPKM))
GEM_sensitive<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
                 "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
                 "DAC-19","DAC-37","DAC-39","DAC-38")
GEM_res_S8495_sen<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
GEM_res_S8495_res<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")
all<-c(GEM_sensitive,GEM_res_S8495_sen,GEM_res_S8495_res)
library(data.table)
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_GEM/2_ATACseq/annotation")
homer_anno_up<-fread("sensitivity-resistance-GEM-AUC-P0.01-logFC0-DESeq2-up-annotation.txt",header=T,data.table=F)
colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
length(unique(Anno_gene_100Kb_up$`Gene Name`))
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_GEM/1_RNAseq")
up_gene<-read.table("sensitivity-resistance-GEM-DESeq2_AUC_RNAseq_P0.05_FC0_up_DEGs.txt",header=T,sep="\t")
DEGs<-as.character(up_gene[,1])
int_all_DEGs<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)

all_peaks_anno<-homer_anno_up
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_GEM/2_ATACseq")
DApeak<-fread("sensitivity-resistance-all-DESeq2_GEM_AUC_ATAC.txt",header=T,data.table=F)
different_FC<-NULL
for(i in 1:length(int_all_DEGs)){
  DEGs_peaks1<-all_peaks_anno[all_peaks_anno$`Gene Name` %in% int_all_DEGs[i],1]
  DEGs_peaks_FC<-as.data.frame(DApeak[match(DEGs_peaks1,DApeak$peak_id),3])
  gene_name<-matrix(rep(int_all_DEGs[i],length(DEGs_peaks1)),ncol=1)
  DEGs_FC<-matrix(rep(up_gene[match(int_all_DEGs[i],DEGs),3],ncol=1))
  different_FC1<-cbind(gene_name,DEGs_FC,DEGs_peaks1,DEGs_peaks_FC)
  different_FC<-rbind(different_FC,different_FC1)
}
sum(different_FC[,2]>0) #up DEGs
sum(different_FC[,4]>0) #up DARs
colnames(different_FC)<-c("GEM_sensitive_DEGs","DEGs_FC","GEM_sensitive_DARs","DARs_FC")
GEM_sensitive_DARs<-as.character(different_FC[which(different_FC[,2]>0 & different_FC[,4]>0),3])

setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495/2_ATACseq/annotation")
homer_anno_up<-fread("sensitivity-resistance-S8495-AUC-P0.01-logFC0-DESeq2-down-annotation.txt",header=T,data.table=F)
colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
length(unique(Anno_gene_100Kb_up$`Gene Name`))
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495/1_RNAseq")
up_gene<-read.table("sensitivity-resistance-S8495-DESeq2_AUC_RNAseq_P0.05_FC0_down_DEGs.txt",header=T,sep="\t")
DEGs<-as.character(up_gene[,1])
int_all_DEGs<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)

all_peaks_anno<-homer_anno_up
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495/2_ATACseq")
DApeak<-fread("sensitivity-resistance-all-DESeq2_S8495_AUC_ATAC.txt",header=T,data.table=F)

different_FC<-NULL
for(i in 1:length(int_all_DEGs)){
  DEGs_peaks1<-all_peaks_anno[all_peaks_anno$`Gene Name` %in% int_all_DEGs[i],1]
  DEGs_peaks_FC<-as.data.frame(DApeak[match(DEGs_peaks1,DApeak$peak_id),3])
  gene_name<-matrix(rep(int_all_DEGs[i],length(DEGs_peaks1)),ncol=1)
  DEGs_FC<-matrix(rep(up_gene[match(int_all_DEGs[i],DEGs),3],ncol=1))
  different_FC1<-cbind(gene_name,DEGs_FC,DEGs_peaks1,DEGs_peaks_FC)
  different_FC<-rbind(different_FC,different_FC1)
}
sum(different_FC[,2]<0) #up DEGs
sum(different_FC[,4]<0) #up DARs
colnames(different_FC)<-c("S8495_resistant_DEGs","DEGs_FC","S8495_resistant_DARs","DARs_FC")
S8495_resistant_DARs<-as.character(different_FC[which(different_FC[,2]<0 & different_FC[,4]<0),3])

setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495/2_ATACseq/annotation")
homer_anno_up<-fread("sensitivity-resistance-S8495-AUC-P0.01-logFC0-DESeq2-up-annotation.txt",header=T,data.table=F)
colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
length(unique(Anno_gene_100Kb_up$`Gene Name`))

setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495/1_RNAseq")
up_gene<-read.table("sensitivity-resistance-S8495-DESeq2_AUC_RNAseq_P0.05_FC0_up_DEGs.txt",header=T,sep="\t")
DEGs<-as.character(up_gene[,1])
int_all_DEGs<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)
all_peaks_anno<-homer_anno_up

different_FC<-NULL
for(i in 1:length(int_all_DEGs)){
  DEGs_peaks1<-all_peaks_anno[all_peaks_anno$`Gene Name` %in% int_all_DEGs[i],1]
  DEGs_peaks_FC<-as.data.frame(DApeak[match(DEGs_peaks1,DApeak$peak_id),3])
  gene_name<-matrix(rep(int_all_DEGs[i],length(DEGs_peaks1)),ncol=1)
  DEGs_FC<-matrix(rep(up_gene[match(int_all_DEGs[i],DEGs),3],ncol=1))
  different_FC1<-cbind(gene_name,DEGs_FC,DEGs_peaks1,DEGs_peaks_FC)
  different_FC<-rbind(different_FC,different_FC1)
}
sum(different_FC[,2]>0) #up DEGs
sum(different_FC[,4]>0) #up DARs
colnames(different_FC)<-c("S8495_sensitive_DEGs","DEGs_FC","S8495_sensitive_DARs","DARs_FC")
S8495_sensitive_DARs<-as.character(different_FC[which(different_FC[,2]>0 & different_FC[,4]>0),3])

intersect(S8495_sensitive_DARs,S8495_resistant_DARs)
a<-intersect(GEM_sensitive_DARs,S8495_sensitive_DARs)
if(length(a)==0){
  GEM_sensitive_only1<-GEM_sensitive_DARs
  S8495_sensitive_only<-S8495_sensitive_DARs
}
if(length(a)!=0){
  GEM_sensitive_only1<-GEM_sensitive_DARs[-match(a,GEM_sensitive_DARs)]
  S8495_sensitive_only<-S8495_sensitive_DARs[-match(a,S8495_sensitive_DARs)]
}

b<-intersect(GEM_sensitive_only1,S8495_resistant_DARs)
if(length(b)==0){
  GEM_sensitive_only<-GEM_sensitive_only1
  S8495_resistant_only<-S8495_resistant_DARs
}
if(length(b)!=0){
  GEM_sensitive_only<-GEM_sensitive_only1[-match(b,GEM_sensitive_only1)]
  S8495_resistant_only<-S8495_resistant_DARs[-match(b,S8495_resistant_DARs)]
}
DARs<-c(GEM_sensitive_only,S8495_sensitive_only,S8495_resistant_only)
#GEM_sensitive<-c("DAC-24","DAC-27","DAC-28","DAC-4","DAC-33","DAC-34","DAC-36","DAC-6","DAC-1","DAC-9","DAC-11","DAC-15","DAC-3","DAC-16","DAC-18","DAC-19","DAC-37","DAC-38")
#GEM_res_S8495_sen<-c("DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-12","DAC-14")
#GEM_res_S8495_res<-c("DAC-20","DAC-21","DAC-22","DAC-23","DAC-25","DAC-29","DAC-30","DAC-32","DAC-35","DAC-5","DAC-13","DAC-39")

GEM_sensitive<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
                 "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
                 "DAC-19","DAC-37","DAC-39","DAC-38")
GEM_res_S8495_sen<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
GEM_res_S8495_res<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")
all<-c(GEM_sensitive,GEM_res_S8495_sen,GEM_res_S8495_res)
sample_order<-c(GEM_sensitive,GEM_res_S8495_sen,GEM_res_S8495_res)

setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_RPKM<-read.table("peak_RPKM.txt",sep="\t",header=T)
peak_RPKM<-read.table("peak_count.txt",sep="\t",header=T)

colnames(peak_RPKM)<-gsub("\\.","-",colnames(peak_RPKM))

peaks<-c(as.character(GEM_sensitive_only),as.character(S8495_sensitive_only),as.character(S8495_resistant_only))
peak_RPKM1<-peak_RPKM[match(peaks,peak_RPKM[,1]),]
peaks1<-matrix(c(rep("Class1",length(GEM_sensitive_only)),rep("Class2",length(S8495_sensitive_only)),rep("Class3",length(S8495_resistant_only))),ncol=1)
peak_RPKM_order<-peak_RPKM1[,match(sample_order,colnames(peak_RPKM1))]

col_cut=c(length(GEM_sensitive),length(GEM_sensitive)+length(GEM_res_S8495_sen))
row_cut=c(length(GEM_sensitive_only),length(GEM_sensitive_only)+length(S8495_sensitive_only))
library(pheatmap)
data0=cbind(peaks1,peak_RPKM_order) ##0值很多

anno1=c(rep("Class1",length(GEM_sensitive)),rep("Class2",length(GEM_res_S8495_sen)),rep("Class3",length(GEM_res_S8495_res)))
anno1=data.frame(anno1)

rownames(anno1)=as.character(t(colnames(data0))[2:38])
ann_colors = list(
  anno1 = c(Class1 = "#4dbbd5", Class2 = "#00a087",Class3 = "#f39b7f"), anno=c("Class1"= "#4dbbd5","Class2"= "#00a087","Class3"= "#f39b7f"))
H=6
W=5
library(RColorBrewer)
color1=colorRampPalette(c("#3363aa","#1384d5","#23b2ae"))(3)
color10=colorRampPalette(c(color1[1],color1[2]))(35)
color11=colorRampPalette(c(color1[2],color1[3]))(15)
color2=colorRampPalette(c( "#23b2ae","#c5bc5e","#faf513"))(3)
color20=colorRampPalette(c(color2[1],color2[3]))(15)
color21=colorRampPalette(c(color2[3],color2[4]))(35)
mycolor=c(color10,color11,color20,color21)
anno=data.frame(data0$peaks1)
rownames(data0)=as.character(1:dim(data0)[1])
rownames(anno)=rownames(data0)
colnames(anno)="anno"
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_GEM/4_GEM_S8495_pheatmap")
p<-pheatmap(data0[2:38],scale="row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize=7,fontsize_row = 70*H/dim(data0)[1], fontsize_col = 7,cluster_cols=FALSE, cluster_rows=FALSE,annotation_col=anno1,gaps_col = col_cut,gaps_row = row_cut, annotation_row=anno,annotation_colors = ann_colors)

p<-pheatmap(data0[2:38],scale="row",color = mycolor,fontsize=7,fontsize_row = 70*H/dim(data0)[1], fontsize_col = 7,cluster_cols=FALSE, cluster_rows=FALSE,annotation_col=anno1,gaps_col = col_cut,gaps_row = row_cut, annotation_row=anno,annotation_colors = ann_colors)
ggsave("ATACseq_int_RNAseq_heatmap.pdf",p,width = 10, height = 8)
dev.off()


###################### 37个样本  单独使用peak作为特征  ##################################
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_GEM/2_ATACseq")
GEM<-read.table(paste("sensitivity-resistance-all-DESeq2_GEM_AUC_ATAC.txt",sep=""),sep = '\t',header=T)##数据输出
resSig<-GEM[which(GEM$pvalue<0.05 & abs(GEM$log2FoldChange)>2),]
#resSig<-GEM[which(GEM$pvalue<0.01),]
resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
cat(sum(resSig$up_down=='up'),"\n")
GEM_sensitive_DARs<-as.character(resSig[which(resSig$log2FoldChange>0),1])

setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_S8495/2_ATACseq")
S8495<-read.table(paste("sensitivity-resistance-all-DESeq2_S8495_AUC_ATAC.txt",sep=""),sep = '\t',header=T)##数据输出
resSig_S8495<-S8495[which(S8495$pvalue<0.05 & abs(S8495$log2FoldChange)>2),]
#resSig<-GEM[which(GEM$pvalue<0.01),]
resSig_S8495[which(resSig_S8495$log2FoldChange>0),'up_down']<-'up'
resSig_S8495[which(resSig_S8495$log2FoldChange<0),'up_down']<-'down'
cat(sum(resSig_S8495$up_down=='up'),"\n")
cat(sum(resSig_S8495$up_down=='down'),"\n")
S8495_sensitive_DARs<-as.character(resSig_S8495[which(resSig_S8495$log2FoldChange>0),1])
S8495_resistant_DARs<-as.character(resSig_S8495[which(resSig_S8495$log2FoldChange<0),1])
a<-intersect(GEM_sensitive_DARs,S8495_sensitive_DARs)
if(length(a)==0){
  GEM_sensitive_only1<-GEM_sensitive_DARs
  S8495_sensitive_only<-S8495_sensitive_DARs
}
if(length(a)!=0){
  GEM_sensitive_only1<-GEM_sensitive_DARs[-match(a,GEM_sensitive_DARs)]
  S8495_sensitive_only<-S8495_sensitive_DARs[-match(a,S8495_sensitive_DARs)]
}

b<-intersect(GEM_sensitive_only1,S8495_resistant_DARs)
if(length(b)==0){
  GEM_sensitive_only<-GEM_sensitive_only1
  S8495_resistant_only<-S8495_resistant_DARs
}
if(length(b)!=0){
  GEM_sensitive_only<-GEM_sensitive_only1[-match(b,GEM_sensitive_only1)]
  S8495_resistant_only<-S8495_resistant_DARs[-match(b,S8495_resistant_DARs)]
}
GEM_sensitive<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
                 "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
                 "DAC-19","DAC-37","DAC-39","DAC-38")
GEM_res_S8495_sen<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
GEM_res_S8495_res<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")
####median
#GEM_sensitive<-c("DAC-24","DAC-27","DAC-28","DAC-4","DAC-33","DAC-34","DAC-36","DAC-6","DAC-1","DAC-9","DAC-11","DAC-15","DAC-3","DAC-16","DAC-18","DAC-19","DAC-37","DAC-38")
#GEM_res_S8495_sen<-c("DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-12","DAC-14")
#GEM_res_S8495_res<-c("DAC-20","DAC-21","DAC-22","DAC-23","DAC-25","DAC-29","DAC-30","DAC-32","DAC-35","DAC-5","DAC-13","DAC-39")

sample_order<-c(GEM_sensitive,GEM_res_S8495_sen,GEM_res_S8495_res)

setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_RPKM<-read.table("peak_RPKM.txt",sep="\t",header=T)
peak_RPKM<-read.table("peak_count.txt",sep="\t",header=T)

colnames(peak_RPKM)<-gsub("\\.","-",colnames(peak_RPKM))

peaks<-c(as.character(GEM_sensitive_only),as.character(S8495_sensitive_only),as.character(S8495_resistant_only))
peak_RPKM1<-peak_RPKM[match(peaks,peak_RPKM[,1]),]
peaks1<-matrix(c(rep("Class1",length(GEM_sensitive_only)),rep("Class2",length(S8495_sensitive_only)),rep("Class3",length(S8495_resistant_only))),ncol=1)
peak_RPKM_order<-peak_RPKM1[,match(sample_order,colnames(peak_RPKM1))]

col_cut=c(length(GEM_sensitive),length(GEM_sensitive)+length(GEM_res_S8495_sen))
row_cut=c(length(GEM_sensitive_only),length(GEM_sensitive_only)+length(S8495_sensitive_only))
library(pheatmap)
data0=cbind(peaks1,peak_RPKM_order) ##0值很多

anno1=c(rep("Class1",length(GEM_sensitive)),rep("Class2",length(GEM_res_S8495_sen)),rep("Class3",length(GEM_res_S8495_res)))
anno1=data.frame(anno1)

rownames(anno1)=as.character(t(colnames(data0))[2:38])
ann_colors = list(
  anno1 = c(Class1 = "#4dbbd5", Class2 = "#00a087",Class3 = "#f39b7f"), anno=c("Class1"= "#4dbbd5","Class2"= "#00a087","Class3"= "#f39b7f"))
H=6
W=5
library(RColorBrewer)
color1=colorRampPalette(c("#3363aa","#1384d5","#23b2ae"))(3)
color10=colorRampPalette(c(color1[1],color1[2]))(35)
color11=colorRampPalette(c(color1[2],color1[3]))(15)
color2=colorRampPalette(c( "#23b2ae","#c5bc5e","#faf513"))(3)
color20=colorRampPalette(c(color2[1],color2[3]))(15)
color21=colorRampPalette(c(color2[3],color2[4]))(35)
mycolor=c(color10,color11,color20,color21)
anno=data.frame(data0$peaks1)
rownames(data0)=as.character(1:dim(data0)[1])
rownames(anno)=rownames(data0)
colnames(anno)="anno"
#setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_GEM/4_GEM_S8495_pheatmap")
p<-pheatmap(data0[2:38],scale="row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize=7,fontsize_row = 70*H/dim(data0)[1], fontsize_col = 7,cluster_cols=FALSE, cluster_rows=FALSE,annotation_col=anno1,gaps_col = col_cut,gaps_row = row_cut, annotation_row=anno,annotation_colors = ann_colors)

p<-pheatmap(data0[2:38],scale="row",color = mycolor,fontsize=7,fontsize_row = 70*H/dim(data0)[1], fontsize_col = 7,cluster_cols=FALSE, cluster_rows=FALSE,annotation_col=anno1,gaps_col = col_cut,gaps_row = row_cut, annotation_row=anno,annotation_colors = ann_colors)
setwd("~/xjj/drug/drug_result/HDACi_chemo618/sen_or_res_to_GEM_median/3_GEM_S8495_pheatmap")
library(ggplot2)
ggsave("ATACseqheatmap.pdf",p,width = 10, height = 8)
dev.off()
#### 注意要对每一行的数据进行标准化




#########################  留一法构造 ########################################################
#SVM 独立训练 独立验证 leave one out
setwd("~/xjj/drug/drug_result/HDACi_chemo618")
first_category_name = list.files("classifier_median_leave_one")  
dir = paste("~/xjj/drug/drug_result/HDACi_chemo618/classifier_median_leave_one/",first_category_name,sep="")  
n = length(dir) 
result<-c()
for(i in 1:n){      
  setwd(dir[i])
  facets<-read.table("text.txt",sep='\t',header=F)
  a<-as.character(facets[1,1])
  result<-c(result,a)
}

rm(list=ls())
library(DESeq2)
library("data.table")
library(dplyr)
#setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
#peak_count_name<-read.table("peak_count.txt",sep="\t",header=T) #没有全0的
#colnames(peak_count_name)<-gsub("\\.","-",colnames(peak_count_name))

#GEM_sensitive<-c("DAC-24","DAC-27","DAC-28","DAC-4","DAC-33","DAC-34","DAC-36","DAC-6","DAC-1","DAC-9","DAC-11","DAC-15","DAC-3","DAC-16","DAC-18","DAC-19","DAC-37","DAC-38")
#GEM_res_S8495_sen<-c("DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-12","DAC-14")
#GEM_res_S8495_res<-c("DAC-20","DAC-21","DAC-22","DAC-23","DAC-25","DAC-29","DAC-30","DAC-32","DAC-35","DAC-5","DAC-13","DAC-39")


GEM_sensitive<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
                 "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
                 "DAC-19","DAC-37","DAC-39","DAC-38")
GEM_res_S8495_sen<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
GEM_res_S8495_res<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")
all<-c(GEM_sensitive,GEM_res_S8495_sen,GEM_res_S8495_res)

setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
library("data.table")
exp11<-fread("star_rsem.GeneSymbol.readscount.xls",header=T,data.table=F)
sensitivity<-c("DAC-24","DAC-27","DAC-28","DAC-4","DAC-33","DAC-34","DAC-36","DAC-6","DAC-1","DAC-9","DAC-11","DAC-15","DAC-3","DAC-16","DAC-18","DAC-19","DAC-37","DAC-38")
resistance<-c("DAC-20","DAC-21","DAC-22","DAC-23","DAC-25","DAC-26","DAC-29","DAC-30","DAC-31","DAC-32","DAC-35","DAC-5","DAC-7","DAC-8","DAC-2","DAC-12","DAC-13","DAC-14","DAC-39")
all<-c(sensitivity,resistance)
exp1<-exp11[,match(all,colnames(exp11))]
exp2<-floor(exp1)
delete<-apply(exp2,1,function(x) mean(x==0))
cou<-which(delete>0.75)####在超过75%样本以上的都是0的基因删去
peak_count_name<-exp2[(-cou),]
geneid<-exp11[(-cou),1]



text_leave_one_out1<-matrix(0,37,2)
for(i in 1:37){
  cat(i,"\n")
  GEM_sensitive<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
                   "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
                   "DAC-19","DAC-37","DAC-39","DAC-38")
  GEM_res_S8495_sen<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
  GEM_res_S8495_res<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")
  all<-c(GEM_sensitive,GEM_res_S8495_sen,GEM_res_S8495_res)
  text<-all[i]
  
  text_leave_one_out1[i,1]<-i
  text_leave_one_out1[i,2]<-text
  if(sum(GEM_sensitive==text)==1){
    text_label<-"GEM_sensitive"
    GEM_sensitive<-GEM_sensitive[-match(text,GEM_sensitive)]
  }
  if(sum(GEM_res_S8495_sen==text)==1){
    text_label<-"GEM_res_S8495_sen"
    GEM_res_S8495_sen<-GEM_res_S8495_sen[-match(text,GEM_res_S8495_sen)]
  }
  if(sum(GEM_res_S8495_res==text)==1){
    text_label<-"GEM_res_S8495_res"
    GEM_res_S8495_res<-GEM_res_S8495_res[-match(text,GEM_res_S8495_res)]
  }
  sensitivity<-GEM_sensitive
  resistance<-c(GEM_res_S8495_sen,GEM_res_S8495_res)
  resistance_peak<-peak_count_name[,match(resistance,colnames(peak_count_name))]
  sensitivity_peak<-peak_count_name[,match(sensitivity,colnames(peak_count_name))]
  #peak_name<-peak_count_name[,1]
  peak_name<-geneid
  colDate<-data.frame(row.names = c(as.vector(resistance),as.vector(sensitivity)),
                      condition=factor(c(rep("resistance",length(resistance)),rep("sensitivity",length(sensitivity))))
  )
  datexpr<-cbind(resistance_peak,sensitivity_peak)##前面的相对于后面的上下调
  counts <- apply(datexpr,2,as.numeric)   ###矩阵中必须是数值
  dds<-DESeqDataSetFromMatrix(countData = counts,colData = colDate,design = ~condition)
  dds<-DESeq(dds)##进行标准化分析
  res<-results(dds)##将结果输出
  res<-as.data.frame(res)
  res<-cbind(peak_name,res)  ##对数据增加一列
  colnames(res)<- c('peak_id',"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj" )##对列名进行修改
  
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier_RNAseq_leave_one/classifier',i,sep='_')
  if (!file.exists(outputPath)){
    dir.create(outputPath, recursive=TRUE)
  }
  write.table(res,paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_GEM_AUC_RNA.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  write.table(text,paste(outputPath,'/',"text.txt",sep=""),sep = '\t',col.names = F,row.names = F,quote = FALSE)##数据输出
  
  sensitivity1<-GEM_res_S8495_sen
  resistance1<-GEM_res_S8495_res
  resistance_peak1<-peak_count_name[,match(resistance1,colnames(peak_count_name))]
  sensitivity_peak1<-peak_count_name[,match(sensitivity1,colnames(peak_count_name))]
  #peak_name<-peak_count_name[,1]
  peak_name<-geneid
  colDate<-data.frame(row.names = c(as.vector(resistance1),as.vector(sensitivity1)),
                      condition=factor(c(rep("resistance1",length(resistance1)),rep("sensitivity1",length(sensitivity1))))
  )
  datexpr1<-cbind(resistance_peak1,sensitivity_peak1)##前面的相对于后面的上下调
  counts1 <- apply(datexpr1,2,as.numeric)   ###矩阵中必须是数值
  dds1<-DESeqDataSetFromMatrix(countData = counts1,colData = colDate,design = ~condition)
  dds1<-DESeq(dds1)##进行标准化分析
  res1<-results(dds1)##将结果输出
  res1<-as.data.frame(res1)
  res1<-cbind(peak_name,res1)  ##对数据增加一列
  colnames(res1)<- c('peak_id',"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj" )##对列名进行修改
  
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier_RNAseq_leave_one/classifier',i,sep='_')
  if (!file.exists(outputPath)){
    dir.create(outputPath, recursive=TRUE)
  }
  write.table(res1,paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_S8495_AUC_RNA.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
}

#######################  三分类-留一法验证 ############################
rm(list=ls())
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_RPKM<-read.table("peak_RPKM.txt",sep="\t",header=T) #没有全0的
colnames(peak_RPKM)<-gsub("\\.","-",colnames(peak_RPKM))

GEM_sensitive<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
                 "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
                 "DAC-19","DAC-37","DAC-39","DAC-38")
GEM_res_S8495_sen<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
GEM_res_S8495_res<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")

svm_result<-NULL
library(doParallel)
library(parallel)
library(iterators)

cl <- makeCluster(20)
registerDoParallel(cl)
foreach(i=1:18,.combine='rbind')%dopar%{
  cat(i,"\n")
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier1/classifier',i,sep='_')
  GEM<-read.table(paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_GEM_AUC_ATAC.txt",sep=""),sep = '\t',header=T)##数据输出
  resSig<-GEM[which(GEM$pvalue<0.05 & abs(GEM$log2FoldChange)>2),]
  #resSig<-GEM[which(GEM$pvalue<0.01),]
  resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
  resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
  cat(sum(resSig$up_down=='up'),"\n")
  GEM_sensitive_DARs<-as.character(resSig[which(resSig$log2FoldChange>0),1])
  
  S8495_sen<-read.table(paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_GEM_res_S8495_sen_AUC_ATAC.txt",sep=""),sep = '\t',header=T)##数据输出
  resSig_S8495<-S8495_sen[which(S8495_sen$pvalue<0.05 & abs(S8495_sen$log2FoldChange)>2),]
  #resSig<-GEM[which(GEM$pvalue<0.01),]
  resSig_S8495[which(resSig_S8495$log2FoldChange>0),'up_down']<-'up'
  cat(sum(resSig_S8495$up_down=='up'),"\n")
  S8495_sensitive_DARs<-as.character(resSig_S8495[which(resSig_S8495$log2FoldChange>0),1])
  
  S8495_res<-read.table(paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_GEM_res_S8495_res_AUC_ATAC.txt",sep=""),sep = '\t',header=T)##数据输出
  resSig_S8495<-S8495_res[which(S8495_res$pvalue<0.05 & abs(S8495_res$log2FoldChange)>2),]
  #resSig<-GEM[which(GEM$pvalue<0.01),]
  resSig_S8495[which(resSig_S8495$log2FoldChange>0),'up_down']<-'up'
  cat(sum(resSig_S8495$up_down=='up'),"\n")
  S8495_resistant_DARs<-as.character(resSig_S8495[which(resSig_S8495$log2FoldChange>0),1])
  
  
  a<-intersect(GEM_sensitive_DARs,S8495_sensitive_DARs)
  if(length(a)==0){
    GEM_sensitive_only1<-GEM_sensitive_DARs
    S8495_sensitive_only<-S8495_sensitive_DARs
  }
  if(length(a)!=0){
    GEM_sensitive_only1<-GEM_sensitive_DARs[-match(a,GEM_sensitive_DARs)]
    S8495_sensitive_only<-S8495_sensitive_DARs[-match(a,S8495_sensitive_DARs)]
  }
  
  b<-intersect(GEM_sensitive_only1,S8495_resistant_DARs)
  if(length(b)==0){
    GEM_sensitive_only<-GEM_sensitive_only1
    S8495_resistant_only<-S8495_resistant_DARs
  }
  if(length(b)!=0){
    GEM_sensitive_only<-GEM_sensitive_only1[-match(b,GEM_sensitive_only1)]
    S8495_resistant_only<-S8495_resistant_DARs[-match(b,S8495_resistant_DARs)]
  }
  DARs<-c(GEM_sensitive_only,S8495_sensitive_only,S8495_resistant_only)
  all<-c(GEM_sensitive,GEM_res_S8495_sen,GEM_res_S8495_res)
  cat(length(DARs),"\n")
  peak_RPKM37<-peak_RPKM[,match(all,colnames(peak_RPKM))]
  peak_RPKM37_DARs<-peak_RPKM37[match(DARs,peak_RPKM[,1]),]
  peak_RPKM37_DARs1<-t(peak_RPKM37_DARs)
  colnames(peak_RPKM37_DARs1)<-DARs
  y<-matrix(c(rep("GEM_sensitive",24),rep("GEM_res_S8495_sen",7),rep("GEM_res_S8495_res",6)),ncol=1)
  colnames(y)<-"y"
  data<-as.data.frame(cbind(peak_RPKM37_DARs1,y))
  text<-read.table(paste(outputPath,'/',"text.txt",sep=""),sep = '\t',header=F)##数据输出
  library(e1071)
  set.seed(12345)
  text_data<-data[match(as.character(text[1,1]),all),]
  train_data<-data[-match(as.character(text[1,1]),all),]
  #SvmFit<-svm(y~.,data=train_data,type="C-classification",kernel="radial",cost=5,gamma=1,scale=FALSE)
  #head(SvmFit$decision.values)
  #yPred<-predict(SvmFit,text_data)
  #ConfM<-table(yPred,text_data$y)
  #Err<-(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  
  #### randomForest
  #library(randomForest)
  #rf_train<-randomForest(as.factor(y)~.,data=train_data,mtry=6,ntree=500,importance=TRUE,proximity=TRUE)
  #importance<-importance(rf_train)
  #yPred<- predict(rf_train,newdata=text_data)
  #ConfM<-table(yPred,text_data$y)
  #Err<-(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  
  ### NaiveBayes
  #library(klaR)
  #fit_Bayes1=NaiveBayes(y~.,data=train_data)
  #yPred=predict(fit_Bayes1,text_data)
  #ConfM<-table(text_data$y,yPred$class)
  #Err=sum(as.numeric(as.numeric(yPred$class)!=as.numeric(text_data$y)))/nrow(text_data) 
  
  library(caret)
  knn.model1 = knn3(y~.,data = train_data, k = 3) 
  yPred = predict(knn.model1,text_data,type = "class") 
  ConfM<-table(yPred,text_data$y)
  Err=(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  
  cat(Err,"\n")
  if(sum(GEM_sensitive==as.character(text[1,1]))==1){
    text_label<-"GEM_sensitive"
  }
  if(sum(GEM_res_S8495_sen==as.character(text[1,1]))==1){
    text_label<-"GEM_res_S8495_sen"
  }
  if(sum(GEM_res_S8495_res==as.character(text[1,1]))==1){
    text_label<-"GEM_res_S8495_res"
  }
  svm_result1<-as.data.frame(matrix(c(i,length(DARs),text_label,as.character(yPred),Err,1-Err),nrow=1))
  svm_result<-rbind(svm_result,svm_result1)
}
stopCluster(cl)
colnames(svm_result)<-c("times","length of biomarks","text_label","predict_label","error","precision")
sum(as.numeric(as.character(svm_result[,6])))/nrow(svm_result)

sum(as.numeric(as.character(svm_result[,6])))
setwd("~/xjj/drug/drug_result/HDACi_chemo618/classifier1")
write.table(svm_result,"rf_result_leave-one-out_new.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出


#######################  留一法验证 ############################
rm(list=ls())
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_RPKM<-read.table("peak_RPKM.txt",sep="\t",header=T) #没有全0的
colnames(peak_RPKM)<-gsub("\\.","-",colnames(peak_RPKM))

#GEM_sensitive<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
#                 "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
#                 "DAC-19","DAC-37","DAC-39","DAC-38")
#GEM_res_S8495_sen<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
#GEM_res_S8495_res<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")

GEM_sensitive<-c("DAC-24","DAC-27","DAC-28","DAC-4","DAC-33","DAC-34","DAC-36","DAC-6","DAC-1","DAC-9","DAC-11","DAC-15","DAC-3","DAC-16","DAC-18","DAC-19","DAC-37","DAC-38")
GEM_res_S8495_sen<-c("DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-12","DAC-14")
GEM_res_S8495_res<-c("DAC-20","DAC-21","DAC-22","DAC-23","DAC-25","DAC-29","DAC-30","DAC-32","DAC-35","DAC-5","DAC-13","DAC-39")


svm_result<-NULL
library(doParallel)
library(parallel)
library(iterators)

cl <- makeCluster(37)
registerDoParallel(cl)
foreach(i=1:37,.combine='rbind')%dopar%{
#for(i in 1:37){
  cat(i,"\n")
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier_median_leave_one/classifier',i,sep='_')
  GEM<-read.table(paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_GEM_AUC_ATAC.txt",sep=""),sep = '\t',header=T)##数据输出
  resSig<-GEM[which(GEM$pvalue<0.05 & abs(GEM$log2FoldChange)>2),]
  #resSig<-GEM[which(GEM$pvalue<0.01),]
  resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
  resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
  cat(sum(resSig$up_down=='up'),"\n")
  GEM_sensitive_DARs<-as.character(resSig[which(resSig$log2FoldChange>0),1])
  
  S8495<-read.table(paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_S8495_AUC_ATAC.txt",sep=""),sep = '\t',header=T)##数据输出
  resSig_S8495<-S8495[which(S8495$pvalue<0.05 & abs(S8495$log2FoldChange)>2),]
  #resSig<-GEM[which(GEM$pvalue<0.01),]
  resSig_S8495[which(resSig_S8495$log2FoldChange>0),'up_down']<-'up'
  resSig_S8495[which(resSig_S8495$log2FoldChange<0),'up_down']<-'down'
  cat(sum(resSig_S8495$up_down=='up'),"\n")
  cat(sum(resSig_S8495$up_down=='down'),"\n")
  S8495_sensitive_DARs<-as.character(resSig_S8495[which(resSig_S8495$log2FoldChange>0),1])
  S8495_resistant_DARs<-as.character(resSig_S8495[which(resSig_S8495$log2FoldChange<0),1])
  a<-intersect(GEM_sensitive_DARs,S8495_sensitive_DARs)
  if(length(a)==0){
    GEM_sensitive_only1<-GEM_sensitive_DARs
    S8495_sensitive_only<-S8495_sensitive_DARs
  }
  if(length(a)!=0){
    GEM_sensitive_only1<-GEM_sensitive_DARs[-match(a,GEM_sensitive_DARs)]
    S8495_sensitive_only<-S8495_sensitive_DARs[-match(a,S8495_sensitive_DARs)]
  }
  
  b<-intersect(GEM_sensitive_only1,S8495_resistant_DARs)
  if(length(b)==0){
    GEM_sensitive_only<-GEM_sensitive_only1
    S8495_resistant_only<-S8495_resistant_DARs
  }
  if(length(b)!=0){
    GEM_sensitive_only<-GEM_sensitive_only1[-match(b,GEM_sensitive_only1)]
    S8495_resistant_only<-S8495_resistant_DARs[-match(b,S8495_resistant_DARs)]
  }
  DARs<-c(GEM_sensitive_only,S8495_sensitive_only,S8495_resistant_only)
  all<-c(GEM_sensitive,GEM_res_S8495_sen,GEM_res_S8495_res)
  cat(length(DARs),"\n")
  peak_RPKM37<-peak_RPKM[,match(all,colnames(peak_RPKM))]
  peak_RPKM37_DARs<-peak_RPKM37[match(DARs,peak_RPKM[,1]),]
  peak_RPKM37_DARs1<-t(peak_RPKM37_DARs)
  colnames(peak_RPKM37_DARs1)<-DARs
  y<-matrix(c(rep("GEM_sensitive",18),rep("GEM_res_S8495_sen",7),rep("GEM_res_S8495_res",12)),ncol=1)
  colnames(y)<-"y"
  data<-as.data.frame(cbind(peak_RPKM37_DARs1,y))
  text<-read.table(paste(outputPath,'/',"text.txt",sep=""),sep = '\t',header=F)##数据输出
  library(e1071)
  set.seed(123)
  
  text_data<-data[match(as.character(text[1,1]),all),]
  train_data<-data[-match(as.character(text[1,1]),all),]
  #wts <- sum(train_data$y=="GEM_sensitive") / table(train_data$y)
  
  #SvmFit<-svm(y~.,data=train_data,type="C-classification",kernel="radial",cost=5,gamma=1,scale=FALSE)
  #head(SvmFit$decision.values)
  #yPred<-predict(SvmFit,text_data)
  #ConfM<-table(yPred,text_data$y)
  #Err<-(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  
  #### randomForest
  #library(randomForest)
  #rf_train<-randomForest(as.factor(y)~.,data=train_data,mtry=6,ntree=500,importance=TRUE,proximity=TRUE)
  #importance<-importance(rf_train)
  #yPred<- predict(rf_train,newdata=text_data)
  #ConfM<-table(yPred,text_data$y)
  #Err<-(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  
  ### NaiveBayes
  #library(klaR)
  #fit_Bayes1=NaiveBayes(y~.,data=train_data)
  #yPred=predict(fit_Bayes1,text_data)
  #ConfM<-table(text_data$y,yPred$class)
  #Err=sum(as.numeric(as.numeric(yPred$class)!=as.numeric(text_data$y)))/nrow(text_data) 
  
  library(caret)
  knn.model1 = knn3(y~.,data = train_data, k = 5) 
  yPred = predict(knn.model1,text_data,type = "class") 
  ConfM<-table(yPred,text_data$y)
  Err=(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  cat(Err,"\n")
  if(sum(GEM_sensitive==as.character(text[1,1]))==1){
    text_label<-"GEM_sensitive"
  }
  if(sum(GEM_res_S8495_sen==as.character(text[1,1]))==1){
    text_label<-"GEM_res_S8495_sen"
  }
  if(sum(GEM_res_S8495_res==as.character(text[1,1]))==1){
    text_label<-"GEM_res_S8495_res"
  }
  svm_result1<-as.data.frame(matrix(c(i,length(DARs),text_label,as.character(yPred),Err,1-Err),nrow=1))
  svm_result<-rbind(svm_result,svm_result1)
}
colnames(svm_result)<-c("times","length of biomarks","text_label","predict_label","error","precision")
sum(as.numeric(as.character(svm_result[,6])))/nrow(svm_result)

sum(as.numeric(as.character(svm_result[,6])))
setwd("~/xjj/drug/drug_result/HDACi_chemo618/classifier_median_leave_one")
write.table(svm_result,"KNN_result_leave-one-out.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
stopCluster(cl)

################  留一法验证  p值最显著特定数量
rm(list=ls())
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_RPKM<-read.table("peak_RPKM.txt",sep="\t",header=T) #没有全0的
colnames(peak_RPKM)<-gsub("\\.","-",colnames(peak_RPKM))

#GEM_sensitive<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
#                 "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
#                 "DAC-19","DAC-37","DAC-39","DAC-38")
#GEM_res_S8495_sen<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
#GEM_res_S8495_res<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")

GEM_sensitive<-c("DAC-24","DAC-27","DAC-28","DAC-4","DAC-33","DAC-34","DAC-36","DAC-6","DAC-1","DAC-9","DAC-11","DAC-15","DAC-3","DAC-16","DAC-18","DAC-19","DAC-37","DAC-38")
GEM_res_S8495_sen<-c("DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-12","DAC-14")
GEM_res_S8495_res<-c("DAC-20","DAC-21","DAC-22","DAC-23","DAC-25","DAC-29","DAC-30","DAC-32","DAC-35","DAC-5","DAC-13","DAC-39")


svm_result<-NULL
library(doParallel)
library(parallel)
library(iterators)

cl <- makeCluster(37)
registerDoParallel(cl)
foreach(i=1:37,.combine='rbind')%dopar%{
#for(i in 1:37){
  cat(i,"\n")
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier_median_leave_one/classifier',i,sep='_')
  GEM<-read.table(paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_GEM_AUC_ATAC.txt",sep=""),sep = '\t',header=T)##数据输出
  resSig<-GEM[order(GEM$pvalue),]
  resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
  resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
  GEM_sensitive_DARs<-as.character(resSig[which(resSig$up_down=="up")[1:100],1])
  
  
  S8495<-read.table(paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_S8495_AUC_ATAC.txt",sep=""),sep = '\t',header=T)##数据输出
  resSig_S8495<-S8495[order(S8495$pvalue),]
  resSig_S8495[which(resSig_S8495$log2FoldChange>0),'up_down']<-'up'
  resSig_S8495[which(resSig_S8495$log2FoldChange<0),'up_down']<-'down'
  S8495_sensitive_DARs<-as.character(resSig_S8495[which(resSig_S8495$up_down=="up")[1:100],1])
  S8495_resistant_DARs<-as.character(resSig_S8495[which(resSig_S8495$up_down=="down")[1:100],1])
  a<-intersect(GEM_sensitive_DARs,S8495_sensitive_DARs)
  if(length(a)==0){
    GEM_sensitive_only1<-GEM_sensitive_DARs
    S8495_sensitive_only<-S8495_sensitive_DARs
  }
  if(length(a)!=0){
    GEM_sensitive_only1<-GEM_sensitive_DARs[-match(a,GEM_sensitive_DARs)]
    S8495_sensitive_only<-S8495_sensitive_DARs[-match(a,S8495_sensitive_DARs)]
  }
  
  b<-intersect(GEM_sensitive_only1,S8495_resistant_DARs)
  if(length(b)==0){
    GEM_sensitive_only<-GEM_sensitive_only1
    S8495_resistant_only<-S8495_resistant_DARs
  }
  if(length(b)!=0){
    GEM_sensitive_only<-GEM_sensitive_only1[-match(b,GEM_sensitive_only1)]
    S8495_resistant_only<-S8495_resistant_DARs[-match(b,S8495_resistant_DARs)]
  }
  DARs<-c(GEM_sensitive_only,S8495_sensitive_only,S8495_resistant_only)
  all<-c(GEM_sensitive,GEM_res_S8495_sen,GEM_res_S8495_res)
  cat(length(DARs),"\n")
  peak_RPKM37<-peak_RPKM[,match(all,colnames(peak_RPKM))]
  peak_RPKM37_DARs<-peak_RPKM37[match(DARs,peak_RPKM[,1]),]
  peak_RPKM37_DARs1<-t(peak_RPKM37_DARs)
  colnames(peak_RPKM37_DARs1)<-DARs
  y<-matrix(c(rep("GEM_sensitive",18),rep("GEM_res_S8495_sen",7),rep("GEM_res_S8495_res",12)),ncol=1)
  colnames(y)<-"y"
  data<-as.data.frame(cbind(peak_RPKM37_DARs1,y))
  text<-read.table(paste(outputPath,'/',"text.txt",sep=""),sep = '\t',header=F)##数据输出
  library(e1071)
  set.seed(123)
  
  text_data<-data[match(as.character(text[1,1]),all),]
  train_data<-data[-match(as.character(text[1,1]),all),]
  #wts <- sum(train_data$y=="GEM_sensitive") / table(train_data$y)
  
  #SvmFit<-svm(y~.,data=train_data,type="C-classification",kernel="radial",cost=5,gamma=1,scale=FALSE)
  #head(SvmFit$decision.values)
  #yPred<-predict(SvmFit,text_data)
  #ConfM<-table(yPred,text_data$y)
  #Err<-(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  
  #### randomForest
  #library(randomForest)
  #rf_train<-randomForest(as.factor(y)~.,data=train_data,mtry=6,ntree=500,importance=TRUE,proximity=TRUE)
  #importance<-importance(rf_train)
  #yPred<- predict(rf_train,newdata=text_data)
  #ConfM<-table(yPred,text_data$y)
  #Err<-(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  
  ### NaiveBayes
  #library(klaR)
  #fit_Bayes1=NaiveBayes(y~.,data=train_data)
  #yPred=predict(fit_Bayes1,text_data)
  #ConfM<-table(text_data$y,yPred$class)
  #Err=sum(as.numeric(as.numeric(yPred$class)!=as.numeric(text_data$y)))/nrow(text_data) 
  
  library(caret)
  knn.model1 = knn3(y~.,data = train_data, k = 5) 
  yPred = predict(knn.model1,text_data,type = "class") 
  ConfM<-table(yPred,text_data$y)
  Err=(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  cat(Err,"\n")
  if(sum(GEM_sensitive==as.character(text[1,1]))==1){
    text_label<-"GEM_sensitive"
  }
  if(sum(GEM_res_S8495_sen==as.character(text[1,1]))==1){
    text_label<-"GEM_res_S8495_sen"
  }
  if(sum(GEM_res_S8495_res==as.character(text[1,1]))==1){
    text_label<-"GEM_res_S8495_res"
  }
  svm_result1<-as.data.frame(matrix(c(i,length(DARs),text_label,as.character(yPred),Err,1-Err),nrow=1))
  svm_result<-rbind(svm_result,svm_result1)
}
colnames(svm_result)<-c("times","length of biomarks","text_label","predict_label","error","precision")
sum(as.numeric(as.character(svm_result[,6])))/nrow(svm_result)

sum(as.numeric(as.character(svm_result[,6])))
setwd("~/xjj/drug/drug_result/HDACi_chemo618/classifier_median_leave_one")
write.table(svm_result,"SVM_100result_leave-one-out.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
stopCluster(cl)



#######################  留一法验证 RNAseq ############################
rm(list=ls())
setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
peak_RPKM<-read.table("star_rsem.GeneSymbol.FPKM.xls",sep="\t",header=T) #没有全0的
colnames(peak_RPKM)<-gsub("\\.","-",colnames(peak_RPKM))

#GEM_sensitive<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
#                 "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
#                 "DAC-19","DAC-37","DAC-39","DAC-38")
#GEM_res_S8495_sen<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
#GEM_res_S8495_res<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")

GEM_sensitive<-c("DAC-24","DAC-27","DAC-28","DAC-4","DAC-33","DAC-34","DAC-36","DAC-6","DAC-1","DAC-9","DAC-11","DAC-15","DAC-3","DAC-16","DAC-18","DAC-19","DAC-37","DAC-38")
GEM_res_S8495_sen<-c("DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-12","DAC-14")
GEM_res_S8495_res<-c("DAC-20","DAC-21","DAC-22","DAC-23","DAC-25","DAC-29","DAC-30","DAC-32","DAC-35","DAC-5","DAC-13","DAC-39")


svm_result<-NULL
for(i in 1:37){
  cat(i,"\n")
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier_median_RNAseq_leave_one/classifier',i,sep='_')
  GEM<-read.table(paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_GEM_AUC_RNA.txt",sep=""),sep = '\t',header=T)##数据输出
  resSig<-GEM[which(GEM$pvalue<0.05 & abs(GEM$log2FoldChange)>2),]
  #resSig<-GEM[which(GEM$pvalue<0.01),]
  resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
  resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
  cat(sum(resSig$up_down=='up'),"\n")
  GEM_sensitive_DARs<-as.character(resSig[which(resSig$log2FoldChange>0),1])
  
  S8495<-read.table(paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_S8495_AUC_RNA.txt",sep=""),sep = '\t',header=T)##数据输出
  resSig_S8495<-S8495[which(S8495$pvalue<0.05 & abs(S8495$log2FoldChange)>2),]
  #resSig<-GEM[which(GEM$pvalue<0.01),]
  resSig_S8495[which(resSig_S8495$log2FoldChange>0),'up_down']<-'up'
  resSig_S8495[which(resSig_S8495$log2FoldChange<0),'up_down']<-'down'
  cat(sum(resSig_S8495$up_down=='up'),"\n")
  cat(sum(resSig_S8495$up_down=='down'),"\n")
  S8495_sensitive_DARs<-as.character(resSig_S8495[which(resSig_S8495$log2FoldChange>0),1])
  S8495_resistant_DARs<-as.character(resSig_S8495[which(resSig_S8495$log2FoldChange<0),1])
  a<-intersect(GEM_sensitive_DARs,S8495_sensitive_DARs)
  if(length(a)==0){
    GEM_sensitive_only1<-GEM_sensitive_DARs
    S8495_sensitive_only<-S8495_sensitive_DARs
  }
  if(length(a)!=0){
    GEM_sensitive_only1<-GEM_sensitive_DARs[-match(a,GEM_sensitive_DARs)]
    S8495_sensitive_only<-S8495_sensitive_DARs[-match(a,S8495_sensitive_DARs)]
  }
  
  b<-intersect(GEM_sensitive_only1,S8495_resistant_DARs)
  if(length(b)==0){
    GEM_sensitive_only<-GEM_sensitive_only1
    S8495_resistant_only<-S8495_resistant_DARs
  }
  if(length(b)!=0){
    GEM_sensitive_only<-GEM_sensitive_only1[-match(b,GEM_sensitive_only1)]
    S8495_resistant_only<-S8495_resistant_DARs[-match(b,S8495_resistant_DARs)]
  }
  DARs<-c(GEM_sensitive_only,S8495_sensitive_only,S8495_resistant_only)
  all<-c(GEM_sensitive,GEM_res_S8495_sen,GEM_res_S8495_res)
  cat(length(DARs),"\n")
  peak_RPKM37<-peak_RPKM[,match(all,colnames(peak_RPKM))]
  peak_RPKM37_DARs<-peak_RPKM37[match(DARs,peak_RPKM[,1]),]
  peak_RPKM37_DARs1<-t(peak_RPKM37_DARs)
  colnames(peak_RPKM37_DARs1)<-DARs
  y<-matrix(c(rep("GEM_sensitive",18),rep("GEM_res_S8495_sen",7),rep("GEM_res_S8495_res",12)),ncol=1)
  colnames(y)<-"y"
  data<-as.data.frame(cbind(peak_RPKM37_DARs1,y))
  count_DARs<-c()
  for(m in 1:length(DARs)){
    count_DARs1<-length(unlist(strsplit(DARs[m],"-")))
    count_DARs<-c(count_DARs,count_DARs1)  
  }
  data<-data[,-(which(count_DARs==2))]
  text<-read.table(paste(outputPath,'/',"text.txt",sep=""),sep = '\t',header=F)##数据输出
  text_data<-data[match(as.character(text[1,1]),all),]
  train_data<-data[-match(as.character(text[1,1]),all),]
  wts <- sum(train_data$y=="GEM_sensitive") / table(train_data$y)
  
  #SvmFit<-svm(y~.,data=train_data,type="C-classification",kernel="radial",cost=5,gamma=1,scale=FALSE)
  #SvmFit<-svm(y~.,data=train_data,type="C-classification",kernel="radial",cost=5,gamma=1,scale=FALSE)
  #head(SvmFit$decision.values)
  #yPred<-predict(SvmFit,text_data)
  #ConfM<-table(yPred,text_data$y)
  #Err<-(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  
  #### randomForest
  #library(randomForest)
  #rf_train<-randomForest(as.factor(y)~.,data=train_data,mtry=6,ntree=500,importance=TRUE,proximity=TRUE)
  #importance<-importance(rf_train)
  #yPred<- predict(rf_train,newdata=text_data)
  #ConfM<-table(yPred,text_data$y)
  #Err<-(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  
  ### NaiveBayes
  #library(klaR)
  #fit_Bayes1=NaiveBayes(y~.,data=train_data)
  #yPred=predict(fit_Bayes1,text_data)
  #ConfM<-table(text_data$y,yPred$class)
  #Err=sum(as.numeric(as.numeric(yPred$class)!=as.numeric(text_data$y)))/nrow(text_data) 
  
  library(caret)
  knn.model1 = knn3(y~.,data = train_data, k = 3) 
  yPred = predict(knn.model1,text_data,type = "class") 
  ConfM<-table(yPred,text_data$y)
  Err=(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  cat(Err,"\n")
  if(sum(GEM_sensitive==as.character(text[1,1]))==1){
    text_label<-"GEM_sensitive"
  }
  if(sum(GEM_res_S8495_sen==as.character(text[1,1]))==1){
    text_label<-"GEM_res_S8495_sen"
  }
  if(sum(GEM_res_S8495_res==as.character(text[1,1]))==1){
    text_label<-"GEM_res_S8495_res"
  }
  svm_result1<-as.data.frame(matrix(c(i,length(DARs),text_label,as.character(yPred),Err,1-Err),nrow=1))
  svm_result<-rbind(svm_result,svm_result1)
}
colnames(svm_result)<-c("times","length of biomarks","text_label","predict_label","error","precision")
sum(as.numeric(as.character(svm_result[,6])))/nrow(svm_result)

sum(as.numeric(as.character(svm_result[,6])))
setwd("~/xjj/drug/drug_result/HDACi_chemo618/classifier_median_RNAseq_leave_one")
write.table(svm_result,"KNN3_result_leave-one-out.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
stopCluster(cl)

#######################  留一法验证 RNAseq 前100个基因############################
rm(list=ls())
setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
peak_RPKM<-read.table("star_rsem.GeneSymbol.FPKM.xls",sep="\t",header=T) #没有全0的
colnames(peak_RPKM)<-gsub("\\.","-",colnames(peak_RPKM))

#GEM_sensitive<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
#                 "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
#                 "DAC-19","DAC-37","DAC-39","DAC-38")
#GEM_res_S8495_sen<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
#GEM_res_S8495_res<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")

GEM_sensitive<-c("DAC-24","DAC-27","DAC-28","DAC-4","DAC-33","DAC-34","DAC-36","DAC-6","DAC-1","DAC-9","DAC-11","DAC-15","DAC-3","DAC-16","DAC-18","DAC-19","DAC-37","DAC-38")
GEM_res_S8495_sen<-c("DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-12","DAC-14")
GEM_res_S8495_res<-c("DAC-20","DAC-21","DAC-22","DAC-23","DAC-25","DAC-29","DAC-30","DAC-32","DAC-35","DAC-5","DAC-13","DAC-39")


svm_result<-NULL
for(i in 1:37){
  cat(i,"\n")
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier_median_RNAseq_leave_one/classifier',i,sep='_')
  GEM<-read.table(paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_GEM_AUC_RNA.txt",sep=""),sep = '\t',header=T)##数据输出
  resSig<-GEM[order(GEM$pvalue),]
  resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
  resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
  GEM_sensitive_DARs<-as.character(resSig[which(resSig$up_down=="up")[1:100],1])
  
  
  S8495<-read.table(paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_S8495_AUC_RNA.txt",sep=""),sep = '\t',header=T)##数据输出
  resSig_S8495<-S8495[order(S8495$pvalue),]
  resSig_S8495[which(resSig_S8495$log2FoldChange>0),'up_down']<-'up'
  resSig_S8495[which(resSig_S8495$log2FoldChange<0),'up_down']<-'down'
  S8495_sensitive_DARs<-as.character(resSig_S8495[which(resSig_S8495$up_down=="up")[1:100],1])
  S8495_resistant_DARs<-as.character(resSig_S8495[which(resSig_S8495$up_down=="down")[1:100],1])
  a<-intersect(GEM_sensitive_DARs,S8495_sensitive_DARs)
  if(length(a)==0){
    GEM_sensitive_only1<-GEM_sensitive_DARs
    S8495_sensitive_only<-S8495_sensitive_DARs
  }
  if(length(a)!=0){
    GEM_sensitive_only1<-GEM_sensitive_DARs[-match(a,GEM_sensitive_DARs)]
    S8495_sensitive_only<-S8495_sensitive_DARs[-match(a,S8495_sensitive_DARs)]
  }
  
  b<-intersect(GEM_sensitive_only1,S8495_resistant_DARs)
  if(length(b)==0){
    GEM_sensitive_only<-GEM_sensitive_only1
    S8495_resistant_only<-S8495_resistant_DARs
  }
  if(length(b)!=0){
    GEM_sensitive_only<-GEM_sensitive_only1[-match(b,GEM_sensitive_only1)]
    S8495_resistant_only<-S8495_resistant_DARs[-match(b,S8495_resistant_DARs)]
  }
  DARs<-c(GEM_sensitive_only,S8495_sensitive_only,S8495_resistant_only)
  all<-c(GEM_sensitive,GEM_res_S8495_sen,GEM_res_S8495_res)
  cat(length(DARs),"\n")
  peak_RPKM37<-peak_RPKM[,match(all,colnames(peak_RPKM))]
  peak_RPKM37_DARs<-peak_RPKM37[match(DARs,peak_RPKM[,1]),]
  peak_RPKM37_DARs1<-t(peak_RPKM37_DARs)
  colnames(peak_RPKM37_DARs1)<-DARs
  y<-matrix(c(rep("GEM_sensitive",18),rep("GEM_res_S8495_sen",7),rep("GEM_res_S8495_res",12)),ncol=1)
  colnames(y)<-"y"
  data<-as.data.frame(cbind(peak_RPKM37_DARs1,y))
  count_DARs<-c()
  for(m in 1:length(DARs)){
    count_DARs1<-length(unlist(strsplit(DARs[m],"-")))
    count_DARs<-c(count_DARs,count_DARs1)  
  }
  if(sum(count_DARs==2)>0){
    data<-data[,-(which(count_DARs==2))]
  }
  if(sum(count_DARs==2)==0){
    data<-data
  }
  
  count_DARs2<-c()
  for(m in 1:length(colnames(data))){
    count_DARs1<-length(unlist(strsplit(colnames(data)[m],"/")))
    count_DARs2<-c(count_DARs2,count_DARs1)  
  }
  if(sum(count_DARs2==2)>0){
    data<-data[,-(which(count_DARs2==2))]
  }
  if(sum(count_DARs2==2)==0){
    data<-data
  }
  
  text<-read.table(paste(outputPath,'/',"text.txt",sep=""),sep = '\t',header=F)##数据输出
  library(e1071)
  set.seed(123)
  
  text_data<-data[match(as.character(text[1,1]),all),]
  train_data<-data[-match(as.character(text[1,1]),all),]
  #wts <- sum(train_data$y=="GEM_sensitive") / table(train_data$y)
  
  #SvmFit<-svm(y~.,data=train_data,type="C-classification",kernel="radial",cost=5,gamma=1,scale=FALSE)
  #head(SvmFit$decision.values)
  #yPred<-predict(SvmFit,text_data)
  #ConfM<-table(yPred,text_data$y)
  #Err<-(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  
  #### randomForest
  #library(randomForest)
  #rf_train<-randomForest(as.factor(y)~.,data=train_data,mtry=6,ntree=500,importance=TRUE,proximity=TRUE)
  #importance<-importance(rf_train)
  #yPred<- predict(rf_train,newdata=text_data)
  #ConfM<-table(yPred,text_data$y)
  #Err<-(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  
  ### NaiveBayes
  library(klaR)
  fit_Bayes1=NaiveBayes(y~.,data=train_data)
  yPred=predict(fit_Bayes1,text_data)
  ConfM<-table(text_data$y,yPred$class)
  Err=sum(as.numeric(as.numeric(yPred$class)!=as.numeric(text_data$y)))/nrow(text_data) 
  
  #library(caret)
  #knn.model1 = knn3(y~.,data = train_data, k = 3) 
  #yPred = predict(knn.model1,text_data,type = "class") 
  #ConfM<-table(yPred,text_data$y)
  #Err=(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  cat(Err,"\n")
  if(sum(GEM_sensitive==as.character(text[1,1]))==1){
    text_label<-"GEM_sensitive"
  }
  if(sum(GEM_res_S8495_sen==as.character(text[1,1]))==1){
    text_label<-"GEM_res_S8495_sen"
  }
  if(sum(GEM_res_S8495_res==as.character(text[1,1]))==1){
    text_label<-"GEM_res_S8495_res"
  }
  svm_result1<-as.data.frame(matrix(c(i,length(DARs),text_label,as.character(yPred),Err,1-Err),nrow=1))
  svm_result<-rbind(svm_result,svm_result1)
}
colnames(svm_result)<-c("times","length of biomarks","text_label","predict_label","error","precision")
sum(as.numeric(as.character(svm_result[,6])))/nrow(svm_result)

sum(as.numeric(as.character(svm_result[,6])))
setwd("~/xjj/drug/drug_result/HDACi_chemo618/classifier_median_RNAseq_leave_one")
write.table(svm_result,"NaiveBayes_result_leave-one-out100.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
stopCluster(cl)



#######################  留一法验证 先找出biomarkers############################
rm(list=ls())
GEM_sensitive<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
                 "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
                 "DAC-19","DAC-37","DAC-39","DAC-38")
GEM_res_S8495_sen<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
GEM_res_S8495_res<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")

DARs_result<-c()
for(i in 1:37){
  cat(i,"\n")
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier1/classifier',i,sep='_')
  GEM<-read.table(paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_GEM_AUC_ATAC.txt",sep=""),sep = '\t',header=T)##数据输出
  resSig<-GEM[which(GEM$pvalue<0.05 & abs(GEM$log2FoldChange)>2),]
  #resSig<-GEM[which(GEM$pvalue<0.01),]
  resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
  resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
  cat(sum(resSig$up_down=='up'),"\n")
  GEM_sensitive_DARs<-as.character(resSig[which(resSig$log2FoldChange>0),1])
  
  S8495<-read.table(paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_S8495_AUC_ATAC.txt",sep=""),sep = '\t',header=T)##数据输出
  resSig_S8495<-S8495[which(S8495$pvalue<0.05 & abs(S8495$log2FoldChange)>2),]
  #resSig<-GEM[which(GEM$pvalue<0.01),]
  resSig_S8495[which(resSig_S8495$log2FoldChange>0),'up_down']<-'up'
  resSig_S8495[which(resSig_S8495$log2FoldChange<0),'up_down']<-'down'
  cat(sum(resSig_S8495$up_down=='up'),"\n")
  cat(sum(resSig_S8495$up_down=='down'),"\n")
  S8495_sensitive_DARs<-as.character(resSig_S8495[which(resSig_S8495$log2FoldChange>0),1])
  S8495_resistant_DARs<-as.character(resSig_S8495[which(resSig_S8495$log2FoldChange<0),1])
  a<-intersect(GEM_sensitive_DARs,S8495_sensitive_DARs)
  if(length(a)==0){
    GEM_sensitive_only1<-GEM_sensitive_DARs
    S8495_sensitive_only<-S8495_sensitive_DARs
  }
  if(length(a)!=0){
    GEM_sensitive_only1<-GEM_sensitive_DARs[-match(a,GEM_sensitive_DARs)]
    S8495_sensitive_only<-S8495_sensitive_DARs[-match(a,S8495_sensitive_DARs)]
  }
  
  b<-intersect(GEM_sensitive_only1,S8495_resistant_DARs)
  if(length(b)==0){
    GEM_sensitive_only<-GEM_sensitive_only1
    S8495_resistant_only<-S8495_resistant_DARs
  }
  if(length(b)!=0){
    GEM_sensitive_only<-GEM_sensitive_only1[-match(b,GEM_sensitive_only1)]
    S8495_resistant_only<-S8495_resistant_DARs[-match(b,S8495_resistant_DARs)]
  }
  DARs<-c(GEM_sensitive_only,S8495_sensitive_only,S8495_resistant_only)
  DARs_result<-c(DARs_result,DARs)
}
setwd("~/xjj/drug/drug_result/HDACi_chemo618/classifier1/biomarker")
write.table(DARs_result,file="all_peaks_P0.05_FC2.txt",sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')

uni_peak<-unique(DARs_result)
peaks_result<-matrix(0,length(uni_peak),2)
for(i in 1:length(uni_peak)){
  peaks_result[i,1]<-uni_peak[i]
  peaks_result[i,2]<-sum(DARs_result==uni_peak[i])
}
biomarker<-peaks_result[which(peaks_result[,2]==37),1]
setwd("~/xjj/drug/drug_result/HDACi_chemo618/classifier1/biomarker")
write.table(biomarker,file="biomarker_P0.05_FC2.txt",sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')
write.table(peaks_result,file="peaks_counts_P0.05_FC2.txt",sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')


##############  使用biomaker进行分类
rm(list=ls())
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_RPKM<-read.table("peak_RPKM.txt",sep="\t",header=T) #没有全0的
colnames(peak_RPKM)<-gsub("\\.","-",colnames(peak_RPKM))

GEM_sensitive<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
                 "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
                 "DAC-19","DAC-37","DAC-39","DAC-38")
GEM_res_S8495_sen<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
GEM_res_S8495_res<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")

#library(doParallel)
#library(parallel)
#library(iterators)
#cl <- makeCluster(37)
#registerDoParallel(cl)
svm_result<-NULL
for(i in 1:37){
  cat(i,"\n")
  setwd("~/xjj/drug/drug_result/HDACi_chemo618/classifier1/biomarker")
  DARs1<-read.table("biomarker_P0.05_FC2.txt",sep="\t",header=T)
  DARs<-as.character(DARs1[,1])
  all<-c(GEM_sensitive,GEM_res_S8495_sen,GEM_res_S8495_res)
  cat(length(DARs),"\n")
  peak_RPKM37<-peak_RPKM[,match(all,colnames(peak_RPKM))]
  peak_RPKM37_DARs<-peak_RPKM37[match(DARs,peak_RPKM[,1]),]
  peak_RPKM37_DARs1<-t(peak_RPKM37_DARs)
  colnames(peak_RPKM37_DARs1)<-DARs
  y<-matrix(c(rep("GEM_sensitive",24),rep("GEM_res_S8495_sen",7),rep("GEM_res_S8495_res",6)),ncol=1)
  colnames(y)<-"y"
  data<-as.data.frame(cbind(peak_RPKM37_DARs1,y))
  
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier1/classifier',i,sep='_')
  text<-read.table(paste(outputPath,'/',"text.txt",sep=""),sep = '\t',header=F)##数据输出
  
  text_data<-data[match(as.character(text[1,1]),all),]
  train_data<-data[-match(as.character(text[1,1]),all),]
  #wts <- sum(train_data$y=="GEM_sensitive") / table(train_data$y)
  
  #SvmFit<-svm(y~.,data=train_data,type="C-classification",kernel="radial",cost=5,gamma=1,scale=FALSE,class.weights=wts)
  library(e1071)
  set.seed(123)
  #SvmFit<-svm(y~.,data=train_data,type="C-classification",kernel="radial",cost=5,gamma=1,scale=FALSE)
  #head(SvmFit$decision.values)
  #yPred<-predict(SvmFit,text_data)
  #ConfM<-table(yPred,text_data$y)
  #Err<-(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  
  #### randomForest
  #library(randomForest)
  #rf_train<-randomForest(as.factor(y)~.,data=train_data,mtry=6,ntree=500,importance=TRUE,proximity=TRUE)
  #importance<-importance(rf_train)
  #yPred<- predict(rf_train,newdata=text_data)
  #ConfM<-table(yPred,text_data$y)
  #Err<-(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  
  ### NaiveBayes
  #library(klaR)
  #fit_Bayes1=NaiveBayes(y~.,data=train_data)
  #yPred=predict(fit_Bayes1,text_data)
  #ConfM<-table(text_data$y,yPred$class)
  #Err=sum(as.numeric(as.numeric(yPred$class)!=as.numeric(text_data$y)))/nrow(text_data) 
  
  library(caret)
  knn.model1 = knn3(y~.,data = train_data, k = 5) 
  yPred = predict(knn.model1,text_data,type = "class") 
  ConfM<-table(yPred,text_data$y)
  Err=(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  cat(Err,"\n")
  if(sum(GEM_sensitive==as.character(text[1,1]))==1){
    text_label<-"GEM_sensitive"
  }
  if(sum(GEM_res_S8495_sen==as.character(text[1,1]))==1){
    text_label<-"GEM_res_S8495_sen"
  }
  if(sum(GEM_res_S8495_res==as.character(text[1,1]))==1){
    text_label<-"GEM_res_S8495_res"
  }
  #svm_result1<-as.data.frame(matrix(c(i,length(DARs),text_label,as.character(yPred$class),Err,1-Err),nrow=1))
  svm_result1<-as.data.frame(matrix(c(i,length(DARs),text_label,as.character(yPred),Err,1-Err),nrow=1))
  svm_result<-rbind(svm_result,svm_result1)
}
colnames(svm_result)<-c("times","length of biomarks","text_label","predict_label","error","precision")
sum(as.numeric(as.character(svm_result[,6])))/nrow(svm_result)

sum(as.numeric(as.character(svm_result[,6])))
setwd("~/xjj/drug/drug_result/HDACi_chemo618/classifier1")
setwd("~/xjj/drug/drug_result/HDACi_chemo618/classifier1/biomarker")
write.table(svm_result,"knn_result.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
stopCluster(cl)



########################  差异peaks做注释
for(i in 1:37){
  cat(i,"\n")
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier_median_leave_one/classifier',i,sep='_')
  GEM<-read.table(paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_GEM_AUC_ATAC.txt",sep=""),sep = '\t',header=T)##数据输出
  #resSig<-GEM[which(GEM$pvalue<0.05 & abs(GEM$log2FoldChange)>2),]
  resSig<-GEM[which(GEM$pvalue<0.01),]
  
  # pvalue padj
  ##新增一列，将log2FoldChange>0标注为up，<0标准为down
  resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
  resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
  up_gene<-as.data.frame(resSig[which(resSig$log2FoldChange>0),1])
  down_gene<-as.data.frame(resSig[which(resSig$log2FoldChange<0),1])
  library(dplyr)
  a_up<-apply(up_gene,1,function(x) unlist(strsplit(as.character(x), "_")))
  d_up<-t(a_up)
  e_up<-data.frame(rep("+",nrow(d_up)))
  peaks_up<-cbind(up_gene,d_up,e_up)
  colnames(peaks_up)<-c("peak_id","chr","start","end","strand")
  setwd(paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier_median_leave_one/classifier',i,sep='_'))
  out_file<-"sensitivity-resistance-GEM-AUC-P0.01-logFC0"
  write.table(peaks_up,file=paste(out_file,"-DESeq2-up.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  a_down<-apply(down_gene,1,function(x) unlist(strsplit(as.character(x), "_")))
  d_down<-t(a_down)
  e_down<-data.frame(rep("+",nrow(d_down)))
  peaks_down<-cbind(down_gene,d_down,e_down)
  colnames(peaks_down)<-c("peak_id","chr","start","end","strand")
  write.table(peaks_down,file=paste(out_file,"-DESeq2-down.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
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
  
  
  S8495<-read.table(paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_S8495_AUC_ATAC.txt",sep=""),sep = '\t',header=T)##数据输出
  resSig1<-S8495[which(S8495$pvalue<0.01),]
  # pvalue padj
  ##新增一列，将log2FoldChange>0标注为up，<0标准为down
  resSig1[which(resSig1$log2FoldChange>0),'up_down']<-'up'
  resSig1[which(resSig1$log2FoldChange<0),'up_down']<-'down'
  up_gene1<-as.data.frame(resSig1[which(resSig1$log2FoldChange>0),1])
  down_gene1<-as.data.frame(resSig1[which(resSig1$log2FoldChange<0),1])
  library(dplyr)
  a_up<-apply(up_gene1,1,function(x) unlist(strsplit(as.character(x), "_")))
  d_up<-t(a_up)
  e_up<-data.frame(rep("+",nrow(d_up)))
  peaks_up1<-cbind(up_gene1,d_up,e_up)
  colnames(peaks_up1)<-c("peak_id","chr","start","end","strand")
  setwd(paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier_median_leave_one/classifier',i,sep='_'))
  out_file1<-"sensitivity-resistance-S8495-AUC-P0.01-logFC0"
  write.table(peaks_up1,file=paste(out_file1,"-DESeq2-up.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  a_down<-apply(down_gene1,1,function(x) unlist(strsplit(as.character(x), "_")))
  d_down<-t(a_down)
  e_down<-data.frame(rep("+",nrow(d_down)))
  peaks_down1<-cbind(down_gene1,d_down,e_down)
  colnames(peaks_down1)<-c("peak_id","chr","start","end","strand")
  write.table(peaks_down1,file=paste(out_file1,"-DESeq2-down.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  #####  bed
  f_up<-data.frame(rep(1,nrow(d_up)))
  peaks_bed_up<-cbind(d_up,up_gene1,f_up,e_up)
  colnames(peaks_bed_up)<-c("chr","start","end","peak_id","score","strand")
  write.table(peaks_bed_up,file=paste(out_file1,"-DESeq2-up.bed",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  f_down<-data.frame(rep(1,nrow(d_down)))
  peaks_bed_down<-cbind(d_down,down_gene1,f_down,e_down)
  colnames(peaks_bed_down)<-c("chr","start","end","peak_id","score","strand")
  write.table(peaks_bed_down,file=paste(out_file1,"-DESeq2-down.bed",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
}

result1<-c()
for(i in 1:30){
  cat(i,"\n")
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier',i,sep='_')
  text<-read.table(paste(outputPath,'/',"text.txt",sep=""),sep = '\t',header=F)##数据输出
  result1<-c(result1,as.character(text[1,1]))
}  

match(unique(result1),result1)
unique(result1)

result2<-c("DAC-31","DAC-35","DAC-13","DAC-14","DAC-29","DAC-18","DAC-25","DAC-7","DAC-27","DAC-22","DAC-32","DAC-9","DAC-37","DAC-16","DAC-26","DAC-1","DAC-36","DAC-28","DAC-12","DAC-19","DAC-21")


########## DEGs ##########################################################
########## DEGs
rm(list=ls())
library(DESeq2)
library("data.table")
library(dplyr)

setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
library("data.table")
exp1<-fread("star_rsem.GeneSymbol.readscount.xls",header=T,data.table=F)
exp2<-floor(exp1[,-1])
delete<-apply(exp2,1,function(x) mean(x==0))
cou<-which(delete>0.75)####在超过75%样本以上的都是0的基因删去
exp<-exp2[(-cou),]
geneid<-exp1[(-cou),1]


text_leave_one_out3<-matrix(0,37,2)
for(i in 22:27){
  cat(i,"\n")
  GEM_sensitive<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
                   "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
                   "DAC-19","DAC-37","DAC-39","DAC-38")
  GEM_res_S8495_sen<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
  GEM_res_S8495_res<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")
  all<-c(GEM_sensitive,GEM_res_S8495_sen,GEM_res_S8495_res)
  
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier',i,sep='_')
  text1<-read.table(paste(outputPath,'/',"text.txt",sep=""),sep = '\t',header=F)##数据输出
  text<-as.character(text1[1,1])
  
  text_leave_one_out3[i,1]<-i
  text_leave_one_out3[i,2]<-text
  if(sum(GEM_sensitive==text)==1){
    text_label<-"GEM_sensitive"
    GEM_sensitive<-GEM_sensitive[-match(text,GEM_sensitive)]
  }
  if(sum(GEM_res_S8495_sen==text)==1){
    text_label<-"GEM_res_S8495_sen"
    GEM_res_S8495_sen<-GEM_res_S8495_sen[-match(text,GEM_res_S8495_sen)]
  }
  if(sum(GEM_res_S8495_res==text)==1){
    text_label<-"GEM_res_S8495_res"
    GEM_res_S8495_res<-GEM_res_S8495_res[-match(text,GEM_res_S8495_res)]
  }
  sensitivity<-GEM_sensitive
  resistance<-c(GEM_res_S8495_sen,GEM_res_S8495_res)
  resistance_exp<-exp[,match(resistance,colnames(exp))]
  sensitivity_exp<-exp[,match(sensitivity,colnames(exp))]
  colDate<-data.frame(row.names = c(as.vector(resistance),as.vector(sensitivity)),
                      condition=factor(c(rep("resistance",length(resistance)),rep("sensitivity",length(sensitivity))))
  )
  datexpr<-cbind(resistance_exp,sensitivity_exp)
  counts <- apply(datexpr,2,as.numeric)###矩阵中必须是数值
  dds<-DESeqDataSetFromMatrix(countData = counts,colData = colDate,design = ~condition)
  dds<-DESeq(dds)##进行标准化分析
  res<-results(dds)##将结果输出
  res<-as.data.frame(res)
  res<-cbind(geneid,res)  ##对数据增加一列
  colnames(res)<- c('gene_id',"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj" )##对列名进行修改
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier',i,sep='_')
  if (!file.exists(outputPath)){
    dir.create(outputPath, recursive=TRUE)
  }
  write.table(res,paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_GEM_AUC_RNAseq_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  resSig<-res[which(res$pvalue<0.05),]
  resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
  resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
  up_gene<-resSig[which(resSig$log2FoldChange>0),]
  down_gene<-resSig[which(resSig$log2FoldChange<0),]
  
  write.table(up_gene,paste(outputPath,'/',"sensitivity-resistance-GEM-DESeq2_AUC_RNAseq_P0.05_FC0_up_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  write.table(down_gene,paste(outputPath,'/',"sensitivity-resistance-GEM-DESeq2_AUC_RNAseq_P0.05_FC0_down_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  
  sensitivity1<-GEM_res_S8495_sen
  resistance1<-GEM_res_S8495_res
  resistance_exp1<-exp[,match(resistance1,colnames(exp))]
  sensitivity_exp1<-exp[,match(sensitivity1,colnames(exp))]
  colDate1<-data.frame(row.names = c(as.vector(resistance1),as.vector(sensitivity1)),
                      condition=factor(c(rep("resistance1",length(resistance1)),rep("sensitivity1",length(sensitivity1))))
  )
  datexpr1<-cbind(resistance_exp1,sensitivity_exp1)
  counts1 <- apply(datexpr1,2,as.numeric)###矩阵中必须是数值
  dds1<-DESeqDataSetFromMatrix(countData = counts1,colData = colDate1,design = ~condition)
  dds1<-DESeq(dds1)##进行标准化分析
  res1<-results(dds1)##将结果输出
  res1<-as.data.frame(res1)
  res1<-cbind(geneid,res1)  ##对数据增加一列
  colnames(res1)<- c('gene_id',"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj" )##对列名进行修改
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier',i,sep='_')
  if (!file.exists(outputPath)){
    dir.create(outputPath, recursive=TRUE)
  }
  write.table(res1,paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_S8495_AUC_RNAseq_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  resSig1<-res1[which(res1$pvalue<0.05),]
  resSig1[which(resSig1$log2FoldChange>0),'up_down']<-'up'
  resSig1[which(resSig1$log2FoldChange<0),'up_down']<-'down'
  up_gene1<-resSig1[which(resSig1$log2FoldChange>0),]
  down_gene1<-resSig1[which(resSig1$log2FoldChange<0),]
  if (!file.exists(outputPath)){
    dir.create(outputPath, recursive=TRUE)
  }
  write.table(up_gene1,paste(outputPath,'/',"sensitivity-resistance-S8495-DESeq2_AUC_RNAseq_P0.05_FC0_up_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  write.table(down_gene1,paste(outputPath,'/',"sensitivity-resistance-S8495-DESeq2_AUC_RNAseq_P0.05_FC0_down_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
}

################## RNAseq ATACseq联合分析筛选biomarker ######
rm(list=ls())
library(data.table)
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_RPKM<-read.table("peak_RPKM.txt",sep="\t",header=T) #没有全0的
colnames(peak_RPKM)<-gsub("\\.","-",colnames(peak_RPKM))

svm_result0<-NULL
library(doParallel)
library(parallel)
library(iterators)
library(data.table)
cl <- makeCluster(37)
registerDoParallel(cl)
#foreach(i=1:37,.combine='rbind')%dopar%{
for(k in 1:37){
  cat(k,"\n")
  library(data.table)
  setwd(paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier_median_leave_one/classifier',k,sep='_'))
  homer_anno_up<-fread("sensitivity-resistance-GEM-AUC-P0.01-logFC0-DESeq2-up-annotation.txt",header=T,data.table=F)
  colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
  Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
  length(unique(Anno_gene_100Kb_up$`Gene Name`))
  
  
  setwd(paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier_median_RNAseq_leave_one/classifier',k,sep='_'))
  df<-read.table("sensitivity-resistance-all-DESeq2_GEM_AUC_RNA.txt",header=T,sep="\t")
  up_gene<-df[which(df$pvalue<0.05&df$log2FoldChange>0),]
  
  DEGs<-as.character(up_gene[,1])
  int_all_DEGs<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)
  all_peaks_anno<-homer_anno_up
  setwd(paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier_median_leave_one/classifier',k,sep='_'))
  DApeak<-fread("sensitivity-resistance-all-DESeq2_GEM_AUC_ATAC.txt",header=T,data.table=F)
  
  
  different_FC<-NULL
  for(i in 1:length(int_all_DEGs)){
    DEGs_peaks1<-all_peaks_anno[all_peaks_anno$`Gene Name` %in% int_all_DEGs[i],1]
    DEGs_peaks_FC<-as.data.frame(DApeak[match(DEGs_peaks1,DApeak$peak_id),3])
    gene_name<-matrix(rep(int_all_DEGs[i],length(DEGs_peaks1)),ncol=1)
    DEGs_FC<-matrix(rep(up_gene[match(int_all_DEGs[i],DEGs),3],ncol=1))
    different_FC1<-cbind(gene_name,DEGs_FC,DEGs_peaks1,DEGs_peaks_FC)
    different_FC<-rbind(different_FC,different_FC1)
  }
  sum(different_FC[,2]>0) #up DEGs
  sum(different_FC[,4]>0) #up DARs
  colnames(different_FC)<-c("GEM_sensitive_DEGs","DEGs_FC","GEM_sensitive_DARs","DARs_FC")
  GEM_sensitive_DARs<-as.character(different_FC[which(different_FC[,2]>0 & different_FC[,4]>0),3])
  
  setwd(paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier_median_leave_one/classifier',k,sep='_'))
  homer_anno_up<-fread("sensitivity-resistance-S8495-AUC-P0.01-logFC0-DESeq2-down-annotation.txt",header=T,data.table=F)
  colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
  Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
  length(unique(Anno_gene_100Kb_up$`Gene Name`))
  all_peaks_anno<-homer_anno_up
  DApeak<-fread("sensitivity-resistance-all-DESeq2_S8495_AUC_ATAC.txt",header=T,data.table=F)
  
  setwd(paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier_median_RNAseq_leave_one/classifier',k,sep='_'))
  df1<-read.table("sensitivity-resistance-all-DESeq2_S8495_AUC_RNA.txt",header=T,sep="\t")
  down_gene<-df1[which(df1$pvalue<0.05&df1$log2FoldChange<0),]
  DEGs<-as.character(down_gene[,1])
  int_all_DEGs<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)
  
  
  different_FC<-NULL
  for(i in 1:length(int_all_DEGs)){
    DEGs_peaks1<-all_peaks_anno[all_peaks_anno$`Gene Name` %in% int_all_DEGs[i],1]
    DEGs_peaks_FC<-as.data.frame(DApeak[match(DEGs_peaks1,DApeak$peak_id),3])
    gene_name<-matrix(rep(int_all_DEGs[i],length(DEGs_peaks1)),ncol=1)
    DEGs_FC<-matrix(rep(up_gene[match(int_all_DEGs[i],DEGs),3],ncol=1))
    different_FC1<-cbind(gene_name,DEGs_FC,DEGs_peaks1,DEGs_peaks_FC)
    different_FC<-rbind(different_FC,different_FC1)
  }
  sum(different_FC[,2]<0) #up DEGs
  sum(different_FC[,4]<0) #up DARs
  colnames(different_FC)<-c("S8495_resistant_DEGs","DEGs_FC","S8495_resistant_DARs","DARs_FC")
  S8495_resistant_DARs<-as.character(different_FC[which(different_FC[,2]<0 & different_FC[,4]<0),3])
  
  setwd(paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier_median_leave_one/classifier',k,sep='_'))
  homer_anno_up<-fread("sensitivity-resistance-S8495-AUC-P0.01-logFC0-DESeq2-up-annotation.txt",header=T,data.table=F)
  colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
  Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
  length(unique(Anno_gene_100Kb_up$`Gene Name`))
  
  setwd(paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier_median_RNAseq_leave_one/classifier',k,sep='_'))
  df1<-read.table("sensitivity-resistance-all-DESeq2_S8495_AUC_RNA.txt",header=T,sep="\t")
  up_gene<-df1[which(df1$pvalue<0.05&df1$log2FoldChange>0),]
  DEGs<-as.character(up_gene[,1])
  int_all_DEGs<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)
  all_peaks_anno<-homer_anno_up
  
  different_FC<-NULL
  for(i in 1:length(int_all_DEGs)){
    DEGs_peaks1<-all_peaks_anno[all_peaks_anno$`Gene Name` %in% int_all_DEGs[i],1]
    DEGs_peaks_FC<-as.data.frame(DApeak[match(DEGs_peaks1,DApeak$peak_id),3])
    gene_name<-matrix(rep(int_all_DEGs[i],length(DEGs_peaks1)),ncol=1)
    DEGs_FC<-matrix(rep(up_gene[match(int_all_DEGs[i],DEGs),3],ncol=1))
    different_FC1<-cbind(gene_name,DEGs_FC,DEGs_peaks1,DEGs_peaks_FC)
    different_FC<-rbind(different_FC,different_FC1)
  }
  sum(different_FC[,2]>0) #up DEGs
  sum(different_FC[,4]>0) #up DARs
  colnames(different_FC)<-c("S8495_sensitive_DEGs","DEGs_FC","S8495_sensitive_DARs","DARs_FC")
  S8495_sensitive_DARs<-as.character(different_FC[which(different_FC[,2]>0 & different_FC[,4]>0),3])
  
  intersect(S8495_sensitive_DARs,S8495_resistant_DARs)
  a<-intersect(GEM_sensitive_DARs,S8495_sensitive_DARs)
  if(length(a)==0){
    GEM_sensitive_only1<-GEM_sensitive_DARs
    S8495_sensitive_only<-S8495_sensitive_DARs
  }
  if(length(a)!=0){
    GEM_sensitive_only1<-GEM_sensitive_DARs[-match(a,GEM_sensitive_DARs)]
    S8495_sensitive_only<-S8495_sensitive_DARs[-match(a,S8495_sensitive_DARs)]
  }

  b<-intersect(GEM_sensitive_only1,S8495_resistant_DARs)
  if(length(b)==0){
    GEM_sensitive_only<-GEM_sensitive_only1
    S8495_resistant_only<-S8495_resistant_DARs
  }
  if(length(b)!=0){
    GEM_sensitive_only<-GEM_sensitive_only1[-match(b,GEM_sensitive_only1)]
    S8495_resistant_only<-S8495_resistant_DARs[-match(b,S8495_resistant_DARs)]
  }
  DARs<-c(GEM_sensitive_only,S8495_sensitive_only,S8495_resistant_only)
  
  GEM_sensitive<-c("DAC-24","DAC-27","DAC-28","DAC-4","DAC-33","DAC-34","DAC-36","DAC-6","DAC-1","DAC-9","DAC-11","DAC-15","DAC-3","DAC-16","DAC-18","DAC-19","DAC-37","DAC-38")
  GEM_res_S8495_sen<-c("DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-12","DAC-14")
  GEM_res_S8495_res<-c("DAC-20","DAC-21","DAC-22","DAC-23","DAC-25","DAC-29","DAC-30","DAC-32","DAC-35","DAC-5","DAC-13","DAC-39")
  all<-c(GEM_sensitive,GEM_res_S8495_sen,GEM_res_S8495_res)
  peak_RPKM37<-peak_RPKM[,match(all,colnames(peak_RPKM))]
  peak_RPKM37_DARs<-peak_RPKM37[match(DARs,peak_RPKM[,1]),]
  peak_RPKM37_DARs1<-t(peak_RPKM37_DARs)
  colnames(peak_RPKM37_DARs1)<-DARs
  y<-matrix(c(rep("GEM_sensitive",18),rep("GEM_res_S8495_sen",7),rep("GEM_res_S8495_res",12)),ncol=1)
  colnames(y)<-"y"
  data<-as.data.frame(cbind(peak_RPKM37_DARs1,y))
  library(e1071)
  set.seed(12345)
  setwd(paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier_median_leave_one/classifier',k,sep='_'))
  text<-read.table("text.txt",sep = '\t',header=F)##数据输出
  
  text_data<-data[match(as.character(text[1,1]),all),]
  train_data<-data[-(match(as.character(text[1,1]),all)),]
  #SvmFit<-svm(y~.,data=train_data,type="C-classification",kernel="radial",cost=5,gamma=1,scale=FALSE)
  #head(SvmFit$decision.values)
  #yPred<-predict(SvmFit,text_data)
  #ConfM<-table(yPred,text_data$y)
  #Err<-(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  
  #### randomForest
  #library(randomForest)
  #rf_train<-randomForest(as.factor(y)~.,data=train_data,mtry=6,ntree=500,importance=TRUE,proximity=TRUE)
  #importance<-importance(rf_train)
  #yPred<- predict(rf_train,newdata=text_data)
  #ConfM<-table(yPred,text_data$y)
  #Err<-(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  
  ### NaiveBayes
  library(klaR)
  fit_Bayes1=NaiveBayes(y~.,data=train_data)
  yPred=predict(fit_Bayes1,text_data)
  ConfM<-table(text_data$y,yPred$class)
  Err=sum(as.numeric(as.numeric(yPred$class)!=as.numeric(text_data$y)))/nrow(text_data) 
  
  #library(caret)
  #knn.model1 = knn3(y~.,data = train_data, k = 5) 
  #yPred = predict(knn.model1,text_data,type = "class") 
  #ConfM<-table(yPred,text_data$y)
  #Err=(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  cat(Err,"\n")
  if(sum(GEM_sensitive==as.character(text[1,1]))==1){
    text_label<-"GEM_sensitive"
  }
  if(sum(GEM_res_S8495_sen==as.character(text[1,1]))==1){
    text_label<-"GEM_res_S8495_sen"
  }
  if(sum(GEM_res_S8495_res==as.character(text[1,1]))==1){
    text_label<-"GEM_res_S8495_res"
  }
  svm_result10<-as.data.frame(matrix(c(i,length(DARs),text_label,as.character(yPred),Err,1-Err),nrow=1))
  svm_result0<-rbind(svm_result0,svm_result10)
  
}
colnames(svm_result0)<-c("Times","length of DARs","text_label","predict_label","error","precision")
sum(as.numeric(as.character(svm_result0[,6])))





########################################################################################
###############独立验证 10 fold cross validation #########################################################################
#SVM 独立训练 独立验证 10 fold cross validation
rm(list=ls())
library(DESeq2)
library("data.table")
library(dplyr)

setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_count_name<-read.table("peak_count.txt",sep="\t",header=T) #没有全0的
colnames(peak_count_name)<-gsub("\\.","-",colnames(peak_count_name))

setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
library("data.table")
exp1<-fread("star_rsem.GeneSymbol.readscount.xls",header=T,data.table=F)
exp2<-floor(exp1[,-1])
delete<-apply(exp2,1,function(x) mean(x==0))
cou<-which(delete>0.75)####在超过75%样本以上的都是0的基因删去
exp<-exp2[(-cou),]
geneid<-exp1[(-cou),1]

##median
GEM_sensitive1<-c("DAC-24","DAC-27","DAC-28","DAC-4","DAC-33","DAC-34","DAC-36","DAC-6","DAC-1","DAC-9","DAC-11","DAC-15","DAC-3","DAC-16","DAC-18","DAC-19","DAC-37","DAC-38")
GEM_res_S8495_sen1<-c("DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-12","DAC-14")
GEM_res_S8495_res1<-c("DAC-20","DAC-21","DAC-22","DAC-23","DAC-25","DAC-29","DAC-30","DAC-32","DAC-35","DAC-5","DAC-13","DAC-39")
##
library(doParallel)
library(parallel)
library(iterators)

cl <- makeCluster(20)
registerDoParallel(cl)

for(i in 76:100){
  cat(i,"\n")
  #GEM_sensitive1<-c("DAC-24","DAC-27","DAC-28","DAC-4","DAC-33","DAC-34","DAC-36","DAC-6","DAC-1","DAC-9","DAC-11","DAC-15","DAC-3","DAC-16","DAC-18","DAC-19","DAC-37","DAC-38")
  #GEM_res_S8495_sen1<-c("DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-12","DAC-14")
  #GEM_res_S8495_res1<-c("DAC-20","DAC-21","DAC-22","DAC-23","DAC-25","DAC-29","DAC-30","DAC-32","DAC-35","DAC-5","DAC-13","DAC-39")
  
  
  GEM_sensitive1<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
                   "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
                   "DAC-19","DAC-37","DAC-39","DAC-38")
  GEM_res_S8495_sen1<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
  GEM_res_S8495_res1<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")
  all1<-c(GEM_sensitive1,GEM_res_S8495_sen1,GEM_res_S8495_res1)
  text<-all1[sample(1:length(all1),3)]
  
  all<-all1[-match(text,all1)]
  GEM_sensitive<-intersect(GEM_sensitive1,all)
  GEM_res_S8495_sen<-intersect(GEM_res_S8495_sen1,all)
  GEM_res_S8495_res<-intersect(GEM_res_S8495_res1,all)
  
  sensitivity<-GEM_sensitive
  resistance<-c(GEM_res_S8495_sen,GEM_res_S8495_res)
  resistance_peak<-peak_count_name[,match(resistance,colnames(peak_count_name))]
  sensitivity_peak<-peak_count_name[,match(sensitivity,colnames(peak_count_name))]
  peak_name<-peak_count_name[,1]
  colDate<-data.frame(row.names = c(as.vector(resistance),as.vector(sensitivity)),
                      condition=factor(c(rep("resistance",length(resistance)),rep("sensitivity",length(sensitivity))))
  )
  datexpr<-cbind(resistance_peak,sensitivity_peak)##前面的相对于后面的上下调
  counts <- apply(datexpr,2,as.numeric)   ###矩阵中必须是数值
  dds<-DESeqDataSetFromMatrix(countData = counts,colData = colDate,design = ~condition)
  dds<-DESeq(dds)##进行标准化分析
  res<-results(dds)##将结果输出
  res<-as.data.frame(res)
  res<-cbind(peak_count_name[,1],res)  ##对数据增加一列
  colnames(res)<- c('peak_id',"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj" )##对列名进行修改
  
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier/CV10',i,sep='_')
  if (!file.exists(outputPath)){
    dir.create(outputPath, recursive=TRUE)
  }
  write.table(res,paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_GEM_AUC_ATAC.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  write.table(text,paste(outputPath,'/',"text.txt",sep=""),sep = '\t',col.names = F,row.names = F,quote = FALSE)##数据输出
  resSig<-res[which(res$pvalue<0.01),]
  # pvalue padj
  ##新增一列，将log2FoldChange>0标注为up，<0标准为down
  resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
  resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
  up_gene<-as.data.frame(resSig[which(resSig$log2FoldChange>0),1])
  down_gene<-as.data.frame(resSig[which(resSig$log2FoldChange<0),1])
  a_up<-apply(up_gene,1,function(x) unlist(strsplit(as.character(x), "_")))
  d_up<-t(a_up)
  e_up<-data.frame(rep("+",nrow(d_up)))
  peaks_up<-cbind(up_gene,d_up,e_up)
  colnames(peaks_up)<-c("peak_id","chr","start","end","strand")
  setwd(paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier/CV10',i,sep='_'))
  out_file<-"sensitivity-resistance-GEM-AUC-P0.01-logFC0"
  write.table(peaks_up,file=paste(out_file,"-DESeq2-up.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  a_down<-apply(down_gene,1,function(x) unlist(strsplit(as.character(x), "_")))
  d_down<-t(a_down)
  e_down<-data.frame(rep("+",nrow(d_down)))
  peaks_down<-cbind(down_gene,d_down,e_down)
  colnames(peaks_down)<-c("peak_id","chr","start","end","strand")
  write.table(peaks_down,file=paste(out_file,"-DESeq2-down.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
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
  
  resistance_exp<-exp[,match(resistance,colnames(exp))]
  sensitivity_exp<-exp[,match(sensitivity,colnames(exp))]
  colDate<-data.frame(row.names = c(as.vector(resistance),as.vector(sensitivity)),
                      condition=factor(c(rep("resistance",length(resistance)),rep("sensitivity",length(sensitivity))))
  )
  datexpr<-cbind(resistance_exp,sensitivity_exp)
  counts <- apply(datexpr,2,as.numeric)###矩阵中必须是数值
  dds<-DESeqDataSetFromMatrix(countData = counts,colData = colDate,design = ~condition)
  dds<-DESeq(dds)##进行标准化分析
  res<-results(dds)##将结果输出
  res<-as.data.frame(res)
  res<-cbind(geneid,res)  ##对数据增加一列
  colnames(res)<- c('gene_id',"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj" )##对列名进行修改
  write.table(res,paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_GEM_AUC_RNAseq_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  resSig<-res[which(res$pvalue<0.05),]
  resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
  resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
  up_gene<-resSig[which(resSig$log2FoldChange>0),]
  down_gene<-resSig[which(resSig$log2FoldChange<0),]
  write.table(up_gene,paste(outputPath,'/',"sensitivity-resistance-GEM-DESeq2_AUC_RNAseq_P0.05_FC0_up_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  write.table(down_gene,paste(outputPath,'/',"sensitivity-resistance-GEM-DESeq2_AUC_RNAseq_P0.05_FC0_down_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  
  
  sensitivity1<-GEM_res_S8495_sen
  resistance1<-GEM_res_S8495_res
  resistance_peak1<-peak_count_name[,match(resistance1,colnames(peak_count_name))]
  sensitivity_peak1<-peak_count_name[,match(sensitivity1,colnames(peak_count_name))]
  peak_name<-peak_count_name[,1]
  colDate1<-data.frame(row.names = c(as.vector(resistance1),as.vector(sensitivity1)),
                      condition=factor(c(rep("resistance1",length(resistance1)),rep("sensitivity1",length(sensitivity1))))
  )
  datexpr1<-cbind(resistance_peak1,sensitivity_peak1)##前面的相对于后面的上下调
  counts1 <- apply(datexpr1,2,as.numeric)   ###矩阵中必须是数值
  dds1<-DESeqDataSetFromMatrix(countData = counts1,colData = colDate1,design = ~condition)
  dds1<-DESeq(dds1)##进行标准化分析
  res1<-results(dds1)##将结果输出
  res1<-as.data.frame(res1)
  res1<-cbind(peak_count_name[,1],res1)  ##对数据增加一列
  colnames(res1)<- c('peak_id',"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj" )##对列名进行修改
  write.table(res1,paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_S8495_AUC_ATAC.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  resSig1<-res1[which(res1$pvalue<0.01),]
  # pvalue padj
  ##新增一列，将log2FoldChange>0标注为up，<0标准为down
  resSig1[which(resSig1$log2FoldChange>0),'up_down']<-'up'
  resSig1[which(resSig1$log2FoldChange<0),'up_down']<-'down'
  up_gene1<-as.data.frame(resSig1[which(resSig1$log2FoldChange>0),1])
  down_gene1<-as.data.frame(resSig1[which(resSig1$log2FoldChange<0),1])
  library(dplyr)
  a_up<-apply(up_gene1,1,function(x) unlist(strsplit(as.character(x), "_")))
  d_up<-t(a_up)
  e_up<-data.frame(rep("+",nrow(d_up)))
  peaks_up1<-cbind(up_gene1,d_up,e_up)
  colnames(peaks_up1)<-c("peak_id","chr","start","end","strand")
  setwd(paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier/CV10',i,sep='_'))
  out_file1<-"sensitivity-resistance-S8495-AUC-P0.01-logFC0"
  write.table(peaks_up1,file=paste(out_file1,"-DESeq2-up.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  a_down<-apply(down_gene1,1,function(x) unlist(strsplit(as.character(x), "_")))
  d_down<-t(a_down)
  e_down<-data.frame(rep("+",nrow(d_down)))
  peaks_down1<-cbind(down_gene1,d_down,e_down)
  colnames(peaks_down1)<-c("peak_id","chr","start","end","strand")
  write.table(peaks_down1,file=paste(out_file1,"-DESeq2-down.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  #####  bed
  f_up<-data.frame(rep(1,nrow(d_up)))
  peaks_bed_up<-cbind(d_up,up_gene1,f_up,e_up)
  colnames(peaks_bed_up)<-c("chr","start","end","peak_id","score","strand")
  write.table(peaks_bed_up,file=paste(out_file1,"-DESeq2-up.bed",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  f_down<-data.frame(rep(1,nrow(d_down)))
  peaks_bed_down<-cbind(d_down,down_gene1,f_down,e_down)
  colnames(peaks_bed_down)<-c("chr","start","end","peak_id","score","strand")
  write.table(peaks_bed_down,file=paste(out_file1,"-DESeq2-down.bed",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  
  resistance_exp1<-exp[,match(resistance1,colnames(exp))]
  sensitivity_exp1<-exp[,match(sensitivity1,colnames(exp))]
  colDate1<-data.frame(row.names = c(as.vector(resistance1),as.vector(sensitivity1)),
                       condition=factor(c(rep("resistance1",length(resistance1)),rep("sensitivity1",length(sensitivity1))))
  )
  datexpr1<-cbind(resistance_exp1,sensitivity_exp1)
  counts1 <- apply(datexpr1,2,as.numeric)###矩阵中必须是数值
  dds1<-DESeqDataSetFromMatrix(countData = counts1,colData = colDate1,design = ~condition)
  dds1<-DESeq(dds1)##进行标准化分析
  res1<-results(dds1)##将结果输出
  res1<-as.data.frame(res1)
  res1<-cbind(geneid,res1)  ##对数据增加一列
  colnames(res1)<- c('gene_id',"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj" )##对列名进行修改
  write.table(res1,paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_S8495_AUC_RNAseq_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  resSig1<-res1[which(res1$pvalue<0.05),]
  resSig1[which(resSig1$log2FoldChange>0),'up_down']<-'up'
  resSig1[which(resSig1$log2FoldChange<0),'up_down']<-'down'
  up_gene1<-resSig1[which(resSig1$log2FoldChange>0),]
  down_gene1<-resSig1[which(resSig1$log2FoldChange<0),]
  write.table(up_gene1,paste(outputPath,'/',"sensitivity-resistance-S8495-DESeq2_AUC_RNAseq_P0.05_FC0_up_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  write.table(down_gene1,paste(outputPath,'/',"sensitivity-resistance-S8495-DESeq2_AUC_RNAseq_P0.05_FC0_down_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

}


###################################################
################## RNAseq ATACseq联合分析筛选biomarker ######
rm(list=ls())
library(data.table)
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_RPKM<-read.table("peak_RPKM.txt",sep="\t",header=T) #没有全0的
colnames(peak_RPKM)<-gsub("\\.","-",colnames(peak_RPKM))

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
svm_result<-matrix(0,100,6)
for(k in 51:100){
  cat(k,"\n")
  setwd(paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier/CV10',k,sep='_'))
  homer_anno_up<-fread("sensitivity-resistance-GEM-up-annotation.txt",header=T,data.table=F)
  colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
  Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
  length(unique(Anno_gene_100Kb_up$`Gene Name`))
  
  up_gene<-read.table("sensitivity-resistance-GEM-DESeq2_AUC_RNAseq_P0.05_FC0_up_DEGs.txt",header=T,sep="\t")
  DEGs<-as.character(up_gene[,1])
  int_all_DEGs<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)
  
  all_peaks_anno<-homer_anno_up
  DApeak<-fread("sensitivity-resistance-all-DESeq2_GEM_AUC_ATAC.txt",header=T,data.table=F)
  
  
  different_FC<-NULL
  for(i in 1:length(int_all_DEGs)){
    DEGs_peaks1<-all_peaks_anno[all_peaks_anno$`Gene Name` %in% int_all_DEGs[i],1]
    DEGs_peaks_FC<-as.data.frame(DApeak[match(DEGs_peaks1,DApeak$peak_id),3])
    gene_name<-matrix(rep(int_all_DEGs[i],length(DEGs_peaks1)),ncol=1)
    DEGs_FC<-matrix(rep(up_gene[match(int_all_DEGs[i],DEGs),3],ncol=1))
    different_FC1<-cbind(gene_name,DEGs_FC,DEGs_peaks1,DEGs_peaks_FC)
    different_FC<-rbind(different_FC,different_FC1)
  }
  sum(different_FC[,2]>0) #up DEGs
  sum(different_FC[,4]>0) #up DARs
  colnames(different_FC)<-c("GEM_sensitive_DEGs","DEGs_FC","GEM_sensitive_DARs","DARs_FC")
  GEM_sensitive_DARs1<-as.character(different_FC[which(different_FC[,2]>0 & different_FC[,4]>0),3])
  
  text<-read.table("text.txt",sep = '\t',header=F)
  peak_RPKM37_GEM<-peak_RPKM37[match(GEM_sensitive_DARs1,peak_RPKM[,1]),]
  peak_RPKM37_GEM1<-peak_RPKM37_GEM[,-match(as.character(text[,1]),colnames(peak_RPKM37_GEM))]
  want_AUC_chemo1<-want_AUC_chemo[,-match(as.character(text[,1]),colnames(want_AUC_chemo))]
  GEM_result<-matrix(0,nrow(peak_RPKM37_GEM1),3)
  for(j in 1:nrow(peak_RPKM37_GEM1)){
    cor_result<-cor.test(as.numeric(peak_RPKM37_GEM1[j,]),as.numeric(want_AUC_chemo1[1,]),alternative = "two.sided",method = "pearson")
    GEM_result[j,1]<-GEM_sensitive_DARs1[j]
    GEM_result[j,2]<-cor_result$p.value
    GEM_result[j,3]<-cor_result$estimate
  }
  GEM_sensitive_DARs<-GEM_result[which(GEM_result[,2]<0.05),1]
  
  
  homer_anno_up<-fread("sensitivity-resistance-S8495-down-annotation.txt",header=T,data.table=F)
  colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
  Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
  length(unique(Anno_gene_100Kb_up$`Gene Name`))
  
  up_gene<-read.table("sensitivity-resistance-S8495-DESeq2_AUC_RNAseq_P0.05_FC0_down_DEGs.txt",header=T,sep="\t")
  DEGs<-as.character(up_gene[,1])
  int_all_DEGs<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)
  
  all_peaks_anno<-homer_anno_up
  DApeak<-fread("sensitivity-resistance-all-DESeq2_S8495_AUC_ATAC.txt",header=T,data.table=F)
  
  
  different_FC<-NULL
  for(i in 1:length(int_all_DEGs)){
    DEGs_peaks1<-all_peaks_anno[all_peaks_anno$`Gene Name` %in% int_all_DEGs[i],1]
    DEGs_peaks_FC<-as.data.frame(DApeak[match(DEGs_peaks1,DApeak$peak_id),3])
    gene_name<-matrix(rep(int_all_DEGs[i],length(DEGs_peaks1)),ncol=1)
    DEGs_FC<-matrix(rep(up_gene[match(int_all_DEGs[i],DEGs),3],ncol=1))
    different_FC1<-cbind(gene_name,DEGs_FC,DEGs_peaks1,DEGs_peaks_FC)
    different_FC<-rbind(different_FC,different_FC1)
  }
  sum(different_FC[,2]<0) #up DEGs
  sum(different_FC[,4]<0) #up DARs
  colnames(different_FC)<-c("S8495_resistant_DEGs","DEGs_FC","S8495_resistant_DARs","DARs_FC")
  S8495_resistant_DARs1<-as.character(different_FC[which(different_FC[,2]<0 & different_FC[,4]<0),3])
  
  text<-read.table("text.txt",sep = '\t',header=F)
  peak_RPKM37_S8495<-peak_RPKM37[match(S8495_resistant_DARs1,peak_RPKM[,1]),]
  peak_RPKM37_S84951<-peak_RPKM37_S8495[,-match(as.character(text[,1]),colnames(peak_RPKM37_S8495))]
  want_AUC_HDAC1<-want_AUC_HDAC[,-match(as.character(text[,1]),colnames(want_AUC_HDAC))]
  S8495_result<-matrix(0,nrow(peak_RPKM37_S84951),3)
  for(j in 1:nrow(peak_RPKM37_S84951)){
    cor_result<-cor.test(as.numeric(peak_RPKM37_S84951[j,]),as.numeric(want_AUC_HDAC1[1,]),alternative = "two.sided",method = "pearson")
    S8495_result[j,1]<-S8495_resistant_DARs1[j]
    S8495_result[j,2]<-cor_result$p.value
    S8495_result[j,3]<-cor_result$estimate
  }
  S8495_resistant_DARs<-S8495_result[which(S8495_result[,2]<0.05),1]
  
  
  homer_anno_up<-fread("sensitivity-resistance-S8495-up-annotation.txt",header=T,data.table=F)
  colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
  Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
  length(unique(Anno_gene_100Kb_up$`Gene Name`))
  
  up_gene<-read.table("sensitivity-resistance-S8495-DESeq2_AUC_RNAseq_P0.05_FC0_up_DEGs.txt",header=T,sep="\t")
  DEGs<-as.character(up_gene[,1])
  int_all_DEGs<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)
  
  all_peaks_anno<-homer_anno_up
  DApeak<-fread("sensitivity-resistance-all-DESeq2_S8495_AUC_ATAC.txt",header=T,data.table=F)
  
  different_FC<-NULL
  for(i in 1:length(int_all_DEGs)){
    DEGs_peaks1<-all_peaks_anno[all_peaks_anno$`Gene Name` %in% int_all_DEGs[i],1]
    DEGs_peaks_FC<-as.data.frame(DApeak[match(DEGs_peaks1,DApeak$peak_id),3])
    gene_name<-matrix(rep(int_all_DEGs[i],length(DEGs_peaks1)),ncol=1)
    DEGs_FC<-matrix(rep(up_gene[match(int_all_DEGs[i],DEGs),3],ncol=1))
    different_FC1<-cbind(gene_name,DEGs_FC,DEGs_peaks1,DEGs_peaks_FC)
    different_FC<-rbind(different_FC,different_FC1)
  }
  sum(different_FC[,2]>0) #up DEGs
  sum(different_FC[,4]>0) #up DARs
  colnames(different_FC)<-c("S8495_sensitive_DEGs","DEGs_FC","S8495_sensitive_DARs","DARs_FC")
  S8495_sensitive_DARs1<-as.character(different_FC[which(different_FC[,2]>0 & different_FC[,4]>0),3])
  
  text<-read.table("text.txt",sep = '\t',header=F)
  peak_RPKM37_S8495<-peak_RPKM37[match(S8495_sensitive_DARs1,peak_RPKM[,1]),]
  peak_RPKM37_S84951<-peak_RPKM37_S8495[,-match(as.character(text[,1]),colnames(peak_RPKM37_S8495))]
  want_AUC_HDAC1<-want_AUC_HDAC[,-match(as.character(text[,1]),colnames(want_AUC_HDAC))]
  S8495_result<-matrix(0,nrow(peak_RPKM37_S84951),3)
  for(j in 1:nrow(peak_RPKM37_S84951)){
    cor_result<-cor.test(as.numeric(peak_RPKM37_S84951[j,]),as.numeric(want_AUC_HDAC1[1,]),alternative = "two.sided",method = "pearson")
    S8495_result[j,1]<-S8495_sensitive_DARs1[j]
    S8495_result[j,2]<-cor_result$p.value
    S8495_result[j,3]<-cor_result$estimate
  }
  S8495_sensitive_DARs<-S8495_result[which(S8495_result[,2]<0.05),1]
  
  
  intersect(S8495_sensitive_DARs,S8495_resistant_DARs)
  a<-intersect(GEM_sensitive_DARs,S8495_sensitive_DARs)
  if(length(a)==0){
    GEM_sensitive_only1<-GEM_sensitive_DARs
    S8495_sensitive_only<-S8495_sensitive_DARs
  }
  if(length(a)!=0){
    GEM_sensitive_only1<-GEM_sensitive_DARs[-match(a,GEM_sensitive_DARs)]
    S8495_sensitive_only<-S8495_sensitive_DARs[-match(a,S8495_sensitive_DARs)]
  }
  
  b<-intersect(GEM_sensitive_only1,S8495_resistant_DARs)
  if(length(b)==0){
    GEM_sensitive_only<-GEM_sensitive_only1
    S8495_resistant_only<-S8495_resistant_DARs
  }
  if(length(b)!=0){
    GEM_sensitive_only<-GEM_sensitive_only1[-match(b,GEM_sensitive_only1)]
    S8495_resistant_only<-S8495_resistant_DARs[-match(b,S8495_resistant_DARs)]
  }
  DARs<-c(GEM_sensitive_only,S8495_sensitive_only,S8495_resistant_only)
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
  peak_RPKM37_DARs<-peak_RPKM37[match(DARs,peak_RPKM[,1]),]
  peak_RPKM37_DARs1<-t(peak_RPKM37_DARs)
  y<-matrix(c(rep("GEM_sensitive",24),rep("GEM_res_S8495_sen",7),rep("GEM_res_S8495_res",6)),ncol=1)
  colnames(y)<-"y"
  data<-as.data.frame(cbind(peak_RPKM37_DARs1,y))
  library(e1071)
  #set.seed(12345)
  text<-read.table("text.txt",sep = '\t',header=F)##数据输出
  
  text_data<-data[match(as.character(text[,1]),all),]
  train_data<-data[-(match(as.character(text[,1]),all)),]
  SvmFit<-svm(y~.,data=train_data,type="C-classification",kernel="radial",cost=5,gamma=1,scale=FALSE)
  head(SvmFit$decision.values)
  yPred<-predict(SvmFit,text_data)
  ConfM<-table(yPred,text_data$y)
  Err<-(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  
  svm_result[k,1]<-k
  svm_result[k,2]<-paste(length(GEM_sensitive_only),length(S8495_sensitive_only),length(S8495_resistant_only),length(DARs),sep=";")
  svm_result[k,3]<-paste0(as.character(text_data[,ncol(text_data)]),collapse =",")
  svm_result[k,4]<-paste0(as.character(yPred),collapse =",")
  svm_result[k,5]<-Err
  svm_result[k,6]<-1-Err
}
colnames(svm_result)<-c("Times","length of DARs","text_label","predict_label","error","precision")
mean(as.numeric(as.character(svm_result[,6])))
setwd('~/xjj/drug/drug_result/HDACi_chemo618/classifier')
write.table(svm_result,"svm_result_50-100fold_peak_gene_drug_correlationP0.05.txt",sep="\t",col.names=T,row.names=F)



#######################  没有和药物有关联
################## RNAseq ATACseq联合分析筛选biomarker ######
rm(list=ls())
library(data.table)
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_RPKM<-read.table("peak_RPKM.txt",sep="\t",header=T) #没有全0的
colnames(peak_RPKM)<-gsub("\\.","-",colnames(peak_RPKM))

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
svm_result<-matrix(0,40,6)
for(k in 1:40){
  cat(k,"\n")
  setwd(paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier4/CV10',k,sep='_'))
  homer_anno_up<-fread("sensitivity-resistance-GEM-up-annotation.txt",header=T,data.table=F)
  colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
  Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
  length(unique(Anno_gene_100Kb_up$`Gene Name`))
  
  up_gene<-read.table("sensitivity-resistance-GEM-DESeq2_AUC_RNAseq_P0.05_FC0_up_DEGs.txt",header=T,sep="\t")
  DEGs<-as.character(up_gene[,1])
  int_all_DEGs<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)
  
  all_peaks_anno<-homer_anno_up
  DApeak<-fread("sensitivity-resistance-all-DESeq2_GEM_AUC_ATAC.txt",header=T,data.table=F)
  
  
  different_FC<-NULL
  for(i in 1:length(int_all_DEGs)){
    DEGs_peaks1<-all_peaks_anno[all_peaks_anno$`Gene Name` %in% int_all_DEGs[i],1]
    DEGs_peaks_FC<-as.data.frame(DApeak[match(DEGs_peaks1,DApeak$peak_id),3])
    gene_name<-matrix(rep(int_all_DEGs[i],length(DEGs_peaks1)),ncol=1)
    DEGs_FC<-matrix(rep(up_gene[match(int_all_DEGs[i],DEGs),3],ncol=1))
    different_FC1<-cbind(gene_name,DEGs_FC,DEGs_peaks1,DEGs_peaks_FC)
    different_FC<-rbind(different_FC,different_FC1)
  }
  sum(different_FC[,2]>0) #up DEGs
  sum(different_FC[,4]>0) #up DARs
  colnames(different_FC)<-c("GEM_sensitive_DEGs","DEGs_FC","GEM_sensitive_DARs","DARs_FC")
  GEM_sensitive_DARs<-as.character(different_FC[which(different_FC[,2]>0 & different_FC[,4]>0),3])
  
  homer_anno_up<-fread("sensitivity-resistance-S8495-down-annotation.txt",header=T,data.table=F)
  colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
  Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
  length(unique(Anno_gene_100Kb_up$`Gene Name`))
  
  up_gene<-read.table("sensitivity-resistance-S8495-DESeq2_AUC_RNAseq_P0.05_FC0_down_DEGs.txt",header=T,sep="\t")
  DEGs<-as.character(up_gene[,1])
  int_all_DEGs<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)
  
  all_peaks_anno<-homer_anno_up
  DApeak<-fread("sensitivity-resistance-all-DESeq2_S8495_AUC_ATAC.txt",header=T,data.table=F)
  
  
  different_FC<-NULL
  for(i in 1:length(int_all_DEGs)){
    DEGs_peaks1<-all_peaks_anno[all_peaks_anno$`Gene Name` %in% int_all_DEGs[i],1]
    DEGs_peaks_FC<-as.data.frame(DApeak[match(DEGs_peaks1,DApeak$peak_id),3])
    gene_name<-matrix(rep(int_all_DEGs[i],length(DEGs_peaks1)),ncol=1)
    DEGs_FC<-matrix(rep(up_gene[match(int_all_DEGs[i],DEGs),3],ncol=1))
    different_FC1<-cbind(gene_name,DEGs_FC,DEGs_peaks1,DEGs_peaks_FC)
    different_FC<-rbind(different_FC,different_FC1)
  }
  sum(different_FC[,2]<0) #up DEGs
  sum(different_FC[,4]<0) #up DARs
  colnames(different_FC)<-c("S8495_resistant_DEGs","DEGs_FC","S8495_resistant_DARs","DARs_FC")
  S8495_resistant_DARs<-as.character(different_FC[which(different_FC[,2]<0 & different_FC[,4]<0),3])
  
  homer_anno_up<-fread("sensitivity-resistance-S8495-up-annotation.txt",header=T,data.table=F)
  colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
  Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
  length(unique(Anno_gene_100Kb_up$`Gene Name`))
  
  up_gene<-read.table("sensitivity-resistance-S8495-DESeq2_AUC_RNAseq_P0.05_FC0_up_DEGs.txt",header=T,sep="\t")
  DEGs<-as.character(up_gene[,1])
  int_all_DEGs<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)
  
  all_peaks_anno<-homer_anno_up
  DApeak<-fread("sensitivity-resistance-all-DESeq2_S8495_AUC_ATAC.txt",header=T,data.table=F)
  
  different_FC<-NULL
  for(i in 1:length(int_all_DEGs)){
    DEGs_peaks1<-all_peaks_anno[all_peaks_anno$`Gene Name` %in% int_all_DEGs[i],1]
    DEGs_peaks_FC<-as.data.frame(DApeak[match(DEGs_peaks1,DApeak$peak_id),3])
    gene_name<-matrix(rep(int_all_DEGs[i],length(DEGs_peaks1)),ncol=1)
    DEGs_FC<-matrix(rep(up_gene[match(int_all_DEGs[i],DEGs),3],ncol=1))
    different_FC1<-cbind(gene_name,DEGs_FC,DEGs_peaks1,DEGs_peaks_FC)
    different_FC<-rbind(different_FC,different_FC1)
  }
  sum(different_FC[,2]>0) #up DEGs
  sum(different_FC[,4]>0) #up DARs
  colnames(different_FC)<-c("S8495_sensitive_DEGs","DEGs_FC","S8495_sensitive_DARs","DARs_FC")
  S8495_sensitive_DARs<-as.character(different_FC[which(different_FC[,2]>0 & different_FC[,4]>0),3])
  
  intersect(S8495_sensitive_DARs,S8495_resistant_DARs)
  a<-intersect(GEM_sensitive_DARs,S8495_sensitive_DARs)
  if(length(a)==0){
    GEM_sensitive_only1<-GEM_sensitive_DARs
    S8495_sensitive_only<-S8495_sensitive_DARs
  }
  if(length(a)!=0){
    GEM_sensitive_only1<-GEM_sensitive_DARs[-match(a,GEM_sensitive_DARs)]
    S8495_sensitive_only<-S8495_sensitive_DARs[-match(a,S8495_sensitive_DARs)]
  }
  
  b<-intersect(GEM_sensitive_only1,S8495_resistant_DARs)
  if(length(b)==0){
    GEM_sensitive_only<-GEM_sensitive_only1
    S8495_resistant_only<-S8495_resistant_DARs
  }
  if(length(b)!=0){
    GEM_sensitive_only<-GEM_sensitive_only1[-match(b,GEM_sensitive_only1)]
    S8495_resistant_only<-S8495_resistant_DARs[-match(b,S8495_resistant_DARs)]
  }
  DARs<-c(GEM_sensitive_only,S8495_sensitive_only,S8495_resistant_only)
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
  peak_RPKM37_DARs<-peak_RPKM37[match(DARs,peak_RPKM[,1]),]
  peak_RPKM37_DARs1<-t(peak_RPKM37_DARs)
  y<-matrix(c(rep("GEM_sensitive",24),rep("GEM_res_S8495_sen",7),rep("GEM_res_S8495_res",6)),ncol=1)
  colnames(y)<-"y"
  data<-as.data.frame(cbind(peak_RPKM37_DARs1,y))
  library(e1071)
  set.seed(12345)
  text<-read.table("text.txt",sep = '\t',header=F)##数据输出
  
  text_data<-data[match(as.character(text[,1]),all),]
  train_data<-data[-(match(as.character(text[,1]),all)),]
  SvmFit<-svm(y~.,data=train_data,type="C-classification",kernel="radial",cost=5,gamma=1,scale=FALSE)
  head(SvmFit$decision.values)
  yPred<-predict(SvmFit,text_data)
  ConfM<-table(yPred,text_data$y)
  Err<-(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  
  svm_result[k,1]<-k
  svm_result[k,2]<-paste(length(GEM_sensitive_only),length(S8495_sensitive_only),length(S8495_resistant_only),length(DARs),sep=";")
  svm_result[k,3]<-paste0(as.character(text_data[,ncol(text_data)]),collapse =",")
  svm_result[k,4]<-paste0(as.character(yPred),collapse =",")
  svm_result[k,5]<-Err
  svm_result[k,6]<-1-Err
}
colnames(svm_result)<-c("Times","length of DARs","text_label","predict_label","error","precision")
mean(as.numeric(as.character(svm_result[,6])))
setwd('~/xjj/drug/drug_result/HDACi_chemo618/classifier4')
write.table(svm_result,"svm_result_40fold_peak_gene.txt",sep="\t",col.names=T,row.names=F)





###########  单纯只用peaks作为特征，没有涉及DEGs
###################################################
rm(list=ls())
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_RPKM<-read.table("peak_RPKM.txt",sep="\t",header=T) #没有全0的
colnames(peak_RPKM)<-gsub("\\.","-",colnames(peak_RPKM))

GEM_sensitive<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
                 "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
                 "DAC-19","DAC-37","DAC-39","DAC-38")
GEM_res_S8495_sen<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
GEM_res_S8495_res<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")

##median
GEM_sensitive<-c("DAC-24","DAC-27","DAC-28","DAC-4","DAC-33","DAC-34","DAC-36","DAC-6","DAC-1","DAC-9","DAC-11","DAC-15","DAC-3","DAC-16","DAC-18","DAC-19","DAC-37","DAC-38")
GEM_res_S8495_sen<-c("DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-12","DAC-14")
GEM_res_S8495_res<-c("DAC-20","DAC-21","DAC-22","DAC-23","DAC-25","DAC-29","DAC-30","DAC-32","DAC-35","DAC-5","DAC-13","DAC-39")
##

svm_result<-NULL
library(doParallel)
library(parallel)
library(iterators)

cl <- makeCluster(40)
registerDoParallel(cl)
foreach(k=41:80,.combine='rbind')%dopar%{
  cat(k,"\n")
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/10fold-CV/CV10',k,sep='_')
  GEM<-read.table(paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_GEM_AUC_ATAC.txt",sep=""),sep = '\t',header=T)##数据输出
  resSig<-GEM[which(GEM$pvalue<0.05 & abs(GEM$log2FoldChange)>2),]
  #resSig<-GEM[which(GEM$pvalue<0.01),]
  resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
  resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
  cat(sum(resSig$up_down=='up'),"\n")
  GEM_sensitive_DARs<-as.character(resSig[which(resSig$log2FoldChange>0),1])
  
  S8495<-read.table(paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_S8495_AUC_ATAC.txt",sep=""),sep = '\t',header=T)##数据输出
  resSig_S8495<-S8495[which(S8495$pvalue<0.05 & abs(S8495$log2FoldChange)>2),]
  #resSig<-GEM[which(GEM$pvalue<0.01),]
  resSig_S8495[which(resSig_S8495$log2FoldChange>0),'up_down']<-'up'
  resSig_S8495[which(resSig_S8495$log2FoldChange<0),'up_down']<-'down'
  cat(sum(resSig_S8495$up_down=='up'),"\n")
  cat(sum(resSig_S8495$up_down=='down'),"\n")
  S8495_sensitive_DARs<-as.character(resSig_S8495[which(resSig_S8495$log2FoldChange>0),1])
  S8495_resistant_DARs<-as.character(resSig_S8495[which(resSig_S8495$log2FoldChange<0),1])
  a<-intersect(GEM_sensitive_DARs,S8495_sensitive_DARs)
  if(length(a)==0){
    GEM_sensitive_only1<-GEM_sensitive_DARs
    S8495_sensitive_only<-S8495_sensitive_DARs
  }
  if(length(a)!=0){
    GEM_sensitive_only1<-GEM_sensitive_DARs[-match(a,GEM_sensitive_DARs)]
    S8495_sensitive_only<-S8495_sensitive_DARs[-match(a,S8495_sensitive_DARs)]
  }
  
  b<-intersect(GEM_sensitive_only1,S8495_resistant_DARs)
  if(length(b)==0){
    GEM_sensitive_only<-GEM_sensitive_only1
    S8495_resistant_only<-S8495_resistant_DARs
  }
  if(length(b)!=0){
    GEM_sensitive_only<-GEM_sensitive_only1[-match(b,GEM_sensitive_only1)]
    S8495_resistant_only<-S8495_resistant_DARs[-match(b,S8495_resistant_DARs)]
  }
  DARs<-c(GEM_sensitive_only,S8495_sensitive_only,S8495_resistant_only)
  all<-c(GEM_sensitive,GEM_res_S8495_sen,GEM_res_S8495_res)
  peak_RPKM37<-peak_RPKM[,match(all,colnames(peak_RPKM))]
  peak_RPKM37_DARs<-peak_RPKM37[match(DARs,peak_RPKM[,1]),]
  peak_RPKM37_DARs1<-t(peak_RPKM37_DARs)
  colnames(peak_RPKM37_DARs1)<-DARs
  y<-matrix(c(rep("GEM_sensitive",24),rep("GEM_res_S8495_sen",7),rep("GEM_res_S8495_res",6)),ncol=1)
  colnames(y)<-"y"
  data<-as.data.frame(cbind(peak_RPKM37_DARs1,y))
  #library(e1071)
  #set.seed(12345)
  text<-read.table(paste(outputPath,'/',"text.txt",sep=""),sep = '\t',header=F)##数据输出
  
  text_data<-data[match(as.character(text[,1]),all),]
  train_data<-data[-(match(as.character(text[,1]),all)),]
  #SvmFit<-svm(y~.,data=train_data,type="C-classification",kernel="radial",cost=10,scale=FALSE)
  #head(SvmFit$decision.values)
  #yPred<-predict(SvmFit,text_data)
  #ConfM<-table(yPred,text_data$y)
  #Err<-(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  
  library(randomForest)
  rf_train<-randomForest(as.factor(y)~.,data=train_data,mtry=6,ntree=500,importance=TRUE,proximity=TRUE)
  importance<-importance(rf_train)
  yPred<- predict(rf_train,newdata=text_data)
  ConfM<-table(yPred,text_data$y)
  Err<-(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  
  #library(klaR)
  #fit_Bayes1=NaiveBayes(y~.,data=train_data)
  #yPred=predict(fit_Bayes1,text_data)
  #ConfM<-table(text_data$y,yPred$class)
  #Err=sum(as.numeric(as.numeric(yPred$class)!=as.numeric(text_data$y)))/nrow(text_data) 
  
  #library(caret)
  #knn.model1 = knn3(y~.,data = train_data, k = 5) 
  #yPred = predict(knn.model1,text_data,type = "class") 
  #ConfM<-table(yPred,text_data$y)
  #Err=(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  
  
  svm_result1<-as.data.frame(matrix(c(k,length(DARs),paste0(as.character(text_data[,ncol(text_data)]),collapse =","),paste0(as.character(yPred),collapse =","),Err,1-Err),nrow=1))
  svm_result<-rbind(svm_result,svm_result1)
}
stopCluster(cl)
colnames(svm_result)<-c("times","length of biomarks","text_label","predict_label","error","precision")
sum(as.numeric(as.character(svm_result[,6])))/nrow(svm_result)
setwd('~/xjj/drug/drug_result/HDACi_chemo618/classifier')
write.table(svm_result,"svm_result_50fold_peak_P0.05_FC2.txt",sep="\t",col.names=T,row.names=F)


#####  随机森林###################################################
#####  朴素贝叶斯 Naive Bayes
#####  GBM
#### randomForest
library(randomForest)
rf_train<-randomForest(as.factor(y)~.,data=train_data,mtry=6,ntree=500,importance=TRUE,proximity=TRUE)
importance<-importance(rf_train)
yPred<- predict(rf_train,newdata=text_data)
ConfM<-table(yPred,text_data$y)
Err<-(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)

### NaiveBayes
library(klaR)
fit_Bayes1=NaiveBayes(y~.,data=train_data)
yPred=predict(fit_Bayes1,text_data)
ConfM<-table(text_data$y,yPred$class)
Err=sum(as.numeric(as.numeric(yPred$class)!=as.numeric(text_data$y)))/nrow(text_data) 

### Decision trees
library(rpart)
# 训练模型
# rpart参考文档
set.seed(42) # 固定交叉验证结果
fit_dt_reg <- rpart(
  y~.,
  data=train_data,
  method = "anova", 
  control = rpart.control(cp = 0.005)
)
# 原始回归树
fit_dt_reg
# 复杂度相关数据
printcp(fit_dt_reg)
plotcp(fit_dt_reg)

# 后剪枝
fit_dt_reg_pruned <- prune(fit_dt_reg, cp = cp1SE)
print(fit_dt_reg_pruned)
summary(fit_dt_reg_pruned)

# 变量重要性数值
fit_dt_reg_pruned$variable.importance
yPred<- predict(fit_dt_reg_pruned,newdata=text_data)
ConfM<-table(yPred,text_data$y)
Err<-(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)

# 变量重要性图示
varimpdata <-
  data.frame(importance = fit_dt_reg_pruned$variable.importance)






### Logistic Regression  还不行
Library(mlogit)
logit.response =mlogit(y~.,data = train_data)


library(stats)
logit = glm(as.factor(y)~.,data = train_data,family = "binomial")
logit = glm(as.factor(y)~.,data = train_data,family = "poisson")
summary(logit)
logit.response = predict(logit,text_data,type = "response")
logit.predict = ifelse(logit.response>0.5,"+","-")
table(logit.predict,test$V16)

#library(nnet)
#mult.cere<-multinom(y~.,data = train_data)

#EW <- glm(as.factor(y)~.,data = train_data,family = "binomial")
#EW <- glm(everwrk~r_maritl+age_p,data=NH11,family=binomial)
#predEW <- with(train_data,
#               expand.grid(y=levels(y)))
#predict(EW,newdata=predEW)

### GBM  还不行
#library(caret)
#ctrl = trainControl(method = "repeatedcv", number = 5, repeats = 5)
#set.seed(300)
#m_gbm = train(V16 ~ ., data=train, method = "gbm",  metric = "Kappa", trControl = ctrl)
#gbm.predict = predict(m_gbm,test)
#table(gbm.predict,test$V16)
#accurancy2 = mean(gbm.predict == test$V16)


### knn
library(caret)
knn.model1 = knn3(y~.,data = train_data, k = 3) 
knn.response1 = predict(knn.model1,text_data,type = "class") 
ConfM<-table(knn.response1,text_data$y)
Err=(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)

### xgboost
#devtools::install_github('dmlc/xgboost',subdir='R-package')
library(xgboost)
require(xgboost)
require(methods)
require(plyr)
library(Matrix)
library(caret)
# 训练集的数据预处理
# 将trainset的1-4列（自变量）转换为矩阵
traindata1 <- as.matrix(train_data[,c(1:(ncol(train_data)-1))])
# 利用Matrix函数，将sparse参数设置为TRUE，转化为稀疏矩阵
traindata2 <- Matrix(traindata1,sparse = T)
# 将因变量转换为numeric类型，-1是为了从0开始计数
train_y <- as.numeric(train_data[,ncol(train_data)])-1
# 将自变量和因变量拼接为list
traindata <- list(data=traindata2,label=train_y)
# 构造模型需要的xgb.DMatrix对象，处理对象为稀疏矩阵
dtrain <- xgb.DMatrix(data = traindata$data, label = traindata$label) 

# 测试集的数据预处理
# 将自变量转化为矩阵
testset1 <- as.matrix(text_data[,c(1:(ncol(text_data)-1))])
# 利用Matrix函数，将sparse参数设置为TRUE，转化为稀疏矩阵
testset2 <- Matrix(testset1,sparse=T) 
# 将因变量转化为numeric
test_y <- as.numeric(text_data[,ncol(text_data)])-1
# 将自变量和因变量拼接为list
testset <- list(data=testset2,label=test_y)
# 构造模型需要的xgb.DMatrix对象，处理对象为稀疏矩阵
dtest <- xgb.DMatrix(data=testset$data,label=testset$label)
# 建立模型
model_xgb <- xgboost(data=dtrain,booster='gbtree',max_depth=6,eta=0.5,objective='multi:softmax',num_class=3,nround=10)
#用测试集预测
pre <- predict(model_xgb,newdata=dtest)
table(pre,test_y)
accurancy4 = mean(pre ==test_y)
Err=1-accurancy4
#模型评估
#xgb.cf <-caret::confusionMatrix(as.factor(pre),as.factor(test_y))
#Err=1-xgb.cf$overall[1]


################## RNAseq ATACseq联合分析筛选biomarker ######
rm(list=ls())
library(data.table)
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_RPKM<-read.table("peak_RPKM.txt",sep="\t",header=T) #没有全0的
colnames(peak_RPKM)<-gsub("\\.","-",colnames(peak_RPKM))

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
rf_result<-matrix(0,100,6)
for(k in 51:100){
  cat(k,"\n")
  setwd(paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier/CV10',k,sep='_'))
  homer_anno_up<-fread("sensitivity-resistance-GEM-up-annotation.txt",header=T,data.table=F)
  colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
  Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
  length(unique(Anno_gene_100Kb_up$`Gene Name`))
  
  up_gene<-read.table("sensitivity-resistance-GEM-DESeq2_AUC_RNAseq_P0.05_FC0_up_DEGs.txt",header=T,sep="\t")
  DEGs<-as.character(up_gene[,1])
  int_all_DEGs<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)
  
  all_peaks_anno<-homer_anno_up
  DApeak<-fread("sensitivity-resistance-all-DESeq2_GEM_AUC_ATAC.txt",header=T,data.table=F)
  
  
  different_FC<-NULL
  for(i in 1:length(int_all_DEGs)){
    DEGs_peaks1<-all_peaks_anno[all_peaks_anno$`Gene Name` %in% int_all_DEGs[i],1]
    DEGs_peaks_FC<-as.data.frame(DApeak[match(DEGs_peaks1,DApeak$peak_id),3])
    gene_name<-matrix(rep(int_all_DEGs[i],length(DEGs_peaks1)),ncol=1)
    DEGs_FC<-matrix(rep(up_gene[match(int_all_DEGs[i],DEGs),3],ncol=1))
    different_FC1<-cbind(gene_name,DEGs_FC,DEGs_peaks1,DEGs_peaks_FC)
    different_FC<-rbind(different_FC,different_FC1)
  }
  sum(different_FC[,2]>0) #up DEGs
  sum(different_FC[,4]>0) #up DARs
  colnames(different_FC)<-c("GEM_sensitive_DEGs","DEGs_FC","GEM_sensitive_DARs","DARs_FC")
  GEM_sensitive_DARs1<-as.character(different_FC[which(different_FC[,2]>0 & different_FC[,4]>0),3])
  
  text<-read.table("text.txt",sep = '\t',header=F)
  peak_RPKM37_GEM<-peak_RPKM37[match(GEM_sensitive_DARs1,peak_RPKM[,1]),]
  peak_RPKM37_GEM1<-peak_RPKM37_GEM[,-match(as.character(text[,1]),colnames(peak_RPKM37_GEM))]
  want_AUC_chemo1<-want_AUC_chemo[,-match(as.character(text[,1]),colnames(want_AUC_chemo))]
  GEM_result<-matrix(0,nrow(peak_RPKM37_GEM1),3)
  for(j in 1:nrow(peak_RPKM37_GEM1)){
    cor_result<-cor.test(as.numeric(peak_RPKM37_GEM1[j,]),as.numeric(want_AUC_chemo1[1,]),alternative = "two.sided",method = "pearson")
    GEM_result[j,1]<-GEM_sensitive_DARs1[j]
    GEM_result[j,2]<-cor_result$p.value
    GEM_result[j,3]<-cor_result$estimate
  }
  GEM_sensitive_DARs<-GEM_result[which(GEM_result[,2]<0.05),1]
  
  
  homer_anno_up<-fread("sensitivity-resistance-S8495-down-annotation.txt",header=T,data.table=F)
  colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
  Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
  length(unique(Anno_gene_100Kb_up$`Gene Name`))
  
  up_gene<-read.table("sensitivity-resistance-S8495-DESeq2_AUC_RNAseq_P0.05_FC0_down_DEGs.txt",header=T,sep="\t")
  DEGs<-as.character(up_gene[,1])
  int_all_DEGs<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)
  
  all_peaks_anno<-homer_anno_up
  DApeak<-fread("sensitivity-resistance-all-DESeq2_S8495_AUC_ATAC.txt",header=T,data.table=F)
  
  
  different_FC<-NULL
  for(i in 1:length(int_all_DEGs)){
    DEGs_peaks1<-all_peaks_anno[all_peaks_anno$`Gene Name` %in% int_all_DEGs[i],1]
    DEGs_peaks_FC<-as.data.frame(DApeak[match(DEGs_peaks1,DApeak$peak_id),3])
    gene_name<-matrix(rep(int_all_DEGs[i],length(DEGs_peaks1)),ncol=1)
    DEGs_FC<-matrix(rep(up_gene[match(int_all_DEGs[i],DEGs),3],ncol=1))
    different_FC1<-cbind(gene_name,DEGs_FC,DEGs_peaks1,DEGs_peaks_FC)
    different_FC<-rbind(different_FC,different_FC1)
  }
  sum(different_FC[,2]<0) #up DEGs
  sum(different_FC[,4]<0) #up DARs
  colnames(different_FC)<-c("S8495_resistant_DEGs","DEGs_FC","S8495_resistant_DARs","DARs_FC")
  S8495_resistant_DARs1<-as.character(different_FC[which(different_FC[,2]<0 & different_FC[,4]<0),3])
  
  text<-read.table("text.txt",sep = '\t',header=F)
  peak_RPKM37_S8495<-peak_RPKM37[match(S8495_resistant_DARs1,peak_RPKM[,1]),]
  peak_RPKM37_S84951<-peak_RPKM37_S8495[,-match(as.character(text[,1]),colnames(peak_RPKM37_S8495))]
  want_AUC_HDAC1<-want_AUC_HDAC[,-match(as.character(text[,1]),colnames(want_AUC_HDAC))]
  S8495_result<-matrix(0,nrow(peak_RPKM37_S84951),3)
  for(j in 1:nrow(peak_RPKM37_S84951)){
    cor_result<-cor.test(as.numeric(peak_RPKM37_S84951[j,]),as.numeric(want_AUC_HDAC1[1,]),alternative = "two.sided",method = "pearson")
    S8495_result[j,1]<-S8495_resistant_DARs1[j]
    S8495_result[j,2]<-cor_result$p.value
    S8495_result[j,3]<-cor_result$estimate
  }
  S8495_resistant_DARs<-S8495_result[which(S8495_result[,2]<0.05),1]
  
  
  homer_anno_up<-fread("sensitivity-resistance-S8495-up-annotation.txt",header=T,data.table=F)
  colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
  Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
  length(unique(Anno_gene_100Kb_up$`Gene Name`))
  
  up_gene<-read.table("sensitivity-resistance-S8495-DESeq2_AUC_RNAseq_P0.05_FC0_up_DEGs.txt",header=T,sep="\t")
  DEGs<-as.character(up_gene[,1])
  int_all_DEGs<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)
  
  all_peaks_anno<-homer_anno_up
  DApeak<-fread("sensitivity-resistance-all-DESeq2_S8495_AUC_ATAC.txt",header=T,data.table=F)
  
  different_FC<-NULL
  for(i in 1:length(int_all_DEGs)){
    DEGs_peaks1<-all_peaks_anno[all_peaks_anno$`Gene Name` %in% int_all_DEGs[i],1]
    DEGs_peaks_FC<-as.data.frame(DApeak[match(DEGs_peaks1,DApeak$peak_id),3])
    gene_name<-matrix(rep(int_all_DEGs[i],length(DEGs_peaks1)),ncol=1)
    DEGs_FC<-matrix(rep(up_gene[match(int_all_DEGs[i],DEGs),3],ncol=1))
    different_FC1<-cbind(gene_name,DEGs_FC,DEGs_peaks1,DEGs_peaks_FC)
    different_FC<-rbind(different_FC,different_FC1)
  }
  sum(different_FC[,2]>0) #up DEGs
  sum(different_FC[,4]>0) #up DARs
  colnames(different_FC)<-c("S8495_sensitive_DEGs","DEGs_FC","S8495_sensitive_DARs","DARs_FC")
  S8495_sensitive_DARs1<-as.character(different_FC[which(different_FC[,2]>0 & different_FC[,4]>0),3])
  
  text<-read.table("text.txt",sep = '\t',header=F)
  peak_RPKM37_S8495<-peak_RPKM37[match(S8495_sensitive_DARs1,peak_RPKM[,1]),]
  peak_RPKM37_S84951<-peak_RPKM37_S8495[,-match(as.character(text[,1]),colnames(peak_RPKM37_S8495))]
  want_AUC_HDAC1<-want_AUC_HDAC[,-match(as.character(text[,1]),colnames(want_AUC_HDAC))]
  S8495_result<-matrix(0,nrow(peak_RPKM37_S84951),3)
  for(j in 1:nrow(peak_RPKM37_S84951)){
    cor_result<-cor.test(as.numeric(peak_RPKM37_S84951[j,]),as.numeric(want_AUC_HDAC1[1,]),alternative = "two.sided",method = "pearson")
    S8495_result[j,1]<-S8495_sensitive_DARs1[j]
    S8495_result[j,2]<-cor_result$p.value
    S8495_result[j,3]<-cor_result$estimate
  }
  S8495_sensitive_DARs<-S8495_result[which(S8495_result[,2]<0.05),1]
  
  
  intersect(S8495_sensitive_DARs,S8495_resistant_DARs)
  a<-intersect(GEM_sensitive_DARs,S8495_sensitive_DARs)
  if(length(a)==0){
    GEM_sensitive_only1<-GEM_sensitive_DARs
    S8495_sensitive_only<-S8495_sensitive_DARs
  }
  if(length(a)!=0){
    GEM_sensitive_only1<-GEM_sensitive_DARs[-match(a,GEM_sensitive_DARs)]
    S8495_sensitive_only<-S8495_sensitive_DARs[-match(a,S8495_sensitive_DARs)]
  }
  
  b<-intersect(GEM_sensitive_only1,S8495_resistant_DARs)
  if(length(b)==0){
    GEM_sensitive_only<-GEM_sensitive_only1
    S8495_resistant_only<-S8495_resistant_DARs
  }
  if(length(b)!=0){
    GEM_sensitive_only<-GEM_sensitive_only1[-match(b,GEM_sensitive_only1)]
    S8495_resistant_only<-S8495_resistant_DARs[-match(b,S8495_resistant_DARs)]
  }
  DARs<-c(GEM_sensitive_only,S8495_sensitive_only,S8495_resistant_only)
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
  peak_RPKM37_DARs<-peak_RPKM37[match(DARs,peak_RPKM[,1]),]
  peak_RPKM37_DARs1<-t(peak_RPKM37_DARs)
  colnames(peak_RPKM37_DARs1)<-DARs
  y<-matrix(c(rep("GEM_sensitive",24),rep("GEM_res_S8495_sen",7),rep("GEM_res_S8495_res",6)),ncol=1)
  colnames(y)<-"y"
  data<-as.data.frame(cbind(peak_RPKM37_DARs1,y))
  text<-read.table("text.txt",sep = '\t',header=F)##数据输出
  text_data<-data[match(as.character(text[,1]),all),]
  train_data<-data[-(match(as.character(text[,1]),all)),]
  library(randomForest)
  rf_train<-randomForest(as.factor(y)~.,data=train_data,mtry=6,ntree=500,importance=TRUE,proximity=TRUE)
  importance<-importance(rf_train)
  yPred<- predict(rf_train,newdata=text_data)
  ConfM<-table(yPred,text_data$y)
  Err<-(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  
  rf_result[k,1]<-k
  rf_result[k,2]<-paste(length(GEM_sensitive_only),length(S8495_sensitive_only),length(S8495_resistant_only),length(DARs),sep=";")
  rf_result[k,3]<-paste0(as.character(text_data[,ncol(text_data)]),collapse =",")
  rf_result[k,4]<-paste0(as.character(yPred),collapse =",")
  rf_result[k,5]<-Err
  rf_result[k,6]<-1-Err
}
colnames(rf_result)<-c("Times","length of DARs","text_label","predict_label","error","precision")
mean(as.numeric(as.character(rf_result[,6])))
setwd('~/xjj/drug/drug_result/HDACi_chemo618/classifier')
write.table(rf_result,"rf_result_51_100fold_peak_gene_drug_correlation_P0.05.txt",sep="\t",col.names=T,row.names=F)


#######################  没有和药物有关联
################## RNAseq ATACseq联合分析筛选biomarker ######
rm(list=ls())
library(data.table)
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_RPKM<-read.table("peak_RPKM.txt",sep="\t",header=T) #没有全0的
colnames(peak_RPKM)<-gsub("\\.","-",colnames(peak_RPKM))

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
rf_result<-matrix(0,50,6)
for(k in 1:50){
  cat(k,"\n")
  setwd(paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier/CV10',k,sep='_'))
  homer_anno_up<-fread("sensitivity-resistance-GEM-up-annotation.txt",header=T,data.table=F)
  colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
  Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
  length(unique(Anno_gene_100Kb_up$`Gene Name`))
  
  up_gene<-read.table("sensitivity-resistance-GEM-DESeq2_AUC_RNAseq_P0.05_FC0_up_DEGs.txt",header=T,sep="\t")
  DEGs<-as.character(up_gene[,1])
  int_all_DEGs<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)
  
  all_peaks_anno<-homer_anno_up
  DApeak<-fread("sensitivity-resistance-all-DESeq2_GEM_AUC_ATAC.txt",header=T,data.table=F)
  
  
  different_FC<-NULL
  for(i in 1:length(int_all_DEGs)){
    DEGs_peaks1<-all_peaks_anno[all_peaks_anno$`Gene Name` %in% int_all_DEGs[i],1]
    DEGs_peaks_FC<-as.data.frame(DApeak[match(DEGs_peaks1,DApeak$peak_id),3])
    gene_name<-matrix(rep(int_all_DEGs[i],length(DEGs_peaks1)),ncol=1)
    DEGs_FC<-matrix(rep(up_gene[match(int_all_DEGs[i],DEGs),3],ncol=1))
    different_FC1<-cbind(gene_name,DEGs_FC,DEGs_peaks1,DEGs_peaks_FC)
    different_FC<-rbind(different_FC,different_FC1)
  }
  sum(different_FC[,2]>0) #up DEGs
  sum(different_FC[,4]>0) #up DARs
  colnames(different_FC)<-c("GEM_sensitive_DEGs","DEGs_FC","GEM_sensitive_DARs","DARs_FC")
  GEM_sensitive_DARs<-as.character(different_FC[which(different_FC[,2]>0 & different_FC[,4]>0),3])
  
  homer_anno_up<-fread("sensitivity-resistance-S8495-down-annotation.txt",header=T,data.table=F)
  colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
  Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
  length(unique(Anno_gene_100Kb_up$`Gene Name`))
  
  up_gene<-read.table("sensitivity-resistance-S8495-DESeq2_AUC_RNAseq_P0.05_FC0_down_DEGs.txt",header=T,sep="\t")
  DEGs<-as.character(up_gene[,1])
  int_all_DEGs<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)
  
  all_peaks_anno<-homer_anno_up
  DApeak<-fread("sensitivity-resistance-all-DESeq2_S8495_AUC_ATAC.txt",header=T,data.table=F)
  
  
  different_FC<-NULL
  for(i in 1:length(int_all_DEGs)){
    DEGs_peaks1<-all_peaks_anno[all_peaks_anno$`Gene Name` %in% int_all_DEGs[i],1]
    DEGs_peaks_FC<-as.data.frame(DApeak[match(DEGs_peaks1,DApeak$peak_id),3])
    gene_name<-matrix(rep(int_all_DEGs[i],length(DEGs_peaks1)),ncol=1)
    DEGs_FC<-matrix(rep(up_gene[match(int_all_DEGs[i],DEGs),3],ncol=1))
    different_FC1<-cbind(gene_name,DEGs_FC,DEGs_peaks1,DEGs_peaks_FC)
    different_FC<-rbind(different_FC,different_FC1)
  }
  sum(different_FC[,2]<0) #up DEGs
  sum(different_FC[,4]<0) #up DARs
  colnames(different_FC)<-c("S8495_resistant_DEGs","DEGs_FC","S8495_resistant_DARs","DARs_FC")
  S8495_resistant_DARs<-as.character(different_FC[which(different_FC[,2]<0 & different_FC[,4]<0),3])
  
  homer_anno_up<-fread("sensitivity-resistance-S8495-up-annotation.txt",header=T,data.table=F)
  colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
  Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
  length(unique(Anno_gene_100Kb_up$`Gene Name`))
  
  up_gene<-read.table("sensitivity-resistance-S8495-DESeq2_AUC_RNAseq_P0.05_FC0_up_DEGs.txt",header=T,sep="\t")
  DEGs<-as.character(up_gene[,1])
  int_all_DEGs<-intersect(Anno_gene_100Kb_up$`Gene Name`,DEGs)
  
  all_peaks_anno<-homer_anno_up
  DApeak<-fread("sensitivity-resistance-all-DESeq2_S8495_AUC_ATAC.txt",header=T,data.table=F)
  
  different_FC<-NULL
  for(i in 1:length(int_all_DEGs)){
    DEGs_peaks1<-all_peaks_anno[all_peaks_anno$`Gene Name` %in% int_all_DEGs[i],1]
    DEGs_peaks_FC<-as.data.frame(DApeak[match(DEGs_peaks1,DApeak$peak_id),3])
    gene_name<-matrix(rep(int_all_DEGs[i],length(DEGs_peaks1)),ncol=1)
    DEGs_FC<-matrix(rep(up_gene[match(int_all_DEGs[i],DEGs),3],ncol=1))
    different_FC1<-cbind(gene_name,DEGs_FC,DEGs_peaks1,DEGs_peaks_FC)
    different_FC<-rbind(different_FC,different_FC1)
  }
  sum(different_FC[,2]>0) #up DEGs
  sum(different_FC[,4]>0) #up DARs
  colnames(different_FC)<-c("S8495_sensitive_DEGs","DEGs_FC","S8495_sensitive_DARs","DARs_FC")
  S8495_sensitive_DARs<-as.character(different_FC[which(different_FC[,2]>0 & different_FC[,4]>0),3])
  
  intersect(S8495_sensitive_DARs,S8495_resistant_DARs)
  a<-intersect(GEM_sensitive_DARs,S8495_sensitive_DARs)
  if(length(a)==0){
    GEM_sensitive_only1<-GEM_sensitive_DARs
    S8495_sensitive_only<-S8495_sensitive_DARs
  }
  if(length(a)!=0){
    GEM_sensitive_only1<-GEM_sensitive_DARs[-match(a,GEM_sensitive_DARs)]
    S8495_sensitive_only<-S8495_sensitive_DARs[-match(a,S8495_sensitive_DARs)]
  }
  
  b<-intersect(GEM_sensitive_only1,S8495_resistant_DARs)
  if(length(b)==0){
    GEM_sensitive_only<-GEM_sensitive_only1
    S8495_resistant_only<-S8495_resistant_DARs
  }
  if(length(b)!=0){
    GEM_sensitive_only<-GEM_sensitive_only1[-match(b,GEM_sensitive_only1)]
    S8495_resistant_only<-S8495_resistant_DARs[-match(b,S8495_resistant_DARs)]
  }
  DARs<-c(GEM_sensitive_only,S8495_sensitive_only,S8495_resistant_only)
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
  peak_RPKM37_DARs<-peak_RPKM37[match(DARs,peak_RPKM[,1]),]
  peak_RPKM37_DARs1<-t(peak_RPKM37_DARs)
  colnames(peak_RPKM37_DARs1)<-DARs
  y<-matrix(c(rep("GEM_sensitive",24),rep("GEM_res_S8495_sen",7),rep("GEM_res_S8495_res",6)),ncol=1)
  colnames(y)<-"y"
  data<-as.data.frame(cbind(peak_RPKM37_DARs1,y))
  text<-read.table("text.txt",sep = '\t',header=F)##数据输出
  text_data<-data[match(as.character(text[,1]),all),]
  train_data<-data[-(match(as.character(text[,1]),all)),]
  
  library(randomForest)
  rf_train<-randomForest(as.factor(y)~.,data=train_data,mtry=6,ntree=500,importance=TRUE,proximity=TRUE)
  importance<-importance(rf_train)
  yPred<- predict(rf_train,newdata=text_data)
  ConfM<-table(yPred,text_data$y)
  Err<-(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  
  #library(klaR)
  #fit_Bayes1=NaiveBayes(y~.,data=train_data)
  #yPred=predict(fit_Bayes1,text_data)
  #ConfM<-table(text_data$y,yPred$class)
  #Err=sum(as.numeric(as.numeric(yPred$class)!=as.numeric(text_data$y)))/nrow(text_data) 
  
  
  rf_result[k,1]<-k
  rf_result[k,2]<-paste(length(GEM_sensitive_only),length(S8495_sensitive_only),length(S8495_resistant_only),length(DARs),sep=";")
  rf_result[k,3]<-paste0(as.character(text_data[,ncol(text_data)]),collapse =",")
  rf_result[k,4]<-paste0(as.character(yPred),collapse =",")
  rf_result[k,5]<-Err
  rf_result[k,6]<-1-Err
}
colnames(rf_result)<-c("Times","length of DARs","text_label","predict_label","error","precision")
mean(as.numeric(as.character(rf_result[,6])))
setwd('~/xjj/drug/drug_result/HDACi_chemo618/classifier')
write.table(rf_result,"rf_result_50fold_peak_gene.txt",sep="\t",col.names=T,row.names=F)



###########  单纯只用peaks作为特征，没有涉及DEGs，两个路径需要改
###################################################
rm(list=ls())
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_RPKM<-read.table("peak_RPKM.txt",sep="\t",header=T) #没有全0的
colnames(peak_RPKM)<-gsub("\\.","-",colnames(peak_RPKM))

GEM_sensitive<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
                 "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
                 "DAC-19","DAC-37","DAC-39","DAC-38")
GEM_res_S8495_sen<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
GEM_res_S8495_res<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")

##median
GEM_sensitive<-c("DAC-24","DAC-27","DAC-28","DAC-4","DAC-33","DAC-34","DAC-36","DAC-6","DAC-1","DAC-9","DAC-11","DAC-15","DAC-3","DAC-16","DAC-18","DAC-19","DAC-37","DAC-38")
GEM_res_S8495_sen<-c("DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-12","DAC-14")
GEM_res_S8495_res<-c("DAC-20","DAC-21","DAC-22","DAC-23","DAC-25","DAC-29","DAC-30","DAC-32","DAC-35","DAC-5","DAC-13","DAC-39")
##


rf_result<-matrix(0,40,6)
for(k in 1:40){
  cat(k,"\n")
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier4/CV10',k,sep='_')
  GEM<-read.table(paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_GEM_AUC_ATAC.txt",sep=""),sep = '\t',header=T)##数据输出
  resSig<-GEM[which(GEM$pvalue<0.05 & abs(GEM$log2FoldChange)>2),]
  #resSig<-GEM[which(GEM$pvalue<0.01),]
  resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
  resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
  cat(sum(resSig$up_down=='up'),"\n")
  GEM_sensitive_DARs<-as.character(resSig[which(resSig$log2FoldChange>0),1])
  
  S8495<-read.table(paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_S8495_AUC_ATAC.txt",sep=""),sep = '\t',header=T)##数据输出
  resSig_S8495<-S8495[which(S8495$pvalue<0.05 & abs(S8495$log2FoldChange)>2),]
  #resSig<-GEM[which(GEM$pvalue<0.01),]
  resSig_S8495[which(resSig_S8495$log2FoldChange>0),'up_down']<-'up'
  resSig_S8495[which(resSig_S8495$log2FoldChange<0),'up_down']<-'down'
  cat(sum(resSig_S8495$up_down=='up'),"\n")
  cat(sum(resSig_S8495$up_down=='down'),"\n")
  S8495_sensitive_DARs<-as.character(resSig_S8495[which(resSig_S8495$log2FoldChange>0),1])
  S8495_resistant_DARs<-as.character(resSig_S8495[which(resSig_S8495$log2FoldChange<0),1])
  a<-intersect(GEM_sensitive_DARs,S8495_sensitive_DARs)
  if(length(a)==0){
    GEM_sensitive_only1<-GEM_sensitive_DARs
    S8495_sensitive_only<-S8495_sensitive_DARs
  }
  if(length(a)!=0){
    GEM_sensitive_only1<-GEM_sensitive_DARs[-match(a,GEM_sensitive_DARs)]
    S8495_sensitive_only<-S8495_sensitive_DARs[-match(a,S8495_sensitive_DARs)]
  }
  
  b<-intersect(GEM_sensitive_only1,S8495_resistant_DARs)
  if(length(b)==0){
    GEM_sensitive_only<-GEM_sensitive_only1
    S8495_resistant_only<-S8495_resistant_DARs
  }
  if(length(b)!=0){
    GEM_sensitive_only<-GEM_sensitive_only1[-match(b,GEM_sensitive_only1)]
    S8495_resistant_only<-S8495_resistant_DARs[-match(b,S8495_resistant_DARs)]
  }
  DARs<-c(GEM_sensitive_only,S8495_sensitive_only,S8495_resistant_only)
  all<-c(GEM_sensitive,GEM_res_S8495_sen,GEM_res_S8495_res)
  peak_RPKM37<-peak_RPKM[,match(all,colnames(peak_RPKM))]
  peak_RPKM37_DARs<-peak_RPKM37[match(DARs,peak_RPKM[,1]),]
  peak_RPKM37_DARs1<-t(peak_RPKM37_DARs)
  colnames(peak_RPKM37_DARs1)<-DARs
  y<-matrix(c(rep("GEM_sensitive",24),rep("GEM_res_S8495_sen",7),rep("GEM_res_S8495_res",6)),ncol=1)
  colnames(y)<-"y"
  data<-as.data.frame(cbind(peak_RPKM37_DARs1,y))
  setwd(paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier4/CV10',k,sep='_'))
  text<-read.table("text.txt",sep = '\t',header=F)##数据输出
  text_data<-data[match(as.character(text[,1]),all),]
  train_data<-data[-(match(as.character(text[,1]),all)),]
  #library(randomForest)
  #rf_train<-randomForest(as.factor(y)~.,data=train_data,mtry=6,ntree=500,importance=TRUE,proximity=TRUE)
  #importance<-importance(rf_train)
  #yPred<- predict(rf_train,newdata=text_data)
  #ConfM<-table(yPred,text_data$y)
  #Err<-(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  
  library(klaR)
  fit_Bayes1=NaiveBayes(y~.,data=train_data)
  yPred=predict(fit_Bayes1,text_data)
  ConfM<-table(text_data$y,yPred$class)
  Err=sum(as.numeric(as.numeric(yPred$class)!=as.numeric(text_data$y)))/nrow(text_data) 
  
  rf_result[k,1]<-k
  rf_result[k,2]<-paste(length(GEM_sensitive_only),length(S8495_sensitive_only),length(S8495_resistant_only),length(DARs),sep=";")
  rf_result[k,3]<-paste0(as.character(text_data[,ncol(text_data)]),collapse =",")
  rf_result[k,4]<-paste0(as.character(yPred$class),collapse =",")
  rf_result[k,5]<-Err
  rf_result[k,6]<-1-Err
}
colnames(rf_result)<-c("times","length of biomarks","text_label","predict_label","error","precision")
sum(as.numeric(as.character(rf_result[,6])))/nrow(rf_result)
setwd('~/xjj/drug/drug_result/HDACi_chemo618/classifier4')
write.table(rf_result,"NaiveBayes_result_40fold_peak_P0.05_FC2.txt",sep="\t",col.names=T,row.names=F)


###########################################################
### knn
library(caret)
knn.model1 = knn3(y~.,data = train_data, k = 3) 
knn.response1 = predict(knn.model1,text_data,type = "class") 
ConfM<-table(knn.response1,text_data$y)
Err=(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)

### xgboost
#devtools::install_github('dmlc/xgboost',subdir='R-package')
library(xgboost)
require(xgboost)
require(methods)
require(plyr)
library(Matrix)
library(caret)
# 训练集的数据预处理
# 将trainset的1-4列（自变量）转换为矩阵
traindata1 <- as.matrix(train_data[,c(1:(ncol(train_data)-1))])
# 利用Matrix函数，将sparse参数设置为TRUE，转化为稀疏矩阵
traindata2 <- Matrix(traindata1,sparse = T)
# 将因变量转换为numeric类型，-1是为了从0开始计数
train_y <- as.numeric(train_data[,ncol(train_data)])-1
# 将自变量和因变量拼接为list
traindata <- list(data=traindata2,label=train_y)
# 构造模型需要的xgb.DMatrix对象，处理对象为稀疏矩阵
dtrain <- xgb.DMatrix(data = traindata$data, label = traindata$label) 

# 测试集的数据预处理
# 将自变量转化为矩阵
testset1 <- as.matrix(text_data[,c(1:(ncol(text_data)-1))])
# 利用Matrix函数，将sparse参数设置为TRUE，转化为稀疏矩阵
testset2 <- Matrix(testset1,sparse=T) 
# 将因变量转化为numeric
test_y <- as.numeric(text_data[,ncol(text_data)])-1
# 将自变量和因变量拼接为list
testset <- list(data=testset2,label=test_y)
# 构造模型需要的xgb.DMatrix对象，处理对象为稀疏矩阵
dtest <- xgb.DMatrix(data=testset$data,label=testset$label)
# 建立模型
model_xgb <- xgboost(data=dtrain,booster='gbtree',max_depth=6,eta=0.5,objective='multi:softmax',num_class=3,nround=10)
#用测试集预测
pre <- predict(model_xgb,newdata=dtest)
table(pre,test_y)
accurancy4 = mean(pre ==test_y)
Err=1-accurancy4


rm(list=ls())
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_RPKM<-read.table("peak_RPKM.txt",sep="\t",header=T) #没有全0的
colnames(peak_RPKM)<-gsub("\\.","-",colnames(peak_RPKM))

GEM_sensitive<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
                 "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
                 "DAC-19","DAC-37","DAC-39","DAC-38")
GEM_res_S8495_sen<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
GEM_res_S8495_res<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")

##median
GEM_sensitive<-c("DAC-24","DAC-27","DAC-28","DAC-4","DAC-33","DAC-34","DAC-36","DAC-6","DAC-1","DAC-9","DAC-11","DAC-15","DAC-3","DAC-16","DAC-18","DAC-19","DAC-37","DAC-38")
GEM_res_S8495_sen<-c("DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-12","DAC-14")
GEM_res_S8495_res<-c("DAC-20","DAC-21","DAC-22","DAC-23","DAC-25","DAC-29","DAC-30","DAC-32","DAC-35","DAC-5","DAC-13","DAC-39")
##


rf_result<-matrix(0,40,6)
for(k in 1:40){
  cat(k,"\n")
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier4/CV10',k,sep='_')
  GEM<-read.table(paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_GEM_AUC_ATAC.txt",sep=""),sep = '\t',header=T)##数据输出
  resSig<-GEM[which(GEM$pvalue<0.05 & abs(GEM$log2FoldChange)>2),]
  #resSig<-GEM[which(GEM$pvalue<0.01),]
  resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
  resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
  cat(sum(resSig$up_down=='up'),"\n")
  GEM_sensitive_DARs<-as.character(resSig[which(resSig$log2FoldChange>0),1])
  
  S8495<-read.table(paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_S8495_AUC_ATAC.txt",sep=""),sep = '\t',header=T)##数据输出
  resSig_S8495<-S8495[which(S8495$pvalue<0.05 & abs(S8495$log2FoldChange)>2),]
  #resSig<-GEM[which(GEM$pvalue<0.01),]
  resSig_S8495[which(resSig_S8495$log2FoldChange>0),'up_down']<-'up'
  resSig_S8495[which(resSig_S8495$log2FoldChange<0),'up_down']<-'down'
  cat(sum(resSig_S8495$up_down=='up'),"\n")
  cat(sum(resSig_S8495$up_down=='down'),"\n")
  S8495_sensitive_DARs<-as.character(resSig_S8495[which(resSig_S8495$log2FoldChange>0),1])
  S8495_resistant_DARs<-as.character(resSig_S8495[which(resSig_S8495$log2FoldChange<0),1])
  a<-intersect(GEM_sensitive_DARs,S8495_sensitive_DARs)
  if(length(a)==0){
    GEM_sensitive_only1<-GEM_sensitive_DARs
    S8495_sensitive_only<-S8495_sensitive_DARs
  }
  if(length(a)!=0){
    GEM_sensitive_only1<-GEM_sensitive_DARs[-match(a,GEM_sensitive_DARs)]
    S8495_sensitive_only<-S8495_sensitive_DARs[-match(a,S8495_sensitive_DARs)]
  }
  
  b<-intersect(GEM_sensitive_only1,S8495_resistant_DARs)
  if(length(b)==0){
    GEM_sensitive_only<-GEM_sensitive_only1
    S8495_resistant_only<-S8495_resistant_DARs
  }
  if(length(b)!=0){
    GEM_sensitive_only<-GEM_sensitive_only1[-match(b,GEM_sensitive_only1)]
    S8495_resistant_only<-S8495_resistant_DARs[-match(b,S8495_resistant_DARs)]
  }
  DARs<-c(GEM_sensitive_only,S8495_sensitive_only,S8495_resistant_only)
  all<-c(GEM_sensitive,GEM_res_S8495_sen,GEM_res_S8495_res)
  peak_RPKM37<-peak_RPKM[,match(all,colnames(peak_RPKM))]
  peak_RPKM37_DARs<-peak_RPKM37[match(DARs,peak_RPKM[,1]),]
  peak_RPKM37_DARs1<-t(peak_RPKM37_DARs)
  colnames(peak_RPKM37_DARs1)<-DARs
  y<-matrix(c(rep("GEM_sensitive",24),rep("GEM_res_S8495_sen",7),rep("GEM_res_S8495_res",6)),ncol=1)
  colnames(y)<-"y"
  data<-as.data.frame(cbind(peak_RPKM37_DARs1,y))
  setwd(paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier4/CV10',k,sep='_'))
  text<-read.table("text.txt",sep = '\t',header=F)##数据输出
  text_data<-data[match(as.character(text[,1]),all),]
  train_data<-data[-(match(as.character(text[,1]),all)),]
  
  library(caret)
  knn.model1 = knn3(y~.,data = train_data, k = 15) 
  knn.response1 = predict(knn.model1,text_data,type = "class") #，type="class"表⽰结果为分类。
  ConfM<-table(knn.response1,text_data$y)
  Err=(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  
  rf_result[k,1]<-k
  rf_result[k,2]<-paste(length(GEM_sensitive_only),length(S8495_sensitive_only),length(S8495_resistant_only),length(DARs),sep=";")
  rf_result[k,3]<-paste0(as.character(text_data[,ncol(text_data)]),collapse =",")
  rf_result[k,4]<-paste0(as.character(knn.response1),collapse =",")
  rf_result[k,5]<-Err
  rf_result[k,6]<-1-Err
}
colnames(rf_result)<-c("times","length of biomarks","text_label","predict_label","error","precision")
sum(as.numeric(as.character(rf_result[,6])))/nrow(rf_result)
setwd('~/xjj/drug/drug_result/HDACi_chemo618/classifier4')
write.table(rf_result,"knn_result_40fold_peak_P0.05_FC2_k15.txt",sep="\t",col.names=T,row.names=F)


#######################################################################################
### xgboost
#devtools::install_github('dmlc/xgboost',subdir='R-package')
library(xgboost)
require(xgboost)
require(methods)
require(plyr)
library(Matrix)
library(caret)
# 训练集的数据预处理
# 将trainset的1-4列（自变量）转换为矩阵
traindata1 <- as.matrix(train_data[,c(1:(ncol(train_data)-1))])
# 利用Matrix函数，将sparse参数设置为TRUE，转化为稀疏矩阵
traindata2 <- Matrix(traindata1,sparse = T)
# 将因变量转换为numeric类型，-1是为了从0开始计数
train_y <- as.numeric(train_data[,ncol(train_data)])-1
# 将自变量和因变量拼接为list
traindata <- list(data=traindata2,label=train_y)
# 构造模型需要的xgb.DMatrix对象，处理对象为稀疏矩阵
dtrain <- xgb.DMatrix(data = traindata$data, label = traindata$label) 

# 测试集的数据预处理
# 将自变量转化为矩阵
testset1 <- as.matrix(text_data[,c(1:(ncol(text_data)-1))])
# 利用Matrix函数，将sparse参数设置为TRUE，转化为稀疏矩阵
testset2 <- Matrix(testset1,sparse=T) 
# 将因变量转化为numeric
test_y <- as.numeric(text_data[,ncol(text_data)])-1
# 将自变量和因变量拼接为list
testset <- list(data=testset2,label=test_y)
# 构造模型需要的xgb.DMatrix对象，处理对象为稀疏矩阵
dtest <- xgb.DMatrix(data=testset$data,label=testset$label)
# 建立模型
model_xgb <- xgboost(data=dtrain,booster='gbtree',max_depth=6,eta=0.5,objective='multi:softmax',num_class=3,nround=10)
#用测试集预测
pre <- predict(model_xgb,newdata=dtest)
table(pre,test_y)
accurancy4 = mean(pre ==test_y)
Err=1-accurancy4


rm(list=ls())
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_RPKM<-read.table("peak_RPKM.txt",sep="\t",header=T) #没有全0的
colnames(peak_RPKM)<-gsub("\\.","-",colnames(peak_RPKM))

GEM_sensitive<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
                 "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
                 "DAC-19","DAC-37","DAC-39","DAC-38")
GEM_res_S8495_sen<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
GEM_res_S8495_res<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")

##median
GEM_sensitive<-c("DAC-24","DAC-27","DAC-28","DAC-4","DAC-33","DAC-34","DAC-36","DAC-6","DAC-1","DAC-9","DAC-11","DAC-15","DAC-3","DAC-16","DAC-18","DAC-19","DAC-37","DAC-38")
GEM_res_S8495_sen<-c("DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-12","DAC-14")
GEM_res_S8495_res<-c("DAC-20","DAC-21","DAC-22","DAC-23","DAC-25","DAC-29","DAC-30","DAC-32","DAC-35","DAC-5","DAC-13","DAC-39")
##


rf_result<-matrix(0,40,6)
for(k in 1:40){
  cat(k,"\n")
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier4/CV10',k,sep='_')
  GEM<-read.table(paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_GEM_AUC_ATAC.txt",sep=""),sep = '\t',header=T)##数据输出
  resSig<-GEM[which(GEM$pvalue<0.05 & abs(GEM$log2FoldChange)>2),]
  #resSig<-GEM[which(GEM$pvalue<0.01),]
  resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
  resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
  cat(sum(resSig$up_down=='up'),"\n")
  GEM_sensitive_DARs<-as.character(resSig[which(resSig$log2FoldChange>0),1])
  
  S8495<-read.table(paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_S8495_AUC_ATAC.txt",sep=""),sep = '\t',header=T)##数据输出
  resSig_S8495<-S8495[which(S8495$pvalue<0.05 & abs(S8495$log2FoldChange)>2),]
  #resSig<-GEM[which(GEM$pvalue<0.01),]
  resSig_S8495[which(resSig_S8495$log2FoldChange>0),'up_down']<-'up'
  resSig_S8495[which(resSig_S8495$log2FoldChange<0),'up_down']<-'down'
  cat(sum(resSig_S8495$up_down=='up'),"\n")
  cat(sum(resSig_S8495$up_down=='down'),"\n")
  S8495_sensitive_DARs<-as.character(resSig_S8495[which(resSig_S8495$log2FoldChange>0),1])
  S8495_resistant_DARs<-as.character(resSig_S8495[which(resSig_S8495$log2FoldChange<0),1])
  a<-intersect(GEM_sensitive_DARs,S8495_sensitive_DARs)
  if(length(a)==0){
    GEM_sensitive_only1<-GEM_sensitive_DARs
    S8495_sensitive_only<-S8495_sensitive_DARs
  }
  if(length(a)!=0){
    GEM_sensitive_only1<-GEM_sensitive_DARs[-match(a,GEM_sensitive_DARs)]
    S8495_sensitive_only<-S8495_sensitive_DARs[-match(a,S8495_sensitive_DARs)]
  }
  
  b<-intersect(GEM_sensitive_only1,S8495_resistant_DARs)
  if(length(b)==0){
    GEM_sensitive_only<-GEM_sensitive_only1
    S8495_resistant_only<-S8495_resistant_DARs
  }
  if(length(b)!=0){
    GEM_sensitive_only<-GEM_sensitive_only1[-match(b,GEM_sensitive_only1)]
    S8495_resistant_only<-S8495_resistant_DARs[-match(b,S8495_resistant_DARs)]
  }
  DARs<-c(GEM_sensitive_only,S8495_sensitive_only,S8495_resistant_only)
  all<-c(GEM_sensitive,GEM_res_S8495_sen,GEM_res_S8495_res)
  peak_RPKM37<-peak_RPKM[,match(all,colnames(peak_RPKM))]
  peak_RPKM37_DARs<-peak_RPKM37[match(DARs,peak_RPKM[,1]),]
  peak_RPKM37_DARs1<-t(peak_RPKM37_DARs)
  colnames(peak_RPKM37_DARs1)<-DARs
  y<-matrix(c(rep("GEM_sensitive",24),rep("GEM_res_S8495_sen",7),rep("GEM_res_S8495_res",6)),ncol=1)
  colnames(y)<-"y"
  data<-as.data.frame(cbind(peak_RPKM37_DARs1,y))
  setwd(paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier4/CV10',k,sep='_'))
  text<-read.table("text.txt",sep = '\t',header=F)##数据输出
  text_data<-data[match(as.character(text[,1]),all),]
  train_data<-data[-(match(as.character(text[,1]),all)),]
  
  library(xgboost)
  require(xgboost)
  require(methods)
  require(plyr)
  library(Matrix)
  library(caret)
  traindata1 <- as.matrix(train_data[,c(1:(ncol(train_data)-1))])
  traindata2 <- Matrix(traindata1,sparse = T)
  train_y <- as.numeric(train_data[,ncol(train_data)])-1
  traindata <- list(data=traindata2,label=train_y)
  dtrain <- xgb.DMatrix(data = traindata$data, label = traindata$label) 
  testset1 <- as.matrix(text_data[,c(1:(ncol(text_data)-1))])
  testset2 <- Matrix(testset1,sparse=T) 
  test_y <- as.numeric(text_data[,ncol(text_data)])-1
  testset <- list(data=testset2,label=test_y)
  dtest <- xgb.DMatrix(data=testset$data,label=testset$label)
  model_xgb <- xgboost(data=dtrain,booster='gbtree',max_depth=6,eta=0.5,objective='multi:softmax',num_class=3,nround=10)
  pre <- predict(model_xgb,newdata=dtest)
  table(pre,test_y)
  accurancy4 = mean(pre ==test_y)
  Err=1-accurancy4
  
  
  rf_result[k,1]<-k
  rf_result[k,2]<-paste(length(GEM_sensitive_only),length(S8495_sensitive_only),length(S8495_resistant_only),length(DARs),sep=";")
  rf_result[k,3]<-paste0(as.character(text_data[,ncol(text_data)]),collapse =",")
  rf_result[k,4]<-paste0(as.character(pre),collapse =",")
  rf_result[k,5]<-Err
  rf_result[k,6]<-accurancy4
}
colnames(rf_result)<-c("times","length of biomarks","text_label","predict_label","error","precision")
sum(as.numeric(as.character(rf_result[,6])))/nrow(rf_result)
setwd('~/xjj/drug/drug_result/HDACi_chemo618/classifier4')
write.table(rf_result,"xgboost_result_40fold_peak_P0.05_FC2.txt",sep="\t",col.names=T,row.names=F)






###################################################################################
##############    三分类     #####################
#########################  留一法构造 ########################################################
#SVM 独立训练 独立验证 leave one out
rm(list=ls())
library(DESeq2)
library("data.table")
library(dplyr)
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_count_name<-read.table("peak_count.txt",sep="\t",header=T) #没有全0的
colnames(peak_count_name)<-gsub("\\.","-",colnames(peak_count_name))
text_leave_one_out3<-matrix(0,40,2)
library(doParallel)
library(parallel)
library(iterators)

for(i in 31:37){
  cat(i,"\n")
  GEM_sensitive<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
                   "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
                   "DAC-19","DAC-37","DAC-39","DAC-38")
  GEM_res_S8495_sen<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
  GEM_res_S8495_res<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")
  all<-c(GEM_sensitive,GEM_res_S8495_sen,GEM_res_S8495_res)
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier1/classifier',i,sep='_')
  text1<-read.table(paste(outputPath,'/',"text.txt",sep=""),header=F,sep="\t")
  text<-text1[1,1]
  
  text_leave_one_out3[i,1]<-i
  text_leave_one_out3[i,2]<-text
  if(sum(GEM_sensitive==text)==1){
    text_label<-"GEM_sensitive"
    GEM_sensitive<-GEM_sensitive[-match(text,GEM_sensitive)]
  }
  if(sum(GEM_res_S8495_sen==text)==1){
    text_label<-"GEM_res_S8495_sen"
    GEM_res_S8495_sen<-GEM_res_S8495_sen[-match(text,GEM_res_S8495_sen)]
  }
  if(sum(GEM_res_S8495_res==text)==1){
    text_label<-"GEM_res_S8495_res"
    GEM_res_S8495_res<-GEM_res_S8495_res[-match(text,GEM_res_S8495_res)]
  }
  sensitivity<-GEM_res_S8495_sen
  resistance<-c(GEM_sensitive,GEM_res_S8495_res)
  resistance_peak<-peak_count_name[,match(resistance,colnames(peak_count_name))]
  sensitivity_peak<-peak_count_name[,match(sensitivity,colnames(peak_count_name))]
  peak_name<-peak_count_name[,1]
  colDate<-data.frame(row.names = c(as.vector(resistance),as.vector(sensitivity)),
                      condition=factor(c(rep("resistance",length(resistance)),rep("sensitivity",length(sensitivity))))
  )
  datexpr<-cbind(resistance_peak,sensitivity_peak)##前面的相对于后面的上下调
  counts <- apply(datexpr,2,as.numeric)   ###矩阵中必须是数值
  dds<-DESeqDataSetFromMatrix(countData = counts,colData = colDate,design = ~condition)
  dds<-DESeq(dds)##进行标准化分析
  res<-results(dds)##将结果输出
  res<-as.data.frame(res)
  res<-cbind(peak_count_name[,1],res)  ##对数据增加一列
  colnames(res)<- c('peak_id',"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj" )##对列名进行修改
  
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier1/classifier',i,sep='_')
  if (!file.exists(outputPath)){
    dir.create(outputPath, recursive=TRUE)
  }
  write.table(res,paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_GEM_res_S8495_sen_AUC_ATAC.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  
  sensitivity1<-GEM_res_S8495_res
  resistance1<-c(GEM_res_S8495_sen,GEM_sensitive)
  resistance_peak1<-peak_count_name[,match(resistance1,colnames(peak_count_name))]
  sensitivity_peak1<-peak_count_name[,match(sensitivity1,colnames(peak_count_name))]
  peak_name<-peak_count_name[,1]
  colDate<-data.frame(row.names = c(as.vector(resistance1),as.vector(sensitivity1)),
                      condition=factor(c(rep("resistance1",length(resistance1)),rep("sensitivity1",length(sensitivity1))))
  )
  datexpr1<-cbind(resistance_peak1,sensitivity_peak1)##前面的相对于后面的上下调
  counts1 <- apply(datexpr1,2,as.numeric)   ###矩阵中必须是数值
  dds1<-DESeqDataSetFromMatrix(countData = counts1,colData = colDate,design = ~condition)
  dds1<-DESeq(dds1)##进行标准化分析
  res1<-results(dds1)##将结果输出
  res1<-as.data.frame(res1)
  res1<-cbind(peak_count_name[,1],res1)  ##对数据增加一列
  colnames(res1)<- c('peak_id',"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj" )##对列名进行修改
  
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier1/classifier',i,sep='_')
  if (!file.exists(outputPath)){
    dir.create(outputPath, recursive=TRUE)
  }
  write.table(res1,paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_GEM_res_S8495_res_AUC_ATAC.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
}
stopCluster(cl)
#######################  留一法验证 ############################
rm(list=ls())
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_RPKM<-read.table("peak_RPKM.txt",sep="\t",header=T) #没有全0的
colnames(peak_RPKM)<-gsub("\\.","-",colnames(peak_RPKM))

GEM_sensitive<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
                 "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
                 "DAC-19","DAC-37","DAC-39","DAC-38")
GEM_res_S8495_sen<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
GEM_res_S8495_res<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")

svm_result<-NULL
for(i in 1:37){
  cat(i,"\n")
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier1/classifier',i,sep='_')
  GEM<-read.table(paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_GEM_AUC_ATAC.txt",sep=""),sep = '\t',header=T)##数据输出
  resSig<-GEM[which(GEM$pvalue<0.05 & abs(GEM$log2FoldChange)>2),]
  #resSig<-GEM[which(GEM$pvalue<0.01),]
  resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
  resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
  cat(sum(resSig$up_down=='up'),"\n")
  GEM_sensitive_DARs<-as.character(resSig[which(resSig$log2FoldChange>0),1])
  
  S8495<-read.table(paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_S8495_AUC_ATAC.txt",sep=""),sep = '\t',header=T)##数据输出
  resSig_S8495<-S8495[which(S8495$pvalue<0.05 & abs(S8495$log2FoldChange)>2),]
  #resSig<-GEM[which(GEM$pvalue<0.01),]
  resSig_S8495[which(resSig_S8495$log2FoldChange>0),'up_down']<-'up'
  resSig_S8495[which(resSig_S8495$log2FoldChange<0),'up_down']<-'down'
  cat(sum(resSig_S8495$up_down=='up'),"\n")
  cat(sum(resSig_S8495$up_down=='down'),"\n")
  S8495_sensitive_DARs<-as.character(resSig_S8495[which(resSig_S8495$log2FoldChange>0),1])
  S8495_resistant_DARs<-as.character(resSig_S8495[which(resSig_S8495$log2FoldChange<0),1])
  a<-intersect(GEM_sensitive_DARs,S8495_sensitive_DARs)
  if(length(a)==0){
    GEM_sensitive_only1<-GEM_sensitive_DARs
    S8495_sensitive_only<-S8495_sensitive_DARs
  }
  if(length(a)!=0){
    GEM_sensitive_only1<-GEM_sensitive_DARs[-match(a,GEM_sensitive_DARs)]
    S8495_sensitive_only<-S8495_sensitive_DARs[-match(a,S8495_sensitive_DARs)]
  }
  
  b<-intersect(GEM_sensitive_only1,S8495_resistant_DARs)
  if(length(b)==0){
    GEM_sensitive_only<-GEM_sensitive_only1
    S8495_resistant_only<-S8495_resistant_DARs
  }
  if(length(b)!=0){
    GEM_sensitive_only<-GEM_sensitive_only1[-match(b,GEM_sensitive_only1)]
    S8495_resistant_only<-S8495_resistant_DARs[-match(b,S8495_resistant_DARs)]
  }
  DARs<-c(GEM_sensitive_only,S8495_sensitive_only,S8495_resistant_only)
  all<-c(GEM_sensitive,GEM_res_S8495_sen,GEM_res_S8495_res)
  cat(length(DARs),"\n")
  peak_RPKM37<-peak_RPKM[,match(all,colnames(peak_RPKM))]
  peak_RPKM37_DARs<-peak_RPKM37[match(DARs,peak_RPKM[,1]),]
  peak_RPKM37_DARs1<-t(peak_RPKM37_DARs)
  y<-matrix(c(rep("GEM_sensitive",24),rep("GEM_res_S8495_sen",7),rep("GEM_res_S8495_res",6)),ncol=1)
  colnames(y)<-"y"
  data<-as.data.frame(cbind(peak_RPKM37_DARs1,y))
  text<-read.table(paste(outputPath,'/',"text.txt",sep=""),sep = '\t',header=F)##数据输出
  library(e1071)
  set.seed(12345)
  text_data<-data[match(as.character(text[1,1]),all),]
  train_data<-data[-match(as.character(text[1,1]),all),]
  SvmFit<-svm(y~.,data=train_data,type="C-classification",kernel="radial",cost=5,gamma=1,scale=FALSE)
  head(SvmFit$decision.values)
  yPred<-predict(SvmFit,text_data)
  ConfM<-table(yPred,text_data$y)
  Err<-(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
  cat(Err,"\n")
  if(sum(GEM_sensitive==as.character(text[1,1]))==1){
    text_label<-"GEM_sensitive"
  }
  if(sum(GEM_res_S8495_sen==as.character(text[1,1]))==1){
    text_label<-"GEM_res_S8495_sen"
  }
  if(sum(GEM_res_S8495_res==as.character(text[1,1]))==1){
    text_label<-"GEM_res_S8495_res"
  }
  svm_result1<-as.data.frame(matrix(c(i,length(DARs),text_label,as.character(yPred),Err,1-Err),nrow=1))
  svm_result<-rbind(svm_result,svm_result1)
}
colnames(svm_result)<-c("times","length of biomarks","text_label","predict_label","error","precision")
sum(as.numeric(as.character(svm_result[,6])))/nrow(svm_result)

sum(as.numeric(as.character(svm_result[,6])))
setwd("~/xjj/drug/drug_result/HDACi_chemo618/classifier1")
write.table(svm_result,"svm_result_leave-one-out.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出


########################  差异peaks做注释
for(i in 28:37){
  cat(i,"\n")
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier',i,sep='_')
  GEM<-read.table(paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_GEM_AUC_ATAC.txt",sep=""),sep = '\t',header=T)##数据输出
  #resSig<-GEM[which(GEM$pvalue<0.05 & abs(GEM$log2FoldChange)>2),]
  resSig<-GEM[which(GEM$pvalue<0.01),]
  
  # pvalue padj
  ##新增一列，将log2FoldChange>0标注为up，<0标准为down
  resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
  resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
  up_gene<-as.data.frame(resSig[which(resSig$log2FoldChange>0),1])
  down_gene<-as.data.frame(resSig[which(resSig$log2FoldChange<0),1])
  library(dplyr)
  a_up<-apply(up_gene,1,function(x) unlist(strsplit(as.character(x), "_")))
  d_up<-t(a_up)
  e_up<-data.frame(rep("+",nrow(d_up)))
  peaks_up<-cbind(up_gene,d_up,e_up)
  colnames(peaks_up)<-c("peak_id","chr","start","end","strand")
  setwd(paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier',i,sep='_'))
  out_file<-"sensitivity-resistance-GEM-AUC-P0.01-logFC0"
  write.table(peaks_up,file=paste(out_file,"-DESeq2-up.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  a_down<-apply(down_gene,1,function(x) unlist(strsplit(as.character(x), "_")))
  d_down<-t(a_down)
  e_down<-data.frame(rep("+",nrow(d_down)))
  peaks_down<-cbind(down_gene,d_down,e_down)
  colnames(peaks_down)<-c("peak_id","chr","start","end","strand")
  write.table(peaks_down,file=paste(out_file,"-DESeq2-down.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
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
  
  
  S8495<-read.table(paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_S8495_AUC_ATAC.txt",sep=""),sep = '\t',header=T)##数据输出
  resSig1<-S8495[which(S8495$pvalue<0.01),]
  # pvalue padj
  ##新增一列，将log2FoldChange>0标注为up，<0标准为down
  resSig1[which(resSig1$log2FoldChange>0),'up_down']<-'up'
  resSig1[which(resSig1$log2FoldChange<0),'up_down']<-'down'
  up_gene1<-as.data.frame(resSig1[which(resSig1$log2FoldChange>0),1])
  down_gene1<-as.data.frame(resSig1[which(resSig1$log2FoldChange<0),1])
  library(dplyr)
  a_up<-apply(up_gene1,1,function(x) unlist(strsplit(as.character(x), "_")))
  d_up<-t(a_up)
  e_up<-data.frame(rep("+",nrow(d_up)))
  peaks_up1<-cbind(up_gene1,d_up,e_up)
  colnames(peaks_up1)<-c("peak_id","chr","start","end","strand")
  setwd(paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier',i,sep='_'))
  out_file1<-"sensitivity-resistance-S8495-AUC-P0.01-logFC0"
  write.table(peaks_up1,file=paste(out_file1,"-DESeq2-up.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  a_down<-apply(down_gene1,1,function(x) unlist(strsplit(as.character(x), "_")))
  d_down<-t(a_down)
  e_down<-data.frame(rep("+",nrow(d_down)))
  peaks_down1<-cbind(down_gene1,d_down,e_down)
  colnames(peaks_down1)<-c("peak_id","chr","start","end","strand")
  write.table(peaks_down1,file=paste(out_file1,"-DESeq2-down.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  #####  bed
  f_up<-data.frame(rep(1,nrow(d_up)))
  peaks_bed_up<-cbind(d_up,up_gene1,f_up,e_up)
  colnames(peaks_bed_up)<-c("chr","start","end","peak_id","score","strand")
  write.table(peaks_bed_up,file=paste(out_file1,"-DESeq2-up.bed",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  f_down<-data.frame(rep(1,nrow(d_down)))
  peaks_bed_down<-cbind(d_down,down_gene1,f_down,e_down)
  colnames(peaks_bed_down)<-c("chr","start","end","peak_id","score","strand")
  write.table(peaks_bed_down,file=paste(out_file1,"-DESeq2-down.bed",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
}

result1<-c()
for(i in 1:30){
  cat(i,"\n")
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier',i,sep='_')
  text<-read.table(paste(outputPath,'/',"text.txt",sep=""),sep = '\t',header=F)##数据输出
  result1<-c(result1,as.character(text[1,1]))
}  

match(unique(result1),result1)
unique(result1)

result2<-c("DAC-31","DAC-35","DAC-13","DAC-14","DAC-29","DAC-18","DAC-25","DAC-7","DAC-27","DAC-22","DAC-32","DAC-9","DAC-37","DAC-16","DAC-26","DAC-1","DAC-36","DAC-28","DAC-12","DAC-19","DAC-21")



#####################三分类十倍交叉验证
#SVM 独立训练 独立验证 10 fold cross validation
rm(list=ls())
library(DESeq2)
library("data.table")
library(dplyr)

setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_count_name<-read.table("peak_count.txt",sep="\t",header=T) #没有全0的
colnames(peak_count_name)<-gsub("\\.","-",colnames(peak_count_name))

setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
library("data.table")
exp1<-fread("star_rsem.GeneSymbol.readscount.xls",header=T,data.table=F)
exp2<-floor(exp1[,-1])
delete<-apply(exp2,1,function(x) mean(x==0))
cou<-which(delete>0.75)####在超过75%样本以上的都是0的基因删去
exp<-exp2[(-cou),]
geneid<-exp1[(-cou),1]

library(doParallel)
library(parallel)
library(iterators)

cl <- makeCluster(15)
registerDoParallel(cl)

for(i in 11:20){
  cat(i,"\n")
  #GEM_sensitive1<-c("DAC-24","DAC-27","DAC-28","DAC-4","DAC-33","DAC-34","DAC-36","DAC-6","DAC-1","DAC-9","DAC-11","DAC-15","DAC-3","DAC-16","DAC-18","DAC-19","DAC-37","DAC-38")
  #GEM_res_S8495_sen1<-c("DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-12","DAC-14")
  #GEM_res_S8495_res1<-c("DAC-20","DAC-21","DAC-22","DAC-23","DAC-25","DAC-29","DAC-30","DAC-32","DAC-35","DAC-5","DAC-13","DAC-39")
  
  
  GEM_sensitive1<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
                    "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
                    "DAC-19","DAC-37","DAC-39","DAC-38")
  GEM_res_S8495_sen1<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
  GEM_res_S8495_res1<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")
  all1<-c(GEM_sensitive1,GEM_res_S8495_sen1,GEM_res_S8495_res1)
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier/CV10',i,sep='_')
  text1<-read.table(paste(outputPath,'/',"text.txt",sep=""),header=F,sep="\t")
  text<-text1[,1]
  all<-all1[-match(text,all1)]
  GEM_sensitive<-intersect(GEM_sensitive1,all)
  GEM_res_S8495_sen<-intersect(GEM_res_S8495_sen1,all)
  GEM_res_S8495_res<-intersect(GEM_res_S8495_res1,all)
  
  sensitivity<-GEM_res_S8495_sen
  resistance<-c(GEM_sensitive,GEM_res_S8495_res)
  resistance_peak<-peak_count_name[,match(resistance,colnames(peak_count_name))]
  sensitivity_peak<-peak_count_name[,match(sensitivity,colnames(peak_count_name))]
  peak_name<-peak_count_name[,1]
  colDate<-data.frame(row.names = c(as.vector(resistance),as.vector(sensitivity)),
                      condition=factor(c(rep("resistance",length(resistance)),rep("sensitivity",length(sensitivity))))
  )
  datexpr<-cbind(resistance_peak,sensitivity_peak)##前面的相对于后面的上下调
  counts <- apply(datexpr,2,as.numeric)   ###矩阵中必须是数值
  library(DESeq2)
  dds<-DESeqDataSetFromMatrix(countData = counts,colData = colDate,design = ~condition)
  dds<-DESeq(dds)##进行标准化分析
  res<-results(dds)##将结果输出
  res<-as.data.frame(res)
  res<-cbind(peak_count_name[,1],res)  ##对数据增加一列
  colnames(res)<- c('peak_id',"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj" )##对列名进行修改
  
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier/CV10',i,sep='_')
  if (!file.exists(outputPath)){
    dir.create(outputPath, recursive=TRUE)
  }
  write.table(res,paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_GEM_res_S8495_sen_AUC_ATAC.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  resSig<-res[which(res$pvalue<0.01),]
  # pvalue padj
  ##新增一列，将log2FoldChange>0标注为up，<0标准为down
  resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
  resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
  up_gene<-as.data.frame(resSig[which(resSig$log2FoldChange>0),1])
  down_gene<-as.data.frame(resSig[which(resSig$log2FoldChange<0),1])
  a_up<-apply(up_gene,1,function(x) unlist(strsplit(as.character(x), "_")))
  d_up<-t(a_up)
  e_up<-data.frame(rep("+",nrow(d_up)))
  peaks_up<-cbind(up_gene,d_up,e_up)
  colnames(peaks_up)<-c("peak_id","chr","start","end","strand")
  setwd(paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier/CV10',i,sep='_'))
  out_file<-"sensitivity-resistance-GEM_res_S8495_sen-AUC-P0.01-logFC0"
  write.table(peaks_up,file=paste(out_file,"-DESeq2-up.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  a_down<-apply(down_gene,1,function(x) unlist(strsplit(as.character(x), "_")))
  d_down<-t(a_down)
  e_down<-data.frame(rep("+",nrow(d_down)))
  peaks_down<-cbind(down_gene,d_down,e_down)
  colnames(peaks_down)<-c("peak_id","chr","start","end","strand")
  write.table(peaks_down,file=paste(out_file,"-DESeq2-down.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
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
  
  resistance_exp<-exp[,match(resistance,colnames(exp))]
  sensitivity_exp<-exp[,match(sensitivity,colnames(exp))]
  colDate<-data.frame(row.names = c(as.vector(resistance),as.vector(sensitivity)),
                      condition=factor(c(rep("resistance",length(resistance)),rep("sensitivity",length(sensitivity))))
  )
  datexpr<-cbind(resistance_exp,sensitivity_exp)
  counts <- apply(datexpr,2,as.numeric)###矩阵中必须是数值
  dds<-DESeqDataSetFromMatrix(countData = counts,colData = colDate,design = ~condition)
  dds<-DESeq(dds)##进行标准化分析
  res<-results(dds)##将结果输出
  res<-as.data.frame(res)
  res<-cbind(geneid,res)  ##对数据增加一列
  colnames(res)<- c('gene_id',"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj" )##对列名进行修改
  write.table(res,paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_GEM_res_S8495_sen_AUC_RNAseq_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  resSig<-res[which(res$pvalue<0.05),]
  resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
  resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
  up_gene<-resSig[which(resSig$log2FoldChange>0),]
  down_gene<-resSig[which(resSig$log2FoldChange<0),]
  write.table(up_gene,paste(outputPath,'/',"sensitivity-resistance-GEM_res_S8495_sen-DESeq2_AUC_RNAseq_P0.05_FC0_up_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  write.table(down_gene,paste(outputPath,'/',"sensitivity-resistance-GEM_res_S8495_sen-DESeq2_AUC_RNAseq_P0.05_FC0_down_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  
  
  sensitivity1<-GEM_res_S8495_res
  resistance1<-c(GEM_res_S8495_sen,GEM_sensitive)
  resistance_peak1<-peak_count_name[,match(resistance1,colnames(peak_count_name))]
  sensitivity_peak1<-peak_count_name[,match(sensitivity1,colnames(peak_count_name))]
  peak_name<-peak_count_name[,1]
  colDate1<-data.frame(row.names = c(as.vector(resistance1),as.vector(sensitivity1)),
                       condition=factor(c(rep("resistance1",length(resistance1)),rep("sensitivity1",length(sensitivity1))))
  )
  datexpr1<-cbind(resistance_peak1,sensitivity_peak1)##前面的相对于后面的上下调
  counts1 <- apply(datexpr1,2,as.numeric)   ###矩阵中必须是数值
  dds1<-DESeqDataSetFromMatrix(countData = counts1,colData = colDate1,design = ~condition)
  dds1<-DESeq(dds1)##进行标准化分析
  res1<-results(dds1)##将结果输出
  res1<-as.data.frame(res1)
  res1<-cbind(peak_count_name[,1],res1)  ##对数据增加一列
  colnames(res1)<- c('peak_id',"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj" )##对列名进行修改
  write.table(res1,paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_GEM_res_S8495_res_AUC_ATAC.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  resSig1<-res1[which(res1$pvalue<0.01),]
  # pvalue padj
  ##新增一列，将log2FoldChange>0标注为up，<0标准为down
  resSig1[which(resSig1$log2FoldChange>0),'up_down']<-'up'
  resSig1[which(resSig1$log2FoldChange<0),'up_down']<-'down'
  up_gene1<-as.data.frame(resSig1[which(resSig1$log2FoldChange>0),1])
  down_gene1<-as.data.frame(resSig1[which(resSig1$log2FoldChange<0),1])
  library(dplyr)
  a_up<-apply(up_gene1,1,function(x) unlist(strsplit(as.character(x), "_")))
  d_up<-t(a_up)
  e_up<-data.frame(rep("+",nrow(d_up)))
  peaks_up1<-cbind(up_gene1,d_up,e_up)
  colnames(peaks_up1)<-c("peak_id","chr","start","end","strand")
  setwd(paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier/CV10',i,sep='_'))
  out_file1<-"sensitivity-resistance-GEM_res_S8495_res-AUC-P0.01-logFC0"
  write.table(peaks_up1,file=paste(out_file1,"-DESeq2-up.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  a_down<-apply(down_gene1,1,function(x) unlist(strsplit(as.character(x), "_")))
  d_down<-t(a_down)
  e_down<-data.frame(rep("+",nrow(d_down)))
  peaks_down1<-cbind(down_gene1,d_down,e_down)
  colnames(peaks_down1)<-c("peak_id","chr","start","end","strand")
  write.table(peaks_down1,file=paste(out_file1,"-DESeq2-down.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  #####  bed
  f_up<-data.frame(rep(1,nrow(d_up)))
  peaks_bed_up<-cbind(d_up,up_gene1,f_up,e_up)
  colnames(peaks_bed_up)<-c("chr","start","end","peak_id","score","strand")
  write.table(peaks_bed_up,file=paste(out_file1,"-DESeq2-up.bed",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  f_down<-data.frame(rep(1,nrow(d_down)))
  peaks_bed_down<-cbind(d_down,down_gene1,f_down,e_down)
  colnames(peaks_bed_down)<-c("chr","start","end","peak_id","score","strand")
  write.table(peaks_bed_down,file=paste(out_file1,"-DESeq2-down.bed",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  
  resistance_exp1<-exp[,match(resistance1,colnames(exp))]
  sensitivity_exp1<-exp[,match(sensitivity1,colnames(exp))]
  colDate1<-data.frame(row.names = c(as.vector(resistance1),as.vector(sensitivity1)),
                       condition=factor(c(rep("resistance1",length(resistance1)),rep("sensitivity1",length(sensitivity1))))
  )
  datexpr1<-cbind(resistance_exp1,sensitivity_exp1)
  counts1 <- apply(datexpr1,2,as.numeric)###矩阵中必须是数值
  dds1<-DESeqDataSetFromMatrix(countData = counts1,colData = colDate1,design = ~condition)
  dds1<-DESeq(dds1)##进行标准化分析
  res1<-results(dds1)##将结果输出
  res1<-as.data.frame(res1)
  res1<-cbind(geneid,res1)  ##对数据增加一列
  colnames(res1)<- c('gene_id',"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj" )##对列名进行修改
  write.table(res1,paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_GEM_res_S8495_res_AUC_RNAseq_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  resSig1<-res1[which(res1$pvalue<0.05),]
  resSig1[which(resSig1$log2FoldChange>0),'up_down']<-'up'
  resSig1[which(resSig1$log2FoldChange<0),'up_down']<-'down'
  up_gene1<-resSig1[which(resSig1$log2FoldChange>0),]
  down_gene1<-resSig1[which(resSig1$log2FoldChange<0),]
  write.table(up_gene1,paste(outputPath,'/',"sensitivity-resistance-GEM_res_S8495_res-DESeq2_AUC_RNAseq_P0.05_FC0_up_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  write.table(down_gene1,paste(outputPath,'/',"sensitivity-resistance-GEM_res_S8495_res-DESeq2_AUC_RNAseq_P0.05_FC0_down_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  
}

stopCluster(cl)

####################################################################
##################### 十倍交叉验证 #############
rm(list=ls())
library(plyr)
CVgroup <- function(k,datasize,seed){
  cvlist <- list()
  set.seed(seed)
  n <- rep(1:k,ceiling(datasize/k))[1:datasize]  #将数据分成K份，并生成的完成数据集n
  temp <- sample(n,datasize)  #把n打乱
  x <- 1:k
  dataseq <- 1:datasize
  cvlist <- lapply(x,function(x) dataseq[temp==x]) #dataseq中随机生成k个随机有序数据列
  return(cvlist)
}
k <- 10
datasize <- 37
cvlist <- CVgroup(k = k,datasize = datasize,seed = 450)#1,1206;2,50;3,100;4,150;5,200；6,250；7,300;8,350;9,400;10,450
cvlist

library(DESeq2)
library("data.table")
library(dplyr)

setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_count_name<-read.table("peak_count.txt",sep="\t",header=T) #没有全0的
colnames(peak_count_name)<-gsub("\\.","-",colnames(peak_count_name))

library(doParallel)
library(parallel)
library(iterators)

cl <- makeCluster(20)
registerDoParallel(cl)

for(i in 91:100){
  cat(i,"\n")
  GEM_sensitive1<-c("DAC-20","DAC-24","DAC-27","DAC-28","DAC-29","DAC-4","DAC-32","DAC-33","DAC-34","DAC-36",
                    "DAC-5","DAC-6","DAC-1","DAC-9","DAC-11","DAC-12","DAC-15","DAC-3","DAC-16","DAC-18",
                    "DAC-19","DAC-37","DAC-39","DAC-38")
  GEM_res_S8495_sen1<-c("DAC-23","DAC-26","DAC-31","DAC-7","DAC-8","DAC-2","DAC-14")
  GEM_res_S8495_res1<-c("DAC-21","DAC-22","DAC-25","DAC-30","DAC-35","DAC-13")
  all1<-c(GEM_sensitive1,GEM_res_S8495_sen1,GEM_res_S8495_res1)
  text<-all1[cvlist[[i-90]]]
  
  all<-all1[-match(text,all1)]
  GEM_sensitive<-intersect(GEM_sensitive1,all)
  GEM_res_S8495_sen<-intersect(GEM_res_S8495_sen1,all)
  GEM_res_S8495_res<-intersect(GEM_res_S8495_res1,all)
  
  sensitivity<-GEM_sensitive
  resistance<-c(GEM_res_S8495_sen,GEM_res_S8495_res)
  resistance_peak<-peak_count_name[,match(resistance,colnames(peak_count_name))]
  sensitivity_peak<-peak_count_name[,match(sensitivity,colnames(peak_count_name))]
  peak_name<-peak_count_name[,1]
  colDate<-data.frame(row.names = c(as.vector(resistance),as.vector(sensitivity)),
                      condition=factor(c(rep("resistance",length(resistance)),rep("sensitivity",length(sensitivity))))
  )
  datexpr<-cbind(resistance_peak,sensitivity_peak)##前面的相对于后面的上下调
  counts <- apply(datexpr,2,as.numeric)   ###矩阵中必须是数值
  dds<-DESeqDataSetFromMatrix(countData = counts,colData = colDate,design = ~condition)
  dds<-DESeq(dds)##进行标准化分析
  res<-results(dds)##将结果输出
  res<-as.data.frame(res)
  res<-cbind(peak_count_name[,1],res)  ##对数据增加一列
  colnames(res)<- c('peak_id',"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj" )##对列名进行修改
  
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/10fold-CV/CV10',i,sep='_')
  if (!file.exists(outputPath)){
    dir.create(outputPath, recursive=TRUE)
  }
  write.table(res,paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_GEM_AUC_ATAC.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  write.table(text,paste(outputPath,'/',"text.txt",sep=""),sep = '\t',col.names = F,row.names = F,quote = FALSE)##数据输出
  sensitivity1<-GEM_res_S8495_sen
  resistance1<-GEM_res_S8495_res
  resistance_peak1<-peak_count_name[,match(resistance1,colnames(peak_count_name))]
  sensitivity_peak1<-peak_count_name[,match(sensitivity1,colnames(peak_count_name))]
  peak_name<-peak_count_name[,1]
  colDate1<-data.frame(row.names = c(as.vector(resistance1),as.vector(sensitivity1)),
                       condition=factor(c(rep("resistance1",length(resistance1)),rep("sensitivity1",length(sensitivity1))))
  )
  datexpr1<-cbind(resistance_peak1,sensitivity_peak1)##前面的相对于后面的上下调
  counts1 <- apply(datexpr1,2,as.numeric)   ###矩阵中必须是数值
  dds1<-DESeqDataSetFromMatrix(countData = counts1,colData = colDate1,design = ~condition)
  dds1<-DESeq(dds1)##进行标准化分析
  res1<-results(dds1)##将结果输出
  res1<-as.data.frame(res1)
  res1<-cbind(peak_count_name[,1],res1)  ##对数据增加一列
  colnames(res1)<- c('peak_id',"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj" )##对列名进行修改
  write.table(res1,paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_S8495_AUC_ATAC.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
}

