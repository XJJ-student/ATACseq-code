###### Drive patients into sen vs res groups base on per drug sensitiviti   ###############
#######  根据药敏AUC数据将样本分为两类
rm(list=ls())
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
setwd("~/xjj/drug")
ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,!(colnames(ic501)%in%c("PC.100","PC.34"))]
colnames(ic50)<-sample_id[match(colnames(ic50),sample_id[,2]),1]
drug_info<-read.csv("compiled_drug_info.csv",header = TRUE,sep = ",")
drug_result<-drug_info[match(rownames(ic50),drug_info$drug_id),]

#没有做标准化
ic2<-ic50
tresult=matrix(0,length(rownames(ic2)),6)
for (i in 1:length(rownames(ic2))){
  median_value<-median(as.numeric(ic2[i,]))
  sen<-ic2[i,which(as.numeric(ic2[i,])>median_value)]
  res<-ic2[i,which(as.numeric(ic2[i,])<=median_value)]
  ttest=t.test(as.numeric(sen),as.numeric(res),alternative = "two.sided")
  sen_names<-colnames(ic2)[which(ic2[i,]>median_value)]
  res_names<-colnames(ic2)[which(ic2[i,]<=median_value)]
  res_names1<-paste0(res_names, collapse = ",")
  sen_names1<-paste0(sen_names, collapse = ",")
  tresult[i,1:5]=c(rownames(ic2)[i],ttest$statistic,ttest$p.value,res_names1,sen_names1)
}
tresult[,6]=p.adjust(tresult[,3],method="BH")
colnames(tresult)<-c("drugid","statistic","p.value","res_names","sen_names","FDR")
sum(as.numeric(tresult[,6])<0.05)
tresult_want<-tresult[as.numeric(tresult[,6])<0.05,]
setwd("~/xjj/drug/drug_result/drug64_AUC/newresult20220515")
write.table(tresult_want,"drug64_sample220515_ttest.txt",sep="\t",col.names=T,row.names=F)
write.table(tresult_want,"drug64_sample220515.txt",sep="\t",col.names=T,row.names=F)

####
wresult=matrix(0,length(rownames(ic2)),6)
for (i in 1:length(rownames(ic2))){
  median_value<-median(as.numeric(ic2[i,]))
  sen<-as.numeric(ic2[i,which(ic2[i,]>median_value)])
  res<-as.numeric(ic2[i,which(ic2[i,]<=median_value)])
  wtest=wilcox.test(as.numeric(sen),as.numeric(res),alternative = "two.sided")#wilcox.test
  sen_names<-colnames(ic2)[which(ic2[i,]>median_value)]
  res_names<-colnames(ic2)[which(ic2[i,]<=median_value)]
  res_names1<-paste0(res_names, collapse = ",")
  sen_names1<-paste0(sen_names, collapse = ",")
  wresult[i,1:5]=c(rownames(ic2)[i],wtest$statistic,wtest$p.value,res_names1,sen_names1)
}
wresult[,6]=p.adjust(wresult[,3],method="BH")
colnames(wresult)<-c("drugid","statistic","p.value","res_names","sen_names","FDR")
sum(as.numeric(wresult[,6])<0.05)
wresult_want<-wresult[as.numeric(wresult[,6])<0.05,]
setwd("~/xjj/drug/drug_result/drug64_AUC/newresult20220515")
write.table(wresult_want,"drug64_sample220515_wilcoxtest.txt",sep="\t",col.names=T,row.names=F)


############################ RNAseq DEGs#################################################
############    RNAseq DEGs  一键式运行，结果整理成一张表
rm(list=ls())
freq<-0.05
freq1<-0.01

setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
library("data.table")
exp1<-fread("star_rsem.GeneSymbol.readscount.xls",header=T,data.table=F)
exp2<-floor(exp1[-1])
delete<-apply(exp2,1,function(x) mean(x==0))
cou<-which(delete>0.75)####在超过75%样本以上都是0的基因删去
exp<-exp2[(-cou),]
setwd("~/xjj/drug/drug_result/drug64_AUC/newresult20220515")
tresult_want<-read.table("drug64_sample220515.txt",sep="\t",header=T)
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
library(DESeq2)
library(dplyr)
library(foreach)

RNAseq_result<-data.frame(matrix(0,nrow(tresult_want),9))
foreach(k=1:nrow(tresult_want),.combine=rbind)%do%{
  G=as.character(tresult_want[k,1])
  cat(G,"\n")
  file_name=as.character(tresult_want[k,1])
  up_name<-paste("P",freq,"_",G,"_up",sep="")
  down_name<-paste("P",freq,"_",G,"_down",sep="")
  up_name1<-paste("P",freq1,"_",G,"_up",sep="")
  down_name1<-paste("P",freq1,"_",G,"_down",sep="")
  
  common_sample<-as.data.frame(tresult_want[match(G,tresult_want[,1]),])
  sensitivity<-unlist(strsplit(as.character(common_sample[1,5]), ","))
  resistance<-unlist(strsplit(as.character(common_sample[1,4]), ","))
  RNAseq_result[k,1]<-G
  RNAseq_result[k,2]<-as.character(common_sample[1,5])
  RNAseq_result[k,3]<-as.character(common_sample[1,4])
  resistance_exp<-exp[,match(resistance,colnames(exp))]
  sensitivity_exp<-exp[,match(sensitivity,colnames(exp))]
  geneid<-exp1[(-cou),1]
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
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/1_RNAseq_DEGs",sep=""))
  write.table(res,paste("sensitivity-resistance-all-DESeq2_",G,"_RNAseq_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  resSig<-res[which(res$pvalue<freq),]
  #pvalue padj
  resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
  resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
  cat(paste("up_DEGs number is",sum(resSig$up_down=='up'),sep=" "),"\n")
  cat(paste("down_DEGs number is",sum(resSig$up_down=='down'),sep=" "),"\n")
  cat(paste("all_DEGs number is",nrow(resSig),sep=" "),"\n")
  RNAseq_result[k,4]<-sum(resSig$up_down=='up')
  RNAseq_result[k,5]<-sum(resSig$up_down=='down')
  RNAseq_result[k,6]<-nrow(resSig)
  up_gene<-resSig[which(resSig$log2FoldChange>0),]
  down_gene<-resSig[which(resSig$log2FoldChange<0),]
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/1_RNAseq_DEGs",sep=""))
  write.table(up_gene,paste("sensitivity-resistance-all-DESeq2_RNAseq_",up_name,"_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  write.table(down_gene,paste("sensitivity-resistance-all-DESeq2_RNAseq_",down_name,"_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  
  resSig1<-res[which(res$pvalue<freq1),]
  #pvalue padj
  resSig1[which(resSig1$log2FoldChange>0),'up_down']<-'up'
  resSig1[which(resSig1$log2FoldChange<0),'up_down']<-'down'
  cat(paste("up_DEGs number is",sum(resSig1$up_down=='up'),sep=" "),"\n")
  cat(paste("down_DEGs number is",sum(resSig1$up_down=='down'),sep=" "),"\n")
  cat(paste("all_DEGs number is",nrow(resSig1),sep=" "),"\n")
  RNAseq_result[k,7]<-sum(resSig1$up_down=='up')
  RNAseq_result[k,8]<-sum(resSig1$up_down=='down')
  RNAseq_result[k,9]<-nrow(resSig1)
  up_gene1<-resSig1[which(resSig1$log2FoldChange>0),]
  down_gene1<-resSig1[which(resSig1$log2FoldChange<0),]
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/1_RNAseq_DEGs",sep=""))
  write.table(up_gene1,paste("sensitivity-resistance-all-DESeq2_RNAseq_",up_name1,"_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  write.table(down_gene1,paste("sensitivity-resistance-all-DESeq2_RNAseq_",down_name1,"_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  
}
colnames(RNAseq_result)<-c("drug","sensitivity","resistance","up-DEGs-p0.05","down-DEGs-p0.05","all-DEGs-p0.05","up-DEGs-p0.01","down-DEGs-p0.01","all-DEGs-p0.01")
setwd("~/xjj/drug/drug_result/drug64_AUC/newresult20220515")
write.table(RNAseq_result,"RNAseq_result_64drug_FC0_P0.05_P0.01.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

############################ RNAseq DEGs Volcano Plot##########################################
#########  RNAseq DEGs Volcano Plot  一键式完成
rm(list=ls())
library(ggrepel) 
library(ggplot2)
freq<-0.05 # pvalue
FC<-0
text_freq<-0.05  #FDR
text_FC<-1
setwd("~/xjj/drug/drug_result/drug64_AUC/newresult20220515")
tresult_want<-read.table("drug64_sample220515.txt",sep="\t",header=T)

foreach(k=1:nrow(tresult_want))%do%{
  file_name=as.character(tresult_want[k,1])
  cat(file_name,"\n")
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/1_RNAseq_DEGs",sep=""))
  df<-read.table(paste("sensitivity-resistance-all-DESeq2_",file_name,"_RNAseq_DEGs.txt",sep=""),sep = '\t',header= T,row.names = 1)
  df$threshold = factor(ifelse(df$pvalue < freq & abs(df$log2FoldChange) >= FC, ifelse(df$log2FoldChange >= FC ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
  df$gene <- row.names(df) #添加一列基因名，以便备注
  p<-ggplot(df,aes(x=log2FoldChange,y= -log10(pvalue),color=threshold))+
    geom_point(data = df[df$pvalue<freq & abs(df$log2FoldChange)>FC,],size = 1)+ 
    geom_point(data = df[df$pvalue>freq | abs(df$log2FoldChange)<FC,],size = 1)+
    scale_color_manual(values=c('blue','grey','red'))+#确定点的颜色
    geom_text_repel(
      data = df[df$padj<text_freq & abs(df$log2FoldChange)>text_FC,],
      aes(label = gene),
      size = 3,
      color = "black",
      segment.color = "black", show.legend = FALSE )+#添加关注的点的基因名
    ylab('-log10 (pvalue)')+#修改y轴名称
    xlab('log2 (FoldChange)')+#修改x轴名称
    ggtitle(file_name)+
    theme(plot.title = element_text(hjust = 0.5))+
    geom_hline(yintercept=-log10(freq),linetype=(text_FC*2))#添加横线|logFoldChange|>0.25
  geom_vline(xintercept=c(-(text_FC),text_FC),linetype=4)#添加竖线padj<0.05
  #plot(p)
  DEGs<-df[which(df$padj<text_freq & abs(df$log2FoldChange)>text_FC),]
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/1_RNAseq_DEGs/Volcano",sep=""))
  ggsave(paste("FDR_",text_freq,"_FC_",text_FC,"_DEGs.pdf",sep=""),p,width = 10, height = 6)
  write.table(DEGs,paste("FDR_",text_freq,"_FC_",text_FC,"_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
}

############################ DEGs与组蛋白修饰酶的交叠 一键式完成  ########################################
################  DEGs与组蛋白修饰酶的交叠 组蛋白修饰为表观遗传修饰
rm(list=ls())
library(gplots)
library(VennDiagram)
library(foreach)
freq<-0.01 #差异基因的阈值
FC<-0
freq_name<-"P0.05"
setwd("~/xjj/drug/drug_result/drug64_AUC/newresult20220515")
tresult_want<-read.table("drug64_sample220515.txt",sep="\t",header=T)
setwd("~/xjj/drug/drug_result/HDAC_frontiers/5_histone-modifying enzymes")
Acetylase<-read.table("histone-modifying-enzymes.txt",header=T,sep="\t")


Acetylase_result<-data.frame(matrix(0,nrow(tresult_want),2))
foreach(k=1:nrow(tresult_want))%do%{
  file_name=as.character(tresult_want[k,1])
  Acetylase_result[k,1]<-file_name
  cat(file_name,"\n")
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/1_RNAseq_DEGs",sep=""))
  df<-read.table(paste("sensitivity-resistance-all-DESeq2_",file_name,"_RNAseq_DEGs.txt",sep=""),sep = '\t',header= T,row.names = 1)
  DEGs<-df[df$pvalue<freq,]
  DEGs[which(DEGs$log2FoldChange>FC),'up_down']<-'up'
  DEGs[which(DEGs$log2FoldChange<FC),'up_down']<-'down'
  #pvalue padj
  Acetylase_DEGs<-intersect(toupper(Acetylase[,1]),toupper(rownames(DEGs)))
  Acetylase_info1<-Acetylase[match(Acetylase_DEGs,toupper(Acetylase[,1])),]
  DEGs$gene<-toupper(rownames(DEGs))
  Acetylase_info2<-DEGs[match(Acetylase_DEGs,toupper(rownames(DEGs))),7:8]
  Acetylase_info<-cbind(Acetylase_info1,Acetylase_info2)
  Acetylase_name<- unique(toupper(Acetylase[,1]))
  DEGs_name <- unique(toupper(rownames(DEGs)))
  input  <-list(Acetylase_name,DEGs_name)
  ven<-venn(input,showSetLogicLabel=TRUE)
  plot(ven)
  tmp <- venn(input)
  int<-attr(tmp, "intersections")
  Acetylase_result[k,2]<-length(int[["A:B"]])
  Acetylase_result[k,3]<-paste0(int[["A:B"]],collapse =",")
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/1_RNAseq_DEGs/Acetylase",sep=""))
  write.table(Acetylase_info,paste("Acetylase_info_",freq_name,".txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  venn.diagram(x=list(Acetylase=Acetylase_name,DEGs=DEGs_name),cex = 1,margin = 0.3,
               paste("Acetylase_info_",freq_name,".png",sep=""),fill=c("red","blue"))
}

colnames(Acetylase_result)<-c("drug","Number of DEAcetylase","Acetylase")
setwd("~/xjj/drug/drug_result/drug64_AUC/newresult20220515")
write.table(Acetylase_result,"DEAcetylase_result_64drug_P0.01.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

############################ RNAseq clusterProfiler#######################################################
################  clusterProfiler
rm(list=ls())
library(org.Hs.eg.db)
library(clusterProfiler)
library(stringr)
library(foreach)
library(ggplot2)
kegg_freq_name<-"FDR0.05"
go_freq_name<-"FDR0.01"
keggfreq<-0.05 #kegg p<0.05
gofreq<-0.01    #go FDR<0.05
setwd("~/xjj/drug/drug_result/drug64_AUC/newresult20220515")
tresult_want<-read.table("drug64_sample220515.txt",sep="\t",header=T)

clusterProfiler_result<-data.frame(matrix(0,nrow(tresult_want),5))
foreach(k=1:nrow(tresult_want))%do%{
  file_name=as.character(tresult_want[k,1])
  clusterProfiler_result[k,1]<-file_name
  cat(file_name,"\n")
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/1_RNAseq_DEGs",sep=""))
  up_gene<-read.table(paste("sensitivity-resistance-all-DESeq2_RNAseq_P0.05_",file_name,"_up_DEGs.txt",sep=""),sep = '\t',header= T)
  geneu=bitr(up_gene[,1],fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
  keggu<-enrichKEGG(geneu[,2],organism = 'hsa',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'BH',
                    minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.05,use_internal_data = FALSE)
  keu<-barplot(keggu,showCategory=10,drop=T,x = "GeneRatio",color = "p.adjust")
  barplot(keggu,showCategory=10,drop=T,x = "GeneRatio",color = "p.adjust")
  
  eeu<-keggu@result
  ee1u<-eeu[which(eeu$p.adjust<keggfreq),]
  cat(paste("The pathway of KEGG is",nrow(ee1u),sep=" "),"\n")
  clusterProfiler_result[k,2]<-nrow(ee1u)
  #p.adjust
  enrich_geneu<-ee1u$geneID
  pathway_geneu<-unique(unlist(strsplit(enrich_geneu,split="/")))
  pathway_gene2u=bitr(pathway_geneu,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/1_RNAseq_DEGs/pathway",sep=""))
  ggsave(paste("kegg_DEG_pathway_",kegg_freq_name,"_up.pdf",sep=""),keu,width = 8, height = 5)
  write.table(ee1u,paste("kegg_DEG_pathway_",kegg_freq_name,"_up.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  write.table(pathway_gene2u,paste("kegg_DEG_pathway_",kegg_freq_name,"_up_gene.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/1_RNAseq_DEGs",sep=""))
  down_gene<-read.table(paste("sensitivity-resistance-all-DESeq2_RNAseq_P0.05_",file_name,"_down_DEGs.txt",sep=""),sep = '\t',header= T)
  gened=bitr(down_gene[,1],fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
  keggd<-enrichKEGG(gened[,2],organism = 'hsa',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'BH',
                    minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.05,use_internal_data = FALSE)
  ked<-barplot(keggd,showCategory=10,drop=T,x = "GeneRatio",color = "p.adjust")
  barplot(keggd,showCategory=10,drop=T,x = "GeneRatio",color = "p.adjust")
  
  eed<-keggd@result
  ee1d<-eed[which(eed$p.adjust<keggfreq),]
  cat(paste("The pathway of KEGG is",nrow(ee1d),sep=" "),"\n")
  clusterProfiler_result[k,3]<-nrow(ee1d)
  #p.adjust
  enrich_gened<-ee1d$geneID
  pathway_gened<-unique(unlist(strsplit(enrich_gened,split="/")))
  pathway_gene2d=bitr(pathway_gened,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/1_RNAseq_DEGs/pathway",sep=""))
  ggsave(paste("kegg_DEG_pathway_",kegg_freq_name,"_down.pdf",sep=""),ked,width = 8, height = 10)
  write.table(ee1d,paste("kegg_DEG_pathway_",kegg_freq_name,"_down.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  write.table(pathway_gene2d,paste("kegg_DEG_pathway_",kegg_freq_name,"_down_gene.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  
  gou<-enrichGO(geneu[,2],OrgDb = org.Hs.eg.db,ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.01,keyType = 'ENTREZID')
  au<-gou@result
  go_BPu<-au[which(au$p.adjust<gofreq),]
  cat(paste("The pathway of GO is",nrow(go_BPu),sep=" "),"\n")
  clusterProfiler_result[k,4]<-nrow(go_BPu)
  gouu<-dotplot(gou,showCategory=10)
  ggsave(paste("GO_DEG_pathway_",go_freq_name,"_up.pdf",sep=""),gouu,width = 8, height = 5)
  enrich_genegou<-go_BPu$geneID
  pathway_genegou<-unique(unlist(strsplit(enrich_genegou,split="/")))
  pathway_genego2u=bitr(pathway_genegou,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
  write.table(pathway_genego2u,paste("GO_DEG_pathway_",go_freq_name,"_up_gene.txt",sep=""),sep = '\t',col.names = F,row.names = F,quote = FALSE)##数据输出
  write.table(go_BPu,paste("GO_DEG_pathway_",go_freq_name,"_up.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  god<-enrichGO(gened[,2],OrgDb = org.Hs.eg.db,ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.01,keyType = 'ENTREZID')
  godd<-dotplot(god,showCategory=10)
  ggsave(paste("GO_DEG_pathway_",go_freq_name,"_down.pdf",sep=""),godd,width = 8, height = 5)
  ad<-god@result
  go_BPd<-ad[which(ad$p.adjust<gofreq),]
  cat(paste("The pathway of GO is",nrow(go_BPd),sep=" "),"\n")
  clusterProfiler_result[k,5]<-nrow(go_BPd)
  enrich_genegod<-go_BPd$geneID
  pathway_genegod<-unique(unlist(strsplit(enrich_genegod,split="/")))
  pathway_genego2d=bitr(pathway_genegod,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
  write.table(pathway_genego2d,paste("GO_DEG_pathway_",go_freq_name,"_down_gene.txt",sep=""),sep = '\t',col.names = F,row.names = F,quote = FALSE)##数据输出
  write.table(go_BPd,paste("GO_DEG_pathway_",go_freq_name,"_down.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
}
colnames(clusterProfiler_result)<-c("drug","KEGG_up_FDR0.05","KEGG_down_FDR0.05","GO_up_FDR0.01","GO_down_FDR0.01")
setwd("~/xjj/drug/drug_result/drug64_AUC/newresult20220515")
write.table(clusterProfiler_result,"RNAseq_clusterProfiler_result_64drug_DEGP0.05_KEGG_FDR0.05_GO_FDR0.01.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
############################ TCGA_PAAD#######################################
#############  TCGA_PAAD   ####### 有符合条件的基因
rm(list=ls())
library(org.Hs.eg.db)
library(stringr)
library(clusterProfiler)
library(survival)
library(survminer)

setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/6_TCGA") #老位置不变
chemotherapy_drug_sample<-read.table("drug_paad1.txt",header=T,sep="\t") #就用这个文件，化疗药物信息
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
setwd("~/xjj/drug/drug_result/drug64_AUC/newresult20220515")
tresult_want<-read.table("drug64_sample220515.txt",sep="\t",header=T)

for(k in 1:nrow(tresult_want)){
  file_name=tresult_want[k,1]
  
}
setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/1_RNAseq_DEGs/Volcano",sep=""))
ggsave(paste("FDR_",text_freq,"_FC_",text_FC,"_DEGs.pdf",sep=""),p,width = 10, height = 6)
write.table(DEGs,paste("FDR_",text_freq,"_FC_",text_FC,"_DEGs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出


setwd("~/xjj/drug/drug_result/chemotherapy/1_RNAseq")
df<-read.table("sensitivity-resistance-all-DESeq2_chemotherapy_AUC_RNAseq_DEGs.txt",sep = '\t',header= T)
#g2<-as.character(df[which(df$pvalue<0.05),1])
g2 = as.character(df[which(df$pvalue<0.05),1])
g2 =as.character(df[which(df$padj<0.05 & abs(df$log2FoldChange)>2),1])

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

setwd("~/xjj/drug/drug_result/chemotherapy/1_RNAseq/TCGA_PAAD")
library(ggplot2)
##手动保存
ggsave(paste(g,"up-regulation in sensitive.pdf",sep=" "),p,width = 5, height = 5)
write.table(a,"DEGs_P0.05_logrank_P0.05_TCGA_survived.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(result,"FDR05_FC2_DEGs_TCGA_survived.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出


############################ ATACseq 一体式   ###############################
#############  DApeaks
rm(list=ls())
freq<-0.05
freq1<-0.01
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_count_name<-read.table("peak_count.txt",sep="\t",header=T)
setwd("~/xjj/drug/drug_result/drug64_AUC/newresult20220515")
tresult_want<-read.table("drug64_sample220515.txt",sep="\t",header=T)

setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
library(DESeq2)
library(dplyr)


ATACseq_result<-data.frame(matrix(0,30,9))
for(k in 35:64){
  G=as.character(tresult_want[k,1])
  cat(G,"\n")
  file_name=as.character(tresult_want[k,1])
  out_file<-paste("sen-res-",G,"-0.05-P",sep="")
  out_file1<-paste("sen-res-",G,"-0.01-P",sep="")
  
  common_sample<-as.data.frame(tresult_want[match(G,tresult_want[,1]),])
  sensitivity<-unlist(strsplit(as.character(common_sample[1,5]), ","))
  resistance<-unlist(strsplit(as.character(common_sample[1,4]), ","))
  ATACseq_result[k,1]<-G
  ATACseq_result[k,2]<-as.character(common_sample[1,5])
  ATACseq_result[k,3]<-as.character(common_sample[1,4])
  colnames(peak_count_name)<-gsub("\\.","-",colnames(peak_count_name))
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
  #sizeFactors(dds)##查看每个主成分的标准化值
  res<-results(dds)##将结果输出
  res<-as.data.frame(res)
  res<-cbind(peak_count_name[,1],res)  ##对数据增加一列
  colnames(res)<- c('peak_id',"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj" )##对列名进行修改
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/2_ATACseq_DApeaks",sep=""))
  write.table(res,paste("sensitivity-resistance-all-DESeq2_",file_name,"_ATAC.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  resSig<-res[which(res$pvalue<freq),]
  # pvalue padj
  resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
  resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
  cat(paste("up DApeaks is",sum(resSig$up_down=='up'),sep=" "),"\n")
  cat(paste("down DApeaks is",sum(resSig$up_down=='down'),sep=" "),"\n")
  cat(paste("all DApeaks is",nrow(resSig),sep=" "),"\n")
  ATACseq_result[k,4]<-sum(resSig$up_down=='up')
  ATACseq_result[k,5]<-sum(resSig$up_down=='down')
  ATACseq_result[k,6]<-nrow(resSig)
  up_gene<-as.data.frame(resSig[which(resSig$log2FoldChange>0),1])
  down_gene<-as.data.frame(resSig[which(resSig$log2FoldChange<0),1])
  up1<-match(up_gene[,1],peak_count_name[,1])
  down1<-match(down_gene[,1],peak_count_name[,1])
  a_up<-apply(up_gene,1,function(x) unlist(strsplit(as.character(x), "_")))
  d_up<-t(a_up)
  e_up<-data.frame(rep("+",nrow(d_up)))
  peaks_up<-cbind(up_gene,d_up,e_up)
  colnames(peaks_up)<-c("peak_id","chr","start","end","strand")
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/2_ATACseq_DApeaks",sep=""))
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
#################  
  resSig1<-res[which(res$pvalue<freq1),]
  # pvalue padj
  resSig1[which(resSig1$log2FoldChange>0),'up_down']<-'up'
  resSig1[which(resSig1$log2FoldChange<0),'up_down']<-'down'
  cat(paste("up DApeaks is",sum(resSig1$up_down=='up'),sep=" "),"\n")
  cat(paste("down DApeaks is",sum(resSig1$up_down=='down'),sep=" "),"\n")
  cat(paste("all DApeaks is",nrow(resSig1),sep=" "),"\n")
  ATACseq_result[k,7]<-sum(resSig1$up_down=='up')
  ATACseq_result[k,8]<-sum(resSig1$up_down=='down')
  ATACseq_result[k,9]<-nrow(resSig1)
  up_gene1<-as.data.frame(resSig1[which(resSig1$log2FoldChange>0),1])
  down_gene1<-as.data.frame(resSig1[which(resSig1$log2FoldChange<0),1])
  up11<-match(up_gene1[,1],peak_count_name[,1])
  down11<-match(down_gene1[,1],peak_count_name[,1])
  a_up1<-apply(up_gene1,1,function(x) unlist(strsplit(as.character(x), "_")))
  d_up1<-t(a_up1)
  e_up1<-data.frame(rep("+",nrow(d_up1)))
  peaks_up1<-cbind(up_gene1,d_up1,e_up1)
  colnames(peaks_up1)<-c("peak_id","chr","start","end","strand")
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/2_ATACseq_DApeaks",sep=""))
  write.table(peaks_up1,file=paste(out_file1,"-logFC0-DESeq2-up.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  write.table(d_up1,file=paste(out_file1,"-logFC0-DESeq2-up.gff",sep=""),sep = '\t',
              col.names = F,row.names = F,quote = FALSE)
  a_down1<-apply(down_gene1,1,function(x) unlist(strsplit(as.character(x), "_")))
  d_down1<-t(a_down1)
  e_down1<-data.frame(rep("+",nrow(d_down1)))
  peaks_down1<-cbind(down_gene1,d_down1,e_down1)
  colnames(peaks_down1)<-c("peak_id","chr","start","end","strand")
  write.table(peaks_down1,file=paste(out_file1,"-logFC0-DESeq2-down.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  write.table(d_down1,file=paste(out_file1,"-logFC0-DESeq2-down.gff",sep=""),sep = '\t',
              col.names = F,row.names = F,quote = FALSE)
  #####  bed
  f_up1<-data.frame(rep(1,nrow(d_up1)))
  peaks_bed_up1<-cbind(d_up1,up_gene1,f_up1,e_up1)
  colnames(peaks_bed_up1)<-c("chr","start","end","peak_id","score","strand")
  write.table(peaks_bed_up1,file=paste(out_file1,"-logFC0-DESeq2-up.bed",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  f_down1<-data.frame(rep(1,nrow(d_down1)))
  peaks_bed_down1<-cbind(d_down1,down_gene1,f_down1,e_down1)
  colnames(peaks_bed_down1)<-c("chr","start","end","peak_id","score","strand")
  write.table(peaks_bed_down1,file=paste(out_file1,"-logFC0-DESeq2-down.bed",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  
}
colnames(ATACseq_result)<-c("drug","sensitivity","resistance","up-peaks-p0.05","down-peaks-p0.05","all-peaks-p0.05","up-peaks-p0.01","down-peaks-p0.01","all-peaks-p0.01")
setwd("~/xjj/drug/drug_result/drug64_AUC/newresult20220515")
write.table(ATACseq_result,"ATACseq_result_64drug_FC0_P0.05_P0.01.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出












############################ peaks associated with drug sensitivity##############################
######   a chromatin accessibility signature associated with drug sensitivity
######   including both positive and negative correlations
rm(list=ls())
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
setwd("~/xjj/drug")
ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,!(colnames(ic501)%in%c("PC.100","PC.34"))]## all of PDAC
colnames(ic50)<-sample_id[match(colnames(ic50),sample_id[,2]),1]

setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_RPKM_name<-read.table("peak_RPKM.txt",sep="\t",header=T)
colnames(peak_RPKM_name)<-gsub("\\.","-",colnames(peak_RPKM_name))
peak_clinical<-peak_RPKM_name[,match(colnames(ic50),colnames(peak_RPKM_name))]
peak_name<-as.character(peak_RPKM_name[,1])

correlation_all<-matrix(0,64,7)
for(i in 1:nrow(ic50)){
  file_name<-rownames(ic50)[i]
  cat(file_name,"\n")
  correlation_all[i,1]<-file_name
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/2_ATACseq_DApeaks",sep=""))
  DApeak<-fread(paste("sensitivity-resistance-all-DESeq2_",file_name,"_ATAC.txt",sep=""),header=T,data.table=F)
  DAR<-as.character(DApeak[which(DApeak$pvalue<0.01),1])
  DApeaks<-peak_name[match(DAR,peak_name)] #DApeaks and position
  DApeak_clinical<-peak_clinical[match(DAR,peak_name),]
  chemotherapy5_29<-ic50[match(file_name,rownames(ic50)),]
  corralation_all1<-matrix(0,nrow(DApeak_clinical),8)
  for(k in 1:nrow(DApeak_clinical)){
    corralation_all1[k,1]<-file_name
    corralation_all1[k,2]<-as.character(DApeaks[k])
    corralation_all1[k,3]<-"spearman"
    cor_result_spearman<-cor.test(as.numeric(chemotherapy5_29),as.numeric(DApeak_clinical[k,]),alternative = "two.sided",method = "spearman")
    corralation_all1[k,4]<-cor_result_spearman$estimate
    corralation_all1[k,5]<-cor_result_spearman$p.value
    corralation_all1[k,6]<-"pearson"
    cor_result_pearson<-cor.test(as.numeric(chemotherapy5_29),as.numeric(DApeak_clinical[k,]),alternative = "two.sided",method = "pearson")
    corralation_all1[k,7]<-cor_result_pearson$estimate
    corralation_all1[k,8]<-cor_result_pearson$p.value
  }
  corralation_all2<-corralation_all1[which(as.numeric(corralation_all1[,5])<0.05),]
  corralation_all3<-corralation_all1[which(as.numeric(corralation_all1[,8])<0.05),]
  colnames(corralation_all2)<-c("Drug","Peak_name","spearman","spearman_estimate","spearman_pvalue","pearson","pearson_estimate","pearson_pvalue")
  colnames(corralation_all3)<-c("Drug","Peak_name","spearman","spearman_estimate","spearman_pvalue","pearson","pearson_estimate","pearson_pvalue")
  correlation_all[i,2]<-nrow(corralation_all2)
  correlation_all[i,3]<-sum(as.numeric(corralation_all2[,4])>0)
  correlation_all[i,4]<-sum(as.numeric(corralation_all2[,4])<0)
  correlation_all[i,5]<-nrow(corralation_all3)
  correlation_all[i,6]<-sum(as.numeric(corralation_all3[,7])>0)
  correlation_all[i,7]<-sum(as.numeric(corralation_all3[,7])<0)
  outputPath = paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/8_ATACseq_drugs",sep='')
  if (!file.exists(outputPath)){
    dir.create(outputPath, recursive=TRUE)
  }
  write.table(corralation_all2,paste(outputPath,'/',file_name,"-spearman-P0.05-DARs-drugs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)
  write.table(corralation_all3,paste(outputPath,'/',file_name,"-pearson-P0.05-DARs-drugs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)
}
colnames(correlation_all)<-c("drug","Spearman-significant-correlation-all","Spearman-positive","Spearman-negative","Pearson-significant-correlation-all","Pearson-positive","Pearson-negative")
setwd("~/xjj/drug/drug_result/drug64_AUC/newresult20220515")
write.table(correlation_all,"ATACseq_64drug_spearman_pearson.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

############################ peaks identify subgroup of patients###########################################
library(data.table)
file_name<-"GEM"
setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/8_ATACseq_drugs",sep=""))
spearman_DARs_drugs<-fread(paste(file_name,"-spearman-P0.05-DARs-drugs.txt",sep=""),header=T,data.table=F)
one_positive<-spearman_DARs_drugs[which(spearman_DARs_drugs$spearman_estimate>0),2]

setwd("~/xjj/drug")
drug_info<-read.csv("compiled_drug_info.csv",header = TRUE,sep = ",")
HDAC1<-drug_info[drug_info$target=="HDAC",]
ccc<-c(as.character(HDAC1[,1]),"S1848","S2759","S1194","S1047","5-FU","GEM","IRI","OXA","PAC")

result<-matrix(0,length(one_positive),(length(ccc)+1))
for(i in 1:length(ccc)){
  file_name<-ccc[i]
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/8_ATACseq_drugs",sep=""))
  spearman_DARs_drugs2<-fread(paste(file_name,"-spearman-P0.05-DARs-drugs.txt",sep=""),header=T,data.table=F)
  two_negative<-spearman_DARs_drugs2[which(spearman_DARs_drugs2$spearman_estimate<0),2]
  int_peaks<-intersect(one_positive,two_negative)
  result[match(int_peaks,one_positive),i+1]<-1
}
a<-apply(result[,-1],1,sum)
result[,1]<-one_positive
result[order(a,decreasing=T)[1:21],1]

############################ DApeaks-DEGs-drug-correlation############################################
#######  DApeaks-DEGs-drug-correlation
rm(list=ls())
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
setwd("~/xjj/drug")
ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,!(colnames(ic501)%in%c("PC.100","PC.34"))]## all of PDAC
colnames(ic50)<-sample_id[match(colnames(ic50),sample_id[,2]),1]
setwd("~/xjj/drug")
drug_info<-read.csv("compiled_drug_info.csv",header = TRUE,sep = ",")
HDAC1<-drug_info[drug_info$target=="HDAC",]
ccc<-c(as.character(HDAC1[,1]),"S1848","S2759","S1194","S1047","5-FU","GEM","IRI","OXA","PAC")
chemotherapy5<-ic50[match(ccc,rownames(ic50)),]

setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_RPKM_name<-read.table("peak_RPKM.txt",sep="\t",header=T)
colnames(peak_RPKM_name)<-gsub("\\.","-",colnames(peak_RPKM_name))
peak_clinical<-peak_RPKM_name[,match(colnames(ic50),colnames(peak_RPKM_name))]
peak_name<-peak_RPKM_name[,1]

correlation_all<-matrix(0,25,7)
for(i in 1:25){
  file_name<-rownames(ic50)[i]
  cat(file_name,"\n")
  correlation_all[i,1]<-file_name
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/2_ATACseq_DApeaks",sep=""))
  library(data.table)
  homer_anno_down<-fread(paste("sen-res-",file_name,"-0.01-P-logFC0-DESeq2-down-annotation.txt",sep=""),header=T,data.table=F)
  colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
  Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]
  homer_anno_up<-fread(paste("sen-res-",file_name,"-0.01-P-logFC0-DESeq2-up-annotation.txt",sep=""),header=T,data.table=F)
  colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
  Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/1_RNAseq_DEGs",sep=""))
  up_gene<-read.table(paste("sensitivity-resistance-all-DESeq2_RNAseq_P0.05_",file_name,"_up_DEGs.txt",sep=""),header=T,sep="\t")
  down_gene<-read.table(paste("sensitivity-resistance-all-DESeq2_RNAseq_P0.05_",file_name,"_down_DEGs.txt",sep=""),header=T,sep="\t")
  downDEGs_peaks<-c()
  for(j in 1:nrow(down_gene)){
    DEGs_peaks1<-Anno_gene_100Kb_down[Anno_gene_100Kb_down$`Gene Name` %in% down_gene[j,1],1]
    downDEGs_peaks<-c(downDEGs_peaks,DEGs_peaks1)
  }
  upDEGs_peaks<-c()
  for(j in 1:nrow(up_gene)){
    DEGs_peaks2<-Anno_gene_100Kb_up[Anno_gene_100Kb_up$`Gene Name` %in% up_gene[j,1],1]
    upDEGs_peaks<-c(upDEGs_peaks,DEGs_peaks2)
  }
  DAR<-unique(c(downDEGs_peaks,upDEGs_peaks))
  DApeaks<-peak_name[match(DAR,peak_name)] #DApeaks and position
  DApeak_clinical<-peak_clinical[match(DAR,peak_name),]

  chemotherapy5_29<-ic50[match(file_name,rownames(ic50)),]
  corralation_all1<-matrix(0,nrow(DApeak_clinical),9)
for(k in 1:nrow(DApeak_clinical)){
  corralation_all1[k,1]<-file_name
  corralation_all1[k,2]<-as.character(DApeaks[k])
  corralation_all1[k,3]<-"spearman"
  cor_result_spearman<-cor.test(as.numeric(chemotherapy5_29),as.numeric(DApeak_clinical[k,]),alternative = "two.sided",method = "spearman")
  corralation_all1[k,4]<-cor_result_spearman$estimate
  corralation_all1[k,5]<-cor_result_spearman$p.value
  corralation_all1[k,6]<-"pearson"
  cor_result_pearson<-cor.test(as.numeric(chemotherapy5_29),as.numeric(DApeak_clinical[k,]),alternative = "two.sided",method = "pearson")
  corralation_all1[k,7]<-cor_result_pearson$estimate
  corralation_all1[k,8]<-cor_result_pearson$p.value
  Anno_gene<-rbind(Anno_gene_100Kb_down,Anno_gene_100Kb_up)
  corralation_all1[k,9]<-Anno_gene$`Gene Name`[match(as.character(DApeaks[k]),Anno_gene[,1])]
}
corralation_all2<-corralation_all1[which(as.numeric(corralation_all1[,5])<0.05),]
corralation_all3<-corralation_all1[which(as.numeric(corralation_all1[,8])<0.05),]
colnames(corralation_all2)<-c("Drug","Peak_name","spearman","spearman_estimate","spearman_pvalue","pearson","pearson_estimate","pearson_pvalue","gene")
colnames(corralation_all3)<-c("Drug","Peak_name","spearman","spearman_estimate","spearman_pvalue","pearson","pearson_estimate","pearson_pvalue","gene")
correlation_all[i,2]<-nrow(corralation_all2)
correlation_all[i,3]<-sum(as.numeric(corralation_all2[,4])>0)
correlation_all[i,4]<-sum(as.numeric(corralation_all2[,4])<0)
correlation_all[i,5]<-nrow(corralation_all3)
correlation_all[i,6]<-sum(as.numeric(corralation_all3[,7])>0)
correlation_all[i,7]<-sum(as.numeric(corralation_all3[,7])<0)
outputPath = paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/8_ATACseq_drugs",sep='')
if (!file.exists(outputPath)){
  dir.create(outputPath, recursive=TRUE)
}
write.table(corralation_all2,paste(outputPath,'/',file_name,"-spearman-P0.05-DARs-DEGs-drugs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)
write.table(corralation_all3,paste(outputPath,'/',file_name,"-pearson-P0.05-DARs-DEGs-drugs.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)
}
colnames(correlation_all)<-c("drug","Spearman-significant-correlation-all","Spearman-positive","Spearman-negative","Pearson-significant-correlation-all","Pearson-positive","Pearson-negative")
setwd("~/xjj/drug/drug_result/drug64_AUC/newresult20220515")
write.table(correlation_all,"ATACseq_64drug_spearman_pearson_DARs_DEGs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出



















############################ 筛选每种药物特异最显著的peaks  ################################





############################ 在化疗药中不敏感的peaks,在HDAC中敏感  try##########################
outputPath = paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/8_ATACseq_drugs",sep='')
spearman_DARs_drugs<-read.table(paste(outputPath,'/',file_name,"-spearman-P0.05-DARs-drugs.txt",sep=""),sep = '\t',header=T)

ccc<-c("5-FU","GEM","IRI","OXA","PAC")
i=1
file_name=ccc[i]
setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/8_ATACseq_drugs",sep=''))
up<-read.table(paste(file_name,"-spearman-P0.05-DARs-drugs.txt",sep=""),sep = '\t',header=T)
cat(dim(up),"\n")
all_chemo_up<-as.character(up[which(as.numeric(up$spearman_pvalue)<0.05 & as.numeric(up$spearman_estimate)>0),2])



ccc<-c("5-FU","GEM","IRI","OXA","PAC")
i=1
file_name=ccc[i]
setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/8_ATACseq_drugs",sep=''))
up<-read.table(paste(file_name,"-spearman-P0.05-DARs-drugs.txt",sep=""),sep = '\t',header=T)
cat(dim(up),"\n")
all_chemo_up<-as.character(up[which(as.numeric(up$spearman_pvalue)<0.05 & as.numeric(up$spearman_estimate)>0),2])
for(i in 2:length(ccc)){
  file_name=ccc[i]
  cat(file_name,"\n")
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/8_ATACseq_drugs",sep=''))
  up<-read.table(paste(file_name,"-spearman-P0.05-DARs-drugs.txt",sep=""),sep = '\t',header=T)
  cat(dim(up),"\n")
  all_chemo_up1<-as.character(up[which(as.numeric(up$spearman_pvalue)<0.05 & as.numeric(up$spearman_estimate)>0),2])
  all_chemo_up<-intersect(all_chemo_up,all_chemo_up1)
}

setwd("~/xjj/drug")
drug_info<-read.csv("compiled_drug_info.csv",header = TRUE,sep = ",")
HDAC1<-drug_info[drug_info$target=="HDAC",]
ccc<-c(as.character(HDAC1[,1]),"S1848","S2759","S1194","S1047")
i=1
file_name=ccc[i]
cat(file_name,"\n")
setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/2_ATACseq_DApeaks",sep=""))
down<-read.table(paste("sen-res-",file_name,"-0.01-P-logFC0-DESeq2-down.txt",sep=""),sep = '\t',header=T)
cat(dim(down),"\n")
all_HDAC_down<-as.character(down[,1])
for(i in 2:length(ccc)){
  file_name=ccc[i]
  cat(file_name,"\n")
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/2_ATACseq_DApeaks",sep=""))
  down<-read.table(paste("sen-res-",file_name,"-0.01-P-logFC0-DESeq2-down.txt",sep=""),sep = '\t',header=T)
  cat(dim(down),"\n")
  all_HDAC_down<-intersect(all_HDAC_down,as.character(down[,1]))
}

##对化疗敏感，所以在化疗中是up，在HDAC中耐药，所以在HDAC中down
int_peak<-intersect(all_chemo_up,all_HDAC_down)
chemo_only<-all_chemo_up[-match(int_peak,all_chemo_up)]
HDAC_only<-all_HDAC_down[-match(int_peak,all_HDAC_down)]


setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
setwd("~/xjj/drug")
ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,!(colnames(ic501)%in%c("PC.100","PC.34"))]## all of PDAC
colnames(ic50)<-sample_id[match(colnames(ic50),sample_id[,2]),1]
ccc<-c("5-FU","GEM","IRI","OXA","PAC")
chemotherapy5<-ic50[match(ccc,rownames(ic50)),]

setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_RPKM_name<-read.table("peak_RPKM.txt",sep="\t",header=T)
colnames(peak_RPKM_name)<-gsub("\\.","-",colnames(peak_RPKM_name))
peak_name<-peak_RPKM_name[,1]
peak_clinical<-peak_RPKM_name[,match(colnames(ic50),colnames(peak_RPKM_name))]

DAR<-int_peak
DApeaks<-peak_name[match(DAR,peak_name)] #DApeaks and position
DApeak_clinical<-peak_clinical[match(DAR,peak_name),]

corralation_all_s1<-NULL
corralation_all_p1<-NULL
for(i in 1:nrow(chemotherapy5)){
  drug_name<-rownames(chemotherapy5)[i]
  cat(drug_name,"\n")
  corralation_all1<-matrix(0,nrow(DApeak_clinical),8)
  for(k in 1:nrow(DApeak_clinical)){
    corralation_all1[k,1]<-drug_name
    corralation_all1[k,2]<-as.character(DApeaks[k])
    corralation_all1[k,3]<-"spearman"
    cor_result_spearman<-cor.test(as.numeric(chemotherapy5[i,]),as.numeric(DApeak_clinical[k,]),alternative = "two.sided",method = "spearman")
    corralation_all1[k,4]<-cor_result_spearman$estimate
    corralation_all1[k,5]<-cor_result_spearman$p.value
    corralation_all1[k,6]<-"pearson"
    cor_result_pearson<-cor.test(as.numeric(chemotherapy5[i,]),as.numeric(DApeak_clinical[k,]),alternative = "two.sided",method = "pearson")
    corralation_all1[k,7]<-cor_result_pearson$estimate
    corralation_all1[k,8]<-cor_result_pearson$p.value
  }
  corralation_all2<-corralation_all1[which(as.numeric(corralation_all1[,5])<0.05),]
  corralation_all3<-corralation_all1[which(as.numeric(corralation_all1[,8])<0.05),]
  corralation_all_s1<-rbind(corralation_all_s1,corralation_all2)
  corralation_all_p1<-rbind(corralation_all_p1,corralation_all3)
}
colnames(corralation_all_s)<-c("Drug","Peak_name","spearman","spearman_estimate","spearman_pvalue","pearson","pearson_estimate","pearson_pvalue")
colnames(corralation_all_p)<-c("Drug","Peak_name","spearman","spearman_estimate","spearman_pvalue","pearson","pearson_estimate","pearson_pvalue")



for(i in 1:nrow(chemotherapy5_29)){
  drug_name<-rownames(chemotherapy5_29)[i]
  cat(drug_name,"\n")
  drug_profile<-corralation_all_s[corralation_all_s[,1] %in% drug_name,]
  Top5_cor<-drug_profile[order(drug_profile[,4],decreasing=TRUE)[1:5],]
}


####   单药对单药
ccc<-c("5-FU","GEM","IRI","OXA","PAC")
i=1
file_name=ccc[i]
setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/2_ATACseq_DApeaks",sep=""))
up<-read.table(paste("sen-res-",file_name,"-0.01-P-logFC0-DESeq2-up.txt",sep=""),sep = '\t',header=T)
cat(dim(up),"\n")
all_chemo_up<-as.character(up[,1])

setwd("~/xjj/drug")
drug_info<-read.csv("compiled_drug_info.csv",header = TRUE,sep = ",")
HDAC1<-drug_info[drug_info$target=="HDAC",]
ccc<-c(as.character(HDAC1[,1]),"S1848","S2759","S1194","S1047")
i=1
file_name=ccc[i]
cat(file_name,"\n")
setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/2_ATACseq_DApeaks",sep=""))
down<-read.table(paste("sen-res-",file_name,"-0.01-P-logFC0-DESeq2-down.txt",sep=""),sep = '\t',header=T)
cat(dim(down),"\n")
all_HDAC_down<-as.character(down[,1])

##对化疗敏感，所以在化疗中是up，在HDAC中耐药，所以在HDAC中down
int_peak<-intersect(all_chemo_up,all_HDAC_down)
chemo_only<-all_chemo_up[-match(int_peak,all_chemo_up)]
HDAC_only<-all_HDAC_down[-match(int_peak,all_HDAC_down)]




######################################################################
rm(list=ls())
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
ccc<-c(as.character(HDAC1[,1]),"S1848","S2759","S1194","S1047")
bbb<-c("5-FU","GEM","IRI","OXA","PAC")

result<-NULL
for(j in 1:5){
  G<-bbb[j]
  common_sample<-as.data.frame(tresult_want[match(G,tresult_want[,1]),])
  sensitivity<-unlist(strsplit(as.character(common_sample[1,5]), ","))
  resistance<-unlist(strsplit(as.character(common_sample[1,4]), ","))
  result1<-matrix(0,length(ccc),6)
for(i in 1:length(ccc)){
  result1[i,1]<-G
  file_name<-ccc[i]
  result1[i,2]<-file_name
  common_sample1<-as.data.frame(tresult_want[match(file_name,tresult_want[,1]),])
  sensitivity1<-unlist(strsplit(as.character(common_sample1[1,5]), ","))
  resistance1<-unlist(strsplit(as.character(common_sample1[1,4]), ","))
  int_sen<-intersect(sensitivity,resistance1)
  int_res<-intersect(resistance,sensitivity1)
  result1[i,3]<-length(sensitivity)
  result1[i,4]<-length(resistance)
  result1[i,5]<-length(int_sen)
  result1[i,6]<-length(int_res)
}
  result<-rbind(result,result1)
}


###########################新想法 从病人角度分对化疗药物和HDACI敏感的病人  ############
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
ccc<-c(as.character(HDAC1[,1]),"S1848","S2759","S1194","S1047")
bbb<-c("5-FU","GEM","IRI","OXA","PAC")

sen_chemo<-NULL
res_chemo<-NULL
for(j in 1:5){
  G<-bbb[j]
  common_sample<-as.data.frame(tresult_want[match(G,tresult_want[,1]),])
  sensitivity<-unlist(strsplit(as.character(common_sample[1,5]), ","))
  resistance<-unlist(strsplit(as.character(common_sample[1,4]), ","))
  sen_chemo<-union(sen_chemo,sensitivity)
  res_chemo<-union(res_chemo,resistance)
}
############################
sen_HDAC<-NULL
res_HDAC<-NULL
ccc<-c(as.character(HDAC1[,1]),"S1848","S2759","S1194","S1047")
for(j in 1:length(ccc)){
  G<-ccc[j]
  common_sample<-as.data.frame(tresult_want[match(G,tresult_want[,1]),])
  sensitivity<-unlist(strsplit(as.character(common_sample[1,5]), ","))
  resistance<-unlist(strsplit(as.character(common_sample[1,4]), ","))
  sen_HDAC<-union(sen_HDAC,sensitivity)
  res_HDAC<-union(res_HDAC,resistance)
}

############3
bbb<-c("5-FU","GEM","IRI","OXA","PAC")
result_chemo<-matrix(0,length(bbb),length(colnames(ic50)))
for(j in 1:length(bbb)){
  G<-bbb[j]
  common_sample<-as.data.frame(tresult_want[match(G,tresult_want[,1]),])
  sensitivity<-unlist(strsplit(as.character(common_sample[1,5]), ","))
  resistance<-unlist(strsplit(as.character(common_sample[1,4]), ","))
  result_chemo[j,match(sensitivity,colnames(ic50))]<-1
}
sen_sum<-apply(result_chemo,2,sum)
colnames(result_chemo)<-colnames(ic50)
sensitive1<-colnames(result_chemo)[which((sen_sum>2))]
resistant1<-colnames(result_chemo)[-match(sensitive1,colnames(result_chemo))]

################################## 在至少一半药物以上都敏感才定义为敏感
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
ccc<-c(as.character(HDAC1[,1]),"S1848","S2759","S1194","S1047")
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
sensitivity<-colnames(result)[which((sen_sum>=10))] #### 19 samples
resistance<-colnames(result)[-match(sensitive,colnames(result))]

################## 筛选HDAC's biomarker ###########################
#######################    ATACseq   ################################################
#############  DApeaks
rm(list=ls())
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_count_name<-read.table("peak_count.txt",sep="\t",header=T) #没有全0的
setwd("~/xjj/drug/drug_result/chemotherapy/0_clinical")

sensitivity<-colnames(result)[which((sen_sum>=10))] #### 19 samples
resistance<-colnames(result)[-match(sensitive,colnames(result))]

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
setwd("~/xjj/drug/drug_result/chemotherapy/2_ATACseq")
write.table(res,"sensitivity-resistance-all-DESeq2_chemotherapy_AUC_ATAC.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
#write.table(exprSet_new,"sensitivity-resistance-HDAC-DESeq2_normlization.txt",sep = '\t',col.names = T,row.names = T,quote = FALSE)##数据输出

#setwd("~/xjj/drug/drug_result/HDAC_drug/heatmap")
res<-read.table("sensitivity-resistance-all-DESeq2_HDAC_chemotherapy_AUC_ATAC.txt",header = TRUE,sep = "\t")
#resSig<-res[which(res$padj<0.05 & abs(res$log2FoldChange>1)),]
resSig<-res[which(res$pvalue<0.05 & abs(res$log2FoldChange)>1),]
resSig<-res[which(res$pvalue<0.001),]

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
out_file<-"sen-res-chemotherapy_AUC-0.05-P"
write.table(peaks_up,file=paste(out_file,"-logFC1-DESeq2-up.txt",sep=""),sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')
write.table(d_up,file=paste(out_file,"-logFC1-DESeq2-up.gff",sep=""),sep = '\t',
            col.names = F,row.names = F,quote = FALSE)
a_down<-apply(down_gene,1,function(x) unlist(strsplit(as.character(x), "_")))
d_down<-t(a_down)
e_down<-data.frame(rep("+",nrow(d_down)))
peaks_down<-cbind(down_gene,d_down,e_down)
colnames(peaks_down)<-c("peak_id","chr","start","end","strand")
write.table(peaks_down,file=paste(out_file,"-logFC1-DESeq2-down.txt",sep=""),sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')
write.table(d_down,file=paste(out_file,"-logFC1-DESeq2-down.gff",sep=""),sep = '\t',
            col.names = F,row.names = F,quote = FALSE)
#####  bed
f_up<-data.frame(rep(1,nrow(d_up)))
peaks_bed_up<-cbind(d_up,up_gene,f_up,e_up)
colnames(peaks_bed_up)<-c("chr","start","end","peak_id","score","strand")
write.table(peaks_bed_up,file=paste(out_file,"-logFC1-DESeq2-up.bed",sep=""),sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')
f_down<-data.frame(rep(1,nrow(d_down)))
peaks_bed_down<-cbind(d_down,down_gene,f_down,e_down)
colnames(peaks_bed_down)<-c("chr","start","end","peak_id","score","strand")
write.table(peaks_bed_down,file=paste(out_file,"-logFC1-DESeq2-down.bed",sep=""),sep = '\t',
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
setwd("~/xjj/drug/drug_result/chemotherapy/2_ATACseq")
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







