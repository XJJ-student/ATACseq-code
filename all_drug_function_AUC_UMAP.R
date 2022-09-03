#### 5-cross-validation
###function of cross-validation
rm(list=ls())
library(plyr)
CVgroup <- function(k,datasize,seed){
  cvlist <- list()
  set.seed(seed)
  n <- rep(1:k,ceiling(datasize/k))[1:datasize]    #将数据分成K份，并生成的完成数据集n
  temp <- sample(n,datasize)   #把n打乱
  x <- 1:k
  dataseq <- 1:datasize
  cvlist <- lapply(x,function(x) dataseq[temp==x])  #dataseq中随机生成k个随机有序数据列
  return(cvlist)
}

setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/6_application/data")
peak_count_name<-read.table("sample37_peaks_count.txt",sep="\t",header=T)
peak_RPKM_name<-read.table("sample37_peaks_RPKM.txt",sep="\t",header=T)

setwd("~/xjj/drug/drug_result/drug64_AUC")
tresult_want<-read.table("drug64_sample.txt",sep="\t",header=T)
library(DESeq2)
library(dplyr)
library("data.table")
library(foreach)

#setwd("~/xjj/drug")
#ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
#ic50<-ic501[,!(colnames(ic501)%in%c("PC.100","PC.34"))]
#drug_info<-read.csv("compiled_drug_info.csv",header = TRUE,sep = ",")
#drug_result<-drug_info[match(rownames(ic50),drug_info$drug_id),]
#HDAC1<-drug_result[drug_result$target=="HDAC",]
#aaa<-c("S1848","S2759","S1194","S1047")
#bbb<-c("5-FU","GEM","IRI","OXA","PAC")
#ccc<-c(aaa,bbb)
#HDAC2<-drug_result[drug_result$drug_id%in%ccc,]
#HDAC3<-rbind(HDAC1,HDAC2)
#HDAC<-HDAC3[order(HDAC3[,1],decreasing = T),]
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/6_application/data/data")
first_category_name = list.files("fold-second5-100") 
aa<-unique(substring(first_category_name,1,5))
b<-aa[-c(1:9)]
a<-c("5-FU",  "GEM","S1030", "S1090", "S1181","S1249","S1451", "S1515", "S1826", "S1848", "S2018", "S2244","S4900")
a<-c("5-FU","GEM","IRI","OXA","PAC",b)

tresult_want2<-tresult_want1[-(match(a,tresult_want1[,1])),]
tresult_want<-tresult_want2[-c(1:10),]

#tresult_want<-tresult_want1[-c(1:13),]


library(doParallel)
library(parallel)
library(iterators)

cl <- makeCluster(40)
registerDoParallel(cl)
#sample_label<-NULL
#foreach(k=1:nrow(tresult_want),.combine='rbind',.packages = 'DESeq2')%dopar%{
foreach(k=1:nrow(tresult_want),.combine='rbind')%dopar%{
  file_name=as.character(tresult_want[k,1])
  cat(file_name,"\n")
  common_sample<-as.data.frame(tresult_want[match(file_name,tresult_want[,1]),])
  sensitivity<-unlist(strsplit(as.character(common_sample[1,5]), ","))
  resistance<-unlist(strsplit(as.character(common_sample[1,4]), ","))
  #cross-validation
  cv_resistance <- CVgroup(k = 5,datasize = length(resistance),seed = 12345)# first 1206, second 100, third 1005, fourth 50.fiveth 12345
  cv_sensitivity <- CVgroup(k = 5,datasize = length(sensitivity),seed = 12345)
  
  sample_label1<-matrix(0,5,5)
  for(i in 1:5){
    resistance_train<-resistance[-cv_resistance[[i]]]
    resistance_validation<-resistance[cv_resistance[[i]]]
    sensitivity_train<-sensitivity[-cv_sensitivity[[i]]]
    sensitivity_validation<-sensitivity[cv_sensitivity[[i]]]
    sample_label1[i,1]<-file_name
    sample_label1[i,2]<-paste0(resistance_train,collapse=",")
    sample_label1[i,3]<-paste0(sensitivity_train,collapse=",")
    sample_label1[i,4]<-paste0(resistance_validation,collapse=",")
    sample_label1[i,5]<-paste0(sensitivity_validation,collapse=",")
    resistance_peak<-peak_count_name[,match(resistance_train,colnames(peak_count_name))]
    sensitivity_peak<-peak_count_name[,match(sensitivity_train,colnames(peak_count_name))]
    peak_name<-peak_count_name[,1]
    library(DESeq2)
    colDate<-data.frame(row.names = c(as.vector(resistance_train),as.vector(sensitivity_train)),
                        condition=factor(c(rep("resistance",length(resistance_train)),rep("sensitivity",length(sensitivity_train)))))
    datexpr<-cbind(resistance_peak,sensitivity_peak)##前面的相对于后面的上下调
    counts <- apply(datexpr,2,as.numeric)   ###矩阵中必须是数值 
    dds<-DESeqDataSetFromMatrix(countData = counts,colData = colDate,design = ~condition)
    dds<-DESeq(dds)##进行标准化分析
    #sizeFactors(dds)##查看每个主成分的标准化值
    res<-cbind(peak_name,as.data.frame(results(dds)))##将结果输出
    colnames(res)<- c('peak_id',"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj" )##对列名进行修改
    resSig<-res[order(res$pvalue)[1:100],1]  ## 这里需要改Top500，改成500
    sensitivity_train_top<-peak_RPKM_name[match(resSig,peak_RPKM_name[,1]),match(sensitivity_train,colnames(peak_RPKM_name))]
    sensitivity_validation_top<-peak_RPKM_name[match(resSig,peak_RPKM_name[,1]),match(sensitivity_validation,colnames(peak_RPKM_name))]
    resistance_train_top<-peak_RPKM_name[match(resSig,peak_RPKM_name[,1]),match(resistance_train,colnames(peak_RPKM_name))]
    resistance_validation_top<-peak_RPKM_name[match(resSig,peak_RPKM_name[,1]),match(resistance_validation,colnames(peak_RPKM_name))]
    
    outputPath = paste('~/xjj/drug/drug_result/HDAC20_chemotherapy5/6_application/data/data/fold-fifth',i,"-100",sep='')
    if (!file.exists(outputPath)){
      dir.create(outputPath, recursive=TRUE)
    }
    write.table(sensitivity_train_top,paste(outputPath,'/',file_name,"-sensitivity-train-",i,"-100.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
    write.table(sensitivity_validation_top,paste(outputPath,'/',file_name,"-sensitivity-validation-",i,"-100.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
    write.table(resistance_train_top,paste(outputPath,'/',file_name,"-resistance-train-",i,"-100.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
    write.table(resistance_validation_top,paste(outputPath,'/',file_name,"-resistance-validation-",i,"-100.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出}
    #write.table(res,paste(outputPath,'/',file_name,"-DARs-",i,".txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  
    resSig5<-res[order(res$pvalue)[1:500],1]  ## 这里需要改Top500，改成500
    sensitivity_train_top5<-peak_RPKM_name[match(resSig5,peak_RPKM_name[,1]),match(sensitivity_train,colnames(peak_RPKM_name))]
    sensitivity_validation_top5<-peak_RPKM_name[match(resSig5,peak_RPKM_name[,1]),match(sensitivity_validation,colnames(peak_RPKM_name))]
    resistance_train_top5<-peak_RPKM_name[match(resSig5,peak_RPKM_name[,1]),match(resistance_train,colnames(peak_RPKM_name))]
    resistance_validation_top5<-peak_RPKM_name[match(resSig5,peak_RPKM_name[,1]),match(resistance_validation,colnames(peak_RPKM_name))]
    
    outputPath5 = paste('~/xjj/drug/drug_result/HDAC20_chemotherapy5/6_application/data/data/fold-fifth',i,"-500",sep='')
    if (!file.exists(outputPath5)){
      dir.create(outputPath5, recursive=TRUE)
    }
    write.table(sensitivity_train_top5,paste(outputPath5,'/',file_name,"-sensitivity-train-",i,"-500.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
    write.table(sensitivity_validation_top5,paste(outputPath5,'/',file_name,"-sensitivity-validation-",i,"-500.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
    write.table(resistance_train_top5,paste(outputPath5,'/',file_name,"-resistance-train-",i,"-500.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
    write.table(resistance_validation_top5,paste(outputPath5,'/',file_name,"-resistance-validation-",i,"-500.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
    
    resSig10<-res[order(res$pvalue)[1:1000],1]  ## 这里需要改Top500，改成500
    sensitivity_train_top10<-peak_RPKM_name[match(resSig10,peak_RPKM_name[,1]),match(sensitivity_train,colnames(peak_RPKM_name))]
    sensitivity_validation_top10<-peak_RPKM_name[match(resSig10,peak_RPKM_name[,1]),match(sensitivity_validation,colnames(peak_RPKM_name))]
    resistance_train_top10<-peak_RPKM_name[match(resSig10,peak_RPKM_name[,1]),match(resistance_train,colnames(peak_RPKM_name))]
    resistance_validation_top10<-peak_RPKM_name[match(resSig10,peak_RPKM_name[,1]),match(resistance_validation,colnames(peak_RPKM_name))]
    
    outputPath10 = paste('~/xjj/drug/drug_result/HDAC20_chemotherapy5/6_application/data/data/fold-fifth',i,"-1000",sep='')
    if (!file.exists(outputPath10)){
      dir.create(outputPath5, recursive=TRUE)
    }
    write.table(sensitivity_train_top10,paste(outputPath5,'/',file_name,"-sensitivity-train-",i,"-1000.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
    write.table(sensitivity_validation_top10,paste(outputPath5,'/',file_name,"-sensitivity-validation-",i,"-1000.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
    write.table(resistance_train_top10,paste(outputPath5,'/',file_name,"-resistance-train-",i,"-1000.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
    write.table(resistance_validation_top10,paste(outputPath5,'/',file_name,"-resistance-validation-",i,"-1000.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
     }
    #sample_label<-rbind(sample_label,sample_label1)
}
stopCluster(cl)
colnames(sample_label)<-c("drug","resistance_train","sensitivity_train","resistance_validation","sensitivity_validation")

######################


setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_count_name<-read.table("peak_count.txt",sep="\t",header=T)
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
colnames(ic50)
resistance_peak111<-peak_count_name[,match(colnames(ic50),gsub("\\.","-",colnames(peak_count_name)))]
resistance_peak112<-cbind(peak_count_name[,1],resistance_peak111)
colnames(resistance_peak111)<-sample_id[match(gsub("\\.","-",colnames(resistance_peak111)),sample_id[,1]),2]
colnames(resistance_peak112)<-c("peaks_name",colnames(resistance_peak111))
write.table(resistance_peak112,"sample37_peaks_count.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出



library(plyr)
CVgroup <- function(k,datasize,seed){
  cvlist <- list()
  set.seed(seed)
  n <- rep(1:k,ceiling(datasize/k))[1:datasize]    #将数据分成K份，并生成的完成数据集n
  temp <- sample(n,datasize)   #把n打乱
  x <- 1:k
  dataseq <- 1:datasize
  cvlist <- lapply(x,function(x) dataseq[temp==x])  #dataseq中随机生成k个随机有序数据列
  return(cvlist)
}
k <- 5
datasize <- 37
cvlist <- CVgroup(k = k,datasize = datasize,seed = 1206)
cvlist


#############   guiding treatment of pancreatic cancer   ################
####################     uwot    AUC渐变     #####################
rm(list=ls())
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_count_name<-read.table("peak_RPKM.txt",sep="\t",header=T)
setwd("~/xjj/drug/drug_result/drug64_AUC")
tresult_want<-read.table("drug64_sample.txt",sep="\t",header=T)
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
library(DESeq2)
library(dplyr)
library("data.table")
library(foreach)
setwd("~/xjj/drug")
ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,-(which(colnames(ic501)==c("PC.100","PC.34")))]
colnames(ic50)<-sample_id[match(colnames(ic50),sample_id[,2]),1]

foreach(k=1:nrow(tresult_want),.combine='rbind')%do%{
  G=as.character(tresult_want[k,1])
  cat(G,"\n")
  file_name=as.character(tresult_want[k,1])
  common_sample<-as.data.frame(tresult_want[match(G,tresult_want[,1]),])
  sensitivity<-unlist(strsplit(as.character(common_sample[1,5]), ","))
  resistance<-unlist(strsplit(as.character(common_sample[1,4]), ","))
  resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
  sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]
  resistance_peak<-peak_count_name[,match(resistance_new,gsub("\\.","-",colnames(peak_count_name)))]
  sensitivity_peak<-peak_count_name[,match(sensitivity_new,gsub("\\.","-",colnames(peak_count_name)))]
  peak_name<-peak_count_name[,1]
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/2_ATACseq_DApeaks",sep=""))
  res<-fread(paste("sensitivity-resistance-all-DESeq2_",file_name,"_ATAC.txt",sep=''),header=T,data.table=F)
  resSig<-res[order(res$pvalue)[1:1000],1]   ##features
  ##前19耐药，后18敏感
  peak<-cbind(resistance_peak,sensitivity_peak)
  ic502<-ic50[rownames(ic50)%in%file_name,]
  AUC<-round(ic502[match(gsub("\\.","-",colnames(peak)),colnames(ic502))],2)
  feature_peak<-peak[match(resSig,peak_name),]
  feature_peak1<-t(feature_peak)
  labelll<-as.data.frame(matrix(c(rep("resistance",ncol(resistance_peak)),rep("sensitivity",ncol(sensitivity_peak))),ncol=1))
  expbbb<-cbind(feature_peak1,labelll)
  colnames(expbbb)<-c(peak_name[match(resSig,peak_name)],"cancertype")
  library(uwot)
  set.seed(100)
  iris_umap <- uwot::umap(expbbb,n_neighbors =10)
  iris_sumap_res <- data.frame(iris_umap,cancertype=expbbb$cancertype,AUC=t(AUC))
  colnames(iris_sumap_res)<-c("UMAP1","UMAP2","cancertype","AUC")
  head(iris_sumap_res)
  library(ggplot2)
  p<-ggplot(iris_sumap_res,aes(UMAP1,UMAP2,color=AUC,shape=cancertype)) + 
    geom_point(size = 3) +theme_bw() + 
    scale_color_gradient(low="green", high="red")+
    theme(plot.title = element_text(hjust = 0.5)) + 
    labs(x="UMAP_1",y="UMAP_2",
         title = paste("A UMAP visualization of the ",file_name,sep=""))+
    geom_text(aes(label = AUC),colour="black",size = 2.5,check_overlap = TRUE,nudge_y = (-0.1))
  setwd("~/xjj/drug/drug_result/drug64_AUC/figure_feature/try")
  ggsave(paste(file_name,"_Top1000_DApeaks.pdf",sep=""),p,width = 6, height = 5)
}

###############  敏感和耐药两个颜色  ##################################
foreach(k=1:nrow(tresult_want),.combine='rbind')%do%{
  G=as.character(tresult_want[k,1])
  cat(G,"\n")
  file_name=as.character(tresult_want[k,1])
  common_sample<-as.data.frame(tresult_want[match(G,tresult_want[,1]),])
  sensitivity<-unlist(strsplit(as.character(common_sample[1,5]), ","))
  resistance<-unlist(strsplit(as.character(common_sample[1,4]), ","))
  resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
  sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]
  resistance_peak<-peak_count_name[,match(resistance_new,gsub("\\.","-",colnames(peak_count_name)))]
  sensitivity_peak<-peak_count_name[,match(sensitivity_new,gsub("\\.","-",colnames(peak_count_name)))]
  peak_name<-peak_count_name[,1]
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/2_ATACseq_DApeaks",sep=""))
  res<-fread(paste("sensitivity-resistance-all-DESeq2_",file_name,"_ATAC.txt",sep=''),header=T,data.table=F)
  resSig<-res[order(res$pvalue)[1:1000],1]   ##features
  ##前19耐药，后18敏感
  peak<-cbind(resistance_peak,sensitivity_peak)
  ic502<-ic50[rownames(ic50)%in%file_name,]
  AUC<-round(ic502[match(gsub("\\.","-",colnames(peak)),colnames(ic502))],2)
  feature_peak<-peak[match(resSig,peak_name),]
  feature_peak1<-t(feature_peak)
  labelll<-as.data.frame(matrix(c(rep("resistance",ncol(resistance_peak)),rep("sensitivity",ncol(sensitivity_peak))),ncol=1))
  expbbb<-cbind(feature_peak1,labelll)
  colnames(expbbb)<-c(peak_name[match(resSig,peak_name)],"cancertype")
  library(uwot)
  set.seed(100)
  iris_umap <- uwot::umap(expbbb)
  iris_sumap_res <- data.frame(iris_umap,cancertype=expbbb$cancertype,AUC=t(AUC))
  colnames(iris_sumap_res)<-c("UMAP1","UMAP2","cancertype","AUC")
  head(iris_sumap_res)
  library(ggplot2)
  p<-ggplot(iris_sumap_res,aes(UMAP1,UMAP2,color=cancertype)) + 
    geom_point(size = 2) +theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    labs(x="UMAP_1",y="UMAP_2",
         title = paste("A UMAP visualization of the ",file_name,sep=""))+
    geom_text(aes(label = AUC),colour="black",size = 2.5,check_overlap = TRUE,nudge_y = (-0.15))
  setwd("~/xjj/drug/drug_result/drug64_AUC/figure_feature/Pmin_Top1000_DApeaks")
  ggsave(paste(file_name,"_Top1000_DApeaks.pdf",sep=""),p,width = 6, height = 5)
}




### log2FC>1

foreach(k=1:nrow(tresult_want),.combine='rbind')%do%{
  G=as.character(tresult_want[k,1])
  cat(G,"\n")
  file_name=as.character(tresult_want[k,1])
  common_sample<-as.data.frame(tresult_want[match(G,tresult_want[,1]),])
  sensitivity<-unlist(strsplit(as.character(common_sample[1,5]), ","))
  resistance<-unlist(strsplit(as.character(common_sample[1,4]), ","))
  resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
  sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]
  resistance_peak<-peak_count_name[,match(resistance_new,gsub("\\.","-",colnames(peak_count_name)))]
  sensitivity_peak<-peak_count_name[,match(sensitivity_new,gsub("\\.","-",colnames(peak_count_name)))]
  peak_name<-peak_count_name[,1]
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/2_ATACseq_DApeaks",sep=""))
  res<-fread(paste("sensitivity-resistance-all-DESeq2_",file_name,"_ATAC.txt",sep=''),header=T,data.table=F)
  res1<-res[abs(res$log2FoldChange)>1,]
  resSig<-res1[order(res1$pvalue)[1:1000],1]   ##features
  ##前19耐药，后18敏感
  peak<-cbind(resistance_peak,sensitivity_peak)
  feature_peak<-peak[match(resSig,peak_name),]
  feature_peak1<-t(feature_peak)
  labelll<-as.data.frame(matrix(c(rep("resistance",ncol(resistance_peak)),rep("sensitivity",ncol(sensitivity_peak))),ncol=1))
  expbbb<-cbind(feature_peak1,labelll)
  colnames(expbbb)<-c(peak_name[match(resSig,peak_name)],"cancertype")
  library(uwot)
  iris_umap <- uwot::umap(expbbb)
  iris_sumap_res <- data.frame(iris_umap,cancertype=expbbb$cancertype)
  head(iris_sumap_res)
  library(ggplot2)
  p<-ggplot(iris_sumap_res,aes(X1,X2,color=cancertype)) + 
    geom_point(size = 1) +theme_bw() + 
    #geom_hline(yintercept = 0,lty=2,col="red") + 
    #geom_vline(xintercept = 0,lty=2,col="blue",lwd=1) +
    theme(plot.title = element_text(hjust = 0.5)) + 
    labs(x="UMAP_1",y="UMAP_2",
         title = paste("A UMAP visualization of the ",file_name,sep=""))
  #geom_text(aes(label = cancertype),size = 3,check_overlap = TRUE,angle = 30)
  setwd("~/xjj/drug/drug_result/drug64_AUC/figure_feature/log2FC_Pmin_Top1000_DApeaks")
  ggsave(paste(file_name,"_log2FC1_P_Top1000_DApeaks.pdf",sep=""),p,width = 6, height = 5)
}
################  umap packages  #############
rm(list=ls())
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_count_name<-read.table("peak_RPKM.txt",sep="\t",header=T)
setwd("~/xjj/drug/drug_result/drug64_AUC")
tresult_want<-read.table("drug64_sample.txt",sep="\t",header=T)
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
library(DESeq2)
library(dplyr)
library("data.table")
library(foreach)

foreach(k=1:nrow(tresult_want),.combine='rbind')%do%{
  G=as.character(tresult_want[k,1])
  cat(G,"\n")
  file_name=as.character(tresult_want[k,1])
  common_sample<-as.data.frame(tresult_want[match(G,tresult_want[,1]),])
  sensitivity<-unlist(strsplit(as.character(common_sample[1,5]), ","))
  resistance<-unlist(strsplit(as.character(common_sample[1,4]), ","))
  resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
  sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]
  resistance_peak<-peak_count_name[,match(resistance_new,gsub("\\.","-",colnames(peak_count_name)))]
  sensitivity_peak<-peak_count_name[,match(sensitivity_new,gsub("\\.","-",colnames(peak_count_name)))]
  peak_name<-peak_count_name[,1]
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/2_ATACseq_DApeaks",sep=""))
  res<-fread(paste("sensitivity-resistance-all-DESeq2_",file_name,"_ATAC.txt",sep=''),header=T,data.table=F)
  resSig<-res[order(res$pvalue)[1:1000],1]   ##features
  ##前19耐药，后18敏感
  peak<-cbind(resistance_peak,sensitivity_peak)
  feature_peak<-peak[match(resSig,peak_name),]
  feature_peak1<-t(feature_peak)
  labelll<-as.data.frame(matrix(c(rep("resistance",ncol(resistance_peak)),rep("sensitivity",ncol(sensitivity_peak))),ncol=1))
  expbbb<-cbind(feature_peak1,labelll)
  colnames(expbbb)<-c(peak_name[match(resSig,peak_name)],"cancertype")
  library(umap)
  iris.data = expbbb[,-(ncol(expbbb))]
  iris.umap = umap::umap(iris.data)
  head(iris.umap$layout)
  iris_sumap_res <- data.frame(iris.umap$layout,cancertype=expbbb$cancertype)
  head(iris_sumap_res)
  library(ggplot2)
  p<-ggplot(iris_sumap_res,aes(X1,X2,color=cancertype)) + 
    geom_point(size = 1) +theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    labs(x="UMAP_1",y="UMAP_2",
         title = paste("A UMAP visualization of the ",file_name,sep=""))
  setwd("~/xjj/drug/drug_result/drug64_AUC/figure_feature/umap_packages_Pmin_Top1000_DApeaks")
  ggsave(paste(file_name,"_umap_packages_Pmin_Top1000_DApeaks.pdf",sep=""),p,width = 6, height = 5)
}



setwd("~/xjj/drug/drug_result/HDACi_chemo618")
iris.data<-read.csv("covid_and_health_matrix.csv",header = FALSE,sep = "\t")
library(umap)
iris.umap = umap::umap(iris.data,n_components=2,input=dist)
head(iris.umap$layout)
iris_sumap_res <- data.frame(iris.umap$layout)
head(iris_sumap_res)
write.table(iris_sumap_res,"iris_sumap_res.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出


##############   t-SNE  #################
foreach(k=1:nrow(tresult_want),.combine='rbind')%do%{
  G=as.character(tresult_want[k,1])
  cat(G,"\n")
  file_name=as.character(tresult_want[k,1])
  common_sample<-as.data.frame(tresult_want[match(G,tresult_want[,1]),])
  sensitivity<-unlist(strsplit(as.character(common_sample[1,5]), ","))
  resistance<-unlist(strsplit(as.character(common_sample[1,4]), ","))
  resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
  sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]
  resistance_peak<-peak_count_name[,match(resistance_new,gsub("\\.","-",colnames(peak_count_name)))]
  sensitivity_peak<-peak_count_name[,match(sensitivity_new,gsub("\\.","-",colnames(peak_count_name)))]
  peak_name<-peak_count_name[,1]
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/2_ATACseq_DApeaks",sep=""))
  res<-fread(paste("sensitivity-resistance-all-DESeq2_",file_name,"_ATAC.txt",sep=''),header=T,data.table=F)
  resSig<-res[order(res$pvalue)[1:1000],1]   ##features
  ##前19耐药，后18敏感
  peak<-cbind(resistance_peak,sensitivity_peak)
  ic502<-ic50[rownames(ic50)%in%file_name,]
  AUC<-round(ic502[match(gsub("\\.","-",colnames(peak)),colnames(ic502))],2)
  feature_peak<-peak[match(resSig,peak_name),]
  feature_peak1<-t(feature_peak)
  labelll<-as.data.frame(matrix(c(rep("resistance",ncol(resistance_peak)),rep("sensitivity",ncol(sensitivity_peak))),ncol=1))
  expbbb<-cbind(feature_peak1,labelll)
  colnames(expbbb)<-c(peak_name[match(resSig,peak_name)],"cancertype")
  set.seed(42)
  library(Rtsne)
  tsne_out <- Rtsne(as.matrix(expbbb[,-(ncol(expbbb))]),pca=FALSE,dims=2,
                    perplexity=10,theta=0.0) # Run TSNE
  iris_sumap_res <- data.frame(tsne_out$Y,cancertype=expbbb$cancertype,AUC=t(AUC))
  colnames(iris_sumap_res)<-c("tSNE1","tSNE2","cancertype","AUC")
  head(iris_sumap_res)
  library(ggplot2)
  p<-ggplot(iris_sumap_res,aes(tSNE1,tSNE2,color=AUC,shape=cancertype)) + 
    geom_point(size = 3) +theme_bw() + 
    scale_color_gradient(low="green", high="red")+
    theme(plot.title = element_text(hjust = 0.5)) + 
    labs(x="tSNE_1",y="tSNE_2",
         title = paste("A tSNE visualization of the ",file_name,sep=""))+
    geom_text(aes(label = AUC),colour="black",size = 2.5,check_overlap = TRUE,nudge_y = (-3))
  setwd("~/xjj/drug/drug_result/drug64_AUC/figure_feature/Rtsne_packages_Pmin_Top1000_DApeaks_Gradual_change")
  ggsave(paste(file_name,"_Rtsne_packages_Pmin_Top1000_DApeaks.pdf",sep=""),p,width = 6, height = 5)
}


########################   PCA   #######################
foreach(k=1:nrow(tresult_want),.combine='rbind')%do%{
  G=as.character(tresult_want[k,1])
  cat(G,"\n")
  file_name=as.character(tresult_want[k,1])
  common_sample<-as.data.frame(tresult_want[match(G,tresult_want[,1]),])
  sensitivity<-unlist(strsplit(as.character(common_sample[1,5]), ","))
  resistance<-unlist(strsplit(as.character(common_sample[1,4]), ","))
  resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
  sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]
  resistance_peak<-peak_count_name[,match(resistance_new,gsub("\\.","-",colnames(peak_count_name)))]
  sensitivity_peak<-peak_count_name[,match(sensitivity_new,gsub("\\.","-",colnames(peak_count_name)))]
  peak_name<-peak_count_name[,1]
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/2_ATACseq_DApeaks",sep=""))
  res<-fread(paste("sensitivity-resistance-all-DESeq2_",file_name,"_ATAC.txt",sep=''),header=T,data.table=F)
  resSig<-res[order(res$pvalue)[1:1000],1]   ##features
  ##前19耐药，后18敏感
  peak<-cbind(resistance_peak,sensitivity_peak)
  ic502<-ic50[rownames(ic50)%in%file_name,]
  AUC<-round(ic502[match(gsub("\\.","-",colnames(peak)),colnames(ic502))],2)
  feature_peak<-peak[match(resSig,peak_name),]
  feature_peak1<-t(feature_peak)
  labelll<-as.data.frame(matrix(c(rep("resistance",ncol(resistance_peak)),rep("sensitivity",ncol(sensitivity_peak))),ncol=1))
  expbbb<-cbind(feature_peak1,labelll)
  colnames(expbbb)<-c(peak_name[match(resSig,peak_name)],"cancertype")
  set.seed(1000)
  df_pca <- prcomp(as.matrix(expbbb[,-(ncol(expbbb))])) # Run PCA
  df_pcs <-data.frame(df_pca$x,cancertype=expbbb$cancertype)
  iris_umap<- data.frame(df_pcs[,c(1,2)])
  iris_sumap_res <- data.frame(df_pcs[,c(1,2,38)],AUC=t(AUC))
  colnames(iris_sumap_res)<-c("PC1","PC2","cancertype","AUC")
  head(iris_sumap_res)
  library(ggplot2)
  p<-ggplot(iris_sumap_res,aes(PC1,PC2,color=AUC,shape=cancertype)) + 
    geom_point(size = 3) +theme_bw() + 
    scale_color_gradient(low="green", high="red")+
    theme(plot.title = element_text(hjust = 0.5)) + 
    labs(x="PCA_1",y="PCA_2",
         title = paste("A PCA visualization of the ",file_name,sep=""))+
    geom_text(aes(label = AUC),colour="black",size = 2.5,check_overlap = TRUE,nudge_y = (-0.5))
  
  setwd("~/xjj/drug/drug_result/drug64_AUC/figure_feature/PCA_Pmin_Top1000_DApeaks_Gradual_change")
  ggsave(paste(file_name,"_PCA_Pmin_Top1000_DApeaks.pdf",sep=""),p,width = 6, height = 5)
}


##################   MDS    ##############################
library(stats)
library(ggplot2)
foreach(k=1:nrow(tresult_want),.combine='rbind')%do%{
  G=as.character(tresult_want[k,1])
  cat(G,"\n")
  file_name=as.character(tresult_want[k,1])
  common_sample<-as.data.frame(tresult_want[match(G,tresult_want[,1]),])
  sensitivity<-unlist(strsplit(as.character(common_sample[1,5]), ","))
  resistance<-unlist(strsplit(as.character(common_sample[1,4]), ","))
  resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
  sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]
  resistance_peak<-peak_count_name[,match(resistance_new,gsub("\\.","-",colnames(peak_count_name)))]
  sensitivity_peak<-peak_count_name[,match(sensitivity_new,gsub("\\.","-",colnames(peak_count_name)))]
  peak_name<-peak_count_name[,1]
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/2_ATACseq_DApeaks",sep=""))
  res<-fread(paste("sensitivity-resistance-all-DESeq2_",file_name,"_ATAC.txt",sep=''),header=T,data.table=F)
  resSig<-res[order(res$pvalue)[1:1000],1]   ##features
  ##前19耐药，后18敏感
  peak<-cbind(resistance_peak,sensitivity_peak)
  ic502<-ic50[rownames(ic50)%in%file_name,]
  AUC<-round(ic502[match(gsub("\\.","-",colnames(peak)),colnames(ic502))],2)
  feature_peak<-peak[match(resSig,peak_name),]
  feature_peak1<-t(feature_peak)
  labelll<-as.data.frame(matrix(c(rep("resistance",ncol(resistance_peak)),rep("sensitivity",ncol(sensitivity_peak))),ncol=1))
  expbbb<-cbind(feature_peak1,labelll)
  colnames(expbbb)<-c(peak_name[match(resSig,peak_name)],"cancertype")
  set.seed(1000)
  dis_iris = dist(expbbb[,-(ncol(expbbb))],p = 2)
  mds_x = cmdscale(dis_iris)
  mds_x = data.frame(mds_x)
  iris_sumap_res <- data.frame(mds_x,cancertype=expbbb$cancertype,AUC=t(AUC))
  colnames(iris_sumap_res)<-c("MDS1","MDS2","cancertype","AUC")
  head(iris_sumap_res)
  library(ggplot2)
  p<-ggplot(iris_sumap_res,aes(MDS1,MDS2,color=AUC,shape=cancertype)) + 
    geom_point(size = 3) +theme_bw() + 
    scale_color_gradient(low="green", high="red")+
    theme(plot.title = element_text(hjust = 0.5)) + 
    labs(x="MDS_1",y="MDS_2",
         title = paste("A MDS visualization of the ",file_name,sep=""))+
    geom_text(aes(label = AUC),colour="black",size = 2,check_overlap = TRUE,nudge_y = (-0.5))
  
  setwd("~/xjj/drug/drug_result/drug64_AUC/figure_feature/MDS_Pmin_Top1000_DApeaks_Gradual_change")
  ggsave(paste(file_name,"_MDS_Pmin_Top1000_DApeaks.pdf",sep=""),p,width = 6, height = 5)
}

ggplot(iris_sumap_res,aes(MDS1,MDS2,color=cancertype)) + 
  geom_point(size = 3) +theme_bw()
##############   两样本之间的相关性和欧氏距离   #######################
rm(list=ls())
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_count_name<-read.table("peak_RPKM.txt",sep="\t",header=T)
setwd("~/xjj/drug/drug_result/drug64_AUC")
tresult_want<-read.table("drug64_sample.txt",sep="\t",header=T)
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
library(dplyr)
library("data.table")
library(foreach)
library(philentropy)

drug_result<-NULL
cor_result<-matrix(0,64,3)
foreach(k=1:nrow(tresult_want),.combine='rbind')%do%{
  G=as.character(tresult_want[k,1])
  cat(G,"\n")
  file_name=as.character(tresult_want[k,1])
  common_sample<-as.data.frame(tresult_want[match(G,tresult_want[,1]),])
  sensitivity<-unlist(strsplit(as.character(common_sample[1,5]), ","))
  resistance<-unlist(strsplit(as.character(common_sample[1,4]), ","))
  resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
  sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]
  resistance_peak<-peak_count_name[,match(resistance_new,gsub("\\.","-",colnames(peak_count_name)))]
  sensitivity_peak<-peak_count_name[,match(sensitivity_new,gsub("\\.","-",colnames(peak_count_name)))]
  peak_name<-peak_count_name[,1]
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/2_ATACseq_DApeaks",sep=""))
  res<-fread(paste("sensitivity-resistance-all-DESeq2_",file_name,"_ATAC.txt",sep=''),header=T,data.table=F)
  resSig<-res[order(res$pvalue)[1:1000],1]   ##features
  
  peak<-cbind(resistance_peak,sensitivity_peak) 
  feature_peak<-peak[match(resSig,peak_name),]  #1000_peak_exp
  feature_peak1<-t(feature_peak)
  labelll<-as.data.frame(matrix(c(rep("resistance",ncol(resistance_peak)),rep("sensitivity",ncol(sensitivity_peak))),ncol=1))
  expbbb<-cbind(feature_peak1,labelll)
  colnames(expbbb)<-c(peak_name[match(resSig,peak_name)],"cancertype")
  
  #set.seed(42)
  #library(Rtsne)
  #tsne_out <- Rtsne(as.matrix(expbbb[,-(ncol(expbbb))]),pca=FALSE,dims=2,perplexity=10,theta=0.0)
  #iris_umap <- data.frame(tsne_out$Y)
  
  library(uwot)
  set.seed(100)
  iris_umap <- uwot::umap(expbbb)
  rownames(iris_umap)<-colnames(feature_peak)
  iris_umap<-as.data.frame(iris_umap)
  
  #set.seed(1000)
  #df_pca <- prcomp(as.matrix(expbbb[,-(ncol(expbbb))])) # Run PCA
  #df_pcs <-data.frame(df_pca$x,cancertype=expbbb$cancertype)
  #iris_umap<- data.frame(df_pcs[,c(1,2)])
  
  #set.seed(1000)
  #dis_iris = dist(expbbb[,-(ncol(expbbb))],p = 2)
  #mds_x = cmdscale(dis_iris)
  #iris_umap = data.frame(mds_x)
  
  
  ###两点之间的欧氏距离
  result<-NULL
  for(j in 1:(ncol(feature_peak)-1)){
    a<-feature_peak[,j]
    x<-iris_umap[j,]
    result1<-NULL
    for(m in (j+1):ncol(feature_peak)){
      b<-feature_peak[,m]
      y<-iris_umap[m,]
      dat<-rbind(x,y)
      d<-distance(dat,method = "euclidean")
      r<-cor.test(a,b,method ="pearson",alternative="two.sided")$estimate
      result2<-matrix(c(r,d),nrow=1)
      result1<-rbind(result1,result2)
    }
    result<-rbind(result,result1)
  }
  dat1<-as.data.frame(result)
  colnames(dat1)<-c("correlation","distance")
  pp<-ggplot(data=dat1, aes(x=correlation, y=distance))+geom_point(color="red",size=0.1)+stat_smooth(method="lm",se=FALSE)+stat_cor(data=dat1, method = "pearson")+
    ggtitle(paste(file_name,"-pearson","-UMAP_correlation_distance",sep="")) +
    theme(plot.title = element_text(hjust = 0.5))
  setwd("~/xjj/drug/drug_result/drug64_AUC/figure_feature/try")
  ggsave(paste("correlation_distance_",file_name,"_UMAP.pdf",sep=""),pp,width = 5, height = 5)
  #ggsave(paste("correlation_distance_",file_name,".pdf",sep=""),pp,width = 5, height = 5)
  
  colnames(result)<-c(paste(file_name,"_correlation",sep=""),paste(file_name,"_distance",sep=""))
  drug_result<-cbind(drug_result,result)
  cor_result[k,1]<-file_name
  cor_result[k,2]<-cor.test(result[,1],result[,2],method ="pearson",alternative="two.sided")$estimate
  cor_result[k,3]<-cor.test(result[,1],result[,2],method ="pearson",alternative="two.sided")$p.value
}
setwd("~/xjj/drug/drug_result/drug64_AUC/figure_feature/try")
colnames(cor_result)<-c("drug","pearson-estimate","pearson-pvalue")
write.table(cor_result,"MDS_correlation_distance.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(drug_result,"MDS_correlation_distance_detailed_information.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出




###################   HDAC     ######################
rm(list=ls())
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_count_name<-read.table("peak_RPKM.txt",sep="\t",header=T)
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
library(dplyr)
library("data.table")
library(foreach)
setwd("~/xjj/drug")
ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,-(which(colnames(ic501)==c("PC.100","PC.34")))]
colnames(ic50)<-sample_id[match(colnames(ic50),sample_id[,2]),1]
sensitivity<-c("PC.101","PC.109","PC.115","PC.116","PC.117","PC.13","PC.130","PC.136","PC.139","PC.16","PC.2","PC.27","PC.64","PC.78","PC.8","PC.81","PC.97","PC.98", "PC.L")
resistance<-c("PC.102","PC.104","PC.105","PC.111","PC.112","PC.119","PC.121","PC.134","PC.135","PC.14","PC.18","PC.22","PC.40","PC.5","PC.52","PC.56","PC.G","PC.I" )
resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]
resistance_peak<-peak_count_name[,match(resistance_new,gsub("\\.","-",colnames(peak_count_name)))]
sensitivity_peak<-peak_count_name[,match(sensitivity_new,gsub("\\.","-",colnames(peak_count_name)))]
peak_name<-peak_count_name[,1]
  
setwd("~/xjj/drug/drug_result/HDAC_frontiers_AUC_correlation/2_ATACseq_DApeaks")
res<-fread("sensitivity-resistance-all-DESeq2_HDAC_AUC_ATAC.txt",header=T,data.table=F)
resSig<-res[order(res$pvalue)[1:1000],1]   ##features
  
peak<-cbind(resistance_peak,sensitivity_peak)
file_name<-c("S8043","S7596","S7569","S7555","S2779","S2759","S2693","S2244","S2170","S1515",
             "S1194","S1122","S1096","S1095","S1090","S1085","S1047","S1030")
ic502<-ic50[rownames(ic50)%in%file_name,]
ic503<-apply(ic502,2,median)
AUC<-matrix(round(ic503[match(gsub("\\.","-",colnames(peak)),colnames(ic502))],2),nrow=1)
feature_peak<-peak[match(resSig,peak_name),]
feature_peak1<-t(feature_peak)
labelll<-as.data.frame(matrix(c(rep("resistance",ncol(resistance_peak)),rep("sensitivity",ncol(sensitivity_peak))),ncol=1))
expbbb<-cbind(feature_peak1,labelll)
colnames(expbbb)<-c(peak_name[match(resSig,peak_name)],"cancertype")

library(uwot)
set.seed(100)
#iris_umap <- uwot::umap(expbbb)
iris_umap <- uwot::umap(expbbb,n_neighbors = 10)
iris_sumap_res <- data.frame(iris_umap,cancertype=expbbb$cancertype,AUC=t(AUC))
colnames(iris_sumap_res)<-c("UMAP1","UMAP2","cancertype","AUC")
head(iris_sumap_res)
library(ggplot2)
p<-ggplot(iris_sumap_res,aes(UMAP1,UMAP2,color=AUC,shape=cancertype)) + 
    geom_point(size = 3) +theme_bw() + 
    scale_color_gradient(low="green", high="red")+
    theme(plot.title = element_text(hjust = 0.5)) + 
    labs(x="UMAP_1",y="UMAP_2",
         title = paste("A UMAP visualization of the ","HDAC",sep=""))+
    geom_text(aes(label = AUC),colour="black",size = 2.5,check_overlap = TRUE,nudge_y = (-0.1))
setwd("~/xjj/drug/drug_result/drug64_AUC/figure_feature/UMAP_Pmin_Top1000_DApeaks_Gradual_change")
ggsave(paste("HDAC","_Top1000_DApeaks.pdf",sep=""),p,width = 6, height = 5)
########  tSNE
set.seed(100)
library(Rtsne)
tsne_out <- Rtsne(as.matrix(expbbb[,-(ncol(expbbb))]),pca=FALSE,dims=2,
                  perplexity=10,theta=0.0) # Run TSNE
iris_umap <- data.frame(tsne_out$Y)
iris_sumap_res <- data.frame(tsne_out$Y,cancertype=expbbb$cancertype,AUC=t(AUC))
colnames(iris_sumap_res)<-c("tSNE1","tSNE2","cancertype","AUC")
head(iris_sumap_res)
library(ggplot2)
p<-ggplot(iris_sumap_res,aes(tSNE1,tSNE2,color=AUC,shape=cancertype)) + 
  geom_point(size = 3) +theme_bw() + 
  scale_color_gradient(low="green", high="red")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x="tSNE_1",y="tSNE_2",
       title = paste("A tSNE visualization of the ","HDAC",sep=""))+
  geom_text(aes(label = AUC),colour="black",size = 2.5,check_overlap = TRUE,nudge_y = (-3))
setwd("~/xjj/drug/drug_result/drug64_AUC/figure_feature/Rtsne_packages_Pmin_Top1000_DApeaks_Gradual_change")
ggsave(paste("HDAC","_Rtsne_packages_Pmin_Top1000_DApeaks.pdf",sep=""),p,width = 6, height = 5)

##########  PCA
set.seed(1000)
df_pca <- prcomp(as.matrix(expbbb[,-(ncol(expbbb))]),cor=T) # Run PCAprincomp(Inx,cor=T)
df_pcs <-data.frame(df_pca$x,cancertype=expbbb$cancertype)
iris_umap<- data.frame(df_pcs[,c(1,2)])
iris_sumap_res <- data.frame(df_pcs[,c(1,2,38)],AUC=t(AUC))
colnames(iris_sumap_res)<-c("PC1","PC2","cancertype","AUC")
head(iris_sumap_res)
library(ggplot2)
p<-ggplot(iris_sumap_res,aes(PC1,PC2,color=AUC,shape=cancertype)) + 
  geom_point(size = 3) +theme_bw() + 
  scale_color_gradient(low="green", high="red")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x="PCA_1",y="PCA_2",
       title = paste("A PCA visualization of the ","HDAC",sep=""))+
  geom_text(aes(label = AUC),colour="black",size = 2.5,check_overlap = TRUE,nudge_y = (-0.5))
  
setwd("~/xjj/drug/drug_result/drug64_AUC/figure_feature/PCA_Pmin_Top1000_DApeaks_Gradual_change")
ggsave(paste("HDAC","_PCA_Pmin_Top1000_DApeaks.pdf",sep=""),p,width = 6, height = 5)


####### MDS
set.seed(1000)
dis_iris = dist(expbbb[,-(ncol(expbbb))],p = 2)
mds_x = cmdscale(dis_iris)
mds_x = data.frame(mds_x)
iris_sumap_res <- data.frame(mds_x,cancertype=expbbb$cancertype,AUC=t(AUC))
colnames(iris_sumap_res)<-c("MDS1","MDS2","cancertype","AUC")
head(iris_sumap_res)
library(ggplot2)
p<-ggplot(iris_sumap_res,aes(MDS1,MDS2,color=AUC,shape=cancertype)) + 
  geom_point(size = 3) +theme_bw() + 
  scale_color_gradient(low="green", high="red")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x="MDS_1",y="MDS_2",
       title = paste("A MDS visualization of the ","HDAC",sep=""))+
  geom_text(aes(label = AUC),colour="black",size = 2,check_overlap = TRUE,nudge_y = (-0.5))

setwd("~/xjj/drug/drug_result/drug64_AUC/figure_feature/MDS_Pmin_Top1000_DApeaks_Gradual_change")
ggsave(paste("HDAC","_MDS_Pmin_Top1000_DApeaks.pdf",sep=""),p,width = 6, height = 5)


########  两点之间的欧氏距离
result<-NULL
for(j in 1:(ncol(feature_peak)-1)){
  a<-feature_peak[,j]
  x<-iris_umap[j,]
  result1<-NULL
  for(m in (j+1):ncol(feature_peak)){
    b<-feature_peak[,m]
    y<-iris_umap[m,]
    dat<-rbind(x,y)
    d<-distance(dat,method = "euclidean")
    r<-cor.test(a,b,method ="pearson",alternative="two.sided")$estimate
    result2<-matrix(c(r,d),nrow=1)
    result1<-rbind(result1,result2)
  }
  result<-rbind(result,result1)
}
dat1<-as.data.frame(result)
colnames(dat1)<-c("correlation","distance")
pp<-ggplot(data=dat1, aes(x=correlation, y=distance))+geom_point(color="red",size=0.1)+stat_smooth(method="lm",se=FALSE)+stat_cor(data=dat1, method = "pearson")+
  ggtitle(paste("HDAC","-pearson","-MDS",sep="")) +
  theme(plot.title = element_text(hjust = 0.5))
setwd("~/xjj/drug/drug_result/drug64_AUC/figure_feature/MDS_correlation_distance")
ggsave(paste("correlation_distance_","HDAC","_MDS.pdf",sep=""),pp,width = 5, height = 5)



###################   HDAC  +  chemotherapy   ######################
rm(list=ls())
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_count_name<-read.table("peak_RPKM.txt",sep="\t",header=T)
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
library(dplyr)
library("data.table")
library(foreach)
setwd("~/xjj/drug")
ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,-(which(colnames(ic501)==c("PC.100","PC.34")))]
colnames(ic50)<-sample_id[match(colnames(ic50),sample_id[,2]),1]
sensitivity<-c("PC.101","PC.109","PC.115","PC.116","PC.117","PC.13","PC.130","PC.136","PC.139","PC.16","PC.18","PC.2","PC.27","PC.64","PC.78","PC.8","PC.81","PC.97","PC.98","PC.L" )
resistance<-c("PC.102","PC.104","PC.105","PC.111","PC.112","PC.119","PC.121","PC.134","PC.135","PC.14","PC.22","PC.40","PC.5","PC.52","PC.56","PC.G","PC.I" )
resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]
resistance_peak<-peak_count_name[,match(resistance_new,gsub("\\.","-",colnames(peak_count_name)))]
sensitivity_peak<-peak_count_name[,match(sensitivity_new,gsub("\\.","-",colnames(peak_count_name)))]
peak_name<-peak_count_name[,1]

setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/2_ATACseq_DApeaks")
res<-fread("sensitivity-resistance-all-DESeq2_HDAC_chemotherapy_AUC_ATAC.txt",header=T,data.table=F)
resSig<-res[order(res$pvalue)[1:1000],1]   ##features

peak<-cbind(resistance_peak,sensitivity_peak)
file_name<-c("S8043","S7596","S7569","S7555","S2779","S2759","S2693","S2244","S2170","S1515",
             "S1194","S1122","S1096","S1095","S1090","S1085","S1047","S1030","5-FU","GEM","IRI","OXA","PAC")
ic502<-ic50[rownames(ic50)%in%file_name,]
ic503<-apply(ic502,2,median)
AUC<-matrix(round(ic503[match(gsub("\\.","-",colnames(peak)),colnames(ic502))],2),nrow=1)
feature_peak<-peak[match(resSig,peak_name),]
feature_peak1<-t(feature_peak)
labelll<-as.data.frame(matrix(c(rep("resistance",ncol(resistance_peak)),rep("sensitivity",ncol(sensitivity_peak))),ncol=1))
expbbb<-cbind(feature_peak1,labelll)
colnames(expbbb)<-c(peak_name[match(resSig,peak_name)],"cancertype")

library(uwot)
set.seed(100)
#iris_umap <- uwot::umap(expbbb)
iris_umap <- uwot::umap(expbbb,n_neighbors = 10)
iris_sumap_res <- data.frame(iris_umap,cancertype=expbbb$cancertype,AUC=t(AUC))
colnames(iris_sumap_res)<-c("UMAP1","UMAP2","cancertype","AUC")
head(iris_sumap_res)
library(ggplot2)
p<-ggplot(iris_sumap_res,aes(UMAP1,UMAP2,color=AUC,shape=cancertype)) + 
  geom_point(size = 3) +theme_bw() + 
  scale_color_gradient(low="green", high="red")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x="UMAP_1",y="UMAP_2",
       title = paste("A UMAP visualization of the ","HDAC + chemotherapy",sep=""))+
  geom_text(aes(label = AUC),colour="black",size = 2.5,check_overlap = TRUE,nudge_y = (-0.1))
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/5_dimensionality_reduction")
ggsave(paste("HDAC + chemotherapy","_Top1000_DApeaks_UMAP.pdf",sep=""),p,width = 6, height = 5)
########  tSNE
set.seed(100)
library(Rtsne)
tsne_out <- Rtsne(as.matrix(expbbb[,-(ncol(expbbb))]),pca=FALSE,dims=2,
                  perplexity=10,theta=0.0) # Run TSNE
iris_umap <- data.frame(tsne_out$Y)
iris_sumap_res <- data.frame(tsne_out$Y,cancertype=expbbb$cancertype,AUC=t(AUC))
colnames(iris_sumap_res)<-c("tSNE1","tSNE2","cancertype","AUC")
head(iris_sumap_res)
library(ggplot2)
p<-ggplot(iris_sumap_res,aes(tSNE1,tSNE2,color=AUC,shape=cancertype)) + 
  geom_point(size = 3) +theme_bw() + 
  scale_color_gradient(low="green", high="red")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x="tSNE_1",y="tSNE_2",
       title = paste("A tSNE visualization of the ","HDAC + chemotherapy",sep=""))+
  geom_text(aes(label = AUC),colour="black",size = 2.5,check_overlap = TRUE,nudge_y = (-3))
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/5_dimensionality_reduction")
ggsave(paste("HDAC + chemotherapy","_Rtsne_packages_Pmin_Top1000_DApeaks.pdf",sep=""),p,width = 6, height = 5)

##########  PCA
set.seed(1000)
df_pca <- prcomp(as.matrix(expbbb[,-(ncol(expbbb))]),cor=T) # Run PCAprincomp(Inx,cor=T)
df_pcs <-data.frame(df_pca$x,cancertype=expbbb$cancertype)
iris_umap<- data.frame(df_pcs[,c(1,2)])
iris_sumap_res <- data.frame(df_pcs[,c(1,2,38)],AUC=t(AUC))
colnames(iris_sumap_res)<-c("PC1","PC2","cancertype","AUC")
head(iris_sumap_res)
library(ggplot2)
p<-ggplot(iris_sumap_res,aes(PC1,PC2,color=AUC,shape=cancertype)) + 
  geom_point(size = 3) +theme_bw() + 
  scale_color_gradient(low="green", high="red")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x="PCA_1",y="PCA_2",
       title = paste("A PCA visualization of the ","HDAC + chemotherapy",sep=""))+
  geom_text(aes(label = AUC),colour="black",size = 2.5,check_overlap = TRUE,nudge_y = (-0.5))

setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/5_dimensionality_reduction")
ggsave(paste("HDAC + chemotherapy","_PCA_Pmin_Top1000_DApeaks.pdf",sep=""),p,width = 6, height = 5)


####### MDS
set.seed(1000)
dis_iris = dist(expbbb[,-(ncol(expbbb))],p = 2)
mds_x = cmdscale(dis_iris)
mds_x = data.frame(mds_x)
iris_sumap_res <- data.frame(mds_x,cancertype=expbbb$cancertype,AUC=t(AUC))
colnames(iris_sumap_res)<-c("MDS1","MDS2","cancertype","AUC")
head(iris_sumap_res)
library(ggplot2)
p<-ggplot(iris_sumap_res,aes(MDS1,MDS2,color=AUC,shape=cancertype)) + 
  geom_point(size = 3) +theme_bw() + 
  scale_color_gradient(low="green", high="red")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x="MDS_1",y="MDS_2",
       title = paste("A MDS visualization of the ","HDAC + chemotherapy",sep=""))+
  geom_text(aes(label = AUC),colour="black",size = 2,check_overlap = TRUE,nudge_y = (-0.5))

setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/5_dimensionality_reduction")
ggsave(paste("HDAC + chemotherapy","_MDS_Pmin_Top1000_DApeaks.pdf",sep=""),p,width = 6, height = 5)


########  两点之间的欧氏距离
result<-NULL
for(j in 1:(ncol(feature_peak)-1)){
  a<-feature_peak[,j]
  x<-iris_umap[j,]
  result1<-NULL
  for(m in (j+1):ncol(feature_peak)){
    b<-feature_peak[,m]
    y<-iris_umap[m,]
    dat<-rbind(x,y)
    d<-distance(dat,method = "euclidean")
    r<-cor.test(a,b,method ="pearson",alternative="two.sided")$estimate
    result2<-matrix(c(r,d),nrow=1)
    result1<-rbind(result1,result2)
  }
  result<-rbind(result,result1)
}
dat1<-as.data.frame(result)
colnames(dat1)<-c("correlation","distance")
pp<-ggplot(data=dat1, aes(x=correlation, y=distance))+geom_point(color="red",size=0.1)+stat_smooth(method="lm",se=FALSE)+stat_cor(data=dat1, method = "pearson")+
  ggtitle(paste("HDAC + chemotherapy","-pearson","-tSNE",sep="")) +
  theme(plot.title = element_text(hjust = 0.5))
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/5_dimensionality_reduction")
ggsave(paste("correlation_distance_","HDAC + chemotherapy","_tSNE.pdf",sep=""),pp,width = 5, height = 5)




############################################################
rm(list=ls())
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_count_name<-read.table("peak_RPKM.txt",sep="\t",header=T)
setwd("~/xjj/drug/drug_result/drug64_AUC")
tresult_want<-read.table("drug64_sample.txt",sep="\t",header=T)
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
library(dplyr)
library("data.table")
library(foreach)
setwd("~/xjj/drug")
ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,-(which(colnames(ic501)==c("PC.100","PC.34")))]
colnames(ic50)<-sample_id[match(colnames(ic50),sample_id[,2]),1]

foreach(k=1:nrow(tresult_want),.combine='rbind')%do%{
  file_name=as.character(tresult_want[k,1])
  cat(file_name,"\n")
  common_sample<-as.data.frame(tresult_want[match(file_name,tresult_want[,1]),])
  sensitivity<-unlist(strsplit(as.character(common_sample[1,5]), ","))
  resistance<-unlist(strsplit(as.character(common_sample[1,4]), ","))
  resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
  sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]
  resistance_peak<-peak_count_name[,match(resistance_new,gsub("\\.","-",colnames(peak_count_name)))]
  sensitivity_peak<-peak_count_name[,match(sensitivity_new,gsub("\\.","-",colnames(peak_count_name)))]
  peak_name<-peak_count_name[,1]
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/2_ATACseq_DApeaks",sep=""))
  res<-fread(paste("sensitivity-resistance-all-DESeq2_",file_name,"_ATAC.txt",sep=''),header=T,data.table=F)
  #res1<-res[abs(res$log2FoldChange)>1,]
  #cat(nrow(res1),"\n")
  resSig<-res[order(res$pvalue)[1:100],1]   ##features
  peak<-cbind(resistance_peak,sensitivity_peak)
  feature_peak<-peak[match(resSig,peak_name),]
  colnames(feature_peak)<-sample_id[match(gsub("\\.","-",colnames(feature_peak)),sample_id[,1]),2]
  feature_peak1<-cbind(peak_name[match(resSig,peak_name)],feature_peak)
  colnames(feature_peak1)<-c("Top 100 peaks",colnames(feature_peak))
  setwd("~/xjj/ATAC/Top 100 peaks")
  write.table(feature_peak1,file=paste(file_name,"_top_100_peaks.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')
  
  }


####################   HDAC + chemotherapy
rm(list=ls())
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_count_name<-read.table("peak_RPKM.txt",sep="\t",header=T)
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
library(dplyr)
library("data.table")
library(foreach)
setwd("~/xjj/drug")
ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,-(which(colnames(ic501)==c("PC.100","PC.34")))]
colnames(ic50)<-sample_id[match(colnames(ic50),sample_id[,2]),1]

  sensitivity<-c("PC.101","PC.109","PC.115","PC.116","PC.117","PC.13","PC.130","PC.136","PC.139","PC.16","PC.18","PC.2","PC.27","PC.64","PC.78","PC.8","PC.81","PC.97","PC.98","PC.L" )
  resistance<-c("PC.102","PC.104","PC.105","PC.111","PC.112","PC.119","PC.121","PC.134","PC.135","PC.14","PC.22","PC.40","PC.5","PC.52","PC.56","PC.G","PC.I" )
  resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
  sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]
  resistance_peak<-peak_count_name[,match(resistance_new,gsub("\\.","-",colnames(peak_count_name)))]
  sensitivity_peak<-peak_count_name[,match(sensitivity_new,gsub("\\.","-",colnames(peak_count_name)))]
  peak_name<-peak_count_name[,1]
  setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/2_ATACseq_DApeaks")
  res<-fread("sensitivity-resistance-all-DESeq2_HDAC_chemotherapy_AUC_ATAC.txt",header=T,data.table=F)
  #res1<-res[abs(res$log2FoldChange)>1,]
  #cat(nrow(res1),"\n")
  resSig<-res[order(res$pvalue)[1:100],1]   ##features
  peak<-cbind(resistance_peak,sensitivity_peak)
  feature_peak<-peak[match(resSig,peak_name),]
  colnames(feature_peak)<-sample_id[match(gsub("\\.","-",colnames(feature_peak)),sample_id[,1]),2]
  feature_peak1<-cbind(peak_name[match(resSig,peak_name)],feature_peak)
  colnames(feature_peak1)<-c("Top 100 peaks",colnames(feature_peak))
  setwd("~/xjj/ATAC/Top 100 peaks")
  write.table(feature_peak1,file=paste("HDAC_chemotherapy","_top_100_peaks.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')
  


  
 ##########  加验证集
  rm(list=ls())
  setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
  peak_count_name<-read.table("peak_RPKM.txt",sep="\t",header=T)
  setwd("~/xjj/drug")
  sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
  library(dplyr)
  library("data.table")
  library(foreach)
  setwd("~/xjj/drug")
  ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
  ic50<-ic501[,-(which(colnames(ic501)==c("PC.100","PC.34")))]
  colnames(ic50)<-sample_id[match(colnames(ic50),sample_id[,2]),1]
  
  va<-ic501[,which(colnames(ic501)==c("PC.100","PC.34"))]
  colnames(va)<-sample_id[match(colnames(va),sample_id[,2]),1]
  
  sensitivity<-c("PC.101","PC.109","PC.115","PC.116","PC.117","PC.13","PC.130","PC.136","PC.139","PC.16","PC.18","PC.2","PC.27","PC.64","PC.78","PC.8","PC.81","PC.97","PC.98","PC.L" )
  resistance<-c("PC.102","PC.104","PC.105","PC.111","PC.112","PC.119","PC.121","PC.134","PC.135","PC.14","PC.22","PC.40","PC.5","PC.52","PC.56","PC.G","PC.I" )
  resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
  sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]
  resistance_peak<-peak_count_name[,match(resistance_new,gsub("\\.","-",colnames(peak_count_name)))]
  sensitivity_peak<-peak_count_name[,match(sensitivity_new,gsub("\\.","-",colnames(peak_count_name)))]
  va_peak<-peak_count_name[,match(colnames(va),gsub("\\.","-",colnames(peak_count_name)))]
  peak_name<-peak_count_name[,1]
  
  setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/2_ATACseq_DApeaks")
  res<-fread("sensitivity-resistance-all-DESeq2_HDAC_chemotherapy_AUC_ATAC.txt",header=T,data.table=F)
  resSig<-res[order(res$pvalue)[1:1000],1]   ##features
  
  peak<-cbind(resistance_peak,sensitivity_peak)
  file_name<-c("S8043","S7596","S7569","S7555","S2779","S2759","S2693","S2244","S2170","S1515","S1848","S8495",
               "S1194","S1122","S1096","S1095","S1090","S1085","S1047","S1030","5-FU","GEM","IRI","OXA","PAC")
  ic502<-ic50[rownames(ic50)%in%file_name,]
  ic503<-apply(ic502,2,median)
  va1<-va[rownames(va)%in%file_name,]
  va2<-apply(va1,2,median)
  AUC1<-matrix(round(ic503[match(gsub("\\.","-",colnames(peak)),colnames(ic502))],2),nrow=1)
  AUC2<-matrix(round(va2[match(gsub("\\.","-",colnames(va_peak)),colnames(va1))],2),nrow=1)
  AUC<-cbind(AUC1,AUC2)
  feature_peak<-peak[match(resSig,peak_name),]
  vali_peak<-va_peak[match(resSig,peak_name),]
  feature_peak1<-t(cbind(feature_peak,vali_peak))
  labelll<-as.data.frame(matrix(c(rep("resistance",ncol(resistance_peak)),rep("sensitivity",ncol(sensitivity_peak)),rep("validation",2)),ncol=1))
  expbbb<-cbind(feature_peak1,labelll)
  colnames(expbbb)<-c(peak_name[match(resSig,peak_name)],"cancertype")
  
  library(uwot)
  set.seed(100)
  #iris_umap <- uwot::umap(expbbb)
  iris_umap <- uwot::umap(expbbb,n_neighbors = 10)
  iris_sumap_res <- data.frame(iris_umap,cancertype=expbbb$cancertype,AUC=t(AUC))
  colnames(iris_sumap_res)<-c("UMAP1","UMAP2","cancertype","AUC")
  head(iris_sumap_res)
  library(ggplot2)
  p<-ggplot(iris_sumap_res,aes(UMAP1,UMAP2,color=AUC,shape=cancertype)) + 
    geom_point(size = 3) +theme_bw() + 
    scale_color_gradient(low="green", high="red")+
    theme(plot.title = element_text(hjust = 0.5)) + 
    labs(x="UMAP_1",y="UMAP_2",
         title = paste("A UMAP visualization of the ","HDAC + chemotherapy",sep=""))+
    geom_text(aes(label = AUC),colour="black",size = 2.5,check_overlap = TRUE,nudge_y = (-0.1))
  setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/5_dimensionality_reduction")
  ggsave(paste("HDAC + chemotherapy","_Top1000_DApeaks_UMAP.pdf",sep=""),p,width = 6, height = 5)
  
  
  
  
  
  setwd("~/xjj/drug/drug_result/chemotherapy")
  feature_peak<-read.csv("covid_and_health_matrix.csv",header=F,sep="\t")
  feature_peak1<-feature_peak
  labelll<-as.data.frame(matrix(c(rep("resistance",ncol(resistance_peak)),rep("sensitivity",ncol(sensitivity_peak))),ncol=1))
  expbbb<-cbind(feature_peak1,labelll)
  colnames(expbbb)<-c(peak_name[match(resSig,peak_name)],"cancertype")
  library(uwot)
  set.seed(100)
  iris_umap <- uwot::umap(feature_peak1,n_neighbors =10,n_components =3) 
  write.table(iris_umap,file="3D-euclidean.txt",sep = '\t',col.names = F,row.names = F,quote = FALSE)
  
