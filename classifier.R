#########################  留一法构造 ########################################################
#SVM 独立训练 独立验证 leave one out
rm(list=ls())
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
setwd("~/xjj/drug")
ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,!(colnames(ic501)%in%c("PC.100","PC.34"))]## all of PDAC
colnames(ic50)<-sample_id[match(colnames(ic50),sample_id[,2]),1]
want_AUC_chemo<-ic50[rownames(ic50)%in%"GEM",]
median<-median(as.numeric(want_AUC_chemo))

sensitivity<-colnames(want_AUC_chemo)[which(as.numeric(want_AUC_chemo)>median)]
resistance<-colnames(want_AUC_chemo)[which(as.numeric(want_AUC_chemo)<=median)]

sensitivity_AUC<-want_AUC_chemo[,match(sensitivity,colnames(want_AUC_chemo))]
GEM_sensitive<-colnames(sensitivity_AUC[,order(sensitivity_AUC,decreasing =TRUE)][1:6])


want_AUC_HDAC<-ic50[rownames(ic50)%in%"S8495",]
GEM_sensitivity1<-colnames(want_AUC_chemo)[which(as.numeric(want_AUC_chemo)>median(as.numeric(want_AUC_chemo)))]
GEM_resistance1<-colnames(want_AUC_chemo)[which(as.numeric(want_AUC_chemo)<=median(as.numeric(want_AUC_chemo)))]
S8495_sensitivity1<-colnames(want_AUC_HDAC)[which(as.numeric(want_AUC_HDAC)>median(as.numeric(want_AUC_HDAC)))]
S8495_resistance1<-colnames(want_AUC_HDAC)[which(as.numeric(want_AUC_HDAC)<=median(as.numeric(want_AUC_HDAC)))]

GEM_res_S8495_sen1<-intersect(GEM_resistance1,S8495_sensitivity1)
GEM_res_S8495_res1<-intersect(GEM_resistance1,S8495_resistance1)

GEM_res_S8495_sen_AUC<-want_AUC_HDAC[,match(GEM_res_S8495_sen1,colnames(want_AUC_HDAC))]
GEM_res_S8495_sen<-colnames(GEM_res_S8495_sen_AUC[,order(GEM_res_S8495_sen_AUC,decreasing =TRUE)][1:6])

GEM_res_S8495_res_AUC<-want_AUC_HDAC[,match(GEM_res_S8495_res1,colnames(want_AUC_HDAC))]
GEM_res_S8495_res<-colnames(GEM_res_S8495_res_AUC[,order(GEM_res_S8495_res_AUC,decreasing =TRUE)][1:6])



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

GEM_sensitive<-c("DAC-1","DAC-18","DAC-38","DAC-24","DAC-19","DAC-6" )
GEM_res_S8495_sen<-c("DAC-14","DAC-7","DAC-8","DAC-12","DAC-31","DAC-26")
GEM_res_S8495_res<-c("DAC-23","DAC-20","DAC-21","DAC-22","DAC-25","DAC-39")
all<-c(GEM_sensitive,GEM_res_S8495_sen,GEM_res_S8495_res)

setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
library("data.table")
exp11<-fread("star_rsem.GeneSymbol.readscount.xls",header=T,data.table=F)
exp1<-exp11[,match(all,colnames(exp11))]
exp2<-floor(exp1)
delete<-apply(exp2,1,function(x) mean(x==0))
cou<-which(delete>0.75)####在超过75%样本以上的都是0的基因删去
peak_count_name<-exp2[(-cou),]
geneid<-exp11[(-cou),1]



text_leave_one_out1<-matrix(0,18,2)
for(i in 1:18){
  cat(i,"\n")
  GEM_sensitive<-c("DAC-1","DAC-18","DAC-38","DAC-24","DAC-19","DAC-6" )
  GEM_res_S8495_sen<-c("DAC-14","DAC-7","DAC-8","DAC-12","DAC-31","DAC-26")
  GEM_res_S8495_res<-c("DAC-23","DAC-20","DAC-21","DAC-22","DAC-25","DAC-39")
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
  
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier_RNAseq_leave_one18/classifier',i,sep='_')
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
  
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier_RNAseq_leave_one18/classifier',i,sep='_')
  if (!file.exists(outputPath)){
    dir.create(outputPath, recursive=TRUE)
  }
  write.table(res1,paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_S8495_AUC_RNA.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
}


###############  ATACseq
rm(list=ls())
library(DESeq2)
library("data.table")
library(dplyr)
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_count_name<-read.table("peak_count.txt",sep="\t",header=T) #没有全0的
colnames(peak_count_name)<-gsub("\\.","-",colnames(peak_count_name))

GEM_sensitive<-c("DAC-1","DAC-18","DAC-38","DAC-24","DAC-19","DAC-6" )
GEM_res_S8495_sen<-c("DAC-14","DAC-7","DAC-8","DAC-12","DAC-31","DAC-26")
GEM_res_S8495_res<-c("DAC-23","DAC-20","DAC-21","DAC-22","DAC-25","DAC-39")
all<-c(GEM_sensitive,GEM_res_S8495_sen,GEM_res_S8495_res)


geneid<-peak_count_name[,1]



text_leave_one_out1<-matrix(0,18,2)
for(i in 1:18){
  cat(i,"\n")
  GEM_sensitive<-c("DAC-1","DAC-18","DAC-38","DAC-24","DAC-19","DAC-6" )
  GEM_res_S8495_sen<-c("DAC-14","DAC-7","DAC-8","DAC-12","DAC-31","DAC-26")
  GEM_res_S8495_res<-c("DAC-23","DAC-20","DAC-21","DAC-22","DAC-25","DAC-39")
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
  res<-cbind(peak_name,res)  ##对数据增加一列
  colnames(res)<- c('peak_id',"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj" )##对列名进行修改
  
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier_RNAseq_leave_one18/classifier',i,sep='_')
  if (!file.exists(outputPath)){
    dir.create(outputPath, recursive=TRUE)
  }
  write.table(res,paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_GEM_AUC_ATAC.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  write.table(text,paste(outputPath,'/',"textATAC.txt",sep=""),sep = '\t',col.names = F,row.names = F,quote = FALSE)##数据输出
  
  sensitivity1<-GEM_res_S8495_sen
  resistance1<-GEM_res_S8495_res
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
  res1<-cbind(peak_name,res1)  ##对数据增加一列
  colnames(res1)<- c('peak_id',"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj" )##对列名进行修改
  
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier_RNAseq_leave_one18/classifier',i,sep='_')
  if (!file.exists(outputPath)){
    dir.create(outputPath, recursive=TRUE)
  }
  write.table(res1,paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_S8495_AUC_ATAC.txt",sep=""),sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
}

############################################################################
#######################  留一法验证 RNAseq ############################
rm(list=ls())
setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
peak_RPKM<-read.table("star_rsem.GeneSymbol.FPKM.xls",sep="\t",header=T) #没有全0的
colnames(peak_RPKM)<-gsub("\\.","-",colnames(peak_RPKM))

GEM_sensitive<-c("DAC-1","DAC-18","DAC-38","DAC-24","DAC-19","DAC-6" )
GEM_res_S8495_sen<-c("DAC-14","DAC-7","DAC-8","DAC-12","DAC-31","DAC-26")
GEM_res_S8495_res<-c("DAC-23","DAC-20","DAC-21","DAC-22","DAC-25","DAC-39")

svm_result<-NULL
for(i in 1:18){
  cat(i,"\n")
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier_RNAseq_leave_one18/classifier',i,sep='_')
  GEM<-read.table(paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_GEM_AUC_RNA.txt",sep=""),sep = '\t',header=T)##数据输出
  resSig<-GEM[which(GEM$padj<0.05 & abs(GEM$log2FoldChange)>2),]
  #resSig<-GEM[which(GEM$pvalue<0.01),]
  resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
  resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
  cat(sum(resSig$up_down=='up'),"\n")
  GEM_sensitive_DARs<-as.character(resSig[which(resSig$log2FoldChange>0),1])
  
  S8495<-read.table(paste(outputPath,'/',"sensitivity-resistance-all-DESeq2_S8495_AUC_RNA.txt",sep=""),sep = '\t',header=T)##数据输出
  resSig_S8495<-S8495[which(S8495$padj<0.05 & abs(S8495$log2FoldChange)>2),]
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
  y<-matrix(c(rep("GEM_sensitive",6),rep("GEM_res_S8495_sen",6),rep("GEM_res_S8495_res",6)),ncol=1)
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
  #SvmFit<-svm(y~.,data=train_data,type="C-classification",kernel="radial",cost=1,gamma=1,scale=FALSE)
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
setwd("~/xjj/drug/drug_result/HDACi_chemo618/classifier_RNAseq_leave_one18")
write.table(svm_result,"KNN3_result_leave-one-out.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
stopCluster(cl)

#######################  留一法验证 RNAseq 前100个基因############################
rm(list=ls())
setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
peak_RPKM<-read.table("star_rsem.GeneSymbol.FPKM.xls",sep="\t",header=T) #没有全0的
colnames(peak_RPKM)<-gsub("\\.","-",colnames(peak_RPKM))

GEM_sensitive<-c("DAC-1","DAC-18","DAC-38","DAC-24","DAC-19","DAC-6" )
GEM_res_S8495_sen<-c("DAC-14","DAC-7","DAC-8","DAC-12","DAC-31","DAC-26")
GEM_res_S8495_res<-c("DAC-23","DAC-20","DAC-21","DAC-22","DAC-25","DAC-39")


svm_result<-NULL
for(i in 1:18){
  cat(i,"\n")
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier_RNAseq_leave_one18/classifier',i,sep='_')
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
  y<-matrix(c(rep("GEM_sensitive",6),rep("GEM_res_S8495_sen",6),rep("GEM_res_S8495_res",6)),ncol=1)
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
setwd("~/xjj/drug/drug_result/HDACi_chemo618/classifier_RNAseq_leave_one18")
write.table(svm_result,"NaiveBayes_result_leave-one-out100.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
stopCluster(cl)


#######################  留一法验证 ############################
rm(list=ls())
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_RPKM<-read.table("peak_RPKM.txt",sep="\t",header=T) #没有全0的
colnames(peak_RPKM)<-gsub("\\.","-",colnames(peak_RPKM))

GEM_sensitive<-c("DAC-1","DAC-18","DAC-38","DAC-24","DAC-19","DAC-6" )
GEM_res_S8495_sen<-c("DAC-14","DAC-7","DAC-8","DAC-12","DAC-31","DAC-26")
GEM_res_S8495_res<-c("DAC-23","DAC-20","DAC-21","DAC-22","DAC-25","DAC-39")


svm_result<-NULL
library(doParallel)
library(parallel)
library(iterators)

cl <- makeCluster(18)
registerDoParallel(cl)
foreach(i=1:18,.combine='rbind')%dopar%{
#for(i in 1:37){
  cat(i,"\n")
  outputPath = paste('~/xjj/drug/drug_result/HDACi_chemo618/classifier_RNAseq_leave_one18/classifier',i,sep='_')
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
  y<-matrix(c(rep("GEM_sensitive",6),rep("GEM_res_S8495_sen",6),rep("GEM_res_S8495_res",6)),ncol=1)
  colnames(y)<-"y"
  data<-as.data.frame(cbind(peak_RPKM37_DARs1,y))
  text<-read.table(paste(outputPath,'/',"textATAC.txt",sep=""),sep = '\t',header=F)##数据输出
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
setwd("~/xjj/drug/drug_result/HDACi_chemo618/classifier_RNAseq_leave_one18")
write.table(svm_result,"KNN_result_leave-one-out.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
stopCluster(cl)


