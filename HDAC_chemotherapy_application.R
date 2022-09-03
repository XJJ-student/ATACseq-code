########  Classifier
###  SVM
rm(list=ls())
library(e1071)
setwd("~/xjj/drug/drug_result/drug64_AUC")
tresult_want<-read.table("drug64_sample.txt",sep="\t",header=T)


result_final<-NULL
for(i in 1:5){
  setwd(paste("~/xjj/drug/drug_result/HDAC20_chemotherapy5/6_application/data/data/fold-second",i,"-500",sep=""))
  cat(i,"\n")
  result_each<-matrix(0,64,3)
  foreach(k=1:nrow(tresult_want),.combine='rbind')%do%{
    file_name=as.character(tresult_want[k,1])
    cat(file_name,"\n")
    res_train<-read.table(paste(file_name,"-resistance-train-",i,"-500.txt",sep=""),sep="\t",header=T)
    sen_train<-read.table(paste(file_name,"-sensitivity-train-",i,"-500.txt",sep=""),sep="\t",header=T)
    train_x<-apply(cbind(res_train[1:50,],sen_train[1:50,]),1,as.numeric)
    train_y<-c(rep("resistance",ncol(res_train)),rep("sensitivity",ncol(sen_train)))
    data_train<-data.frame(train_x,train_y) 
    res_val<-read.table(paste(file_name,"-resistance-validation-",i,"-500.txt",sep=""),sep="\t",header=T)
    sen_val<-read.table(paste(file_name,"-sensitivity-validation-",i,"-500.txt",sep=""),sep="\t",header=T)
    val_x<-apply(cbind(res_val[1:50,],sen_val[1:50,]),1,as.numeric)
    val_y<-c(rep("resistance",ncol(res_val)),rep("sensitivity",ncol(sen_val)))
    data_val<-data.frame(val_x,val_y)
    #model1 <- svm(train_y ~ ., data = data_train) 
    #summary(model1)
    #model3 = svm(as.data.frame(train_x),as.factor(train_y),kernel="radial")
    #pred<-predict(model3,train_x)
    #table(pred,train_y)
    SvmFit<-svm(train_y~.,data=train_x,type="C-classification",kernel="linear",cost=0.1,scale=FALSE)
    summary(SvmFit)
    set.seed(12345)
    tObj<-tune.svm(train_y~.,data = data_train,type="C-classification",kernel="linear",cost=c(0.001,0.01,0.1,1,5,10,100,1000),scale=FALSE)
    summary(tObj)
    BestSvm<-tObj$best.model
    summary(BestSvm)
    yPred<-predict(BestSvm,data_val)
    ConfM<-table(yPred,val_y)
    #预测误差率
    Err<-(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
    result_each[k,1]<-file_name
    result_each[k,2]<-i
    result_each[k,3]<-1-Err
  }
  result_final<-cbind(result_final,result_each)
}
colnames(result_final)<-rep(c("Drug","cross-validation","correct"),times=5)
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/6_application/data/data")
write.table(result_final,"SVM_second_50.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)
data<-as.data.frame(rbind(result_final[,2:3],result_final[,5:6],result_final[,8:9],result_final[,11:12],result_final[,14:15]))
colnames(data)<-c("group","value")
library(ggplot2)
p<-ggplot(as.data.frame(data),aes(x=as.factor(group),y=as.numeric(as.character(value)),fill=group))+
  #geom_boxplot(width=0.5)+
  geom_violin()+
  geom_jitter(aes(color=group),width=.2,size=0.5)+
  ylab('Correct rate')+#修改y轴名称
  xlab('Cross-validation')+
  labs(title = "SVM_second_50")+
  theme(plot.title = element_text(hjust = 0.5))
plot(p) 
ggsave("SVM_second_50.pdf",p,width = 10, height = 6)



###  randomForest
library(randomForest)
result_rf<-NULL
for(i in 1:5){
  setwd(paste("~/xjj/drug/drug_result/HDAC20_chemotherapy5/6_application/data/data/fold-second",i,"-500",sep=""))
  cat(i,"\n")
  result_each<-matrix(0,64,3)
  foreach(k=1:nrow(tresult_want),.combine='rbind')%do%{
    file_name=as.character(tresult_want[k,1])
    cat(file_name,"\n")
    res_train<-read.table(paste(file_name,"-resistance-train-",i,"-500.txt",sep=""),sep="\t",header=T)
    sen_train<-read.table(paste(file_name,"-sensitivity-train-",i,"-500.txt",sep=""),sep="\t",header=T)
    train_x<-apply(cbind(res_train,sen_train),1,as.numeric)
    train_y<-c(rep("resistance",ncol(res_train)),rep("sensitivity",ncol(sen_train)))
    data_train<-data.frame(train_x,train_y) 
    res_val<-read.table(paste(file_name,"-resistance-validation-",i,"-500.txt",sep=""),sep="\t",header=T)
    sen_val<-read.table(paste(file_name,"-sensitivity-validation-",i,"-500.txt",sep=""),sep="\t",header=T)
    val_x<-apply(cbind(res_val,sen_val),1,as.numeric)
    val_y<-c(rep("resistance",ncol(res_val)),rep("sensitivity",ncol(sen_val)))
    data_val<-data.frame(val_x,val_y)
    set.seed(12345)
    rf.train<-randomForest(as.factor(train_y)~.,data = data_train,importance=TRUE,na.action = na.pass)
    #需要的是类别，所以type="class"
    rf.test<-predict(rf.train,newdata=data_val,type="class")
    #ROC需要一个概率，所以type=“prob”
    rf.test2<-predict(rf.train,newdata=data_val,type="prob")
    ConfM<-table(rf.test,val_y)
    #预测误差率
    Err<-(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
    result_each[k,1]<-file_name
    result_each[k,2]<-i
    result_each[k,3]<-1-Err
  }
  result_rf<-cbind(result_rf,result_each)
}

colnames(result_rf)<-rep(c("Drug","cross-validation","correct"),times=5)
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/6_application/data/data")
write.table(result_rf,"randomForest_second_500.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)
data<-as.data.frame(rbind(result_rf[,2:3],result_rf[,5:6],result_rf[,8:9],result_rf[,11:12],result_rf[,14:15]))
colnames(data)<-c("group","value")
library(ggplot2)
p<-ggplot(as.data.frame(data),aes(x=as.factor(group),y=as.numeric(as.character(value)),fill=group))+
  #geom_boxplot(width=0.5)+
  geom_violin()+
  geom_jitter(aes(color=group),width=.2,size=0.5)+
  ylab('Correct rate')+#修改y轴名称
  xlab('Cross-validation')+
  labs(title = "randomForest_second_500")+
  theme(plot.title = element_text(hjust = 0.5))
plot(p) 
ggsave("randomForest_second_500.pdf",p,width = 10, height = 6)



### KNN
library(class)
result_knn<-NULL
for(i in 1:5){
  setwd(paste("~/xjj/drug/drug_result/HDAC20_chemotherapy5/6_application/data/data/fold-second",i,"-500",sep=""))
  cat(i,"\n")
  result_each<-matrix(0,64,3)
  foreach(k=1:nrow(tresult_want),.combine='rbind')%do%{
    file_name=as.character(tresult_want[k,1])
    cat(file_name,"\n")
    res_train<-read.table(paste(file_name,"-resistance-train-",i,"-500.txt",sep=""),sep="\t",header=T)
    sen_train<-read.table(paste(file_name,"-sensitivity-train-",i,"-500.txt",sep=""),sep="\t",header=T)
    train_x<-apply(cbind(res_train,sen_train),1,as.numeric)
    train_y<-c(rep("resistance",ncol(res_train)),rep("sensitivity",ncol(sen_train)))
    data_train<-data.frame(train_x,train_y) 
    res_val<-read.table(paste(file_name,"-resistance-validation-",i,"-500.txt",sep=""),sep="\t",header=T)
    sen_val<-read.table(paste(file_name,"-sensitivity-validation-",i,"-500.txt",sep=""),sep="\t",header=T)
    val_x<-apply(cbind(res_val,sen_val),1,as.numeric)
    val_y<-c(rep("resistance",ncol(res_val)),rep("sensitivity",ncol(sen_val)))
    data_val<-data.frame(val_x,val_y)
    set.seed(12345)
    cl <- factor(train_y)
    knn.test<-knn(train_x, val_x, cl, k = 3, prob = FALSE)
    ConfM<-table(knn.test,val_y)
    #预测误差率
    Err<-(sum(ConfM)-sum(diag(ConfM)))/sum(ConfM)
    result_each[k,1]<-file_name
    result_each[k,2]<-i
    result_each[k,3]<-1-Err
  }
  result_knn<-cbind(result_knn,result_each)
}

colnames(result_knn)<-rep(c("Drug","cross-validation","correct"),times=5)
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/6_application/data/data")
write.table(result_knn,"knn_second_500.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)
data<-as.data.frame(rbind(result_rf[,2:3],result_rf[,5:6],result_rf[,8:9],result_rf[,11:12],result_rf[,14:15]))
colnames(data)<-c("group","value")
library(ggplot2)
p<-ggplot(as.data.frame(data),aes(x=as.factor(group),y=as.numeric(as.character(value)),fill=group))+
  #geom_boxplot(width=0.5)+
  geom_violin()+
  geom_jitter(aes(color=group),width=.2,size=0.5)+
  ylab('Correct rate')+#修改y轴名称
  xlab('Cross-validation')+
  labs(title = "knn_second_500")+
  theme(plot.title = element_text(hjust = 0.5))
plot(p) 
ggsave("knn_second_500.pdf",p,width = 10, height = 6)







############### clinical PFS  ################
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/7_clinical")
clinical<-read.table("clinical.txt",header = TRUE,sep = "\t")
setwd("~/xjj/drug/drug_result/drug64_AUC")
tresult_want<-read.table("drug64_sample.txt",sep="\t",header=T)
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
setwd("~/xjj/drug")
ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,-(which(colnames(ic501)==c("PC.100","PC.34")))]
colnames(ic50)<-sample_id[match(colnames(ic50),sample_id[,2]),1]
library(dplyr)
library("data.table")
library(foreach)

### all
clinical_want<-clinical[match(intersect(clinical[,1],colnames(ic50)),clinical[,1]),]
ic50_want<-ic50[,match(intersect(clinical[,1],colnames(ic50)),colnames(ic50))]

corralation_all<-matrix(0,64,7)
for(k in 1:nrow(ic50_want)){
  file_name=as.character(rownames(ic50_want)[k])
  cat(file_name,"\n")
  corralation_all[k,1]<-file_name
  corralation_all[k,2]<-"spearman"
  cor_result_spearman<-cor.test(as.numeric(ic50_want[k,]),clinical_want[,3],alternative = "two.sided",method = "spearman")
  corralation_all[k,3]<-cor_result_spearman$estimate
  corralation_all[k,4]<-cor_result_spearman$p.value
  corralation_all[k,5]<-"pearson"
  cor_result_pearson<-cor.test(as.numeric(ic50_want[k,]),clinical_want[,3],alternative = "two.sided",method = "pearson")
  corralation_all[k,6]<-cor_result_pearson$estimate
  corralation_all[k,7]<-cor_result_pearson$p.value
}
colnames(corralation_all)<-c("Drug","spearman","spearman-estimate","spearman-pvalue","pearson","pearson-estimate","pearson-pvalue")
corralation_all[which(corralation_all[,4]<0.05),] ##only "S7958"

### sen vs res ##########
corralation_all<-matrix(0,64,13)
for(k in 1:nrow(tresult_want)){
  file_name=as.character(tresult_want[k,1])
  cat(file_name,"\n")
  common_sample<-as.data.frame(tresult_want[match(file_name,tresult_want[,1]),])
  sensitivity<-unlist(strsplit(as.character(common_sample[1,5]), ","))
  resistance<-unlist(strsplit(as.character(common_sample[1,4]), ","))
  resistance_new<-as.character(sample_id[match(resistance,sample_id[,2]),1])
  sensitivity_new<-as.character(sample_id[match(sensitivity,sample_id[,2]),1])
  
  corralation_all[k,1]<-file_name
  clinical_want_res<-clinical[match(intersect(clinical[,1],resistance_new),clinical[,1]),]
  ic50_want_res<-ic50[,match(intersect(clinical[,1],resistance_new),colnames(ic50))]
  corralation_all[k,2]<-"spearman"
  cor_result_spearman<-cor.test(as.numeric(ic50_want_res[k,]),clinical_want_res[,3],alternative = "two.sided",method = "spearman")
  corralation_all[k,3]<-cor_result_spearman$estimate
  corralation_all[k,4]<-cor_result_spearman$p.value
  corralation_all[k,5]<-"pearson"
  cor_result_pearson<-cor.test(as.numeric(ic50_want_res[k,]),clinical_want_res[,3],alternative = "two.sided",method = "pearson")
  corralation_all[k,6]<-cor_result_pearson$estimate
  corralation_all[k,7]<-cor_result_pearson$p.value
  
  clinical_want_sen<-clinical[match(intersect(clinical[,1],sensitivity_new),clinical[,1]),]
  ic50_want_sen<-ic50[,match(intersect(clinical[,1],sensitivity_new),colnames(ic50))]
  corralation_all[k,8]<-"spearman"
  cor_result_spearman_sen<-cor.test(as.numeric(ic50_want_sen[k,]),clinical_want_sen[,3],alternative = "two.sided",method = "spearman")
  corralation_all[k,9]<-cor_result_spearman_sen$estimate
  corralation_all[k,10]<-cor_result_spearman_sen$p.value
  corralation_all[k,11]<-"pearson"
  cor_result_pearson_sen<-cor.test(as.numeric(ic50_want_sen[k,]),clinical_want_sen[,3],alternative = "two.sided",method = "pearson")
  corralation_all[k,12]<-cor_result_pearson_sen$estimate
  corralation_all[k,13]<-cor_result_pearson_sen$p.value
  
}
colnames(corralation_all)<-c("Drug","spearman-res","spearman-estimate-res","spearman-pvalue-res","pearson-res","pearson-estimate-res","pearson-pvalue-res","spearman-sen","spearman-estimate-sen","spearman-pvalue-sen","pearson-sen","pearson-estimate-sen","pearson-pvalue-sen")
corralation_all[which(corralation_all[,4]<0.05),1] ## res only "S2718" "S7596"
corralation_all[which(corralation_all[,7]<0.05),1] ## res only "S7575"
corralation_all[which(corralation_all[,10]<0.05),1] ##sen only "S1047" "S2244"
corralation_all[which(corralation_all[,13]<0.05),1] ##sen only "S7958"

##################################################################################3333
################ 37 patient   
rm(list=ls())
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
setwd("~/xjj/drug")
ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
ic5011<-ic501[,!(colnames(ic501)%in%c("PC.100","PC.34"))]
colnames(ic5011)<-sample_id[match(colnames(ic5011),sample_id[,2]),1]
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/7_clinical")
NC_AUC<-read.table("NC-AUC.txt",header = TRUE,sep = "\t")
colnames(NC_AUC)<-gsub("\\.","-",gsub("CAS.","",colnames(NC_AUC)))
ic50<-NC_AUC[,match(colnames(ic5011),colnames(NC_AUC))]
rownames(ic50)<-NC_AUC[,1]
rownames(ic50)[5]<-"PAC"
setwd("~/xjj/drug")
drug_info<-read.csv("compiled_drug_info.csv",header = TRUE,sep = ",")
drug_result<-drug_info[match(rownames(ic50),drug_info$drug_id),]
HDAC1<-drug_result[drug_result$target=="HDAC",]
ccc<-c("S1848","S2759","S1194","S1047")
HDAC2<-drug_result[drug_result$drug_id%in%ccc,]
HDAC3<-rbind(HDAC1,HDAC2)
HDAC<-HDAC3[order(HDAC3[,1],decreasing = T),]
HDAC_ic50<-ic50[match(as.character(HDAC[,1]),rownames(ic50)),]
HDAC_ic50<-HDAC_ic50[-which(rownames(HDAC_ic50)=="S7555"|rownames(HDAC_ic50)=="S1122"|rownames(HDAC_ic50)=="S1848"|rownames(HDAC_ic50)=="S8495"),]

HDAC_ic50<-HDAC_ic50[-which(rownames(HDAC_ic50)=="S7555"|rownames(HDAC_ic50)=="S1122"|rownames(HDAC_ic50)=="S1848"|rownames(HDAC_ic50)=="S8495"|rownames(HDAC_ic50)=="S7569"|rownames(HDAC_ic50)=="S2779"),]

ccc<-c("S1096","S1515","S1095","S1030","S2170","S1090")

library(pheatmap)
a=cor(t(HDAC_ic50))
result=pheatmap(a,scale = "none",main = "The correlation coefficients between HDACi and chemotherapy' AUC",show_rownames=T,show_colnames=T,
                clustering_distance_rows = "correlation",
                clustering_distance_cols = "correlation",clustering_method = "complete",cutree_rows=2,cutree_cols=2)

##########   20种HDAC+5种化疗药
rm(list=ls())
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
setwd("~/xjj/drug")
ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,!(colnames(ic501)%in%c("PC.100","PC.34"))]
colnames(ic50)<-sample_id[match(colnames(ic50),sample_id[,2]),1]

drug_info<-read.csv("compiled_drug_info.csv",header = TRUE,sep = ",")
drug_result<-drug_info[match(rownames(ic50),drug_info$drug_id),]
HDAC1<-drug_result[drug_result$target=="HDAC",]
ccc<-c("S1848","S2759","S1194","S1047")
ccc<-c("5-FU","GEM","IRI","OXA","PAC")
HDAC2<-drug_result[drug_result$drug_id%in%ccc,]
HDAC3<-rbind(HDAC1,HDAC2)
HDAC<-HDAC3[order(HDAC3[,1],decreasing = T),]
HDAC_ic50<-ic50[match(as.character(HDAC[,1]),rownames(ic50)),]

bb<-t(apply(HDAC_ic50,1,function(x) scale(x)))
colnames(bb)<-colnames(HDAC_ic50)
norma_result<-pheatmap(bb,scale = "none",main = "The AUC of HDACi and chemotherapy are normalized",show_rownames=T,show_colnames=T,
                       clustering_distance_rows = "correlation",
                       clustering_distance_cols = "euclidean",clustering_method = "ward.D2")


norma_result<-pheatmap(HDAC_ic50,scale = "none",main = "The AUC of HDACi and chemotherapy are normalized",show_rownames=T,show_colnames=T,
                       clustering_distance_rows = "correlation",
                       clustering_distance_cols = "euclidean",clustering_method = "ward.D2")


clusters_filter<-cutree(norma_result$tree_col,k=2)
sensitivity=colnames(HDAC_ic50[,clusters_filter==1])
resistance=colnames(HDAC_ic50[,clusters_filter==2])
#resistance=c("PC.2","PC.22","PC.13","PC.23")
#medianl<-resistance[-match(resistance,resistance1)]


library(org.Hs.eg.db)
library(stringr)
library(clusterProfiler)
library(survival)
library(survminer)

setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/7_clinical")
clinical<-read.table("clinical.txt",header = TRUE,sep = "\t")
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
library("data.table")
resistance_cli<-clinical[match(intersect(clinical[,1],resistance),clinical[,1]),]
sensitivity_cli<-clinical[match(intersect(clinical[,1],sensitivity),clinical[,1]),]
all<-rbind(resistance_cli,sensitivity_cli)
label<-as.matrix(c(rep("resistance",nrow(resistance_cli)),rep("sensitivity",nrow(sensitivity_cli))),ncol=1)
TTP=as.numeric(all[,3])  ##生存时间
status_TTP=as.matrix(all[,2])  ##生存状态
status_TTP[which(status_TTP=="Disease Free")]=0
status_TTP[which(status_TTP=="Recurred")]=1
status_TTP<-as.numeric(status_TTP)
surv_info1<-as.data.frame(cbind(TTP,status_TTP))

surv_TTP<-survfit(Surv(TTP, status_TTP) ~ label,data=surv_info1)
ggsurvplot(surv_TTP,
           pval = TRUE, conf.int = TRUE, #是否显示p值/显示风险区域
           risk.table = TRUE,# 将风险表显示在生存曲线下面
           ggtheme = theme_bw(),
           tables.theme = theme_cleantable(),
           title="HDACIs combined with Chemotherapeutic drugs")

######
clusters_filter<-cutree(norma_result$tree_col,k=2)
sensitivity=colnames(HDAC_ic50[,clusters_filter==1])
resistance1=colnames(HDAC_ic50[,clusters_filter==2])
resistance=c("DAC-22","DAC-13","DAC-23")
medianl<-resistance1[-match(resistance,resistance1)]

setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/7_clinical")
clinical<-read.table("clinical.txt",header = TRUE,sep = "\t")

resistance_cli<-clinical[match(intersect(clinical[,1],resistance),clinical[,1]),]
sensitivity_cli<-clinical[match(intersect(clinical[,1],sensitivity),clinical[,1]),]
medianl_cli<-clinical[match(intersect(clinical[,1],medianl),clinical[,1]),]

all<-rbind(resistance_cli,sensitivity_cli,medianl_cli)
label<-as.matrix(c(rep("resistance",nrow(resistance_cli)),rep("sensitivity",nrow(sensitivity_cli)),rep("medianl",nrow(medianl_cli))),ncol=1)
TTP=as.numeric(all[,3])  ##生存时间
status_TTP=as.matrix(all[,2])  ##生存状态
status_TTP[which(status_TTP=="Disease Free")]=0
status_TTP[which(status_TTP=="Recurred")]=1
status_TTP<-as.numeric(status_TTP)
surv_info1<-as.data.frame(cbind(TTP,status_TTP))

surv_TTP<-survfit(Surv(TTP, status_TTP) ~ label,data=surv_info1)
ggsurvplot(surv_TTP,
           pval = TRUE, conf.int = TRUE, #是否显示p值/显示风险区域
           risk.table = TRUE,# 将风险表显示在生存曲线下面
           ggtheme = theme_bw(),
           tables.theme = theme_cleantable(),
           title="HDACIs combined with Chemotherapeutic drugs")




########  31 patient
clusters_filter<-cutree(norma_result$tree_col,k=2)
sensitivity=colnames(HDAC_ic50[,clusters_filter==1])
resistance=colnames(HDAC_ic50[,clusters_filter==2])

setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/7_clinical")
clinical<-read.table("31patient_clinical.txt",header = TRUE,sep = "\t")
clinical[,1]<-gsub("CAS-","",clinical[,1])
resistance_cli<-clinical[match(intersect(clinical[,1],resistance),clinical[,1]),]
sensitivity_cli<-clinical[match(intersect(clinical[,1],sensitivity),clinical[,1]),]
all<-rbind(resistance_cli,sensitivity_cli)
label<-as.matrix(c(rep("resistance",nrow(resistance_cli)),rep("sensitivity",nrow(sensitivity_cli))),ncol=1)
TTP=as.numeric(all[,4])  ##生存时间
status_TTP=as.matrix(all[,3])  ##生存状态
status_TTP[which(status_TTP=="Disease Free")]=0
status_TTP[which(status_TTP=="Recurred")]=1
status_TTP<-as.numeric(status_TTP)
surv_info1<-as.data.frame(cbind(TTP,status_TTP))

surv_TTP<-survfit(Surv(TTP, status_TTP) ~ label,data=surv_info1)
ggsurvplot(surv_TTP,
           pval = TRUE, conf.int = TRUE, #是否显示p值/显示风险区域
           risk.table = TRUE,# 将风险表显示在生存曲线下面
           ggtheme = theme_bw(),
           tables.theme = theme_cleantable(),
           title="HDACIs combined with Chemotherapeutic drugs")

######
clusters_filter<-cutree(norma_result$tree_col,k=2)
sensitivity=colnames(HDAC_ic50[,clusters_filter==1])
resistance1=colnames(HDAC_ic50[,clusters_filter==2])
resistance=c("DAC-2","DAC-22","DAC-13","DAC-23")
medianl<-resistance1[-match(resistance,resistance1)]
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/7_clinical")
clinical<-read.table("31patient_clinical.txt",header = TRUE,sep = "\t")
clinical[,1]<-gsub("CAS-","",clinical[,1])


resistance_cli<-clinical[match(intersect(clinical[,1],resistance),clinical[,1]),]
sensitivity_cli<-clinical[match(intersect(clinical[,1],sensitivity),clinical[,1]),]
medianl_cli<-clinical[match(intersect(clinical[,1],medianl),clinical[,1]),]

all<-rbind(resistance_cli,sensitivity_cli,medianl_cli)
label<-as.matrix(c(rep("resistance",nrow(resistance_cli)),rep("sensitivity",nrow(sensitivity_cli)),rep("medianl",nrow(medianl_cli))),ncol=1)
TTP=as.numeric(all[,4])  ##生存时间
status_TTP=as.matrix(all[,3])  ##生存状态
status_TTP[which(status_TTP=="Disease Free")]=0
status_TTP[which(status_TTP=="Recurred")]=1
status_TTP<-as.numeric(status_TTP)
surv_info1<-as.data.frame(cbind(TTP,status_TTP))

surv_TTP<-survfit(Surv(TTP, status_TTP) ~ label,data=surv_info1)
ggsurvplot(surv_TTP,
           pval = TRUE, conf.int = TRUE, #是否显示p值/显示风险区域
           risk.table = TRUE,# 将风险表显示在生存曲线下面
           ggtheme = theme_bw(),
           tables.theme = theme_cleantable(),
           title="HDACIs combined with Chemotherapeutic drugs")


### NC文章
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/7_clinical")
clinical<-read.table("31patient_clinical.txt",header = TRUE,sep = "\t")


resistance_cli<-clinical[which(clinical[,2]=="Resistant"),]
sensitivity_cli<-clinical[which(clinical[,2]=="Sensitive"),]
medianl_cli<-clinical[which(clinical[,2]=="Intermediate"),]
all<-rbind(resistance_cli,sensitivity_cli)
label<-as.matrix(c(rep("resistance",nrow(resistance_cli)),rep("sensitivity",nrow(sensitivity_cli))),ncol=1)

all<-rbind(resistance_cli,sensitivity_cli,medianl_cli)
label<-as.matrix(c(rep("resistance",nrow(resistance_cli)),rep("sensitivity",nrow(sensitivity_cli)),rep("medianl",nrow(medianl_cli))),ncol=1)

all<-rbind(sensitivity_cli,medianl_cli)
label<-as.matrix(c(rep("sensitivity",nrow(sensitivity_cli)),rep("medianl",nrow(medianl_cli))),ncol=1)

TTP=as.numeric(all[,4])  ##生存时间
status_TTP=as.matrix(all[,3])  ##生存状态
status_TTP[which(status_TTP=="Disease Free")]=0
status_TTP[which(status_TTP=="Recurred")]=1
status_TTP<-as.numeric(status_TTP)
surv_info1<-as.data.frame(cbind(TTP,status_TTP))

surv_TTP<-survfit(Surv(TTP, status_TTP) ~ label,data=surv_info1)
#surv_TTP<-coxph(Surv(TTP, status_TTP) ~ label,data=surv_info1)
ggsurvplot(surv_TTP,
           pval = TRUE, conf.int = TRUE, #是否显示p值/显示风险区域
           risk.table = TRUE,# 将风险表显示在生存曲线下面
           ggtheme = theme_bw(),
           tables.theme = theme_cleantable(),
           title="HDACIs combined with Chemotherapeutic drugs")

###### 单药二等分
###单药
rm(list=ls())
setwd("~/xjj/drug/drug_result/drug64_AUC")
tresult_want<-read.table("drug64_sample.txt",sep="\t",header=T)
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/7_clinical")
clinical<-read.table("clinical.txt",header = TRUE,sep = "\t")
library(dplyr)
library(tidyr)
library(broom)
result<-data.frame(matrix(0,nrow(tresult_want),3))
foreach(k=1:nrow(tresult_want),.combine=rbind)%do%{
  file_name=as.character(tresult_want[k,1])
  cat(file_name,"\n")
  common_sample<-as.data.frame(tresult_want[match(file_name,tresult_want[,1]),])
  sensitivity<-unlist(strsplit(as.character(common_sample[1,5]), ","))
  resistance<-unlist(strsplit(as.character(common_sample[1,4]), ","))
  resistance_new<-as.character(sample_id[match(resistance,sample_id[,2]),1])
  sensitivity_new<-as.character(sample_id[match(sensitivity,sample_id[,2]),1])
  
  resistance_cli<-clinical[match(intersect(clinical[,1],resistance_new),clinical[,1]),]
  sensitivity_cli<-clinical[match(intersect(clinical[,1],sensitivity_new),clinical[,1]),]
  all<-rbind(resistance_cli,sensitivity_cli)
  label<-as.matrix(c(rep("resistance",nrow(resistance_cli)),rep("sensitivity",nrow(sensitivity_cli))),ncol=1)
  TTP=as.numeric(all[,3])  ##生存时间
  status_TTP=as.matrix(all[,2])  ##生存状态
  status_TTP[which(status_TTP=="Disease Free")]=0
  status_TTP[which(status_TTP=="Recurred")]=1
  status_TTP<-as.numeric(status_TTP)
  surv_info1<-as.data.frame(cbind(TTP,status_TTP))
  surv_TTP<-survfit(Surv(TTP, status_TTP) ~ label,data=surv_info1)
  coxp<-coxph(Surv(TTP, status_TTP) ~ label,data=surv_info1)
  result[k,1]<-file_name
  result[k,2]<-tidy(coxp)$estimate
  result[k,3]<-tidy(coxp)$p.value
  
  library(ggplot2)
  p<-ggsurvplot(surv_TTP,
                pval = TRUE, conf.int = TRUE, #是否显示p值/显示风险区域
                risk.table = TRUE,# 将风险表显示在生存曲线下面
                ggtheme = theme_bw(),
                tables.theme = theme_cleantable(),
                title=paste("    ",file_name,"resistance-sensitive",sep=" "))
  setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/7_clinical")
  ggsave(paste(file_name,"resistance-sensitive.pdf",sep=" "),p,width = 5, height = 5)
}  
colnames(result)<-c("geneid","estimate","p.value")
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/7_clinical")
write.table(result,file="each_drug_median.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')

###### 单药三等分
###单药
rm(list=ls())
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/7_clinical")
clinical<-read.table("31patient_clinical.txt",header = TRUE,sep = "\t")
clinical[,1]<-gsub("CAS-","",clinical[,1])

setwd("~/xjj/drug")
ic50<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
#ic50<-ic501[,-(which(colnames(ic501)==c("PC.100","PC.34")))]
colnames(ic50)<-sample_id[match(colnames(ic50),sample_id[,2]),1]

library(dplyr)
library(tidyr)
library(broom)
result<-data.frame(matrix(0,nrow(ic50),3))
foreach(k=1:nrow(ic50),.combine=rbind)%do%{
  file_name=as.character(rownames(ic50)[k])
  cat(file_name,"\n")
  ic50_sample<-ic50[match(file_name,rownames(ic50)),]
  resistance<-colnames(ic50_sample)[order(ic50_sample)[1:13]]
  sensitivity<-colnames(ic50_sample)[order(ic50_sample)[27:39]]
  resistance_cli<-clinical[match(intersect(clinical[,1],resistance),clinical[,1]),]
  sensitivity_cli<-clinical[match(intersect(clinical[,1],sensitivity),clinical[,1]),]
  all<-rbind(resistance_cli,sensitivity_cli)
  label<-as.matrix(c(rep("resistance",nrow(resistance_cli)),rep("sensitivity",nrow(sensitivity_cli))),ncol=1)
  TTP=as.numeric(all[,4])  ##生存时间
  status_TTP=as.matrix(all[,3])  ##生存状态
  status_TTP[which(status_TTP=="Disease Free")]=0
  status_TTP[which(status_TTP=="Recurred")]=1
  status_TTP<-as.numeric(status_TTP)
  surv_info1<-as.data.frame(cbind(TTP,status_TTP))
  surv_TTP<-survfit(Surv(TTP, status_TTP) ~ label,data=surv_info1)
  coxp<-coxph(Surv(TTP, status_TTP) ~ label,data=surv_info1)
  result[k,1]<-file_name
  result[k,2]<-tidy(coxp)$estimate
  result[k,3]<-tidy(coxp)$p.value
  #library(ggplot2)
  #p<-ggsurvplot(surv_TTP,
                #pval = TRUE, conf.int = TRUE, #是否显示p值/显示风险区域
                #risk.table = TRUE,# 将风险表显示在生存曲线下面
                #ggtheme = theme_bw(),
                #tables.theme = theme_cleantable(),
                #title=paste("    ",file_name,"resistance-sensitive",sep=" "))
  #setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/7_clinical")
  #ggsave(paste(file_name,"resistance-sensitive-30%.pdf",sep=" "),p,width = 5, height = 5)
}  
colnames(result)<-c("geneid","estimate","p.value")
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/7_clinical")
write.table(result,file="each_drug_30%.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')

##39个样本都放下去的时候，有5个药的预后有显著差别"GEM"   "S1129" "S1134" "S1648" "S7575"
##37个PDAC样本都放下去的时候，有3个药的预后有显著差别"GEM"   "S1134" "S1648"
##########  比较两份AUC
rm(list=ls())
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/7_clinical")
NC_AUC<-read.table("NC-AUC.txt",header = TRUE,sep = "\t")
colnames(NC_AUC)<-gsub("CAS-","",colnames(NC_AUC))

setwd("~/xjj/drug")
ic50<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
colnames(ic50)<-sample_id[match(colnames(ic50),sample_id[,2]),1]
library(sqldf)

result<-NULL
for(i in 1:64){
  NC<-NC_AUC[match(rownames(ic50)[i],NC_AUC[,1]),-1]
  NC_name<-gsub("\\.","-",gsub("CAS.","",colnames(NC)[order(NC)]))
  ic<-ic50[match(rownames(ic50)[i],rownames(ic50)),]
  ic_name<-colnames(ic)[order(ic,decreasing = TRUE)]
  result1<-rbind(NC_name,ic_name)
  result<-rbind(result,result1)
}
colnames(result)<-c(1:39)
length(which(!is.na(match(result[1,1:13],result[2,1:13]))))
length(which(!is.na(match(result[3,1:13],result[4,1:13]))))
d<-c()
for(i in 1:64){
  c1<-length(which(!is.na(match(result[(2*i-1),1:13],result[(2*i),1:13]))))
  d<-c(d,c1)
}
z<-c()
for(i in 1:64){
  c1<-length(which(!is.na(match(result[(2*i-1),14:26],result[(2*i),14:26]))))
  z<-c(z,c1)
}
g<-c()
for(i in 1:64){
  c1<-length(which(!is.na(match(result[(2*i-1),27:39],result[(2*i),27:39]))))
  g<-c(g,c1)
}

#################################################################
########## 结合AUC值和临床信息将五种化疗药物分成两组
rm(list=ls())
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
setwd("~/xjj/drug")
ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,!(colnames(ic501)%in%c("PC.100","PC.34"))]## all of PDAC
colnames(ic50)<-sample_id[match(colnames(ic50),sample_id[,2]),1]

setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/7_clinical")
clinical<-read.table("31patient_clinical.txt",header = TRUE,sep = "\t")
clinical[,1]<-gsub("CAS-","",clinical[,1])

clinical_sample<-intersect(colnames(ic50),clinical[,1]) ##PDAC with clinical information

ccc<-c("5-FU","GEM","IRI","OXA","PAC")
chemotherapy5<-ic50[match(ccc,rownames(ic50)),]
result<-matrix(0,5,38)
result_sen_res<-matrix(0,5,29)
for(i in 1:5){
  sensitive<-colnames(chemotherapy5)[which(as.numeric(chemotherapy5[i,])>median(as.numeric(chemotherapy5[i,])))]
  resistant<-colnames(chemotherapy5)[which(as.numeric(chemotherapy5[i,])<=median(as.numeric(chemotherapy5[i,])))]
  result[i,1]<-rownames(chemotherapy5)[i]
  result[i,2:19]<-sensitive
  result[i,20:38]<-resistant
  result_sen_res[i,match(intersect(sensitive,clinical_sample),clinical_sample)]<-"sen"
  result_sen_res[i,match(intersect(resistant,clinical_sample),clinical_sample)]<-"res"

}
colnames(result)<-c("drug",rep("sensitive",18),rep("resistant",19))
colnames(result_sen_res)<-clinical_sample
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/7_clinical")
write.table(result,file="chemotherapy5.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')
write.table(result_sen_res,file="chemotherapy5_29.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')

### 29个病人根据5个化疗药二分
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/7_clinical")
clinical<-read.table("31patient_clinical.txt",header = TRUE,sep = "\t")
clinical[,1]<-gsub("CAS-","",clinical[,1])

setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/7_clinical")
two<-read.table("29patient_clinical.txt",header = TRUE,sep = "\t")
resistance<-as.character(two[which(two[,2]=="res"),1])
sensitivity<-as.character(two[which(two[,2]=="sen"),1])
resistance_cli<-clinical[match(intersect(clinical[,1],resistance),clinical[,1]),]
sensitivity_cli<-clinical[match(intersect(clinical[,1],sensitivity),clinical[,1]),]
all<-rbind(resistance_cli,sensitivity_cli)
label<-as.matrix(c(rep("resistance",nrow(resistance_cli)),rep("sensitivity",nrow(sensitivity_cli))),ncol=1)
TTP=as.numeric(all[,4])  ##生存时间
status_TTP=as.matrix(all[,3])  ##生存状态
status_TTP[which(status_TTP=="Disease Free")]=0
status_TTP[which(status_TTP=="Recurred")]=1
status_TTP<-as.numeric(status_TTP)
surv_info1<-as.data.frame(cbind(TTP,status_TTP))

surv_TTP<-survfit(Surv(TTP, status_TTP) ~ label,data=surv_info1)
ggsurvplot(surv_TTP,
           pval = TRUE, conf.int = TRUE, #是否显示p值/显示风险区域
           risk.table = TRUE,# 将风险表显示在生存曲线下面
           ggtheme = theme_bw(),
           tables.theme = theme_cleantable(),
           title="5 Chemotherapeutic drugs")

