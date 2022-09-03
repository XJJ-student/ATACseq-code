library(MASS)
lm.ridge(formula, data, subset, na.action, lambda = 0, model = FALSE,
         x = FALSE, y = FALSE, contrasts = NULL, ...)

rm(list=ls())
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/2_ATACseq_DApeaks/FIMO_result/up0.01")
library(data.table)
fimo1<-fread("fimo.tsv",header=T,data.table=F)
colnames(fimo1)<-c(colnames(fimo1)[1:8],"FDR","matched_sequence")
fimo<-fimo1[which(fimo1$FDR<0.05),]
library(dplyr)
#TF_fimo<-fimo %>% distinct(motif_id,motif_alt_id,FDR, .keep_all = TRUE)
peak<-unique(fimo[,3])
TF<-unique(fimo[,2])
X_matrix<-matrix(0,length(peak),length(TF))
colnames(X_matrix)<-TF
rownames(X_matrix)<-peak

for(i in 1:length(TF)){
  TF_peak<-fimo[fimo[,2] %in% TF[i],3]
  X_matrix[match(unique(TF_peak),peak),i]<-1
}
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/2_ATACseq_DApeaks")
DApeak<-fread("sensitivity-resistance-all-DESeq2_HDAC_chemotherapy_AUC_ATAC.txt",header=T,data.table=F)
#DApeak<-fread(paste("sensitivity-resistance-all-DESeq2_",file_name,"_ATAC.txt",sep=""),header=T,data.table=F)
a<-gsub(":","_",peak)
b<-gsub("-","_",a)
Y<-DApeak[match(b,DApeak[,1]),3]

##### ridge regression
library(glmnet)
y<-Y
x<-X_matrix
lambdas <- 10^seq(2, -2, by = -.1)
#lambdas <- seq(0,5, length.out = 200)
fit = glmnet(x,y,alpha = 0,family = "gaussian",standardize=TRUE)
print(fit)
plot(fit)
summary(fit)

cv.fit <- cv.glmnet(x,y,alpha = 0,family = 'gaussian',standardize=TRUE,grouped=FALSE,nfolds = 5,lambda = lambdas)
plot(cv.fit)
coef(cv.fit,s = "lambda.min")


##check performance
#bestlam1=cv.fit$lambda.min
#bestlam2=cv.fit$lambda.1se
#y_predicted <- predict(cv.fit,s=c(bestlam1,bestlam2),newx=x)  #这样子会输出两列y'，分别来自于两种λ
y_predicted <- predict(cv.fit,s=cv.fit$lambda.min,newx=x)

# Sum of Squares Total and Error
sst <- sum((y - mean(y))^2)
sse <- sum((y_predicted - y)^2)
# R squared
rsq <- round(1 - sse / sst,3)
rsq


### 检查summary，看输出的数据是否正确
tmp_coeffs <- coef(cv.fit, s = "lambda.min")
output_coef=data.frame(features = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
output_coef<-output_coef[order(output_coef[,"coefficient"],decreasing=TRUE),]
output_coef





############################################################################
#############  一键完成  #####################
rm(list=ls())
library(dplyr)
library(data.table)
library(glmnet)
library(foreach)
library(ggplot2)
setwd("~/xjj/drug/drug_result/drug64_AUC")
tresult_want<-read.table("drug64_sample.txt",sep="\t",header=T)

FIMO_result<-data.frame(matrix(0,nrow(tresult_want),13))
foreach(k=1:nrow(tresult_want))%do%{
  file_name=as.character(tresult_want[k,1])
  cat(file_name,"\n")
  FIMO_result[k,1]<-file_name
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/4_TFs/FIMO_result/down",sep=""))
  fimo_down<-fread("fimo.tsv",header=T,data.table=F)
  colnames(fimo_down)<-c(colnames(fimo_down)[1:8],"FDR","matched_sequence")
  fimo_d<-fimo_down[which(fimo_down$FDR<0.05),]
  peak_d<-unique(fimo_d[,3])
  TF_d<-unique(fimo_d[,2])
  X_matrix_d<-matrix(0,length(peak_d),length(TF_d))
  colnames(X_matrix_d)<-TF_d
  rownames(X_matrix_d)<-peak_d
  for(i in 1:length(TF_d)){
  TF_peak_d<-fimo_d[fimo_d[,2] %in% TF_d[i],3]
  X_matrix_d[match(unique(TF_peak_d),peak_d),i]<-1
}
setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/2_ATACseq_DApeaks",sep=""))
DApeak<-fread(paste("sensitivity-resistance-all-DESeq2_",file_name,"_ATAC.txt",sep=""),header=T,data.table=F)
setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/1_RNAseq_DEGs",sep=""))
DE_gene<-read.table(paste("sensitivity-resistance-all-DESeq2_",file_name,"_RNAseq_DEGs.txt",sep=""),header=T,sep="\t")

a_d<-gsub(":","_",peak_d)
b_d<-gsub("-","_",a_d)
Y_d<-DApeak[match(b_d,DApeak[,1]),3]

##### ridge regression
y<-Y_d
x<-X_matrix_d
lambdas <- 10^seq(2, -2, by = -.1)
fit = glmnet(x,y,alpha = 0,family = "gaussian",standardize=TRUE)

cv.fit <- cv.glmnet(x,y,alpha = 0,family = 'gaussian',standardize=TRUE,grouped=FALSE,nfolds = 5,lambda = lambdas)
pd<-plot(cv.fit)
setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/4_TFs/FIMO_result/down",sep=""))
ggsave(paste("ridge_regression_",file_name,"lambda_down.pdf",sep=""),pd,width = 10, height = 4)

y_predicted <- predict(cv.fit,s=cv.fit$lambda.min,newx=x)

# Sum of Squares Total and Error
sst <- sum((y - mean(y))^2)
sse <- sum((y_predicted - y)^2)
# R squared
rsq_d <- round(1 - sse / sst,3)
FIMO_result[k,2]<-rsq_d #R squared sen vs ren down
### 检查summary，看输出的数据是否正确
tmp_coeffs <- coef(cv.fit, s = "lambda.min")
output_coef=data.frame(features = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
output_coef<-output_coef[order(output_coef[,"coefficient"],decreasing=TRUE),]
setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/4_TFs/FIMO_result/down",sep=""))
write.table(output_coef,file=paste(file_name,"_ridge_regression_down.txt",sep=""),sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')
res_positive1<-output_coef[which(output_coef$coefficient>0),1]
FIMO_result[k,3]<-length(res_positive1)
FIMO_result[k,4]<-paste0(res_positive1,collapse =",")
resSig<-DE_gene[which(DE_gene$pvalue<0.05),1]
down_int<-intersect(res_positive1,resSig)
FIMO_result[k,5]<-length(down_int)
FIMO_result[k,6]<-paste0(down_int,collapse =",")

setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/4_TFs/FIMO_result/up",sep=""))
fimo_up<-fread("fimo.tsv",header=T,data.table=F)
colnames(fimo_up)<-c(colnames(fimo_up)[1:8],"FDR","matched_sequence")
fimo_u<-fimo_up[which(fimo_up$FDR<0.05),]
peak_u<-unique(fimo_u[,3])
TF_u<-unique(fimo_u[,2])
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
setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/4_TFs/FIMO_result/up",sep=""))
ggsave(paste("ridge_regression_",file_name,"lambda_up.pdf",sep=""),pu,width = 10, height = 4)
y_predicted_u <- predict(cv.fit_u,s=cv.fit_u$lambda.min,newx=x)
# Sum of Squares Total and Error
sst <- sum((y - mean(y))^2)
sse <- sum((y_predicted_u - y)^2)
# R squared
rsq_u <- round(1 - sse / sst,3)
FIMO_result[k,7]<-rsq_u #R squared sen vs ren down
### 检查summary，看输出的数据是否正确
tmp_coeffs_u <- coef(cv.fit_u, s = "lambda.min")
output_coef_u=data.frame(features = tmp_coeffs_u@Dimnames[[1]][tmp_coeffs_u@i + 1], coefficient = tmp_coeffs_u@x)
output_coef_u<-output_coef_u[order(output_coef_u[,"coefficient"],decreasing=TRUE),]
setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/4_TFs/FIMO_result/up",sep=""))
write.table(output_coef_u,file=paste(file_name,"_ridge_regression_up.txt",sep=""),sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')
res_positive1_u<-output_coef_u[which(output_coef_u$coefficient>0),1]
FIMO_result[k,8]<-length(res_positive1_u)
FIMO_result[k,9]<-paste0(res_positive1_u,collapse =",")
up_int<-intersect(res_positive1_u,resSig)
FIMO_result[k,10]<-length(up_int)
FIMO_result[k,11]<-paste0(up_int,collapse =",")
int_DApeak<-intersect(res_positive1_u,res_positive1)
FIMO_result[k,12]<-length(int_DApeak)
FIMO_result[k,13]<-paste0(int_DApeak,collapse =",")
}

colnames(FIMO_result)<-c("drug","down_R_squared","length_of_down_positive_TF","down_positive_TF","int_down_length","int_down","up_R_squared","length_of_up_positive_TF","up_positive_TF","int_up_length","int_up","length_of_int_up_down","int_up_down")
setwd("~/xjj/drug/drug_result/drug64_AUC")
write.table(FIMO_result,"FIMO_result_64drug.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出



##################   画柱状图  #####################
rm(list=ls())
library(dplyr)
library(data.table)
library(glmnet)
library(foreach)
library(ggplot2)
setwd("~/xjj/drug/drug_result/drug64_AUC")
tresult_want<-read.table("drug64_sample.txt",sep="\t",header=T)

foreach(k=1:nrow(tresult_want))%do%{
  file_name=as.character(tresult_want[k,1])
  cat(file_name,"\n")
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/4_TFs/FIMO_result/up",sep=""))
  up_TF<-read.table(paste(file_name,"_ridge_regression_up.txt",sep=""),sep = '\t',header=T)
  up_TFs<-up_TF[which(up_TF$coefficient>0),]

setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/4_TFs/FIMO_result/down",sep=""))
down_TF<-read.table(paste(file_name,"_ridge_regression_down.txt",sep=""),sep = '\t',header=T)
down_TFs<-down_TF[which(down_TF$coefficient>0),]

int_TF1<-intersect(up_TFs[,1],down_TFs[,1]) ## the intersect of up and down TFs
int_TF<-c(int_TF1,"(Intercept)")

up_TFs1<-up_TFs[-(match(int_TF,up_TFs[,1])[which(!is.na(match(int_TF,up_TFs[,1])))]),] #up TFs
down_TFs1<-down_TFs[-(match(int_TF,down_TFs[,1])[which(!is.na(match(int_TF,down_TFs[,1])))]),] #down TFs

## 指定绘图顺序
up_TFs1$features <- factor(up_TFs1$features,
                          levels = rev(unique(up_TFs1$features)),
                          ordered = T)
down_TFs1$features <- factor(down_TFs1$features,
                           levels = rev(unique(down_TFs1$features)),
                           ordered = T)
## 绘制右侧的条形图
left_1 <- ggplot(up_TFs1,aes(x=features,y=coefficient*-1))+ #*-1 change direction
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
setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/4_TFs/FIMO_result",sep=""))
ggsave(paste("ridge_regression_",file_name,"TF_coefficient.pdf",sep=""),pu,width = 8, height = 6)
}


####################   上下调TFs中去掉重叠的TFs   #####################
rm(list=ls())
library(dplyr)
library(data.table)
library(glmnet)
library(foreach)
library(ggplot2)
setwd("~/xjj/drug/drug_result/drug64_AUC")
tresult_want<-read.table("drug64_sample.txt",sep="\t",header=T)

FIMO_result_clean<-data.frame(matrix(0,nrow(tresult_want),5))
foreach(k=1:nrow(tresult_want))%do%{
  file_name=as.character(tresult_want[k,1])
  cat(file_name,"\n")
  FIMO_result_clean[k,1]<-file_name
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/4_TFs/FIMO_result/up",sep=""))
  up_TF<-read.table(paste(file_name,"_ridge_regression_up.txt",sep=""),sep = '\t',header=T)
  up_TFs<-up_TF[which(up_TF$coefficient>0),]
  
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/4_TFs/FIMO_result/down",sep=""))
  down_TF<-read.table(paste(file_name,"_ridge_regression_down.txt",sep=""),sep = '\t',header=T)
  down_TFs<-down_TF[which(down_TF$coefficient>0),]
  
  int_TF1<-intersect(up_TFs[,1],down_TFs[,1])
  int_TF<-c(int_TF1,"(Intercept)")
  
  up_TFs1<-up_TFs[-(match(int_TF,up_TFs[,1])[which(!is.na(match(int_TF,up_TFs[,1])))]),]
  down_TFs1<-down_TFs[-(match(int_TF,down_TFs[,1])[which(!is.na(match(int_TF,down_TFs[,1])))]),]
  FIMO_result_clean[k,2]<-nrow(up_TFs1)
  FIMO_result_clean[k,3]<-paste0(up_TFs1[,1],collapse =",")
  FIMO_result_clean[k,4]<-nrow(down_TFs1)
  FIMO_result_clean[k,5]<-paste0(down_TFs1[,1],collapse =",")

  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/4_TFs/FIMO_result/up",sep=""))
  write.table(up_TFs1,file=paste(file_name,"_ridge_regression_up_clean.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/4_TFs/FIMO_result/down",sep=""))
  write.table(down_TFs1,file=paste(file_name,"_ridge_regression_down_clean.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
}
colnames(FIMO_result_clean)<-c("drug","length of up TFs","up TFs","length of down TFs","down TFs")
setwd("~/xjj/drug/drug_result/drug64_AUC")
write.table(FIMO_result_clean,"FIMO_result_clean_64drug.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出





############ 单独一种药
rm(list=ls())
library(dplyr)
library(data.table)
library(glmnet)
library(foreach)
library(ggplot2)

setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/2_ATACseq_DApeaks/FIMO_result/down0.01")
fimo_down<-fread("fimo.tsv",header=T,data.table=F)
colnames(fimo_down)<-c(colnames(fimo_down)[1:8],"FDR","matched_sequence")
fimo_d<-fimo_down[which(fimo_down$FDR<0.05),]
peak_d<-unique(fimo_d[,3])
TF_d<-unique(fimo_d[,2])
X_matrix_d<-matrix(0,length(peak_d),length(TF_d))
  colnames(X_matrix_d)<-TF_d
  rownames(X_matrix_d)<-peak_d
  for(i in 1:length(TF_d)){
    TF_peak_d<-fimo_d[fimo_d[,2] %in% TF_d[i],3]
    X_matrix_d[match(unique(TF_peak_d),peak_d),i]<-1
  }
  setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/2_ATACseq_DApeaks")
  DApeak<-fread("sensitivity-resistance-all-DESeq2_HDAC_chemotherapy_AUC_ATAC.txt",header=T,data.table=F)
  setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/1_RNAseq")
  DE_gene<-read.table("sensitivity-resistance-all-DESeq2_HDAC_chemotherapy_AUC_RNAseq_DEGs.txt",header=T,sep="\t")
  
  a_d<-gsub(":","_",peak_d)
  b_d<-gsub("-","_",a_d)
  Y_d<-DApeak[match(b_d,DApeak[,1]),3]
  
  ##### ridge regression
  y<-Y_d
  x<-X_matrix_d
  lambdas <- 10^seq(2, -2, by = -.1)
  fit = glmnet(x,y,alpha = 0,family = "gaussian",standardize=TRUE)
  
  cv.fit <- cv.glmnet(x,y,alpha = 0,family = 'gaussian',standardize=TRUE,grouped=FALSE,nfolds = 5,lambda = lambdas)
  pd<-plot(cv.fit)
  setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/2_ATACseq_DApeaks/FIMO_result/down0.01")
  ggsave(paste("ridge_regression_","0.01","_lambda_down.pdf",sep=""),pd,width = 10, height = 4)
  
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
  setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/2_ATACseq_DApeaks/FIMO_result/down0.01")
  write.table(output_coef,file="ridge_regression_down0.01.txt",sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  res_positive1<-output_coef[which(output_coef$coefficient>0),1]
  
  resSig<-DE_gene[which(DE_gene$pvalue<0.05),1]
  down_int<-intersect(res_positive1,resSig)
  
  
  setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/2_ATACseq_DApeaks/FIMO_result/up0.01")
  fimo_up<-fread("fimo.tsv",header=T,data.table=F)
  colnames(fimo_up)<-c(colnames(fimo_up)[1:8],"FDR","matched_sequence")
  fimo_u<-fimo_up[which(fimo_up$FDR<0.05),]
  peak_u<-unique(fimo_u[,3])
  TF_u<-unique(fimo_u[,2])
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
  setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/2_ATACseq_DApeaks/FIMO_result/up0.01")
  ggsave(paste("ridge_regression_","0.01","_lambda_up.pdf",sep=""),pu,width = 10, height = 4)
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
  setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/2_ATACseq_DApeaks/FIMO_result/up0.01")
  write.table(output_coef_u,file="ridge_regression_up0.01.txt",sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  res_positive1_u<-output_coef_u[which(output_coef_u$coefficient>0),1]
  up_int<-intersect(res_positive1_u,resSig)
  int_DApeak<-intersect(res_positive1_u,res_positive1)
  

##################   画柱状图  #####################
rm(list=ls())
library(dplyr)
library(data.table)
library(glmnet)
library(foreach)
library(ggplot2)
  setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/2_ATACseq_DApeaks/FIMO_result/up0.01")
  up_TF<-read.table("ridge_regression_up0.01.txt",sep = '\t',header=T)
  up_TFs<-up_TF[which(up_TF$coefficient>0),]
  
  setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/2_ATACseq_DApeaks/FIMO_result/down0.01")
  down_TF<-read.table("ridge_regression_down0.01.txt",sep = '\t',header=T)
  down_TFs<-down_TF[which(down_TF$coefficient>0),]
  
  int_TF1<-intersect(up_TFs[,1],down_TFs[,1]) ## the intersect of up and down TFs
  int_TF<-c(int_TF1,"(Intercept)")
  
  up_TFs1<-up_TFs[-(match(int_TF,up_TFs[,1])[which(!is.na(match(int_TF,up_TFs[,1])))]),] #up TFs
  down_TFs1<-down_TFs[-(match(int_TF,down_TFs[,1])[which(!is.na(match(int_TF,down_TFs[,1])))]),] #down TFs
  down_TFs1[5,1]<-"TFAP2C"
  down_TFs1[6,1]<-"Rarb"
  ## 指定绘图顺序
  up_TFs1$features <- factor(up_TFs1$features,
                             levels = rev(unique(up_TFs1$features)),
                             ordered = T)
  down_TFs1$features <- factor(down_TFs1$features,
                               levels = rev(unique(down_TFs1$features)),
                               ordered = T)
  ## 绘制右侧的条形图
  left_1 <- ggplot(up_TFs1,aes(x=features,y=coefficient*-1))+ #*-1 change direction
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
  setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/2_ATACseq_DApeaks/FIMO_result")
  ggsave(paste("ridge_regression_","DApeakP0.01","TF_coefficient.pdf",sep=""),pu,width = 8, height = 6)



####################   上下调TFs中去掉重叠的TFs   #####################
rm(list=ls())
library(dplyr)
library(data.table)
library(glmnet)
library(foreach)
library(ggplot2)
setwd("~/xjj/drug/drug_result/drug64_AUC")
tresult_want<-read.table("drug64_sample.txt",sep="\t",header=T)

FIMO_result_clean<-data.frame(matrix(0,nrow(tresult_want),5))
foreach(k=1:nrow(tresult_want))%do%{
  file_name=as.character(tresult_want[k,1])
  cat(file_name,"\n")
  FIMO_result_clean[k,1]<-file_name
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/4_TFs/FIMO_result/up",sep=""))
  up_TF<-read.table(paste(file_name,"_ridge_regression_up.txt",sep=""),sep = '\t',header=T)
  up_TFs<-up_TF[which(up_TF$coefficient>0),]
  
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/4_TFs/FIMO_result/down",sep=""))
  down_TF<-read.table(paste(file_name,"_ridge_regression_down.txt",sep=""),sep = '\t',header=T)
  down_TFs<-down_TF[which(down_TF$coefficient>0),]
  
  int_TF1<-intersect(up_TFs[,1],down_TFs[,1])
  int_TF<-c(int_TF1,"(Intercept)")
  
  up_TFs1<-up_TFs[-(match(int_TF,up_TFs[,1])[which(!is.na(match(int_TF,up_TFs[,1])))]),]
  down_TFs1<-down_TFs[-(match(int_TF,down_TFs[,1])[which(!is.na(match(int_TF,down_TFs[,1])))]),]
  FIMO_result_clean[k,2]<-nrow(up_TFs1)
  FIMO_result_clean[k,3]<-paste0(up_TFs1[,1],collapse =",")
  FIMO_result_clean[k,4]<-nrow(down_TFs1)
  FIMO_result_clean[k,5]<-paste0(down_TFs1[,1],collapse =",")
  
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/4_TFs/FIMO_result/up",sep=""))
  write.table(up_TFs1,file=paste(file_name,"_ridge_regression_up_clean.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/4_TFs/FIMO_result/down",sep=""))
  write.table(down_TFs1,file=paste(file_name,"_ridge_regression_down_clean.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
}
colnames(FIMO_result_clean)<-c("drug","length of up TFs","up TFs","length of down TFs","down TFs")
setwd("~/xjj/drug/drug_result/drug64_AUC")
write.table(FIMO_result_clean,"FIMO_result_clean_64drug.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出





