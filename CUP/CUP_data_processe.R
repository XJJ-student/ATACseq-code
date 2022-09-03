#############     TCGA   Fresh    ###########################
rm(list=ls())
setwd("~/xjj/CUP")
first_category_name = list.files("TCGA")  
dir = paste("~/xjj/CUP/TCGA/",first_category_name,sep="")  
n = length(dir)
library(foreach)
library("data.table")

result<-data.frame(matrix(0,n,7))
foreach(i=1:n,.combine='rbind')%do%{
  setwd(dir[i])
  clinical<-fread(paste("TCGA-",first_category_name[i],".GDC_phenotype.tsv",sep=""),header=T,data.table=F)
  sampletype<-substring(clinical[,1],16)
  clinical_want<-clinical[sampletype%in%"A",]##The clinical samples of fresh
  exp<-fread("convert_exp.txt",header=T,data.table=F)
  exp1<-exp[,-1] #exp matrix
  geneid<-exp[,1] #exp geneid
  clincal_exp<-intersect(substring(clinical_want[,1],1,15),colnames(exp1))##the exp of fresh sample
  int_exp<-exp1[,match(clincal_exp,colnames(exp1))]######the exp of fresh sample
  colnames_int_exp<-substring(colnames(int_exp),14,15)
  primary_sample<-colnames(int_exp)[which(colnames_int_exp=="01"|colnames_int_exp=="03")]
  clinical_want_primary<-clinical_want[match(primary_sample,substring(clinical_want[,1],1,15)),]##在编码是01和03的条件下看看type是否为"primary"
  primary<-int_exp[,which(colnames_int_exp=="01"|colnames_int_exp=="03")]
  primary1<-cbind(geneid,primary)
  colnames(primary1)<-c("geneid",colnames(primary))
  write.table(primary1,file=paste(first_category_name[i],"_primary_exp.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  if(sum(colnames_int_exp=="06")==1){
    metastatic<-as.matrix(exp1[,which(colnames_int_exp=="06")])
    metastatic1<-cbind(geneid,metastatic)
    colnames(metastatic1)<-c("geneid",colnames(exp1)[which(colnames_int_exp=="06")])
    write.table(metastatic1,file=paste(first_category_name[i],"_metastatic_exp.txt",sep=""),sep = '\t',
                col.names = T,row.names = F,quote = FALSE,na='')
  }
  if(sum(colnames_int_exp=="06")!=1){
    metastatic<-exp1[,which(colnames_int_exp=="06")]
    metastatic1<-cbind(geneid,metastatic)
    colnames(metastatic1)<-c("geneid",colnames(metastatic))
    write.table(metastatic1,file=paste(first_category_name[i],"_metastatic_exp.txt",sep=""),sep = '\t',
                col.names = T,row.names = F,quote = FALSE,na='')
  }
  cat(first_category_name[i],"\n")
  cat(ncol(exp1),"\n")
  cat(length(clincal_exp),"\n")
  cat(sum(colnames_int_exp=="01"|colnames_int_exp=="03"),"\n")
  cat(sum(colnames_int_exp=="06"),"\n")
  cat(sum(clinical_want_primary$sample_type.samples=="Primary Tumor"|clinical_want_primary$sample_type.samples=="Primary Blood Derived Cancer - Peripheral Blood"),"\n")
  
  result[i,1]<-first_category_name[i]
  result[i,2]<-ncol(exp1) ## all case
  result[i,3]<-length(clincal_exp)  ## all case of clincal
  result[i,4]<-sum(colnames_int_exp=="01"|colnames_int_exp=="03") #primary
  result[i,5]<-sum(colnames_int_exp=="06") #metastatic
  result[i,6]<-sum(clinical_want_primary$sample_type.samples=="Primary Tumor"|clinical_want_primary$sample_type.samples=="Primary Blood Derived Cancer - Peripheral Blood")
  result[i,7]<-length(unique(geneid))
}  
colnames(result)<-c("cancer","the number of all case","the number of clinical case","primary","metastatic","Numbers of Primary Tumor","Numbers of genes")
setwd("~/xjj/CUP")
write.table(result,"TCGA_primary_metastatic_sample_fresh.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出


#############     TCGA   FFPE    ###########################
rm(list=ls())
setwd("~/xjj/CUP")
first_category_name = list.files("TCGA")  
dir = paste("~/xjj/CUP/TCGA/",first_category_name,sep="")  
n = length(dir)
library(foreach)
library("data.table")

result<-data.frame(matrix(0,n,7))
foreach(i=1:n,.combine='rbind')%do%{
  setwd(dir[i])
  clinical<-fread(paste("TCGA-",first_category_name[i],".GDC_phenotype.tsv",sep=""),header=T,data.table=F)
  sampletype<-substring(clinical[,1],16)
  clinical_want<-clinical[sampletype%in%"B",]##The clinical samples of fresh
  exp<-fread("convert_exp.txt",header=T,data.table=F)
  exp1<-exp[,-1] #exp matrix
  geneid<-exp[,1] #exp geneid
  clincal_exp<-intersect(substring(clinical_want[,1],1,15),colnames(exp1))##the exp of fresh sample
  
  if(length(clincal_exp)==1){
    int_exp<-as.matrix(exp1[,match(clincal_exp,colnames(exp1))])######the exp of fresh sample
    colnames(int_exp)<-colnames(exp1)[match(clincal_exp,colnames(exp1))]
    colnames_int_exp<-substring(colnames(int_exp),14,15)
    primary_sample<-colnames(int_exp)[which(colnames_int_exp=="01"|colnames_int_exp=="03")]
    clinical_want_primary<-clinical_want[match(primary_sample,substring(clinical_want[,1],1,15)),]##在编码是01和03的条件下看看type是否为"primary"
    primary<-int_exp[,which(colnames_int_exp=="01"|colnames_int_exp=="03")]
    primary1<-cbind(geneid,primary)
    colnames(primary1)<-c("geneid",colnames(int_exp))
    write.table(primary1,file=paste(first_category_name[i],"_primary_exp_FFPE.txt",sep=""),sep = '\t',
                col.names = T,row.names = F,quote = FALSE,na='')
  }
  if(length(clincal_exp)!=1){
  int_exp<-exp1[,match(clincal_exp,colnames(exp1))]######the exp of fresh sample
  colnames_int_exp<-substring(colnames(int_exp),14,15)
  primary_sample<-colnames(int_exp)[which(colnames_int_exp=="01"|colnames_int_exp=="03")]
  clinical_want_primary<-clinical_want[match(primary_sample,substring(clinical_want[,1],1,15)),]##在编码是01和03的条件下看看type是否为"primary"
  primary<-int_exp[,which(colnames_int_exp=="01"|colnames_int_exp=="03")]
  primary1<-cbind(geneid,primary)
  colnames(primary1)<-c("geneid",colnames(primary))
  write.table(primary1,file=paste(first_category_name[i],"_primary_exp_FFPE.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  }
  if(sum(colnames_int_exp=="06")==1){
    metastatic<-as.matrix(exp1[,which(colnames_int_exp=="06")])
    metastatic1<-cbind(geneid,metastatic)
    colnames(metastatic1)<-c("geneid",colnames(exp1)[which(colnames_int_exp=="06")])
    write.table(metastatic1,file=paste(first_category_name[i],"_metastatic_exp_FFPE.txt",sep=""),sep = '\t',
                col.names = T,row.names = F,quote = FALSE,na='')
  }
  if(sum(colnames_int_exp=="06")!=1){
    metastatic<-exp1[,which(colnames_int_exp=="06")]
    metastatic1<-cbind(geneid,metastatic)
    colnames(metastatic1)<-c("geneid",colnames(metastatic))
    write.table(metastatic1,file=paste(first_category_name[i],"_metastatic_exp_FFPE.txt",sep=""),sep = '\t',
                col.names = T,row.names = F,quote = FALSE,na='')
  }
  cat(first_category_name[i],"\n")
  cat(ncol(exp1),"\n")
  cat(length(clincal_exp),"\n")
  cat(sum(colnames_int_exp=="01"|colnames_int_exp=="03"),"\n")
  cat(sum(colnames_int_exp=="06"),"\n")
  cat(sum(clinical_want_primary$sample_type.samples=="Primary Tumor"|clinical_want_primary$sample_type.samples=="Primary Blood Derived Cancer - Peripheral Blood"),"\n")
  
  result[i,1]<-first_category_name[i]
  result[i,2]<-ncol(exp1) ## all case
  result[i,3]<-length(clincal_exp)  ## all case of clincal
  result[i,4]<-sum(colnames_int_exp=="01"|colnames_int_exp=="03") #primary
  result[i,5]<-sum(colnames_int_exp=="06") #metastatic
  result[i,6]<-sum(clinical_want_primary$sample_type.samples=="Primary Tumor"|clinical_want_primary$sample_type.samples=="Primary Blood Derived Cancer - Peripheral Blood")
  result[i,7]<-length(unique(geneid))
}  
colnames(result)<-c("cancer","the number of all case","the number of clinical case","primary","metastatic","Numbers of Primary Tumor","Numbers of genes")
setwd("~/xjj/CUP")
write.table(result,"TCGA_primary_metastatic_sample_FFPE.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出


#############     这一步先不进行 TCGA COAD+READ   #############################
rm(list=ls())
setwd("~/xjj/CUP/TCGA/COAD")
COAD<-fread("COAD_primary_exp.txt",header=T,data.table=F)
COAD_exp<-COAD[,-1]
COAD_gene<-COAD[,1]
setwd("~/xjj/CUP/TCGA/READ")
READ<-fread("READ_primary_exp.txt",header=T,data.table=F)
READ_exp<-READ[,-1]
READ_gene<-READ[,1]
gene<-intersect(READ_gene,COAD_gene)
exp<-cbind(gene,COAD_exp,READ_exp)
colnames(exp)<-c("geneid",colnames(COAD_exp),colnames(READ_exp))
setwd("~/xjj/CUP/TCGA/COADREAD")
write.table(exp,"COADREAD_primary_exp.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

setwd("~/xjj/CUP/TCGA/COADREAD")
COAD<-fread("TCGA-COAD.GDC_phenotype.tsv",header=T,data.table=F)
READ<-fread("TCGA-READ.GDC_phenotype.tsv",header=T,data.table=F)
int_sam<-intersect(colnames(COAD),colnames(READ))
COAD_exp<-COAD[,match(int_sam,colnames(COAD))]
READ_exp<-READ[,match(int_sam,colnames(READ))]

exp<-rbind(COAD_exp,READ_exp)
colnames(exp)<-colnames(READ_exp)
write.table(exp,"TCGA-COADREAD.GDC_phenotype.tsv",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

###########################################################################
##########################     only exp    ################################
#############     TCGA COAD+READ  convert_exp #############################
rm(list=ls())
setwd("~/xjj/CUP/TCGA/COAD")
COAD<-fread("convert_exp.txt",header=T,data.table=F)
COAD_exp<-COAD[,-1]
COAD_gene<-COAD[,1]
setwd("~/xjj/CUP/TCGA/READ")
READ<-fread("convert_exp.txt",header=T,data.table=F)
READ_exp<-READ[,-1]
READ_gene<-READ[,1]
gene<-intersect(READ_gene,COAD_gene)
exp<-cbind(gene,COAD_exp,READ_exp)
colnames(exp)<-c("geneid",colnames(COAD_exp),colnames(READ_exp))
setwd("~/xjj/CUP/TCGA/COADREAD")
write.table(exp,"convert_exp.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出


#################  单单从表达谱的标签去分组，没有结合临床数据  ######
rm(list=ls())
setwd("~/xjj/CUP")
first_category_name = list.files("TCGA")  
dir = paste("~/xjj/CUP/TCGA/",first_category_name,sep="")  
n = length(dir)
library(foreach)
library("data.table")

result<-data.frame(matrix(0,n,5))
foreach(i=1:n,.combine='rbind')%do%{
  setwd(dir[i])
  exp<-fread("convert_exp.txt",header=T,data.table=F)
  exp1<-exp[,-1] #exp matrix
  geneid<-exp[,1] #exp geneid
  colnames_int_exp<-substring(colnames(exp1),14,15)
  primary<-exp1[,which(colnames_int_exp=="01"|colnames_int_exp=="03")]
  primary1<-cbind(geneid,primary)
  colnames(primary1)<-c("geneid",colnames(primary))
  write.table(primary1,file=paste(first_category_name[i],"_primary_exp_noclinical.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  if(sum(colnames_int_exp=="06")==1){
    metastatic<-as.matrix(exp1[,which(colnames_int_exp=="06")])
    metastatic1<-cbind(geneid,metastatic)
    colnames(metastatic1)<-c("geneid",colnames(exp1)[which(colnames_int_exp=="06")])
    write.table(metastatic1,file=paste(first_category_name[i],"_metastatic_exp.txt",sep=""),sep = '\t',
                col.names = T,row.names = F,quote = FALSE,na='')
  }
  if(sum(colnames_int_exp=="06")!=1){
    metastatic<-exp1[,which(colnames_int_exp=="06")]
    metastatic1<-cbind(geneid,metastatic)
    colnames(metastatic1)<-c("geneid",colnames(metastatic))
    write.table(metastatic1,file=paste(first_category_name[i],"_metastatic_exp.txt",sep=""),sep = '\t',
                col.names = T,row.names = F,quote = FALSE,na='')
  }
  cat(first_category_name[i],"\n")
  cat(ncol(exp1),"\n")
  cat(sum(colnames_int_exp=="01"|colnames_int_exp=="03"),"\n")
  cat(sum(colnames_int_exp=="06"),"\n")
  
  result[i,1]<-first_category_name[i]
  result[i,2]<-ncol(exp1)  ## all case
  result[i,3]<-sum(colnames_int_exp=="01"|colnames_int_exp=="03") #primary
  result[i,4]<-sum(colnames_int_exp=="06") #metastatic
  result[i,5]<-length(unique(geneid))
}  
colnames(result)<-c("cancer","the number of all case","primary","metastatic","Numbers of genes")
setwd("~/xjj/CUP")
write.table(result,"TCGA_primary_metastatic_sample_noclinical.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出




#############   TCGA删除0值   ########
rm(list=ls())
setwd("~/xjj/CUP")
first_category_name = list.files("TCGA")  
dir = paste("~/xjj/CUP/TCGA/",first_category_name,sep="")  
n = length(dir)
library(foreach)
library("data.table")

result<-data.frame(matrix(0,n,4))
foreach(i=1:n,.combine='rbind')%do%{
  setwd(dir[i])
  primary<-fread(paste(first_category_name[i],"_primary_exp.txt",sep=""),header=T,data.table=F)
  exp2<-primary[,-1]
  delete<-apply(exp2,1,function(x) mean(x==0))
  cou<-which(delete>0.5)####在超过一半样本以上的基因删去
  exp<-exp2[(-cou),]
  geneid<-primary[(-cou),1]
  exp3<-cbind(geneid,exp)
  colnames(exp3)<-c("geneid",colnames(exp))
  write.table(exp3,file=paste(first_category_name[i],"_primary_exp_delete50.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  
  cat(first_category_name[i],"\n")
  cat(length(cou),"\n")
  cat(nrow(exp),"\n")
  
  result[i,1]<-first_category_name[i]
  result[i,2]<-nrow(primary)  ##origenal
  result[i,3]<-length(cou)    ##delete
  result[i,4]<-nrow(exp)      ##residue
}
colnames(result)<-c("cancer","the number of origenal","Numbers of delete","Numbers of residue")
setwd("~/xjj/CUP")
write.table(result,"TCGA_delete50_primary_sample.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

############     TCGA所有的表达谱merge   #######################
rm(list=ls())
setwd("~/xjj/CUP")
first_category_name = list.files("TCGA")  
dir = paste("~/xjj/CUP/TCGA/",first_category_name,sep="")  
n = length(dir)
library(foreach)
library("data.table")
i=1
setwd(dir[i])
primary<-fread(paste(first_category_name[i],"_primary_exp_delete50.txt",sep=""),header=T,data.table=F)
geneid<-primary[,1]

gene<-geneid
for(i in 1:(n-1)){
  i=i+1
  setwd(dir[i])
  primary1<-fread(paste(first_category_name[i],"_primary_exp_delete50.txt",sep=""),header=T,data.table=F)
  geneid1<-primary1[,1]
  gene<-intersect(gene,geneid1)
}
int_gene<-data.frame(gene)
colnames(int_gene)<-"TCGA_all_cancer_geneid"
setwd("~/xjj/CUP")
write.table(int_gene,"TCGA_all_cancer_geneid.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

#############    TCGA exp    ###########################
rm(list=ls())
setwd("~/xjj/CUP")
int_gene<-as.matrix(read.delim("TCGA_all_cancer_geneid.txt",sep = '\t',header = T))##数据输出
setwd("~/xjj/CUP")
first_category_name = list.files("TCGA")  
dir = paste("~/xjj/CUP/TCGA/",first_category_name,sep="")  
n = length(dir)
library(foreach)
library("data.table")

foreach(i=1:n,.combine='rbind')%do%{
  setwd(dir[i])
  primary1<-fread(paste(first_category_name[i],"_primary_exp_delete50.txt",sep=""),header=T,data.table=F)
  exp<-primary1[match(int_gene,primary1[,1]),]
  write.table(exp,file=paste(first_category_name[i],"_primary_exp_delete50_consistent_gene.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
}

#############   使用uwot包对TCGA数据进行UMAP降维可视化分析   ##########
rm(list=ls())
setwd("~/xjj/CUP")
first_category_name = list.files("TCGA")  
dir = paste("~/xjj/CUP/TCGA/",first_category_name,sep="")  
n = length(dir)
library(foreach)
library("data.table")
i=1
setwd(dir[i])
primary<-fread(paste(first_category_name[i],"_primary_exp_delete50_consistent_gene.txt",sep=""),header=T,data.table=F)
exp<-primary[,-1]
gene<-as.matrix(primary[,1])
label<-rep(first_category_name[i],ncol(exp))

expp<-exp
labell<-label
for(i in 1:(n-1)){
  i=i+1
  setwd(dir[i])
  primary1<-fread(paste(first_category_name[i],"_primary_exp_delete50_consistent_gene.txt",sep=""),header=T,data.table=F)
  exp1<-primary1[,-1]
  label1<-rep(first_category_name[i],ncol(exp1))
  expp<-cbind(expp,exp1)
  labell<-c(labell,label1)
}

expp_sum<-apply(expp,1,sum)
exppp<-expp[order(expp_sum,decreasing=TRUE)[1:floor(length(expp_sum)*0.25)],]
gene_sum<-gene[order(expp_sum,decreasing=TRUE)[1:floor(length(expp_sum)*0.25)],1]
expb<-data.frame(t(exppp))
labelll<-as.data.frame(matrix(labell,ncol=1))
expbbb<-cbind(expb,labelll)
#expbbb<-expbb[-(which(labelll=="COAD"|labelll=="READ")),]
colnames(expbbb)<-c(gene_sum,"cancertype")
colnames(expbbb)<-c(gene,"cancertype")

library(umap)
iris.data = expbbb[,-(ncol(expbbb))]
iris.umap = umap::umap(iris.data)
head(iris.umap$layout)
iris_sumap_res <- data.frame(iris.umap$layout,cancertype=expbbb$cancertype)
head(iris_sumap_res)

library(uwot)
set.seed(100)
iris_umap <- uwot::umap(expbbb,min_dist = 0.1)
iris_sumap_res <- data.frame(iris_umap,cancertype=expbbb$cancertype)
colnames(iris_sumap_res)<-c("UMAP1","UMAP2","cancertype")
head(iris_sumap_res)
library(ggplot2)
ggplot(iris_sumap_res,aes(UMAP1,UMAP2,color=cancertype)) + 
  geom_point(size = 0.1) +theme_bw() + 
  #geom_hline(yintercept = 0,lty=2,col="red") + 
  #geom_vline(xintercept = 0,lty=2,col="blue",lwd=1) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x="UMAP_1",y="UMAP_2",
       title = "A UMAP visualization of the TCGA dataset")
geom_text(aes(label = cancertype),size = 3,check_overlap = TRUE,nudge_y = (-0.1))

p<-ggplot(iris_sumap_res,aes(UMAP1,UMAP2,color=cancertype)) + 
  geom_point(size = 3) +theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x="UMAP_1",y="UMAP_2",
       title = paste("A UMAP visualization of the ","HDAC",sep=""))+
  geom_text(aes(label = AUC),colour="black",size = 2.5,check_overlap = TRUE,nudge_y = (-0.1))
setwd("~/xjj/drug/drug_result/drug64_AUC/figure_feature/UMAP_Pmin_Top1000_DApeaks_Gradual_change")
ggsave(paste("HDAC","_Top1000_DApeaks.pdf",sep=""),p,width = 6, height = 5)




library(uwot)
iris_sumap <- uwot::umap(expbbb, n_neighbors = 15, min_dist = 0.01,
                         y = expbbb$cancertype, target_weight = 0.5)
head(iris_sumap)
iris_sumap_res <- data.frame(iris_sumap,cancertype=expbbb$cancertype)
head(iris_sumap_res)
library(ggplot2)

ggplot(iris_sumap_res,aes(X1,X2,color=cancertype)) + 
  geom_point(size = 0.1) +theme_bw() + 
  #geom_hline(yintercept = 0,lty=2,col="red") + 
  #geom_vline(xintercept = 0,lty=2,col="blue",lwd=1) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x="UMAP_1",y="UMAP_2",
       title = "A UMAP visualization of the TCGA dataset")
  geom_text(aes(label = cancertype),size = 3,check_overlap = TRUE,angle = 30)

p<-ggplot(iris_sumap_res,aes(X1,X2,color=cancertype)) + 
  geom_point(size = 0.1) +theme_bw() + 
  #geom_hline(yintercept = 0,lty=2,col="red") + 
  #geom_vline(xintercept = 0,lty=2,col="blue",lwd=1) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x="UMAP_1",y="UMAP_2",
       title = "A UMAP visualization of the TCGA dataset")+
  geom_text(aes(label = cancertype),size = 3,check_overlap = TRUE,angle = 30)
setwd("~/xjj/CUP")
ggsave("mean025_0.1.pdf",p,width = 12, height = 10)



##############   PCA
rm(list=ls())
setwd("~/xjj/CUP")
first_category_name = list.files("TCGA")  
dir = paste("~/xjj/CUP/TCGA/",first_category_name,sep="")  
n = length(dir)
library(foreach)
library("data.table")
i=1
setwd(dir[i])
primary<-fread(paste(first_category_name[i],"_primary_exp_delete50_consistent_gene.txt",sep=""),header=T,data.table=F)
exp<-primary[,-1]
gene<-as.matrix(primary[,1])
label<-rep(first_category_name[i],ncol(exp))

expp<-exp
labell<-label
for(i in 1:(n-1)){
  i=i+1
  setwd(dir[i])
  primary1<-fread(paste(first_category_name[i],"_primary_exp_delete50_consistent_gene.txt",sep=""),header=T,data.table=F)
  exp1<-primary1[,-1]
  label1<-rep(first_category_name[i],ncol(exp1))
  expp<-cbind(expp,exp1)
  labell<-c(labell,label1)
}

expp_sum<-apply(expp,1,sum)
exppp<-expp[order(expp_sum,decreasing=TRUE)[1:floor(length(expp_sum)*0.25)],]
gene_sum<-gene[order(expp_sum,decreasing=TRUE)[1:floor(length(expp_sum)*0.25)],1]
expb<-data.frame(t(exppp))
#expb<-data.frame(t(expp))
labelll<-as.data.frame(matrix(labell,ncol=1))
expbbb<-cbind(expb,labelll)
colnames(expbbb)<-c(gene_sum,"cancertype")
colnames(expbbb)<-c(gene,"cancertype")



set.seed(42)
library(Rtsne)
tsne_out <- Rtsne(as.matrix(expbbb[,-(ncol(expbbb))]),pca=FALSE,dims=2,
                  perplexity=10,theta=0.0) # Run TSNE
iris_sumap_res <- data.frame(tsne_out$Y,cancertype=expbbb$cancertype)
colnames(iris_sumap_res)<-c("tSNE1","tSNE2","cancertype")
head(iris_sumap_res)
library(ggplot2)
p<-ggplot(iris_sumap_res,aes(tSNE1,tSNE2,color=cancertype)) + 
  geom_point(size = 0.5) +theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x="tSNE_1",y="tSNE_2",
       title = "A tSNE visualization of the TCGA")
  geom_text(aes(label = cancertype),colour="black",size = 2.5,check_overlap = TRUE,nudge_y = (-3))


set.seed(1000)
df_pca <- prcomp(as.matrix(expbbb[,-(ncol(expbbb))])) # Run PCA princomp(Inx,cor=T)
df_pcs <-data.frame(df_pca$x,cancertype=expbbb$cancertype)
iris_umap<- data.frame(df_pcs[,c(1,2)])
iris_sumap_res <- data.frame(df_pcs[,c(1,2,38)])
colnames(iris_sumap_res)<-c("PC1","PC2","cancertype")
head(iris_sumap_res)
library(ggplot2)
p<-ggplot(iris_sumap_res,aes(PC1,PC2,color=cancertype)) + 
  geom_point(size = 3) +theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x="PCA_1",y="PCA_2",
       title = "A PCA visualization of the TCGA dataset")+
  geom_text(aes(label = cancertype),colour="black",size = 2.5,check_overlap = TRUE,nudge_y = (-0.5))



#############     Feature gene selection    #############################
rm(list=ls())
setwd("~/xjj/CUP")
first_category_name = list.files("TCGA")  
dir = paste("~/xjj/CUP/TCGA/",first_category_name,sep="")  
n = length(dir)
library(foreach)
library("data.table")
i=1
setwd(dir[i])
primary<-fread(paste(first_category_name[i],"_primary_exp_delete50_consistent_gene.txt",sep=""),header=T,data.table=F)
exp<-primary[,-1]
gene<-as.matrix(primary[,1])
label<-rep(first_category_name[i],ncol(exp))

expp<-exp
labell<-label
for(i in 1:(n-1)){
  i=i+1
  setwd(dir[i])
  primary1<-fread(paste(first_category_name[i],"_primary_exp_delete50_consistent_gene.txt",sep=""),header=T,data.table=F)
  exp1<-primary1[,-1]
  label1<-rep(first_category_name[i],ncol(exp1))
  expp<-cbind(expp,exp1)
  labell<-c(labell,label1)
}
colnames(expp)<-labell
#####  RankComp
freq<-0.99
result<-NULL
for(j in 1:n){
  case_exp<-as.matrix(apply(expp[,colnames(expp)%in%first_category_name[j]],1,median))
  control_exp<-as.matrix(apply(expp[,-(which(colnames(expp)==first_category_name[j]))],1,median))
  outlier_dir<-NULL;
  outlier_pvalue<-NULL;
  Lgene<-dim(gene)[1];
  for (k in 1:Lgene){###开始循环每个基因
    print(k)
    Nnorm=control_exp[k,]##Nnorm是第k个基因的正常的表达谱数据
    colN<-dim(control_exp)[2]##正常样本的个数
    colC<-dim(case_exp)[2]###case样本个数
    N_tmp=matrix(rep(Nnorm,Lgene),ncol=colN,byrow=T)-control_exp##matrix(1:12, ncol=4, byrow=T)以四个为一行，构建矩阵
    Nloc_up=which((rowSums(N_tmp>0)/colN)>freq)
    Nloc_down=which((rowSums(N_tmp<0)/colN)>freq)
    reverse=matrix(0,colC,4)###colc是case样本个数,matrix产生一个colc行4列的矩阵，矩阵元素全为0
    reverse[,1]=rep(length(Nloc_up),colC)##fisher检验的norm里面上调的总数
    reverse[,2]=rep(length(Nloc_down),colC)###fisher检验的norm里面下调的总数
    Tcanc=case_exp[k,]	###case的tumor里面的第k个基因	
    if (length(Nloc_up)>0){
      N_tmp=matrix(rep(Tcanc,length(Nloc_up)),ncol=colC,byrow=T)-case_exp[Nloc_up,]###看norm里面》的稳定对在case里面的保持情况
      case_p=colSums(N_tmp<0)##基因k的每个case相对于norm反转的对c
      reverse[,3]=case_p
    }
    if (length(Nloc_down)>0){
      N_tmpp=matrix(rep(Tcanc,length(Nloc_down)),ncol=colC,byrow=T)-case_exp[Nloc_down,]
      case_pp=colSums(N_tmpp>0)###case相对于normal反转的对d
      reverse[,4]=case_pp; 
    }		
    GenePair_sig=NULL
    GenePair=rep(0,colC)
    GenePair[which(reverse[,3]>reverse[,4])]<--1###-1代表基因是下调的
    GenePair[which(reverse[,3]<reverse[,4])]<-1###1代表基因是上调的
    tmp=matrix(c(reverse[,1],reverse[,2],reverse[,1]-reverse[,3]+reverse[,4], reverse[,2]-reverse[,4]+reverse[,3]),ncol=4)
    GenePair_sig<-apply(tmp,1,function(x) fisher.test(matrix(x,ncol=2,byrow=T))$p.value)
    outlier_dir=rbind(outlier_dir,GenePair)
    outlier_pvalue=rbind(outlier_pvalue,GenePair_sig)
  }
  fdr<-apply(outlier_pvalue,2,function(x) p.adjust(x,method="fdr",length(x)))
  Methout<-cbind(gene,outlier_dir,outlier_pvalue,fdr);
  colnames(Methout)<-c(first_category_name[j],"direction","pvalue","fdr")
  result<-cbind(result,Methout)
}  

m=1
position<-result[order(result[,m*4])[1:40],m*4-3]
for(m in 2:n){
  position1<-result[order(result[,m*4])[1:40],m*4-3]
  position<-union(position,position1)
}
feature_RankComp<-as.matrix(position)
colnames(feature_RankComp)<-"feature_gene_RankComp"
setwd("~/xjj/CUP")
write.table(feature_RankComp,file="TCGA_RankComp_40_feature_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')



#####  FC
result<-matrix(0,40,n)
for(j in 1:n){
  cat(first_category_name[j],"\n")
  case_exp<-as.matrix(apply(expp[,colnames(expp)%in%first_category_name[j]],1,median))
  control_exp<-as.matrix(apply(expp[,-(which(colnames(expp)==first_category_name[j]))],1,median))
  FC=case_exp/control_exp
  result[,j]<-gene[order(FC,decreasing =TRUE)[1:40]]
}
colnames(result)<-first_category_name

one<-result[,1]
for(k in 1:(ncol(result)-1)){
  k=k+1
  two<-result[,k]
  one<-union(one,two)
}
feature<-as.matrix(one)
colnames(feature)<-"feature_gene_FC"
setwd("~/xjj/CUP")
write.table(feature,file="TCGA_FC_40_feature_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')


#############     ICGC       ##############00000000000000000000000000#############
rm(list=ls())
setwd("~/xjj/CUP")
first_category_name = list.files("ICGC")  
dir = paste("~/xjj/CUP/ICGC/",first_category_name,sep="")  
n = length(dir)
library(foreach)
library("data.table")

result<-data.frame(matrix(0,n,4))
foreach(i=1:n,.combine='rbind')%do%{
  setwd(dir[i])
  clinical<-fread("Merge_RNAseq_clinical.txt",header=T,data.table=F)
  exp<-fread("Merge_RNAseq_Count.txt",header=T,data.table=F)
  exp1<-exp[,-1] #exp matrix
  geneid<-exp[,1] #exp geneid
  clincal_exp<-intersect(clinical[,1],colnames(exp1)) ## int samples
  int_clinical<-clinical[match(clincal_exp,clinical[,1]),]##all int
  exp2<-exp1[,match(int_clinical[int_clinical$specimen_type=="Primary tumour - solid tissue",1],colnames(exp1))]
  Primary<-cbind(geneid,exp2)
  colnames(Primary)<-c("geneid",colnames(exp2))
  write.table(Primary,file=paste(first_category_name[i],"_primary_exp.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
  
  cat(first_category_name[i],"\n")
  cat(nrow(int_clinical),"\n")
  cat(sum(int_clinical$specimen_type=="Primary tumour - solid tissue"),"\n")
  cat(sum(int_clinical$specimen_type!="Primary tumour - solid tissue"),"\n")
  
  result[i,1]<-first_category_name[i]
  result[i,2]<-nrow(int_clinical)  ## all case
  result[i,3]<-sum(int_clinical$specimen_type=="Primary tumour - solid tissue") #primary
  result[i,4]<-sum(int_clinical$specimen_type!="Primary tumour - solid tissue") #no primary
}  
colnames(result)<-c("cancer","the number of all case","primary","no primary")
setwd("~/xjj/CUP")
write.table(result,"ICGC_primary_no_primary_sample.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出



#################   GEO    #########################################
rm(list=ls())
#setwd("~/xjj/CUP/GEO")
#first_category_name1 = list.files("Lung") 
first_category_name<-c("GSE40367")
dir = paste("~/xjj/CUP/GEO/Liver/",first_category_name,sep="")  
n = length(dir)
library(foreach)
library("data.table")

foreach(i=1:n,.combine='rbind')%do%{
  setwd(dir[i])
  cat(first_category_name[i],"\n")
  exp<-read.delim("prob_exp1.txt",sep="\t",row.names=1)
  probe_gene<-read.table("GPL570.txt",sep="\t",header=T)	##经过预处理，删除单一探针对应多基因和探针无对应基因的行 ，保留下的列表第一列为探针，第二列为基因。
  #genelist<-factor(probe_gene[,2])##取探针基因表中的第二列改为因子型，因为因子型的水平只有单独一个
  exp<-t(exp)
  rownames(exp)<-gsub("X","",rownames(exp))
  exp<-exp[match(probe_gene[,1],rownames(exp)),]	##探针id和表达谱列名匹配，取只在列表中出现的基因，返回那个位置的所有列
  colnames(exp)<-sub("(_\\w+)?(\\.\\w+)+","",colnames(exp)) ##第一个括号内是去掉下划线后面的字符，第二个括号是去掉点后面的字符。
  exp[1:5,1:2]
  drop1=grep(pattern='///',probe_gene[,2])##对应到多个基因上的探针的行
  drop2=which(is.na(probe_gene[,2]))##没有对应到基因上的探针的行
  geneid_drop = c(drop1, drop2)##去掉这两行
  probe_gene1=probe_gene[-geneid_drop,]##合并去掉后的数据
  int_probe<-intersect(probe_gene1[,1],rownames(exp))##共有的探针
  probe_gene2<-probe_gene1[match(int_probe,probe_gene1[,1]),]
  genelist<-factor(probe_gene2[,2])
  exp2<-exp[match(int_probe,rownames(exp)),] ###在同样位置去掉探针对应的表达谱
  dim(exp2)
  exp3<-apply(exp2,2,function(x) tapply(x,genelist,mean))##多探针对应一个基因取均值
  expp<-cbind(rownames(exp3),exp3)
  colnames(expp)<-c("geneid",colnames(exp3))
  exppp<-na.omit(as.matrix(expp))
  exppp[1:5,1:2]
  dim(exppp)
  write.table(exppp,file=paste("exp",first_category_name[i],".txt",sep=""),col.names=T,sep="\t",row.names = F)
  geneid<-as.matrix(rownames(exppp))##geneid为表达谱的行名
  colnames(geneid)<-"geneid"
  write.table(geneid,file=paste("geneid",first_category_name[i],".txt",sep=""),sep="\t",col.names=T,row.names = F)
}

first_category_name<-c("GSE19249","GSE5851")

rm(list=ls())
library("data.table")
clinical<-fread("GSE22932_clinical_sample.txt",header=T,data.table=F)
exppp<-fread("expGSE40367.txt",header=T,data.table=F)
dim(exppp)

library(clusterProfiler)
gene=bitr(exppp[,1],fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
ge<-intersect(gene[,1],exppp[,1])
exppp<-exppp[match(ge,exppp[,1]),]
symbol<-gene[match(ge,gene[,1]),2]
exppp[,1]<-symbol

unique(clinical[,3])
a<-clinical[clinical[,3]%in%"Low-passage culture of visceral metastasis from cutaneous melanoma",1]
#a<-clinical[-which(clinical[,3]=="Frozen tissue of non tumoral colorectal mucosa"),2]
#primary<-exppp[,match(a,substring(colnames(exppp),1,9))]
primary<-as.data.frame(exppp[,match(a,colnames(exppp))])
#primary1<-cbind(rownames(primary),primary)
primary1<-cbind(exppp[,1],primary)
colnames(primary1)<-c("geneid",a)
head(primary1)[1:5,1:2]
write.table(primary1,file="GSE46141_breast_cancer_lung_metastasis.txt",sep="\t",col.names=T,row.names = F)

primary1<-exppp[,c(1,2)]
write.table(primary1,file="SRP003173_skin_normal.txt",sep="\t",col.names=T,row.names = F)
write.table(primary1,file="GSE33532_normal_lung.txt",sep="\t",col.names=T,row.names = F)
write.table(exppp,file="GSE116174_primary_liver.txt",sep="\t",col.names=T,row.names = F)



######################     METABRIC     ############
library("data.table")
setwd("~/xjj/CUP/GEO/Breast/METABRIC/brca_metabric")
exp<-fread("data_mrna_agilent_microarray_zscores_ref_all_samples.txt",header=T,data.table=F)
#exp<-fread("data_mrna_agilent_microarray.txt",header=T,data.table=F)
clinical<-fread("data_clinical_sample.txt",header=T,data.table=F)
unique(clinical[,10])
a<-clinical[clinical[,10]%in%"Primary",2]
b<-intersect(a,colnames(exp))
primary<-exp[,match(b,colnames(exp))]
write.table(cbind(exp[,1:2],primary),file="METABRIC_primary_breast_zscores.txt",sep="\t",col.names=T,row.names = F)

exp1<-fread("METABRIC_primary_breast_zscores.txt",header=T,data.table=F)




####################  FFPE  ####################
#############   TCGA删除0值   ########
rm(list=ls())
setwd("~/xjj/CUP")
first_category_name = list.files("TCGA")  
dir = paste("~/xjj/CUP/TCGA/",first_category_name,sep="")  
n = length(dir)
library(foreach)
library("data.table")
i=1
setwd(dir[i])
#clinical<-fread("Merge_clinical.txt",header=T,data.table=F)
primary<-fread(paste(first_category_name[i],"_primary_exp_FFPE.txt",sep=""),header=T,data.table=F)
clinical<-fread(paste("TCGA-",first_category_name[i],".GDC_phenotype.tsv",sep=""),header=T,data.table=F)
exp1<-primary[,-1] #exp matrix
geneid<-primary[,1] #exp geneid
if(ncol(exp1)>0){
sampletype<-substring(clinical[,1],16)
clinical_want<-clinical[sampletype%in%"B",]
clincal_exp<-intersect(substring(clinical_want[,1],1,15),colnames(primary))##the exp of fresh sample
int_exp<-exp1[,match(clincal_exp,colnames(exp1))]######the exp of fresh sample
clinical1<-clinical_want[match(clincal_exp,substring(clinical_want[,1],1,15)),]
sample1<-matrix(rep(first_category_name[i],nrow(clinical1)),ncol=1)
result1<-cbind(sample1,clinical1$submitter_id.samples,clinical1$batch_number)
}
if(ncol(exp1)==0){
  result11<-matrix(0,1,3)
  result11[1,1]<-first_category_name[i]
  result1<-result11
}


result<-result1
foreach(i=1:n,.combine='rbind')%do%{
  i=i+1
  setwd(dir[i])
  primary<-fread(paste(first_category_name[i],"_primary_exp_FFPE.txt",sep=""),header=T,data.table=F)
  clinical<-fread(paste("TCGA-",first_category_name[i],".GDC_phenotype.tsv",sep=""),header=T,data.table=F)
  exp1<-primary[,-1] #exp matrix
  geneid<-primary[,1] #exp geneid
  if(ncol(exp1)>0){
  sampletype<-substring(clinical[,1],16)
  clinical_want<-clinical[sampletype%in%"B",]
  clincal_exp<-intersect(substring(clinical_want[,1],1,15),colnames(primary))##the exp of fresh sample
  int_exp<-exp1[,match(clincal_exp,colnames(exp1))]######the exp of fresh sample
  clinical2<-clinical_want[match(clincal_exp,substring(clinical_want[,1],1,15)),]
  sample2<-matrix(rep(first_category_name[i],nrow(clinical2)),ncol=1)
  result2<-cbind(sample2,clinical2$submitter_id.samples,clinical2$batch_number)
  result<-rbind(result,result2)
  }
  if(ncol(exp1)==0){
    result11<-matrix(0,1,3)
    result11[1,1]<-first_category_name[i]
    result2<-result11
    result<-rbind(result,result2)
  }
}


colnames(result)<-c("cancer type","sample","batch_number")
setwd("~/xjj/CUP")
write.table(result,"TCGA_FFPE_batch_number.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

############     TCGA所有的表达谱merge   #######################
rm(list=ls())
setwd("~/xjj/CUP")
first_category_name = list.files("TCGA")  
dir = paste("~/xjj/CUP/TCGA/",first_category_name,sep="")  
n = length(dir)
library(foreach)
library("data.table")
i=1
setwd(dir[i])
primary<-fread(paste(first_category_name[i],"_primary_exp_delete50_FFPE.txt",sep=""),header=T,data.table=F)
geneid<-primary[,1]

gene<-geneid
for(i in 1:(n-1)){
  i=i+1
  setwd(dir[i])
  primary1<-fread(paste(first_category_name[i],"_primary_exp_delete50_FFPE.txt",sep=""),header=T,data.table=F)
  geneid1<-primary1[,1]
  gene<-intersect(gene,geneid1)
}
int_gene<-data.frame(gene)
colnames(int_gene)<-"TCGA_all_cancer_geneid"
setwd("~/xjj/CUP")
write.table(int_gene,"TCGA_all_cancer_geneid_FFPE.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

#############    TCGA exp    ###########################
rm(list=ls())
setwd("~/xjj/CUP")
int_gene<-as.matrix(read.delim("TCGA_all_cancer_geneid.txt",sep = '\t',header = T))##数据输出
setwd("~/xjj/CUP")
first_category_name = list.files("TCGA")  
dir = paste("~/xjj/CUP/TCGA/",first_category_name,sep="")  
n = length(dir)
library(foreach)
library("data.table")

foreach(i=1:n,.combine='rbind')%do%{
  setwd(dir[i])
  primary1<-fread(paste(first_category_name[i],"_primary_exp_FFPE.txt",sep=""),header=T,data.table=F)
  exp<-primary1[match(int_gene,primary1[,1]),]
  write.table(exp,file=paste(first_category_name[i],"_primary_exp_delete50_consistent_gene_FFPE.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
}


foreach(i=1:n,.combine='rbind')%do%{
  setwd(dir[i])
  primary1<-fread(paste(first_category_name[i],"_metastatic_exp.txt",sep=""),header=T,data.table=F)
  exp<-primary1[match(int_gene,primary1[,1]),]
  write.table(exp,file=paste(first_category_name[i],"_metastatic_exp_delete50_consistent_gene.txt",sep=""),sep = '\t',
              col.names = T,row.names = F,quote = FALSE,na='')
}

######################  biomarker 预处理   #########################
setwd("~/xjj/CUP/biomarker")
library("data.table")
#genetics<-fread("all_sequence_variants.tsv",header=T,data.table=F)
genetics<-fread("all_diagnostic_proteins.tsv",header=F,data.table=F)
head(genetics)
disease1<-unique(genetics[,4])
disease<-disease1[order(disease1)]
disease_gene<-data.frame()
for(i in 1:length(disease)){
  gene<-unique(genetics[genetics[,4]%in%disease[i],2])
  count<-length(gene)
  gene1<-paste0(gene,collapse=",")
  disease_gene[i,1]<-disease[i]
  disease_gene[i,2]<-count
  disease_gene[i,3]<-gene1
}
colnames(disease_gene)<-c("Disease","Number of protein","proteins")

setwd("~/xjj/CUP/biomarker")
write.table(disease_gene,file="disease_only_one_biomarker.txt",sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')

## each biomarker match disease
genetics1<-fread("all_diagnostic_proteins.tsv",header=F,data.table=F)
genetics<-genetics1[-which(genetics1[,4]=="Normal"),]
head(genetics)
biomarker<-unique(genetics[,2])
biomarker_disease<-data.frame()
for(i in 1:length(biomarker)){
  gene<-unique(genetics[genetics[,2]%in%biomarker[i],4])
  count<-length(gene)
  gene1<-paste0(gene,collapse=",")
  biomarker_disease[i,1]<-biomarker[i]
  biomarker_disease[i,2]<-count
  biomarker_disease[i,3]<-gene1
}
colnames(biomarker_disease)<-c("biomarker","Number of disease","disease")
biomarker_disease_gene11<-biomarker_disease[order(biomarker_disease[,2]),]
biomarker_disease_gene12<-biomarker_disease_gene11[biomarker_disease_gene11[,2]%in%"1",]
setwd("~/xjj/CUP/biomarker")
write.table(biomarker_disease_gene12,file="each_protein_biomarker_match_one_disease.txt",sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')

bcd<-unique(biomarker_disease_gene12[,3])
biomarker_disease56<-data.frame()
for(i in 1:length(bcd)){
  gene<-unique(biomarker_disease_gene12[biomarker_disease_gene12[,3]%in%bcd[i],1])
  count<-length(gene)
  gene1<-paste0(gene,collapse=",")
  biomarker_disease56[i,1]<-bcd[i]
  biomarker_disease56[i,2]<-count
  biomarker_disease56[i,3]<-gene1
}
colnames(biomarker_disease56)<-c("disease","Number of biomarker","biomarker")
setwd("~/xjj/CUP/biomarker")
write.table(biomarker_disease56,file="each_disease_match_special_biomarker.txt",sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')

####################################################################
cancer_type<-read.table("cancer_type.txt",header=F,sep="\t")
cancer_result<-disease_gene[disease_gene[,1]%in%cancer_type[,1],]

result2<-NULL
for(i in 1:nrow(cancer_type)){
  gene<-matrix(unique(genetics[genetics[,4]%in%cancer_type[i,1],2]),ncol=1)
  type<-matrix(rep(cancer_type[i,1],length(gene)),ncol=1)
  result1<-cbind(type,gene)
  result2<-rbind(result2,result1)
}
colnames(disease_gene)<-c("Disease","Number of protein","proteins")


write.table(cancer_result,file="cancer_type_Diagnostic_Protein_markers_disease_gene.txt",sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')







