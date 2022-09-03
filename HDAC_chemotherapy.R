#######  根据药敏数据将样本分为两类
################################
#############   HDAC+chemotherapy  #######################
rm(list=ls())
setwd("~/xjj/drug")
ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,!(colnames(ic501)%in%c("PC.100","PC.34"))]

drug_info<-read.csv("compiled_drug_info.csv",header = TRUE,sep = ",")
drug_result<-drug_info[match(rownames(ic50),drug_info$drug_id),]
HDAC1<-drug_result[drug_result$target=="HDAC",]
aaa<-c("S1848","S2759","S1194","S1047")
bbb<-c("5-FU","GEM","IRI","OXA","PAC")
ccc<-c(aaa,bbb)
HDAC2<-drug_result[drug_result$drug_id%in%ccc,]
HDAC3<-rbind(HDAC1,HDAC2)
HDAC<-HDAC3[order(HDAC3[,1],decreasing = T),]
HDAC_ic50<-ic50[match(as.character(HDAC[,1]),rownames(ic50)),]
library(pheatmap)
a=cor(t(HDAC_ic50))
result=pheatmap(a,scale = "none",main = "The correlation coefficients between HDACi and chemotherapy' AUC",show_rownames=T,show_colnames=T,
                clustering_distance_rows = "correlation",
                clustering_distance_cols = "correlation",clustering_method = "complete",cutree_rows=2,cutree_cols=2)

##########   20种HDAC+5种化疗药
bb<-t(apply(HDAC_ic50,1,function(x) scale(x)))
colnames(bb)<-colnames(HDAC_ic50)
norma_result<-pheatmap(bb,scale = "none",main = "The AUC of HDACi and chemotherapy are normalized",show_rownames=T,show_colnames=T,
                       clustering_distance_rows = "correlation",
                       clustering_distance_cols = "euclidean",clustering_method = "ward.D2")


clusters_filter<-cutree(norma_result$tree_col,k=2)
sensitivity=colnames(HDAC_ic50[,clusters_filter==1])
resistance=colnames(HDAC_ic50[,clusters_filter==2])
#resistance1=c("PC.104","PC.105","PC.5","PC.56")
#medianl<-resistance[-match(resistance1,resistance)]


sensitivity_AUC<-HDAC_ic50[,match(sensitivity,colnames(HDAC_ic50))]
resistance_AUC<-HDAC_ic50[,match(resistance,colnames(HDAC_ic50))]

geneid<-rownames(HDAC_ic50)
tresult=matrix(0,length(geneid),4)
for (i in 1:length(geneid)){
  ttest=t.test(sensitivity_AUC[i,],resistance_AUC[i,])
  tresult[i,1:3]=c(geneid[i],ttest$statistic,ttest$p.value)
}
tresult[,4]=p.adjust(tresult[,3],method="BH")
colnames(tresult)<-c("geneid","statistic","p.value","FDR")
which(tresult[,4]<0.05)

tresult=matrix(0,length(geneid),4)
for (i in 1:length(geneid)){
  ttest=wilcox.test(as.numeric(sensitivity_AUC[i,]),as.numeric(resistance_AUC[i,]),alternative ="two.sided")
  tresult[i,1:3]=c(geneid[i],ttest$statistic,ttest$p.value)
}
tresult[,4]=p.adjust(tresult[,3],method="BH")


rm(list=ls())
setwd("~/xjj/drug")
ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,!(colnames(ic501)%in%c("PC.100","PC.34"))]

drug_info<-read.csv("compiled_drug_info.csv",header = TRUE,sep = ",")
drug_result<-drug_info[match(rownames(ic50),drug_info$drug_id),]
HDAC1<-drug_result[drug_result$target=="HDAC",]
aaa<-c("S2759","S1194","S1047")
bbb<-c("GEM","IRI","PAC")

ccc<-c(aaa,bbb)
HDAC2<-drug_result[drug_result$drug_id%in%ccc,]
HDAC3<-rbind(HDAC1,HDAC2)
HDAC<-HDAC3[order(HDAC3[,1],decreasing = T),]
HDAC_ic50<-ic50[match(as.character(HDAC[,1]),rownames(ic50)),]
library(pheatmap)
a=cor(t(HDAC_ic50))
result=pheatmap(a,scale = "none",main = "The correlation coefficients between HDACi and chemotherapy' AUC",show_rownames=T,show_colnames=T,
                clustering_distance_rows = "correlation",
                clustering_distance_cols = "correlation",clustering_method = "complete",cutree_rows=2,cutree_cols=2)

##########   20种HDAC+5种化疗药
bb<-t(apply(HDAC_ic50,1,function(x) scale(x)))
colnames(bb)<-colnames(HDAC_ic50)
norma_result<-pheatmap(bb,scale = "none",main = "The AUC of HDACi and chemotherapy are normalized",show_rownames=T,show_colnames=T,
                       clustering_distance_rows = "correlation",
                       clustering_distance_cols = "euclidean",clustering_method = "ward.D2")


clusters_filter<-cutree(norma_result$tree_col,k=2)
sensitivity=colnames(HDAC_ic50[,clusters_filter==1])
resistance=colnames(HDAC_ic50[,clusters_filter==2])



##############
sensitivity<-c("PC.101","PC.109","PC.115","PC.116","PC.117","PC.13","PC.130","PC.136","PC.139","PC.16","PC.18","PC.2","PC.27","PC.64","PC.78","PC.8","PC.81","PC.97","PC.98","PC.L" )
resistance<-c("PC.102","PC.104","PC.105","PC.111","PC.112","PC.119","PC.121","PC.134","PC.135","PC.14","PC.22","PC.40","PC.5","PC.52","PC.56","PC.G","PC.I" )

########  clinical
#new clinical
rm(list=ls())
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
library("data.table")
clinical<-fread("sample_information.txt",header=T,data.table=F)
sensitivity<-c("PC.101","PC.109","PC.115","PC.116","PC.117","PC.13","PC.130","PC.136","PC.139","PC.16","PC.18","PC.2","PC.27","PC.64","PC.78","PC.8","PC.81","PC.97","PC.98","PC.L" )
resistance<-c("PC.102","PC.104","PC.105","PC.111","PC.112","PC.119","PC.121","PC.134","PC.135","PC.14","PC.22","PC.40","PC.5","PC.52","PC.56","PC.G","PC.I" )
resistance_new<-as.character(sample_id[match(resistance,sample_id[,2]),1])
sensitivity_new<-as.character(sample_id[match(sensitivity,sample_id[,2]),1])
sen<-clinical[match(sensitivity_new,clinical[,1]),]
res<-clinical[match(resistance_new,clinical[,1]),]
names<-matrix(c(rep("sensitivity",nrow(sen)),rep("resistance",nrow(res))),ncol=1)
clinical_table<-cbind(names,rbind(sen,res))
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5")
write.table(clinical_table,"clinical_table.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

##old
rm(list=ls())
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
library("data.table")
clinical<-fread("HDAC_IC50_information.csv",header=T,data.table=F)
  clinical_result<-data.frame(matrix(0,7,3))
  sensitivity<-c("PC.101","PC.109","PC.115","PC.116","PC.117","PC.13","PC.130","PC.136","PC.139","PC.16","PC.18","PC.2","PC.27","PC.64","PC.78","PC.8","PC.81","PC.97","PC.98","PC.L" )
  resistance<-c("PC.102","PC.104","PC.105","PC.111","PC.112","PC.119","PC.121","PC.134","PC.135","PC.14","PC.22","PC.40","PC.5","PC.52","PC.56","PC.G","PC.I" )
  resistance_new<-as.character(sample_id[match(resistance,sample_id[,2]),1])
  sensitivity_new<-as.character(sample_id[match(sensitivity,sample_id[,2]),1])
  sen<-clinical[match(sensitivity_new,clinical$sample_new),]
  res<-clinical[match(resistance_new,clinical$sample_new),]
  clinical_result[,1]<-as.matrix(c("Age","Sex","Pathological type","location","Degree of differentiation","Neoadjuvant therapy","Tumor stage"),ncol=1)
  clinical_result[1,2]<-paste("median",median(sen$age),"Range",min(sen$age),"-",max(sen$age),sep=",")
  clinical_result[1,3]<-paste("median",median(res$age),"Range",min(res$age),"-",max(res$age),sep=",")
  clinical_result[2,2]<-paste("female",sum(sen$sex==1),"Ratio",round(sum(sen$sex==1)/length(sen$sex),digits = 2),"male",sum(sen$sex==2),"Ratio",round(sum(sen$sex==2)/length(sen$sex),digits = 2),sep=",")
  clinical_result[2,3]<-paste("female",sum(res$sex==1),"Ratio",round(sum(res$sex==1)/length(res$sex),digits = 2),"male",sum(res$sex==2),"Ratio",round(sum(res$sex==2)/length(res$sex),digits = 2),sep=",")
  clinical_result[3,2]<-paste("PDAC",length(sen$sex),"100%",sep=",")
  clinical_result[3,3]<-paste("PDAC",length(res$sex),"100%",sep=",")
  clinical_result[4,2]<-paste("Pancreatic head",sum(sen$location==1),"Pancreatic tail",sum(sen$location==2),"liver metastasis",sum(sen$location==3),sep=",")
  clinical_result[4,3]<-paste("Pancreatic head",sum(res$location==1),"Pancreatic tail",sum(res$location==2),"liver metastasis",sum(res$location==3),sep=",")
  clinical_result[5,2]<-paste("Low differentiation",sum(sen$`Degree of differentiation`==1),"Medium differentiation",sum(sen$`Degree of differentiation`==2),"none",sum(sen$`Degree of differentiation`==0),sep=",")
  clinical_result[5,3]<-paste("Low differentiation",sum(res$`Degree of differentiation`==1),"Medium differentiation",sum(res$`Degree of differentiation`==2),"none",sum(res$`Degree of differentiation`==0),sep=",")
  clinical_result[6,2]<-paste("Yes",sum(sen$`Neoadjuvant therapy`==1),"No",sum(sen$`Neoadjuvant therapy`==0),sep=",")
  clinical_result[6,3]<-paste("Yes",sum(res$`Neoadjuvant therapy`==1),"No",sum(res$`Neoadjuvant therapy`==0),sep=",")
  clinical_result[7,2]<-paste("IB",sum(sen$`Tumor stage`=="IB"),"IIA",sum(sen$`Tumor stage`=="IIA"),"IIB",sum(sen$`Tumor stage`=="IIB"),"III",sum(sen$`Tumor stage`=="III"),"IV",sum(sen$`Tumor stage`=="IV"),sep=",")
  clinical_result[7,3]<-paste("IB",sum(res$`Tumor stage`=="IB"),"IIA",sum(res$`Tumor stage`=="IIA"),"IIB",sum(res$`Tumor stage`=="IIB"),"III",sum(res$`Tumor stage`=="III"),"IV",sum(res$`Tumor stage`=="IV"),sep=",")
  colnames(clinical_result)<-c("Characteristic","Sensitive","Resistance")
  setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5")
  write.table(clinical_result,"clinical_information.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

#############   mutation landscape Fisher's exact test
#new
rm(list=ls())
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
sensitivity<-c("PC.101","PC.109","PC.115","PC.116","PC.117","PC.13","PC.130","PC.136","PC.139","PC.16","PC.18","PC.2","PC.27","PC.64","PC.78","PC.8","PC.81","PC.97","PC.98","PC.L" )
resistance<-c("PC.102","PC.104","PC.105","PC.111","PC.112","PC.119","PC.121","PC.134","PC.135","PC.14","PC.22","PC.40","PC.5","PC.52","PC.56","PC.G","PC.I" )
resistance_new<-as.character(sample_id[match(resistance,sample_id[,2]),1])
sensitivity_new<-as.character(sample_id[match(sensitivity,sample_id[,2]),1])
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/3_mutation")
mut<-read.table("nocoding_intergenic_binary.txt",header = TRUE,sep = "\t",row.names=1)
colnames(mut)<-gsub("\\.","-",gsub("CAS.","",colnames(mut)))
res_int<-intersect(resistance_new,colnames(mut))
resistance_mut<-mut[,match(res_int,colnames(mut))]
sen_int<-intersect(sensitivity_new,colnames(mut))
sensitivity_mut<-mut[,match(sen_int,colnames(mut))]
muty<-cbind(resistance_mut,sensitivity_mut)

geneid<-rownames(mut)
#geneid[which(delete>=10)]
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
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/3_mutation")
write.table(data_result,"nocoding_intergenic_fisher_all.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
#plot
rm(list=ls())
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
sensitivity<-c("PC.101","PC.109","PC.115","PC.116","PC.117","PC.13","PC.130","PC.136","PC.139","PC.16","PC.18","PC.2","PC.27","PC.64","PC.78","PC.8","PC.81","PC.97","PC.98","PC.L" )
resistance<-c("PC.102","PC.104","PC.105","PC.111","PC.112","PC.119","PC.121","PC.134","PC.135","PC.14","PC.22","PC.40","PC.5","PC.52","PC.56","PC.G","PC.I" )
resistance_new<-as.character(sample_id[match(resistance,sample_id[,2]),1])
sensitivity_new<-as.character(sample_id[match(sensitivity,sample_id[,2]),1])
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/3_mutation")
mut1<-read.table("nocoding_UTR_binary.txt",header = TRUE,sep = "\t",row.names=1)
mut2<-read.table("nocoding_up-downstream_binary.txt",header = TRUE,sep = "\t",row.names=1)
mut3<-read.table("nocoding_ncRNA_binary.txt",header = TRUE,sep = "\t",row.names=1)
mut4<-read.table("nocoding_intronic_binary.txt",header = TRUE,sep = "\t",row.names=1)
mut5<-read.table("nocoding_intergenic_binary.txt",header = TRUE,sep = "\t",row.names=1)
setwd("~/xjj/WGS_CNV/mutation")
library("data.table")
mut<-fread("Organoid_mut_binary.txt",header=T,data.table=F)
mut6<-mut[,-1]
colnames(mut1)<-gsub("\\.","-",gsub("CAS.","",colnames(mut1)))
colnames(mut2)<-gsub("\\.","-",gsub("CAS.","",colnames(mut2)))
colnames(mut3)<-gsub("\\.","-",gsub("CAS.","",colnames(mut3)))
colnames(mut4)<-gsub("\\.","-",gsub("CAS.","",colnames(mut4)))
colnames(mut5)<-gsub("\\.","-",gsub("CAS.","",colnames(mut5)))
colname<-intersect(intersect(intersect(intersect(intersect(colnames(mut1),colnames(mut2)),colnames(mut3)),colnames(mut4)),colnames(mut5)),colnames(mut6))
res_int<-intersect(resistance_new,colname)
sen_int<-intersect(sensitivity_new,colname)
all_auc<-c(res_int,sen_int)
all_mut1<-mut1[,match(all_auc,colnames(mut1))]
all_mut2<-mut2[,match(all_auc,colnames(mut2))]
all_mut3<-mut3[,match(all_auc,colnames(mut3))]
all_mut4<-mut4[,match(all_auc,colnames(mut4))]
all_mut5<-mut5[,match(all_auc,colnames(mut5))]
all_mut6<-mut6[,match(all_auc,colnames(mut6))]

name1<-rownames(all_mut1)[order(apply(all_mut1,1,sum),decreasing = TRUE)]
name2<-rownames(all_mut2)[order(apply(all_mut2,1,sum),decreasing = TRUE)]
name3<-rownames(all_mut3)[order(apply(all_mut3,1,sum),decreasing = TRUE)]
name4<-rownames(all_mut4)[order(apply(all_mut4,1,sum),decreasing = TRUE)]
name5<-rownames(all_mut5)[order(apply(all_mut5,1,sum),decreasing = TRUE)]
name6<-mut[order(apply(all_mut6,1,sum),decreasing = TRUE),1]
#hotplot


#old
  rm(list=ls())
  library("data.table")
  library(dplyr)
  library(stringr)
  library(foreach)
  setwd("~/xjj/WGS_CNV/mutation")
  library("data.table")
  mut<-fread("Organoid_mut_binary.txt",header=T,data.table=F)
  mut1<-mut[,-1]
  setwd("~/xjj/drug")
  sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
  
  mutation_result<-data.frame(matrix(0,1,3))
    sensitivity<-c("PC.101","PC.109","PC.115","PC.116","PC.117","PC.13","PC.130","PC.136","PC.139","PC.16","PC.18","PC.2","PC.27","PC.64","PC.78","PC.8","PC.81","PC.97","PC.98","PC.L" )
    resistance<-c("PC.102","PC.104","PC.105","PC.111","PC.112","PC.119","PC.121","PC.134","PC.135","PC.14","PC.22","PC.40","PC.5","PC.52","PC.56","PC.G","PC.I" )
    resistance_new<-as.character(sample_id[match(resistance,sample_id[,2]),1])
    sensitivity_new<-as.character(sample_id[match(sensitivity,sample_id[,2]),1])
    res_int<-intersect(resistance_new,colnames(mut))
    resistance_mut<-mut[,match(res_int,colnames(mut))]
    sen_int<-intersect(sensitivity_new,colnames(mut))
    sensitivity_mut<-mut[,match(sen_int,colnames(mut))]
    muty<-cbind(resistance_mut,sensitivity_mut)
    #delete<-apply(muty,1,function(x) sum(x==1))
    #length(delete[which(delete>=10)])
    geneid<-mut[,1]
    #geneid[which(delete>=10)]
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

colnames(mutation_result)<-c("drug","The number of gene mutation","mutation gene")
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/3_mutation")
write.table(data,"HDAC_chemotherapy_mutation_all.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

###TCGA-mutation

  
##################   CNV  
##TCGA-CNV
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/0_CNV")
TCGA_CNV<-read.table("Merge_GeneLevelCopyNumber.txt",header = TRUE,sep = "\t")
ENSG_TCGA<-substring(TCGA_CNV$Gene.Symbol,1,15)
setwd("~/xjj/CUP/gtf")
ENSG<-read.table("ENSG.txt",sep="\t",header=T)
id<-read.table("2id.txt",sep="\t",header=T)
int_gene<-intersect(ENSG[,2],id[,2])
id1<-id[match(int_gene,id[,2]),]
ENSG1<-ENSG[match(int_gene,ENSG[,2]),]
ENSG_id_symbol<-cbind(ENSG1,id1)
int_ENSG<-intersect(ENSG_TCGA,ENSG_id_symbol[,1])
TCGA_CNV1<-TCGA_CNV[match(int_ENSG,ENSG_TCGA),]
symbol<-ENSG_id_symbol[match(int_ENSG,ENSG_id_symbol[,1]),3]
TCGA_CNV_symbol<-cbind(symbol,TCGA_CNV1)
#接下去直接做fisher
TCGA_CNV_symbol<-TCGA_CNV_symbol[,-c(2:4)]
CNV<-apply(TCGA_CNV_symbol[,-1],1,sum)
CNV_gene<-TCGA_CNV_symbol[order(CNV,decreasing = TRUE),1]


rm(list=ls())
  library("data.table")
  library(dplyr)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(stringr)
  library(foreach)
  setwd("~/xjj/WGS_CNV/CNV_FACETS")
  mut<-fread("CNV_FACETS_median_CNA_matrix.txt",header=T,data.table=F)# -2,2 在
  #mut<-fread("CNV_FACETS_median_CNA_matrix_del_loss_gain_amp.txt",header=T,data.table=F)##-2,-1,1,2都在
  setwd("~/xjj/drug")
  sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
  
  CNV_result<-data.frame(matrix(0,1,9))
  sensitivity<-c("PC.101","PC.109","PC.115","PC.116","PC.117","PC.13","PC.130","PC.136","PC.139","PC.16","PC.18","PC.2","PC.27","PC.64","PC.78","PC.8","PC.81","PC.97","PC.98","PC.L" )
  resistance<-c("PC.102","PC.104","PC.105","PC.111","PC.112","PC.119","PC.121","PC.134","PC.135","PC.14","PC.22","PC.40","PC.5","PC.52","PC.56","PC.G","PC.I" )
  resistance_new<-as.character(sample_id[match(resistance,sample_id[,2]),1])
  sensitivity_new<-as.character(sample_id[match(sensitivity,sample_id[,2]),1])
    res_int<-intersect(resistance_new,colnames(mut))
    resistance_mut<-mut[,match(res_int,colnames(mut))]
    sen_int<-intersect(sensitivity_new,colnames(mut))
    sensitivity_mut<-mut[,match(sen_int,colnames(mut))]
    geneid<-mut[,1]
    data_amp<-data.frame()
    for(i in 1:length(geneid)){
      a<-sum(sensitivity_mut[i,]==1)
      c<-sum(sensitivity_mut[i,]!=1)
      b<-sum(resistance_mut[i,]==1)
      d<-sum(resistance_mut[i,]!=1)
      tmp<-matrix(c(a, c,b, d),
                  nrow = 2,
                  dimnames = list(Truth = c("amplification", "wild"),
                                  Guess = c("sensitivity", "resistance")))
      P<-fisher.test(tmp,alternative = "two.sided",conf.level = 0.95)$p.value
      data_amp[i,1]<-geneid[i]
      data_amp[i,2]<-P
    }
    data_amp[,3]<-p.adjust(data_amp[,2],method="BH")
    amp_gene<-data_amp[which(data_amp[,2]<0.1),]  
    colnames(amp_gene)<-c("geneid","Pvalue","FDR")
    cat(paste("The number of amp gene is ",nrow(amp_gene),sep=""),"\n")
    cat(paste(amp_gene[,1],sep=""),"\n")
    if(nrow(amp_gene)>0){
      CNV_result[1,2]<-nrow(amp_gene)
      CNV_result[1,3]<-paste0(amp_gene[,1],collapse =",")
    }
    if(nrow(amp_gene)==0){
      CNV_result[1,3]<-0
    }
    
    data_censored<-data.frame()
    for(j in 1:length(geneid)){
      a<-sum(sensitivity_mut[j,]==(-1))
      c<-sum(sensitivity_mut[j,]!=(-1))
      b<-sum(resistance_mut[j,]==(-1))
      d<-sum(resistance_mut[j,]!=(-1))
      tmp<-matrix(c(a, c,b, d),
                  nrow = 2,
                  dimnames = list(Truth = c("censored", "wild"),
                                  Guess = c("sensitivity", "resistance")))
      P<-fisher.test(tmp,alternative = "two.sided",conf.level = 0.95)$p.value
      data_censored[j,1]<-geneid[i]
      data_censored[j,2]<-P
    }
    data_censored[,3]<-p.adjust(data_censored[,2],method="BH")
    censored_gene<-data_censored[which(data_censored[,2]<0.1),]  ###控制P值
    colnames(censored_gene)<-c("geneid","Pvalue","FDR")
    cat(paste("The number of censored gene is ",nrow(censored_gene),sep=""),"\n")
    cat(paste(censored_gene[,1],sep=""),"\n")
    CNV_result[1,4]<-nrow(censored_gene)
    if(nrow(censored_gene)>0){
      CNV_result[1,5]<-paste0(censored_gene[,1],collapse =",")
    }
    if(nrow(censored_gene)==0){
      CNV_result[1,5]<-0
    }
    
    int_CNV<-intersect(censored_gene[,1],amp_gene[,1])
    censored_gene1<-censored_gene#去掉既是扩增也是删失的基因
    amp_gene1<-amp_gene
    cat(paste("The number of amp int censored is ",length(int_CNV),sep=""),"\n")
    cat(paste("The number of final amp is ",nrow(amp_gene1),sep=""),"\n")
    cat(paste("The number of final censored is ",nrow(censored_gene1),sep=""),"\n")
    CNV_result[1,6]<-length(int_CNV)
    if(length(int_CNV)>0){
      CNV_result[1,7]<-paste0(int_CNV,collapse =",")
    }
    if(length(int_CNV)==0){
      CNV_result[1,7]<-0
    }
    
    if (nrow(amp_gene1)>0){
      setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/0_CNV/fisher")
      write.table(amp_gene1,"sensitivity-resistance-amp-CNV-gene-p0.1.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
      
      genea=bitr(amp_gene1[,1],fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
      kegga<-enrichKEGG(genea[,2],organism = 'hsa',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'none',
                        minGSSize = 10,maxGSSize = 500,qvalueCutoff = 1,use_internal_data = FALSE)
      ka<-barplot(kegga,showCategory=10,drop=T,x = "GeneRatio",color = "pvalue")
      eea<-kegga@result
      eea1<-eea[which(eea$pvalue<0.01),]
      CNV_result[1,8]<-nrow(eea1)
      setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/0_CNV/fisher")
      write.table(eea1,"kegg_CNV_amp_p0.1_pathway_P0.1.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
      ggsave(paste("kegg_CNV_amp_pathway_P0.1_DEGs_P0.1.pdf",sep=""),ka,width = 8, height = 6)
      
    }
    if (nrow(censored_gene1)>0){
      setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/0_CNV/fisher")
      write.table(censored_gene1,"sensitivity-resistance-censored-CNV-gene-p0.05.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
      genec=bitr(censored_gene1[,1],fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
      keggc<-enrichKEGG(genec[,2],organism = 'hsa',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'none',
                        minGSSize = 10,maxGSSize = 500,qvalueCutoff = 1,use_internal_data = FALSE)
      kc<-barplot(keggc,showCategory=10,drop=T,x = "GeneRatio",color = "pvalue")
      eec<-keggc@result
      eec1<-eec[which(eec$pvalue<0.05),]
      CNV_result[1,9]<-nrow(eec1)
      setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/0_CNV/fisher")
      write.table(eec1,"kegg_CNV_censored_P0.05_pathway_P0.05.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
      ggsave(paste("kegg_CNV_censored_pathway_P0.05_DEGs_P0.05.pdf",sep=""),kc,width = 8, height = 4)
}
    
  colnames(CNV_result)<-c("drug","Number of amp","amp_name","number of delete","delete_name","Number of int gene","int_gene_name","amp pathway","delete pathway")
  setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/0_CNV")
  write.table(CNV_result,"CNV_result_DEG_P0.1.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  
  
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
  setwd("~/xjj/drug/drug_result/drug64_AUC")
  tresult_want<-read.table("drug64_sample.txt",sep="\t",header=T)
  setwd("~/xjj/drug")
  sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
  
  cnv_t_result<-data.frame(matrix(0,nrow(tresult_want),3))
    sensitivity<-c("PC.101","PC.109","PC.115","PC.116","PC.117","PC.13","PC.130","PC.136","PC.139","PC.16","PC.18","PC.2","PC.27","PC.64","PC.78","PC.8","PC.81","PC.97","PC.98","PC.L" )
    resistance<-c("PC.102","PC.104","PC.105","PC.111","PC.112","PC.119","PC.121","PC.134","PC.135","PC.14","PC.22","PC.40","PC.5","PC.52","PC.56","PC.G","PC.I" )
    resistance_new<-as.character(sample_id[match(resistance,sample_id[,2]),1])
    sensitivity_new<-as.character(sample_id[match(sensitivity,sample_id[,2]),1])
    res_int<-intersect(resistance_new,colnames(mut))
    resistance_mut<-mut[,match(res_int,colnames(mut))]
    sen_int<-intersect(sensitivity_new,colnames(mut))
    sensitivity_mut<-mut[,match(sen_int,colnames(mut))]
    for(k in 1:ncol(sensitivity_mut)){
      sensitivity_mut[which(sensitivity_mut[,k]!=1),k]<-0
    }
    for(l in 1:ncol(resistance_mut)){
      resistance_mut[which(resistance_mut[,l]!=1),l]<-0
    }
    #res<-apply(resistance_mut,2,sum)
    #sen<-apply(sensitivity_mut,2,sum)
    #ttest=t.test(sen,res)
    geneid<-mut[,1]
    tresult=matrix(0,length(geneid),4)
    for (i in 1:length(geneid)){
      ttest=t.test(sensitivity_mut[i,],resistance_mut[i,])
      tresult[i,1:3]=c(geneid[i],ttest$statistic,ttest$p.value)
    }
    tresult[,4]=p.adjust(tresult[,3],method="BH")
    colnames(tresult)<-c("geneid","statistic","p.value","FDR")
    a1<-tresult[which(tresult[,3]<0.05),]
    up<-a1[which(a1[,2]>0),1] #sensitivity
    down<-a1[which(a1[,2]<0),1] #resistance
    
    a2<-tresult[which(tresult[,3]<0.05),]
    up2<-a2[which(a2[,2]>0),1] #resistance
    down2<-a2[which(a2[,2]<0),1] #sensitivity
    
    setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/0_CNV/T-test")
  write.table(a1,"cnv_result_T_test_sensitivity_amp.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  
  want_drug<-cnv_t_result[cnv_t_result$`T-p.value`<0.05,]
  
  
  setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/0_CNV/T-test")
  amp_gene1<-read.table("cnv_result_T_test_sensitivity_amp.txt",sep = '\t',header = T)
  genea=bitr(amp_gene1[,1],fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
  kegga<-enrichKEGG(genea[,2],organism = 'hsa',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'none',
                    minGSSize = 10,maxGSSize = 500,qvalueCutoff = 1,use_internal_data = FALSE)
  ka<-barplot(kegga,showCategory=10,drop=T,x = "GeneRatio",color = "pvalue")
  eea<-kegga@result
  eea1<-eea[which(eea$pvalue<0.05),]
  setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/0_CNV/T-test")
  write.table(eea1,"kegg_CNV_amp_p0.05_pathway_P0.05.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
  ggsave(paste("kegg_CNV_amp_pathway_P0.05_Ttest_DEGs_P0.05.pdf",sep=""),ka,width = 8, height = 6)
  
  
  dhyper(x, m, n, k)
############    RNAseq DEGs
rm(list=ls())
setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
library("data.table")
exp1<-fread("star_rsem.GeneSymbol.readscount.xls",header=T,data.table=F)
exp2<-floor(exp1[-1])
delete<-apply(exp2,1,function(x) mean(x==0))
cou<-which(delete>0.5)####在超过一半样本以上的基因删去
exp<-exp2[(-cou),]
sensitivity<-c("PC.101","PC.109","PC.115","PC.116","PC.117","PC.13","PC.130","PC.136","PC.139","PC.16","PC.18","PC.2","PC.27","PC.64","PC.78","PC.8","PC.81","PC.97","PC.98","PC.L" )
resistance<-c("PC.102","PC.104","PC.105","PC.111","PC.112","PC.119","PC.121","PC.134","PC.135","PC.14","PC.22","PC.40","PC.5","PC.52","PC.56","PC.G","PC.I" )
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]
resistance_exp<-exp[,match(resistance_new,colnames(exp))]
sensitivity_exp<-exp[,match(sensitivity_new,colnames(exp))]
geneid<-exp1[(-cou),1]
library(DESeq2)
colDate<-data.frame(row.names = c(as.vector(resistance_new),as.vector(sensitivity_new)),
                    condition=factor(c(rep("resistance",length(resistance_new)),rep("sensitivity",length(sensitivity_new))))
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
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/1_RNAseq")
write.table(res,"sensitivity-resistance-all-DESeq2_HDAC_chemotherapy_AUC_RNAseq_DEGs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
#resSig<-res[which(res$pvalue<0.05 & abs(res$log2FoldChange>1)),]
resSig<-res[which(res$padj<0.05),]
#pvalue padj
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

setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/1_RNAseq")
write.table(up_gene,"sensitivity-resistance-all-DESeq2_AUC_RNAseq_FDR0.05_FC0_up_DEGs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(down_gene,"sensitivity-resistance-all-DESeq2_AUC_RNAseq_FRD0.05_FC0_down_DEGs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

#########  RNAseq DEGs Volcano Plot
rm(list=ls())
library(ggrepel) 
library(ggplot2)
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/1_RNAseq")
df<-read.table("sensitivity-resistance-all-DESeq2_HDAC_chemotherapy_AUC_RNAseq_DEGs.txt",sep = '\t',header= T,row.names = 1)
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
data = df[df$padj<0.05&abs(df$log2FoldChange)>2,]
DEGs<-na.omit(as.matrix(data))

DEGs<-df[which(df$padj<0.01 & abs(df$log2FoldChange)>2),]
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/1_RNAseq/Volcano")
write.table(DEGs,"FDR05_FC2_DEGs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
ggsave("Volcano_FDR0.05_FC2.pdf",p,width = 8, height = 6)


#############  TCGA_PAAD   ####### 有符合条件的基因
rm(list=ls())
library(org.Hs.eg.db)
library(stringr)
library(clusterProfiler)
library(survival)
library(survminer)


setwd("~/xjj/drug/drug_result/chemotherapy_frontiers/6_TCGA")
chemotherapy_drug_sample<-read.table("drug_paad1.txt",header=T,sep="\t") #就用这个文件，化疗药物信息
#chemotherapy_drug_sample<-chemotherapy_drug_sample1[which(chemotherapy_drug_sample1$pharmaceutical_therapy_drug_name=="Oxaliplatin"),]
#chemotherapy_drug_sample<-chemotherapy_drug_sample1[which(chemotherapy_drug_sample1$pharmaceutical_therapy_drug_name %in% c("5-FU","5 FU","5FU","5-Fluorouracil","5-fu","5-fluorouracil","5-Fluorouracil?")),]
#chemotherapy_drug_sample<-read.table("drug_paad1.txt",header=T,sep="\t") #化疗药物信息
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
g2<-c("ZNF684","NR2C2","RFX5","ZNF528","ZNF384","SPI1","ZNF382","ZNF140","BCL6B","MEF2A","Foxd3","TFAP2A","RFX4","RFX2")
g2<-c("ESR2","RARA::RXRA","IRF2","MTF1","TFAP2C","Rarb","STAT1::STAT2","NR6A1","FOXK1","SRF","Plagl1","E2F6"
      ,"PBX3","ONECUT3","Nr1h3::Rxra","ZNF263","ZNF460","KLF4","ZNF16","EWSR1-FLI1","MAF","ZNF135","ZNF148","RREB1","Nr5a2","KLF2","TP63","EHF","KLF9","NFIB")

g1<-intersect(g2,rownames(exp1))
result<-matrix(0,13,3)
for(i in 1:length(g1)){
  g<-g1[i]
  result[i,1]<-g1[i]
  median<-median(as.numeric(exp1[match(g,rownames(exp1)),]))###中位数，基因需要换
  up<-which(as.numeric(exp1[match(g,rownames(exp1)),])>=median)
  down<-which(as.numeric(exp1[match(g,rownames(exp1)),])<median)
  class_ind<-matrix(0,ncol(exp1),1)
  class_ind[up,1]<-"high"   #基因表达更高，越敏感，理论上生存越好
  class_ind[down,1]<-"low" #基因表达更低，越耐药，理论上生存越差
  exp1_colname1<-substring(colnames(exp1),1,12)
  exp1_colname<-gsub("\\.","-",exp1_colname1) #提取符合条件的样本名
  surv_info<-state_time[match(exp1_colname,state_time[,1]),c(2,4)]#生存时间和生存状态
  label<-as.matrix(class_ind)
  TTP=as.numeric(surv_info[,2])  ##生存时间
  status_TTP=as.matrix(surv_info[,1])  ##生存状态
  status_TTP[which(status_TTP=="Alive")]=0
  status_TTP[which(status_TTP=="Dead")]=1
  status_TTP<-as.numeric(status_TTP)
  surv_info1<-as.data.frame(cbind(TTP,status_TTP))
  
  #surv_TTP<-survfit(Surv(TTP, status_TTP) ~ label,data=surv_info1)
  #ggsurvplot(surv_TTP,
  #pval = TRUE, conf.int = TRUE, #是否显示p值/显示风险区域
  #risk.table = TRUE,# 将风险表显示在生存曲线下面
  #ggtheme = theme_bw(),
  #tables.theme = theme_cleantable(),
  #title=paste("    ",g,"up-regulation in resistance",sep=" "))
             
  coxp<-coxph(Surv(TTP, status_TTP) ~ label,data=surv_info1)
  result[i,1]<-g1[i]
  result[i,2]<-tidy(coxp)$estimate
  result[i,3]<-tidy(coxp)$p.value
}
colnames(result)<-c("geneid","estimate","p.value")

summary(coxph(Surv(TTP, status_TTP) ~ label,data=surv_info1)) 
surv_TTP<-survfit(Surv(TTP, status_TTP) ~ label,data=surv_info1)
library(ggplot2)
ggsurvplot(surv_TTP,
           pval = TRUE, conf.int = TRUE, #是否显示p值/显示风险区域
           risk.table = TRUE,# 将风险表显示在生存曲线下面
           ggtheme = theme_bw(),
           tables.theme = theme_cleantable(),
           title=paste("    ",g,"up-regulation in sensitive",sep=" "))
plot(P)
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/4_TCGA_PAAD/TF_surviver")
library(ggplot2)
##手动保存
ggsave(paste(g,"up-regulation in sensitive.pdf",sep=" "),p,width = 5, height = 5)
write.table(result,"TFs0.01_down_TCGA_survived.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出


#################################################################################
################  DEGs与组蛋白修饰酶的交叠
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDAC_frontiers/5_histone-modifying enzymes")
Acetylase<-read.table("histone-modifying-enzymes.txt",header=T,sep="\t")
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/1_RNAseq")
df<-read.table("sensitivity-resistance-all-DESeq2_HDAC_chemotherapy_AUC_RNAseq_DEGs.txt",sep = '\t',header= T,row.names = 1)
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
input  <-list(Acetylase_name,DEGs_name)
venn(input,showSetLogicLabel=TRUE)
tmp <- venn(input)
int<-attr(tmp, "intersections")
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/1_RNAseq/Volcano")
write.table(Acetylase_info,"Acetylase_info_P0.05.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
library(VennDiagram)
venn.diagram(x=list(Acetylase=Acetylase_name,DEGs=DEGs_name),cex = 1,margin = 0.1, "Acetylase_info_P005.png",fill=c("red","blue"))

setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
library("data.table")
exp<-fread("star_rsem.GeneSymbol.TPM.xls",header=T,data.table=F)
sensitivity<-c("PC.101","PC.109","PC.115","PC.116","PC.117","PC.13","PC.130","PC.136","PC.139","PC.16","PC.18","PC.2","PC.27","PC.64","PC.78","PC.8","PC.81","PC.97","PC.98","PC.L" )
resistance<-c("PC.102","PC.104","PC.105","PC.111","PC.112","PC.119","PC.121","PC.134","PC.135","PC.14","PC.22","PC.40","PC.5","PC.52","PC.56","PC.G","PC.I" )
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]
resistance_exp<-exp[,match(resistance_new,colnames(exp))]
sensitivity_exp<-exp[,match(sensitivity_new,colnames(exp))]
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
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/1_RNAseq/Volcano")
write.table(Acetylase_info_value,"Acetylase_info_P0.05_value.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

library(ggplot2)
df2<-cbind(rbind(as.matrix(Acetylase_info_value[,1]),as.matrix(Acetylase_info_value[,1])),
           rbind(as.matrix(sensitivity_exp11),as.matrix(resistance_exp11)),
           rbind(as.matrix(rep("sensitive",17)),as.matrix(rep("resistant",17))),
           rbind(as.matrix(Acetylase_info_value$sensitive+sensitivity_sd),as.matrix(Acetylase_info_value$resistant+resistance_sd)),
           rbind(as.matrix(Acetylase_info_value$sensitive-sensitivity_sd),as.matrix(Acetylase_info_value$resistant-resistance_sd)))
colnames(df2)<-c("gene","value","type","valuesd1","valuesd2")
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/1_RNAseq/Volcano")
write.table(df2,"Acetylase_info_P0.05_value_plot.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

df3<-read.table("Acetylase_info_P0.05_value_plot.txt",sep = '\t',header= T)
temp=df3[df3[,3]=="sensitive",1]
df3[,1]=factor(df3[,1],levels=temp)
df3[which(df3$pvalue>0.01),'pva']<-"*"
df3[which(df3$pvalue<0.01 & df3$pvalue>0.001),'pva']<-"**"
df3[18:34,7]=NA
df3$valuesd1<-df3$value+1
p3 <- ggplot(df3, aes(x = gene, y = value, fill = type)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge(),width=0.4)+ 
  geom_errorbar(aes(ymin = valuesd1, ymax = value), width = 0.2, position = position_dodge(0.4))+
  labs(x = "Differently expression of histone-modifying enzymes", y = "mRNA average expression levels (TPM)")+
  geom_text(aes(y=60,label = pva),size = 3)
plot(p3)
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/1_RNAseq/Volcano")
ggsave("Acetylase_info_P0.05_value_plot.pdf",p3,width = 10, height = 6)


#################################################################################
################  clusterProfiler
rm(list=ls())
library(org.Hs.eg.db)
library(clusterProfiler)
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/1_RNAseq")
up_gene<-read.table("sensitivity-resistance-all-DESeq2_AUC_RNAseq_P0.05_FC0_down_DEGs.txt",header=T,sep="\t")
library(stringr)
gene=bitr(up_gene[,1],fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
kegg<-enrichKEGG(gene[,2],organism = 'hsa',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'BH',
                 minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.05,use_internal_data = FALSE)
keu<-barplot(kegg,showCategory=13,drop=T,x = "GeneRatio",color = "p.adjust")

ee<-kegg@result
ee1<-ee[which(ee$p.adjust<0.05),]
#p.adjust
#dotplot(kegg,showCategory=8,x = "GeneRatio",color = "p.adjust",title = "DEGs-P05-KEGG")
#enrichplot::gseaplot2(ee1,1,pvalue_table=T,color="#086538")
enrich_gene<-ee1$geneID
pathway_gene<-unique(unlist(strsplit(enrich_gene,split="/")))
pathway_gene2=bitr(pathway_gene,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/1_RNAseq/pathway")
write.table(ee1,"kegg_DEG_pathway_FDR0.05_down_drug_target.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(pathway_gene2,"kegg_DEG_pathway_FDR0.05_down_drug_target_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
ggsave("kegg_DEG_pathway_FDR0.05_down_drug_target_gene.pdf",keu,width = 8, height = 6)

#dotplot(kegg,showCategory=20)

go<-enrichGO(gene[,2],OrgDb = org.Hs.eg.db,ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.01,keyType = 'ENTREZID')
ked<-dotplot(go,showCategory=10)
a<-go@result
go_BP<-a[which(a$p.adjust<0.05),]

enrich_genego<-go_BP$geneID
pathway_genego<-unique(unlist(strsplit(enrich_genego,split="/")))
pathway_genego2=bitr(pathway_genego,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/1_RNAseq/pathway")
write.table(pathway_genego2,"GO_DEG_pathway_FDR0.05_down_drug_target_gene.txt",sep = '\t',col.names = F,row.names = F,quote = FALSE)##数据输出
write.table(go_BP,"GO_DEG_pathway_FDR0.05_down_drug_target.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
ggsave("GO_DEG_pathway_FDR0.01_up_drug_target_gene.pdf",ked,width = 8, height = 6)



#######################    ATACseq   ################################################
#############  DApeaks
rm(list=ls())
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_count_name<-read.table("peak_count.txt",sep="\t",header=T)
sensitivity<-c("PC.101","PC.109","PC.115","PC.116","PC.117","PC.13","PC.130","PC.136","PC.139","PC.16","PC.18","PC.2","PC.27","PC.64","PC.78","PC.8","PC.81","PC.97","PC.98","PC.L" )
resistance<-c("PC.102","PC.104","PC.105","PC.111","PC.112","PC.119","PC.121","PC.134","PC.135","PC.14","PC.22","PC.40","PC.5","PC.52","PC.56","PC.G","PC.I" )
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]
gsub("\\.","-",colnames(peak_count_name))
resistance_peak<-peak_count_name[,match(resistance_new,gsub("\\.","-",colnames(peak_count_name)))]
sensitivity_peak<-peak_count_name[,match(sensitivity_new,gsub("\\.","-",colnames(peak_count_name)))]
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
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/2_ATACseq_DApeaks")
write.table(res,"sensitivity-resistance-all-DESeq2_HDAC_chemotherapy_AUC_ATAC.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
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
out_file<-"sen-res-HDAC_chemotherapy_AUC-0.05-P"
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
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/2_ATACseq_DApeaks/annotation")
write.table(peaks_bed_up,file="all_peaks.bed",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')


#########  暂时运行
rm(list=ls())
freq<-0.05
setwd("~/176/Changhai_ATAC/results/bwa/mergedReplicate/macs/peaks_merge/txt_result")
peak_count_name<-read.table("peak_count.txt",sep="\t",header=T)
setwd("~/xjj/drug/drug_result/drug64_AUC")
tresult_want<-read.table("drug64_sample.txt",sep="\t",header=T)
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
library(DESeq2)
library(dplyr)

for(k in 6:20){
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
  colDate<-data.frame(row.names = c(as.vector(resistance),as.vector(sensitivity)),
                      condition=factor(c(rep("resistance",length(resistance)),rep("sensitivity",length(sensitivity))))
  )
  datexpr<-cbind(resistance_peak,sensitivity_peak)##前面的相对于后面的上下调
  counts <- apply(datexpr,2,as.numeric)   ###矩阵中必须是数值
  dds<-DESeqDataSetFromMatrix(countData = counts,colData = colDate,design = ~condition)
  dds<-DESeq(dds)##进行标准化分析
  rld <- rlogTransformation(dds)  ## 得到经过DESeq2软件normlization的表达矩阵！
  exprSet_new<-assay(rld)
  rownames(exprSet_new)<-peak_name
  cat(exprSet_new[1,],"\n")
  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/2_ATACseq_DApeaks",sep=""))
  write.table(exprSet_new,paste("sensitivity-resistance-",file_name,"-DESeq2_normlization.txt",sep=""),sep = '\t',col.names = T,row.names = T,quote = FALSE)##数据输出
}


#######################################################################################
####  homer注释完的信息,画饼图
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/2_ATACseq_DApeaks/annotation")
library(data.table)
homer_anno<-fread("sen-res-HDAC_chemotherapy_AUC-0.01-P-logFC0-DESeq2-up-annotation.txt",header=T,data.table=F)
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
  result[i,2]<-sum(Annotation %in% position[i])
  result[i,3]<-(sum(Annotation %in% position[i])/(length(Annotation)))*100
}

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
  labs(title = "down0.05-FC1-peaks")+
  theme(plot.title = element_text(hjust = 0.5))
#geom_text(aes(y = A/2 + c(0, cumsum(A)[-length(A)]), x = sum(A)/20, label = myLabel), size = 5)   # 在图中加上百分比：x 调节标签到圆心的距离, y 调节标签的左右位置
print(p) #显示饼图
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/2_ATACseq_DApeaks/annotation")
ggsave("down0.05_FC1_annotation.pdf",p,width = 6, height = 5)



###############  peaks都处于什么位置  ######################    
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/2_ATACseq_DApeaks/annotation")
library(data.table)
homer_anno1<-fread("all_peaks-annotation.txt",header=T,data.table=F)

setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/2_ATACseq_DApeaks")
res<-read.table("sensitivity-resistance-all-DESeq2_HDAC_chemotherapy_AUC_ATAC.txt",header = TRUE,sep = "\t")
resSig<-res[which(res$pvalue<0.01),]
int_peak<-intersect(resSig[,1],homer_anno1[,1])

homer_anno<-homer_anno1[match(int_peak,homer_anno1[,1]),]
one<-sum(abs(homer_anno$`Distance to TSS`)<=100)
two<-length(which(abs(homer_anno$`Distance to TSS`)>=100 & abs(homer_anno$`Distance to TSS`)<=1000))
three<-length(which(abs(homer_anno$`Distance to TSS`)>=1000 & abs(homer_anno$`Distance to TSS`)<=10000))
four<-length(which(abs(homer_anno$`Distance to TSS`)>=10000 & abs(homer_anno$`Distance to TSS`)<=100000))
five<-sum(abs(homer_anno$`Distance to TSS`)>=100000)
ma<-matrix(rep("all-peaks",5),ncol=1)
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
  labs(x = "", y = "", title = "DApeaks-position") +  # 将横纵坐标的标签设为空
  theme(axis.ticks = element_blank()) +  # 将左上角边框的刻度去掉
  theme(legend.title = element_blank(), legend.position = "left")+   ## 将图例标题设为空，并把图例方放在左边位置
  scale_fill_discrete(breaks = dt$B, labels = myLabel)+   # 将原来的图例标签换成现在的myLabel
  theme(axis.text.x = element_blank())+   ## 去掉饼图的外框上的数值，即去除原柱状图的X轴，把X轴的刻度文字去掉
  theme(plot.title = element_text(hjust = 0.5))
  #geom_text(aes(y = A/2 + c(0, cumsum(A)[-length(A)]), x = sum(A)/20, label = myLabel), size = 5)   # 在图中加上百分比：x 调节标签到圆心的距离, y 调节标签的左右位置
print(p) #显示饼图
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/2_ATACseq_DApeaks/annotation")
ggsave("DApeaks-0.01-position.pdf",p,width = 6, height = 5)


#####################  符合在TSS100KB之内条件的注释基因
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/2_ATACseq_DApeaks/annotation")
library(data.table)
homer_anno_down<-fread("sen-res-HDAC_chemotherapy_AUC-0.01-P-logFC0-DESeq2-down-annotation.txt",header=T,data.table=F)
colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]

homer_anno_up<-fread("sen-res-HDAC_chemotherapy_AUC-0.01-P-logFC0-DESeq2-up-annotation.txt",header=T,data.table=F)
colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]

setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/1_RNAseq")
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

setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/1_2_intersect/venn")
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

setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/1_2_intersect/venn")
write.table(int_up,"Up-peak0.01-Up-DEG0.05.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(int_down,"Down-peak0.01-Down-DEG0.05.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(int_up_peak,"int_up_peak-Up-peak0.01-Up-DEG0.05.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(int_down_peak,"int_down_peak-Down-peak0.01-Down-DEG0.05.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出



#################################################################
##############    ATAC-seq CNV
setwd("~/xjj/drug/drug_result/HDAC_frontiers_AUC/6_CNV")
amp_CNV<-read.table("sensitivity-resistance-amp-CNV-gene.txt",sep = '\t',header=T)##数据输出
censored_CNV<-read.table("sensitivity-resistance-censored-CNV-gene.txt",sep = '\t',header=T)##数据输出
setwd("~/xjj/WGS_CNV/CNV_FACETS")
gene_position<-read.table("gene_position.txt",sep = '\t',header=T)##数据输出
amp_position<-gene_position[match(amp_CNV[,1],gene_position[,5]),]
censored_position<-gene_position[match(censored_CNV[,1],gene_position[,5]),]
setwd("~/xjj/drug/drug_result/HDAC_frontiers_AUC/2_ATACseq_DApeaks")
library(data.table)
homer_anno_down<-fread("sen-res-HDAC_AUC-0.05-P-down-annotation.txt",header=T,data.table=F)
colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]

homer_anno_up<-fread("sen-res-HDAC_AUC-0.05-P-up-annotation.txt",header=T,data.table=F)
colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]
amp_CNV_DApeak<-intersect(amp_position[,5],Anno_gene_100Kb_up$`Gene Name`)
amp_CNV_DApeak1<-intersect(amp_position[,5],Anno_gene_100Kb_down$`Gene Name`)

censored_CNV_DApeak1<-intersect(censored_position[,5],Anno_gene_100Kb_down$`Gene Name`)



##################################################################################
##############  DEpeaks在不同样本变化的倍数与DEGs在不同样本变化的倍数之间的相关性
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/2_ATACseq_DApeaks/annotation")
library(data.table)
homer_anno_down<-fread("sen-res-HDAC_chemotherapy_AUC-0.01-P-logFC0-DESeq2-down-annotation.txt",header=T,data.table=F)
colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]

homer_anno_up<-fread("sen-res-HDAC_chemotherapy_AUC-0.01-P-logFC0-DESeq2-up-annotation.txt",header=T,data.table=F)
colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]

setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/1_RNAseq")
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

setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/2_ATACseq_DApeaks")
DApeak<-fread("sensitivity-resistance-all-DESeq2_HDAC_chemotherapy_AUC_ATAC.txt",header=T,data.table=F)
DEGs_peaks_FC<-DApeak[match(DEGs_peaks,DApeak$peak_id),3]
DEGs_FC<-all[match(int_all_DEGs,DEGs),3]

#pearson", "kendall", "spearman
cor_result<-cor.test(DEGs_peaks_FC, DEGs_FC,alternative = "two.sided",method = "spearman")
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
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/1_2_intersect/correlation")
write.table(dat1,"DEGs_FC_DEGs_peaks_FC_correlation.txt",sep = '\t',col.names = T,row.names = T,quote = FALSE)##数据输出

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
    aes(label = rownames(want_dat)),
    size = 3,
    color = "black",
    segment.color = "black", show.legend = FALSE )
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/1_2_intersect/correlation")
ggsave("correlation-DApeaksP0.01-DEGsP0.05-pearson.pdf",ps,width = 10, height = 10)

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


#############################  对高低可及peaks-DEGs进行KEGG通路富集分析
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/2_ATACseq_DApeaks/annotation")
library(data.table)
homer_anno_down<-fread("sen-res-HDAC_chemotherapy_AUC-0.01-P-logFC0-DESeq2-down-annotation.txt",header=T,data.table=F)
colnames(homer_anno_down)<-c("peakid",colnames(homer_anno_down)[-1])
Anno_gene_100Kb_down<-homer_anno_down[which(abs(homer_anno_down$`Distance to TSS`)<=100000),]

homer_anno_up<-fread("sen-res-HDAC_chemotherapy_AUC-0.01-P-logFC0-DESeq2-up-annotation.txt",header=T,data.table=F)
colnames(homer_anno_up)<-c("peakid",colnames(homer_anno_up)[-1])
Anno_gene_100Kb_up<-homer_anno_up[which(abs(homer_anno_up$`Distance to TSS`)<=100000),]

setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/1_RNAseq")
up_gene<-read.table("sensitivity-resistance-all-DESeq2_AUC_RNAseq_P0.05_FC0_up_DEGs.txt",header=T,sep="\t")
down_gene<-read.table("sensitivity-resistance-all-DESeq2_AUC_RNAseq_P0.05_FC0_down_DEGs.txt",header=T,sep="\t")
#all<-rbind(up_gene,down_gene)
#DEGs<-as.character(all[,1])
down_DEGs_DApeaks<-intersect(Anno_gene_100Kb_down$`Gene Name`,down_gene[,1])
up_DEGs_DApeaks<-intersect(Anno_gene_100Kb_up$`Gene Name`,up_gene[,1])

library(org.Hs.eg.db)
library(clusterProfiler)
up_gene<-down_DEGs_DApeaks  ## 修改
gene=bitr(up_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
## 去重
library(stringr)
kegg<-enrichKEGG(gene[,2],organism = 'hsa',keyType = 'kegg',pvalueCutoff = 0.01,pAdjustMethod = 'none',
                 minGSSize = 10,maxGSSize = 500,qvalueCutoff = 1,use_internal_data = FALSE)
ee<-kegg@result
ee1<-ee[which(ee$pvalue<0.01),]
#ee1<-ee[which(ee$p.adjust<0.05),]

p<-barplot(kegg,showCategory=15,drop=T,x = "GeneRatio",color = "pvalue")
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


setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/1_2_intersect/pathway")
write.table(ee1,"DEGP0.05_DApeaksP0.01_kegg_pathway_P0.01_down.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(pathway_gene,"DEGP0.05_DApeaksP0.01_kegg_pathway_P0.01_down_gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
ggsave("DEGP0.05_DApeaksP0.01_kegg_pathway_P0.01_down.pdf",p,width = 5, height = 3)
write.table(result2,"DEGP0.05_DApeaksP0.01_kegg_pathway_P0.01_up_gene_symbol.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(a4,"DEGP0.05_DApeaksP0.01_kegg_pathway_P0.01_up_pathway_symbol.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出


#############################################################################
######################  使用homer进行motif富集，并对应到TFs
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDAC_frontiers_AUC/2_ATACseq_DApeaks/homer_result/sen-res-HDAC_AUC-0.05-P-logFC0-DESeq2-down")
knownResults<-read.delim("knownResults.txt",header=T,sep='\t')
FDR05<-knownResults[knownResults$q.value..Benjamini.<0.05,]
#FDR05<-knownResults[knownResults$P.value<0.01,]
library(dplyr)
homer_TF1<-data.frame(lapply(strsplit(as.character(FDR05[,1]),'/'), function(x) x[3])%>%unlist())
FDR051<-FDR05[which(homer_TF1=="Homer"),]
homer_TF2<-data.frame(lapply(strsplit(as.character(FDR051[,1]),'/'), function(x) x[1])%>%unlist())
library(stringr)
homer_TF3<-apply(homer_TF2,1,function(x) str_replace_all(as.character(x),"\\((.*?)\\)",""))##取出括号前的字符
homer_TF4<-unique(homer_TF3)
length(homer_TF4)
new_homer_result<-cbind(homer_TF3,FDR051[,-1])
colnames(new_homer_result)<-c("TF_name",colnames(new_homer_result)[-1])
setwd("~/xjj/drug/drug_result/HDAC_frontiers_AUC/4_TF/homer_result")
write.table(new_homer_result,"DApeaks_P0.05-logFC0-DESeq2-up_FDR0.05_TFs.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

#####  画柱状图
data1<-cbind(as.data.frame(new_homer_result[,1]),new_homer_result$P.value)
data1[,2]<-(-log2(data1[,2]))
colnames(data1)<-c("TFname","Pvalue")
data11<-data1[match(unique(data1[,1]),data1[,1]),]##重叠的TFqu取第一个
data2<-data11[1:10,]
data2$TFname=factor(data2$TFname,levels = data2[,1])
new_homer_result1<-new_homer_result[match(unique(data1[,1]),data1[,1]),]
val<-new_homer_result1$X..of.Target.Sequences.with.Motif[1:10]
library(ggplot2)
ggplot(data=data2,mapping=aes(x=TFname,y=as.numeric(Pvalue)))+
  geom_bar(stat="identity",fill="blue4")+
  #geom_text(aes(label = val, vjust = -0.8, hjust = 0.5), show.legend = TRUE)+
  ylab("-log2(Pvalue)")+
  ggtitle("DApeaks_P0.05-logFC0-DESeq2-down_FDR0.05_TFs") +
  theme(plot.title = element_text(hjust = 0.5))


#############################################################################
######################  使用FIMO进行motif富集，并对应到TFs
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDAC_frontiers/4_TFs/FIMO_result")
library(data.table)
fimo<-fread("up05fimo.tsv",header=T,data.table=F)
colnames(fimo)<-c(colnames(fimo)[1:8],"FDR","matched_sequence")
library(dplyr)
TF_fimo<-fimo %>% distinct(motif_id,motif_alt_id,FDR, .keep_all = TRUE)
df<-TF_fimo[which(TF_fimo$FDR<0.05),]
length(unique(df$motif_alt_id))

TF_motif<-data.frame()
for(i in 1:length(unique(df[,2]))){
  TF_motif[i,1]<-unique(df[,2])[i]
  TF_motif[i,2]<-sum(df[,2] %in% unique(df[,2])[i])
}
TF_motif_new<-TF_motif[order(TF_motif[,2],decreasing = T),]
colnames(TF_motif_new)<-c("TF","number")

TF_motif_new$TF=factor(TF_motif_new$TF,levels = TF_motif_new[,1])
ggplot(data = TF_motif_new[1:20,],mapping = aes(x = TF[1:20], y =number[1:20]))+ 
  geom_bar(stat = 'identity',position="dodge")+
  geom_text(aes(label = number, vjust = -0.8, hjust = 0.5), show.legend = TRUE)+
  ylab("Top20-TFs")+
  xlab("TFs-names")+
  ggtitle("up05fimo_TFs") +
  theme(plot.title = element_text(hjust = 0.5))

######################################################################################################
##################  ridge regresion 差异表达转录因子(DETF)及其与DApeaks相关的目标DEGs  上下调的转录因子是分开看的
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/1_RNAseq")
up_gene<-read.table("sensitivity-resistance-all-DESeq2_AUC_RNAseq_P0.05_FC0_down_DEGs.txt",header=T,sep="\t")

setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/2_ATACseq_DApeaks/FIMO_result")
knownResults<-read.delim("DApeaksP0.01_ridge_regression_down_clean.txt",header=T,sep='\t')
int_TF<-intersect(as.character(knownResults[,1]),as.character(up_gene[,1]))



#FDR05<-knownResults[knownResults$P.value<0.05,]
FDR05<-knownResults[knownResults$q.value..Benjamini.<0.05,]

library(dplyr)
homer_TF1<-data.frame(lapply(strsplit(as.character(FDR05[,1]),'/'), function(x) x[3])%>%unlist())
FDR051<-FDR05[which(homer_TF1=="Homer"),]
homer_TF2<-data.frame(lapply(strsplit(as.character(FDR051[,1]),'/'), function(x) x[1])%>%unlist())
library(stringr)
homer_TF3<-apply(homer_TF2,1,function(x) str_replace_all(as.character(x),"\\((.*?)\\)",""))##取出括号前的字符
homer_TF4<-unique(homer_TF3)
length(homer_TF4)
new_homer_result<-cbind(homer_TF3,FDR051[,-1])
colnames(new_homer_result)<-c("TF_name",colnames(new_homer_result)[-1])
length(intersect(toupper(new_homer_result$TF_name),toupper(up_gene$gene_id)))
DETFs<-intersect(toupper(new_homer_result$TF_name),toupper(up_gene$gene_id))
setwd("~/xjj/drug/drug_result/HDAC_frontiers_AUC/3_intersect/DETFs_venn")
library(VennDiagram)
venn.diagram(x=list(TF_enrich_hyper_DApeaks=toupper(new_homer_result$TF_name),up_DEGs=toupper(up_gene$gene_id)), "up_DETFs_P005.png",fill=c("red","blue"),margin = 0.3)

###############   TF找到对应的靶基因
DETFs<-c("ZNF684","NR2C2","RFX5","ZNF528","ZNF384","SPI1","ZNF382","ZNF140","BCL6B","MEF2A","Foxd3","TFAP2A","RFX4","RFX2")
DETFs<-c("ESR2","RARA::RXRA","IRF2","MTF1","TFAP2C","Rarb","STAT1::STAT2","NR6A1","FOXK1","SRF","Plagl1","E2F6"
      ,"PBX3","ONECUT3","Nr1h3::Rxra","ZNF263","ZNF460","KLF4","ZNF16","EWSR1-FLI1","MAF","ZNF135","ZNF148","RREB1","Nr5a2","KLF2","TP63","EHF","KLF9","NFIB")

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
result<-NULL
for(i in 1:length(DETF_want)){
  result1<-target_TF_gene[target_TF_gene[,2] %in% DETF_want[i],]
  result<-rbind(result,result1)
}
### 查看靶标基因中有多少是DEGs
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/1_RNAseq")
up_gene<-read.table("sensitivity-resistance-all-DESeq2_AUC_RNAseq_P0.05_FC0_up_DEGs.txt",header=T,sep="\t")
down_gene<-read.table("sensitivity-resistance-all-DESeq2_AUC_RNAseq_P0.05_FC0_down_DEGs.txt",header=T,sep="\t")

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
colnames(data)<-c("DETFs","target-DEGs","target-DEGs-dysregulation")

aa<-unique(data[,1])
re<-data.frame()
for(m in 1:length(unique(data[,1]))){
  dat<-data[data[,1]%in%aa[m],]
  up<-sum(dat[,3]=="up")
  down<-sum(dat[,3]=="down")
  re[m,1]<-aa[m]
  re[m,2]<-up
  re[m,3]<-down
}
colnames(re)<-c("geneid","up","down")
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/2_ATACseq_DApeaks/TF_target")
write.table(data,"sen_res_down_TFs_TRANSFAC-gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(data,"sen_res_down_TFs_JASPAR-gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(data,"sen_res_down_TFs_MotifMap-gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(data,"sen_res_down_TFs_ENCODE-gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
write.table(re,"number_of_down_JASPAR-gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出

###################################################################################
################   分类器 randomForest
rm(list=ls())
setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
library("data.table")
exp1<-fread("star_rsem.GeneSymbol.TPM.xls",header=T,data.table=F)
exp2<-log2(exp1[-1]+1)
delete<-apply(exp2,1,function(x) mean(x==0))
cou<-which(delete>0.5)####在超过一半样本以上的基因删去
exp<-exp2
resistance<-c("PC.104","PC.5","PC.105","PC.56","PC.119","PC.135","PC.G","PC.102","PC.I","PC.14","PC.40","PC.134","PC.111","PC.121","PC.112","PC.22")
sensitivity<-c("PC.97","PC.130","PC.136","PC.2","PC.L","PC.117","PC.101","PC.27","PC.64","PC.81","PC.115","PC.116","PC.16","PC.78","PC.98","PC.18","PC.109","PC.139","PC.13","PC.52","PC.8")
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]
resistance_exp<-exp[,match(resistance_new,colnames(exp))]
sensitivity_exp<-exp[,match(sensitivity_new,colnames(exp))]
geneid<-exp1[,1]

data<-cbind(geneid,resistance_exp,sensitivity_exp) ####37个样本的表达谱
symbol<-as.matrix(data[,1])####基因名
glist<-symbol
exp<-t(data[,-1])  ###纯表达谱
#colnames(exp)<-paste("g",1:length(symbol),sep="")
colnames(exp)<-symbol
exp<-as.data.frame(exp)
gt<-rep(c(0,1),c(length(resistance),length(sensitivity))) #ground_truth 耐药17 敏感20

setwd("~/xjj/drug/drug_result/HDAC_frontiers_AUC/1_RNAseq_DEGs/Acetylase")
Acetylase_info<-read.table("Acetylase_info_P0.05.txt",sep = '\t',header = T)


k=100
loc<-matrix(1,37,k)  ###构建一张37行100列的矩阵，用于做测试集的样本位置
rownames(loc)<-rownames(exp) ##37个样本名
set.seed(123)
for(i in 1:k){
  rr<-sample(1:length(resistance),4)  ##耐药中随机选取2个作为作为验证集
  sr<-sample((length(resistance)+1):37,4) ##敏感中随机选取2个样本作为验证集
  loc[c(rr,sr),i]<-2 #测试集
}
IMP<-matrix(0,length(symbol),2*k)  ###基因名长度，有200列
colnames(IMP)[seq(2,2*k,2)]<-"MeanDecreaseGini"
colnames(IMP)[seq(1,2*k,2)]<-"DRG"
rownames(IMP)<-symbol
ability<-matrix(NA,6,k)
rownames(ability)<-c("sen","spe","acc","val_sen","val_sep","val_acc")

####
for(j in 1:k)
{
  print(j)
  train.data<-exp[which(loc[,j]==1),]
  train.gt<-gt[which(loc[,j]==1)]
  test.data<-exp[which(loc[,j]==2),]
  test.gt<-gt[which(loc[,j]==2)]
  
  dat<-cbind(train.gt,train.data)
  Acetylase_DEGs<-Acetylase_info[,1]
  IMP[match(Acetylase_DEGs,colnames(exp)),2*j-1]<-1

  dat1<-dat[,c(1,match(Acetylase_DEGs,colnames(dat)))]
  dat1$train.gt[dat1$train.gt==0]="resistant"
  dat1$train.gt[dat1$train.gt==1]="sensitive"
  dat1$train.gt<-as.factor(dat1$train.gt)
  
  library(randomForest)
  set.seed(123)
  dat1.rf <- randomForest(train.gt ~ ., data=dat1, importance=T,proximity=T)
  print(dat1.rf) #展示随机森林模型简要信息
  ability[1,j]<-dat1.rf$ confusion[1,1]/(length(resistance)-4)
  ability[2,j]<-dat1.rf$ confusion[2,2]/(length(sensitivity)-4)
  ability[3,j]<-(dat1.rf$ confusion[1,1]+dat1.rf$ confusion[2,2])/29
  imp<-dat1.rf$importance
  IMP[match(rownames(imp),colnames(exp)),2*j]<-imp[,4]
  #提取随机森林模型中以准确率递减方法得到维度重要性值。type=2为基尼系数方法
  varImpPlot(dat1.rf,sort=TRUE,n.var=nrow(dat1.rf$importance), main = "variable importance")
  hist(treesize(dat1.rf))   #展示随机森林模型中每棵决策树的节点数
  max(treesize(dat1.rf));min(treesize(dat1.rf))
  MDSplot(dat1.rf, dat1$train.gt, palette=rep(1, 2), pch=as.numeric(dat1$train.gt)) #展示数据集在二维情况下各类别的具体分布情况
  pre<- predict(dat1.rf,newdata=test.data)
  print(pre)
  pred_out_1<-predict(object=dat1.rf,newdata=test.data,type="prob")
  print(pred_out_1)
  ability[4,j]<-length(grep("res",pre[1:4]))/4
  ability[5,j]<-length(grep("sen",pre[5:8]))/4
  ability[6,j]<-(length(grep("res",pre[1:4]))+length(grep("sen",pre[5:8])))/8
}
IMP1<-IMP[match(Acetylase_info[,1],colnames(exp)),]
setwd("~/xjj/drug/drug_result/HDAC_frontiers_AUC/5_classifier/randomForest/Acetylase_DEGs")
write.table(loc,"Acetylase_DEGs_train_test4.txt",sep="\t",col.names=F)
write.table(IMP1,"Acetylase_DEGs_MeanDecreaseGini4.txt",sep="\t",col.names=F)
write.table(ability,"Acetylase_DEGs_train_test_ability4.txt",sep="\t",col.names=F)

###############    分类器  SVM
rm(list=ls())
setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
library("data.table")
exp1<-fread("star_rsem.GeneSymbol.TPM.xls",header=T,data.table=F)
exp2<-log2(exp1[-1]+1)
delete<-apply(exp2,1,function(x) mean(x==0))
cou<-which(delete>0.5)####在超过一半样本以上的基因删去
exp<-exp2
resistance<-c("PC.104","PC.5","PC.105","PC.56","PC.119","PC.135","PC.G","PC.102","PC.I","PC.14","PC.40","PC.134","PC.111","PC.121","PC.112","PC.22")
sensitivity<-c("PC.97","PC.130","PC.136","PC.2","PC.L","PC.117","PC.101","PC.27","PC.64","PC.81","PC.115","PC.116","PC.16","PC.78","PC.98","PC.18","PC.109","PC.139","PC.13","PC.52","PC.8")
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
resistance_new<-sample_id[match(resistance,sample_id[,2]),1]
sensitivity_new<-sample_id[match(sensitivity,sample_id[,2]),1]
resistance_exp<-exp[,match(resistance_new,colnames(exp))]
sensitivity_exp<-exp[,match(sensitivity_new,colnames(exp))]
geneid<-exp1[,1]

data<-cbind(geneid,resistance_exp,sensitivity_exp) ####37个样本的表达谱
symbol<-as.matrix(data[,1])####基因名
glist<-symbol
exp<-t(data[,-1])  ###纯表达谱
#colnames(exp)<-paste("g",1:length(symbol),sep="")
colnames(exp)<-symbol
exp<-as.data.frame(exp)
gt<-rep(c(0,1),c(length(resistance),length(sensitivity))) #ground_truth 耐药17 敏感20

setwd("~/xjj/drug/drug_result/HDAC_frontiers_AUC/1_RNAseq_DEGs/Acetylase")
Acetylase_info<-read.table("Acetylase_info_P0.05.txt",sep = '\t',header = T)

#随机100次
k=100
loc<-matrix(1,37,k)
rownames(loc)<-rownames(exp)
set.seed(123)
for(i in 1:k){
  rr<-sample(1:length(resistance),4)
  sr<-sample((length(resistance)+1):37,4)
  loc[c(rr,sr),i]<-2 #测试集
}

IMP<-matrix(0,length(symbol),k)
rownames(IMP)<-symbol
ability<-matrix(NA,6,k)
rownames(ability)<-c("sen","spe","acc","val_sen","val_sep","val_acc")

library(e1071)
for(j in 1:k){
  print(j)
  train.data<-exp[which(loc[,j]==1),]
  train.gt<-gt[which(loc[,j]==1)]
  test.data<-exp[which(loc[,j]==2),]
  test.gt<-gt[which(loc[,j]==2)]
  
  dat<-cbind(train.gt,train.data)
  Acetylase_DEGs<-Acetylase_info[,1]
  IMP[match(Acetylase_DEGs,colnames(exp)),j]<-1
  dat1<-dat[,c(1,match(Acetylase_DEGs,colnames(dat)))]
  dat1$train.gt<-as.factor(dat1$train.gt)
  
  set.seed(100) # for reproducing results
  tuned<-tune.svm(train.gt~.,data = dat1,gamma = 10^(-6:-1),cost = 10^(0:4))
  summary(tuned)
  model.tuned<-svm(train.gt~.,data = dat1,gamma=tuned$best.parameters$gamma,cost=tuned$best.parameters$cost)
  print(summary(model.tuned))
  pred<-predict(model.tuned,dat1)
  
  ability[1,j]<-length(grep(0,pred[1:(length(resistance)-4)]))/(length(resistance)-4)
  ability[2,j]<-length(grep(1,pred[(length(resistance)-4+1):29]))/(29-length(resistance)+4-1)
  ability[3,j]<-(length(grep(0,pred[1:(length(resistance)-4)]))+length(grep(1,pred[(length(resistance)-4+1):29])))/29
  
  pre<- predict(model.tuned,newdata=test.data)
  print(pre)
  ability[4,j]<-length(grep(0,pre[1:4]))/4
  ability[5,j]<-length(grep(1,pre[5:8]))/4
  ability[6,j]<-(length(grep(0,pre[1:4]))+length(grep(1,pre[5:8])))/8
}
IMP1<-IMP[match(Acetylase_info[,1],colnames(exp)),]
setwd("~/xjj/drug/drug_result/HDAC_frontiers_AUC/5_classifier/SVM/Acetylase_DEGs")
write.table(loc,"Acetylase_DEGs_train_test4_SVM.txt",sep="\t",col.names=T)
write.table(IMP1,"Acetylase_DEGs_MeanDecreaseGini4_SVM.txt",sep="\t",col.names=T)
write.table(ability,"Acetylase_DEGs_train_test_ability4_SVM.txt",sep="\t",col.names=T)

###########  分类器   Logistic regression model


#############   HDAC+chemotherapy  #######################
rm(list=ls())
setwd("~/xjj/drug")
ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,!(colnames(ic501)%in%c("PC.100","PC.34"))]

drug_info<-read.csv("compiled_drug_info.csv",header = TRUE,sep = ",")
drug_result<-drug_info[match(rownames(ic50),drug_info$drug_id),]
HDAC1<-drug_result[drug_result$target=="HDAC",]
aaa<-c("S1848","S2759","S1194","S1047")
bbb<-c("5-FU","GEM","IRI","OXA","PAC")
ccc<-c(aaa,bbb)
HDAC2<-drug_result[drug_result$drug_id%in%ccc,]
HDAC3<-rbind(HDAC1,HDAC2)
HDAC<-HDAC3[order(HDAC3[,1],decreasing = T),]
HDAC_ic50<-ic50[match(as.character(HDAC[,1]),rownames(ic50)),]
library(pheatmap)
a=cor(t(HDAC_ic50))
result=pheatmap(a,scale = "none",main = "The correlation coefficients between HDACi and chemotherapy' AUC",show_rownames=T,show_colnames=T,
                clustering_distance_rows = "correlation",
                clustering_distance_cols = "correlation",clustering_method = "complete",cutree_rows=2,cutree_cols=2)

###########  HDACi20种药物过滤不相关的2个药物
clusters<-cutree(result$tree_row,k=3)
#a=table(clusters)
ic50_filter=HDAC_ic50[clusters==1,]

pheatmap(ic50_filter,scale = "none",main = "The AUC of HDACi is not normalized",show_rownames=T,show_colnames=T,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",clustering_method = "ward.D2")

bb<-t(apply(ic50_filter,1,function(x) scale(x)))
colnames(bb)<-colnames(ic50_filter)
norma_result<-pheatmap(bb,scale = "none",main = "The AUC of HDACi and chemotherapy are normalized",show_rownames=T,show_colnames=T,
                       clustering_distance_rows = "correlation",
                       clustering_distance_cols = "euclidean",clustering_method = "ward.D2")

bb<-t(apply(HDAC_ic50,1,function(x) scale(x)))
colnames(bb)<-colnames(HDAC_ic50)
norma_result<-pheatmap(bb,scale = "none",main = "The AUC of HDACi and chemotherapy are normalized",show_rownames=T,show_colnames=T,
                       clustering_distance_rows = "correlation",
                       clustering_distance_cols = "euclidean",clustering_method = "ward.D2")


###########
norma_result<-pheatmap(t(scale(t(ic50_filter))),scale = "row",main = "The AUC of HDACi is normalized",show_rownames=T,show_colnames=T,
                       clustering_distance_rows = "euclidean",
                       clustering_distance_cols = "euclidean",clustering_method = "ward.D2")

clusters_filter<-cutree(norma_result$tree_col,k=2)
sensitivity=colnames(ic50_filter[,clusters_filter==1])
resistance=colnames(ic50_filter[,clusters_filter==2])
resistance1=c("PC.104","PC.105","PC.5","PC.56")
medianl<-resistance[-match(resistance1,resistance)]


#########   clinical
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
resistance<-c("PC.117","PC.119","PC.121","PC.134","PC.135","PC.14","PC.16","PC.2","PC.27","PC.40","PC.78","PC.8","PC.97","PC.98","PC.G")
sensitivity<-c("PC.101","PC.102","PC.104","PC.105","PC.109","PC.111","PC.112","PC.115","PC.116","PC.13","PC.130","PC.136","PC.139","PC.18","PC.22","PC.5","PC.52","PC.56","PC.64","PC.81","PC.I","PC.L")

sensitivity<-c("PC.101","PC.109","PC.115","PC.116","PC.117","PC.13","PC.130","PC.136","PC.139","PC.16","PC.18","PC.2","PC.27","PC.64","PC.78","PC.8","PC.81","PC.97","PC.98","PC.L" )
resistance<-c("PC.102","PC.104","PC.105","PC.111","PC.112","PC.119","PC.121","PC.134","PC.135","PC.14","PC.22","PC.40","PC.5","PC.52","PC.56","PC.G","PC.I" )
resistance_new<-as.character(sample_id[match(resistance,sample_id[,2]),1])
sensitivity_new<-as.character(sample_id[match(sensitivity,sample_id[,2]),1])
resistance_cli<-clinical[match(intersect(clinical[,1],resistance_new),clinical[,1]),]
sensitivity_cli<-clinical[match(intersect(clinical[,1],sensitivity_new),clinical[,1]),]
all<-rbind(resistance_cli,sensitivity_cli)


  class_ind<-matrix(0,nrow(resistance_cli)+nrow(sensitivity_cli),1)
  class_ind[1:nrow(resistance_cli),1]<-"resistance"   #基因表达更高，越敏感，理论上生存越好
  class_ind[(nrow(resistance_cli)+1):nrow(class_ind),1]<-"sensitivity" #基因表达更低，越耐药，理论上生存越差
  
  label<-as.matrix(class_ind)
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
  
##2
  rm(list=ls())
  setwd("~/xjj/drug")
  ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
  ic50<-ic501[,!(colnames(ic501)%in%c("PC.100","PC.34"))]
  
  drug_info<-read.csv("compiled_drug_info.csv",header = TRUE,sep = ",")
  drug_result<-drug_info[match(rownames(ic50),drug_info$drug_id),]
  HDAC1<-drug_result[drug_result$target=="HDAC",]
  aaa<-c("S1848","S2759","S1194","S1047")
  bbb<-c("5-FU","GEM","IRI","OXA","PAC")
  ccc<-c(aaa,bbb)
  HDAC2<-drug_result[drug_result$drug_id%in%ccc,]
  HDAC3<-rbind(HDAC1,HDAC2)
  HDAC<-HDAC3[order(HDAC3[,1],decreasing = T),]
  HDAC_ic50<-ic50[match(as.character(HDAC[,1]),rownames(ic50)),]
  library(pheatmap)
  a=cor(t(HDAC_ic50))
  result=pheatmap(a,scale = "none",main = "The correlation coefficients between HDACi and chemotherapy' AUC",show_rownames=T,show_colnames=T,
                  clustering_distance_rows = "correlation",
                  clustering_distance_cols = "correlation",clustering_method = "complete",cutree_rows=2,cutree_cols=2)
  
  ##########   20种HDAC+5种化疗药
  bb<-t(apply(HDAC_ic50,1,function(x) scale(x)))
  colnames(bb)<-colnames(HDAC_ic50)
  norma_result<-pheatmap(bb,scale = "none",main = "The AUC of HDACi and chemotherapy are normalized",show_rownames=T,show_colnames=T,
                         clustering_distance_rows = "correlation",
                         clustering_distance_cols = "euclidean",clustering_method = "ward.D2")
  
  
  clusters_filter<-cutree(norma_result$tree_col,k=2)
  sensitivity=colnames(HDAC_ic50[,clusters_filter==1])
  resistance1=colnames(HDAC_ic50[,clusters_filter==2])
  resistance=c("PC.104","PC.105","PC.56")
  medianl<-resistance1[-match(resistance,resistance1)]
  
  
  setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/7_clinical")
  clinical<-read.table("clinical.txt",header = TRUE,sep = "\t")
  setwd("~/xjj/drug")
  sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
  library("data.table")
  resistance_new<-as.character(sample_id[match(resistance,sample_id[,2]),1])
  sensitivity_new<-as.character(sample_id[match(sensitivity,sample_id[,2]),1])
  medianl_new<-as.character(sample_id[match(medianl,sample_id[,2]),1])
  
  resistance_cli<-clinical[match(intersect(clinical[,1],resistance_new),clinical[,1]),]
  sensitivity_cli<-clinical[match(intersect(clinical[,1],sensitivity_new),clinical[,1]),]
  medianl_cli<-clinical[match(intersect(clinical[,1],medianl_new),clinical[,1]),]
  
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

  
### NC文章
  setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/7_clinical")
  clinical<-read.table("31patient_clinical.txt",header = TRUE,sep = "\t")
  
  
  resistance_cli<-clinical[which(clinical[,2]=="Resistant"),]
  sensitivity_cli<-clinical[which(clinical[,2]=="Sensitive"),]
  medianl_cli<-clinical[which(clinical[,2]=="Intermediate"),]
  
  all<-rbind(resistance_cli,sensitivity_cli,medianl_cli)
  label<-as.matrix(c(rep("resistance",nrow(resistance_cli)),rep("sensitivity",nrow(sensitivity_cli)),rep("medianl",nrow(medianl_cli))),ncol=1)
  TTP=as.numeric(all[,4])  ##生存时间
  status_TTP=as.matrix(all[,3])  ##生存状态
  status_TTP[which(status_TTP=="Disease Free")]=0
  status_TTP[which(status_TTP=="Recurred")]=1
  status_TTP<-as.numeric(status_TTP)
  surv_info1<-as.data.frame(cbind(TTP,status_TTP))
  
  surv_TTP<-survfit(Surv(TTP, status_TTP) ~ label,data=surv_info1)
  surv_TTP<-coxph(Surv(TTP, status_TTP) ~ label,data=surv_info1)
  ggsurvplot(surv_TTP,
             pval = TRUE, conf.int = TRUE, #是否显示p值/显示风险区域
             risk.table = TRUE,# 将风险表显示在生存曲线下面
             ggtheme = theme_bw(),
             tables.theme = theme_cleantable(),
             title="HDACIs combined with Chemotherapeutic drugs")


###单药
rm(list=ls())
setwd("~/xjj/drug/drug_result/drug64_AUC")
tresult_want<-read.table("drug64_sample.txt",sep="\t",header=T)
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/7_clinical")
clinical<-read.table("clinical.txt",header = TRUE,sep = "\t")

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

######### application
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/7_clinical")
clinical<-read.table("clinical.txt",header = TRUE,sep = "\t")
setwd("~/xjj/drug/drug_result/drug64_AUC")
tresult_want<-read.table("drug64_sample.txt",sep="\t",header=T)
setwd("~/xjj/drug")
sample_id<-read.csv("sample_id.csv",header = TRUE,sep = ",")
library(dplyr)
library("data.table")
library(foreach)
clinical[,4]<-sample_id[match(clinical[,1],sample_id[,1]),2]

setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/6_application/data")
peak_RPKM_name<-read.table("sample37_peaks_RPKM.txt",sep="\t",header=T)

result_final<-NULL
foreach(k=1:nrow(tresult_want),.combine='rbind')%do%{
  file_name=as.character(tresult_want[k,1])
  cat(file_name,"\n")
  common_sample<-as.data.frame(tresult_want[match(file_name,tresult_want[,1]),])
  sensitivity<-unlist(strsplit(as.character(common_sample[1,5]), ","))
  resistance<-unlist(strsplit(as.character(common_sample[1,4]), ","))
  sen_int<-intersect(sensitivity,clinical[,4])
  res_int<-intersect(resistance,clinical[,4])
  sensitivity_cli<-clinical[match(sen_int,clinical[,4]),] ##sensitivity samples recurrence state
  resistance_cli<-clinical[match(res_int,clinical[,4]),]
  sen_peak<-peak_RPKM_name[,match(sen_int,colnames(peak_RPKM_name))]
  res_peak<-peak_RPKM_name[,match(res_int,colnames(peak_RPKM_name))]
  peak_name<-peak_RPKM_name[,1]

  setwd(paste("~/xjj/drug/drug_result/drug64_AUC/",file_name,"/2_ATACseq_DApeaks",sep=""))
  res<-fread(paste("sensitivity-resistance-all-DESeq2_",file_name,"_ATAC.txt",sep=''),header=T,data.table=F)
  resSig<-res[order(format(res$pvalue,scientific = FALSE))[1:100],1]   ##features
  peak<-cbind(sen_peak,res_peak)
  state<-cbind(rbind(sensitivity_cli,resistance_cli),matrix(c(rep("sensitivity",nrow(sensitivity_cli)),rep("resistance",nrow(resistance_cli))),ncol=1))
  feature_peak<-peak[match(resSig,peak_name),]
  
  result<-NULL
  resultt<-NULL
  for(j in 1:ncol(feature_peak)){
    a<-feature_peak[,j]
    state_one<-state[match(colnames(feature_peak)[j],state[,4]),2]
    label_one<-state[match(colnames(feature_peak)[j],state[,4]),5]
    name_one<-colnames(feature_peak)[j]
    result1<-NULL
    for(m in 1:ncol(feature_peak)){
      b<-feature_peak[,m]
      name_two<-colnames(feature_peak)[m]
      r<-cor.test(a,b,method ="pearson",alternative="two.sided")$estimate
      p<-format(cor.test(a,b,method ="pearson",alternative="two.sided")$p.value,scientific = FALSE)
      state_two<-state[match(colnames(feature_peak)[m],state[,4]),2]
      label_two<-state[match(colnames(feature_peak)[m],state[,4]),5]
      result2<-matrix(c(colnames(feature_peak)[j],as.character(state_one),as.character(label_one),colnames(feature_peak)[m],as.character(state_two),as.character(label_two),r,p),nrow=1)
      result1<-rbind(result1,result2)
    }
    result3<-result1[which(result1[,8]<0.05),]
    if(class(result3)=="matrix"){
    result4<-as.matrix(result3[-which(result3[,4]==name_one),])
    result5<-result4[order(as.numeric(result4[,7]),decreasing=TRUE)[1:5],]
    state_sum<-sum(result5[,5]==state_one)
    label_sum<-sum(result5[,6]==label_one)
    result_sum<-as.data.frame(matrix(c(file_name,name_one,as.character(state_one),as.numeric(state_sum),as.numeric(state_sum/5),as.character(label_one),as.numeric(label_sum),as.numeric(label_sum/5)),nrow=1))
    result<-rbind(result,result_sum)
    resultt<-rbind(resultt,result1)
    }
  }
  colnames(resultt)<-c("sample1","state1","lable1","sample2","state2","lable2","r","p")
  colnames(result)<-c("Drug","sample","state","state_around","state_ratio","label","label_around","label_ratio")
  result_final<-rbind(result_final,result)
}  
setwd("~/xjj/drug/drug_result/HDAC20_chemotherapy5/6_application")
write.table(result_final,"result_final_500.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出
result_final<-read.table("result_final_500.txt",sep="\t",header=T)
class(result_final[,4])
result_final[which(as.numeric(result_final[,4])>=3),9]<-"Yes"
result_final[which(as.numeric(result_final[,4])<3),9]<-"No"

result111<-NULL
foreach(k=1:nrow(tresult_want),.combine='rbind')%do%{
  file_name=as.character(tresult_want[k,1])
  cat(file_name,"\n")
  yes_sum<-sum(result_final[result_final[,1]%in%file_name,9]=="Yes")
  no_sum<-sum(result_final[result_final[,1]%in%file_name,9]=="No")
  each_drug<-matrix(c(file_name,file_name,"Yes","No",yes_sum,no_sum),nrow=2)
  result111<-rbind(result111,each_drug)
}
colnames(result111)<-c("group","condictions","values")
result112<-as.data.frame(result111)
A<-result112[-which(is.na(result112[,3])),]
library(ggplot2)
library(forcats)
A$condictions <- factor(A$condictions,levels=c("Yes","No"))
#A$condictions <- fct_inorder(A$condictions)
A$group <- as.factor(A$group)
#A$group <- fct_inorder(A$group)
A$values<-as.numeric(as.character(A$values))

ggplot(A, aes(fill=condictions, y=values, x=group))+
  geom_bar(position=position_dodge(),stat="summary",colour = "black",size=0.5)+
  theme_classic(base_size = 12)+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.line = element_line(colour = "#000000",size = 0.5),
        axis.text = element_text(colour = "#000000" ,size = 10),
        axis.text.x = element_text(angle = 45,vjust = 0.85,hjust = 0.75), ##就是这里
        axis.ticks = element_line(colour = "#000000" ,size = 1) ,
        axis.ticks.length = unit(1,'mm'),
        plot.margin = unit(c(0.05,0,0,0),"cm"),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

setwd("~/xjj/drug")
ic501<-read.csv("changhai_auc.csv",header = TRUE,sep = ",",row.names=1)
ic50<-ic501[,!(colnames(ic501)%in%c("PC.100","PC.34"))]

drug_info<-read.csv("compiled_drug_info.csv",header = TRUE,sep = ",")
drug_result<-drug_info[match(rownames(ic50),drug_info$drug_id),]
HDAC1<-drug_result[drug_result$target=="HDAC",]
aaa<-c("S1848","S2759","S1194","S1047")
bbb<-c("5-FU","GEM","IRI","OXA","PAC")
ccc<-c(aaa,bbb)
HDAC2<-drug_result[drug_result$drug_id%in%ccc,]
HDAC3<-rbind(HDAC1,HDAC2)
HDAC<-HDAC3[order(HDAC3[,1],decreasing = T),]
B<-A[A[,1]%in%HDAC[,1],]
B$condictions <- factor(B$condictions,levels=c("Yes","No"))
B$group <- as.factor(B$group)
B$values<-as.numeric(as.character(B$values))
ggplot(B, aes(fill=condictions, y=values, x=group))+
  geom_bar(position=position_dodge(),stat="summary",colour = "black",size=0.5)+
  theme_classic(base_size = 12)+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.line = element_line(colour = "#000000",size = 0.5),
        axis.text = element_text(colour = "#000000" ,size = 10),
        axis.text.x = element_text(angle = 45,vjust = 0.85,hjust = 0.75), ##就是这里
        axis.ticks = element_line(colour = "#000000" ,size = 1) ,
        axis.ticks.length = unit(1,'mm'),
        plot.margin = unit(c(0.05,0,0,0),"cm"),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
