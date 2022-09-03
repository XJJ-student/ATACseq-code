setwd("F:\\Organoid\\突变数据")
mut_original<-read.csv("mut.csv",header = TRUE,sep = ",")
a<-unique(mut_original$Tumor_Sample_Barcode)
b<-unique(mut_original$Hugo_Symbol)
result<-matrix(0,length(b),length(a))
for(i in 1:length(a)){
  want_sample<-mut_original[mut_original$Tumor_Sample_Barcode%in%a[i],]
  result[match(want_sample$Hugo_Symbol,b),i]<-1
}
colnames(result)<-a
rownames(result)<-b
write.table(result,"Organoid_mut_binary.txt",sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE);


###read.maf,maftools包将无意义突变给滤掉了
BiocManager::install("maftools")
library(maftools)
setwd("F:\\Organoid\\突变数据")
laml<-read.maf(maf="WGS.maf")
write.table(laml@data,"WGS_mut_all_gene.txt",sep="\t",row.names=F,col.names=TRUE)


### 读取最初的没有过滤的所有数据,csv文件和xlsx文件的数据是一样的，使用csv文件进行输入
setwd("F:\\Organoid\\突变数据")
mut_original<-read.csv("mut_all_gene.csv",header = TRUE,sep = ",")
b<-as.matrix(unique(mut_original$Hugo_Symbol))
write.table(b,"WGS_mut_all_gene.txt",sep="\t",row.names=F,col.names=TRUE)


###  oncoplot
library(maftools)
setwd("F:\\Organoid\\突变数据")
var_maf<-read.maf(maf="WGS.maf")
plotmafSummary(maf = var_maf, rmOutlier = TRUE, addStat = 'median',dashboard = TRUE,titvRaw = FALSE)
oncoplot(maf = var_maf, top = 40, writeMatrix=T,removeNonMutated = F,showTumorSampleBarcodes=T)

#显示特定基因
oncostrip(maf = var_maf, genes = c('KRAS','TP53', 'VSX1'), writeMatrix=T,removeNonMutated = F,showTumorSampleBarcodes=T)

#titv函数将SNP分类为Transitions_vs_Transversions，并以各种方式返回汇总表的列表。 汇总数据也可以显示为一个箱线图，显示六种不同转换的总体分布，并作为堆积条形图显示每个样本中的转换比例。
titv(var_maf, useSyn = FALSE, plot = TRUE, file = NULL)

#Hotspot mutation,TAD——拓扑相关结构域
lollipopPlot
lollipopPlot(maf = var_maf, gene = 'TP53', AACol = 'AAChange.refGene')
lollipopPlot(maf = var_maf, gene = 'KRAS', AACol = 'AAChange.refGene')


#共突变和突变互斥
somaticInteractions(maf = var_maf, top = 25, genes = NULL,pvalue = c(0.05, 0.1),
                    returnAll = TRUE,geneOrder = NULL,fontSize = 0.6,showSigSymbols = TRUE,
                    showCounts = FALSE,countStats = "all",countType = "all",countsFontSize = 0.5,
                    countsFontColor = "black",colPal = "BrBG",showSum = TRUE,
                    colNC = 6,nShiftSymbols = 5,sigSymbolsSize = 2,sigSymbolsFontSize = 0.9,
                    pvSymbols = c(46, 42),limitColorBreaks = TRUE)


##计算TMB
library(maftools)
setwd("F:\\Organoid\\突变数据")
var_maf<-read.maf(maf="WGS.maf")
x = tmb(maf = var_maf)#not log (logScale = F)
head(x)
quantile(x$total_perMB)
write.table(x,"WGS_mut_TMB.txt",sep="\t",row.names=TRUE,col.names=TRUE)

#使用sequenza软件判定肿瘤纯度

#判断拷贝数是扩增还是删失的原则
#TCN_EM	Integer	Total copy number. 2 for normal diploid，LCN_EM	Integer	Lesser (minor) copy number. 1 for normal diploid
#TCN_EM要是不等于2就一定是拷贝数变异，只有当TCN_EM==2和LCN_EM==1时才是正常的，要是LCN_EM==0时，则为LOH(杂合性缺失)
rm(list=ls())
install.packages("vcfR")
library(vcfR)
setwd("G:\\类器官\\WGS_CNV") 
first_category_name = list.files("CNV_FACETS")  
dir = paste("G:\\类器官\\WGS_CNV\\CNV_FACETS\\",first_category_name,sep="")  
n = length(dir) 
result<-NULL
for(i in 1:n){      
  setwd(dir[i])
  facets<-read.vcfR(paste(first_category_name[i],"facets.vcf.gz",sep='.'))
  a<-as.matrix(facets@fix)
  b<-as.matrix(a[,8])
  exp<-apply(b,1,function(x) unlist(strsplit(x, "[;]")))
  TCN<-matrix(exp[12,],ncol=1)
  LCN<-matrix(exp[13,],ncol=1)
  one_col<-matrix(rep(first_category_name[i],nrow(a)),ncol=1)
  result1<-cbind(one_col,a,TCN,LCN)
  result<-rbind(result,result1)
}
colnames(result)<-c("samples",colnames(a),"TCN_EM","LCN_EM")
setwd("G:\\类器官\\WGS_CNV\\CNV_FACETS")
write.table(result,"WGS_CNV_all.txt",sep="\t",row.names=FALSE,col.names=TRUE)

rm(list=ls())
setwd("G:\\类器官\\WGS_CNV\\CNV_FACETS")
exp<-read.table("WGS_CNV_all.txt",header=T,sep='\t',stringsAsFactors = FALSE)
no_cnv1<-which(exp[,10]=="TCN_EM=2")
no_cnv2<-which(exp[,11]=="LCN_EM=1")
no_cnv<-intersect(no_cnv1,no_cnv2)
exp[no_cnv,12]<-0
exp[-(no_cnv),12]<-1
write.table(exp,"WGS_CNV_all_01.txt",sep="\t",row.names=FALSE,col.names=TRUE)

BiocManager::install("GenomicFeatures")
library(GenomicFeatures)
setwd("G:\\类器官")
txdb<-makeTxDbFromGFF(file="hg19_gene.gff3",format="gff3")

#############
file=list.files(dir[i])
grep("vcf.gz", file, value = TRUE)
facets<-read.vcfR("file\\*.vcf.gz")
file[grepl("vcf.gz", file)]
##############

#######################
######  sequenza
rm(list=ls())
library(data.table)
setwd("~/176/Changhai_WGS/results/2_Variants/SomaticVcfs")
ENSG <- fread("Homo_sapiens_protein_coding.GRCh37.75.gtf")
ENSG_gene<-ENSG[which(ENSG[,3]=="gene"),]
ENSG1<-as.matrix(ENSG_gene)
setwd("~/176/Changhai_WGS/results/2_Variants/SomaticVcfs")
first_category_name = list.files("scarHRD")  
dir = paste("~/176/Changhai_WGS/results/2_Variants/SomaticVcfs/scarHRD/",first_category_name,sep="")  
n = length(dir)-1 

out<-foreach(i=1:n,.combine='rbind')%do%{
  setwd(dir[i])
  facets<-read.table(paste(first_category_name[i],"segments.txt"),sep="\t",header=T)
  a<-facets$CNt
  facets[which(facets$CNt==2),14]<-0
  facets[which(facets$CNt!=2),14]<-1
  gene_matrix<-NULL
  for(j in 1:nrow(facets)){
    d<-which(ENSG1[,1]==facets[j,1] & ENSG1[,4]>facets[j,2] & ENSG1[,5]<facets[j,3])
    f<-ENSG1[d,9]
    if(length(d)==0){
      gene_matrix1<-matrix(c(facets[j,],"null"),nrow=1)
      gene_matrix<-rbind(gene_matrix,gene_matrix1)
    }
    if(length(d)==1){
      gene_matrix1<-matrix(c(facets[j,],f),nrow=1)
      gene_matrix<-rbind(gene_matrix,gene_matrix1)
    }
    if(length(d)>1){
      gene_matrix1<-cbind(as.matrix(facets[rep(j,length(f)),]),f)
      gene_matrix<-rbind(gene_matrix,gene_matrix1)
    }
  }
  one_col<-matrix(rep(first_category_name[i],nrow(gene_matrix)),ncol=1)
  result<-data.frame(one_col,gene_matrix)
  return(result)
}
setwd("~/xjj/WGS_CNV/CNV_sequenza")
save(out,  file = "CNV_sequenza_gene.RData")

out1<-out[-(which(out[,16]=="null")),]
out2<-out1[which(out1[,15]==1),]
CNV_gene_origial<-out2[,c(1,15,16)]
colnames(CNV_gene_origial)<-c("sample","CNV","gene")
b<-as.matrix(CNV_gene_origial[,3])
exp<-apply(b,1,function(x) unlist(strsplit(as.character(x), "[;]")))
exp1<-data.frame(matrix(unlist(exp), nrow=nrow(CNV_gene_origial), byrow=T),stringsAsFactors=FALSE)
library(dplyr)
gene<-data.frame(lapply(strsplit(exp1[,2],'\"'), function(x) x[2])%>%unlist()) ## get gene_names
CNV_gene_origial[,4]<-gene
a<-unique(CNV_gene_origial[,1])
b<-unique(CNV_gene_origial[,4])
result<-matrix(0,length(b),length(a))
for(i in 1:length(a)){
  want_sample<-CNV_gene_origial[CNV_gene_origial[,1]%in%a[i],]
  result[match(want_sample[,4],b),i]<-1
}
colnames(result)<-a
rownames(result)<-b
write.table(result,"CNV_sequenza_binary.txt",sep="\t",row.names=TRUE,col.names=TRUE)


#######################    简化版本 ensamble数据库的参考基因组
######  sequenza
rm(list=ls())
library(data.table)
setwd("~/xjj/Reference genome")
ENSG <- fread("Homo_sapiens.GRCh37.75.gtf")
ENSG_gene<-ENSG[which(ENSG[,3]=="gene"),]
ENSG2<-apply(ENSG_gene[,1],1,function(x) paste("chr",x,sep=''))
ENSG_gene[,1]<-ENSG2
ENSG1<-as.matrix(ENSG_gene)
setwd("~/176/Changhai_WGS/results/2_Variants/SomaticVcfs")
first_category_name = list.files("scarHRD")  
dir = paste("~/176/Changhai_WGS/results/2_Variants/SomaticVcfs/scarHRD/",first_category_name,sep="")  
n = length(dir)-1 

out<-foreach(i=1:n,.combine='rbind')%do%{
  setwd(dir[i])
  facets1<-read.table(paste(first_category_name[i],"segments.txt",sep='_'),sep="\t",header=T)
  facets<-facets1[-(which(facets1$CNt==2 & facets1$A==1 & facets1$B==1)),]
  gene_matrix<-NULL
  for(j in 1:nrow(facets)){
    d<-which(ENSG1[,1]==facets[j,1] & ENSG1[,4]>facets[j,2] & ENSG1[,5]<facets[j,3])
    f<-ENSG1[d,9]
    if(length(d)==0){
      gene_matrix1<-matrix(c(facets[j,],"null"),nrow=1)
      gene_matrix<-rbind(gene_matrix,gene_matrix1)
    }
    if(length(d)==1){
      gene_matrix1<-matrix(c(facets[j,],f),nrow=1)
      gene_matrix<-rbind(gene_matrix,gene_matrix1)
    }
    if(length(d)>1){
      gene_matrix1<-cbind(as.matrix(facets[rep(j,nrow(f)),]),f)
      gene_matrix<-rbind(gene_matrix,gene_matrix1)
    }
  }
  one_col<-matrix(rep(first_category_name[i],nrow(gene_matrix)),ncol=1)
  result<-data.frame(one_col,gene_matrix)
  return(result)
}
setwd("~/xjj/WGS_CNV/CNV_sequenza")
save(out,  file = "CNV_sequenza_gene_newhg19.RData")

out2<-out[-(which(out[,15]=="null")),]  #facets$CNt!=2  CNn = 2 number of alleles in the normal sample
CNV_gene_origial<-out2[,c(1,14,15)]
b<-as.matrix(CNV_gene_origial[,3])
exp<-apply(b,1,function(x) unlist(strsplit(as.character(x), "[;]")))
exp1<-data.frame(matrix(unlist(exp), nrow=nrow(CNV_gene_origial), byrow=T),stringsAsFactors=FALSE)
library(dplyr)
gene<-data.frame(lapply(strsplit(exp1[,2],'\"'), function(x) x[2])%>%unlist()) ## get gene_names
b<-as.matrix(gene)
gene1<-apply(b,1,function(x) unlist(strsplit(as.character(x), "-"))[1])
gene2<-data.frame(matrix(unlist(gene1), nrow=nrow(CNV_gene_origial), byrow=F),stringsAsFactors=FALSE)
w=gene2[1:nrow(gene2),]
gene3<-sub("(\\.\\w+)+","",w)

CNV_gene_origial[,4]<-gene3
a<-unique(CNV_gene_origial[,1])
b<-unique(CNV_gene_origial[,4])
result<-matrix(0,length(b),length(a))
for(i in 1:length(a)){
  want_sample<-CNV_gene_origial[CNV_gene_origial[,1]%in%a[i],]
  result[match(want_sample[,4],b),i]<-1
}
colnames(result)<-a
rownames(result)<-b
write.table(result,"CNV_sequenza_binary_easily_newhg19.txt",sep="\t",row.names=TRUE,col.names=TRUE)




#######################    简化版本 ANNOVAR .gtf文件
######  sequenza 
rm(list=ls())
library(data.table)
setwd("~/176/Changhai_WGS")
ENSG <- fread("hg19.gtf")
colnames(ENSG)<-c("chr","dev","region","start","end","no1","no2","no3","gene")
ENSG13<-ENSG[which(ENSG[,3]=="exon"),]
b<-as.matrix(ENSG13[,9])
exp<-apply(b,1,function(x) unlist(strsplit(as.character(x), "[;]")))
exp1<-data.frame(matrix(unlist(exp), nrow=nrow(ENSG13), byrow=T),stringsAsFactors=FALSE)
library(dplyr)
gene1<-data.frame(lapply(strsplit(exp1[,5],'\"'), function(x) x[2])%>%unlist()) ## get gene_names
b<-as.matrix(gene1)
gene2<-apply(b,1,function(x) unlist(strsplit(as.character(x), "-"))[1])
gene<-data.frame(matrix(unlist(gene2), nrow=length(gene2), byrow=F),stringsAsFactors=FALSE)
ENSG12<-cbind(ENSG13,gene)
colnames(ENSG12)<-c(colnames(ENSG13),"gene_alone")
library(dplyr, warn.conflicts = FALSE)
ENSG11 <- ENSG12 %>% distinct(chr,start,end,gene_alone, .keep_all = T)
ENSG1<-as.matrix(ENSG11)

setwd("~/176/Changhai_WGS/results/2_Variants/SomaticVcfs")
first_category_name = list.files("scarHRD")  
dir = paste("~/176/Changhai_WGS/results/2_Variants/SomaticVcfs/scarHRD/",first_category_name,sep="")  
n = length(dir)-1 
library(foreach)
out2<-foreach(i=51:99,.combine='rbind')%do%{
  setwd(dir[i])
  facets<-read.table(paste(first_category_name[i],"segments.txt",sep='_'),sep="\t",header=T)
  gene_matrix<-data.frame()
  for(j in 1:nrow(facets)){
    d<-which(ENSG1[,1]==facets[j,1] & ENSG1[,4]>facets[j,2] & ENSG1[,5]<facets[j,3])
    f<-ENSG11[d,c(3,10)]
    if(length(d)==0){
      gene_matrix1<-data.frame(c(facets[j,c(1,2,3,10,11,12)],c("null","null")))
      colnames(gene_matrix1)<-c(colnames(facets)[c(1,2,3,10,11,12)],"type","gene_alone")
      gene_matrix<-rbind(gene_matrix,gene_matrix1)
    }
    if(length(d)!=0){
      gene_matrix1<-data.frame(cbind(facets[rep(j,nrow(f)),c(1,2,3,10,11,12)]),f)
      colnames(gene_matrix1)<-c(colnames(facets)[c(1,2,3,10,11,12)],"type","gene_alone")
      gene_matrix<-rbind(gene_matrix,gene_matrix1)
    }
  }
  gene_matrix2 <- gene_matrix %>% distinct(chromosome,start.pos,end.pos,gene_alone, .keep_all = T)
  write.table(gene_matrix2,file=paste(first_category_name[i],"exon_gene.txt",sep='_'),sep="\t",row.names=FALSE,col.names=TRUE)
  one_col<-matrix(rep(first_category_name[i],nrow(gene_matrix2)),ncol=1)
  result<-data.frame(one_col,gene_matrix2)
  return(result)
}
out111<-rbind(out1,out2)
out<-rbind(out111,out3)
setwd("~/xjj/WGS_CNV/CNV_sequenza")
save(out,  file = "CNV_sequenza_gene_ANNOVAR_exon.RData")

outt<-out[-(which(out$CNt==2 & out$A==1 & out$B==1)),]
out2<-outt[-(which(outt$gene_alone=="null")),] 
CNV_gene_origial<-out2[,c(1,9)]####还没修改
a<-unique(CNV_gene_origial[,1])
b<-unique(CNV_gene_origial[,2])
result<-matrix(0,length(b),length(a))
for(i in 1:length(a)){
  want_sample<-CNV_gene_origial[CNV_gene_origial[,1]%in%a[i],]
  result[match(want_sample[,2],b),i]<-1
}
colnames(result)<-a
rownames(result)<-b
write.table(result,"CNV_sequenza_binary_ANNOVAR_exon.txt",sep="\t",row.names=TRUE,col.names=TRUE)

##### TMB数据处理
setwd("~/176/Changhai_WGS/results/2_Variants/SomaticVcfs")
first_category_name = list.files("strelka2_somatic_TMB")  
dir = paste("~/176/Changhai_WGS/results/2_Variants/SomaticVcfs/strelka2_somatic_TMB/",first_category_name,sep="")  
n = length(dir)-3
result<-data.frame()
for(i in 1:n){
    setwd(dir[i])
  facets<-read.table(paste(first_category_name[i],"TMB.results.txt",sep='.'),sep="\t",header=T)
  result<-rbind(result,facets)
}  
setwd("~/xjj/WGS_CNV/TMB_result")
write.table(result,"TMB_result.txt",sep="\t",row.names=TRUE,col.names=TRUE)


##### ATAC_stat










