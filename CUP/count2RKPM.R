
##计算外显子并集
len_sum=function(test4){
  #计算
  if(is.null(dim(test4)[1])){
    return(test4[2]-test4[1])
  }
  
  test4=test4[order(test4[,1]),]
  exon_len=0
  exon_range=test4[1,]
  for (i in 2:dim(test4)[1]) {
    if(test4[i,1]<test4[i-1,2]){
      exon_range=c(min(test4[i,1],test4[i-1,1]),max(test4[i,2],test4[i-1,2]))
    }else{
      exon_len=exon_len+(exon_range[2]-exon_range[1])
      exon_range=test4[i,]
    }
  }
  exon_len=exon_len+(exon_range[2]-exon_range[1])
  return(exon_len)
}

#计算外显子长度，并导出文件
rm(list=ls())
#setwd("~/xjj/RNAseq/gtf")
setwd("~/xjj/CUP/gtf")
library(data.table)
library(stringr)
ENSG <- fread("Homo_sapiens.GRCh37.75.gtf")
ENSG=ENSG[which(as.matrix(ENSG[,3])=="exon"),]

len=as.matrix(ENSG[,c(4,5)]) #后续计算外显子长度
ENSG=as.matrix(ENSG[,9])     #名称

#提取名称
ENSG2=apply(as.matrix(ENSG), 1,function(x) strsplit(str_extract(x,"gene_id \"(\\w)*\""),"\"")[[1]][2])#如果是基因名gene_name，gene_id

ENSG1=apply(as.matrix(ENSG), 1,function(x) strsplit(str_extract(x,"gene_name \"(\\w)*\""),"\"")[[1]][2])#如果是基因名gene_name，gene_id
len=len[!is.na(ENSG2),]
ENSG2=ENSG2[!is.na(ENSG2)]

len_all=c()
ENSG_uni=unique(ENSG2)
for (i in 1:length(ENSG_uni)) {
  len_all=c(len_all,len_sum(len[ENSG2==ENSG_uni[i],]))
}

ENSG1_uni=ENSG1[match(ENSG_uni,ENSG2)]
write.table(cbind(ENSG_uni,ENSG1_uni),"exon_name_37.75.txt",sep="\t",col.names = T,row.names = F)

write.table(cbind(ENSG1_uni,len_all),"exon_symbol_len_37.75.txt",sep="\t",col.names = T,row.names = F)
colnames(result)<-c("ENSG","genename")
result<-as.data.frame(result)
library(dplyr)
result1<-result %>% distinct(ENSG,genename, .keep_all = T)
#读取外显子长度文件，并处理文件
temp<-read.table("G:\\pre\\data\\ref\\exon_len_id_37.75.txt",sep="\t",header=F)
temp<-read.table("G:\\pre\\data\\ref\\exon_len_name_37.75.txt",sep="\t",header=F)
temp<-read.table("G:\\pre\\data\\ref\\exon_len_id_38.94.txt",sep="\t",header=F)
setwd("~/xjj/CUP/gtf")
temp<-read.table("exon_len_37.75.txt",sep="\t",header=T)
ENSG_uni=temp[,1]
len_all=temp[,2]
#(1000000*C)/(N*L/1000)
#设C 为比对到 gene A 的 reads数（read count），N为比对到所有 gene 的总 reads 数，L 为 gene A 的碱基数。
ENSGid<-read.table("exon_name_37.75.txt",sep="\t",header=T)
id<-as.character(ENSGid[match(intersect(ENSG_uni,ENSGid[,1]),ENSGid[,1]),2])
len_id<-len_all[match(intersect(ENSG_uni,ENSGid[,1]),ENSG_uni)]

setwd("~/176/Changhai_mRNA/2_Expression/Summary/star_rsem")
exp<-read.csv("GSE179979_Raw_counts.csv",sep=",",header=T)
exp<-read.table("star_rsem.GeneSymbol.readscount.xls",sep="\t",header=T)
#exp<-read.table("merge.txt",sep="\t",header=T)
gene = exp[,1]
exp = floor(exp[,-1])
#temp=intersect(as.matrix(gene),as.matrix(ENSG_uni))
temp=intersect(as.matrix(gene),as.matrix(id))

exp_0=exp[match(temp,as.matrix(gene)),]
len_all_0=len_id[match(temp,as.matrix(id))]
#len_all_0=len_all[match(temp,ENSG_uni)]

exp_RPKM=apply(exp_0,2,function(x) 1000000000*x/sum(x)/len_all_0)

write.table(cbind(as.matrix(temp),exp_RPKM),"gene_exp_RPKM_million.txt",col.names=T,row.names=F,sep="\t")
ENSG<-read.table("ENSG.txt",sep="\t",header=T)
id<-read.table("2id.txt",sep="\t",header=T)
int_gene<-intersect(ENSG[,2],id[,2])
id1<-id[match(int_gene,id[,2]),]
ENSG1<-ENSG[match(int_gene,ENSG[,2]),]
ENSG_id_symbol<-cbind(ENSG1,id1)

int_ENSG<-intersect(temp,ENSG[,1])
final_exp<-exp_RPKM[match(int_ENSG,temp),]
final_exp_geneid<-cbind(as.matrix(int_ENSG),final_exp)
colnames(final_exp_geneid)<-c("geneid",colnames(exp_RPKM))
write.table(final_exp_geneid,"geneid_exp_RPKM.txt",col.names=T,row.names=F,sep="\t")



probe_to_gene("probe_exp.txt","G:\\pre\\data\\GPL\\ENSG.txt","gene_exp.txt")#探针到基因
label<-read.table("label.txt",sep="\t",header=T)
label=label[,c(18,45)]
exp<-as.matrix(read.table("gene_exp.txt",sep="\t",header=T))
gene = exp[,1]
exp = exp[,-1]

exp=exp[,-match(label[143:159,1],unlist(lapply(strsplit(colnames(exp),".fastq"), function(x) x[1])))]

ENSG_uni=unlist(lapply(strsplit(colnames(exp),"."), function(x) x[1]))


############     TCGA所有的表达谱merge   #######################
setwd("~/xjj/CUP/GEO/SRP069243")
first_category_name = list.files("GEO_download")  
 
n = length(first_category_name)
library(foreach)
library("data.table")
setwd("~/xjj/CUP/GEO/SRP069243/GEO_download")
i=1
primary<-fread(first_category_name[i],header=F,data.table=F)
geneid<-primary[,1]

gene<-geneid
for(i in 1:(n-1)){
  i=i+1
  setwd("~/xjj/CUP/GEO/SRP069243/GEO_download")
  primary1<-fread(first_category_name[i],header=F,data.table=F)
  geneid1<-primary1[,1]
  gene<-intersect(gene,geneid1)
}
int_gene<-data.frame(gene)
setwd("~/xjj/CUP/GEO/SRP069243/GEO_download")
write.table(int_gene,"all_geneid.txt",sep = '\t',col.names = F,row.names = F,quote = FALSE)##数据输出

#############    TCGA exp    ###########################
rm(list=ls())
setwd("~/xjj/CUP/GEO/SRP069243")
int_gene<-as.matrix(read.delim("all_geneid.txt",sep = '\t',header = F))##数据输出
setwd("~/xjj/CUP/GEO/SRP069243")
first_category_name = list.files("GEO_download")  

n = length(first_category_name)
library(foreach)
library("data.table")
i=1
setwd("~/xjj/CUP/GEO/SRP069243/GEO_download")
primary<-fread(first_category_name[i],header=F,data.table=F)
exp<-primary[,-1]
gene<-as.matrix(primary[,1])
exp1<-primary[match(int_gene,primary[,1]),]

expp<-exp1
foreach(i=1:(n-1),.combine='cbind')%do%{
  setwd("~/xjj/CUP/GEO/SRP069243/GEO_download")
  i=i+1
  primary1<-fread(first_category_name[i],header=F,data.table=F)
  exp2<-primary1[match(int_gene,primary1[,1]),2]
  expp<-cbind(expp,exp2)
}
colnames(expp)<-c("geneid",first_category_name)
write.table(expp,"gene_exp.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)##数据输出


