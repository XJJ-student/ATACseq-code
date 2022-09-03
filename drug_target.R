#####  drug-target
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDAC_drug/HDAC13_heatmap/cancer_modules_gene")
cancer_modules_genes<-read.table("heatmap13HDAC_P0.01_up_anno_DEGs_GSEA_cancer_modules.txt",header=T,sep='\t')
nrow(cancer_modules_genes)
#symbol_id<-read.table("2id.txt",header=T,sep='\t')
#cancer_modules_genes<-symbol_id[match(cancer_modules_genes1[,1],symbol_id[,2]),1]
setwd("~/xjj/DESeq_DEGs/DESeq2_P05/DEGs_intersect_Diffpeakanno_gene")
library("data.table")
interactions<-fread("interactions.tsv",header=T,data.table=F)
int_gene<-intersect(cancer_modules_genes[,1],interactions[,1])
int_gene_drug1<-interactions[interactions[,1]%in%int_gene,]
int_gene_drug2<-int_gene_drug1[order(int_gene_drug1[,1],decreasing =F),]
int_gene_time<-cancer_modules_genes[match(int_gene,cancer_modules_genes[,1]),]

gene_drug<-NULL
a<-unique(int_gene_drug2[,1])
for(i in 1:length(a)){
  max_positon1<-int_gene_drug2[which(int_gene_drug2[,1]==a[i]),]
  max_positon2<-max_positon1[order(max_positon1[,10],decreasing =T),]
  max_positon3<-int_gene_time[which(int_gene_time[,1]==a[i]),2]
  max_positon4<-rep(max_positon3,nrow(max_positon2))
  max_positon5<-cbind(max_positon4,max_positon2)
  gene_drug1<-max_positon5[which(!is.na(max_positon5[,11])),]
  gene_drug<-rbind(gene_drug,gene_drug1)
}
name<-data.frame(rep(toupper("HDAC"),nrow(gene_drug)))
gene_drug_final<-cbind(name,gene_drug)
colnames(gene_drug_final)<-c("drug","times",colnames(int_gene_drug2))
gene_drug_final1<-gene_drug_final[order(gene_drug_final[,2],decreasing =T),]
gene_drug_final2<-gene_drug_final1[,c(1,2,3,10,12,4,5,6,7,8,9,11,13)]
#want_target<-gene_drug_final2[which(toupper(gene_drug_final2$drug_name)==toupper("Paclitaxel")),c(3,4)]
#want_target
length(unique(gene_drug_final2$gene_name))
setwd("~/xjj/drug/drug_result/HDAC_drug/HDAC13_heatmap/cancer_associated_gene_drug_target")
write.table(gene_drug_final2,"heatmap13HDAC_P0.01_up_anno_DEGs_GSEA_cancer_modules_drug_target.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)

type<-unique(gene_drug_final2$interaction_types)
genename<-unique(gene_drug_final2$gene_name)
gene_drug_final_type<-data.frame(gene_drug_final2[match(genename,gene_drug_final2$gene_name),c(1:3)])
for(i in 1:length(genename)){
  onegene<-gene_drug_final2[(gene_drug_final2$gene_name) %in% genename[i],]
  drugtype1<-paste0(onegene$interaction_types,collapse = ",")
  gene_drug_final_type[i,4]<-drugtype1
  gene_drug_final_type[i,5]<-nrow(onegene)
}
colnames(gene_drug_final_type)<-c("drug","time","gene_name","drug_type","drug_number")
gene_drug_final_type<-gene_drug_final_type[order(gene_drug_final_type[,5],decreasing =T),]

for(i in 1:nrow(gene_drug_final_type)){
  type<-unlist(strsplit(as.character(gene_drug_final_type[i,4]), ","))
  uniquetype<-unique(type)
  group<-c()
  for(j in 1:length(uniquetype)){
    count<-sum(type==uniquetype[j])
    group1<-paste(uniquetype[j],count,sep = "-")
    group<-c(group,group1)
  }
  gene_drug_final_type[i,6]<-paste0(group,collapse = ",")
}
colnames(gene_drug_final_type)<-c("drug","time","gene_name","drug_type","drug_number","drug_type_number")
setwd("~/xjj/drug/drug_result/HDAC_drug/HDAC13_heatmap/cancer_associated_gene_drug_target")
write.table(gene_drug_final_type,"heatmap13HDAC_P0.01_up_anno_DEGs_GSEA_cancer_modules_drug_target_summary.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE)


################# TF-drug
rm(list=ls())
setwd("~/xjj/DESeq_DEGs/DESeq2_P05/DEGs_intersect_Diffpeakanno_gene")
library("data.table")
interactions<-fread("interactions.tsv",header=T,data.table=F)

TF_target<-interactions[which(toupper(interactions$gene_name)=="ELF5"),]
TF_drug1<-TF_target[which(!is.na(TF_target[,10])),]
library(dplyr)
TF_drug<-TF_drug1 %>% distinct(gene_name,drug_name, .keep_all = TRUE)
nrow(TF_drug)
TF_drug[which(toupper(TF_drug$drug_name)==toupper("Oxaliplatin")),]


########   homer-TF-drug
setwd("~/xjj/drug/drug_result/HDAC_drug/S2693/homer_result/S2693-30-p-0.05-logFC0-DESeq2-up")
knownResults<-read.delim("knownResults.txt",header=T,sep='\t')
FDR05<-knownResults[knownResults$q.value..Benjamini.<0.05,]
library(dplyr)
homer_TF1<-data.frame(lapply(strsplit(as.character(FDR05[,1]),'/'), function(x) x[3])%>%unlist())
FDR051<-FDR05[which(homer_TF1=="Homer"),]
homer_TF2<-data.frame(lapply(strsplit(as.character(FDR051[,1]),'/'), function(x) x[1])%>%unlist())
library(stringr)
homer_TF3<-apply(homer_TF2,1,function(x) str_replace_all(as.character(x),"\\((.*?)\\)",""))##取出括号前的字符
homer_TF<-toupper(homer_TF3)
length(homer_TF)


setwd("~/xjj/DESeq_DEGs/DESeq2_P05/DEGs_intersect_Diffpeakanno_gene")
library("data.table")
interactions<-fread("interactions.tsv",header=T,data.table=F)
int_gene<-intersect(interactions[,1],homer_TF)
length(int_gene)
int_gene_drug1<-interactions[interactions[,1]%in%int_gene,]
int_gene_drug2<-int_gene_drug1[which(!is.na(int_gene_drug1[,10])),]
int_gene_drug3<-int_gene_drug2[order(int_gene_drug2[,1],decreasing =F),]

library(dplyr)
TF_drug<-int_gene_drug3 %>% distinct(gene_name,drug_name, .keep_all = TRUE)
nrow(TF_drug)
TF_drug[which(toupper(TF_drug$drug_name)==toupper("Resminostat")),c(1,8)]
as.matrix(table(TF_drug$gene_name))
length(unique(TF_drug$gene_name))

which(toupper(interactions$drug_name)==toupper("CAY10603"))

#################################################################
#####  pathway-drug-target
rm(list=ls())
setwd("~/xjj/drug/drug_result/HDAC_drug/HDAC13_heatmap/pathway_gene")
cancer_modules_genes1<-read.table("kegg_pathway_P0.05_up_gene.txt",header=F,sep='\t')
nrow(cancer_modules_genes1)
symbol_id<-read.table("2id.txt",header=T,sep='\t')
HDAC_drugname<-read.table("HDAC_drugname.txt",header=F,sep='\t')
HDAC<-read.table("HDAC.txt",header=T,sep='\t')

cancer_modules_genes<-as.data.frame(symbol_id[match(cancer_modules_genes1[,1],symbol_id[,2]),1])
setwd("~/xjj/DESeq_DEGs/DESeq2_P05/DEGs_intersect_Diffpeakanno_gene")
library("data.table")
interactions<-fread("interactions.tsv",header=T,data.table=F)
int_gene<-intersect(cancer_modules_genes[,1],interactions[,1])
int_gene_drug1<-interactions[interactions[,1]%in%int_gene,]
int_gene_drug2<-int_gene_drug1[order(int_gene_drug1[,1],decreasing =F),]
int_gene_drug3<-int_gene_drug2[which(!is.na(int_gene_drug2[,10])),]

aaa<-intersect(toupper(int_gene_drug3$drug_name),toupper(HDAC_drugname[,1]))
unique(int_gene_drug3[toupper(int_gene_drug3$drug_name)%in%aaa,1])

HDAC[match(aaa,toupper(HDAC$drug_name)),2]







