setwd("/scratch/jibzhang/project/DHL_WES/")
gc(rm(list=ls()))
source("/home/zgu_labs/bin/R_script/functions.R")
library(dplyr)
library(reshape2)
#mutect2 mutations
file_list=system("ls /scratch/jibzhang/project/DHL_WES/*/EXOME/search_dbsnp_human_Mutect2noNormal/*.search_dbsnp.tsv",intern = T)
file_one=file_list[1]
data_mutation=bind_rows(lapply(file_list,function(file_one){
  data_one=read.table(file_one,sep="\t",stringsAsFactors = F,header = T)
  data_one$AF_gnomADGenome=as.character(data_one$AF_gnomADGenome)
  data_one$AF_gnomADExome=as.character(data_one$AF_gnomADExome)
  data_one$AF_dbsnpAll=as.character(data_one$AF_dbsnpAll)
  data_one$COH_sample=stringr::word(basename(file_one),1,sep="[.]")
  data_one
}))

data_mutation$new_var="database"
data_mutation1=data_mutation %>% select(names(data_mutation)[c(1:25,29)]) %>% distinct()

data_mutation_1=dcast(data_mutation,snp_i+COH_sample~new_var,function(x){paste0(x,collapse = ",")},value.var = "AF_gnomADGenome") %>% rename(AF_gnomADGenome=database)
data_mutation_2=dcast(data_mutation,snp_i+COH_sample~new_var,function(x){paste0(x,collapse = ",")},value.var = "AF_gnomADExome") %>% rename(AF_gnomADExome=database)
data_mutation_3=dcast(data_mutation,snp_i+COH_sample~new_var,function(x){paste0(x,collapse = ",")},value.var = "AF_dbsnpAll") %>% rename(AF_dbsnpAll=database)

data_mutation2=data_mutation1 %>% left_join(data_mutation_1) %>% left_join(data_mutation_2) %>% left_join(data_mutation_3)

bam_dir="/scratch/jibzhang/project/DHL_WES/"
ref_genome="/home/zgu_labs/refData/genome/human/GRCh38/bwa/GRCh38.fa"
head(data_mutation)
data_mutation3=data_mutation2 %>% mutate(samtools=paste0("samtools tview -p ",chr,":",start," ",bam_dir,COH_sample,"/EXOME/bam/",COH_sample,".bam ",ref_genome),
                                       snp=paste0(COH_sample,":",chr,":",start)
                                       )

length(unique(data_mutation3$snp_i))

write_tsv(data_mutation3,"id/data_mutation_exome_noNormal.tsv")

#mutect2 mutations with matched normal

file_list=system("ls /scratch/jibzhang/project/DHL_WES/*/EXOME/search_dbsnp_human_Mutect2withNormal/*.search_dbsnp.tsv",intern = T)
file_one=file_list[1]
data_mutation=bind_rows(lapply(file_list,function(file_one){
  data_one=read.table(file_one,sep="\t",stringsAsFactors = F,header = T)
  data_one$AF_gnomADGenome=as.character(data_one$AF_gnomADGenome)
  data_one$AF_gnomADExome=as.character(data_one$AF_gnomADExome)
  data_one$AF_dbsnpAll=as.character(data_one$AF_dbsnpAll)
  data_one$COH_sample=stringr::word(basename(file_one),1,sep="[.]")
  data_one
}))

data_mutation$new_var="database"
data_mutation1=data_mutation %>% select(names(data_mutation)[c(1:8,11:27,31)]) %>% distinct()

data_mutation_1=dcast(data_mutation,snp_i+COH_sample~new_var,function(x){paste0(x,collapse = ",")},value.var = "AF_gnomADGenome") %>% rename(AF_gnomADGenome=database)
data_mutation_2=dcast(data_mutation,snp_i+COH_sample~new_var,function(x){paste0(x,collapse = ",")},value.var = "AF_gnomADExome") %>% rename(AF_gnomADExome=database)
data_mutation_3=dcast(data_mutation,snp_i+COH_sample~new_var,function(x){paste0(x,collapse = ",")},value.var = "AF_dbsnpAll") %>% rename(AF_dbsnpAll=database)

data_mutation2=data_mutation1 %>% left_join(data_mutation_1) %>% left_join(data_mutation_2) %>% left_join(data_mutation_3)

bam_dir="/scratch/jibzhang/project/DHL_WES/"
ref_genome="/home/zgu_labs/refData/genome/human/GRCh38/bwa/GRCh38.fa"
head(data_mutation)
data_mutation3=data_mutation2 %>% mutate(samtools=paste0("samtools tview -p ",chr,":",start," ",bam_dir,COH_sample,"/EXOME/bam/",COH_sample,".bam ",ref_genome),
                                       snp=paste0(COH_sample,":",chr,":",start)
                                       )

length(unique(data_mutation3$snp_i))

write_tsv(data_mutation3,"id/data_mutation_exome_withNormal.tsv")




