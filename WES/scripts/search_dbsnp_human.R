######################################## load ##################################
setwd("/scratch/jibzhang/project/DHL_WES/")

library(dplyr)
library(reshape2)
gc(rm(list=ls()))
options(scipen = 200)
source("/home/jibzhang/bin/R_script/functions.R")
##################### set file names
#tsv_snv="SUDHL6/EXOME/search_gnomAD_Mutect2noNormal/SUDHL6.search_gnomAD.snv.tsv"
#tsv_indel="SUDHL6/EXOME/search_gnomAD_Mutect2noNormal/SUDHL6.search_gnomAD.indel.tsv"

#df_dbsnp="/ref_genomes/dbSNP/human/153_hg38/dbSnp153.bb"
#bigBedToBed="/home/zgu_labs/anaconda3/bin/bigBedToBed"

#out_snv="xx.snv"
#out_indel="xx.indel"

args <- commandArgs(trailingOnly = TRUE)
print(args)

tsv_snv=args[1]
tsv_indel=args[2]

df_dbsnp=args[3]
bigBedToBed=args[4]

out_snv=args[5]
out_indel=args[6]


# read data -------------
data_snv=read_tsv(tsv_snv) %>% mutate(snp_i=paste0(chr,":",start,":",refBase,varBase))
data_indel=read_tsv(tsv_indel)


#search snv -------------------
#indata=data_snv[data_snv$chr %in% c("chr21","chr22"),] 


search_bb_snv_Rare_database=function(indata,df_dbsnp,bigBedToBed,database_name){
  search_out_snv=search_bb_dataset(indata,df_dbsnp,bigBedToBed)
  search_out_snv_keep=search_out_snv %>% 
    filter(is.na(AF_database) | (AF_database<0.01 & chr==chr_database & refBase==ref_database & varBase==alt_database)) %>% 
    mutate(newvar="AF_database")
  search_out_snv_keep1=dcast(search_out_snv_keep,snp_i~newvar,function(x){paste0(round(as.numeric(x),4),collapse = ",")},value.var = "AF_database")
  search_out_snv_keep2=search_out_snv_keep %>% select(-c(AF_database,newvar)) %>% distinct() %>% left_join(search_out_snv_keep1)
  search_out_snv_keep2=search_out_snv_keep2 %>% select(c(names(indata),"AF_database"))
  
  names(search_out_snv_keep2)[length(names(search_out_snv_keep2))]=paste0("AF_",database_name)
  search_out_snv_keep2
}


search_bb_indel_Rare_database=function(indata,df_dbsnp,bigBedToBed,database_name){
  search_out_indel=search_bb_dataset(indata,df_dbsnp,bigBedToBed)

  search_out_indel_exclude=search_out_indel %>% filter(AF_database>=0.01 & chr==chr_database) 
  
  search_out_indel_keep=search_out_indel %>% filter(!snp_i %in% search_out_indel_exclude$snp_i) %>% 
    mutate(newvar="AF_database",AF_database=round(as.numeric(AF_database),4))
  
  formula_in=formula(paste0(paste0(names(indata),collapse = "+"),"~","newvar"))
  
  search_out_indel=dcast(search_out_indel_keep,formula_in,function(x){paste0(x,collapse = ",")},value.var = "AF_database")
  names(search_out_indel)[ncol(search_out_indel)]=paste0("AF_",database_name)
  
  search_out_indel
}

filter_out_snv=search_bb_snv_Rare_database(data_snv,df_dbsnp,bigBedToBed,"dbsnpAll")
if (nrow(data_indel)>=1) {
data_indel=data_indel%>% mutate(snp_i=paste0(chr,":",start,":",refBase,varBase))
filter_out_indel=search_bb_indel_Rare_database(data_indel,df_dbsnp,bigBedToBed,"dbsnpAll")
}else{
filter_out_indel=data_indel
}
filter_out_snv=unique(filter_out_snv)
filter_out_indel=unique(filter_out_indel)


write_tsv(filter_out_snv,out_snv)

write_tsv(filter_out_indel,out_indel)























