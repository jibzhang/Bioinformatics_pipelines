######################################## load ##################################
#setwd("/scratch/jibzhang/project/SW_FL_DLBCL/")

library(dplyr)
library(reshape2)
gc(rm(list=ls()))
options(scipen = 200)
source("/home/jibzhang/bin/R_script/functions.R")
##################### set file names
#tsv_snv="SUDHL6/EXOME/perl_annotation_Mutect2noNormal/SUDHL6.annotation.nonsilent.snv.tsv"
#tsv_indel="SUDHL6/EXOME/perl_annotation_Mutect2noNormal/SUDHL6.annotation.nonsilent.indel.tsv"

#vcf_dir_genome="/ref_genomes/gnomad/GRCh38/Genomes/v3.1/"

#vcf_gnomADExome="/ref_genomes/gnomad/GRCh38/Exomes/r2.1.1/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz"

#bcftools="/home/jibzhang/miniconda3/envs/bv/bin/bcftools"

#out_snv="xx.snv"
#out_indel="xx.indel"

args <- commandArgs(trailingOnly = TRUE)
print(args)

tsv_snv=args[1]
tsv_indel=args[2]
vcf_dir_genome=args[3]
vcf_gnomADExome=args[4]
bcftools=args[5]

out_snv=args[6]
out_indel=args[7]


# read data -------------
data_snv=read_tsv(tsv_snv) %>% mutate(snp_i=paste0(chr,":",start,":",refBase,varBase))
data_indel=read_tsv(tsv_indel)
#search snv -------------------
#indata=data_snv[data_snv$chr %in% c("chr21","chr22"),] 

filter_snv_gnomad_PassRare_onedatabase=function(indata,bcftools,vcf_file,database_name){
  search_out_snv=search_vcf_dataset(indata,vcf_file,bcftools)
  search_out_snv_af=search_out_snv %>% 
    filter((chr==chr_database & refBase==ref_database & varBase== alt_database & start==start_database)) %>% 
    select(snp_i,AF_database) %>% filter(!is.na(AF_database)) %>% distinct() 
  
  search_out_snv1=search_out_snv%>% 
    filter((chr==chr_database & refBase==ref_database & varBase== alt_database & start==start_database)) %>%
    filter((!filter_database=="PASS" )| as.numeric(AF_database) >0.01)
  search_out_snv2=search_out_snv %>% filter(!snp_i %in% search_out_snv1$snp_i) %>% select(names(indata)) %>% distinct() %>% left_join(search_out_snv_af)
  names(search_out_snv2)[length(names(search_out_snv2))]=paste0("AF_",database_name)
  search_out_snv2
}


filter_snv_gnomad_PassRare_onedatabase_genome=function(indata,bcftools,vcf_dir_genome,database_name){
  chr_list=unique(indata$chr)
  
  bind_rows(lapply(chr_list, function(chr_one){
    data_one=indata[indata$chr==chr_one,]
    print(nrow(data_one))
    vcf_file=paste0(vcf_dir_genome,"gnomad.genomes.v3.1.sites.",chr_one,".vcf.bgz")
    filter_snv_gnomad_PassRare_onedatabase(data_one,bcftools,vcf_file,database_name)
  }))
}

filter_out_snv_gnomADGenome=filter_snv_gnomad_PassRare_onedatabase_genome(data_snv,bcftools,vcf_dir_genome,"gnomADGenome")
filter_out_snv_gnomADExome=filter_snv_gnomad_PassRare_onedatabase(data_snv,bcftools,vcf_gnomADExome,"gnomADExome") %>%
  select(snp_i,AF_gnomADExome)
  
filter_out_snv=filter_out_snv_gnomADGenome %>% inner_join(filter_out_snv_gnomADExome)

write_tsv(filter_out_snv,out_snv)
#search indel -------------------

#indata=data_indel[data_indel$chr %in% c("chr21","chr22"),] 
#vcf_file="/scratch/zuhu/df_tmp/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz"

filter_indel_gnomad_PassRare_onedatabase=function(indata,bcftools,vcf_file,database_name){
  search_out_indel=search_vcf_dataset(indata,vcf_file,bcftools)
  
  search_out_indel_keep=search_out_indel %>% filter((chr==chr_database & filter_database=="PASS" & as.numeric(AF_database) <0.01) | is.na(AF_database)) %>%
    select(c(names(indata),"AF_database"))  %>% mutate(newvar="AF_database")
  
  search_out_indel_keep1=dcast(search_out_indel_keep,snp_i~newvar,function(x){paste0(round(as.numeric(x),4),collapse = ",")},value.var = "AF_database")
  
  search_out_indel_keep2=search_out_indel_keep %>% select(-c(AF_database,newvar)) %>% distinct() %>% left_join(search_out_indel_keep1)
  names(search_out_indel_keep2)[length(names(search_out_indel_keep2))]=paste0("AF_",database_name)
  search_out_indel_keep2
}


filter_indel_gnomad_PassRare_onedatabase_genome=function(indata,bcftools,vcf_dir_genome,database_name){
  chr_list=unique(indata$chr) 
  bind_rows(lapply(chr_list, function(chr_one){
    data_one=indata[indata$chr==chr_one,]
    print(nrow(data_one))
    vcf_file=paste0(vcf_dir_genome,"gnomad.genomes.v3.1.sites.",chr_one,".vcf.bgz")
    filter_indel_gnomad_PassRare_onedatabase(data_one,bcftools,vcf_file,database_name)
  }))
}

if (nrow(data_indel)>1) {
data_indel=data_indel%>% mutate(snp_i=paste0(chr,":",start,":",refBase,varBase))
filter_out_indel_gnomADGenome=filter_indel_gnomad_PassRare_onedatabase_genome(data_indel,bcftools,vcf_dir_genome,"gnomADGenome")
filter_out_indel_gnomADExome=filter_indel_gnomad_PassRare_onedatabase(data_indel,bcftools,vcf_gnomADExome,"gnomADExome") %>%
  select(snp_i,AF_gnomADExome)
filter_out_indel=filter_out_indel_gnomADGenome %>% inner_join(filter_out_indel_gnomADExome)
} else {
filter_out_indel=data_indel
}

write_tsv(filter_out_indel,out_indel)























