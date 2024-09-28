######################################## load ##################################
#setwd("/scratch/jibzhang/project/DHL_WES")

library(dplyr)
library(stringr)
library(tidyr)
library(reshape2)
#gc(rm(list=ls()))
options(scipen = 200)
##################### set file names
#file_out_variants_sample="Jejoye/EXOME/funcotator/Jejoye_SNV.tsv"

#bigBedToBed="/home/zgu_labs/bin/software/bigbed/bigBedToBed"
#bcftools="bcftools"

#dbsnp_common="/ref_genomes/dbSNP/human/153_hg38/dbSnp153Common.bb"
#dbsnp_all="/ref_genomes/dbSNP/human/153_hg38/dbSnp153.bb"

#gnomad_genome="/ref_genomes/gnomad/gnomad.genomes.r3.0.sites.vcf.bgz"
#gnomad_exome= "/ref_genomes/gnomad/GRCh38/Exomes/r2.1.1/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz"

#file_out="Jejoye_SNV.tsv"

args <- commandArgs(trailingOnly = TRUE)
print(args)

file_out_variants_sample=args[1]

bigBedToBed=args[2]
bcftools=args[3]

dbsnp_common=args[4]
dbsnp_all=args[5]

gnomad_genome= args[6]
gnomad_exome=args[7]
file_out=args[8]

##################### read data
variants_sample=read.table(file_out_variants_sample,header=F,stringsAsFactors=F,sep="\t")
names(variants_sample)=c("CHROM","POS","REF","ALT","TYPE","gene","mutation_type","ensembl_ID","AAC","POPAF","DP","AD","AF")

variants_table=variants_sample[c(1,2,3,4,6,7,10,11,12,13,9)] %>% distinct() %>%
  mutate(type_new=ifelse(nchar(REF)==1 & nchar(ALT)==1,"SNV",
                         ifelse(nchar(REF)>1 | nchar(ALT)>1,"INDEL",NA)),
         length_ref=nchar(REF),
         length_alt=nchar(ALT),
         start=POS,
         end=POS-1+ifelse(length_ref>=length_alt,length_ref,length_alt),
         length_ref=NULL,length_alt=NULL) %>% distinct()

#################### search databases

#######search_function:
#the input position of the function get_search_dbsnp is 1-based. Not 0-based.
get_max_maf=function(x){
  if(x ==""){x1=NA}
  if(x!=""){x1=max(as.numeric(unlist(strsplit(x,","))))}
  x1
}


#chr='chrY'
#start=11334004
#end=11334011
#dbsnp_in=dbsnp_all
get_search_dbsnp_one=function(chr,start,end,dbsnp_in){
  start=start-1
  command=paste0(bigBedToBed," -chrom=",chr," -start=",start," -end=",end," ",dbsnp_in," stdout")
  search_out=system(command,intern=T)
  data_search_out=data.frame(search_out=search_out,stringsAsFactors = F)
  
  get_max_maf=function(x){
    if(x ==""){x1=NA}
    if(x!=""){x1=max(as.numeric(unlist(strsplit(x,","))))}
    x1
  }
  
  data_search_out=data_search_out %>% 
    mutate(CHROM=stringr::word(search_out,1,sep="\t"),
           POS=as.numeric(stringr::word(search_out,2,sep="\t"))+1,
           REF=stringr::word(search_out,5,sep="\t"),
           ALT=stringr::word(search_out,7,sep="\t"),
           V10=stringr::word(search_out,10,sep="\t"),
    )
  
  data_search_out=tidyr::separate_rows(data_search_out,ALT,sep=",") %>% filter(!ALT=="")
  data_search_out$AF=unlist(lapply(data_search_out$V10, get_max_maf)) 
  data_search_out = data_search_out %>% mutate(start=POS,
                                               length_ref=nchar(REF),
                                               length_alt=nchar(ALT),
                                               end=POS-1+ifelse(length_ref>=length_alt,length_ref,length_alt),
                                               length_ref=NULL,length_alt=NULL,V10=NULL,search_out=NULL)
  if(nrow(data_search_out)==0){data_search_out=data.frame(CHROM=chr,POS=NA,REF="",ALT="",AF=NA,start=NA,end=NA,stringsAsFactors = F)}
  if(nrow(data_search_out)>=1){data_search_out=data_search_out}
  data_search_out
}


#the input position of the function get_search_gnomad is 1-based. Not 0-based.
#gnomad_in=gnomad_genome
get_search_gnomad_one=function(chr,start,end,gnomad_in){
  para_f="'%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t%INFO/AF\\t%INFO/AC\\t%INFO/AN\\n'"
  search_region=paste0(chr,":",start,"-",end)
  
  command=paste0(bcftools," query -f ",para_f," ",gnomad_in," -r ",search_region," 2>/dev/null ")
  search_out=system(command,intern=T)
  
  data_search_out=data.frame(
    CHROM=stringr::word(search_out,1,sep="\t"),
    POS=as.numeric(stringr::word(search_out,2,sep="\t")),
    REF=stringr::word(search_out,3,sep="\t"),
    ALT=stringr::word(search_out,4,sep="\t"),
    QUAL=stringr::word(search_out,5,sep="\t"),
    FILTER=stringr::word(search_out,6,sep="\t"),
    AF=stringr::word(search_out,7,sep="\t"),
    AC=stringr::word(search_out,8,sep="\t"),
    AN=stringr::word(search_out,9,sep="\t"),
    stringsAsFactors = F
  )
  data_search_out=tidyr::separate_rows(data_search_out,ALT,sep=",") %>% filter(!ALT=="")
  
  data_search_out = data_search_out %>% mutate(start=POS,
                                               length_ref=nchar(REF),
                                               length_alt=nchar(ALT),
                                               end=POS-1+ifelse(length_ref>=length_alt,length_ref,length_alt),
                                               length_ref=NULL,length_alt=NULL)
  if(nrow(data_search_out)==0){data_search_out=data.frame(CHROM=chr,POS=NA,REF="",ALT="",QUAL="",FILTER="",AF="",AC="",AN="",start=NA,end=NA,stringsAsFactors = F)}
  if(nrow(data_search_out)>=1){data_search_out=data_search_out}
  data_search_out
}

#search_n=100
#i=5
#start=23198227
#end=23198227
#i=28
#indata=variants_table
get_search_dataset=function(indata){
  bind_rows(lapply(1:nrow(indata), function(i){
    #print(i)
    for_search_one=indata[i,]
    
    chr=for_search_one$CHROM
    start=for_search_one$start
    end=for_search_one$end
    
    out1=get_search_dbsnp_one(chr,start,end,dbsnp_common)
    out2=get_search_dbsnp_one(chr,start,end,dbsnp_all)

    names(out1)[5]="AF_dbsnp_common"
    names(out2)[5]="AF_dbsnp_all"
    
    out3=get_search_gnomad_one(chr,start,end,gnomad_exome)
    out4=get_search_gnomad_one(chr,start,end,gnomad_genome)
    
    names(out3)[5:9]=paste0(names(out3)[5:9],"_gnomad_exome")
    names(out4)[5:9]=paste0(names(out4)[5:9],"_gnomad_genome")
    
    search_out=merge(
      merge(out4, out3, by = c("CHROM","POS","REF","ALT","start","end"), all = TRUE),
      merge(out1, out2, by = c("CHROM","POS","REF","ALT","start","end"), all = TRUE),
      by = c("CHROM","POS","REF","ALT","start","end"), all = TRUE) %>% filter(!is.na(POS)) %>%
      select(CHROM,POS,REF,ALT,start,end,
             QUAL_gnomad_genome,FILTER_gnomad_genome,AF_gnomad_genome,AC_gnomad_genome,AN_gnomad_genome,
             QUAL_gnomad_exome,FILTER_gnomad_exome,AF_gnomad_exome,AC_gnomad_exome,AN_gnomad_exome,
             AF_dbsnp_all,AF_dbsnp_common)

    names(search_out)[2:6]=paste0(names(search_out)[2:6],"_database")
    merge_=merge(for_search_one,search_out,by="CHROM",all=TRUE)
    merge_
  }))
}

search_out=get_search_dataset(variants_table)
search_out_new=search_out %>% distinct() %>%  filter(REF==REF_database & ALT==ALT_database & start==start_database)

##### filtering based on PASS in database
search_out1=search_out %>% filter(!(type_new=="SNV" & (REF!=REF_database | ALT!=ALT_database))) %>% mutate(newvar="vartmp")


pass_gnomad=dcast(search_out1,CHROM+POS+REF+ALT+type_new+start+end~newvar,function(x){paste0(x,collapse = ",")},value.var = "FILTER_gnomad_genome") %>% rename(FILTER_gnomad_genome1=vartmp)
pass_exome=dcast(search_out1,CHROM+POS+REF+ALT+type_new+start+end~newvar,function(x){paste0(x,collapse = ",")},value.var = "FILTER_gnomad_exome") %>% rename(FILTER_gnomad_exome1=vartmp)

AF_gnomad=dcast(search_out1,CHROM+POS+REF+ALT+type_new+start+end~newvar,function(x){paste0(x,collapse = ",")},value.var = "AF_gnomad_genome") %>% rename(AF_gnomad_genome1=vartmp)
AF_exome=dcast(search_out1,CHROM+POS+REF+ALT+type_new+start+end~newvar,function(x){paste0(x,collapse = ",")},value.var = "AF_gnomad_exome") %>% rename(AF_gnomad_exome1=vartmp)
AF_common=dcast(search_out1,CHROM+POS+REF+ALT+type_new+start+end~newvar,function(x){paste0(x,collapse = ",")},value.var = "AF_dbsnp_common") %>% rename(AF_dbsnp_common1=vartmp)
AF_all=dcast(search_out1,CHROM+POS+REF+ALT+type_new+start+end~newvar,function(x){paste0(x,collapse = ",")},value.var = "AF_dbsnp_all") %>% rename(AF_dbsnp_all1=vartmp)

filters_merge=pass_gnomad %>% left_join(pass_exome) %>% left_join(AF_gnomad) %>% left_join(AF_exome) %>% left_join(AF_common) %>% left_join(AF_all)

merged_=variants_table %>% left_join(search_out_new) %>% left_join(filters_merge)

var_af=c("AF_gnomad_genome1","AF_gnomad_exome1","AF_dbsnp_common1","AF_dbsnp_all1")

var_filter=c('FILTER_gnomad_genome1','FILTER_gnomad_exome1')

merged_$rare_variants=apply(merged_[var_af],1,function(x){
  value=as.numeric(unlist(strsplit(paste0(x,collapse = ","),",")))
  value1=value[!is.na(value)]
  ifelse(length(value1)==0,NA,ifelse(all(value1<0.01),"Rare","NonRare"))})

merged_$pass_gnomad=apply(merged_[var_filter],1,function(x){
  value=unlist(strsplit(paste0(x,collapse = ","),","))
  value1=value[!value=="NA"]
  ifelse(length(value1)==0,NA,ifelse(all(value1=="PASS"),"PASS_gnomad","NotPass_gnomad"))
})


write.table(merged_,file_out,col.names = T,row.names=F,quote = F,na="",sep="\t")






