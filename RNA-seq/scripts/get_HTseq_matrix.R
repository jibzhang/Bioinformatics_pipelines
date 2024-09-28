######################################## set ##################################
setwd("/scratch/jibzhang/project/Rat_heart_RNAseq")

library(dplyr)
library(stringr)
library(reshape2)

gc(rm(list=ls()))

#Get HTSeq count matrix --------------------
get_htseq_matrix=function(files_htseq){
  for (file_htseq in files_htseq){
    #df=read.table(file_htseq,header = F,sep = "\t")
    value=readLines(file_htseq)
    
    feature=word(value,1,sep="\t")[substr(value,1,3)=="ENS"]
    count=word(value,2,sep="\t")[substr(value,1,3)=="ENS"]
    df=data.frame(
      feature=feature,
      count=count
    )
    
    names(df)=c("feature",word(basename(file_htseq),1,sep="[.]"))
    print(match(file_htseq,files_htseq))
    if (match(file_htseq,files_htseq)==1){count_all=df}
    if (match(file_htseq,files_htseq)>1){count_all=count_all %>% left_join(df)}
  }
  count_all=count_all %>% arrange(feature)
  count_all
}
#htseq--------------------
files_htseq=readLines("id/HTSeq.list")

count_htseq=get_htseq_matrix(files_htseq)

write.table(count_htseq,"counts/HTSeq_count_matrix.tsv",quote = F,col.names = T,row.names = F,sep = "\t")

