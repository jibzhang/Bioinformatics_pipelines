#set----------------------------------------------------------------------------------------------------------------------------------
setwd("/scratch/jibzhang/project/DHL_RNAseq/THIT-1/TRANSCRIPTOME/trust4_bam/")
gc(rm(list=ls()))

library(stringr)
library(dplyr)
source("/home/zgu_labs/bin/R/functions.R")

# in_file="COHMM000042_D1_report.tsv"
# top_n=50
# out_prefix="COHMM000042_D1"

args <- commandArgs(trailingOnly = TRUE)
print(args)

in_file=args[1]
top_n=as.numeric(args[2])
out_prefix=args[3]

# in_file="/scratch/zuhu/project/ZhaohuiGu/B-ALL_subtyping/out/COH000304_D1/TRANSCRIPTOME/trust4_fq/COH000304_D1_report.tsv"
# top_n=50
# out_file="/scratch/zuhu/project/ZhaohuiGu/B-ALL_subtyping/out/COH000304_D1/TRANSCRIPTOME/trust4_fq/COH000304_D1_report.topClone.pdf"

#read data ----------------
df_trust=read.table(in_file,skip=1,header = F,stringsAsFactors = F) 

df_n=nrow(df_trust)
n_in=ifelse(df_n>top_n,top_n,df_n)
print(paste0("n_in=",n_in))


df_trust_=df_trust %>%
  mutate(
    v_type=substr(V5,1,3),
    v=word(V5,1,sep="[*]"),
    d=word(V6,1,sep="[*]"),
    j=word(V7,1,sep="[*]"),
    c=substr(V8,1,4),
    in_frame=ifelse(V4=="out_of_frame","",ifelse(V4!=".","*","")),
    vdjc=paste0(v,"_",d,"_",j,"_",c," ",in_frame),
    n=1:n()
         )

df_v=data.frame(
  v_type=c('IGH','IGK','IGL','TRA','TRB','TRD','TRG',"."),
  col=c("seagreen2","gold2","skyblue","#3E9F32","#CCCC33","#1F78B5","#808000",'grey40'),
  stringsAsFactors = F
)

df_trust_top=df_trust_ %>% filter(n<=n_in) %>% left_join(df_v) %>% arrange(V1)

df_trust_top_IGH=df_trust_ %>% filter(v_type=="IGH") %>% mutate(n=1:n()) %>%filter(n<=n_in) %>% left_join(df_v) %>% arrange(V1)

print(paste0("Row N=",nrow(df_trust_top)))

draw_plot=function(indata){
  par(mar=c(5,8,2,2))
  barplot(indata$V1,col=indata$col,horiz=T,names.arg=indata$vdjc,
          las=2,cex.names=.3,
          xlab="Count")
  
  df_v1=df_v %>% filter(v_type %in% indata$v_type)
  legend("bottomright",legend = df_v1$v_type,col=df_v1$col,pch=15,title="V type",bty="n")
}

print(paste0(out_prefix,".trust4Top50.png"))

png(paste0(out_prefix,".trust4Top50.png"),width=1200,height=1600,res=300)
draw_plot(df_trust_top)
dev.off()

print(paste0(out_prefix,".trust4Top50IGH.png"))

png(paste0(out_prefix,".trust4Top50IGH.png"),width=1200,height=1600,res=300)
draw_plot(df_trust_top_IGH)
dev.off()











