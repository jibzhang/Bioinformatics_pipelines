setwd("/coh_labs/jochan/PRDM1_NK/")

library(dplyr)
library(stringr)

gc(rm(list=ls()))
#Featurecounts --------------------
files_featurecounts=readLines("id/featureCounts.list")
file_featurecounts=files_featurecounts[1]
rm(count_all_featurecounts)
for (file_featurecounts in files_featurecounts){
df=read.table(file_featurecounts,header = T,sep = "\t")
df=df[c(1,ncol(df))]

names(df)=c("feature",word(basename(file_featurecounts),1,sep="[.]"))
print(match(file_featurecounts,files_featurecounts))
if (match(file_featurecounts,files_featurecounts)==1){count_all_featurecounts=df}
if (match(file_featurecounts,files_featurecounts)>1){count_all_featurecounts=count_all_featurecounts %>% left_join(df)}
}
count_all_featurecounts=count_all_featurecounts %>% arrange(feature)
write.table(count_all_featurecounts,"id/featureCounts_all.tsv",quote = F,col.names = T,row.names = F,sep = "\t")
