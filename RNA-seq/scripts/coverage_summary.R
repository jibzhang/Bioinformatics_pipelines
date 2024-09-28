#gc(rm(list=ls()))
library(dplyr)


file_coverage_list = "/scratch/jibzhang/project/DHL_RNAseq/sample_infor/coverage.list"

coverage_list=readLines(file_coverage_list)


file_one=coverage_list[1]


data_coverage=bind_rows(lapply(coverage_list, function(file_one){
  data_one=read.csv(file_one,stringsAsFactors = F)
  data_one
})) %>% distinct()

#names(data_coverage)=c(names(data_coverage)[1:3],paste0("gte_",1:30,"X"))

#openxlsx::write.xlsx(data_coverage,"/scratch/zuhu/project/BinghuiShen/ProstateCancer/sample_id/data_coverage_WES.xlsx")
write.csv(data_coverage,"/scratch/jibzhang/project/DHL_RNAseq/sample_infor/data_coverage_WES.csv",row.names=F)















