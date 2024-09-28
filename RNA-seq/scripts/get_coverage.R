#gc(rm(list=ls()))

args <- commandArgs(trailingOnly = TRUE)
print(args)

#sample_id="COH0000001_G1"
#file_dup_rate = "/scratch/zuhu/project/BinghuiShen/ProstateCancer/COH0000001_G1/markDup/COH0000001_G1.MarkDup.metrics.txt"
#file_count = "/scratch/zuhu/project/BinghuiShen/ProstateCancer/COH0000001_G1/coverage/COH0000001_G1.sample_cumulative_coverage_counts"
#file_coverage="/scratch/zuhu/project/BinghuiShen/ProstateCancer/COH0000001_G1/coverage/COH0000001_G1.sample_cumulative_coverage_proportions"
#log="/scratch/zuhu/project/BinghuiShen/ProstateCancer/COH0000001_G1/coverage/COH0000001_G1_DepthOfCoverage.log"
#out="/scratch/zuhu/project/BinghuiShen/ProstateCancer/COH0000001_G1/coverage/COH0000001_G1_coverage_dupRate.csv"


sample_id=args[1]
file_dup_rate = args[2]
file_count = args[3]
file_coverage=args[4]
#reads=args[5]
log=args[5]
out=args[6]


data1=readLines(file_dup_rate)

dup_rate=as.numeric(unlist(strsplit(data1[8],"\t"))[9])

data3=read.csv(file_coverage,stringsAsFactors = F)
coverage=unlist(data3[3:32])


data_log=readLines(log)
total_bp=as.numeric(unlist(strsplit(data_log[grep("IntervalArgumentCollection",data_log)]," "))[7])
 
data2=read.csv(file_count,stringsAsFactors = F)

#i=2
#if(i==1){counts=unlist(data2[i+2])*i}
#if(i>1){counts=(unlist(data2[i+2])-unlist(data2[i+1]))*i}


count_2_max=unlist(data2[4:ncol(data2)])
count_1_max_1=unlist(data2[3:(ncol(data2)-1)])
count_value=count_1_max_1-count_2_max

average_coverage=(sum(unlist(lapply(count_value,function(x){x*match(x,count_value)}))))/total_bp


#reads_count=readLines(reads)
#reads_count1=reads_count[grepl("0 read1",reads_count,fixed=F)]
#reads_count2=unlist(strsplit(reads_count1," "))[1]

data_out_one=data.frame(matrix(
  c(sample_id,
    round(dup_rate,6),
    round(average_coverage,3),
    coverage),
  nrow=1),stringsAsFactors = F)

names(data_out_one)=c('sample_id','dup_rate','average_depth',paste0("gte_",1:30,"X"))

write.csv(data_out_one,out,quote = F,row.names = F)




