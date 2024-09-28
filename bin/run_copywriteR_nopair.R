#setwd("/coh_labs/jochan/DHL_WES/CNV_analysis")

library(CopywriteR)

source("/home/zgu_labs/bin/R/functions.R")
# GRch38 unpaired --------------------
threads=4
dir_helper="/home/zgu_labs/bin/software/copywriteR/mm10/mm10_20kb_chr"
file_bam="/home/zgu_labs/pipeline/WES_WGS/test/out/test1_G1/EXOME/bam/test1_G1.bam"
dir_out="/home/zgu_labs/pipeline/WES_WGS/test/out/test1_G1/EXOME/copywriteR"


# arguments --------------------
args <- commandArgs(trailingOnly = TRUE)
print(args)

threads=as.numeric(args[1])
dir_helper=args[2]
file_bam=args[3]
dir_out=args[4]

#parellel parameters --------------------
bp.param <- SnowParam(workers = threads, type = "SOCK")



# run copywriteR ----------------


sample.control = data.frame(samples=file_bam,controls=file_bam)

if(dir.exists(dir_out)){system(paste0("rm -rf ",dir_out))}

if(! dir.exists(dir_out)){dir.create(dir_out)}

CopywriteR(sample.control = sample.control,
           destination.folder = file.path(dir_out),
           reference.folder = file.path(dir_helper),
           bp.param = bp.param)
plotCNA(destination.folder = file.path(dir_out))

load(paste0(dir_out,"/CNAprofiles/segment.Rdata"))

output=segment.CNA.object$output
write_tsv(output,paste0(dir_out,"/results_output.tsv"))


















