library(pheatmap)
library(tidyverse)
gc(rm(list=ls()))
args <- commandArgs(trailingOnly = TRUE)
print(args)
file=args[1]
output=args[2]

setwd("/coh_labs/jochan/CD4_TET2KO_H3K27ac/DiffBind")
PM_counts=read.table(file, sep="\t", header=TRUE, stringsAsFactor=FALSE)
nrow(PM_counts)
PM_counts=PM_counts[,4:15]
PM_counts = PM_counts[rowSums(PM_counts != 0)> 0,]
nrow(PM_counts)
PM_counts <- PM_counts[!apply(PM_counts, 1, function(row) length(unique(row)) == 1), ]
nrow(PM_counts)
heatmap <- pheatmap(t(scale(t(PM_counts))),fontsize = 6,cluster_rows=TRUE,cluster_cols=TRUE, clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation", show_rownames = F, clustering_method = "ward.D2")
save_pheatmap_png <- function(x, filename, width, height) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  png(filename, width=width, height=height, units="in",res=300)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(heatmap,paste0(output,"_PercMethMatrix_1kbbin_btcenhancer.png"),5,9)
