library(pheatmap)
library(tidyverse)

setwd("/coh_labs/jochan/CD4_TET2KO_RRBS/methylkit/")
gc(rm(list=ls()))
args <- commandArgs(trailingOnly = TRUE)
print(args)

group1=args[1]
group2=args[2]
outdir=paste0(group2,"vs",group1,"/")
prefix=paste0(group2,"vs",group1)

diffmeth <- read.table(paste0(outdir,prefix,"_whole_diffmethylation_result.txt"),sep="\t",head=T) %>% 
  filter(Gene.name !="" & diffmeth != "NO" & qvalue <.01) %>% 
  slice_max(order_by = abs(meth.diff), n = 100, with_ties = FALSE) %>% 
  mutate(DMR=paste0(chr,".",start,".",end), feature = paste0(Gene.name,"_",dist.to.feature))

percentage <- read.table(paste0(outdir,prefix,"methylation_percentage.txt"),sep="\t",head=T)
overlap <- merge(diffmeth,percentage,by.x="DMR",by.y="row.names")
overlap <- overlap[,c(15:19)]
heatmap_data <- overlap %>% column_to_rownames(var="feature")
heatmap <- pheatmap(t(scale(t(heatmap_data))),fontsize = 6,cluster_rows=TRUE,cluster_cols=TRUE, clustering_distance_rows = "correlation",
                 clustering_distance_cols = "correlation", show_rownames = T, clustering_method = "ward.D2")
save_pheatmap_png <- function(x, filename, width, height) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  png(filename, width=width, height=height, units="in",res=300)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(heatmap,paste0(outdir,prefix,"_heatmap of top 100 DMRs.png"),4,10)
