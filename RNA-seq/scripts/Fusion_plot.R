library(tidyverse)
library(ggplot2)
setwd("/coh_labs/jochan/HIV_DLBCL_RNAseq/fusion")
lis <- list.files(pattern="*_candidate-fusion-genes.txt")
positive <- lis[c(7:10,12,17,18,20,22,24:26,30,33,34,36,39,40)]
negative <- lis[c(1:6,11,13:16,19,21,23,27:29,31,32,35,37,38)]

plot_fusion <- function(lis, name) {
  datalist <- list()
  for (i in 1:length(lis)) {
    data <- read.table(file=lis[i],header = T,sep = "\t")
    data[, c("Gene_1_symbol.5end_fusion_partner.", "Gene_2_symbol.3end_fusion_partner.")] <- t(apply(data[, c("Gene_1_symbol.5end_fusion_partner.", "Gene_2_symbol.3end_fusion_partner.")], 1, function(x) sort(x)))
    data <- data %>% mutate(fusion=paste0(Gene_1_symbol.5end_fusion_partner.,"_",Gene_2_symbol.3end_fusion_partner.))
    datalist[[i]] = data[!duplicated(data$fusion),]
  }
  combined <- do.call(rbind, datalist)
  head(combined)
  fusion_count <- as.data.frame(table(combined$fusion)) %>% top_n(20,Freq)
  fusion_count
  ggplot(fusion_count, aes(x=reorder(Var1,-Freq),y=Freq)) + geom_col(position='dodge') + xlab("Fusion") + ylab("Frequency") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(paste(name,"top 20 fusion frequency.png"), dpi=300,width=8, height = 5, bg="white")
}

plot_fusion(positive,"positive")
plot_fusion(negative,"negative")