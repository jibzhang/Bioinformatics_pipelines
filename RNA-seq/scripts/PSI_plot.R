setwd("C:/Users/jibzhang/PycharmProjects/Mouse_RNAseq_splice")

library(ggbio)
library(DEXSeq)
countData<- read.table("selectgene_ranges.txt",header = T)
RYR2 <- countData %>% filter(grepl("ENSMUSG00000021313",exon_ID))
TITIN <- countData %>% filter(grepl("ENSMUSG00000051747",exon_ID))
CAMK2D <-  countData %>% filter(grepl("ENSMUSG00000053819",exon_ID))
TPM2 <-  countData %>% filter(grepl("ENSMUSG00000028464",exon_ID))

exon.RYR2<-makeGRangesFromDataFrame(RYR2, keep.extra.columns = TRUE)
exon.TITIN<-makeGRangesFromDataFrame(TITIN, keep.extra.columns = TRUE)
exon.CAMK2D<-makeGRangesFromDataFrame(CAMK2D, keep.extra.columns = TRUE)
exon.TPM2<-makeGRangesFromDataFrame(TPM2, keep.extra.columns = TRUE)

plotRangesLinkedToData(exon.RYR2, stat.y = c("RSdel","WT"), stat.ylab="PSI", linetype = 5)
ggsave("RYR2_PSI.pdf",width=15,height = 8, dpi=300)
plotRangesLinkedToData(exon.TITIN, stat.y = c("RSdel","WT"), stat.ylab="PSI", linetype = 5)
ggsave("TITIN_PSI.pdf",width=15,height = 8, dpi=300)
plotRangesLinkedToData(exon.CAMK2D, stat.y = c("RSdel","WT"), stat.ylab="PSI", linetype = 5)
ggsave("CAMK2D_PSI.pdf",width=15,height = 8, dpi=300)
plot<-plotRangesLinkedToData(exon.TPM2, stat.y = c("RSdel","WT"), stat.ylab="PSI", linetype = 5)
ggsave("TPM2_PSI.pdf",width=15,height = 8, dpi=300)

coeff <- 
ylimn <- c(0, max(coeff, na.rm = TRUE))
coeff <- vst(coeff, object)
numexons <- nrow(count)
rt <- which(TPM2$exon_ID == "ENSMUSG00000028464:001")
rcolorlines <- ifelse(each <= FDR, "#F219ED60", "#B3B3B360")
rango <- seq(along = rt)
intervals <- (0:nrow(count))/nrow(count)
drawPlot(matr = coeff, ylimn, object, intervals, rango, 
         textAxis = "Exon usage", rt = rt, color = rep(color[colnames(coeff)], 
                                                       each = numexons), colorlines = colorlines, ...)
