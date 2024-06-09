setwd("C:/Users/jibzhang/Documents/PRDM1_ATAC")
d6nofeeder<- read.table("D6_nofeeder_consensus.annotation.bed",sep="\t",header=T)
d2nofeeder<- read.table("D2_nofeeder_consensus.annotation.bed",sep="\t",header=T)
d13feeder<- read.table("D13_feeder_consensus.annotation.bed",sep="\t",header=T)
d28feeder<- read.table("D28_feeder_consensus.annotation.bed",sep="\t",header=T)

library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(4, "Pastel2")

# Chart
venn.diagram(
  x = list(d6nofeeder$Gene.Name, d2nofeeder$Gene.Name, d13feeder$Gene.Name,d28feeder$Gene.Name),
  category.names = c("D6_nofeeder" , "D2_nofeeder" , "D13_feeder", "D28_feeder"),
  filename = 'venndiagramm of gene overlap of consensus peak.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 780 , 
  width = 880 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",

)

