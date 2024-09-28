library(tidyverse)
library(ggplot2)
library(reshape2)
setwd("/scratch/jibzhang/project/DHL_WES/id/")
qc <- read.table("qc_EXOME.tsv", header = TRUE, sep = "\t",check.names=FALSE)
tiff("reads_number_per_sample.tiff", units="in", width=10, height=8, res=400)
ggplot(qc, aes(x=Sample, y=ReadsNumber)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,20000000),labels = scales::comma) +
  geom_bar(stat="identity") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x=element_text(angle = 90, colour = "black",vjust = 0.5),
        axis.text.y=element_text(colour = "black",size = 18))
dev.off()

qc <- qc %>% mutate(mappercent=as.numeric(substr(MappingPercentage,1,5)))
tiff("mapping_percent_per_sample.tiff", units="in", width=10, height=8, res=400)
ggplot(qc, aes(x=Sample, y=mappercent)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,100)) +
  geom_bar(stat="identity") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x=element_text(angle = 90, colour = "black",vjust = 0.5),
        axis.text.y=element_text(colour = "black",size = 18))
dev.off()

tiff("dup_rate_per_sample.tiff", units="in", width=10, height=8, res=400)
ggplot(qc, aes(x=Sample, y=DupRates)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,0.8)) +
  geom_bar(stat="identity") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x=element_text(angle = 90, colour = "black",vjust = 0.5),
        axis.text.y=element_text(colour = "black",size = 18))
dev.off()

tiff("average_depth_per_sample.tiff", units="in", width=10, height=8, res=400)
ggplot(qc, aes(x=Sample, y=Depth)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,15)) +
  geom_bar(stat="identity") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x=element_text(angle = 90, colour = "black",vjust = 0.5),
        axis.text.y=element_text(colour = "black",size = 18))
dev.off()

qc2 <- qc %>% gather(key, value, -c("Sample","Depth","ReadsNumber","MappingPercentage","DupRates","mappercent"))
qc2$value <- as.numeric(gsub("%","",as.character(qc2$value)))
qc2$key <- factor(qc2$key,
               levels=c(">=30X",">=29X",">=28X",">=27X",">=26X",">=25X",">=24X",">=23X",">=22X",">=21X",">=20X",">=19X",">=18X",">=17X",">=16X",">=15X",">=14X",">=13X",">=12X",">=11X",">=10X",">=9X",">=8X",">=7X",">=6X",">=5X",">=4X",">=3X",">=2X",">=1X"))

tiff("coverage_depth_per_sample.tiff", units="in", width=10, height=8, res=400)
ggplot(qc2, aes(x=Sample, y=value, fill=key)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,260)) +
  geom_bar(position="stack", stat="identity") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x=element_text(angle = 90, colour = "black",vjust = 0.5),
        axis.text.y=element_text(colour = "black",size = 18))
dev.off()