library(CopywriteR)

# preCopywriteR(output.folder = file.path("./"),
#              bin.size = 100000,  # dependent on coverage
#              ref.genome = "hg19",
#              prefix = "chr")

bp.param <- SnowParam(workers = 8, type = "SOCK")
samples <- list.files(path = "./BAM", pattern = "*recal.bam$", full.names = TRUE)
ctrlIdx <- c()
for (i in seq(2, 50, by = 2)) {
  ctrlIdx <- c(ctrlIdx, i, i)
}
controls <- samples[ctrlIdx]
sample.control <- data.frame(samples, controls)
CopywriteR(sample.control = sample.control,
           destination.folder = file.path("."),
           reference.folder = file.path(".", "hg19_100kb_chr"),
           bp.param = bp.param)
plotCNA(destination.folder = file.path("."))

## Extract segment data and identify CNA regions

load("CNA.segment.Rdata")
cn <- segment.CNA.object$output
cn <- cn[grep("none", cn$ID),]
cn$ID <- gsub("log2.", "", cn$ID)
cn$ID <- gsub(".recal.bam.vs.none", "", cn$ID)

ids <- unique(cn$ID)
cn$CNA <- "Normal"
cn$loc.start <- round(cn$loc.start)
cn$loc.end <- round(cn$loc.end)

# determine sex based on mean of chrX and chrY
yMean <- c()
xMean <- c()
for (i in ids) {
  yMean[i] <- 0
  xMean[i] <- 0
}
for (i in 1:nrow(cn)) {
  if (cn$chrom[i] == 23) {
    xMean[cn$ID[i]] = xMean[cn$ID[i]] + cn$num.mark[i] * cn$seg.mean[i]
  } else if (cn$chrom[i] == 24) {
    yMean[cn$ID[i]] = yMean[cn$ID[i]] + cn$num.mark[i] * cn$seg.mean[i]  
  }
}
plot(xMean, yMean)
male <- c()
for (i in ids) {
  if (yMean[i] < -510) { # value set based on the previous plot
    male[i] <- FALSE
  } else{
    male[i] <- TRUE
  }
}
table(male)
names(male[which(!male)])

# correct the sex information based on clinical data (if any)
# ...
# ...

cn <- cn[which(!(cn$chrom > 23 & !male[cn$ID])), ]

# Assuming tumor purity >= 20%
# Will not determine CNV for samples with <30% tumor purity
attach(cn)
cn[which(
  num.mark >= 20 &
    ((chrom <= 22 & seg.mean >= log2(2.2 / 2)) |
       (chrom >= 23 & male[ID] & seg.mean >= log2(1.2 / 2)) |
       (chrom == 23 & !male[ID] & seg.mean >= log2(2.2 / 2)))),
  'CNA'] <- 'Gain'
cn[which(
  num.mark >= 20 &
    ((chrom <= 22 & seg.mean <= log2(1.8 / 2)) |
       (chrom >= 23 & male[ID] & seg.mean <= log2(0.8 / 2)) |
       (chrom == 23 & !male[ID] & seg.mean <= log2(1.8 / 2)))),
  'CNA'] <- 'Loss'
# remove the two regions with 100% Loss
cn[which(chrom == 6 & loc.end >= 28050000 & loc.start <= 32850000), 'CNA'] <- 'Normal'
cn[which(chrom == 23 & loc.start <= 2250000), 'CNA'] <- 'Normal'
detach(cn)

cna <- cn[which(cn$CNA != "Normal"), ]

write.csv(cn,"NKTCL.CN.csv",row.names = F)
write.csv(cna,"NKTCL.CNA.csv",row.names = F)


# calculate frequencies of CNA along the genome
gsize <- read.table("genome.fa.fai")
head(gsize)
gsize <- gsize[,c(1:2)]
colnames(gsize) <- c('chr', 'size')
gsize$chr <- gsub("chr", "", gsize$chr)
lossFreq <- c()
gainFreq <- c()
for (i in 1:nrow(gsize)) {
  if (gsize$chr[i] == "M") next
  if (gsize$chr[i] == "X"){
    gsize$chr[i] <- 23
  } else if (gsize$chr[i] == "Y") {
    gsize$chr[i] <- 24
  }
  for (j in seq(5E4, gsize$size[i], 1E5)) {
    key <- paste(gsize$chr[i], j, sep = ":")
    lossFreq[key] <- 0
    gainFreq[key] <- 0
  }
}

for (i in 1:nrow(cna)) {
  for (j in seq(cna$loc.start[i], cna$loc.end[i], 1E5)) {
    key <- paste(cna$chrom[i], j, sep = ":")
    if (cna$CNA[i] == "Loss") {
      lossFreq[key] <- lossFreq[key] + 1
    }else{
      gainFreq[key] <- gainFreq[key] + 1
    }
  }
}

quantile(lossFreq)
quantile(gainFreq)
cnaFreq <- data.frame(lossFreq, gainFreq)
head(cnaFreq)
for (i in 1:nrow(cnaFreq)) {
  chrpos <- unlist(strsplit(rownames(cnaFreq)[i], split = ":"))
  cnaFreq$chr[i] <- chrpos[1]
  cnaFreq$pos[i] <- chrpos[2]
} 
cnaFreq$sampleSize <- 209
table(male)  # 155 T 54 F
cnaFreq[which(cnaFreq$chr == 24), 'sampleSize'] <- 155

# smooth the curve to prevent abnormal values
# window = 1 Mb; step = 100 Kb
minChrRow <- 0
maxChrRow <- 0
for (i in 1:24) {
  maxChrRow <- maxChrRow + nrow(cnaFreq[which(cnaFreq$chr == i),])
  for (j in seq(minChrRow + 1, maxChrRow - 9)) {
    cnaFreq[j, 'lossFreq'] <- median(cnaFreq[j : (j + 9), 'lossFreq'])
    cnaFreq[j, 'gainFreq'] <- median(cnaFreq[j : (j + 9), 'gainFreq'])
  }
  minChrRow <- maxChrRow
}
write.csv(cnaFreq, "NKTCL.cnaFreq.csv", row.names = F)
# plot(table(cnaFreq$lossFreq))
# plot(table(cnaFreq$gainFreq))

cnaFreq <- melt(cnaFreq, id.vars = c('chr', 'pos', 'sampleSize'))
colnames(cnaFreq)[4:5] <- c('Type', 'Frequency')
cnaFreq$Frequency <- cnaFreq$Frequency / cnaFreq$sampleSize

# plot the freq
library(ggplot2)
pdf("NKTCL.CnaFreq.pdf")
for(i in 1:24) {
  d <- cnaFreq[which(cnaFreq$chr == i),]
  p <- ggplot(d, aes(as.numeric(pos), Frequency, group = Type, color = Type))
  p <- p + geom_line(size = 2) + ggtitle(paste0("chr", i)) + ylim(0, 0.4)
  print(p)
}
dev.off()

cnaFreq[which(cnaFreq$Frequency>0.5),]

# save.image("NktclCnv.Rdata")
load("NktclCnv.Rdata")
geneInfo <- read.table("NKTCL.genes.withCNA.exp.bed")
colnames(geneInfo) <- c('chr', 'start', 'end', 'LossFreq', 'GainFreq', 'Gene')
geneInfo$pos <- (geneInfo$start + geneInfo$end - 1) / 2
head(geneInfo)
geneInfo <- geneInfo[, c('Gene', 'chr', 'pos')]
cnaFreq[which(cnaFreq$chr == 23), 'chr'] <- 'X'
cnaFreq[which(cnaFreq$chr == 24), 'chr'] <- 'Y'
cnaFreq$chr <- paste0("chr", cnaFreq$chr)
table(cnaFreq$chr)
cnaFreq <- merge(cnaFreq, geneInfo, by = c('chr', 'pos'), all.x = TRUE)
chromosomes <- unique(cnaFreq$chr)
cnaFreq[which(!is.na(cnaFreq$Gene) & cnaFreq$Frequency < 0.1), 'Gene'] <- NA

pdf("NKTCL.CnaFreq.geneMarked.pdf")
for(i in chromosomes) {
  d <- cnaFreq[which(cnaFreq$chr == i),]
  p <- ggplot(d, aes(as.numeric(pos), Frequency, group = Type, color = Type, label = Gene))
  p <- p + geom_line(size = 2) + ggtitle(i) + ylim(0, 0.4) + 
    geom_text(angle = 90, size = 2, color = "black")
  print(p)
}
dev.off()

