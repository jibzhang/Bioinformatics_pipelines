library(methylKit)
library(genomation)
library(tidyverse)
library(fgsea)
library(ggplot2)
library(ggrepel)
library(annotatr)
library(GenomicRanges)

setwd ("/coh_labs/jochan/CD4_TET2KO_RRBS/methylkit/")

gc(rm(list=ls()))
args <- commandArgs(trailingOnly = TRUE)
print(args)
group1=args[1]
group2=args[2]
outdir=paste0(group2,"vs",group1,"/")
prefix=paste0(group2,"vs",group1)

dir.create(prefix)

get_substring <- function(string) {
     parts <- unlist(strsplit(string, "_"))
     substring <- paste(parts[2:4], collapse = "_")
    return(substring)
}

get_group<- function(string) {
  parts <- unlist(strsplit(string, "_"))
  substring <- paste(parts[3:4], collapse = "_")
  return(substring)
}

gene.obj=readTranscriptFeatures(location = "hg38_refseq_gene.bed",up.flank = 2000, down.flank = 2000, remove.unusual = T)
cpg.obj=readFeatureFlank(location = "hg38_cpgi.bed",feature.flank.name=c("CpGi","shores"),flank=2000, remove.unusual = T)

flist <- list.files(pattern="*_pe.bismark.cov")
flist <- grep(paste0(group1,"|",group2),flist,value = T)
names <- sapply(flist,get_substring)
myobj=methRead(location=as.list(flist),
               sample.id=as.list(names),
               assembly="hg38",
               treatment=c(0,1,0,1),
               context="CpG",
               header=FALSE, mincov=5,
               pipeline = "bismarkCoverage")
#myobj2=sapply(myobj,function(x) x[x$chr %in% c(paste0('chr', seq(1,22))),], simplify =T)
tiles = tileMethylCounts(myobj,win.size=1000,step.size=1000)

filtered.myobj=filterByCoverage(tiles,lo.count=5,lo.perc=NULL, hi.count=NULL, hi.perc=99.9)
normalized.myobj=normalizeCoverage(filtered.myobj,method="median")
meth=methylKit::unite(normalized.myobj, destrand=FALSE)
sampleAnnotation=data.frame(names)
groups <- sapply(flist,get_group)
sampleAnnotation$group=groups
sampleAnnotation$donor=c(rep("F27",2),rep("F34",2))
as=assocComp(mBase=meth,sampleAnnotation)
#remove component 1 which correspond to donor
meth <- removeComp(meth,comp=c(1))

# pca <- PCASamples(newobj,obj.return=T)
# components <- data.frame(pca[["x"]])
# components$PC2 <- -components$PC2
# components$PC3 <- -components$PC3
# tot_explained_variance_ratio <- summary(pca)[["importance"]]['Proportion of Variance',][1:3]
# tot_explained_variance_ratio <- 100 * sum(tot_explained_variance_ratio)
# fig <- plot_ly(components, x = ~PC1, y = ~PC2, z = ~PC3,
#                color = ~groups, colors = c('purple','orange',"brown","red","green","blue") ) %>%
#   add_markers()
# 
# fig <- fig %>%
#   layout(
#     title = paste("Total Explained Variance =",tot_explained_variance_ratio),
#     scene = list(bgcolor = "white")
#   ) 
# htmlwidgets::saveWidget(as_widget(fig), "RRBS 3D PCA.html")
# fig %>% config(
#   toImageButtonOptions = list(
#     format = "png",
#     filename = "RRBS 3D PCA",
#     width = 900,
#     height = 700,
#     scale=6
#   ))
#clusterSamples(newobj, dist="correlation", method="ward.D2", plot=TRUE)

my.diffMeth<-calculateDiffMeth(meth,mc.cores = 6)
# get hyper methylated bases
myDiff25p.hyper=getMethylDiff(my.diffMeth,difference=25,qvalue=0.05,type="hyper")
# get hypo methylated bases
myDiff25p.hypo=getMethylDiff(my.diffMeth,difference=25,qvalue=0.05,type="hypo")
# get all differentially methylated bases
myDiff25p=getMethylDiff(my.diffMeth,difference=25,qvalue=0.05)

#make volcano plot with all results
my.diffMeth <- data.frame(my.diffMeth)
my.diffMeth$diffmeth <- "NO"
# if quotlog > 0.5 and pvalue < 0.05, set as "UP" 
my.diffMeth$diffmeth[my.diffMeth$meth.diff > 25 & my.diffMeth$qvalue < 0.05] <- "UP"
# if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"
my.diffMeth$diffmeth[my.diffMeth$meth.diff < -25 & my.diffMeth$qvalue < 0.05]  <- "DOWN"

my.diffMeth <- filter(my.diffMeth, !grepl("random|chrUn|chrM", chr)) %>% rownames_to_column(var="target.row")
my.diffMeth$target.row=seq(1,dim(my.diffMeth)[1],by=1)
class(my.diffMeth$target.row)="integer"
Ann=annotateWithGeneParts(as(my.diffMeth,"GRanges"),gene.obj)
TSSannotaton=getAssociationWithTSS(Ann)
DMRdata <- left_join(my.diffMeth,TSSannotaton,by="target.row")
DMRdata$feature.name<- gsub("\\..*","",DMRdata$feature.name)
Genename <- read.table("Genename with refseqID.txt",sep="\t",header = T)
DMRgenename <- merge(DMRdata,Genename,by.x="feature.name",by.y="ID") %>% distinct()
DMRgenename <- na.omit(DMRgenename)
write.table(DMRgenename,paste0(outdir,prefix,"_whole_diffmethylation_result.txt"),sep="\t", quote=F,row.names = F)

pathways.hallmark <- gmtPathways("h.all.v7.5.1.symbols.gmt")
rankgenes <- DMRgenename %>% select(Gene.name,meth.diff) %>% na.omit() %>% group_by(Gene.name) %>% summarize(stat=mean(meth.diff))
ranks <- deframe(rankgenes)
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks) %>% as.data.frame()
fgseaRes <- fgseaRes %>% dplyr::select(-leadingEdge) %>% arrange(padj)
fgseaRes$pathway <- gsub("HALLMARK_","",fgseaRes$pathway)
write.table(fgseaRes,paste0(outdir,prefix,"GSEA enrichment.txt"),sep="\t",quote=F,row.names = F)
fgseaRes <- fgseaRes %>% filter(pval!="NA")
fgseaRes$pathway <- factor(fgseaRes$pathway, levels = fgseaRes[order(fgseaRes$NES,decreasing = F), "pathway"])
ggplot(fgseaRes,aes(x= NES, y=pathway, colour=-log(pval), size=size)) +
  geom_point(stat = "identity") +
  theme(text = element_text(size=10)) +
  scale_color_gradient(low="blue", high="red") + scale_size() +
  guides(colour = guide_colourbar(order=1), size = guide_legend(order=2)) +
  labs(x="NES", y="Pathways", colour="-logPval", size="Size")
ggsave(paste0(outdir,prefix,"GSEA enrichment.png"),dpi=300,width=7, height=10, bg="white")

ggplot(data=DMRgenename, aes(x=meth.diff, y=-log10(qvalue+10^-310), col=diffmeth)) +
  geom_point(size=0.2) + 
  theme_minimal() +
  theme(legend.position = "none") +
  ggtitle(prefix) +
  geom_text_repel(aes(label = ifelse(abs(meth.diff) > 25 & qvalue < 0.05, as.character(DMRgenename$Gene.name), "")),
                  size=2, max.overlaps = 10,seed =1, point.padding = 0.1,box.padding = 0.1) +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-25, 25), col="red") +
  geom_hline(yintercept=1.3, col="red")
ggsave(paste0(outdir,"volcano plot of diffmeth",prefix,".png"),dpi=300,width=9, height = 9, bg="white")

DMRpromoter <- DMRgenename %>% filter(abs(dist.to.feature)<2000)
ggplot(data=DMRpromoter, aes(x=meth.diff, y=-log10(qvalue+10^-310), col=diffmeth)) +
  geom_point(size=0.2) + 
  theme_minimal() +
  theme(legend.position = "none") +
  ggtitle(prefix) +
  geom_text_repel(aes(label = ifelse(abs(meth.diff) > 25 & qvalue < 0.05, as.character(DMRpromoter$Gene.name), "")),
                  size=2, max.overlaps = 10,seed =1, point.padding = 0.1,box.padding = 0.1) +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-25, 25), col="red") +
  geom_hline(yintercept=1.3, col="red")
ggsave(paste0(outdir,"volcano plot of diffmeth promoter",prefix,".png"),dpi=300,width=9, height = 9, bg="white")

hyperCpGann=annotateWithFeatureFlank(as(myDiff25p.hyper,"GRanges"),cpg.obj$CpGi,cpg.obj$shores,feature.name="CpGi",flank.name="shores")
hypoCpGann=annotateWithFeatureFlank(as(myDiff25p.hypo,"GRanges"),cpg.obj$CpGi,cpg.obj$shores,feature.name="CpGi",flank.name="shores")

bedgraph(myDiff25p,col.name = "meth.diff",file.name = paste0(outdir,prefix,"_diff_cpg_25p.bed"))
write.table(myDiff25p,paste0(outdir,prefix,"_Diffmethyl_25p.txt"),sep="\t",quote=F,row.names=F)
write.table(myDiff25p.hyper,paste0(outdir,prefix,"_Diffmethyl_25p_hyper.txt"),sep="\t",quote=F)
write.table(myDiff25p.hypo,paste0(outdir,prefix,"_Diffmethyl_25p_hypo.txt"),sep="\t",quote=F)

tiff(paste0(outdir,prefix, "percentage of hypermethylation in cpg islands.tiff"),units="in", width=7, height=5, res=300)
plotTargetAnnotation(hyperCpGann,col=c("green","gray","white"),main="differential methylation annotation")
dev.off()
tiff(paste0(outdir,prefix, "percentage of hypomethylation in cpg islands.tiff"),units="in", width=7, height=5, res=300)
plotTargetAnnotation(hypoCpGann,col=c("green","gray","white"),main="differential methylation annotation")
dev.off()

cpgs <- read.table("hg38_cpgi.bed",sep="\t",header=F)
tiff(paste0(outdir,prefix, "_diff_cpg_distribution_across_chromosomes.tiff"),units="in", width=7, height=6, res=300)
diffMethPerChr(myDiff25p,plot=TRUE,qvalue.cutoff=0.05,meth.cutoff=25,exclude=levels(factor(cpgs$V1))[grep("_alt|_random|_fix|chrUn_", levels(factor(cpgs$V1)))],keep.empty.chrom = FALSE)
dev.off()
perc.meth=percMethylation(meth,rowids = T)
write.table(perc.meth,paste0(outdir,prefix,"methylation_percentage.txt"),sep="\t", quote=F)

promoter_analysis <- function(hyperDMRtable, change) {
  hyperDMRtable$target.row=seq(1,dim(hyperDMRtable)[1],by=1)
  class(hyperDMRtable$target.row)="integer"
  hyperAnn=annotateWithGeneParts(as(hyperDMRtable,"GRanges"),gene.obj)
  tiff(paste0(outdir,prefix,"percentage of",change,"in different regions.tiff"),units="in", width=7, height=5, res=300)
  plotTargetAnnotation(hyperAnn,precedence=TRUE,main=paste(change,"annotation"))
  dev.off()
  hyperTSSannotaton=getAssociationWithTSS(hyperAnn)
  write.table(hyperTSSannotaton,paste0(outdir,prefix,"25p TSS annotation.txt"),sep="\t",quote=F)
  hyperDMRpromoter <- hyperTSSannotaton %>% filter(abs(dist.to.feature)<2000)
  hyperDMRdata <- left_join(hyperDMRpromoter,hyperDMRtable,by="target.row")
  hyperDMRdata$feature.name<- gsub("\\..*","",hyperDMRdata$feature.name)
  #Genename <- read.table("Refseq gene name.txt",sep="\t",fill = T,header = T,na.strings=c("","NA")) %>% mutate(ID==ifelse(RefSeq.mRNA.ID=="NA",RefSeq.mRNA.predicted.ID,RefSeq.mRNA.ID)) %>% select(Gene.name,ID)
  #write.table(Genename,"Genename with refseqID.txt", sep = "\t",row.names = FALSE,quote = F)
  hyperDMRgenename <- merge(hyperDMRdata,Genename,by.x="feature.name",by.y="ID") %>% distinct()
  hyperDMRgenename <- na.omit(hyperDMRgenename)
  write.table(hyperDMRgenename,paste0(outdir,prefix,"promoter", change, "gene names.txt"),sep = "\t", quote=F, row.names = F)
}

hypertable <- data.frame(myDiff25p.hyper)
promoter_analysis(hypertable,"hypermethylation")
hypotable <- data.frame(myDiff25p.hypo)
promoter_analysis(hypotable,"hypomethylation")

#combine hyper and hypo methylation regions and annotate with annotator
hypertable$Diffmeth = "hyper"
hypotable$Diffmeth = "hypo"
diff <- rbind(hypertable,hypotable)
dm_regions = makeGRangesFromDataFrame(diff, keep.extra.columns = TRUE)
annots = c('hg38_cpgs', 'hg38_basicgenes')
annotations = build_annotations(genome = 'hg38', annotations = annots)
dm_annotated = annotate_regions(regions = dm_regions,annotations = annotations,ignore.strand = TRUE,quiet = FALSE)
dm_annsum = summarize_annotations(annotated_regions = dm_annotated,quiet = TRUE)
dm_catsum = summarize_categorical(annotated_regions = dm_annotated,by = c('annot.type', 'Diffmeth'),quiet = TRUE)
df_dm_annotated = data.frame(dm_annotated)
write.table(dm_catsum,paste0(outdir,prefix,"_peak_summary_across_regions.txt"),sep="\t",quote=F,row.names = F)
write.table(df_dm_annotated,paste0(outdir,prefix,"_annotatr_annotation.txt"),sep="\t",quote=F,row.names = F)

cpg = dm_catsum %>% dplyr::filter(grepl("cpg",annot.type))
cpg = ddply(cpg, .(Diffmeth), transform, percent = n/sum(n) * 100)
ggplot(cpg, aes(fill=annot.type, y=percent, x=Diffmeth)) + geom_bar(position="stack", stat="identity") +
  ylab("Diffmeth region counts") + xlab("category") + geom_text(aes(label=n),position = position_stack(vjust = 0.5), size = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste0(outdir,prefix," count of diff methylation cpg regions", ".png"), width=3, height=3, dpi=300)

genes = dm_catsum %>% dplyr::filter(grepl("genes",annot.type))
genes = ddply(genes, .(Diffmeth), transform, percent = n/sum(n) * 100)
ggplot(genes, aes(fill=annot.type, y=percent, x=Diffmeth)) + geom_bar(position="stack", stat="identity") +
  ylab("Diffmeth region counts") + xlab("category") + geom_text(aes(label=n),position = position_stack(vjust = 0.5), size = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste0(outdir,prefix," count of diff methylation gene regions", ".png"), width=3, height=3, dpi=300)