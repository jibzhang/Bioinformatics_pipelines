library(methylKit)
library(genomation)
library(tidyverse)
library(fgsea)
library(annotatr)

setwd ("/scratch/jibzhang/project/CD4_TET2KO_RRBS/methylkit_pairwise/")
#setwd ("D:/City_of_Hope/CD4_TET2_RRBS")
gc(rm(list=ls()))
args <- commandArgs(trailingOnly = TRUE)
print(args)

#file1="COH000337_M40D40WT_pe.bismark.cov.out"
#file2="COH000338_M40D40KO_pe.bismark.cov.out"
#name1="M40D40WT"
#name2="M40D40KO"
#prefix="M40D40_KOvsWT"

file1=args[1]
file2=args[2]
name1=args[3]
name2=args[4]
ourdir=args[5]
prefix=args[6]
hmedip_hypo=args[7]

gene.obj=readTranscriptFeatures(location = "hg38_refseq_gene.bed",up.flank = 2000, down.flank = 2000, remove.unusual = T)
cpg.obj=readFeatureFlank(location = "hg38_cpgi.bed",feature.flank.name=c("CpGi","shores"),flank=2000, remove.unusual = T)

#flist <- as.list(file1,file2)
#sublist = flist[c(3,4)]
myobj=methRead(location=list(file1,file2),
               sample.id=list(name1,name2),
               assembly="hg38",
               treatment=c(0,1),
               context="CpG",
               mincov = 5,
               header=FALSE,
               pipeline = "bismarkCoverage")
#myobj=lapply(myobj,function(x) x[x$chr %in% c(paste0('chr', seq(1,22)),'chrX','chrY'),])
tiles = tileMethylCounts(myobj,win.size=300,step.size=300,cov.bases = 10)

filtered.myobj=filterByCoverage(tiles,lo.count=5,lo.perc=NULL, hi.count=NULL, hi.perc=99.9)
normalized.myobj=normalizeCoverage(filtered.myobj,method="median")
meth=methylKit::unite(normalized.myobj,destrand=FALSE)
my.diffMeth<-calculateDiffMeth(meth,mc.cores = 6)

# get hyper methylated bases
myDiff25p.hyper=getMethylDiff(my.diffMeth,difference=25,qvalue=0.01,type="hyper")
# get hypo methylated bases
myDiff25p.hypo=getMethylDiff(my.diffMeth,difference=25,qvalue=0.01,type="hypo")
# get all differentially methylated bases
myDiff25p=getMethylDiff(my.diffMeth,difference=25,qvalue=0.01)

DMRtable <- data.frame(myDiff25p.hyper) 
DMRtable$target.row=seq(1,dim(DMRtable)[1],by=1)
class(DMRtable$target.row)="integer"
diffAnn=annotateWithGeneParts(as(DMRtable,"GRanges"),gene.obj)
TSSannotaton=getAssociationWithTSS(diffAnn)

diffCpGann=annotateWithFeatureFlank(as(myDiff25p,"GRanges"),cpg.obj$CpGi,cpg.obj$shores,feature.name="CpGi",flank.name="shores")
perc.meth=percMethylation(meth,rowids = T)

bedgraph(myDiff25p,col.name = "meth.diff",file.name = paste0(ourdir,prefix,"_diff_cpg_25p.bed"))
write.table(myDiff25p,paste0(ourdir,prefix,"_Diffmethyl_25p.txt"),sep="\t",quote=F)
write.table(myDiff25p.hyper,paste0(ourdir,prefix,"_Diffmethyl_25p_hyper.txt"),sep="\t",quote=F)
write.table(myDiff25p.hypo,paste0(ourdir,prefix,"_Diffmethyl_25p_hypo.txt"),sep="\t",quote=F)
write.table(TSSannotaton,paste0(ourdir,prefix,"_Diffmethyl_25p_TSS_annotation.txt"),sep="\t",quote=F)

tiff(paste0(ourdir,prefix, "_diff_cpg_distribution_across_chromosomes.tiff"),units="in", width=7, height=6, res=300)
diffMethPerChr(myDiff25p,plot=TRUE,qvalue.cutoff=0.05, meth.cutoff=25)
dev.off()
write.table(perc.meth,paste0(ourdir,prefix,"_methylation_percentage.txt"),sep = "\t", quote=F)

DMRpromoter <- TSSannotaton %>% filter(abs(dist.to.feature)<2000)
DMRdata <- left_join(DMRpromoter,DMRtable,by="target.row")
DMRdata$feature.name<- gsub("\\..*","",DMRdata$feature.name)
#Genename <- read.table("Refseq gene name.txt",sep="\t",fill = T,header = T,na.strings=c("","NA")) %>% mutate(ID==ifelse(RefSeq.mRNA.ID=="NA",RefSeq.mRNA.predicted.ID,RefSeq.mRNA.ID)) %>% select(Gene.name,ID)
#write.table(Genename,"Genename with refseqID.txt", sep = "\t",row.names = FALSE,quote = F)
Genename <- read.table("Genename with refseqID.txt",sep="\t",header = T)
DMRgenename <- merge(DMRdata,Genename,by.x="feature.name",by.y="ID") %>% distinct()
write.table(DMRgenename,paste0(ourdir,prefix,"_promoter_hypermethylation_gene_names.txt"),sep = "\t", quote=F, row.names = F)

pathways.hallmark <- gmtPathways("h.all.v7.5.1.symbols.gmt")
rankgenes <- DMRgenename %>% select(Gene.name,meth.diff) %>% na.omit() %>% group_by(Gene.name) %>% summarize(stat=mean(meth.diff))
ranks <- deframe(rankgenes)
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks) %>% as.data.frame()

fgseaRes <- fgseaRes %>% dplyr::select(-leadingEdge) %>% arrange(padj)
write.table(fgseaRes,paste0(ourdir,prefix,"_Diffmethyl_25p_hypermethyl_promoter_GSEA_enrichment.txt"),sep="\t",quote=F)

tiff(paste0(ourdir,prefix,"_GSEA_enrichment_of_25p_hypermethyl_promoter.tiff"),units="in", width=7, height=6, res=300)
ggplot(fgseaRes, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=pval<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
dev.off()

read_annotations(con = "/home/jibzhang/reference/hg38_region/Enhancer_atlas_CD4T_hg38.bed", genome = 'hg38', name = 'CD4Tenhancer', format = 'bed')
read_annotations(con = "/home/jibzhang/reference/hg38_region/GSE138767_TADS.bed", genome = 'hg38', name = 'CD4TTADS', format = 'bed')
read_annotations(con = "/home/jibzhang/reference/hg38_region/Tfh_AllEnhancers.sort.bed", genome = 'hg38', name = 'CD4TfhH3K27', format = 'bed')
read_annotations(con = "/home/jibzhang/reference/hg38_region/Ensembl_TF_binding_sites.txt", genome = 'hg38', name = 'TFbinding', format = 'bed')

#read_annotations(con = "D:/City_of_Hope/CD4_TET2_hMeDIP/M40D40_KO_vs_WT_c3.0_cond1.bed", genome = 'hg38', name = 'hyper5hmc', format = 'bed')
read_annotations(con = hmedip_hypo, genome = 'hg38', name = 'hypo5hmc', format = 'bed')
annots = c('hg38_cpgs', 'hg38_basicgenes', 'hg38_custom_CD4Tenhancer', 'hg38_custom_CD4TTADS', 'hg38_custom_CD4TfhH3K27','hg38_custom_TFbinding','hg38_custom_hypo5hmc')
annotations = build_annotations(genome = 'hg38', annotations = annots)

diffmeth=data.frame(myDiff25p) %>% mutate(stat=ifelse(meth.diff>0, "hyper", "hypo")) 
dm_regions = makeGRangesFromDataFrame(diffmeth, keep.extra.columns = TRUE)
dm_annotated = annotate_regions(regions = dm_regions,annotations = annotations,ignore.strand = TRUE,quiet = FALSE)
dm_annsum = summarize_annotations(annotated_regions = dm_annotated,quiet = TRUE)
dm_catsum = summarize_categorical(annotated_regions = dm_annotated,by = c('annot.type', 'stat'),quiet = TRUE)
write.table(dm_catsum,paste0(ourdir,prefix,"_peak_summary_across_regions.txt"),sep="\t",quote=F)
df_dm_annotated = data.frame(dm_annotated)
write.table(df_dm_annotated,paste0(ourdir,prefix,"_annotatr_annotation.txt"),sep="\t",quote=F)

annots2<-c('hg38_cpgs','hg38_basicgenes','hg38_custom_CD4Tenhancer','hg38_custom_CD4TTADS','hg38_custom_CD4TfhH3K27','hg38_custom_TFbinding')
annotation2 = build_annotations(genome = 'hg38', annotation = annots2)
#overlap<-df_dm_annotated %>% dplyr::filter(annot.type=="hg38_custom_hypo5hmc") %>% dplyr::select(seqnames:stat)
#overlap<-makeGRangesFromDataFrame(overlap, keep.extra.columns = TRUE)
#overlap_annotated= annotate_regions(regions=overlap,annotations = annotation2,ignore.strand = TRUE,quiet = FALSE)
#dm_catsum2 = summarize_categorical(annotated_regions = overlap_annotated,by = c('annot.type', 'stat'),quiet = TRUE)
#write.table(dm_catsum2,paste0(ourdir,prefix,"_overlap_summary_across_regions.txt"),sep="\t",quote=F)
#df_overlap_annotated = data.frame(overlap_annotated)
#write.table(df_overlap_annotated,paste0(ourdir,prefix,"_overlap_annotatr_annotation.txt"),sep="\t",quote=F)


x_order = c('hyper','hypo')
fill_order = c(
  'hg38_custom_hypo5hmc',
  'hg38_custom_CD4Tenhancer',
  'hg38_custom_CD4TTADS',
  'hg38_custom_CD4TfhH3K27',
  'hg38_genes_promoters',
  'hg38_cpg_islands'
  )
plot_categorical(
  annotated_regions = dm_annotated, x='stat', fill='annot.type',
  x_order = x_order, fill_order = fill_order, position='stack',
  plot_title = 'DM Status by annotation Counts',
  legend_title = 'Annotations',
  x_label = 'DM status',
  y_label = 'Counts')
ggsave(paste0(ourdir,prefix,"_DM_Annotation_stackplot",".png"),width=10,height = 7)

#fill_order2 = c('hg38_custom_CD4TTADS','hg38_custom_CD4TfhH3K27','hg38_genes_promoters','hg38_cpg_islands')
#plot_categorical(annotated_regions = overlap_annotated, x='stat', fill='annot.type',x_order = x_order, fill_order = fill_order2, position='stack',plot_title = 'DM Status by annotation Counts',legend_title = 'Annotations',x_label = 'DM status',y_label = 'Counts')
#ggsave(paste0(ourdir,prefix,"_overlap_Annotation_stackplot",".png"),width=10,height = 7)