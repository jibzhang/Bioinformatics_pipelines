library(annotatr)
library(GenomicRanges)
library(ggplot2)
library(plyr)
library(tidyverse)

setwd("/coh_labs/jochan/CD4_TET2KO_RRBS/methylkit/")
gc(rm(list=ls()))
args <- commandArgs(trailingOnly = TRUE)
print(args)

group1=args[1]
group2=args[2]
outdir=paste0(group2,"vs",group1,"/")
prefix=paste0(group2,"vs",group1)

hypertable <- read.table(paste0(outdir,prefix,"_Diffmethyl_25p_hyper.txt"))
hypotable <- read.table(paste0(outdir,prefix,"_Diffmethyl_25p_hypo.txt"))

read_annotations(con = "/coh_labs/jochan/CD4_TET2KO_H3K27ac/consensus_peak/CD4T_btcEnhs.bed", genome = 'hg38', name = 'btcenhancers', format = 'bed')
read_annotations(con = "/coh_labs/jochan/CD4_TET2KO_H3K27ac/consensus_peak/F34_FANTOM5_hg38_enhancers_NOToverlap_with_CD4T_btcEnhs.bed", genome = 'hg38', name = 'otherenhancers', format = 'bed')
# read_annotations(con = "/home/jibzhang/reference/hg38_region/Ensembl_TF_binding_sites.txt", genome = 'hg38', name = 'TFbinding', format = 'bed')

#combine hyper and hypo methylation regions and annotate with annotator
hypertable$Diffmeth = "hyper"
hypotable$Diffmeth = "hypo"
diff <- rbind(hypertable,hypotable)
dm_regions = makeGRangesFromDataFrame(diff, keep.extra.columns = TRUE)
annots = c('hg38_cpgs', 'hg38_basicgenes','hg38_custom_btcenhancers','hg38_custom_otherenhancers')
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
  ylab("Percentage") + xlab("category") + geom_text(aes(label=n),position = position_stack(vjust = 0.5), size = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste0(outdir,prefix," count of diff methylation cpg regions", ".png"), width=3, height=3, dpi=300)

genes = dm_catsum %>% dplyr::filter(grepl("genes",annot.type))
genes = ddply(genes, .(Diffmeth), transform, percent = n/sum(n) * 100)
ggplot(genes, aes(fill=annot.type, y=percent, x=Diffmeth)) + geom_bar(position="stack", stat="identity") +
  ylab("Percentage") + xlab("category") + geom_text(aes(label=n),position = position_stack(vjust = 0.5), size = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste0(outdir,prefix," count of diff methylation gene regions", ".png"), width=4, height=3, dpi=300)

enhancs = dm_catsum %>% dplyr::filter(grepl("enhancers",annot.type))
enhancs = ddply(enhancs, .(Diffmeth), transform, percent = n/sum(n) * 100)
ggplot(enhancs, aes(fill=annot.type, y=percent, x=Diffmeth)) + geom_bar(position="stack", stat="identity") +
  ylab("Percentage") + xlab("category") + geom_text(aes(label=n),position = position_stack(vjust = 0.5), size = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste0(outdir,prefix," count of diff methylation gene enhancers", ".png"), width=4, height=3, dpi=300)
