######################################## set ################################
#setwd("/scratch/jibzhang/project/DHL_WES/")

library(dplyr)
library(reshape2)
gc(rm(list=ls()))
options(scipen = 200)

#read cosmic data ######################
cosmic_snv=read.table("/ref_genomes/COSMIC/V94/Cosmic.all.V94.snv",header = F,sep="\t",stringsAsFactors = F)
cosmic_indel=read.table("/ref_genomes/COSMIC/V94/Cosmic.all.V94.indel",header = F,sep="\t",stringsAsFactors = F)

names(cosmic_snv)=c('chr','cosmic_start','cosmic_refBase','cosmic_varBase','cosmic_end','cosmicCount','cosmic_status','gene','cosmic_aaPos')
names(cosmic_indel)=c('chr','cosmic_start','cosmic_refBase','cosmic_varBase','cosmic_end','cosmicCount','cosmic_status','gene','cosmic_aaPos')

#cosmic_snv$gene=tolower(cosmic_snv$gene)
#cosmic_indel$gene=tolower(cosmic_indel$gene)
#cancer_gene=tolower(readLines("/ref_genomes/COSMIC/V94/cancerGenes.txt"))
#cosmic_snv = cosmic_snv %>% mutate(snp=paste0(gene,"_",cosmic_aaPos)

#save(cancer_gene,cosmic_snv,cosmic_indel,file="/coh_labs/jochan/reference/cosmic/cosmic_human.rdata")

args <- commandArgs(trailingOnly = TRUE)
print(args)
tsv_snv=args[1]
tsv_indel=args[2]
sample_name=args[3]
strelka_SNV=args[4]
strelka_INDEL=args[5]
filter_SNV=args[6]
filter_INDEL=args[7]

load("/coh_labs/jochan/reference/cosmic/cosmic_human.rdata")

#read in DLBCL mutated genes and B cell expressed genes
Bcellgene <- scan("/coh_labs/jochan/pipeline/WES/B_cell_genes.txt", what=character(), sep="\n")
normal_mutations <- scan("/coh_labs/jochan/pipeline/WES/hg38_normal_mutations.txt", what=character(), sep="\n")
COSMIC_mutation <- read.table("/coh_labs/jochan/pipeline/WES/Genes_with_mutation_COSMIC_DLBCL.tsv", header = TRUE, sep = "\t")
COSMIC_mutation <- COSMIC_mutation[!grepl("_", COSMIC_mutation$Gene.name), ]
COSMIC_mutation$Gene.name = recode(COSMIC_mutation$Gene.name,'1-Dec'='DELEC1',"1-Mar"="MARCHF1","10-Mar"="MARCHF10",
                                   "10-Sep"="SEPTIN10","11-Mar"="MARCHF11","14-Sep"="SEPTIN14","3-Mar"="MARCHF3",
                                   "3-Sep"="SEPTIN3","4-Mar"="MARCHF4","4-Sep"="SEPTIN4","5-Mar"="MARCHF5","6-Sep"="SEPTIN6",
                                   "7-Mar"="MARCHF7","7-Sep"="SEPTIN7","8-Mar"="MARCHF8","9-Sep"="SEPTIN9")
COSMIC_mutation <- aggregate(. ~ Gene.name, COSMIC_mutation, sum)
COSMIC_mutation <- COSMIC_mutation %>% mutate(frequency=Mutated.samples/Samples.tested) %>% filter(frequency > 0.1 & Samples.tested > 1000)

#read data ----------------

data_snv=read.table(tsv_snv,sep="\t",header=T) %>% filter(MAF>0.03 & MAF<0.05)
data_indel=read.table(tsv_indel,sep="\t",header=T) %>% filter(MAF>0.03 & MAF<0.05)

data_snv$Gene=tolower(data_snv$Gene)
data_indel$Gene=tolower(data_indel$Gene)

data_snv=data_snv %>% filter(tolower(mutType) %in% c("snp","snv")) %>% mutate(anchor=paste0(chr,":",start,"-",end))
data_indel=data_indel %>% filter(tolower(mutType) %in% c("del","ins","indel")) %>% mutate(anchor=paste0(chr,":",start,"-",end))

cosmic_snv <- cosmic_snv %>% mutate(anchor=paste0(chr,":",cosmic_start,"-",cosmic_end))
cosmic_indel <- cosmic_indel %>% mutate(anchor=paste0(chr,":",cosmic_start,"-",cosmic_end))

data_snv1=data_snv %>% left_join(cosmic_snv,by="anchor") %>% filter(!is.na(cosmic_status)) %>% filter(AAC==cosmic_aaPos)
data_indel1=data_indel %>% left_join(cosmic_indel,by="anchor") %>% filter(!(end < cosmic_start | start > cosmic_end))
if (nrow(data_indel1) > 0) {
cosmic_filtered=bind_rows(data_snv1,data_indel1)
} else {
cosmic_filtered=data_snv1
}
cosmic_filtered$chr=cosmic_filtered$chr.x
cosmic_filtered=cosmic_filtered %>% select(-chr.x,-chr.y,-anchor,-gene,-cosmic_start,-cosmic_end)

strlka_snv=read.table(strelka_SNV,sep="\t",header=F) %>% mutate(location=paste0(V1,":",V2,V3,"_",V4))
filter_snv= read.table(filter_SNV,sep="\t",header=F,skip=1) %>% mutate(location=paste0(V1,":",V2,V3,"_",V4))
data_snv2=read.table(tsv_snv,sep="\t",header=T) %>% mutate(location=paste0(chr,":",start,refBase,"_",varBase)) %>% filter(location %in% unique(c(strlka_snv$location,filter_snv$location)) & MAF>=0.05) %>% select(-location)
strlka_indel=read.table(strelka_INDEL,sep="\t",header=F) %>% mutate(location=paste0(V1,":",V2,V3,"_",V4))
#filter_indel= read.table(filter_INDEL,sep="\t",header=F,skip=1) %>% mutate(location=paste0(V1,":",V2,V3,"_",V4))
data_indel2=read.table(tsv_indel,sep="\t",header=T) %>% mutate(location=paste0(chr,":",start,refBase,"_",varBase)) %>% filter(location %in% strlka_indel$location & MAF>=0.05) %>% select(-location)
if (nrow(data_indel2)>0) {
  largeVAF=bind_rows(data_snv2,data_indel2)
} else {
  largeVAF=data_snv2
}
cosmic_filtered$Gene=toupper(cosmic_filtered$Gene)
cancer_gene <- toupper(cancer_gene)
cosmic_filtered = cosmic_filtered %>% filter(Gene %in% Bcellgene)

all = bind_rows(largeVAF,cosmic_filtered)
all = all %>% mutate(Sample=sample_name) %>% select(Sample,Gene,chr,start,end,strand,mutClass,mutType,refBase,varBase,AAC,MAF) %>% filter(Gene %in% Bcellgene) %>% distinct()
all = all %>% mutate(mutations=paste0(chr,"|",start,"|",refBase,"|",varBase)) %>% filter(!mutations %in% normal_mutations) %>% select(-mutations)

colnames(all) = c("Sample","Hugo_Symbol","Chr","Start_position","End_position","Strand","Variant_Classification","Variant_Type","Reference","Tumor_Seq","Protein_Change","VAF")

#select genes expressed in B cells
all_tumor = all %>% filter(Hugo_Symbol %in% COSMIC_mutation$Gene.name)

largeVAF = largeVAF %>% mutate(Sample=sample_name) 

#data_searched=bind_rows(data_snv1,data_indel1) %>% select(snp,cosmicCount,cosmic_status,cosmic_aaPos) %>%
#mutate(var1="cosmicCount",var2="cosmic_status",var3="cosmic_aaPos")

#data_searched1=dcast(data_searched,snp~var1,function(x){paste0(x,collapse = ";")},value.var = "cosmicCount")
#data_searched2=dcast(data_searched,snp~var2,function(x){paste0(x,collapse = ";")},value.var = "cosmic_status")
#data_searched3=dcast(data_searched,snp~var3,function(x){paste0(x,collapse = ";")},value.var = "cosmic_aaPos")

#bam_dir="/scratch/jibzhang/project/DHL_WES/"
#ref_genome="$GRCH38_FA"

#data_variants_out=data_variants[data_variants$gene %in% cancer_gene,]

#data_variants_out=data_variants[data_variants$gene %in% c(data_snv1$gene,data_indel1$gene) | data_variants$gene %in% cancer_gene,] %>%
#group_by(gene) %>% arrange(start) %>% mutate(x=1,gene_count=sum(x),x=NULL) %>% 
#left_join(data_searched1)  %>% left_join(data_searched2)  %>% left_join(data_searched3) %>%
#arrange(desc(gene_count))

write.table(cosmic_filtered,paste0("/coh_labs/jochan/HIV_DLBCL/",sample_name,"/EXOME/perl_filter_Mutect2noNormal/",sample_name,".cosmicfiltered.tsv"),sep="\t",row.names=F,quote=F)
write.table(all,paste0("/coh_labs/jochan/HIV_DLBCL/",sample_name,"/EXOME/perl_filter_Mutect2noNormal/",sample_name,".allmut.tsv"),sep="\t",row.names=F,quote=F)
write.table(largeVAF,paste0("/coh_labs/jochan/HIV_DLBCL/",sample_name,"/EXOME/perl_filter_Mutect2noNormal/",sample_name,".beforerescue.tsv"),sep="\t",row.names=F,quote=F)












