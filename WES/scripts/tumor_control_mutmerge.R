library ("tidyverse")
options(scipen = 200)

args <- commandArgs(trailingOnly = TRUE)
print(args)
tumor_mut=args[1]
parent_mut=args[2]
out_mut=args[3]

parent <- read_delim(parent_mut, delim = "\t")
tumor <- read_delim(tumor_mut, delim = "\t")
parent<-unite(parent, mut_info, gene, start, refAllele, varAllele, sep="_", remove = FALSE)
tumor<-unite(tumor, mut_info,gene, start, refAllele, varAllele, sep="_", remove = FALSE)
sharemut<-merge(tumor, parent, by.x="mut_info", by.y="mut_info")
sharemut<-sharemut %>% mutate(VAFratio=MAF.x/MAF.y) %>% filter (VAFratio<2) 
  #write.table(sharemut, file = paste("U2932 shared mutation list.csv"), sep=",", row.names=F, quote=F)
#uniquemut<-subset(tumor, !(tumor$mut_info %in% sharemut$mut_info))
#write.table(uniquemut, file = "U2932 unique mutation list.tsv", sep="\t",row.names=F,quote=F)
mergemut<-merge(tumor, parent, by.x="mut_info", by.y="mut_info", all.x=TRUE, no.dups = TRUE)
mutation <- mergemut %>% filter(!mut_info %in% sharemut$mut_info) %>%
 select(Chromosome=`#chrom..x`, start=start.x, REF=refAllele.x, ALT=varAllele.x, end=end.x, 
                                Gene=gene.x, Class=mutClass.x, Type=VARIANT_CLASS.x,  Qual=quality.x, Depth=depth.x, REF_num=refDepth.x, 
                                ALT_num=varDepth.x, mapQ=mapQual.x, VAF=MAF.x, ECNT=ECNT.x, MBQ=MBQ.x, MPOS=MPOS.x, 
                                POPAF=POPAF.x, ControlDepth=depth.y, ControlALTnum=varDepth.y, controlVAF=MAF.y)
mutation <- mutation %>% filter(Class %in% c("missense","nonsense","frame_shift_del","frame_shift_ins","in-frame_del","in-frame_ins","splice","startloss","stoploss") & ALT_num >= 3 & VAF>=0.01)
write.table(mutation, out_mut, sep="\t", row.names=F, col.names=F, quote=F)

