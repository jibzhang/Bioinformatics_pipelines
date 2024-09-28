library(methylKit)
library(genomation)
library(rtacklayer)
setwd ("/scratch/jibzhang/project/CD4_TET2KO_RRBS/methylkit/")
flist <- as.list(list.files(pattern="*.report.txt"))
myobj=methRead(location=flist,
               sample.id=list("M40D90WT","M40D90KO","M40D230KO","F25D90KO_plmd","F25D90WT","F25D90KO_eRNP"),
               assembly="hg38",
               treatment=c(0,1,1,1,0,1),
               context="CpG",
               mincov = 10,
               header=FALSE,
               pipeline = "bismarkCytosineReport"
)
filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)
normalized.myobj=normalizeCoverage(filtered.myobj)
meth=unite(normalized.myobj, destrand=TRUE)
#meth.min=unite(myobj,min.per.group=1L)
getCorrelation(meth,plot=TRUE)
tiff("cluster_dendrogram.tiff",units="in", width=6, height=6, res=300)
clusterSamples(meth, dist="correlation", method="ward.D2", plot=TRUE)
dev.off()
tiff("PCAplot.tiff",units="in", width=7, height=5, res=300)
PCASamples(meth)
dev.off()
tiff("PCA_scree.tiff",units="in", width=5, height=5, res=300)
PCASamples(meth, screeplot=TRUE)
dev.off()

sampleAnnotation=data.frame(sex=c("M","M","M","F","F","F"), time=c(90,90,230,90,90,90))
as=assocComp(mBase=meth,sampleAnnotation)
newObj=removeComp(meth,comp=5)
covariates=data.frame(sex=c(1,1,1,0,0,0), time=c(90,90,230,90,90,90))
my.diffMeth<-calculateDiffMeth(meth,covariates=covariates,overdispersion="MN",test="Chisq",mc.cores=2)

# get hyper methylated bases
myDiff25p.hyper=getMethylDiff(my.diffMeth,difference=25,qvalue=0.01,type="hyper")

# get hypo methylated bases
myDiff25p.hypo=getMethylDiff(my.diffMeth,difference=25,qvalue=0.01,type="hypo")

# get all differentially methylated bases
myDiff25p=getMethylDiff(my.diffMeth,difference=25,qvalue=0.01)

gene.obj=readTranscriptFeatures(location = "hg38_refseq_gene.bed",up.flank = 2000, down.flank = 2000, remove.unusual = T)

diffAnn=annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)

methperchr=diffMethPerChr(myDiff25p,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)

#head(getAssociationWithTSS(diffAnn))
#getTargetAnnotationStats(diffAnn,percentage=TRUE,precedence=TRUE)

tiff("percentage of differential methylation in different regions.tiff",units="in", width=7, height=5, res=300)
plotTargetAnnotation(diffAnn,precedence=TRUE,main="differential methylation annotation")
dev.off()
#getFeatsWithTargetsStats(diffAnn,percentage=TRUE)

cpg.obj=readFeatureFlank(location = "hg38_cpgi.bed",feature.flank.name=c("CpGi","shores"), remove.unusual = T)
diffCpGann=annotateWithFeatureFlank(as(myDiff25p,"GRanges"),cpg.obj$CpGi,cpg.obj$shores,feature.name="CpGi",flank.name="shores")
tiff("percentage of differential methylation in cpg islands.tiff",units="in", width=7, height=5, res=300)
plotTargetAnnotation(diffCpGann,col=c("green","gray","white"),main="differential methylation annotation")
dev.off()

promoters=regionCounts(myobj,gene.obj$promoters,cov.bases = 10)
perc.meth=percMethylation(meth,rowids = T)
TSSannotaton=getAssociationWithTSS(diffAnn)

write.table(methperchr,"Diffmethyl 25p per chromosome.txt",sep="\t",quote=F)
write.table(myDiff25p,"Diffmethyl 25p.txt",sep="\t",quote=F)
write.table(myDiff25p.hyper,"Diffmethyl 25p hyper.txt",sep="\t",quote=F)
write.table(myDiff25p.hypo,"Diffmethyl 25p hypo.txt",sep="\t",quote=F)
write.table(TSSannotaton,"Diffmethyl 25p TSS annotation.txt",sep="\t",quote=F)

### Female sample analysis
female = flist[c(4,5,6)]
myobj=methRead(location=female,
               sample.id=list("F25D90KO_plmd","F25D90WT","F25D90KO_eRNP"),
               assembly="hg38",
               treatment=c(1,0,1),
               context="CpG",
               mincov = 10,
               header=FALSE,
               pipeline = "bismarkCytosineReport"
)
#getMethylationStats(myobj[[1]],plot=TRUE,both.strands=FALSE)
#getCoverageStats(myobj[[1]],plot=TRUE,both.strands=FALSE)
filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)
normalized.myobj=normalizeCoverage(filtered.myobj)
meth=unite(normalized.myobj, destrand=TRUE)
tiff("cluster_dendrogram_female.tiff",units="in", width=5, height=5, res=300)
clusterSamples(meth, dist="correlation", method="ward.D2", plot=TRUE)
dev.off()
tiff("PCAplot_female.tiff",units="in", width=7, height=5, res=300)
PCASamples(meth)
dev.off()
tiff("PCA_scree_female.tiff",units="in", width=5, height=5, res=300)
PCASamples(meth, screeplot=TRUE)
dev.off()
my.diffMeth<-calculateDiffMeth(meth,overdispersion="MN",test="Chisq",mc.cores=2)
# get hyper methylated bases
myDiff25p.hyper=getMethylDiff(my.diffMeth,difference=25,qvalue=0.01,type="hyper")
# get hypo methylated bases
myDiff25p.hypo=getMethylDiff(my.diffMeth,difference=25,qvalue=0.01,type="hypo")
# get all differentially methylated bases
myDiff25p=getMethylDiff(my.diffMeth,difference=25,qvalue=0.01)

diffAnn=annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)
diffCpGann=annotateWithFeatureFlank(as(myDiff25p,"GRanges"),cpg.obj$CpGi,cpg.obj$shores,feature.name="CpGi",flank.name="shores")
perc.meth=percMethylation(meth,rowids = T)
TSSannotaton=getAssociationWithTSS(diffAnn)
tiff("female percentage of differential methylation in different regions.tiff",units="in", width=7, height=5, res=300)
plotTargetAnnotation(diffAnn,precedence=TRUE,main="differential methylation annotation")
dev.off()
tiff("female percentage of differential methylation in cpg islands.tiff",units="in", width=7, height=5, res=300)
plotTargetAnnotation(diffCpGann,col=c("green","gray","white"),main="differential methylation annotation")
dev.off()

write.table(myDiff25p,"female Diffmethyl 25p.txt",sep="\t",quote=F)
write.table(myDiff25p.hyper,"female Diffmethyl 25p hyper.txt",sep="\t",quote=F)
write.table(myDiff25p.hypo,"female Diffmethyl 25p hypo.txt",sep="\t",quote=F)
write.table(TSSannotaton,"female Diffmethyl 25p TSS annotation.txt",sep="\t",quote=F)

#male analysis
male = flist[c(1,2)]
myobj=methRead(location=male,
               sample.id=list("M40D90WT","M40D90KO"),
               assembly="hg38",
               treatment=c(0,1),
               context="CpG",
               mincov = 10,
               header=FALSE,
               pipeline = "bismarkCytosineReport"
)
filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)
normalized.myobj=normalizeCoverage(filtered.myobj)
meth=unite(normalized.myobj, destrand=TRUE)
my.diffMeth<-calculateDiffMeth(meth,overdispersion="MN",test="fast.fisher",mc.cores=2)
# get hyper methylated bases
myDiff25p.hyper=getMethylDiff(my.diffMeth,difference=25,qvalue=0.01,type="hyper")
# get hypo methylated bases
myDiff25p.hypo=getMethylDiff(my.diffMeth,difference=25,qvalue=0.01,type="hypo")
# get all differentially methylated bases
myDiff25p=getMethylDiff(my.diffMeth,difference=25,qvalue=0.01)
methperchr=diffMethPerChr(my.diffMeth,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)
diffAnn=annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)
diffCpGann=annotateWithFeatureFlank(as(myDiff25p,"GRanges"),cpg.obj$CpGi,cpg.obj$shores,feature.name="CpGi",flank.name="shores")
perc.meth=percMethylation(meth,rowids = T)
TSSannotaton=getAssociationWithTSS(diffAnn)

promoters=regionCounts(normalized.myobj,gene.obj$promoters)
meth_promoter=unite(promoters, destrand=TRUE)
promoter.diffMeth<-calculateDiffMeth(meth_promoter,overdispersion="MN",test="fast.fisher",mc.cores=2)
# get all differentially methylated bases
promoter_myDiff25p=getMethylDiff(promoter.diffMeth,difference=25,qvalue=0.05)
meth_diffAnn=annotateWithGeneParts(as(promoter_myDiff25p,"GRanges"),gene.obj)
promoter_TSSannotaton=getAssociationWithTSS(meth_diffAnn)
write.table(promoter_TSSannotaton,"male D90 Diffmethyl 25p promoter TSS annotation.txt",sep="\t",quote=F)

write.table(myDiff25p,"male D90 Diffmethyl 25p.txt",sep="\t",quote=F)
write.table(myDiff25p.hyper,"male D90 Diffmethyl 25p hyper.txt",sep="\t",quote=F)
write.table(myDiff25p.hypo,"male D90 Diffmethyl 25p hypo.txt",sep="\t",quote=F)
write.table(methperchr,"male Diffmethyl 25p per chromosome.txt",sep="\t",quote=F)
write.table(TSSannotaton,"male D90 Diffmethyl 25p TSS annotation.txt",sep="\t",quote=F)