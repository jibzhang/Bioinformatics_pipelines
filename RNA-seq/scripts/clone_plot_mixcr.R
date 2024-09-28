library(data.table)
library(ggplot2)
library(dplyr)

setwd("/coh_labs/jochan/bulkRNA_CD4T/clonal_analysis")

count_TCR <- function(pattern,name) {
  list_file <- list.files(pattern = pattern)
  data <- list()
  for (i in 1:length(list_file)){
	results <- read.table(list_file[i], sep = "\t",header = T)
	results <- results[,c(1:3,6:8)]
	colnames(results)<- c("CloneID","Counts","Fraction","Vchain","Dchain","Jchain")
	results$Vchain <- sub("\\*.*", "", results$Vchain)
	results$Dchain <- sub("\\*.*", "", results$Dchain)
	results$Jchain <- sub("\\*.*", "", results$Jchain)
	results <- results %>% mutate(clone = paste0(Vchain,"_",Jchain)) %>% dplyr::select(Counts,Fraction,clone)
	result <-  results %>% arrange(desc(Counts)) %>% dplyr::mutate(ID=row_number()) 
	result$top_clone <- "other"
	result$top_clone[result$ID<=10]=result$ID[result$ID<=10]
	data[[i]] <- result
  	}
   split_strings <- strsplit(list_file,"_")
   names(data) <- sapply(split_strings, function(x) paste(x[1:4], collapse = "_"))
   TCR <- rbindlist(data,idcol=TRUE) %>% dplyr::select(.id,Counts,Fraction,top_clone,clone)
   write.table(TCR,paste(name,"mixcr clones data replicate2.tsv"),sep="\t",row.names=F)
   return(TCR)
  }
  
TCRB <- count_TCR("*_r2.mixcr_TRB_exportClones.tsv","TCRB")
TCRB_sum <- aggregate(cbind(Counts,Fraction)~.id+top_clone,data=TCRB,FUN=sum)
TCRA <- count_TCR("*_r2.mixcr_TRA_exportClones.tsv","TCRA")
TCRA_sum <- aggregate(cbind(Counts,Fraction)~.id+top_clone,data=TCRA,FUN=sum)
TCRB_total <- aggregate(Counts~.id,data=TCRB_sum,FUN=sum)
TCRA_total <- aggregate(Counts~.id,data=TCRA_sum,FUN=sum)
TCR_total <- rbind(TCRA_total,TCRB_total)
TCR_total$TCR <- c(rep("TCRA",nrow(TCRA_total)),rep("TCRB",nrow(TCRB_total)))
TCR_total$ID <- factor(TCR_total$.id,levels = c("F27_S4_3d_WT","F34_S4_3d_WT","F27_S4_3d_TET2KO","F34_S4_3d_TET2KO","F27_S4_3d_DKO","F34_S4_3d_DKO",
						"F27_S4_10d_WT","F34_S4_10d_WT","F27_S4_10d_TET2KO","F34_S4_10d_TET2KO","F27_S4_10d_DKO","F34_S4_10d_DKO",
						"F27_S8_10d_WT","F34_S8_10d_WT","F27_S8_10d_TET2KO","F34_S8_10d_TET2KO","F27_S8_10d_DKO","F34_S8_10d_DKO"))
colfunc<-colorRampPalette(c("royalblue","springgreen","yellow","red"))
ggplot(TCR_total, aes(x=TCR, y=ID)) + 
	geom_point(aes(fill=log10(Counts),size=log10(Counts)), alpha = 0.8, pch=21) +labs(y = "Samples", x = "Total") +theme_bw() +
	scale_fill_gradientn(colours = colfunc(5))+ theme(axis.text.x = element_text(angle = 90)) + scale_size(range=c(1,9))
ggsave("TCR Dotplot of total mixcr clones replicate2.png",dpi=300,width=4, height = 5)

stack_plot <- function(TCR,name) {
	TCR$ID <- factor(TCR$.id,levels = c("F27_S4_3d_WT","F34_S4_3d_WT","F27_S4_3d_TET2KO","F34_S4_3d_TET2KO","F27_S4_3d_DKO","F34_S4_3d_DKO",
		         "F27_S4_10d_WT","F34_S4_10d_WT","F27_S4_10d_TET2KO","F34_S4_10d_TET2KO","F27_S4_10d_DKO","F34_S4_10d_DKO",
                         "F27_S8_10d_WT","F34_S8_10d_WT","F27_S8_10d_TET2KO","F34_S8_10d_TET2KO","F27_S8_10d_DKO","F34_S8_10d_DKO"))
	TCR$top_clone <- factor(TCR$top_clone,levels = c("1","2","3","4","5","6","7", "8","9","10","other"))
	TCR_sum <- TCR %>% dplyr::group_by(ID) %>% dplyr::summarise(total=sum(Counts),percent=sum(Fraction))
	ggplot(TCR) + aes(fill=forcats::fct_rev(top_clone), y=Fraction, x=ID) + geom_bar(position="stack", stat="identity") +
		ylab("Clone Fraction") + xlab("Samples") + geom_text(data=TCR_sum,aes(x=ID, y=percent-0.5, label=total, fill=NULL),angle = 90, hjust = -0.1, colour="white") +
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
		scale_fill_manual(values=c("grey","#8b0000","red","orange","pink","yellow","green","#90ee90","blue",'purple','black'),"Top Clones")
	ggsave(paste(name,"Stacked barplot of top10 mixcr clones replicate2.png"),dpi=300,width=8, height = 6)
	}

stack_plot(TCRA_sum,"TCRA")
stack_plot(TCRB_sum,"TCRB")

dot_plot <- function(TCR,name) {
	TopTCR <- TCR %>% filter(top_clone == "1")
	TCR <- aggregate(cbind(Counts,Fraction)~.id+clone,data=TCR,FUN=sum) %>% filter(clone %in% TopTCR$clone)
    	TCR$ID <- factor(TCR$.id,levels = c("F27_S4_3d_WT","F34_S4_3d_WT","F27_S4_3d_TET2KO","F34_S4_3d_TET2KO","F27_S4_3d_DKO","F34_S4_3d_DKO",
					"F27_S4_10d_WT","F34_S4_10d_WT","F27_S4_10d_TET2KO","F34_S4_10d_TET2KO","F27_S4_10d_DKO","F34_S4_10d_DKO",
					"F27_S8_10d_WT","F34_S8_10d_WT","F27_S8_10d_TET2KO","F34_S8_10d_TET2KO","F27_S8_10d_DKO","F34_S8_10d_DKO"))
        ggplot(TCR, aes(x=clone, y=ID)) +
	  geom_point(aes(fill=log2(Counts),size=Fraction), alpha = 0.8, pch=21) +labs(y = "Samples", x = "TCR") +theme_bw() +
	  scale_fill_gradientn(colours = viridis::plasma(5))+ theme(axis.text.x = element_text(angle = 90)) + scale_size(range=c(1,9))
	ggsave(paste(name,"Dotplot of top mixcr clones replicate2.png"),dpi=300,width=7, height = 6)
	}
dot_plot(TCRA,"TCRA")
dot_plot(TCRB,"TCRB")
