import pandas as pd

#sample list
intake=pd.read_table(config["intake"])
intake["rg"]="ID:"+intake['ID']+" SM:"+intake['SM']+" PU:"+intake["PU"]+" PL:"+intake["PL"] + " LB:"+intake["LB"]
intake['BAM_file']=intake['BAM_ID']+".bam"

samplelist=list(intake.sample_id.drop_duplicates())
bamTmplist=list(intake.BAM_ID.drop_duplicates())
#dict for fastq files and read groups
fq={}
rg={}
for i in bamTmplist:
    fq[i]=list(intake[intake["BAM_ID"]==i].fq_name.drop_duplicates())
    rg[i]=list(intake[intake["BAM_ID"]==i].rg.drop_duplicates())

#dict for sample to temp bam files
fq_sample={}
sample_tmp_dict={}
for i in set(samplelist):
    sample_tmp_dict[i]=list(intake[intake["sample_id"]==i].BAM_file.drop_duplicates())
    fq_sample[i]=list(intake[intake["sample_id"]==i].fq_name.drop_duplicates())

dir_out=config["dir_out"]
dir_in=config["dir_in"]
squence_type="TRANSCRIPTOME"
multiqc_name="multiqc_DHL_RNAseq"