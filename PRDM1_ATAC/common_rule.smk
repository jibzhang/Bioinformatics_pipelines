import pandas as pd

#sample list
intake=pd.read_table(config["intake"])
intake["index"] = range(len(intake["COH_ID"]))
intake["rg_sm"] = "SM:" + intake['SM']
intake["rg_pu"] = "PU:" + intake["PU"]

#generate BAM file name with tmp id and fq names
intake["tmpid"] = intake.ID.str.split("B",n=2,expand=False).str[1]
sample_tmpid = intake[["COH_sample", "tmpid"]].drop_duplicates()
sample_count = sample_tmpid.COH_sample.value_counts().to_dict()
sample_with_tmp = [k for k, v in sample_count.items() if v >= 2]
intake["BAM_TMP"] = intake.COH_sample + "-TMP" + intake.tmpid
x = lambda i: intake.BAM_TMP[i] if intake.COH_sample[i] in sample_with_tmp else intake.COH_sample[i]
intake["BAM_ID"] = intake["index"].apply(x)
intake['BAM_file']=intake['BAM_ID']+".bam"

#sample id and bam id
samplelist=list(intake.sample_id.drop_duplicates())
bamTmplist=list(intake.BAM_ID.drop_duplicates())
fqnamelist=list(intake.fq_name.drop_duplicates())
grouplist=list(intake.group.drop_duplicates())

#dict for fastq files and read groups
fq={}
rg_sm={}
rg_pu={}
for i in bamTmplist:
    fq[i]=list(intake[intake["BAM_ID"]==i].fq_name.drop_duplicates())
    rg_sm[i]=list(intake[intake["BAM_ID"]==i].rg_sm.drop_duplicates())

contrast={}
for i in set(grouplist):
    contrast[i]=list(intake[intake["group"]==i].sample_id.drop_duplicates())

#dict for sample to temp bam files
fq_sample = {}
sample_tmp_dict={}
for i in set(samplelist):
    sample_tmp_dict[i]=list(intake[intake["sample_id"]==i].BAM_file.drop_duplicates())
    fq_sample[i] = list(intake[intake["sample_id"] == i].fq_name.drop_duplicates())

dir_out=config["dir_out"]
dir_in=config["dir_in"]
squence_type="ATACseq"
multiqc_name="PRDM1_ATAC"