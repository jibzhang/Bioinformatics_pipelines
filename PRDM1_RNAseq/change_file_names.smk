#dir_out="/scratch/jibzhang/project/CD4T_TET2_methylation/fastq"
dir_out="/coh_labs/jochan/bulkRNA_CD4T/fastq"
intake="/coh_labs/jochan/bulkRNA_CD4T/id/intake_CD4T_transcriptome.tsv"

import pandas as pd

intake=pd.read_table(intake)
newname=list(intake.fq_prefix)

name={}
for i in newname:
    name[i]=intake[intake["fq_prefix"]==i].fq_name

rule all:
    input:
        expand(dir_out+"/{link}.fq.gz",link=newname)

rule rename:
    input:
        lambda wildcards: dir_out + '/' + name[wildcards.link]
    output:
        dir_out + "/{link}.fq.gz"
    shell:
        '''
        ln -s {input} {output}
        '''
