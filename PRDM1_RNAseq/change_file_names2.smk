#dir_out="/scratch/jibzhang/project/CD4T_TET2_methylation/fastq"
dir_out="/coh_labs/jochan/bulkRNA_CD4T"
intake="/coh_labs/jochan/bulkRNA_CD4T/id/intake_CD4T_transcriptome2.tsv"
squence_type="TRANSCRIPTOME"

import pandas as pd

intake=pd.read_table(intake)
samplelist=list(intake.sample_id.drop_duplicates())

rule all:
    input:
        expand(dir_out + "/{sample}/" + squence_type + "/FusionCatcher/{sample}_summary.txt",sample=samplelist)

rule rename:
    input:
        dir_out + "/{sample}/" + squence_type + "/FusionCatcher/summary_candidate_fusions.txt"
    output:
        dir_out + "/{sample}/" + squence_type + "/FusionCatcher/{sample}_summary.txt"
    shell:
        '''
        mv {input} {output}
        '''
