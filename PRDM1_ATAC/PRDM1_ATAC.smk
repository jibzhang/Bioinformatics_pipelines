configfile: 'config_sample_PRDM1_ATAC.yaml'
configfile: 'config_software_human.yaml'

include: "common_rule.smk"
include: "GenerateAnalysisReadyBam.smk"
include: "TOBIAS.smk"

individual_mapping= expand(dir_out + "/{sample}/" + squence_type + "/log/{sample}.filter_bed.bmk", sample=samplelist)
group_analysis= expand(dir_out + "/TOBIAS/{group}/log/TOBIAS_score.bmk", group=grouplist)

TASKS = []
TASKS.extend(individual_mapping)
TASKS.extend(group_analysis)

rule all:
    input: TASKS