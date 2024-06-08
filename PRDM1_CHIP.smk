configfile: 'config_sample_PRDM1_ChIPseq.yaml'
configfile: 'config_software_human.yaml'

include: "common_rule.smk"
include: "GenerateAnalysisReadyBam.smk"
#include: "GenerateBam_singleend.smk"
include: "peakcalling.smk"
include: "consensus_analysis.smk"

ALL_mapping = expand(dir_out + "/{sample}/" + squence_type + "/log/done_mapping.txt", sample=samplelist)
ALL_analysis = expand(dir_out + "/{sample}/" + squence_type + "/log/done_analysis.txt", sample=CHIP)
group_analysis = expand(dir_out + "/log/{group}.done_consensus.txt", group=grouplist)

TASKS = []
TASKS.extend(ALL_mapping)
TASKS.extend(ALL_analysis)
TASKS.extend(group_analysis)

rule all:
    input: TASKS