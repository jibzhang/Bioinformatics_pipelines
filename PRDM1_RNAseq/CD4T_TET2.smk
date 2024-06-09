configfile: 'config_software_human.yaml'
configfile: 'config_sample_CD4RNAseq.yaml'

include: "common_rules.smk"
include: "Generate_bam.smk"
include: "Exon_analysis.smk"
include: "VDJ_RNA.smk"
include: "Fusion.smk"

rule all:
    input:
        feature=dir_out+"/{sample}/"+squence_type+"/log/{sample}_featureCounts.bmk",
        VDJ=dir_out + "/{sample}/" + squence_type + "/log/mixcr/{sample}.mixcr_exportClones.bmk",
        bw=dir_out + "/{sample}/" + squence_type + "/log/{sample}.bamcoverage.bmk",
        TE=dir_out + "/{sample}/" + squence_type + "/log/{sample}.TEcount.bmk",
	Cicero=dir_out + "/{sample}/" + squence_type + "/log/{sample}.Cicero_50M.bmk",
	Fusioncather=dir_out + "/{sample}/" + squence_type + "/log/{sample}.FusionCatcher.bmk"
    output:
        dir_out + "/{sample}/" + squence_type + "/log/done_analysis.txt"
    shell:
        '''
        touch {output}
        '''
