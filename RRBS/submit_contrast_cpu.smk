#Each Group/Core is limited to 216 cores and 4.6TB in defq, and 192 Cores, 2.2TB and 24 V100 GPU in gpu
configfile: 'config_sample_RRBS.yaml'
dir_out=config["dir_out"]
contrast=['M40D40_KO_vs_WT', 'M40D90_KO_vs_WT', 'F25eRNPD40_KO_vs_WT', 'F25eRNPD90_KO_vs_WT', 'M40WT_D90_vs_D40',
'F25WT_D90_vs_D40', 'M40KO_D90_vs_D40', 'F25eRNP_D90_vs_D40', 'M40KO_D230_vs_D40', 'M40KO_D230_vs_D90',
'KOD90vsD40','WTD90vsD40','D40KOvsWT','D90KOvsWT','male','female','female_eRNP','all']

rule all:
    input:
        expand(dir_out+"/smk_code/{sample}.sh",sample=contrast)

rule write_code:
    input:
	    key_smk="Pairwise_Contrast.smk"
    output:
        code=dir_out+"/smk_code/{sample}.sh"
    params:
        cores=4,
        mem="60G",
        log=dir_out+"/smk_log/{sample}_homer_%j.log",
        key_file="/scratch/jibzhang/project/CD4_TET2KO_RRBS/methylkit/log/{sample}.homer.bmk"
    shell:
        '''
        echo "#!/bin/bash" > {output.code}
        echo "#SBATCH --job-name=methyl_homer_{wildcards.sample}        # Job name" >> {output.code}
        echo "#SBATCH -n {params.cores}                        # Number of cores" >> {output.code}
        echo "#SBATCH -N 1-1                        # Min - Max Nodes" >> {output.code}
        echo "#SBATCH -p defq                        # gpu queue" >> {output.code}
        echo "#SBATCH --mem={params.mem}            # Amount of memory in GB" >> {output.code}
        echo "#SBATCH --time=96:00:00               # Time limit hrs:min:sec" >> {output.code}
        echo "#SBATCH --output={params.log}    # Standard output and error log" >> {output.code}
        echo "snakemake -s {input.key_smk} -p {params.key_file}  -j{params.cores}"  >> {output.code}
        #sbatch {output.code}
        #scontrol update jobid=179841 partition=gpu-scavenge
        '''

















