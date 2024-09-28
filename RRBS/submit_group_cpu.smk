#Each Group/Core is limited to 216 cores and 4.6TB in defq, and 192 Cores, 2.2TB and 24 V100 GPU in gpu
include: "common_rule.smk"

group1=['s4_TET2KO','s4_DKO','s4_DKO','s8_TET2KO','s8_DKO','s8_DKO']
group2=['s4_WT','s4_WT','s4_TET2KO','s8_WT','s8_WT','s8_TET2KO']
contrast = [f"{a}vs{b}" for a, b in zip(group1, group2)]
contra_dic = dict((z[0],list(z[1:])) for z in zip (contrast,group1,group2))

rule all:
    input:
        expand(dir_out + "/smk_code/{contrast}.sh", contrast=contra_dic.keys())

rule write_code:
    input:
	    key_smk="group_contrast.smk"
    output:
        code=dir_out+"/smk_code/{contrast}.sh"
    params:
        cores=4,
        mem="60G",
        log=dir_out+"/smk_log/{contrast}_DMRs_%j.log",
        key_file="/coh_labs/jochan/CD4_TET2KO_RRBS/group_contr/log/{contrast}.methylkit.bmk"
    shell:
        '''
        echo "#!/bin/bash" > {output.code}
        echo "#SBATCH --job-name=RRBS_DMRs_{wildcards.contrast}         # Job name" >> {output.code}
        echo "#SBATCH -n {params.cores}                          # Number of cores" >> {output.code}
        echo "#SBATCH -N 1-1                        # Min - Max Nodes" >> {output.code}
        echo "#SBATCH -p compute                        # gpu queue" >> {output.code}
        echo "#SBATCH --mem={params.mem}            # Amount of memory in GB" >> {output.code}
        echo "#SBATCH --time=6:00:00               # Time limit hrs:min:sec" >> {output.code}
        echo "#SBATCH --output={params.log}    # Standard output and error log" >> {output.code}
        echo "snakemake -s {input.key_smk} -p {params.key_file}  -j{params.cores}"  >> {output.code}
        #sbatch {output.code}
        #scontrol update jobid=179841 partition=gpu-scavenge
        '''

















