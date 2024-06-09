#Each Group/Core is limited to 216 cores and 4.6TB in defq, and 192 Cores, 2.2TB and 24 V100 GPU in gpu
configfile: 'config_sample_PRDM1_ATAC.yaml'
include: "common_rule.smk"

rule all:
    input:
        expand(dir_out+"/smk_code/{group}.sh",group=samplelist)

rule write_code:
    input:
	    key_smk="PRDM1_ATAC.smk"
    output:
        code=dir_out+"/smk_code/{group}.sh"
    params:
        cores=4,
        mem="60G",
        log=dir_out+"/smk_log/{group}_peakcall_%j.log",
        key_file="/coh_labs/jochan/PRDM1_ATAC/{group}/ATACseq/log/{group}.filter_bed.bmk"
    shell:
        '''
        echo "#!/bin/bash" > {output.code}
        echo "#SBATCH --job-name=peakcall_{wildcards.group}         # Job name" >> {output.code}
        echo "#SBATCH -n {params.cores}                          # Number of cores" >> {output.code}
        echo "#SBATCH -N 1-1                        # Min - Max Nodes" >> {output.code}
        echo "#SBATCH -p defq                        # gpu queue" >> {output.code}
        echo "#SBATCH --mem={params.mem}            # Amount of memory in GB" >> {output.code}
        echo "#SBATCH --time=96:00:00               # Time limit hrs:min:sec" >> {output.code}
        echo "#SBATCH --output={params.log}    # Standard output and error log" >> {output.code}
        echo "snakemake -s {input.key_smk} -p {params.key_file}  -j{params.cores}"  >> {output.code}
        #sbatch {output.code}
        #scontrol update jobid=179841 partition=gpu-scavenge
        '''