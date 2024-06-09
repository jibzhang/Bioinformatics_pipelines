configfile: 'config_sample_HIVDLBCL.yaml'

include: "common_rules.smk"

rule all:
    input:
        expand(dir_out+"/smk_code/{sample}.sh",sample=samplelist)

rule write_code:
    input:
	    key_smk="HIV_DLBCL_RNAseq.smk"
    output:
        code=dir_out+"/smk_code/{sample}.sh"
    params:
        cores=4,
        mem="60G",
        log=dir_out+"/smk_log/{sample}__%j.log",
        key_file="/coh_labs/jochan/HIV_DLBCL_RNAseq/{sample}/TRANSCRIPTOME/log/done_analysis.txt"
    shell:
        '''
        echo "#!/bin/bash" > {output.code}
        echo "#SBATCH --job-name=analysis_{wildcards.sample}         # Job name" >> {output.code}
        echo "#SBATCH -n {params.cores}                          # Number of cores" >> {output.code}
        echo "#SBATCH -N 1-1                        # Min - Max Nodes" >> {output.code}
        echo "#SBATCH -p compute                        # gpu queue" >> {output.code}
        echo "#SBATCH --mem={params.mem}            # Amount of memory in GB" >> {output.code}
        echo "#SBATCH --time=36:00:00               # Time limit hrs:min:sec" >> {output.code}
        echo "#SBATCH --output={params.log}    # Standard output and error log" >> {output.code}
        echo "snakemake -s {input.key_smk} -p {params.key_file} -j{params.cores}"  >> {output.code}
        #sbatch {output.code}
        #scontrol update jobid=179841 partition=gpu-scavenge
        '''