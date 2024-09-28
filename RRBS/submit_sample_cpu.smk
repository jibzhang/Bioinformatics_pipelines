#Each Group/Core is limited to 216 cores and 4.6TB in defq, and 192 Cores, 2.2TB and 24 V100 GPU in gpu
configfile: 'config_sample_RRBS.yaml'

import pandas as pd

#sample list
intake=pd.read_table(config["intake"])
intake["index"] = range(len(intake["COH_ID"]))
#generate BAM file name with tmp id and fq names
intake["tmpid"] = intake.ID.str.split("B",n=2,expand=False).str[1]
sample_tmpid = intake[["COH_sample", "tmpid"]].drop_duplicates()
sample_count = sample_tmpid.COH_sample.value_counts().to_dict()
sample_with_tmp = [k for k, v in sample_count.items() if v >= 2]
intake["BAM_TMP"] = intake.COH_sample + "-TMP" + intake.tmpid
x = lambda i: intake.BAM_TMP[i] if intake.COH_sample[i] in sample_with_tmp else intake.COH_sample[i]
intake["BAM_ID"] = intake["index"].apply(x)
#sample id and bam id
samplelist=list(intake.sample_id.drop_duplicates())
bamTmplist=list(intake.BAM_ID.drop_duplicates())
dir_out=config["dir_out"]

rule all:
    input:
        expand(dir_out+"/smk_code/{sample}.sh",sample=bamTmplist)

rule write_code:
    input:
	    key_smk="GenerateAnalysisReadyBam_RRBS.smk"
    output:
        code=dir_out+"/smk_code/{sample}.sh"
    params:
        cores=4,
        mem="60G",
        log=dir_out+"/smk_log/{sample}_mapping_%j.log",
        key_file="/coh_labs/jochan/CD4_TET2KO_RRBS/{sample}/RRBS/log/done_mapping.txt"
    shell:
        '''
        echo "#!/bin/bash" > {output.code}
        echo "#SBATCH --job-name=RRBS_mapping_{wildcards.sample}         # Job name" >> {output.code}
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

















