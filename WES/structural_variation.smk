rule SvABA_nopair:
    input:
        bam = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam",
        bam_G1 = dir_out + "/Tonsil_23/EXOME/bam/Tonsil_23.bam",
        bam_G2 = dir_out + "/Tonsil_24/EXOME/bam/Tonsil_24.bam",
        bam_G3 = dir_out + "/PBMC_12/EXOME/bam/PBMC_12.bam",
        bam_G4 = dir_out + "/PBMC_13/EXOME/bam/PBMC_13.bam",
        bam_G5 = dir_out + "/DLBCL_Song_N10/EXOME/bam/DLBCL_Song_N10.bam",
        bam_G6 = dir_out + "/DLBCL_Song_N11/EXOME/bam/DLBCL_Song_N11.bam",
        bam_G7 = dir_out + "/DLBCL_Song_N12/EXOME/bam/DLBCL_Song_N12.bam",
        bam_G8 = dir_out + "/DLBCL_Song_N1/EXOME/bam/DLBCL_Song_N1.bam",
        bam_G9 = dir_out + "/DLBCL_Song_N6/EXOME/bam/DLBCL_Song_N6.bam"
    output:
        vcf = dir_out + "/{sample}/" + squence_type + "/SvABA_nopair/{sample}.SvABA_nopair.bps.txt.gz"
    log:        dir_out + "/{sample}/" + squence_type + "/log/{sample}.SvABA_nopair.log"
    benchmark:  dir_out + "/{sample}/" + squence_type + "/log/{sample}.SvABA_nopair.bmk"
    params:
        ref=config["ref_bwa"],
        out_prefix = dir_out + "/{sample}/" + squence_type + "/SvABA_nopair/{sample}.SvABA_nopair",
        bed=config["svaba_region"],
        blacklist=config["blacklist"]
    threads: 4
    shell:
        '''
        svaba run -t {input.bam} -n {input.bam_G1} {input.bam_G2} {input.bam_G3} {input.bam_G4} {input.bam_G5} {input.bam_G6} {input.bam_G7} {input.bam_G8} {input.bam_G9} --threads {threads} -L 6 --id-string {params.out_prefix} -k {params.bed} -B {params.blacklist} --reference-genome  {params.ref}  &> {log}
        '''

rule SvABA_vep:
    input:
        svaba = dir_out + "/{sample}/" + squence_type + "/log/{sample}.SvABA_nopair.bmk"
    output:
        annot = dir_out + "/{sample}/" + squence_type + "/SvABA_nopair/{sample}.SvABA_nopair.sv_vep.vcf"
    benchmark:  dir_out + "/{sample}/" + squence_type + "/log/{sample}.SvABA_vep.bmk"
    params:
        cache=config["vep_cache"],
        vcf = dir_out + "/{sample}/" + squence_type + "/SvABA_nopair/{sample}.SvABA_nopair.svaba.somatic.sv.vcf"
    shell:
        '''
        PATH="/home/zgu_labs/anaconda3/envs/vep/bin:$PATH"
        vep -cache --species homo_sapiens -a GRCh38 --force --pick -e --vcf --filter_common --fork 1 --format vcf --dir_cache {params.cache} -i {params.vcf} -o {output.annot}
        '''

rule Manta:
    input:
        bam = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam",
        bam_G=dir_out + "/Tonsil_23/EXOME/bam/Tonsil_23.bam"
    output:
        vcf=dir_out + "/{sample}/" + squence_type + "/Manta/results/variants/candidateSmallIndels.vcf.gz",
        vcf_sv=dir_out + "/{sample}/" + squence_type + "/Manta/results/variants/somaticSV.vcf.gz",
        vcf_tbi=dir_out + "/{sample}/" + squence_type + "/Manta/results/variants/candidateSmallIndels.vcf.gz.tbi",
        file_runworkflow=dir_out + "/{sample}/" + squence_type + "/Manta/runWorkflow.py"
    log:        dir_out + "/{sample}/" + squence_type + "/log/{sample}.Manta.log"
    benchmark:  dir_out + "/{sample}/" + squence_type + "/log/{sample}.Manta.bmk"
    params:
        ref=config["ref_gatk"],
        bed=config["manta_region"],
        dir_out=dir_out + "/{sample}/" + squence_type + "/Manta"
    threads: 4
    shell:
        '''
        if [ -d {params.dir_out} ] ;then rm -rf {params.dir_out}; fi
        configManta.py \
            --normalBam {input.bam_G} \
            --tumorBam {input.bam} \
            --referenceFasta {params.ref} \
            --exome --callRegions {params.bed} \
            --runDir {params.dir_out} &> {log}
        {output.file_runworkflow}  -m local -j {threads} &>> {log}
        '''

rule Manta_vep:
    input:
        vcf = rules.Manta.output.vcf_sv
    output:
        annot = dir_out + "/{sample}/" + squence_type + "/Manta/results/variants/{sample}_somaticSV_vep.txt"
    benchmark: dir_out + "/{sample}/" + squence_type + "/log/{sample}.Manta_annotation.bmk"
    params:
        gatk = config["gatk"],
        cache=config["vep_cache"],
        vep_vcf = dir_out + "/{sample}/" + squence_type + "/Manta/results/variants/{sample}_somaticSV_vep.vcf",
        vep_table = dir_out + "/{sample}/" + squence_type + "/Manta/results/variants/{sample}_somaticSV_vep.table"
    shell:
        '''
        PATH="/home/zgu_labs/anaconda3/envs/vep/bin:$PATH"
        vep -cache --species homo_sapiens -a GRCh38 --force --pick -e --vcf --filter_common --fork 1 --format vcf --dir_cache {params.cache} -i {input.vcf} -o {params.vep_vcf}
        {params.gatk} VariantsToTable -V {input.vcf} -F CHROM -F POS -F REF -F ALT -F ID -F CSQ -O {params.vep_table}
        cat {params.vep_table}|perl -lane '@ID=split(":",$F[4]); @vep=split("[|]",$F[5]); print join "\\t", $F[0],$F[1],$F[2],$F[3],$ID[0],$vep[1],$vep[2],$vep[3]' > {output.annot}
        '''

rule Strelka:
    input:
        bam = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam",
        #bam_G = dir_out + "/DLBCL_Song_N12/EXOME/bam/DLBCL_Song_N12.bam",
        bam_G = dir_out + "/Tonsil_23/EXOME/bam/Tonsil_23.bam",
        vcf_indel=rules.Manta.output.vcf
    output:
        snv = dir_out + "/{sample}/" + squence_type + "/Strelka/results/variants/filtered.somatic.snvs.txt",
        indel = dir_out + "/{sample}/" + squence_type + "/Strelka/results/variants/filtered.somatic.indel.txt"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.Strelka.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.Strelka.bmk"
    params:
        ref=config["ref_gatk"],
        dir_out=dir_out + "/{sample}/" + squence_type + "/Strelka/",
        file_runworkflow=dir_out + "/{sample}/" + squence_type + "/Strelka/runWorkflow.py",
        vcf_snv=dir_out + "/{sample}/" + squence_type + "/Strelka/results/variants/somatic.snvs.vcf.gz",
        vcf_indel=dir_out + "/{sample}/" + squence_type + "/Strelka/results/variants/somatic.indels.vcf.gz"
    threads: 6
    shell:
        '''
        if [ -d {params.dir_out} ] ;then rm -rf {params.dir_out}; fi
        /home/zgu_labs/bin/software/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
            --normalBam {input.bam_G} \
            --tumorBam {input.bam} \
            --referenceFasta {params.ref} \
            --indelCandidates {input.vcf_indel} --exome \
            --runDir {params.dir_out} &> {log}
        {params.file_runworkflow} -m local -j {threads} &>> {log}
        zcat {params.vcf_snv}|grep PASS|cut -f 1,2,4,5 > {output.snv}
        zcat {params.vcf_indel}|grep PASS|cut -f 1,2,4,5 > {output.indel}
        '''

# rule lumpy_nopair:
#     input:
#         bam=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam",
#         bai=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam.bai",
#         bmk=dir_out + "/{sample}/" + squence_type + "/log/{sample}_leftAlignIndels.bmk",
#     output:
#         vcf=        dir_out + "/{sample}/" + squence_type + "/lumpy_nopair/{sample}.lumpy_nopair.vcf",
#         vcf_svtype= dir_out + "/{sample}/" + squence_type + "/lumpy_nopair/{sample}.lumpy_nopair.svtype.vcf",
#     log:        dir_out + "/{sample}/" + squence_type + "/log/{sample}.lumpy_nopair.log"
#     benchmark:  dir_out + "/{sample}/" + squence_type + "/log/{sample}.lumpy_nopair.bmk"
#     params:
#         ref=config["ref_bwa"],
#         outdir=dir_out + "/{sample}/" + squence_type + "/lumpy_nopair/",
#         discordant_D_unsort=dir_out + "/{sample}/" + squence_type + "/lumpy_nopair/{sample}.discordants.unsorted.bam",
#         splitters_D_unsorted=dir_out + "/{sample}/" + squence_type + "/lumpy_nopair/{sample}.splitters.unsorted.bam",
#         discordant_D=dir_out + "/{sample}/" + squence_type + "/lumpy_nopair/{sample}.discordants.bam",
#         splitters_D=dir_out + "/{sample}/" + squence_type + "/lumpy_nopair/{sample}.splitters.bam",
#     threads: 4
#     shell:
#         '''
#         if [ ! -d {params.outdir} ] ;then mkdir {params.outdir}; fi
#         samtools view -b -F 1294 {input.bam} > {params.discordant_D_unsort}
#
#         samtools view -h {input.bam} |extractSplitReads_BwaMem -i stdin | samtools view -Sb - > {params.splitters_D_unsorted}
#
#         samtools sort -m 10G -@ {threads} -O bam -T {params.outdir} -o {params.discordant_D} {params.discordant_D_unsort}
#
#         samtools sort -m 10G -@ {threads} -O bam -T {params.outdir} -o {params.splitters_D} {params.splitters_D_unsorted}
#
#         lumpyexpress \
#             -B {input.bam} \
#             -S {params.splitters_D} \
#             -D {params.discordant_D} \
#             -o {output.vcf}
#
#         svtyper -B {input.bam} -i {output.vcf} -o {output.vcf_svtype}
#         rm {params.splitters_D};rm {params.discordant_D}
#         rm {params.splitters_D_unsorted};rm {params.discordant_D_unsort};
#         '''

rule done:
    input:
        strelka=rules.Strelka.benchmark,
        svaba=rules.SvABA_vep.benchmark,
        #manta=rules.Manta_vep.benchmark
    output:
        dir_out + "/{sample}/" + squence_type + "/log/done_sv.txt"
    shell:
        '''
        touch {output}
        '''
