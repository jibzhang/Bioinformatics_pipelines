rule FusionCatcher:
    input:
        fq= lambda wildcards: expand(dir_in + '/{i}',i=fq[wildcards.sample])
    output:
        dir= directory(dir_out + "/{sample}/" + squence_type + "/FusionCatcher")
    log:        dir_out + "/{sample}/" + squence_type + "/log/{sample}.FusionCatcher.log"
    benchmark:  dir_out + "/{sample}/" + squence_type + "/log/{sample}.FusionCatcher.bmk"
    params:
        ref=  config["ref_fusioncatcher"],
        fastq= lambda wildcards,input: ','.join(input)
    threads: 16
    shell:
        '''
        export PATH="/home/zgu_labs/anaconda3/envs/fusioncatcher/bin:$PATH"        
        fusioncatcher.py -p {threads} -d {params.ref} -i {params.fastq} -o {output.dir} &>{log}
        '''

# rule get_fq_50M:
#     input:
#         fq1 = dir_in + "/{sample}.R{i}.fq.gz",
#     output:
#         fq1 =       temp(dir_out + "/{sample}/" + squence_type + "/fusion/get_50M_fq/{sample}.R{i}.fq"),
#     benchmark:      dir_out + "/{sample}/" + squence_type + "/fusion/get_50M_fq/{sample}.R{i}.bmk"
#     threads: 16
#     shell:
#         '''
#         seqtk sample -s100 {input.fq1} 50000000 > {output.fq1}
#         '''
#
# rule Star_50M:
#     input:
#         fq1 = dir_out + "/{sample}/" + squence_type + "/fusion/get_50M_fq/{sample}.R1.fq",
#         fq2 = dir_out + "/{sample}/" + squence_type + "/fusion/get_50M_fq/{sample}.R2.fq"
#     output:
#         bam=    temp(dir_out + "/{sample}/" + squence_type + "/fusion/Star_50M/{sample}.bam"),
#     log:        dir_out + "/{sample}/" + squence_type + "/fusion/Star_50M/{sample}.log",
#     benchmark:  dir_out + "/{sample}/" + squence_type + "/fusion/log/{sample}.bmk",
#     params:
#         star=config["star"],
#         star_ref=config["ref_star"],
#         rg=lambda wildcards: rg[wildcards.sample],
#         dir_out=dir_out + "/{sample}/" + squence_type + "/fusion/Star_50M/",
#         bam_tmp=dir_out + "/{sample}/" + squence_type + "/fusion/Star_50M/Aligned.out.bam",
#     threads: 20
#     shell:
#         '''
#         {params.star} \
#         --twopassMode Basic \
#         --runThreadN {threads} \
#         --limitOutSAMoneReadBytes 90000000 \
#         --outSAMtype BAM Unsorted \
#         --outSAMstrandField intronMotif \
#         --outSAMunmapped Within \
#         --outSAMmapqUnique 60 \
#         --outSAMmultNmax 1 \
#         --outReadsUnmapped None \
#         --outMultimapperOrder Random \
#         --outSAMattributes NH HI AS nM NM MD \
#         --outSAMattrRGline {params.rg} \
#         --outFilterType BySJout \
#         --outFilterMultimapNmax 100 \
#         --outFilterMismatchNmax 999 \
#         --outFilterMismatchNoverReadLmax 0.04 \
#         --alignIntronMin 20 \
#         --alignSJstitchMismatchNmax 5 -1 5 5 \
#         --alignIntronMax 1000000 \
#         --alignMatesGapMax 1000000 \
#         --alignSJoverhangMin 8 \
#         --alignSJDBoverhangMin 4 \
#         --alignInsertionFlush Right \
#         --chimOutType Junctions \
#         --chimSegmentMin 12 \
#         --chimJunctionOverhangMin 12 \
#         --chimSegmentReadGapMax 3 \
#         --chimMultimapNmax 100 \
#         --chimMultimapScoreRange 10 \
#         --chimNonchimScoreDropMin 10 \
#         --chimOutJunctionFormat 1 \
#         --peOverlapNbasesMin 12 \
#         --peOverlapMMp 0.1 \
#         --outFileNamePrefix {params.dir_out} \
#         --genomeDir {params.star_ref} \
#         --genomeLoad NoSharedMemory \
#         --readFilesIn {input.fq1} {input.fq2}\
#         --sjdbOverhang 100 &> {log}
#
#         mv {params.bam_tmp} {output.bam}
#         '''
# #--readFilesCommand zcat \
#
# rule samtools_sort_50M:
#     input:
#         bam_unosrt= dir_out + "/{sample}/" + squence_type + "/fusion/Star_50M/{sample}.bam",
#         bmk=        dir_out + "/{sample}/" + squence_type + "/fusion/log/{sample}.bmk",
#     output:
#         bam=        temp(dir_out + "/{sample}/" + squence_type + "/fusion/samtools_sort_50M/{sample}.bam"),
#         bai=        temp(dir_out + "/{sample}/" + squence_type + "/fusion/samtools_sort_50M/{sample}.bam.bai"),
#     log:            dir_out + "/{sample}/" + squence_type + "/fusion/samtools_sort_50M/{sample}.log",
#     benchmark:      dir_out + "/{sample}/" + squence_type + "/fusion/samtools_sort_50M/{sample}.bmk",
#     params:
#         outdir=     dir_out + "/{sample}/" + squence_type + "/fusion/samtools_sort_50M/"
#     threads: 8
#     shell:
#         '''
#         samtools sort -m 5G -@ {threads} -O bam -T {params.outdir} -o {output.bam} {input.bam_unosrt} &> {log}
#         samtools index {output.bam} {output.bai} &>> {log}
#         '''
#
# rule MarkDup_50M:
#     input:
#         bam=        dir_out + "/{sample}/" + squence_type + "/fusion/samtools_sort_50M/{sample}.bam",
#         bmk=        dir_out + "/{sample}/" + squence_type + "/fusion/samtools_sort_50M/{sample}.bmk",
#     output:
#         bam=        temp(dir_out + "/{sample}/" + squence_type + "/fusion/MarkDup_50M/{sample}.bam"),
#         bai=        dir_out + "/{sample}/" + squence_type + "/fusion/MarkDup_50M/{sample}.bam.bai",
#         metrics=    dir_out + "/{sample}/" + squence_type + "/fusion/MarkDup_50M/{sample}.MarkDup.metrics.txt",
#     log:            dir_out + "/{sample}/" + squence_type + "/fusion/MarkDup_50M/{sample}.log",
#     benchmark:      dir_out + "/{sample}/" + squence_type + "/fusion/MarkDup_50M/{sample}.bmk",
#     params:
#         gatk=config["gatk"]
#     shell:
#         '''
#         {params.gatk} MarkDuplicates -I {input.bam} -O {output.bam} -M {output.metrics} &> {log}
#         samtools index {output.bam} {output.bai} &>> {log}
#         '''
#
# rule Samtools_Flagstat_50M:
#     input:
#         bam=dir_out + "/{sample}/" + squence_type + "/fusion/MarkDup_50M/{sample}.bam",
#         bmk=dir_out + "/{sample}/" + squence_type + "/fusion/MarkDup_50M/{sample}.bmk",
#     output:
#         stats=  dir_out + "/{sample}/" + squence_type + "/fusion/MarkDup_50M/{sample}.mappingStats",
#     benchmark:  dir_out + "/{sample}/" + squence_type + "/fusion/MarkDup_50M/{sample}.samtools_flagstat.bmk"
#     shell:
#         '''
#         samtools flagstat {input.bam} > {output.stats}
#         '''
#
# #module load singularity
# rule RNApeg_50M:
#     input:
#         bam=        dir_out + "/{sample}/" + squence_type + "/fusion/MarkDup_50M/{sample}.bam",
#         bai=        dir_out + "/{sample}/" + squence_type + "/fusion/MarkDup_50M/{sample}.bam.bai",
#         bmk=        dir_out + "/{sample}/" + squence_type + "/fusion/MarkDup_50M/{sample}.samtools_flagstat.bmk"
#     output:
#         RNApeg=     dir_out + "/{sample}/" + squence_type + "/fusion/RNApeg_50M/done",
#         junctions=  dir_out + "/{sample}/" + squence_type + "/fusion/RNApeg_50M/{sample}.bam.junctions.tab.shifted.tab"
#     log:            dir_out + "/{sample}/" + squence_type + "/fusion/RNApeg_50M/{sample}.RNApeg.log"
#     benchmark:      dir_out + "/{sample}/" + squence_type + "/fusion/RNApeg_50M/{sample}.RNApeg.bmk"
#     params:
#         ref1=   config["ref1_RNApeg"],
#         ref2=   config["ref2_RNApeg"],
#         dir_out=dir_out + "/{sample}/" + squence_type + "/fusion/RNApeg_50M",
#         dir_bam=dir_out + "/{sample}/" + squence_type + "/fusion/MarkDup_50M",
#         dir_ref=config["dir_ref_RNApeg"],
#     threads: 16
#     shell:
#         '''
#         singularity run --containall \
#         --bind {params.dir_ref},{params.dir_bam},{params.dir_out}:/results docker://mnedmonson/public:rnapeg RNApeg.sh \
#         -b {input.bam} -f {params.ref1} -r {params.ref2}  &>{log}
#         touch {output.RNApeg}
#         '''
#
# rule Cicero_50M:
#     input:
#         bam=        dir_out + "/{sample}/" + squence_type + "/fusion/MarkDup_50M/{sample}.bam",
#         bai=        dir_out + "/{sample}/" + squence_type + "/fusion/MarkDup_50M/{sample}.bam.bai",
#         junctions=  dir_out + "/{sample}/" + squence_type + "/RNApeg/{sample}.bam.junctions.tab.shifted.tab",
#         bmk=        dir_out + "/{sample}/" + squence_type + "/log/{sample}.RNApeg.bmk"
#     output:     dir_out + "/{sample}/" + squence_type + "/log/{sample}.Cicero_50M.bmk"
#     log:        dir_out + "/{sample}/" + squence_type + "/Cicero_50M/{sample}.Cicero.log"
#     benchmark:  dir_out + "/{sample}/" + squence_type + "/log/{sample}.Cicero_50M.bmk"
#     params:
#         ref=        config["ref_Cicero"],
#         dir_cicero= dir_out + "/{sample}/" + squence_type + "/Cicero_50M/",
#         dir_temp=   dir_out + "/{sample}/" + squence_type + "/fusion/Cicero_temp",
#         dir_fusion= dir_out + "/{sample}/" + squence_type + "/fusion/Cicero_temp/CICERO_DATADIR/{sample}",
#         dir_star=dir_out + "/{sample}/" + squence_type + "/fusion/Star_50M/",
#         dir_remove= dir_out + "/{sample}/" + squence_type + "/fusion"
#     threads: 32
#     shell:
#         '''
#         if [ -d {params.dir_temp} ]; then rm -rf {params.dir_temp}; fi
#         singularity exec --bind {params.ref} /packages/singularity-images/cicero_0.3.0p2.sif Cicero.sh \
#         -n {threads} -b {input.bam} -g GRCh38_no_alt -r {params.ref} -j {input.junctions} -s 2 -c 10 -o {params.dir_temp} &>{log}
#         mv {params.dir_fusion}/*.txt {params.dir_cicero}
#         rm -rf {params.dir_remove}
#         '''

