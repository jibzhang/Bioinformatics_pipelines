rule trimmomatic:
    input:
        fq=lambda wildcards: expand(dir_in + '/{i}',i=fq[wildcards.bamTmp])
    output:
        fq1 = dir_out + "/temp/fq_trimed/{bamTmp}_1.fq.gz",
        fq2 = dir_out + "/temp/fq_trimed/{bamTmp}_2.fq.gz"
    log:
        dir_out+"/temp/"+squence_type+"/log/trimmomatic/{bamTmp}_trimmomatic.log"
    benchmark:
        dir_out+"/temp/"+squence_type+"/log/trimmomatic/{bamTmp}_trimmomatic.bmk"
    params:
        trimpara1="ILLUMINACLIP:/home/zgu_labs/bin/software/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10",
        trimpara2="SLIDINGWINDOW:4:15 MINLEN:36",
        dir_tem=dir_out + "/temp/" + squence_type + "/trimmomatic",
        usetrimfq="yes",
        paired1=dir_out + "/temp/" + squence_type + "/trimmomatic/{bamTmp}_R1.paired.fq.gz",
        paired2 = dir_out + "/temp/" + squence_type + "/trimmomatic/{bamTmp}_R2.paired.fq.gz",
        unpaired1 = dir_out + "/temp/" + squence_type + "/trimmomatic/{bamTmp}_R1.unpaired.fq.gz",
        unpaired2 = dir_out + "/temp/" + squence_type + "/trimmomatic/{bamTmp}_R2.unpaired.fq.gz",
        fq1_raw=dir_in + "/{bamTmp}_R1.fq.gz",
        fq2_raw=dir_in + "/{bamTmp}_R2.fq.gz"
    threads:
        config["threads_star"]
    shell:
        '''
        if [ ! -d "{params.dir_tem}" ]; then mkdir {params.dir_tem}; fi
        if [ "{params.usetrimfq}" == "yes" ] ; then
        trimmomatic PE -threads {threads} {input} {params.paired1} {params.unpaired1} {params.paired2} {params.unpaired2} {params.trimpara1} {params.trimpara2} &> {log};
        ln -s {params.paired1} {output.fq1};ln -s {params.paired2} {output.fq2};
        else
        ln -s {params.fq1_raw} {output.fq1};ln -s {params.fq2_raw} {output.fq2};
        fi
        '''

rule Star:
    input:
        fq1=rules.trimmomatic.output.fq1,
        fq2=rules.trimmomatic.output.fq2,
        bmk=dir_out + "/temp/" + squence_type + "/log/trimmomatic/{bamTmp}_trimmomatic.bmk"
    output:
        bam=dir_out+"/temp/"+squence_type+"/star/{bamTmp}.bam",
        junction=dir_out+"/temp/"+squence_type+"/{bamTmp}/star/SJ.out.tab"
    log:
        dir_out+"/temp/"+squence_type+"/log/star/{bamTmp}_star.log"
    benchmark:
        dir_out+"/temp/"+squence_type+"/log/star/{bamTmp}_star.bmk"
    params:
        star=config["star"],
        star_ref=config["ref_star"],
        rg=lambda wildcards: rg[wildcards.bamTmp],
        dir_out=dir_out+"/temp/"+squence_type+"/{bamTmp}/star/",
        bam_tmp=dir_out+"/temp/"+squence_type+"/{bamTmp}/star/Aligned.out.bam",
    threads:
        config["threads_star"]
    shell:
        '''
        {params.star} \
        --twopassMode Basic \
        --runThreadN {threads} \
        --readFilesCommand zcat \
        --limitOutSAMoneReadBytes 90000000 \
        --outSAMtype BAM Unsorted \
        --outSAMstrandField intronMotif \
        --outSAMunmapped Within \
        --outSAMmapqUnique 60 \
        --outSAMmultNmax 1 \
        --outReadsUnmapped None \
        --outMultimapperOrder Random \
        --outSAMattributes NH HI AS nM NM MD \
        --outSAMattrRGline {params.rg} \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --winAnchorMultimapNmax 200 \
        --outFilterMismatchNmax 4 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 4 \
        --alignInsertionFlush Right \
        --chimOutType Junctions \
        --chimSegmentMin 12 \
        --chimJunctionOverhangMin 12 \
        --chimSegmentReadGapMax 3 \
        --chimMultimapNmax 100 \
        --chimMultimapScoreRange 10 \
        --chimNonchimScoreDropMin 10 \
        --chimOutJunctionFormat 1 \
        --peOverlapNbasesMin 12 \
        --peOverlapMMp 0.1 \
        --outFileNamePrefix {params.dir_out} \
        --genomeDir {params.star_ref} \
        --genomeLoad NoSharedMemory \
        --readFilesIn {input.fq1} {input.fq2}\
        --sjdbOverhang 100 &> {log}
        mv {params.bam_tmp} {output.bam}
        '''

rule samtools_merge:
    input:
        lambda wildcards: [dir_out+"/temp/"+squence_type+"/star/"+bam for bam in sample_tmp_dict[wildcards.sample]]
    output:
        bam_unsort=dir_out+"/{sample}/"+squence_type+"/temp/samtools_merge/{sample}_unsort.bam"
    log:
        dir_out+"/{sample}/"+squence_type+"/log/samtools_merge.log",
    benchmark:
        dir_out+"/{sample}/"+squence_type+"/log/samtools_merge.bmk",
    params:
        tmp_n=lambda wildcards: len(sample_tmp_dict[wildcards.sample]),
    threads:
        config["threads_star"]
    shell:
        '''
        if [ {params.tmp_n} -gt 1 ] ; then
        samtools merge {output.bam_unsort} {input} &>{log}
        else
        ln -s {input} {output.bam_unsort}
        fi
        '''

rule samtools_sort:
    input:
        bam_unsort = rules.samtools_merge.output.bam_unsort
    output:
        bam = dir_out + "/{sample}/" + squence_type + "/temp/samtools_sort/{sample}.bam",
        bai = dir_out + "/{sample}/" + squence_type + "/temp/samtools_sort/{sample}.bam.bai",
    log:
        dir_out + "/{sample}/" + squence_type + "/log/samtools_sort.log",
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/samtools_sort.bmk",
    params:
        outdir = dir_out + "/{sample}/" + squence_type + "/temp/samtools_sort/"
    threads:
        config["threads_star"]
    shell:
        '''
        samtools sort -m 10G -@ {threads} -O bam -T {params.outdir} -o {output.bam} {input.bam_unsort}
        samtools index {output.bam} {output.bai} &>> {log}
        '''

rule MarkDup:
    input:
        bam=rules.samtools_sort.output.bam,
        bmk=rules.samtools_sort.benchmark
    output:
        bam=dir_out+"/{sample}/"+squence_type+"/bam/{sample}.bam",
        bai=dir_out+"/{sample}/"+squence_type+"/bam/{sample}.bam.bai",
        metrics=dir_out+"/{sample}/"+squence_type+"/bam/{sample}.MarkDup.metrics.txt",
        stats= dir_out + "/{sample}/" + squence_type + "/bam/{sample}.chr.readcounts.txt"
    log:
        dir_out+"/{sample}/"+squence_type+"/log/{sample}_markDup.log"
    benchmark:
        dir_out+"/{sample}/"+squence_type+"/log/{sample}_markDup.bmk"
    params:
        gatk=config["gatk"]
    shell:
        '''
        {params.gatk} MarkDuplicates -I {input.bam} -O {output.bam} -M {output.metrics} &> {log}
        samtools index {output.bam} {output.bai} &>> {log}
        samtools idxstats {output.bam} | cut -f 1,3 > {output.stats}
        '''

# rule pathseq:
#     input:
#         bam =  dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam",
#         human_image = rules.creatimage.output.human_ima,
#         ebv_image = rules.creatimage.output.ebv_ima,
#         ebv_dict = config["EBV_dict"],
#         kmer = rules.buildkmer.output,
#         ebv_ref= config["EBV_genome"],
#         taxonomy = rules.Taxonomy.output
#     output:
#         bam = dir_out + "/{sample}/" + squence_type + "/pathoseq/{sample}.bam",
#         score = dir_out + "/{sample}/" + squence_type + "/pathoseq/{sample}.score.txt",
#         filter_metrics = dir_out + "/{sample}/" + squence_type + "/pathoseq/{sample}.filter_metrics.txt",
#         score_metrics = dir_out + "/{sample}/" + squence_type + "/pathoseq/{sample}.score_metrics.txt"
#     log:
#         dir_out + "/{sample}/" + squence_type + "/log/{sample}_pathoseq.log"
#     benchmark:
#         dir_out + "/{sample}/" + squence_type + "/log/{sample}_pathoseq.bmk"
#     shell:
#         '''
#         gatk.422 PathSeqPipelineSpark --input {input.bam} --kmer-file {input.kmer} --filter-bwa-image {input.human_image} --microbe-bwa-image {input.ebv_image} \
#         --microbe-dict {input.ebv_dict} --taxonomy-file {input.taxonomy} --min-clipped-read-length 50 --min-score-identity 0.90 --identity-margin 0.02 \
#         --scores-output {output.score} --output {output.bam} --filter-metrics {output.filter_metrics} --score-metrics {output.score_metrics}
#         '''

rule make_bigwig:
    input:
        bam=rules.MarkDup.output.bam,
        bai=rules.MarkDup.output.bai
    output:
        dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bw"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.bamcoverage.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.bamcoverage.bmk"
    shell:
         '''
         bamCoverage -b {input.bam} --normalizeUsing RPKM --binSize 50 --smoothLength 200 -p 5 -o {output} 2> {log}
         '''

rule TE_count:
    input:
        bam=rules.MarkDup.output.bam,
        bai=rules.MarkDup.output.bai
    output:
        dir_out + "/{sample}/" + squence_type + "/TEcount/TEcount_out.cntTable"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.TEcount.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.TEcount.bmk"
    params:
        gtf = config["ref_gtf"],
        TE_gtf = "/coh_labs/jochan/reference/hg38_region/hg38_rmsk_TE_20200804.gtf",
        outdir = dir_out + "/{sample}/" + squence_type + "/TEcount"
    shell:
        '''
        TEcount --sortByPos --format BAM --mode multi -b {input.bam} --GTF {params.gtf} --TE {params.TE_gtf} --outdir {params.outdir} 2> {log}
        '''