
rule bwa:
    input:
        ref=config["ref_bwa"],
        fq=lambda wildcards: expand(dir_in + '/{i}',i=fq[wildcards.bamTmp])
    output:
        bam=dir_out+ "/temp/" + squence_type + "/bwa/{bamTmp}.bam",
    log:
        dir_out + "/temp/" +squence_type+ "/log/{bamTmp}_bwa.log",
    benchmark:
        dir_out + "/temp/" + squence_type+ "/log/{bamTmp}_bwa.bmk",
    params:
        rg=lambda wildcards: rg[wildcards.bamTmp],
    threads:
        config["threads_bwa"]
    shell:
        '''
        bwa mem -t {threads} -R '{params.rg}' {input.ref} {input.fq} 2>{log}| samtools view -@ {threads} -Sb - > {output} 2>>{log}
        '''

rule samtools_merge:
    input:
        lambda wildcards: [dir_out + "/temp/" + squence_type + "/bwa/" + bam for bam in
                            sample_tmp_dict[wildcards.sample]]
    output:
        bam=dir_out + "/{sample}/" + squence_type + "/temp/samtools_merge/{sample}.bam",
    log:
        dir_out + "/{sample}/" + squence_type + "/log/samtools_merge.log",
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/samtools_merge.bmk",
    params:
        tmp_n=lambda wildcards: len(sample_tmp_dict[wildcards.sample]),
    shell:
        '''
        if [ {params.tmp_n} -gt 1 ] ; then
        samtools merge {output.bam} {input} &>{log}
        else
        ln -s {input} {output.bam}
        fi
        '''

rule samtools_sort:
    input:
        bam=rules.samtools_merge.output.bam,
    output:
        bam=dir_out + "/{sample}/" + squence_type + "/temp/samtools_sort/{sample}.bam",
        bai=dir_out + "/{sample}/" + squence_type + "/temp/samtools_sort/{sample}.bai"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/samtools_sort.log",
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/samtools_sort.bmk",
    params:
        outdir=dir_out + "/{sample}/" + squence_type + "/temp/samtools_sort/"
    threads:
        config["threads_bwa"]
    shell:
        '''
        samtools sort -m 10G -@ {threads} -O bam -T {params.outdir} -o {output.bam} {input.bam}
        samtools index {output.bam} {output.bai}
        '''

rule MarkDup:
    input:
        bam=rules.samtools_sort.output.bam,
        bmk=rules.samtools_sort.benchmark
    output:
        bam=dir_out + "/{sample}/" + squence_type + "/temp/markDup/{sample}.bam",
        bai=dir_out + "/{sample}/" + squence_type + "/temp/markDup/{sample}.bai",
        metrics=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.MarkDup.metrics.txt",
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}_markDup.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}_markDup.bmk"
    shell:
        '''
        gatk.422 MarkDuplicates -I {input.bam} -O {output.bam} -M {output.metrics} &> {log}
        samtools index {output.bam} {output.bai} &>> {log}
        '''

rule RealignerTargetCreator:
    input:
        bam=rules.MarkDup.output.bam,
        bai=rules.MarkDup.output.bai,
        ref=config["ref_gatk"],
        known=config["known_indelReligner"],
        bmk=rules.MarkDup.benchmark
    output:
        intervals=dir_out + "/{sample}/" + squence_type + "/temp/indelRealigner/{sample}_indelRealigner.intervals",
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}_RealignerTargetCreator.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}_RealignerTargetCreator.bmk"
    params:
        gatk=config["gatk_indelRealigner"]
    shell:
         '''
         java-1.8 -jar {params.gatk} \
         -T RealignerTargetCreator \
         -R {input.ref} \
         -known {input.known} \
         -I {input.bam} \
         -o {output.intervals} &> {log}
         '''

rule IndelRealigner:
    input:
        bam=rules.MarkDup.output.bam,
        bai=rules.MarkDup.output.bai,
        ref=config["ref_gatk"],
        known=config["known_indelReligner"],
        realigner_intervals=rules.RealignerTargetCreator.output.intervals,
        bmk=rules.RealignerTargetCreator.benchmark
    output:
        bam=dir_out + "/{sample}/" + squence_type + "/temp/indelRealigner/{sample}.bam",
        bai=dir_out + "/{sample}/" + squence_type + "/temp/indelRealigner/{sample}.bai"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}_IndelRealigner.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}_IndelRealigner.bmk"
    params:
        gatk=config["gatk_indelRealigner"]
    shell:
         '''
         java-1.8 -jar {params.gatk} \
         -T IndelRealigner \
         -R {input.ref} \
         -I {input.bam} \
         --targetIntervals {input.realigner_intervals} \
         -o {output.bam} &> {log}
         '''

rule BaseRecalibrate:
    input:
        bam=rules.IndelRealigner.output.bam,
        bai=rules.IndelRealigner.output.bai,
        ref_fa=config["ref_gatk"],
        bmk=rules.IndelRealigner.benchmark
    output:
        table=dir_out + "/{sample}/" + squence_type + "/temp/baseRecalibrate/{sample}.baseRecalibrate.table",
    log:
        dir_out+"/{sample}/" + squence_type + "/log/{sample}_BaseRecalibrator.log"
    benchmark:
        dir_out+"/{sample}/" + squence_type + "/log/{sample}_BaseRecalibrator.bmk"
    params:
        known=config["known_baseRecalibrator"]
    shell:
        '''
        gatk.422 BaseRecalibrator -I {input.bam} \
        -R {input.ref_fa}  \
        {params.known} \
        -O {output.table} &> {log}
        '''

rule ApplyBQSR:
    input:
        bam=rules.IndelRealigner.output.bam,
        bai=rules.IndelRealigner.output.bai,
        table=rules.BaseRecalibrate.output.table,
        ref_fa=config["ref_gatk"],
        bmk=rules.BaseRecalibrate.benchmark
    output:
        bam=dir_out + "/{sample}/" + squence_type + "/temp/baseRecalibrate/{sample}.bam",
        bai=dir_out + "/{sample}/" + squence_type + "/temp/baseRecalibrate/{sample}.bai",
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}_ApplyBQSR.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}_ApplyBQSR.bmk"
    shell:
        '''
        gatk.422 ApplyBQSR \
        -R {input.ref_fa} \
        -I {input.bam} \
        --bqsr-recal-file {input.table} -O {output.bam} &> {log}
        '''

rule LeftAlignIndels:
    input:
        bam=rules.ApplyBQSR.output.bam,
        bai=rules.ApplyBQSR.output.bai,
        ref_fa=config["ref_gatk"],
        bmk=rules.ApplyBQSR.benchmark
    output:
        bam=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam",
        bai=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam.bai",
        stats = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.chr.readcounts.txt"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}_leftAlignIndels.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}_leftAlignIndels.bmk"
    params:
        gatk=config["gatk"],
        bai=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bai"
    shell:
        '''
        gatk.422 LeftAlignIndels \
        -R {input.ref_fa} \
        -I {input.bam} \
        -O {output.bam} &>{log}
        mv {params.bai} {output.bai}
        samtools idxstats {output.bam} | cut -f 1,3 > {output.stats}
        '''

rule get_bw:
    input:
        bam = rules.LeftAlignIndels.output.bam,
        bmk = rules.LeftAlignIndels.benchmark
    output:
        bw=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bw",
    log: dir_out + "/{sample}/" + squence_type + "/log/{sample}.get_bw.log"
    benchmark: dir_out + "/{sample}/" + squence_type + "/log/{sample}.get_bw.bmk"
    shell:
        '''
        bamCoverage -bs 1 -b {input.bam} -o {output.bw} &> {log}
        '''

rule Samtools_Flagstat:
    input:
        bam=rules.LeftAlignIndels.output.bam,
        bai=rules.LeftAlignIndels.output.bai,
        bmk=rules.LeftAlignIndels.benchmark
    output:
        stats=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.mappingStats",
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}_samtools_flagstat.bmk"
    shell:
         '''
         samtools flagstat {input.bam} > {output.stats}
         '''

rule md5:
    input:
        bam=rules.LeftAlignIndels.output.bam,
        bai=rules.LeftAlignIndels.output.bai
    output:
        md5=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam.md5",
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}_md5.bmk"
    shell:
        '''
        md5sum {input.bam} > {output.md5}
        '''

rule cdsDepth:
    input:
        bam=rules.LeftAlignIndels.output.bam,
        bai=rules.LeftAlignIndels.output.bai,
        cds_list=config["cds_bed"]
    output:
        dp=dir_out +"/{sample}/" + squence_type + "/bam/{sample}.dp"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}_cdsDepth.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}_cdsDepth.bmk"
    params:
        cdsDepth=config["cdsDepth"],
    shell:
        '''
        {params.cdsDepth} -f {input.cds_list} -b {input.bam} -o {output.dp} &>{log}
        '''

rule Getpileup:
    input:
        bam=rules.LeftAlignIndels.output.bam,
        bai=rules.LeftAlignIndels.output.bai,
        common_germline_variant=config["common_germline_variant"]
    output:
        table=dir_out + "/{sample}/" + squence_type + "/getpileup/{sample}_pileups.table"
    log:
        dir_out+"/{sample}/"+squence_type+"/log/{sample}_pileupsummaries.log"
    benchmark:
        dir_out+"/{sample}/"+squence_type+"/log/{sample}_pileupsummaries.bmk"
    shell:
        '''
        gatk.422 GetPileupSummaries  \
        -I {input.bam} \
        -V {input.common_germline_variant} \
        -L {input.common_germline_variant} \
        -O {output.table} &> {log}
        '''

rule doneMapping:
    input:
        stats = rules.Samtools_Flagstat.output.stats,
        markdup_m=rules.MarkDup.output.metrics,
        md5=rules.md5.benchmark,
        dp=rules.cdsDepth.benchmark,
        getpileup=rules.Getpileup.output.table,
        get_bw=rules.get_bw.benchmark
    output:
        done=dir_out + "/{sample}/" + squence_type + "/log/doneMapping.txt"
    params:
        dir_tmp=dir_out + "/{sample}/" + squence_type + "/temp",
        mapping_file_prefix=dir_out + "/temp/" + squence_type + "/log/{sample}",
        log_dir=dir_out + "/{sample}/" + squence_type + "/log/"
    shell:
        '''
        mv {params.mapping_file_prefix}* {params.log_dir}
        rm -rf {params.dir_tmp}
        touch {output.done}
        '''

rule qc:
    input:
        expand(dir_out + "/{sample}/" + squence_type + "/log/{sample}_cdsDepth.bmk",sample=samplelist[1])
    output:
        coverage=dir_out + "/id/coverage_" + squence_type + ".tsv",
        qc=dir_out + "/id/qc_" + squence_type + ".tsv",
    benchmark:
        dir_out + "/id/qc_" + squence_type + ".bmk",
    params:
        cds_length=config["cds_length"],
        prefix_dp=dir_out + "/*/" + squence_type + "/bam/*.dp",
        dir_db_out=dir_out + "/depth_temp_" + squence_type,
        prefix_duprate=dir_out + "/*/" + squence_type + "/bam/*.MarkDup.metrics.txt",
        prefix_mappingpercentage=dir_out + "/*/" + squence_type + "/bam/*.mappingStats",
        reads=dir_out + "/id/reads_" + squence_type + ".tsv",
        mapping_percentage=dir_out + "/id/mapping_percentage_" + squence_type + ".tsv",
        duprates=dir_out + "/id/duprates_" + squence_type + ".tsv",
    shell:
        '''
        if [ -d "{params.dir_db_out}" ]; then rm -rf {params.dir_db_out};fi

        mkdir {params.dir_db_out}
        paste <(ls {params.prefix_dp}) \
        <(ls {params.prefix_dp} |awk -F/ '{{print "{params.dir_db_out}/"$NF}}')| awk '$1="ln -s "$1'|bash

        dpStatistics.pl {params.dir_db_out} {output.coverage} {params.cds_length}

        #rm -rf {params.dir_db_out}

        paste <(ls {params.prefix_mappingpercentage} |awk -F/ '{{print $NF}}' |awk -F. '{{print $1}}') \
         <(ls {params.prefix_mappingpercentage}|xargs -i echo " cat {{}}|grep read1|cut -d' ' -f1"|bash) > {params.reads}

        paste <(ls {params.prefix_mappingpercentage} |awk -F/ '{{print $NF}}' |awk -F. '{{print $1}}') \
         <(ls {params.prefix_mappingpercentage}|xargs -i echo "cat {{}} |grep mapped|grep %|cut -d'(' -f2|cut -d' ' -f1"|bash) > {params.mapping_percentage}

          paste <(ls {params.prefix_duprate} |awk -F/ '{{print $NF}}' |awk -F. '{{print $1}}') \
         <(ls {params.prefix_duprate} | xargs  -i echo 'cat {{}} |head -n8|tail -n1|cut -f9'|bash) > {params.duprates}


         cat <(head -n1 {output.coverage}|perl -lane  ' $,="\t";print "Sample\\tReadsNumber\\tMappingPercentage\\tDupRates",@F[1..$#F]') \
         <(join <(join <(join {params.reads} {params.mapping_percentage}) {params.duprates}) \
         <(cat {output.coverage} |grep -v Depth|grep -v Length ) | sed 's/ /\\t/g') > {output.qc}
        '''
