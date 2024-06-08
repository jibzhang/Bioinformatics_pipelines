rule fastqc:
    input:
        fq=dir_in + '/{fqname}.fastq.gz'
    output:
        zip=dir_out + "/temp/" + squence_type + "/fastqc/{fqname}_fastqc.zip"
    log:
        dir_out + "/temp/" + squence_type + "/log/fastqc/{fqname}_fastqc.log"
    benchmark:
        dir_out + "/temp/" + squence_type + "/log/fastqc/{fqname}_fastqc.bmk"
    params:
        outdir=dir_out + "/temp/" + squence_type + "/fastqc/"
    shell:
        """
        if [ ! -d "{params.outdir}" ]; then mkdir {params.outdir}; fi
        fastqc {input} -o {params.outdir} &> {log}
        """

rule multiqc:
    input:
        expand(dir_out + "/temp/" + squence_type + "/log/fastqc/{fqname}_fastqc.bmk",fqname=fqnamelist)
    output:
        dir_out + "/fastqc_multiqc/"+multiqc_name+".html"
    params:
        name=multiqc_name,
        data_dir=dir_out + "/temp/" + squence_type + "/fastqc/",
        outdir=dir_out + "/fastqc_multiqc"
    shell:
        '''
        multiqc {params.data_dir} --filename {params.name} --outdir {params.outdir}
        '''

rule trimmomatic:
    input:
        fq=lambda wildcards: expand(dir_in + '/{i}.fastq.gz',i=fq[wildcards.bamTmp])
    output:
        fq=dir_out + "/fq_trimed/" + squence_type + "/{bamTmp}.fq.gz"
    log:
        dir_out+"/temp/"+squence_type+"/log/trimmomatic/{bamTmp}_trimmomatic.log"
    benchmark:
        dir_out+"/temp/"+squence_type+"/log/trimmomatic/{bamTmp}_trimmomatic.bmk"
    params:
        trimpara="ILLUMINACLIP:/home/zgu_labs/bin/software/trimmomatic/adapters/TruSeq2-PE.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:30",
        dir_tem=dir_out + "/temp/" + squence_type + "/trimmomatic",
        usetrimfq="yes",
        fq=dir_out + "/temp/" + squence_type + "/trimmomatic/{bamTmp}.fq.gz"
    threads:
        config["threads_trimmomatic"]
    shell:
        '''
        if [ ! -d "{params.dir_tem}" ]; then mkdir {params.dir_tem}; fi
        if [ "{params.usetrimfq}" == "yes" ] ; then
        trimmomatic SE -threads {threads} {input} {params.fq} {params.trimpara} &> {log};
        ln -s {params.fq} {output.fq}
        else
        ln -s {input.fq} {output.fq}
        fi
        '''

rule bowtie2:
    input:
        fq=rules.trimmomatic.output.fq
    output:
        bam=dir_out+ "/temp/" + squence_type + "/bowtie2/{bamTmp}.bam",
    log:
        dir_out + "/temp/" +squence_type+ "/log/{bamTmp}_bowtie2.log",
    benchmark:
        dir_out + "/temp/" + squence_type+ "/log/{bamTmp}_bowtie2.bmk",
    params:
        ref=config["bowtie2_index"],
    threads:
        config["threads_bowtie2"]
    shell:
        '''
        bowtie2 --mm -p {threads} -x {params.ref} -U {input.fq} 2>{log}| samtools view -@ {threads} -Sb - > {output.bam} 2>>{log}
        '''

rule samtools_merge:
    input:
        lambda wildcards: [dir_out + "/temp/" + squence_type + "/bowtie2/" + bam for bam in
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
        config["threads_bowtie2"]
    shell:
        '''
        samtools sort -m 10G -@ {threads} -O bam -T {params.outdir} -o {output.bam} {input.bam}
        samtools index {output.bam} {output.bai}
        '''

rule Initial_Stats:
    input:
        bam = rules.samtools_sort.output.bam
    output:
        stats=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.initial.stats",
        idxstats=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.initial.idxstats",
        flagstats=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.initial.flagstats"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.initial_stats.bmk"
    threads:
        config["threads_bowtie2"]
    params:
        ref=config["ref_bwa"]
    shell:
        '''
        samtools stats -@ {threads} -r {params.ref} {input.bam} > {output.stats}
        samtools flagstat -@ {threads} {input.bam} > {output.flagstats}
        samtools idxstats {input.bam} > {output.idxstats}
        '''

rule MarkDup:
    input:
        bam=rules.samtools_sort.output.bam,
        bai=rules.samtools_sort.output.bai
    output:
        bam = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.MarkDup.bam",
        metrics = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.MarkDup.metrics",
        index = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.MarkDup.bai"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}_markDup.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}_markDup.bmk"
    threads:
        config["threads_bowtie2"]
    params:
        ref=config["ref_bwa"]
    shell:
        '''
        gatk.422 MarkDuplicates -I {input.bam} -O {output.bam} -M {output.metrics} --VALIDATION_STRINGENCY LENIENT &> {log}
        samtools index {output.bam} {output.index}
        '''

rule Markdup_Stats:
    input:
        bam = rules.MarkDup.output.bam
    output:
        stats=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.MarkDup.stats",
        idxstats=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.MarkDup.idxstats",
        flagstats=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.MarkDup.flagstats"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.markdup_stats.bmk"
    threads:
        config["threads_bowtie2"]
    params:
        ref=config["ref_bwa"]
    shell:
        '''
        samtools stats -@ {threads} -r {params.ref} {input.bam} > {output.stats}
        samtools flagstat -@ {threads} {input.bam} > {output.flagstats}
        samtools idxstats {input.bam} > {output.idxstats}
        '''

rule samtools_filter:
      input:
        bam=rules.MarkDup.output.bam
      output:
        bam=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.samtools_filter.bam"
      log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.samtools_filter.log",
      benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.samtools_filter.bmk",
      params:
        noblacklist = "/home/jibzhang/reference/hg38_region/hg38_noblacklist.v3.bed",
        json = "/coh_labs/jochan/pipeline/ChIP_seq/bamtools_filter_se.json"
      shell:
        '''
        module load BamTools
        samtools view -F 1804 -q 30 -b {input.bam} -L {params.noblacklist} | bamtools filter -out {output.bam} -script {params.json}
        '''

rule final_bam:
    input:
        bam=rules.samtools_filter.output.bam
    output:
        bam=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam",
        bai=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bai"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.final_bam.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.final_bam.bmk"
    threads:
        config["threads_bowtie2"]
    params:
        outdir = dir_out + "/{sample}/" + squence_type + "/bam/"
    shell:
        '''
        samtools sort -@ {threads} -o {output.bam} -T {params.outdir} {input.bam} 2>{log}
        samtools index -@ {threads} {output.bam} {output.bai} 2>>{log}
        rm {input.bam}
        '''

rule Final_Stats:
    input:
        bam = rules.final_bam.output.bam
    output:
        stats=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.final.stats",
        idxstats=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.final.idxstats",
        flagstats=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.final.flagstats"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.final_stats.bmk"
    threads:
        config["threads_bowtie2"]
    params:
        ref=config["ref_bwa"]
    shell:
         '''
        samtools stats -@ {threads} -r {params.ref} {input.bam} > {output.stats}
        samtools flagstat -@ {threads} {input.bam} > {output.flagstats}
        samtools idxstats {input.bam} > {output.idxstats}
         '''

rule lib_complexity:
    input:
        bam=rules.MarkDup.output.bam,
    output:
        qc=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.pbc.qc",
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.lib_complex.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.lib_complex.bmk"
    params:
        bam = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.srt_temp.bam",
        outdir = dir_out + "/{sample}/" + squence_type + "/bam/"
    threads:
        config["threads_bowtie2"]
    shell:
        '''
        samtools sort -n -@ {threads} -o {params.bam} -T {params.outdir} {input.bam} 2>{log}
        bedtools bamtobed -i {params.bam} | awk 'BEGIN{{OFS="\\t"}}{{print $1,$2,$3,$6}}' | grep -v 'chrM\|chrUn' | sort | uniq -c | awk 'BEGIN{{mt=0;m0=0;m1=0;m2=0}} ($1==1){{m1=m1+1}} ($1==2){{m2=m2+1}} {{m0=m0+1}} {{mt=mt+$1}} END{{printf "%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}}' > {output.qc}
        rm {params.bam}
        '''

rule multimetric:
    input:
        bam = rules.final_bam.output.bam,
        ref = config["ref_bwa"]
    output:
        dir_out + "/{sample}/" + squence_type + "/bam/{sample}.CollectMultipleMetrics.alignment_summary_metrics",
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.multimetric.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.multimetric.bmk"
    params:
        prefix = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.CollectMultipleMetrics",
    shell:
        '''
        gatk.422 CollectMultipleMetrics -I {input.bam} -O {params.prefix} -R {input.ref} &> {log}
        '''

rule phantompeak:
    input:
        bam = rules.final_bam.output.bam
    output:
        plot=dir_out + "/{sample}/" + squence_type + "/phantompeak/{sample}.spp.pdf",
        rdata=dir_out + "/{sample}/" + squence_type + "/phantompeak/{sample}.spp.Rdata",
        result=dir_out + "/{sample}/" + squence_type + "/phantompeak/{sample}.run_spp.out"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.phantompeak.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.phantompeak.bmk"
    params:
        dir=dir_out + "/{sample}/" + squence_type + "/phantompeak/"
    shell:
        '''
        Rscript /home/jibzhang/miniconda3/bin/run_spp.R -c={input.bam} -tmpdir={params.dir} -savp={output.plot} -savd={output.rdata} -rf -p=4 -out={output.result} 2> {log}
        '''

rule genomecov:
    input:
        bam = rules.final_bam.output.bam,
        flagstat = rules.Final_Stats.output.flagstats
    output:
        scale_factor = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.scale_factor.txt",
        bedgraph = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bedGraph"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.genomecov.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.genomecov.bmk"
    shell:
         '''
         SCALE_FACTOR=$(grep '[0-9] mapped (' {input.flagstat} | awk '{{print 1000000/$1}}')
         echo $SCALE_FACTOR > {output.scale_factor}
         bedtools genomecov -ibam {input.bam} -bg -scale $SCALE_FACTOR -fs 200 | sort -T '.' -k1,1 -k2,2n > {output.bedgraph}
         '''

rule bedgraph_to_bigwig:
    input:
        bedgraph = rules.genomecov.output.bedgraph,
        chrom_size = config["chrom_sizes"]
    output:
        dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bigWig"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.bigwig.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.bigwig.bmk"
    shell:
         '''
         bedGraphToBigWig {input.bedgraph} {input.chrom_size} {output}
         rm {input.bedgraph}
         '''

rule make_bigwig:
    input:
        bam=rules.final_bam.output.bam,
        bai=rules.final_bam.output.bai
    output:
        dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bw",
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.bamcoverage.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.bamcoverage.bmk"
    shell:
         '''
         bamCoverage -b {input.bam} --normalizeUsing RPKM --binSize 50 --smoothLength 200 -p 5 -o {output} 2> {log}
         '''

rule done_mapping:
    input:
        phantompeak = rules.phantompeak.benchmark,
        initialstat = rules.Initial_Stats.benchmark,
        finalstat = rules.Final_Stats.benchmark,
        bigwig = rules.make_bigwig.benchmark,
        genomecov = rules.bedgraph_to_bigwig.benchmark,
        complexity = rules.lib_complexity.benchmark,
        multimetric = rules.multimetric.benchmark
    output:
        dir_out + "/{sample}/" + squence_type + "/log/done_mapping.txt"
    shell:
        '''
        touch {output}
        '''