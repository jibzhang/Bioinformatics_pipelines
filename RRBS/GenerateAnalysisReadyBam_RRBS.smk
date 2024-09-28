include: "common_rule.smk"

rule all:
    input:
        expand(dir_out + "/{bamTmp}/" + squence_type + "/log/done_mapping.txt",
               bamTmp=bamTmplist[0]),
        #dir_out + "/id/qc_" + squence_type + ".bmk"

rule trimgalore:
    input:
        fq=lambda wildcards: expand(dir_in + '/{i}',i=fq[wildcards.bamTmp])
    output:
        fq1=dir_out + "/fq_trimed/" + squence_type + "/{bamTmp}_1.fq.gz",
        fq2=dir_out + "/fq_trimed/" + squence_type + "/{bamTmp}_2.fq.gz"
    log:
        dir_out+"/temp/"+squence_type+"/log/{bamTmp}_trimgalore.log"
    benchmark:
        dir_out+"/temp/"+squence_type+"/log/{bamTmp}_trimgalore.bmk"
    params:
        dir_tem=dir_out + "/temp/" + squence_type + "/trimgalore",
        usetrimfq="yes",
        prefix="{bamTmp}",
        paired1=dir_out + "/temp/" + squence_type + "/trimgalore/{bamTmp}_val_1.fq.gz",
        paired2 = dir_out + "/temp/" + squence_type + "/trimgalore/{bamTmp}_val_2.fq.gz",
        fq1_raw=dir_in + "/{bamTmp}_R1.fq.gz",
        fq2_raw=dir_in + "/{bamTmp}_R2.fq.gz"
    shell:
        '''
        if [ ! -d "{params.dir_tem}" ]; then mkdir {params.dir_tem}; fi
        if [ "{params.usetrimfq}" == "yes" ]; then
        trim_galore -q 25 --phred33 --rrbs --paired --non_directional --basename {params.prefix} \
        --length 50  -j 4 --fastqc --gzip --max_n 20 -o {params.dir_tem} {input} &> {log};
        ln -s {params.paired1} {output.fq1};ln -s {params.paired2} {output.fq2};
        else
        ln -s {params.fq1_raw} {output.fq1};ln -s {params.fq2_raw} {output.fq2};
        fi
        '''

rule bismark:
    input:
        fq1=dir_out + "/fq_trimed/" + squence_type + "/{bamTmp}_1.fq.gz",
        fq2=dir_out + "/fq_trimed/" + squence_type + "/{bamTmp}_2.fq.gz"
    output:
        dir_out+ "/temp/" + squence_type + "/bismark/{bamTmp}_pe.bam"
    log:
        dir_out+"/temp/"+squence_type+"/log/{bamTmp}_bismark.log"
    benchmark:
        dir_out+"/temp/"+squence_type+"/log/{bamTmp}_bismark.bmk"
    params:
        prefix="{bamTmp}",
        ref=config["bismark_index"],
        outdir=dir_out+ "/temp/" + squence_type + "/bismark"
    threads:
        config["threads_bowtie2"]
    shell:
        '''
        bismark --genome_folder {params.ref} -1 {input.fq1} -2 {input.fq2} --gzip \
        --basename {params.prefix} -p {threads} -N 0 --bowtie2 \
        --phred33-quals -q --non_directional -o {params.outdir} 2> {log}
        '''

rule samtools_sort:
    input:
        bam=rules.bismark.output
    output:
        bam=dir_out + "/{bamTmp}/" + squence_type + "/samtools_sort/{bamTmp}.bam",
        bai=dir_out + "/{bamTmp}/" + squence_type + "/samtools_sort/{bamTmp}.bai"
    log:
        dir_out + "/{bamTmp}/" + squence_type + "/log/samtools_sort.log",
    benchmark:
        dir_out + "/{bamTmp}/" + squence_type + "/log/samtools_sort.bmk",
    params:
        outdir=dir_out + "/{bamTmp}/" + squence_type + "/samtools_sort/"
    threads:
        6
    shell:
        '''
        samtools sort -m 10G -@ {threads} -O bam -T {params.outdir} -o {output.bam} {input.bam}
        samtools index {output.bam} {output.bai}
        '''

rule bismark_methyextract:
    input:
        bam=rules.bismark.output
    output:
        cov = dir_out+ "/{bamTmp}/" + squence_type + "/bismark_methylextr/{bamTmp}_pe.bismark.cov.gz",
        bedgraph = dir_out+ "/{bamTmp}/" + squence_type + "/bismark_methylextr/{bamTmp}_pe.bedGraph.gz"
    log:
        dir_out + "/{bamTmp}/" + squence_type + "/log/{bamTmp}.bismark_methylextr.log"
    benchmark:
        dir_out + "/{bamTmp}/" + squence_type + "/log/{bamTmp}.bismark_methylextr.bmk"
    params:
        ref=config["bismark_index"],
        outdir=dir_out+ "/{bamTmp}/" + squence_type + "/bismark_methylextr"
    threads:
        config["threads_bismark"]
    shell:
        '''
        bismark_methylation_extractor {input.bam} -p --bedgraph --counts --comprehensive --gzip --report \
        --genome_folder {params.ref} --cytosine_report --remove_spaces --no_overlap --buffer_size 10G --ucsc\
        --parallel {threads} -o {params.outdir} 2> {log}
        '''

rule cover_freq:
    input:
        dir_out+ "/{bamTmp}/" + squence_type + "/bismark_methylextr/{bamTmp}_pe.bismark.cov.gz"
    output:
        dir_out+ "/{bamTmp}/" + squence_type + "/bismark_methylextr/{bamTmp}.freq"
    benchmark:
        dir_out + "/{bamTmp}/" + squence_type + "/log/{bamTmp}.cover_freq.bmk"
    params:
        prefix="{bamTmp}",
        unzip=dir_out+ "/{bamTmp}/" + squence_type + "/bismark_methylextr/{bamTmp}_pe.bismark.cov"
    shell:
        '''
        gunzip {input}
        /coh_labs/jochan/bin/covfreq.sh {params.prefix} {params.unzip} {output}
        '''

rule bdgtobw:
    input:
        bedgraph = rules.bismark_methyextract.output.bedgraph,
        cov = dir_out+ "/{bamTmp}/" + squence_type + "/bismark_methylextr/{bamTmp}_pe.bismark.cov"
    output:
        dir_out+ "/{bamTmp}/" + squence_type + "/bismark_methylextr/{bamTmp}_pe.filtered.bw"
    log:
        dir_out + "/{bamTmp}/" + squence_type + "/log/{bamTmp}.bedgraphtobw.log"
    benchmark:
        dir_out + "/{bamTmp}/" + squence_type + "/log/{bamTmp}.bedgraphtobw.bmk"
    params:
        filter_cov = dir_out + "/{bamTmp}/" + squence_type + "/bismark_methylextr/{bamTmp}_pe.over5.cov",
        sort_bgh = dir_out+ "/{bamTmp}/" + squence_type + "/bismark_methylextr/{bamTmp}_pe.sort.bedGraph",
        filter_bgh= dir_out+ "/{bamTmp}/" + squence_type + "/bismark_methylextr/{bamTmp}_pe.filtered.bedGraph",
        chrom_size = "/ref_genomes/genome/human/GRCh38/GRCh38.chrom.size"
    shell:
        '''
        awk '($5+$6)>=5' {input.cov} | awk 'BEGIN{{OFS="\\t"}} {{ $2=$2-1; print }}' > {params.filter_cov}
        zcat {input.bedgraph} | sort -k1,1 -k2,2n | tail -n +2 > {params.sort_bgh}
        bedtools intersect -a {params.sort_bgh} -b {params.filter_cov} -wa | uniq > {params.filter_bgh}
        bedGraphToBigWig {params.filter_bgh} {params.chrom_size} {output} 2> {log}
        rm {params.sort_bgh}
        '''

rule done:
    input:
        bedgraph=rules.bdgtobw.output,
        biscov=rules.cover_freq.output,
        bam=rules.samtools_sort.output.bam
    output:
        dir_out + "/{bamTmp}/" + squence_type + "/log/done_mapping.txt"
    shell:
        '''
        touch {output}
        '''