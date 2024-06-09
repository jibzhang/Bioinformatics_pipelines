rule cutadapter:
    input:
        fq=lambda wildcards: expand(dir_in + '/{i}.fastq.gz',i=fq[wildcards.bamTmp])
    output:
        fq1=dir_out + "/fq_trimed/" + squence_type + "/{bamTmp}_1.fq.gz",
        fq2=dir_out + "/fq_trimed/" + squence_type + "/{bamTmp}_2.fq.gz"
    log:
        dir_out+"/temp/"+squence_type+"/log/{bamTmp}_cutadapter.log"
    benchmark:
        dir_out+"/temp/"+squence_type+"/log/{bamTmp}_cutadapter.bmk"
    params:
        dir_tem=dir_out + "/temp/" + squence_type + "/cutadapt",
        usetrimfq="yes",
        paired1=dir_out + "/temp/" + squence_type + "/cutadapt/{bamTmp}_R1.paired.fq.gz",
        paired2 = dir_out + "/temp/" + squence_type + "/cutadapt/{bamTmp}_R2.paired.fq.gz",
        fq1_raw=dir_in + "/{bamTmp}_R1.fq.gz",
        fq2_raw=dir_in + "/{bamTmp}_R2.fq.gz"
    shell:
        '''
        if [ ! -d "{params.dir_tem}" ]; then mkdir {params.dir_tem}; fi
        if [ "{params.usetrimfq}" == "yes" ] ; then
        cutadapt -q 25 --minimum-length=5 --pair-filter=any \
        -a CTGTCTCTTATA -A CTGTCTCTTATA -e 0.2 \
        -o {params.paired1} -p {params.paired2} {input} &> {log};
        ln -s {params.paired1} {output.fq1};ln -s {params.paired2} {output.fq2};
        else
        ln -s {params.fq1_raw} {output.fq1};ln -s {params.fq2_raw} {output.fq2};
        fi
        '''

rule bowtie2:
    input:
        fq1=rules.cutadapter.output.fq1,
        fq2=rules.cutadapter.output.fq2
    output:
        bam=dir_out+ "/temp/" + squence_type + "/bowtie2/{bamTmp}.bam",
    log:
        dir_out + "/temp/" +squence_type+ "/log/{bamTmp}_bowtie2.log",
    benchmark:
        dir_out + "/temp/" + squence_type+ "/log/{bamTmp}_bowtie2.bmk",
    params:
        ref=config["bowtie2_index"],
        para_bowtie="-X 2000"
    threads:
        config["threads_bowtie2"]
    shell:
        '''
        bowtie2 {params.para_bowtie} --mm -p {threads} -x {params.ref} -1 {input.fq1} -2 {input.fq2} 2>{log}| samtools view -@ {threads} -Sb - > {output} 2>>{log}
        '''

rule samtools_merge:
    input:
        lambda wildcards: [dir_out + "/temp/" + squence_type + "/bowtie2/" + bam for bam in
                            sample_tmp_dict[wildcards.sample]]
    output:
        bam=dir_out + "/{sample}/" + squence_type + "/temp/samtools_merge/{sample}.bam",
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/samtools_merge.bmk",
    params:
        tmp_n=lambda wildcards: len(sample_tmp_dict[wildcards.sample]),
    shell:
        '''
        if [ {params.tmp_n} -gt 1 ] ; then
        samtools merge {output.bam} {input}
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
        bam=rules.MarkDup.output.bam,
        idx=rules.Markdup_Stats.output.idxstats
      output:
        bam=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.samtools_filter.bam"
      benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.samtools_filter.bmk",
      params:
        noblacklist = "/home/jibzhang/reference/hg38_region/hg38_noblacklist.v3.bed",
        json = "/coh_labs/jochan/pipeline/ChIP_seq/bamtools_filter_pe.json",
        bam_clean = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.chr_rm.bam"
      shell:
        '''
        module load BamTools
        cat {input.idx} | cut -f 1 | grep -v "random$\|^chrUn\|^chrM$\|^chrX$|^chrY$|^chrEBV$" | xargs samtools view {input.bam} -b > {params.bam_clean}
        samtools view -F 1804 -q 30 -f 2 -b {params.bam_clean} -L {params.noblacklist} | bamtools filter -out {output.bam} -script {params.json}
        rm {params.bam_clean}
        '''

rule remove_orphan:
      input:
        bam=rules.samtools_filter.output.bam
      output:
        namesrt_bam = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.name_sorted.bam",
        bam=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.orphan_removed.bam"
      log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.remove_orphan.log",
      benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.remove_orphan.bmk",
      threads:
        config["threads_bowtie2"]
      params:
        dir = dir_out + "/{sample}/" + squence_type + "/bam/"
      shell:
        '''
        samtools sort -n -@ {threads} -o {output.namesrt_bam} -T {params.dir} {input.bam} 2>{log}
        /home/jibzhang/bin/bampe_rm_orphan.py {output.namesrt_bam} {output.bam} 2>>{log}
        '''

rule format:
      input:
        bam = rules.remove_orphan.output.namesrt_bam
      output:
        minibed = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.minimal.bedpe"
      benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.format.bmk"
      params:
        bam = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.fixmate.bam",
        bedpe = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.fixed.bedpe",
        t5bed = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.tn5.bedpe"
      shell:
        '''
        samtools fixmate {input.bam} {params.bam}
        samtools view -bf 0x2 {params.bam} | bedtools bamtobed -i stdin -bedpe > {params.bedpe}
        bash /home/jibzhang/bin/bedpeTn5shift.sh {params.bedpe} > {params.t5bed}
        bash /home/jibzhang/bin/bedpeMinimalConvert.sh {params.t5bed} > {output.minibed}
        rm {params.bam} {params.bedpe} {params.t5bed}
        '''

rule final_bam:
    input:
        bam=rules.remove_orphan.output.bam
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
         bedtools genomecov -ibam {input.bam} -bg -scale $SCALE_FACTOR -pc | sort -T '.' -k1,1 -k2,2n > {output.bedgraph}
         '''

rule bedgraph_to_bigwig:
    input:
        bedgraph = rules.genomecov.output.bedgraph,
        chrom_size = config["chrom_sizes"]
    output:
        dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bigWig"
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

rule Compute_matrix:
    input:
        bw = rules.bedgraph_to_bigwig.output,
        bed = config["gtf_bed"],
        tss = config["tss_bed"]
    output:
        scale_plotmatrix = dir_out + "/{sample}/" + squence_type + "/deeptools/{sample}.scale.matrix.gz",
        tss_plotmatrix = dir_out + "/{sample}/" + squence_type + "/deeptools/{sample}.tss.matrix.gz"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.computmatrix.bmk",
    params:
        scale_matrix = dir_out + "/{sample}/" + squence_type + "/deeptools/{sample}.scale.mat.tab",
        tss_matrix = dir_out + "/{sample}/" + squence_type + "/deeptools/{sample}.tss.mat.tab"
    threads:
        8
    shell:
        '''
        computeMatrix reference-point -S {input.bw} -R {input.bed} --skipZeros\
        -o {output.scale_plotmatrix} --outFileNameMatrix {params.scale_matrix} -p {threads}

        computeMatrix reference-point -S {input.bw} -R {input.tss} --skipZeros\
        -o {output.tss_plotmatrix} --outFileNameMatrix {params.tss_matrix} -p {threads}
        '''

rule plot_profile:
    input:
        matrix = rules.Compute_matrix.output.scale_plotmatrix
    output:
        plot = dir_out+ "/{sample}/" + squence_type + "/deeptools/{sample}.plotProfile.png",
        tab = dir_out + "/{sample}/" + squence_type + "/deeptools/{sample}.plotProfile.tab"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.plotProfile.bmk",
    shell:
        '''
        plotProfile -m {input.matrix} -o {output.plot} --outFileNameData {output.tab}
        '''

rule plot_heatmap:
    input:
        matrix = rules.Compute_matrix.output.tss_plotmatrix
    output:
        plot = dir_out+ "/{sample}/" + squence_type + "/deeptools/{sample}.plotHeatmap.png",
        tab = dir_out + "/{sample}/" + squence_type + "/deeptools/{sample}.plotHeatmap.tab"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.heatmap.bmk",
    shell:
        '''
        plotHeatmap -m {input.matrix} -o {output.plot} --outFileNameMatrix {output.tab}
        '''

rule MACS2_broad:
    input:
        bed = rules.format.output.minibed
    output:
        peak=dir_out + "/{sample}/" + squence_type + "/macs2_broad/{sample}_peaks.broadPeak"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.macs2_broad.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.macs2_broad.bmk"
    params:
        prefix = "{sample}",
        dir = dir_out + "/{sample}/" + squence_type + "/macs2_broad"
    shell:
        '''
        macs2 callpeak -f BEDPE -t {input} -n {params.prefix}  --broad -g 2805665311 --broad-cutoff 0.05 --keep-dup all --outdir {params.dir} 2> {log}
        '''

rule MACS2:
    input:
        bed = rules.format.output.minibed
    output:
        peak=dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_peaks.narrowPeak",
        summit=dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_summits.bed"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.macs2.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.macs2.bmk"
    params:
        prefix = "{sample}",
        dir = dir_out + "/{sample}/" + squence_type + "/macs2",
        fq_n = lambda wildcards:len(fq[wildcards.sample])
    shell:
        '''
        if [ {params.fq_n} -gt 1 ] ; then
        macs2 callpeak -g 2805665311 -f BEDPE -t {input.bed} -n {params.prefix} --keep-dup all --outdir {params.dir} 2> {log}
        else
        macs2 callpeak -g 2805665311 -f BED -t {input.bed} -n {params.prefix} --keep-dup all --outdir {params.dir} 2> {log}
        fi
        '''

rule filter:
    input:
        peaks=rules.MACS2.output.peak,
        blacklist="/home/jibzhang/reference/hg38_region/hg38-ChIPblacklist.v3.bed"
    output:
        filtered = dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_filtered_peaks.bed"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.filter_bed.bmk"
    params:
        blacklist_filtered = dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_peaks.blaclist_filtered.bed"
    shell:
        '''
        bedtools intersect -v -a {input.peaks} -b {input.blacklist} > {params.blacklist_filtered}
        awk -F "\\t" '{{if(($7 > 4) && ($5 > 110)) {{print}} }}' {params.blacklist_filtered} > {output.filtered}
        '''

rule HOMER:
    input:
        narrow = rules.filter.output.filtered,
        broad = rules.MACS2_broad.output.peak,
        fasta = config["ref_bwa"],
        gtf = config["ref_gtf"]
    output:
        narrow=dir_out + "/{sample}/" + squence_type + "/Homer/{sample}.annotation.txt",
        broad=dir_out + "/{sample}/" + squence_type + "/Homer/{sample}.broad.annotation.txt",
        motif=dir_out + "/{sample}/" + squence_type + "/Homer/motif/homerResults/motif1.motif"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.homer.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.homer.bmk"
    params:
        motifdir = dir_out + "/{sample}/" + squence_type + "/Homer/motif",
        KORNAseq = "/coh_labs/jochan/PRDM1_NK/PRDM1KO_vs_WT.txt",
        overexpr = "/coh_labs/jochan/PRDM1_NK/PRDM1_overexpres_vs_low.txt"
    shell:
        '''
        annotatePeaks.pl {input.narrow} {input.fasta} -gtf {input.gtf} -gene {params.KORNAseq} {params.overexpr} -cpu 4 > {output.narrow}
        annotatePeaks.pl {input.broad} hg38 -gene {params.KORNAseq} {params.overexpr} -cpu 4 > {output.broad}
        findMotifsGenome.pl {input.narrow} hg38 {params.motifdir} -size 50 &> {log}
        '''

rule done_individual:
    input:
        initialstat = rules.Initial_Stats.benchmark,
        markdupstat = rules.Markdup_Stats.benchmark,
        bigwig = rules.make_bigwig.benchmark,
        macs2 = rules.MACS2.benchmark,
        multimetric = rules.multimetric.benchmark,
        plotheatmap = rules.plot_heatmap.benchmark,
        plotprofile = rules.plot_profile.benchmark,
        Homer = rules.HOMER.benchmark
    output:
        dir_out + "/{sample}/" + squence_type + "/log/done_mapping.txt"
    shell:
        '''
        touch {output}
        '''

rule consensus_peak:
    input:
        expand(rules.filter.output.filtered,sample=samplelist)
    output:
        bed = dir_out + "/concensus_peak/concensus.bed",
        saf = dir_out + "/concensus_peak/concensus.saf"
    params:
        D2nofeeder_peaks = dir_out + "/*/" + squence_type + "/macs2/*D2nofeeder_peaks.narrowPeak",
        D6nofeeder_peaks = dir_out + "/*/" + squence_type + "/macs2/*D6nofeeder_peaks.narrowPeak",
        D13feeder_peaks = dir_out + "/*/" + squence_type + "/macs2/*D13feeder_peaks.narrowPeak",
        D28feeder_peaks = dir_out + "/*/" + squence_type + "/macs2/*D28feeder_peaks.narrowPeak",
        merged = dir_out + "/concensus_peak/merged_peaks.txt",
        boolean = dir_out + "/concensus_peak/boolean.txt",
        intersect = dir_out + "/concensus_peak/boolean.intersect.txt",
        plot = dir_out + "/concensus_peak/boolean.intersect.plot.pdf"
    benchmark:
        dir_out + "/log/concensus_peak.bmk"
    shell:
        '''
        sort -k1,1 -k2,2n {params.D2nofeeder_peaks} {params.D6nofeeder_peaks} {params.D13feeder_peaks} {params.D28feeder_peaks} | mergeBed -c 2,3,4,5,6,7,8,9,10 -o collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse > {params.merged}
        /home/jibzhang/bin/macs2_merged_expand.py {params.merged} F31_D2nofeeder,M29_D2nofeeder,F31_D6nofeeder,M29_D6nofeeder,F37_D6nofeeder,M37_D6nofeeder,M25_D6nofeeder,F31_D13feeder,F26_D13feeder,M29_D13feeder,M37_D13feeder,F31_D28feeder,M37_D28feeder,M29_D28feeder_peaks {params.boolean} --min_replicates 2 --is_narrow_peak
        awk -v FS='\\t' -v OFS='\\t' 'FNR > 1 {{ print $1, $2, $3, $4, "0", "+" }}' {params.boolean} > {output.bed}
        echo -e "GeneID\\tChr\\tStart\\tEnd\\tStrand" > {output.saf}
        awk -v FS='\\t' -v OFS='\\t' 'FNR > 1 {{ print $4, $1, $2, $3, "+" }}' {params.boolean} >> {output.saf}
        /home/jibzhang/bin/R_script/plot_peak_intersect.r -i {params.intersect} -o {params.plot}
        '''

rule featureCounts:
    input:
        annotation = rules.consensus_peak.output.saf,
        bam = expand(rules.final_bam.output.bam,sample=samplelist)
    output:
        dir_out + "/featureCount/featureCounts.txt"
    benchmark:
        dir_out + "/log/featurecounts.bmk"
    log:
        dir_out + "/log/featurecounts.log"
    threads:
        6
    shell:
        '''
        featureCounts -F SAF -p -O --fracOverlap 0.2 -T {threads} -a {input.annotation} -o {output} {input.bam} &> {log}
        '''

rule plotfingerprint:
    input:
        expand(rules.final_bam.output.bam,sample=samplelist)
    output:
        plot = dir_out + "/deeptools/plotFingerprint.pdf",
        rawcount = dir_out + "/deeptools/plotFingerprint.raw.txt",
        qcmetrics = dir_out + "/deeptools/plotFingerprint.qcmetrics.txt"
    benchmark:
        dir_out + "/log/fingerprint.bmk"
    threads:
        6
    shell:
        '''
        plotFingerprint --bamfiles {input} --minMappingQuality 30 --skipZeros -o {output.plot} --outRawCounts {output.rawcount} --outQualityMetrics {output.qcmetrics} -p {threads}
        '''

rule heatmap_consensus:
    input:
        bw = expand(rules.bedgraph_to_bigwig.output,sample=samplelist),
        peak = rules.consensus_peak.output.bed
    output:
        dir_out+"/deeptools/heatmap/heatmap_consensus.png"
    benchmark:
        dir_out+"/deeptools/heatmap/heatmap_Allpeak_merge.bmk"
    params:
        gene_matrix = dir_out+"/deeptools/heatmap/heatmap_Allpeak_merge.matrix.gz",
        labels = lambda wildcards, input: [bw.split("/")[4] for bw in input["bw"]]
    threads:
        16
    shell:
        '''
        computeMatrix reference-point -S {input.bw} -R {input.peak} \
        --samplesLabel {params.labels} \
        --referencePoint center --skipZeros -p {threads} \
        -o {params.gene_matrix}

        plotHeatmap -m {params.gene_matrix} -out {output} \
        --sortUsing sum --startLabel "Peak Start" --endLabel "Peak End" --xAxisLabel "" --regionsLabel "PeaksInAllSamples"
        '''

rule correlation:
    input:
        bw = expand(rules.make_bigwig.output,sample=samplelist),
        peak = rules.consensus_peak.output.bed
    output:
        pca=dir_out+"/deeptools/correlation/PCA_readcounts.png",
        correlation = dir_out+"/deeptools/correlation/Pearson_corr_heatmap.png",
        mtrx= dir_out+"/deeptools/correlation/Pearson_corr_readCounts.tab"
    benchmark:
        dir_out+"/deeptools/correlation/PCA_correlation.bmk"
    params:
        npz=dir_out+"/deeptools/correlation/results.npz",
    shell:
        '''
        multiBigwigSummary BED-file -b {input.bw} -o {params.npz} --BED {input.peak}
        plotPCA -in {params.npz} -o {output.pca} -T "PCA of read counts"
        plotCorrelation -in {params.npz} --corMethod pearson --skipZeros \
        --plotTitle "Spearman Correlation of Read Counts" \
        --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
        -o {output.correlation} --outFileCorMatrix {output.mtrx}
        '''

rule macs2_qc_plot:
    input:
        expand(dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_peaks.narrowPeak",sample=samplelist)
    output:
        dir_out + "/macs2_qc_plot/macs2_peakqc.plots.pdf"
    benchmark:
        dir_out + "/log/macs2qc.bmk"
    params:
        peaks = lambda wildcards, input: ','.join(input),
        names = ','.join(samplelist),
        outdir = dir_out + "/macs2_qc_plot"
    shell:
        '''
        /home/jibzhang/bin/R_script/plot_macs2_qc.r -i {params.peaks} -s {params.names} -o {params.outdir}
        '''

rule homer_annotate_plot:
    input:
        expand(dir_out + "/{sample}/" + squence_type + "/Homer/{sample}.annotation.txt",sample=samplelist),
    output:
        dir_out + "/homer_qc_plot/homer_annotation.plots.pdf"
    benchmark:
        dir_out + "/log/homerqc.bmk"
    params:
        annos = lambda wildcards, input: ','.join(input),
        names = ','.join(samplelist),
        outdir = dir_out + "/homer_qc_plot"
    shell:
        '''
        /home/jibzhang/bin/R_script/plot_homer_annotatepeaks.r -i {params.annos} -s {params.names} -o {params.outdir}
        '''

rule done_all:
    input:
        featurecount = rules.featureCounts.benchmark,
        fingerprint = rules.plotfingerprint.benchmark,
        heatmap = rules.heatmap_consensus.benchmark,
        correaltion = rules.correlation.benchmark,
        macs2qc = rules.macs2_qc_plot.benchmark,
        homerqc = rules.homer_annotate_plot.benchmark
    output:
        dir_out + "/log/done_analysis.txt"
    shell:
        '''
        touch {output}
        '''