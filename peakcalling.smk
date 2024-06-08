include: "common_rule.smk"

rule MACS2_cutoff:
    input:
        bam=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam",
        control_bam=lambda wildcards: "{dir_out}/{id}/{squence_type}/bam/{id}.bam".format(dir_out=dir_out,squence_type=squence_type,id= pair[wildcards.sample])
    output:
        peak = dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_cutoff_analysis.txt"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.macs2_cutoff.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.macs2_cutoff.bmk"
    params:
        prefix="{sample}",
        para_macs2="-g 2805665311 -cutoff-analysis",
        dir=dir_out + "/{sample}/" + squence_type + "/macs2",
        fq_n=lambda wildcards: len(fq[wildcards.sample])
    shell:
        '''
        if [ {params.fq_n} -gt 1 ] ; then
        macs2 callpeak {params.para_macs2} -f BAMPE -t {input.bam} -c {input.control_bam} --cutoff-analysis --min-length 1000 --nomodel --extsize 150 -n {params.prefix} --outdir {params.dir} 2> {log}
        else
        macs2 callpeak {params.para_macs2} -f BAM -t {input.bam} -c {input.control_bam} --cutoff-analysis --min-length 1000 --nomodel --extsize 150 -n {params.prefix} --outdir {params.dir} 2> {log}
        fi
        '''

rule MACS2:
    input:
        bam = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam",
        control_bam = lambda wildcards: "{dir_out}/{id}/{squence_type}/bam/{id}.bam".format(dir_out=dir_out,squence_type=squence_type,id=pair[wildcards.sample])
    output:
        peak = dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_peaks.narrowPeak",
        summit = dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_summits.bed"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.macs2.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.macs2.bmk"
    params:
        prefix="{sample}",
        para_macs2="-g 2805665311 --call-summits",
        dir=dir_out + "/{sample}/" + squence_type + "/macs2",
        fq_n= lambda wildcards:len(fq[wildcards.sample])
    shell:
        '''
        if [ {params.fq_n} -gt 1 ] ; then
        macs2 callpeak {params.para_macs2} -f BAMPE -t {input.bam} -c {input.control_bam} -n {params.prefix} --outdir {params.dir} 2> {log}
        else
        macs2 callpeak {params.para_macs2} -f BAM -t {input.bam} -c {input.control_bam} -n {params.prefix} --outdir {params.dir} 2> {log}
        fi
        '''

rule peak_summit:
    input:
        summit=rules.MACS2.output.summit
    output:
        range = dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_summits_range.bed"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.range_macs2.bmk"
    shell:
        '''
        cat {input.summit} |awk -v OFS="\t" '{{print $1,$2-100,$3+100,$4,$5}}' > {output.range}
        '''

rule plotfingerprint:
    input:
        bam = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam",
        control_bam = lambda wildcards: "{dir_out}/{id}/{squence_type}/bam/{id}.bam".format(dir_out=dir_out,squence_type=squence_type,id=pair[wildcards.sample])
    output:
        plot = dir_out + "/{sample}/" + squence_type + "/deeptools/plotFingerprint.pdf",
        rawcount = dir_out + "/{sample}/" + squence_type + "/deeptools/plotFingerprint.raw.txt",
        qcmetrics = dir_out + "/{sample}/" + squence_type + "/deeptools/plotFingerprint.qcmetrics.txt"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.plotfingerprint.bmk"
    threads:
        6
    shell:
        '''
        plotFingerprint --bamfiles {input.bam} {input.control_bam} --skipZeros -o {output.plot} --outRawCounts {output.rawcount} --outQualityMetrics {output.qcmetrics} -p {threads}
        '''

rule featureCounts:
    input:
        annotation = dir_out + "/consensus_peak/feeder_nofeeder_consensus.saf",
        bam = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam",
        control_bam = lambda wildcards: "{dir_out}/{id}/{squence_type}/bam/{id}.bam".format(dir_out=dir_out,squence_type=squence_type,id=pair[wildcards.sample])
    output:
        dir_out + "/{sample}/" + squence_type + "/{sample}.featureCounts.txt"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.featurecounts.bmk"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.featurecounts.log"
    threads:
        6
    shell:
        '''
        featureCounts -F SAF -p -O --fracOverlap 0.2 -T {threads} -a {input.annotation} -o {output} {input.bam} {input.control_bam} &> {log}
        '''

rule counts:
    input:
        bam = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam",
        peaks = dir_out + "/consensus_peak/feeder_nofeeder_consensus.bed"
    output:
        dir_out + "/{sample}/" + squence_type + "/bam/{sample}.counts.txt"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.counts.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.counts.bmk"
    params:
        bed = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bed",
        sorted_bam = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.name_sorted.bam"
    shell:
        '''
        samtools sort -n -o {params.sorted_bam} -O bam {input.bam}
        bedtools bamtobed -i {params.sorted_bam} > {params.bed}
        bedtools coverage -a {input.peaks} -b {params.bed} -counts > {output}
        '''

rule FRiP:
    input:
        flagstat = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.final.flagstats",
        counts = rules.counts.output
    output:
        dir_out + "/{sample}/" + squence_type + "/macs2/{sample}.FRiP.txt"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.FRiP.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.FRiP.bmk"
    params:
        prefix="{sample}"
    shell:
        '''
        READS_IN_PEAKS=$(cat {input.counts} | awk -F '\\t' '{{sum += $NF}} END {{print sum}}')
        grep 'mapped (' {input.flagstat} | awk -v a="$READS_IN_PEAKS" -v OFS='\\t' '{{print "{params.prefix}", a/$1}}' > {output}
        '''

rule multiqc_peaks:
    input:
        peaks = rules.MACS2.output.peak,
        FRiP = rules.FRiP.output
    output:
        peak_mqc = dir_out + "/{sample}/" + squence_type + "/macs2/{sample}.peak_count_mqc.tsv",
        FRiP_mqc = dir_out + "/{sample}/" + squence_type + "/macs2/{sample}.FRiP_mqc.tsv"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.multiqc_peaks.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.multiqc_peaks.bmk"
    params:
        prefix = "{sample}",
        peak_header = "/coh_labs/jochan/pipeline/ChIP_seq/peak_count_header.txt",
        frip_header = "/coh_labs/jochan/pipeline/ChIP_seq/frip_score_header.txt"
    shell:
        '''
        cat {input.peaks} | wc -l | awk -v OFS='\\t' '{{ print "{params.prefix}", $1 }}' | cat {params.peak_header} - > {output.peak_mqc}
        cat {params.frip_header} {input.FRiP} > {output.FRiP_mqc}
        '''

rule filter:
    input:
        peaks = rules.MACS2.output.peak,
        blacklist = "/home/jibzhang/reference/hg38_region/hg38-ChIPblacklist.v3.bed",
        summits = rules.peak_summit.output.range
    output:
        filter_peak = dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_filtered_peaks.bed",
        blacklist_filtered = dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_peaks.blaclist_filtered.bed",
        filter_summit = dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_filtered_summits.bed",
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.filter_bed.bmk"
    shell:
        '''
        bedtools intersect -v -a {input.peaks} -b {input.blacklist} > {output.blacklist_filtered}
        awk -F "\\t" '{{if(($7 > 4) && ($5 > 110)) {{print}} }}' {output.blacklist_filtered} > {output.filter_peak}
        bedtools intersect -wa -a {input.summits} -b {output.filter_peak} > {output.filter_summit}
        '''

rule meme_chip:
    input:
        peak_bed=rules.filter.output.filter_peak,
        summit_bed=rules.filter.output.filter_summit
    output:
        peak_motif = dir_out + "/{sample}/" + squence_type + "/memechip/fimo_out_1/fimo.gff",
        summit_motif = dir_out + "/{sample}/" + squence_type + "/memechip_summit/fimo_out_1/fimo.gff",
        peak_fasta = dir_out + "/{sample}/" + squence_type + "/memechip/{sample}_peak.fasta",
        summit_fasta= dir_out + "/{sample}/" + squence_type + "/memechip/{sample}_summit.fasta"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.memechip.bmk"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.memechip.log"
    params:
        ref=config["ref_gatk"],
        HOCO = "/home/jibzhang/motif_databases/HUMAN/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme",
        JASPAR = "/home/jibzhang/motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme",
        out_dir_peak=dir_out + "/{sample}/" + squence_type + "/memechip",
        out_dir_summit=dir_out + "/{sample}/" + squence_type + "/memechip_summit"
    shell:
        '''
        bedtools getfasta -fi {params.ref} -bed {input.peak_bed} > {output.peak_fasta}
        bedtools getfasta -fi {params.ref} -bed {input.summit_bed} > {output.summit_fasta}
        meme-chip -db {params.HOCO} -db {params.JASPAR} {output.peak_fasta} -oc {params.out_dir_peak} -meme-mod anr \
        -minw 6 -maxw 20 -dna -meme-nmotifs 10 -centrimo-score 5.0 -centrimo-ethresh 10.0 -centrimo-local \
        -streme-pvt 0.05 -spamo-skip &> {log}
        meme-chip -db {params.HOCO} -db {params.JASPAR} {output.summit_fasta} -oc {params.out_dir_summit} -meme-mod anr \
        -minw 6 -maxw 20 -dna -meme-nmotifs 10 -centrimo-score 5.0 -centrimo-ethresh 10.0 -centrimo-local \
        -streme-pvt 0.05 -spamo-skip &> {log}
        '''

rule meme_motif_summit:
    input:
        gff_summit = dir_out + "/{sample}/" + squence_type + "/memechip_summit/fimo_out_1/fimo.gff",
        peak_summit = rules.filter.output.filter_summit
    output:
        motif = dir_out + "/{sample}/" + squence_type + "/memechip_summit/fimo_out_1/PRDM1_peaks.bed"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.mememotif_summit.bmk"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.mememotif_summit.log"
    params:
        bed = dir_out + "/{sample}/" + squence_type + "/memechip_summit/fimo_out_1/PRDM1_motif.bed"
    shell:
        '''
        tail -n +2 {input.gff_summit} | awk -v OFS="\t" '{{print $1,$4,$5}}' > {params.bed}
        bedtools intersect -wa -a {input.peak_summit} -b {params.bed} | uniq > {output.motif}
        '''

rule meme_motif:
    input:
        gff = dir_out + "/{sample}/" + squence_type + "/memechip/fimo_out_1/fimo.gff",
        peak = rules.filter.output.filter_peak
    output:
        motif = dir_out + "/{sample}/" + squence_type + "/memechip/fimo_out_1/PRDM1_peaks.bed"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.mememotif.bmk"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.mememotif.log"
    params:
        bed = dir_out + "/{sample}/" + squence_type + "/memechip/fimo_out_1/PRDM1_motif.bed"
    shell:
        '''
        tail -n +2 {input.gff} | awk -v OFS="\t" '{{print $1,$4,$5}}' > {params.bed}
        bedtools intersect -wa -a {input.peak} -b {params.bed} | uniq > {output.motif}
        '''

rule bdgcmp:
    input:
        case=dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_treat_pileup.bdg",
        control=dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_control_lambda.bdg"
    output:
        dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_bdgcmp.bdg"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.macs2_bdgcmp.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.macs2_bdgcmp.bmk"
    shell:
        '''
        macs2 bdgcmp -t {input.case} -c {input.control} -m ppois -o {output} 2> {log}
        '''

rule bamcompare:
    input:
        bam=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam",
        control_bam=lambda wildcards: "{dir_out}/{id}/{squence_type}/bam/{id}.bam".format(dir_out=dir_out,squence_type=squence_type,id=pair[wildcards.sample])
    output:
        dir_out + "/{sample}/" + squence_type + "/deeptools/{sample}.bamcompare.bw"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.bamcompare.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.bamcompare.bmk"
    params:
        para_bamcompare="--binSize 50 --scaleFactorsMethod None --normalizeUsing BPM --smoothLength 200 --centerReads -p 6"
    shell:
        '''
        bamCompare -b1 {input.bam} -b2 {input.control_bam} -o {output} {params.para_bamcompare} 2> {log}
        '''

rule nearbygene:
    input:
        bed=rules.filter.output.filter_peak
    output:
        dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_peaks.stepwiseAnno.bed"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.nearbygene.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.nearbygene.bmk"
    params:
        promoter="/home/jibzhang/reference/hg38_region/promoter5k_hg38.bed",
        upstream="/home/jibzhang/reference/hg38_region/upstream50k_hg38.bed",
        downstream="/home/jibzhang/reference/hg38_region/downstream50k_hg38.bed",
        genebody="/home/jibzhang/reference/hg38_region/genebody2500plus_hg38.bed",
        promoter_overlap= dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_peaks.promoter5k.bed",
        genebody_overlap= dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_peaks.genebody2500plus.bed",
        downstream_overlap= dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_peaks.downstream50k.bed",
        upstream_overlap= dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_peaks.upstream50k.bed",
        prefix=dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_peaks",
        stepwise=dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_peaks.stepwiseAnno.bed.1"
    shell:
        '''
        bedtools intersect -loj -a {input.bed} -b {params.promoter} > {params.promoter_overlap}
        bedtools intersect -loj -a {input.bed} -b {params.genebody} > {params.genebody_overlap}
        bedtools intersect -loj -a {input.bed} -b {params.downstream} > {params.downstream_overlap}
        bedtools intersect -loj -a {input.bed} -b {params.upstream} > {params.upstream_overlap}
        /home/jibzhang/bin/combine.annoPeaks.pl {params.prefix} > {output}
        cat {output} | sort -k1,1 -k2,2n > {params.stepwise}
        mv {params.stepwise} {output}
        '''

rule heatmap_peak:
    input:
        bw=rules.bamcompare.output,
        peak=rules.filter.output.filter_peak
    output:
        dir_out + "/{sample}/" + squence_type + "/heatmap/{sample}_peak.heatmap.png"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}_peak.heatmap.bmk",
    params:
        gene_matrix=dir_out + "/{sample}/" + squence_type + "/heatmap/{sample}_peak.matrix.gz",
    threads:
        8
    shell:
        '''
        computeMatrix reference-point -S {input.bw} -R {input.peak} \
        --referencePoint center --skipZeros -p {threads} -a 1000 -b 1000 \
        -o {params.gene_matrix}

        plotHeatmap -m {params.gene_matrix} -out {output} \
        --sortUsing sum --startLabel "Peak Start" --endLabel "Peak End" --xAxisLabel "" --regionsLabel "Peaks"
        '''

rule HOMER:
    input:
        bed = rules.filter.output.filter_peak,
        fasta = config["ref_bwa"],
        gtf = config["ref_gtf"]
    output:
        anno = dir_out + "/{sample}/" + squence_type + "/Homer/{sample}.annotation.txt"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.homer.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.homer.bmk"
    params:
        motifdir = dir_out + "/{sample}/" + squence_type + "/Homer/motif"
    shell:
        '''
        annotatePeaks.pl {input.bed} {input.fasta} -gtf {input.gtf} -cpu 4 > {output.anno}
        findMotifsGenome.pl {input.bed} hg38 {params.motifdir} -size 200 &> {log}
        '''

rule Homer_motif:
    input:
        fasta = rules.meme_chip.output.summit_fasta,
        PRDM1IRF_motif = dir_out + "/consensus_peak/feeder_d13/Homer/motif/knownResults/PRDM1_IRF.motif",
        PRDM1_motif = dir_out + "/consensus_peak/feeder_d13/Homer/motif/knownResults/known1.motif",
        RUNX1_motif = dir_out + "/consensus_peak/feeder_d13/Homer/motif/knownResults/known9.motif",
        RUNX_motif = dir_out + "/consensus_peak/feeder_d13/Homer/motif/knownResults/RUNX.motif",
        peak = rules.filter.output.filter_summit
    output:
        PRDM1_peak = dir_out + "/{sample}/" + squence_type + "/Homer/motif/PRDM1_peaks.bed",
        RUNX1_peak = dir_out + "/{sample}/" + squence_type + "/Homer/motif/RUNX1_peaks.bed",
        PRDM1IRF_peak = dir_out + "/{sample}/" + squence_type + "/Homer/motif/PRDM1IRF_peaks.bed",
        RUNX_peak = dir_out + "/{sample}/" + squence_type + "/Homer/motif/RUNX_peaks.bed"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.homer.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.homer_motif.bmk"
    params:
        PRDM1_bed = dir_out + "/{sample}/" + squence_type + "/Homer/motif/PRDM1.bed",
        RUNX1_bed = dir_out + "/{sample}/" + squence_type + "/Homer/motif/RUNX1.bed",
        RUNX_bed = dir_out + "/{sample}/" + squence_type + "/Homer/motif/RUNX.bed",
        PRDM1IRF_bed = dir_out + "/{sample}/" + squence_type + "/Homer/motif/PRDM1IRF.bed",
        PRDM1IRF_region = dir_out + "/{sample}/" + squence_type + "/Homer/motif/PRDM1IRF_region.bed",
        PRDM1_region = dir_out + "/{sample}/" + squence_type + "/Homer/motif/PRDM1_region.bed",
        RUNX1_region = dir_out + "/{sample}/" + squence_type + "/Homer/motif/RUNX1_region.bed",
        RUNX_region = dir_out + "/{sample}/" + squence_type + "/Homer/motif/RUNX_region.bed"
    shell:
        '''
        scanMotifGenomeWide.pl {input.PRDM1IRF_motif} {input.fasta} > {params.PRDM1IRF_bed} &> {log}
        cat {params.PRDM1IRF_bed} | awk '{{split($2,a,"[:-]"); printf("%s\\t%s\\t%s\\n",a[1],a[2],a[3]);}}' > {params.PRDM1IRF_region}
        bedtools intersect -wa -a {input.peak} -b {params.PRDM1IRF_region} | uniq > {output.PRDM1IRF_peak}
        scanMotifGenomeWide.pl {input.PRDM1_motif} {input.fasta} > {params.PRDM1_bed} >> {log}
        cat {params.PRDM1_bed} | awk '{{split($2,a,"[:-]"); printf("%s\\t%s\\t%s\\n",a[1],a[2],a[3]);}}' > {params.PRDM1_region}
        bedtools intersect -wa -a {input.peak} -b {params.PRDM1_region} | uniq > {output.PRDM1_peak}
        scanMotifGenomeWide.pl {input.RUNX1_motif} {input.fasta} > {params.RUNX1_bed} >> {log}
        cat {params.RUNX1_bed} | awk '{{split($2,a,"[:-]"); printf("%s\\t%s\\t%s\\n",a[1],a[2],a[3]);}}' > {params.RUNX1_region}
        bedtools intersect -wa -a {input.peak} -b {params.RUNX1_region} | uniq > {output.RUNX1_peak}
        scanMotifGenomeWide.pl {input.RUNX_motif} {input.fasta} > {params.RUNX_bed} >> {log}
        cat {params.RUNX_bed} | awk '{{split($2,a,"[:-]"); printf("%s\\t%s\\t%s\\n",a[1],a[2],a[3]);}}' > {params.RUNX_region}
        bedtools intersect -wa -a {input.peak} -b {params.RUNX_region} | uniq > {output.RUNX_peak}
        '''

rule heatmap_Allpeak:
    input:
        bw=rules.bamcompare.output,
        peak=dir_out + "/consensus_peak/consensus.bed",
    output:
        dir_out + "/{sample}/" + squence_type + "/heatmap/{sample}.heatmap_Allpeak.png"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.heatmap_Allpeak.bmk"
    params:
        gene_matrix=dir_out + "/{sample}/" + squence_type + "/heatmap/{sample}.heatmap_Allpeak.matrix.gz"
    threads:
        8
    shell:
        '''
        computeMatrix reference-point -S {input.bw} -R {input.peak} \
        --referencePoint center --skipZeros -p {threads} -a 3000 -b 3000 \
        -o {params.gene_matrix}

        plotHeatmap -m {params.gene_matrix} -out {output} \
        --sortUsing sum --startLabel "Peak Start" --endLabel "Peak End" --xAxisLabel "" --regionsLabel "PeaksInAllSamples"
        '''

rule done:
    input:
        multiqc = rules.multiqc_peaks.benchmark,
        bdgcmp = rules.bamcompare.benchmark,
        neargene = rules.nearbygene.benchmark,
        fingerprint = rules.plotfingerprint.benchmark,
        filter = rules.filter.benchmark,
        mememotif = rules.meme_motif.benchmark,
        mememotif_summit = rules.meme_motif_summit.benchmark,
        homer_motif = rules.Homer_motif.benchmark,
        featurecount = rules.featureCounts.benchmark,
        counts = rules.counts.benchmark
    output:
        dir_out + "/{sample}/" + squence_type + "/log/done_analysis.txt"
    shell:
        '''
        touch {output}
        '''

rule macs2_qc_plot:
    input:
        expand(rules.MACS2.output.peak,sample=CHIP),
    output:
        dir_out + "/macs2_qc_plot/macs2_peakqc.plots.pdf"
    benchmark:
        dir_out + "/log/macs2qc.bmk"
    params:
        peaks = lambda wildcards, input: ','.join(input),
        names = ','.join(CHIP),
        outdir = dir_out + "/macs2_qc_plot"
    shell:
        '''
        /home/jibzhang/bin/R_script/plot_macs2_qc.r -i {params.peaks} -s {params.names} -o {params.outdir}
        '''

rule homer_annotate_plot:
    input:
        expand(rules.HOMER.output.anno,sample=CHIP),
    output:
        dir_out + "/homer_qc_plot/homer_annotation.plots.pdf"
    benchmark:
        dir_out + "/log/homerqc.bmk"
    params:
        annos = lambda wildcards, input: ','.join(input),
        names = ','.join(CHIP),
        outdir = dir_out + "/homer_qc_plot"
    shell:
        '''
        /home/jibzhang/bin/R_script/plot_homer_annotatepeaks.r -i {params.annos} -s {params.names} -o {params.outdir}
        '''

rule heatmap_Allpeak_merge:
    input:
        bw = expand(rules.bamcompare.output,sample=CHIP),
        peak = dir_out + "/consensus_peak/consensus.bed"
    output:
        dir_out+"/stats/heatmap/heatmap_Allpeak_merge.png",
    benchmark:
        dir_out+"/stats/heatmap/heatmap_Allpeak_merge.bmk",
    params:
        gene_matrix=dir_out+"/stats/heatmap/heatmap_Allpeak_merge.matrix.gz",
    threads:
        16
    shell:
        '''
        computeMatrix reference-point -S {input.bw} -R {input.peak} \
        --samplesLabel F26feeder_d13 F28feeder_d13 F28nofeeder_d6 F29feeder_d28 F29feeder_d31 F23feeder_d13 F23nofeeder_d6 M37_nofeeder_d6 M37_feeder_d10 M37_feeder_d13 M37_feeder_d28 \
        --referencePoint center --skipZeros -p {threads} -a 2000 -b 2000 \
        -o {params.gene_matrix}

        plotHeatmap -m {params.gene_matrix} -out {output} \
        --sortUsing sum --startLabel "Peak Start" --endLabel "Peak End" --xAxisLabel "" --regionsLabel "PeaksInAllSamples"
        '''

rule heatmap_gene_allSample:
    input:
        bw=expand(rules.bamcompare.output,sample=CHIP)
    output:
        heatmap=dir_out+"/stats/heatmap/heatmap_gene_allSample.png",
    benchmark:
        dir_out+"/stats/heatmap/heatmap_gene_allSample.bmk",
    params:
        ref=config["ref_gtf"],
        gene_matrix=dir_out+"/stats/heatmap/heatmap_gene_allSample.matrix.gz",
        blacklist="/home/jibzhang/reference/hg38_region/hg38-ChIPblacklist.v3.bed"
    threads:
        16
    shell:
        '''
        computeMatrix scale-regions -p {threads} -S {input.bw} \
        --samplesLabel F26feeder_d13 F28feeder_d13 F28nofeeder_d6 F29feeder_d28 F29feeder_d31 F23feeder_d13 F23nofeeder_d6 \
        -R {params.ref} -bl {params.blacklist} \
        --skipZeros \
        --beforeRegionStartLength 2500 \
        --regionBodyLength 5000 \
        --afterRegionStartLength 2500 \
         -o {params.gene_matrix}
        plotHeatmap -m {params.gene_matrix} -out {output.heatmap} --sortUsing sum
        '''

rule correlation:
    input:
        bw = expand(rules.bamcompare.output,sample=CHIP),
        peak = dir_out + "/consensus_peak/consensus.bed",
        heatmap = rules.heatmap_Allpeak_merge.output
    output:
        pca=dir_out+"/stats/correlation/PCA_readcounts.png",
        correlation = dir_out+"/stats/correlation/Pearson_corr_heatmap.png",
        mtrx= dir_out+"/stats/correlation/Pearson_corr_readCounts.tab"
    benchmark:
        dir_out+"/stats/correlation/PCA_correlation.bmk"
    params:
        npz=dir_out+"/stats/correlation/results.npz",
    shell:
        '''
        multiBigwigSummary BED-file -b {input.bw} -o {params.npz} --BED {input.peak}
        plotPCA -in {params.npz} -o {output.pca} -T "PCA of read counts"
        plotCorrelation -in {params.npz} --corMethod pearson --skipZeros \
        --plotTitle "Spearman Correlation of Read Counts" \
        --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
        -o {output.correlation} --outFileCorMatrix {output.mtrx}
        '''

rule done_all:
    input:
        heatmap_merge = rules.heatmap_Allpeak_merge.benchmark,
        heatmap_Allgene = rules.heatmap_gene_allSample.benchmark,
        correlation = rules.correlation.benchmark,
        macs2qc = rules.macs2_qc_plot.benchmark,
        homerqc = rules.homer_annotate_plot.benchmark
    output:
        dir_out + "/log/done_all.txt"
    shell:
        '''
        touch {output}
        '''