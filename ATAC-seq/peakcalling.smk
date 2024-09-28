include: "common_rule.smk"

rule MACS2:
    input:
        bam=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam",
        control_bam=lambda wildcards: "{dir_out}/{id}/{squence_type}/bam/{id}.bam".format(dir_out=dir_out,squence_type=squence_type,id=pair[wildcards.sample])
    output:
        peak=dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_peaks.narrowPeak",
        summit=dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_summits.bed"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.macs2.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.macs2.bmk"
    params:
        prefix="{sample}",
        para_macs2="-B -g 2805665311",
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

rule FRiP:
    input:
        bam = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam",
        flagstat = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.final.flagstats",
        peaks = rules.MACS2.output.peak
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
        READS_IN_PEAKS=$(intersectBed -a {input.bam} -b {input.peaks}| awk -F '\\t' '{{sum += $NF}} END {{print sum}}')
        grep 'mapped (' {input.flagstat} | grep -v "primary" | awk -v a="$READS_IN_PEAKS" -v OFS='\\t' '{{print "{params.prefix}", a/$1}}' > {output}
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
        peak_header = "/scratch/jibzhang/pipeline/ChIP_seq/peak_count_header.txt",
        frip_header = "/scratch/jibzhang/pipeline/ChIP_seq/frip_score_header.txt"
    shell:
        '''
        cat {input.peaks} | wc -l | awk -v OFS='\\t' '{{ print "{params.prefix}", $1 }}' | cat {params.peak_header} - > {output.peak_mqc}
        cat {params.frip_header} {input.FRiP} > {output.FRiP_mqc}
        '''

rule featureCounts:
    input:
        saf="/scratch/jibzhang/project/PRDM1_NK/merge_bed/nofeeder.saf",
        bam=lambda wildcards: "{dir_out}/{id}/{squence_type}/bam/{id}.bam".format(dir_out=dir_out,squence_type=squence_type,id=pair[wildcards.sample])
    output:
        tab=dir_out+"/{sample}/"+squence_type+"/featureCounts/{sample}.featureCounts",
    log:
        dir_out+"/{sample}/"+squence_type+"/log/{sample}_featureCounts.log"
    benchmark:
        dir_out+"/{sample}/"+squence_type+"/log/{sample}_featureCounts.bmk"
    threads:
        config["threads_featureCounts"]
    params:
        dir=dir_out+"/{sample}/"+squence_type+"/featureCounts/"
    shell:
        '''
        featureCounts -p -O -donotsort -F SAF --fracOverlap 0.2 --tmpDir {params.dir} \
        -a {input.saf} -o {output.tab} {input.bam} &> {log}
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

rule meme_chip:
    input:
        peak_bed=dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_peaks.blaclist_filtered.bed",
        summit_bed=dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_summits_range.bed"
    output:
        peak_motif = dir_out + "/{sample}/" + squence_type + "/memechip/fimo_out_1/fimo.gff",
        summit_motif = dir_out + "/{sample}/" + squence_type + "/memechip_summit/fimo_out_1/fimo.gff",
        peak_fasta = dir_out + "/{sample}/" + squence_type + "/memechip/{sample}_peak.fasta"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.memechip.bmk"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.memechip.log"
    params:
        ref=config["ref_gatk"],
        summit_fasta= dir_out + "/{sample}/" + squence_type + "/memechip/{sample}_summit.fasta",
        HOCO = "/home/jibzhang/motif_databases/HUMAN/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme",
        JASPAR = "/home/jibzhang/motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme",
        out_dir_peak=dir_out + "/{sample}/" + squence_type + "/memechip",
        out_dir_summit=dir_out + "/{sample}/" + squence_type + "/memechip_summit"
    shell:
        '''
        bedtools getfasta -fi {params.ref} -bed {input.peak_bed} > {output.peak_fasta}
        bedtools getfasta -fi {params.ref} -bed {input.summit_bed} > {params.summit_fasta}
        meme-chip -db {params.HOCO} -db {params.JASPAR} {output.peak_fasta} -oc {params.out_dir_peak} -meme-mod anr \
        -minw 6 -maxw 20 -dna -meme-nmotifs 10 -centrimo-score 5.0 -centrimo-ethresh 10.0 -centrimo-local \
        -streme-pvt 0.05 -spamo-skip &> {log}
        meme-chip -db {params.HOCO} -db {params.JASPAR} {params.summit_fasta} -oc {params.out_dir_summit} -meme-mod anr \
        -minw 6 -maxw 20 -dna -meme-nmotifs 10 -centrimo-score 5.0 -centrimo-ethresh 10.0 -centrimo-local \
        -streme-pvt 0.05 -spamo-skip &> {log}
        '''

rule meme_motif_summit:
    input:
        gff_top500 = dir_out + "/{sample}/" + squence_type + "/memechip_top500/fimo_out_1/fimo.gff",
        peak_top500 = dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_summits_top500.bed",
        gff_summit = dir_out + "/{sample}/" + squence_type + "/memechip_summit/fimo_out_1/fimo.gff",
        peak_summit = dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_summits_range.bed"
    output:
        motif_top500 = dir_out + "/{sample}/" + squence_type + "/memechip_top500/fimo_out_1/PRDM1_peaks.bed",
        motif = dir_out + "/{sample}/" + squence_type + "/memechip_summit/fimo_out_1/PRDM1_peaks.bed"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.mememotif_summit.bmk"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.mememotif_summit.log"
    params:
        bed_top500 = dir_out + "/{sample}/" + squence_type + "/memechip_top500/fimo_out_1/PRDM1_motif.bed",
        bed = dir_out + "/{sample}/" + squence_type + "/memechip_summit/fimo_out_1/PRDM1_motif.bed"
    shell:
        '''
        tail -n +2 {input.gff_top500} | awk -v OFS="\t" '{{print $1,$4,$5}}' > {params.bed_top500}
        bedtools intersect -wa -a {input.peak_top500} -b {params.bed_top500} | uniq > {output.motif_top500}
        tail -n +2 {input.gff_summit} | awk -v OFS="\t" '{{print $1,$4,$5}}' > {params.bed}
        bedtools intersect -wa -a {input.peak_summit} -b {params.bed} | uniq > {output.motif}
        '''

rule peak_summit:
    input:
        summit=rules.MACS2.output.summit
    output:
        fasta = dir_out + "/{sample}/" + squence_type + "/memechip/{sample}_summits.fasta",
        range = dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_summits_range.bed"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.range_macs2.bmk"
    params:
        ref=config["ref_gatk"],
        outdir=dir_out + "/{sample}/" + squence_type + "/memechip"
    shell:
        '''
        if [ ! -d "{params.outdir}" ]; then mkdir {params.outdir}; fi
        cat {input.summit} |awk -v OFS="\t" '{{print $1,$2-100,$3+100,$4,$5}}' > {output.range}
        bedtools getfasta -fi {params.ref} -bed {output.top} > {output.fasta}
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
        bed=dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_peaks.blaclist_filtered.bed",
        bmk=dir_out + "/{sample}/" + squence_type + "/log/{sample}.filter_bed.bmk"
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
        peak=rules.filter.output.filtered
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

rule heatmap_gene:
    input:
        bw=rules.bamcompare.output
    output:
        dir_out + "/{sample}/" + squence_type + "/heatmap/{sample}_gene.heatmap.png",
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}_gene.heatmap.bmk"
    params:
        ref=config["ref_gtf"],
        gene_matrix=dir_out + "/{sample}/" + squence_type + "/heatmap/{sample}_gene.matrix.gz",
    threads:
        8
    shell:
        '''
        computeMatrix scale-regions -p {threads} -S {input.bw} \
        -R {params.ref} \
        --skipZeros \
        --beforeRegionStartLength 3000 \
        --regionBodyLength 5000 \
        --afterRegionStartLength 3000 \
         -o {params.gene_matrix}
        plotHeatmap -m {params.gene_matrix} -out {output} --sortUsing sum
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
        --samplesLabel F26feeder_d13 F28feeder_d13 F28nofeeder_d6 F29feeder_d28 F29feeder_d31 F37M25nofeeder_d6 F23feeder_d13 F23nofeeder_d6 \
        -R {params.ref} -bl {params.blacklist} \
        --skipZeros \
        --beforeRegionStartLength 2500 \
        --regionBodyLength 5000 \
        --afterRegionStartLength 2500 \
         -o {params.gene_matrix}
        plotHeatmap -m {params.gene_matrix} -out {output.heatmap} --sortUsing sum
        '''

rule peak_merge:
    input:
        peak_all=expand(dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_peaks.blaclist_filtered.bed",sample=CHIP)
    output:
        dir_out+"/stats/peak_all.bed",
    shell:
        '''
        cat {input.peak_all}|sort -k1,1 -k2,2n |bedtools merge > {output}
        '''

rule heatmap_Allpeak:
    input:
        bw=rules.bamcompare.output,
        peak=rules.peak_merge.output
    output:
        dir_out + "/{sample}/" + squence_type + "/heatmap/{sample}.heatmap_Allpeak.png",
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.heatmap_Allpeak.bmk",
    params:
        gene_matrix=dir_out + "/{sample}/" + squence_type + "/heatmap/{sample}.heatmap_Allpeak.matrix.gz",
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

rule heatmap_Allpeak_merge:
    input:
        bw=expand(rules.bamcompare.output,sample=CHIP),
        peak=rules.peak_merge.output
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
        --samplesLabel F26feeder_d13 F28feeder_d13 F28nofeeder_d6 F29feeder_d28 F29feeder_d31 F37M25nofeeder_d6 F23feeder_d13 F23nofeeder_d6 M37_nofeeder_d6 M37_feeder_d10 M37_feeder_d13 M37_feeder_d28 \
        --referencePoint center --skipZeros -p {threads} -a 2000 -b 2000 \
        -o {params.gene_matrix}

        plotHeatmap -m {params.gene_matrix} -out {output} \
        --sortUsing sum --startLabel "Peak Start" --endLabel "Peak End" --xAxisLabel "" --regionsLabel "PeaksInAllSamples"
        '''

rule correlation:
    input:
        bw=expand(rules.bamcompare.output,sample=CHIP),
        peak=rules.peak_merge.output,
        heatmap=rules.heatmap_Allpeak_merge.output
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

rule meme:
    input:
        peak=rules.filter.output.filtered
    output:
        meme=dir_out + "/{sample}/" + squence_type + "/meme_top500/out/meme.txt",
        fasta=dir_out + "/{sample}/" + squence_type + "/meme_top500/{sample}.peaks.fasta"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.meme.log",
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.meme.bmk",
    params:
        ref=config["ref_gatk"],
        outdir=dir_out + "/{sample}/" + squence_type + "/meme_top500",
        out_dir=dir_out + "/{sample}/" + squence_type + "/meme_top500/out"
    threads:
        8
    shell:
        '''
        if [ ! -d "{params.outdir}" ]; then mkdir {params.outdir}; fi
        bedtools getfasta -fi {params.ref} -bed {input.peak} > {output.fasta}
        n_bed=$(cat {output.fasta}|wc -l)
        if [ $n_bed -gt 1 ] ; then
        meme -p {threads} -mod anr -nmotifs 10 -minw 6 -maxw 20 -searchsize 100000 -time 2919 -oc {params.out_dir} {output.fasta} &> {log}
        else
        touch {output.meme}
        fi
        '''

rule HOMER:
    input:
        bed = dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_peaks.blaclist_filtered.bed",
        fasta = config["ref_bwa"],
        gtf = config["ref_gtf"]
    output:
        anno=dir_out + "/{sample}/" + squence_type + "/Homer/{sample}.annotation.txt",
        motif=dir_out + "/{sample}/" + squence_type + "/Homer/motif/homerResults/motif1.motif"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.homer.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.homer.bmk"
    params:
        dir=dir_out + "/{sample}/" + squence_type + "/Homer/GO",
        motifdir=dir_out + "/{sample}/" + squence_type + "/Homer/motif",
        KORNAseq = "/scratch/jibzhang/project/PRDM1_NK/PRDM1KO_vs_WT.txt",
        overexpr = "/scratch/jibzhang/project/PRDM1_NK/PRDM1_overexpres_vs_low.txt"
    shell:
        '''
        annotatePeaks.pl {input.bed} {input.fasta} -gtf {input.gtf} -go {params.dir} -gene {params.KORNAseq} {params.overexpr} -cpu 4 > {output.anno}
        findMotifsGenome.pl {input.bed} hg38 {params.motifdir} -size 50 &> {log}
        '''

rule Homer_motif:
    input:
        fasta = dir_out + "/{sample}/" + squence_type + "/memechip/{sample}_peak.fasta",
        motif = dir_out + "/{sample}/" + squence_type + "/Homer/motif/homerResults/motif1.motif",
        peak = dir_out + "/{sample}/" + squence_type + "/macs2/{sample}_peaks.blaclist_filtered.bed"
    output:
        dir_out + "/{sample}/" + squence_type + "/Homer/motif/motif1_peaks.bed"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.homer_motif.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.homer_motif.bmk"
    params:
        bed = dir_out + "/{sample}/" + squence_type + "/Homer/motif/motif1.bed"
    shell:
        '''
        scanMotifGenomeWide.pl {input.motif} {input.fasta} > {params.bed}
        bedtools intersect -wa -a {input.peak} -b {params.bed} | uniq > {output}
        '''

rule counts:
    input:
        dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam"
    output:
        dir_out + "/{sample}/" + squence_type + "/bam/{sample}.counts.txt"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.counts.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.counts.bmk"
    params:
        peaks = "/scratch/jibzhang/project/PRDM1_NK/merged_peaks.bed",
        bed = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bed",
        sorted_bam = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.name_sorted.bam"
    shell:
        '''
        samtools sort -n -o {params.sorted_bam} -O bam {input}
        bedtools bamtobed -i {params.sorted_bam} -bedpe > {params.bed}
        bedtools coverage -a {params.peaks} -b {params.bed} -counts > {output}
        '''

rule done:
    input:
        multiqc = rules.multiqc_peaks.benchmark,
        bdgcmp = rules.bamcompare.benchmark,
        neargene = rules.nearbygene.benchmark,
        heatmap_peak = rules.heatmap_peak.benchmark,
        #heatmap_merge = rules.heatmap_Allpeak_merge.benchmark,
        #heatmap_Allgene = rules.heatmap_gene_allSample.benchmark,
        #correlation = rules.correlation.benchmark,
        meme = rules.meme.benchmark,
        mememotif = rules.meme_motif.benchmark,
        mememotif_summit = rules.meme_motif_summit.benchmark,
        homer_motif = rules.Homer_motif.benchmark,
        #counts=rules.counts.benchmark
    output:
        dir_out + "/{sample}/" + squence_type + "/log/done_analysis.txt"
    shell:
        '''
        touch {output}
        '''
