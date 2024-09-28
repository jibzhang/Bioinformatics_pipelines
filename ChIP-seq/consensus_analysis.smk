rule consensus_peak:
    input:
        lambda wildcards: [dir_out + "/" + rep + "/" + squence_type + "/macs2/" + rep +"_filtered_peaks.bed" for rep in comparison[wildcards.group]]
    output:
        bed = dir_out + "/consensus_peak/{group}/{group}_consensus.bed",
        saf = dir_out + "/consensus_peak/{group}/{group}_consensus.saf",
        boolean = dir_out + "/consensus_peak/{group}/{group}_boolean.txt"
    params:
        merged = dir_out + "/consensus_peak/{group}/{group}_merged.txt",
        intersect = dir_out + "/consensus_peak/{group}/{group}_boolean.intersect.txt",
        plot = dir_out + "/consensus_peak/{group}/{group}_boolean.intersect.plot.pdf",
        names= lambda wildcards,input: ','.join(input.split("/")[4] for input in input)
    benchmark:
        dir_out + "/log/{group}_consensus_peak.bmk"
    shell:
        '''
        sort -k1,1 -k2,2n {input} | mergeBed -c 2,3,4,5,6,7,8,9,10 -o collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse >{params.merged}
        /home/jibzhang/bin/macs2_merged_expand.py {params.merged} {params.names} {output.boolean} --min_replicates 2 --is_narrow_peak
        awk -v FS='\\t' -v OFS='\\t' 'FNR > 1 {{ print $1, $2, $3, $4, "0", "+" }}' {output.boolean} > {output.bed}
        echo -e "GeneID\\tChr\\tStart\\tEnd\\tStrand" > {output.saf}
        awk -v FS='\\t' -v OFS='\\t' 'FNR > 1 {{ print $4, $1, $2, $3, "+" }}' {output.boolean} >> {output.saf}
        /home/jibzhang/bin/R_script/plot_peak_intersect.r -i {params.intersect} -o {params.plot}
        '''

rule consensus_promoter:
    input:
        lambda wildcards: [dir_out + "/" + rep + "/" + squence_type + "/macs2/" + rep +"_peaks.promoter5k.bed" for rep in comparison[wildcards.group]]
    output:
        bed = dir_out + "/consensus_promoter/{group}/{group}_consensus.bed",
        saf = dir_out + "/consensus_promoter/{group}/{group}_consensus.saf",
        boolean = dir_out + "/consensus_promoter/{group}/{group}_boolean.txt"
    params:
        merged = dir_out + "/consensus_promoter/{group}/{group}_merged.txt",
        intersect = dir_out + "/consensus_promoter/{group}/{group}_boolean.intersect.txt",
        plot = dir_out + "/consensus_promoter/{group}/{group}_boolean.intersect.plot.pdf",
        names= lambda wildcards,input: ','.join(input.split("/")[4] for input in input)
    benchmark:
        dir_out + "/log/{group}_consensus_promoter.bmk"
    shell:
        '''

        echo -e "GeneID\\tChr\\tStart\\tEnd\\tStrand" > {output.saf}
        awk -v FS='\\t' -v OFS='\\t' 'FNR > 1 {{ print $4, $1, $2, $3, "+" }}' {output.boolean} >> {output.saf}
        /home/jibzhang/bin/R_script/plot_peak_intersect.r -i {params.intersect} -o {params.plot}
        '''

rule consensus_summit:
    input:
        lambda wildcards: [dir_out + "/" + rep + "/" + squence_type + "/macs2/" + rep +"_filtered_summits.bed" for rep in comparison[wildcards.group]]
    output:
        bed = dir_out + "/consensus_peak/{group}/{group}_consensus.summits",
        summit = dir_out + "/consensus_peak/{group}/{group}_consensus.summits.bed"
    params:
        consensus_peak = dir_out + "/consensus_peak/{group}/{group}_consensus.bed"
    benchmark:
        dir_out + "/log/{group}_consensus_summits.bmk"
    shell:
        '''
        cat {input} | bedtools sort | bedtools merge > {output.bed}
        bedtools intersect -wa -a {output.bed} -b {params.consensus_peak} > {output.summit}
        '''

rule consensus_annotate:
    input:
        bed = rules.consensus_peak.output.bed
    output:
        annotation = dir_out + "/consensus_peak/{group}/{group}_consensus.annotatePeaks.txt",
    benchmark:
        dir_out + "/log/{group}.consensus_annotate.bmk"
    log:
        dir_out + "/log/{group}.consensus_annotate.log"
    params:
        dir = dir_out + "/consensus_peak/{group}/Homer/GO",
        motifdir = dir_out + "/consensus_peak/{group}/Homer/motif",
        KORNAseq = "/coh_labs/jochan/PRDM1_NK/PRDM1KO_vs_WT.txt",
        overexpr6h = "/coh_labs/jochan/PRDM1_NK/PRDM1_overexpres_vs_control_6h.txt",
	    overexpr12h = "/coh_labs/jochan/PRDM1_NK/PRDM1_overexpres_vs_control_12h.txt",
	    feedervsIL2 = "/coh_labs/jochan/PRDM1_NK/D13feeder_vs_D6nofeeder_Deseq2_DEGs.txt"
    shell:
        '''
        annotatePeaks.pl {input.bed} hg38 -go {params.dir} -gene {params.KORNAseq} {params.overexpr6h} {params.overexpr12h} {params.feedervsIL2} -cpu 4 > {output.annotation} &> {log}
        findMotifsGenome.pl {input.bed} hg38 {params.motifdir} -size 200 >> {log}
        '''

rule promoter_annotate:
    input:
        bed = rules.consensus_peak.output.bed
    output:
        annotation = dir_out + "/consensus_promoter/{group}/{group}_consensus.annotatePeaks.txt",
    benchmark:
        dir_out + "/log/{group}.consensus_annotate.bmk"
    log:
        dir_out + "/log/{group}.consensus_annotate.log"
    params:
        dir = dir_out + "/consensus_promoter/{group}/Homer/GO",
        motifdir = dir_out + "/consensus_promoter/{group}/Homer/motif",
        KORNAseq = "/coh_labs/jochan/PRDM1_NK/PRDM1KO_vs_WT.txt",
        overexpr6h = "/coh_labs/jochan/PRDM1_NK/PRDM1_overexpres_vs_control_6h.txt",
	    overexpr12h = "/coh_labs/jochan/PRDM1_NK/PRDM1_overexpres_vs_control_12h.txt",
	    feedervsIL2 = "/coh_labs/jochan/PRDM1_NK/D13feeder_vs_D6nofeeder_Deseq2_DEGs.txt"
    shell:
        '''
        annotatePeaks.pl {input.bed} hg38 -go {params.dir} -gene {params.KORNAseq} {params.overexpr6h} {params.overexpr12h} {params.feedervsIL2} -cpu 4 > {output.annotation} &> {log}
        findMotifsGenome.pl {input.bed} hg38 {params.motifdir} -size 200 >> {log}
        '''

rule annotate_boolean:
    input:
        bed = rules.consensus_annotate.output,
        boolean = rules.consensus_peak.output.boolean
    output:
        dir_out + "/consensus_peak/{group}_boolean.annotatePeaks.txt"
    benchmark:
        dir_out + "/log/{group}_consensus_boolean.bmk"
    params:
        trimmed = dir_out + "/consensus_peak/{group}_consensus.tmp.txt"
    shell:
        '''
        cut -f2- {input.bed} | awk 'NR==1; NR > 1 {{print $0 | "sort -k1,1 -k2,2n"}}' | cut -f6- > {params.trimmed}
        paste {input.boolean} {params.trimmed} > {output}
        rm {params.trimmed}
        '''

rule consensus_memechip:
    input:
        summit_bed = rules.consensus_summit.output.summit
    output:
        PRDM1motif = dir_out + "/consensus_peak/{group}/memechip/fimo_out_1/fimo.gff",
        RUNX1motif = dir_out + "/consensus_peak/{group}/memechip/fimo_out_2/fimo.gff",
        summit_fasta = dir_out + "/consensus_peak/{group}/memechip/{group}_summit.fasta"
    benchmark:
        dir_out + "/log/{group}.memechip.bmk"
    log:
        dir_out + "/log/{group}.memechip.log"
    params:
        ref = config["ref_gatk"],
        HOCO = "/home/jibzhang/motif_databases/HUMAN/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme",
        JASPAR = "/home/jibzhang/motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme",
        out_dir_peak = dir_out + "/consensus_peak/{group}/memechip"
    shell:
        '''
        bedtools getfasta -fi {params.ref} -bed {input.summit_bed} > {output.summit_fasta}
        meme-chip -db {params.HOCO} -db {params.JASPAR} {output.summit_fasta} -oc {params.out_dir_peak} -meme-mod anr \
        -minw 6 -maxw 20 -dna -meme-nmotifs 10 -centrimo-score 5.0 -centrimo-ethresh 10.0 -centrimo-local \
        -streme-pvt 0.05 -spamo-skip &> {log}
        '''

rule meme_motif_consensus:
    input:
        PRDM1_gff = rules.consensus_memechip.output.PRDM1motif,
        RUNX1_gff = rules.consensus_memechip.output.RUNX1motif,
        summit = rules.consensus_summit.output.summit
    output:
        PRDM1 = dir_out + "/consensus_peak/{group}/memechip/fimo_out_1/PRDM1_peaks.bed",
        RUNX1 = dir_out + "/consensus_peak/{group}/memechip/fimo_out_2/RUNX1_peaks.bed"
    benchmark:
        dir_out + "/log/{group}.mememotif_consensus.bmk"
    log:
        dir_out + "/log/{group}.mememotif_consensus.log"
    params:
        PRDM1_bed = dir_out + "/consensus_peak/{group}/memechip/fimo_out_1/PRDM1_motif.bed",
        RUNX1_bed = dir_out + "/consensus_peak/{group}/memechip/fimo_out_2/RUNX1_motif.bed",
    shell:
        '''
        tail -n +2 {input.PRDM1_gff} | awk -v OFS="\\t" '{{print $1,$4,$5}}' > {params.PRDM1_bed}
        bedtools intersect -wa -a {input.summit} -b {params.PRDM1_bed} | uniq > {output.PRDM1}
        tail -n +2 {input.RUNX1_gff} | awk -v OFS="\\t" '{{print $1,$4,$5}}' > {params.RUNX1_bed}
        bedtools intersect -wa -a {input.summit} -b {params.RUNX1_bed} | uniq > {output.RUNX1}
        '''

rule meme_overlap:
    input:
        PRDM1 = rules.meme_motif_consensus.output.PRDM1,
        RUNX1 = rules.meme_motif_consensus.output.RUNX1
    output:
        overlap = dir_out + "/consensus_peak/{group}/memechip/{group}_PRDM1RUNX1_peaks.bed",
        PRDM1_uniq = dir_out + "/consensus_peak/{group}/memechip/{group}_PRDM1uniq_peaks.bed",
        RUNX1_uniq = dir_out + "/consensus_peak/{group}/memechip/{group}_RUNX1uniq_peaks.bed",
        PRDM1_anno = dir_out + "/consensus_peak/{group}_PRDM1.annotatePeaks.txt",
        RUNX1_anno = dir_out + "/consensus_peak/{group}_RUNX1.annotatePeaks.txt",
        overlap_anno = dir_out + "/consensus_peak/{group}_PRDM1RUNX1.annotatePeaks.txt",
        PRDM1uniq_anno = dir_out + "/consensus_peak/{group}_PRDM1uniq.annotatePeaks.txt",
        RUNX1uniq_anno = dir_out + "/consensus_peak/{group}_RUNX1uniq.annotatePeaks.txt"
    benchmark:
        dir_out + "/log/{group}.motif_overlap.bmk"
    log:
        dir_out + "/log/{group}.motif_overlap.log"
    params:
        PRDM1_go = dir_out + "/consensus_peak/{group}/Homer/PRDM1_GO",
        RUNX1_go = dir_out + "/consensus_peak/{group}/Homer/RUNX1_GO",
        overlap_go = dir_out + "/consensus_peak/{group}/Homer/PRDM1RUNX1_GO",
        PRDM1uniq_go = dir_out + "/consensus_peak/{group}/Homer/PRDM1uniq_GO",
        RUNX1uniq_go = dir_out + "/consensus_peak/{group}/Homer/RUNX1uniq_GO",
        KORNAseq = "/coh_labs/jochan/PRDM1_NK/PRDM1KO_vs_WT.txt",
        overexpr6h="/coh_labs/jochan/PRDM1_NK/PRDM1_overexpres_vs_control_6h.txt",
        overexpr12h="/coh_labs/jochan/PRDM1_NK/PRDM1_overexpres_vs_control_12h.txt",
        fvsnof = "/coh_labs/jochan/PRDM1_NK/D13feeder_vs_D6nofeeder_Deseq2_DEGs.txt"
    shell:
        '''
        annotatePeaks.pl {input.PRDM1} hg38 -go {params.PRDM1_go} -gene {params.KORNAseq} {params.overexpr6h} {params.overexpr12h} {params.fvsnof} -cpu 4 > {output.PRDM1_anno} &> {log}
        annotatePeaks.pl {input.RUNX1} hg38 -go {params.RUNX1_go} -gene {params.KORNAseq} {params.overexpr6h} {params.overexpr12h} {params.fvsnof} -cpu 4 > {output.RUNX1_anno} &>> {log}
        bedtools intersect -a {input.PRDM1} -b {input.RUNX1} > {output.overlap}
        bedtools intersect -a {input.PRDM1} -b {input.RUNX1} -v > {output.PRDM1_uniq}
        bedtools intersect -a {input.RUNX1} -b {input.PRDM1} -v > {output.RUNX1_uniq}
        annotatePeaks.pl {output.overlap} hg38 -go {params.overlap_go} -gene {params.KORNAseq} {params.overexpr6h} {params.overexpr12h} {params.fvsnof} -cpu 4 > {output.overlap_anno} &>> {log}
        annotatePeaks.pl {output.PRDM1_uniq} hg38 -go {params.PRDM1uniq_go} -gene {params.KORNAseq} {params.overexpr6h} {params.overexpr12h}  {params.fvsnof} -cpu 4 > {output.PRDM1uniq_anno} &>> {log}
        annotatePeaks.pl {output.RUNX1_uniq} hg38 -go {params.RUNX1uniq_go} -gene {params.KORNAseq} {params.overexpr6h} {params.overexpr12h}  {params.fvsnof} -cpu 4 > {output.RUNX1uniq_anno} &>> {log}
        '''

rule ATAC_overlap:
    input:
        summit = rules.consensus_summit.output.summit,
        PRDM1 = rules.meme_overlap.output.PRDM1_uniq,
        RUNX1 = rules.meme_overlap.output.RUNX1_uniq,
        PRDM1RUNX = rules.meme_overlap.output.overlap
    output:
        OCR = dir_out + "/consensus_peak/{group}/{group}_ChIP_ATAC_overlap.bed",
        CCR = dir_out + "/consensus_peak/{group}/{group}_ChIP_ATAC_nooverlap.bed",
        PRDM1_OCR = dir_out + "/consensus_peak/{group}/memechip/{group}_PRDM1uniq_overlap_ATAC.bed",
        RUNX1_OCR = dir_out + "/consensus_peak/{group}/memechip/{group}_RUNX1uniq_overlap_ATAC.bed",
        PRDM1RUNX_OCR = dir_out + "/consensus_peak/{group}/memechip/{group}_PRDM1RUNX1_overlap_ATAC.bed",
        PRDM1_CCR = dir_out + "/consensus_peak/{group}/memechip/{group}_PRDM1uniq_nooverlap_ATAC.bed",
        RUNX1_CCR = dir_out + "/consensus_peak/{group}/memechip/{group}_RUNX1uniq_nooverlap_ATAC.bed",
        PRDM1RUNX_CCR = dir_out + "/consensus_peak/{group}/memechip/{group}_PRDM1RUNX1_nooverlap_ATAC.bed",
        PRDM1RUNX_OCR_anno = dir_out + "/homer/{group}_PRDM1RUNX1_overlap_ATAC.annotation.txt",
        PRDM1RUNX_CCR_anno = dir_out + "/homer/{group}_PRDM1RUNX1_nooverlap_ATAC.annotation.txt",
        PRDM1_OCR_anno=dir_out + "/homer/{group}_PRDM1uniq_overlap_ATAC.annotation.txt",
        PRDM1_CCR_anno=dir_out + "/homer/{group}_PRDM1uniq_nooverlap_ATAC.annotation.txt",
        RUNX1_OCR_anno=dir_out + "/homer/{group}_RUNX1uniq_overlap_ATAC.annotation.txt",
        RUNX1_CCR_anno=dir_out + "/homer/{group}_RUNX1uniq_nooverlap_ATAC.annotation.txt"
    benchmark:
        dir_out + "/log/{group}.ATAC_overlap.bmk"
    log:
        dir_out + "/log/{group}.ATAC_overlap.log"
    params:
        ATAC = "/coh_labs/jochan/PRDM1_ATAC/TOBIAS/{group}/macs2/{group}_consensus.bed",
        PRDM1RUNX_OCR_go= dir_out + "/homer/{group}_PRDM1RUNX1_overlap_ATAC/GO",
        PRDM1RUNX_CCR_go= dir_out + "/homer/{group}_PRDM1RUNX1_nooverlap_ATAC/GO",
        PRDM1_OCR_go=dir_out + "/homer/{group}_PRDM1uniq_overlap_ATAC/GO",
        PRDM1_CCR_go=dir_out + "/homer/{group}_PRDM1uniq_nooverlap_ATAC/GO",
        RUNX1_OCR_go=dir_out + "/homer/{group}_RUNX1uniq_overlap_ATAC/GO",
        RUNX1_CCR_go=dir_out + "/homer/{group}_RUNX1uniq_nooverlap_ATAC/GO",
        KORNAseq= "/coh_labs/jochan/PRDM1_NK/PRDM1KO_vs_WT.txt",
        overexpr6h="/coh_labs/jochan/PRDM1_NK/PRDM1_overexpres_vs_control_6h.txt",
        overexpr12h="/coh_labs/jochan/PRDM1_NK/PRDM1_overexpres_vs_control_12h.txt",
        fvsnof="/coh_labs/jochan/PRDM1_NK/D13feeder_vs_D6nofeeder_Deseq2_DEGs.txt"
    shell:
        '''
        bedtools intersect -a {input.summit} -b {params.ATAC} -wa > {output.OCR}
        bedtools intersect -a {input.summit} -b {params.ATAC} -v > {output.CCR}
        bedtools intersect -a {input.PRDM1} -b {params.ATAC} -wa > {output.PRDM1_OCR}
        bedtools intersect -a {input.PRDM1} -b {params.ATAC} -v > {output.PRDM1_CCR}
        bedtools intersect -a {input.RUNX1} -b {params.ATAC} -wa > {output.RUNX1_OCR}
        bedtools intersect -a {input.RUNX1} -b {params.ATAC} -v > {output.RUNX1_CCR}
        bedtools intersect -a {input.PRDM1RUNX} -b {params.ATAC} -wa > {output.PRDM1RUNX_OCR}
        bedtools intersect -a {input.PRDM1RUNX} -b {params.ATAC} -v > {output.PRDM1RUNX_CCR}
        annotatePeaks.pl {output.PRDM1_OCR} hg38 -go {params.PRDM1_OCR_go} -gene {params.KORNAseq} {params.overexpr6h} {params.overexpr12h} {params.fvsnof} -cpu 4 > {output.PRDM1_OCR_anno} &> {log}
        annotatePeaks.pl {output.PRDM1_CCR} hg38 -go {params.PRDM1_CCR_go} -gene {params.KORNAseq} {params.overexpr6h} {params.overexpr12h} {params.fvsnof} -cpu 4 > {output.PRDM1_CCR_anno} &>> {log}
        annotatePeaks.pl {output.RUNX1_OCR} hg38 -go {params.RUNX1_OCR_go} -gene {params.KORNAseq} {params.overexpr6h} {params.overexpr12h} {params.fvsnof} -cpu 4 > {output.RUNX1_OCR_anno} &>> {log}
        annotatePeaks.pl {output.RUNX1_CCR} hg38 -go {params.RUNX1_CCR_go} -gene {params.KORNAseq} {params.overexpr6h} {params.overexpr12h} {params.fvsnof} -cpu 4 > {output.RUNX1_CCR_anno} &>> {log}
        annotatePeaks.pl {output.PRDM1RUNX_OCR} hg38 -go {params.PRDM1RUNX_OCR_go} -gene {params.KORNAseq} {params.overexpr6h} {params.overexpr12h} {params.fvsnof} -cpu 4 > {output.PRDM1RUNX_OCR_anno} &>> {log}
        annotatePeaks.pl {output.PRDM1RUNX_CCR} hg38 -go {params.PRDM1RUNX_CCR_go} -gene {params.KORNAseq} {params.overexpr6h} {params.overexpr12h} {params.fvsnof} -cpu 4 > {output.PRDM1RUNX_CCR_anno} &>> {log}
        '''

rule Homer_motif_consensus:
    input:
        fasta = rules.consensus_memechip.output.summit_fasta,
        summits = rules.consensus_summit.output.summit,
        PRDM1 = dir_out + "/consensus_peak/nofeeder/Homer/motif/knownResults/known1.motif",
        RUNX1 = dir_out + "/consensus_peak/nofeeder/Homer/motif/knownResults/known9.motif"
    output:
        PRDM1 = dir_out + "/consensus_peak/{group}/Homer/motif/{group}_PRDM1_summits.bed",
        RUNX1 = dir_out + "/consensus_peak/{group}/Homer/motif/{group}_RUNX1_summits.bed"
    benchmark:
        dir_out + "/log/{group}.homer_motif.bmk"
    params:
        PRDM1_bed = dir_out + "/consensus_peak/{group}/Homer/motif/motif1.bed",
        PRDM1_region = dir_out + "/consensus_peak/{group}/Homer/motif/PRDM1.bed",
        RUNX1_bed = dir_out + "/consensus_peak/{group}/Homer/motif/motif2.bed",
        RUNX1_region = dir_out + "/consensus_peak/{group}/Homer/motif/RUNX1.bed"
    shell:
        '''
        scanMotifGenomeWide.pl {input.PRDM1} {input.fasta} > {params.PRDM1_bed}
        cat {params.PRDM1_bed} | awk '{{split($2,a,"[:-]"); printf("%s\\t%s\\t%s\\n",a[1],a[2],a[3]);}}' > {params.PRDM1_region}
        bedtools intersect -wa -a {input.summits} -b {params.PRDM1_region} | uniq > {output.PRDM1}

        scanMotifGenomeWide.pl {input.RUNX1} {input.fasta} > {params.RUNX1_bed}
        cat {params.RUNX1_bed} | awk '{{split($2,a,"[:-]"); printf("%s\\t%s\\t%s\\n",a[1],a[2],a[3]);}}' > {params.RUNX1_region}
        bedtools intersect -wa -a {input.summits} -b {params.RUNX1_region} | uniq > {output.RUNX1}
        '''

rule Homer_overlap:
    input:
        PRDM1 = rules.Homer_motif_consensus.output.PRDM1,
        RUNX1 = rules.Homer_motif_consensus.output.RUNX1
    output:
        overlap = dir_out + "/consensus_peak/{group}/Homer/motif/{group}_PRDM1RUNX1_summits.bed",
        PRDM1_uniq= dir_out + "/consensus_peak/{group}/memechip/{group}_PRDM1uniq_peaks.bed",
        RUNX1_uniq=dir_out + "/consensus_peak/{group}/memechip/{group}_RUNX1uniq_peaks.bed"
    benchmark:
        dir_out + "/log/{group}.homer_overlap.bmk"
    shell:
        '''
        bedtools intersect -a {input.PRDM1} -b {input.RUNX1} > {output.overlap}
        bedtools intersect -a {input.PRDM1} -b {input.RUNX1} -v > {output.PRDM1_uniq}
        bedtools intersect -a {input.RUNX1} -b {input.PRDM1} -v > {output.RUNX1_uniq}
        '''

rule done_consensus:
    input:
        annotate = rules.annotate_boolean.benchmark,
        meme = rules.ATAC_overlap.benchmark
    output:
        dir_out + "/log/{group}.done_consensus.txt"
    shell:
        '''
        touch {output}
        '''
