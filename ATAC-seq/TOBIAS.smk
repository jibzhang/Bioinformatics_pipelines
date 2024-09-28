rule merge_bam:
    input:
        lambda wildcards: [dir_out + "/" + rep + "/" + squence_type + "/bam/" + rep +".bam" for rep in contrast[wildcards.group]]
    output:
        bam=dir_out + "/TOBIAS/{group}/bam/{group}.bam"
    log:
        dir_out + "/TOBIAS/{group}/log/{group}.merge.log"
    benchmark:
        dir_out + "/TOBIAS/{group}/log/{group}.merge.bmk"
    params:
        tmp_n=lambda wildcards: len(contrast[wildcards.group])
    shell:
        '''
        if [ {params.tmp_n} -gt 1 ] ; then
        samtools merge {output.bam} {input}
        else
        ln -s {input} {output.bam}
        fi
        '''

rule consensus_group:
    input:
        lambda wildcards: [dir_out + "/" + rep + "/" + squence_type + "/macs2/" + rep +"_filtered_peaks.bed" for rep in contrast[wildcards.group]]
    output:
        consensus=dir_out + "/TOBIAS/{group}/macs2/{group}_consensus.bed"
    log:
        dir_out + "/TOBIAS/{group}/log/{group}.consensus.log"
    benchmark:
        dir_out + "/TOBIAS/{group}/log/{group}.consensus.bmk"
    params:
        merged = dir_out + "/TOBIAS/{group}/macs2/{group}_merged_peaks.txt",
        boolean = dir_out + "/TOBIAS/{group}/macs2/{group}_boolean.txt",
        intersect= dir_out + "/consensus_peak/{group}/{group}_boolean.intersect.txt",
        plot= dir_out + "/consensus_peak/{group}/{group}_boolean.intersect.plot.pdf",
        names = lambda wildcards, input: ','.join(input.split("/")[4] for input in input)
    shell:
        '''
        sort -k1,1 -k2,2n {input} | mergeBed -c 2,3,4,5,6,7,8,9,10 -o collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse > {params.merged}
        /home/jibzhang/bin/macs2_merged_expand.py {params.merged} {params.names} {params.boolean} --min_replicates 2 --is_narrow_peak
        awk -v FS='\\t' -v OFS='\\t' 'FNR > 1 {{ print $1, $2, $3, $4, "0", "+" }}' {params.boolean} > {output.consensus}
        /home/jibzhang/bin/R_script/plot_peak_intersect.r -i {params.intersect} -o {params.plot}
        '''

rule MACS2_merged_broad:
    input:
        lambda wildcards: [dir_out + "/" + rep + "/" + squence_type + "/bam/" + rep +".minimal.bedpe" for rep in contrast[wildcards.group]]
    output:
        peak = dir_out + "/TOBIAS/{group}/macs2_broad/{group}_peaks.broadPeak"
    log:
        dir_out + "/TOBIAS/{group}/log/{group}.macs2_broad.log"
    benchmark:
        dir_out + "/TOBIAS/{group}/log/{group}.macs2_broad.bmk"
    params:
        prefix = "{group}",
        dir = dir_out + "/TOBIAS/{group}/macs2_broad"
    shell:
        '''
        macs2 callpeak -t {input} -f BEDPE --keep-dup all -n {params.prefix} -g 2805665311 --broad --broad-cutoff 0.05 --outdir {params.dir} 2> {log}
        '''

rule TOBIAS_correct:
    input:
        bam = rules.merge_bam.output.bam,
        peak = dir_out + "/consensus_peak/D6nofeeder_D13feeder_consensus.bed",
        broad = rules.MACS2_merged_broad.output.peak
    output:
        dir_out + "/TOBIAS/{group}/correct/{group}_corrected.bw"
    benchmark:
        dir_out + "/TOBIAS/{group}/log/TOBIAS_correct.bmk",
    params:
        dir = dir_out + "/TOBIAS/{group}/correct",
        ref = config["ref_bwa"],
        prefix = "{group}"
    threads:
        6
    shell:
        '''
        TOBIAS ATACorrect --bam {input.bam} --genome {params.ref} --peaks {input.peak} --cores {threads} --prefix {params.prefix} --outdir {params.dir}
        '''

rule TOBIAS_score:
    input:
        signal = rules.TOBIAS_correct.output,
        peak = dir_out + "/consensus_peak/D6nofeeder_D13feeder_consensus.bed"
    output:
        dir_out + "/TOBIAS/{group}/{group}_footprints.bw"
    benchmark:
        dir_out + "/TOBIAS/{group}/log/TOBIAS_score.bmk"
    threads:
        6
    shell:
        '''
        TOBIAS FootprintScores --signal {input.signal} --regions {input.peak} --cores {threads} --output {output}
        '''

rule TOBIAS_bindetect:
    input:
        D6nofeeder_signals = dir_out + "/TOBIAS/CD56Dim/CD56Dim_footprints.bw",
        D13feeder_signals = dir_out + "/TOBIAS/nofeeder/nofeeder_footprints.bw",
        peaks = dir_out + "/consensus_peak/D6nofeeder_D13feeder_consensus_anno.bed",
        header = dir_out + "/peaks.annotation.header"
    output:
        result = dir_out + "/TOBIAS/BINDdetect/bindetect_results.txt",
        PRDM1_feeder_bound = dir_out + "/TOBIAS/BINDdetect/_PRDM1_HUMAN.H11MO.0.A/beds/_PRDM1_HUMAN.H11MO.0.A_nofeeder_bound.bed",
        PRDM1_nofeeder_bound = dir_out + "/TOBIAS/BINDdetect/_PRDM1_HUMAN.H11MO.0.A/beds/_PRDM1_HUMAN.H11MO.0.A_CD56Dim_bound.bed",
        PRDM1_feeder_unbound = dir_out + "/TOBIAS/BINDdetect/_PRDM1_HUMAN.H11MO.0.A/beds/_PRDM1_HUMAN.H11MO.0.A_nofeeder_unbound.bed",
        PRDM1_nofeeder_unbound = dir_out + "/TOBIAS/BINDdetect/_PRDM1_HUMAN.H11MO.0.A/beds/_PRDM1_HUMAN.H11MO.0.A_CD56Dim_unbound.bed",
        RUNX1_all = dir_out + "/TOBIAS/BINDdetect/_RUNX1_HUMAN.H11MO.0.A/beds/_RUNX1_HUMAN.H11MO.0.A_all.bed",
        PRDM1_all = dir_out + "/TOBIAS/BINDdetect/_PRDM1_HUMAN.H11MO.0.A/beds/_PRDM1_HUMAN.H11MO.0.A_all.bed"
    benchmark:
        dir_out + "/log/TOBIAS_bindetect.bmk"
    params:
        dir = dir_out + "/TOBIAS/BINDdetect",
        ref = config["ref_bwa"],
        motif = "/home/jibzhang/motif_databases/HUMAN/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme"
    threads:
        6
    shell:
        '''
        TOBIAS BINDetect --motifs {params.motif} --signals  {input.D13feeder_signals} {input.D6nofeeder_signals} --genome {params.ref} --peaks {input.peaks} --peak-header {input.header} --outdir {params.dir} --cond-names nofeeder CD56Dim --cores {threads}
        '''

rule TOBIAS_aggregate:
    input:
        D6nofeeder_signals = dir_out + "/TOBIAS/nofeeder/correct/nofeeder_corrected.bw",
        D13feeder_signals = dir_out + "/TOBIAS/feeder_d13/correct/feeder_d13_corrected.bw",
        D2nofeeder_signals = dir_out + "/TOBIAS/D2_nofeeder/correct/D2_nofeeder_corrected.bw",
        D28feeder_signals = dir_out + "/TOBIAS/D28_feeder/correct/D28_feeder_corrected.bw",
        PRDM1_all = rules.TOBIAS_bindetect.output.PRDM1_all,
        RUNX1_all = rules.TOBIAS_bindetect.output.RUNX1_all
    output:
        PRDM1 = dir_out + "/TOBIAS/TFBS_D13feedervsD6nofeeder_PRDM1_plots.aggregate.png",
        RUNX1 = dir_out + "/TOBIAS/TFBS_D13feedervsD6nofeeder_RUNX1_plots.aggregate.png"
    benchmark:
        dir_out + "/log/TOBIAS_aggregate.bmk"
    shell:
        '''
        TOBIAS PlotAggregate --TFBS {input.PRDM1_all} --signals {input.D2nofeeder_signals} {input.D6nofeeder_signals} {input.D13feeder_signals} {input.D28feeder_signals} --signal_labels D2_nofeeder D6_nofeeder D13_feeder D28_feeder --output {output.PRDM1} --share_y both --plot_boundaries --signal-on-x --title "PRDM1"
        TOBIAS PlotAggregate --TFBS {input.RUNX1_all} --signals {input.D2nofeeder_signals} {input.D6nofeeder_signals} {input.D13feeder_signals} {input.D28feeder_signals} --signal_labels D2_nofeeder D6_nofeeder D13_feeder D28_feeder --output {output.RUNX1} --share_y both --plot_boundaries --signal-on-x --title "RUNX1"
        '''

rule TOBIAS_plotHeatmap:
    input:
        feeder_bound = rules.TOBIAS_bindetect.output.PRDM1_feeder_bound,
        feeder_unbound = rules.TOBIAS_bindetect.output.PRDM1_feeder_unbound,
        feeder_signals = dir_out + "/TOBIAS/feeder_d13/correct/feeder_d13_corrected.bw",
        nofeeder_bound = rules.TOBIAS_bindetect.output.PRDM1_nofeeder_bound,
        nofeeder_unbound = rules.TOBIAS_bindetect.output.PRDM1_nofeeder_unbound,
        nofeeder_signals = dir_out + "/TOBIAS/D6_nofeeder/correct/D6_nofeeder_corrected.bw"
    output:
        dir_out + "/TOBIAS/TFBS_D13feedervsD6nofeeder_PRDM1_plots.heatmap.pdf"
    benchmark:
        dir_out + "/log/TOBIAS_plotheatmap.bmk"
    threads:
        6
    shell:
        '''
        TOBIAS PlotHeatmap --TFBS {input.feeder_bound} {input.feeder_unbound} --TFBS-labels "bound" "unbound" --TFBS {input.nofeeder_bound} {input.nofeeder_unbound} \
        --TFBS-labels "bound" "unbound" --signals {input.feeder_signals} {input.nofeeder_signals} --output {output} --signal_labels feeder nofeeder --share_colorbar \
        --sort_by -1  --title "PRDM1_d13feeder_vs_d6nofeeder_heatmap"
        '''

rule TOBIAS_done:
    input:
        plotheatmap = rules.TOBIAS_plotHeatmap.benchmark,
        aggregate = rules.TOBIAS_aggregate.benchmark
    output:
        dir_out + "/log/done_TOBIAS.txt"
    shell:
        '''
        touch {output}
        '''
