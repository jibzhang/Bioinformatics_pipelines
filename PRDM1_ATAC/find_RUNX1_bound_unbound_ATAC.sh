#!/bin/sh

ATAC_PRDM1="/scratch/jibzhang/project/PRDM1_ATAC/TOBIAS/BINDdetect/_RUNX1_HUMAN.H11MO.0.A/beds/_RUNX1_HUMAN.H11MO.0.A_all.bed"
bw_feeder="/scratch/jibzhang/project/PRDM1_ATAC/TOBIAS/feeder/correct/feeder_corrected.bw"
bw_nofeeder="/scratch/jibzhang/project/PRDM1_ATAC/TOBIAS/nofeeder/correct/nofeeder_corrected.bw"

nofeeder_ChIP_PRDM1="/scratch/jibzhang/project/PRDM1_NK/consensus_peak/nofeeder/Homer/motif/nofeeder_PRDM1uniq_peaks.bed"
nofeeder_ChIP_RUNX1="/scratch/jibzhang/project/PRDM1_NK/consensus_peak/nofeeder/Homer/motif/nofeeder_RUNX1uniq_peaks.bed"
nofeeder_ChIP_PRDM1RUNX1="/scratch/jibzhang/project/PRDM1_NK/consensus_peak/nofeeder/Homer/motif/nofeeder_PRDM1RUNX1_peaks.bed"
nofeeder_PRDM1_bound="/scratch/jibzhang/project/PRDM1_ATAC/TOBIAS/BINDdetect/_RUNX1_HUMAN.H11MO.0.A/beds/nofeeder_PRDM1only_bound.bed"
nofeeder_PRDM1_unbound="/scratch/jibzhang/project/PRDM1_ATAC/TOBIAS/BINDdetect/_RUNX1_HUMAN.H11MO.0.A/beds/nofeeder_PRDM1only_unbound.bed"
nofeeder_RUNX1_bound="/scratch/jibzhang/project/PRDM1_ATAC/TOBIAS/BINDdetect/_RUNX1_HUMAN.H11MO.0.A/beds/nofeeder_RUNX1only_bound.bed"
nofeeder_RUNX1_unbound="/scratch/jibzhang/project/PRDM1_ATAC/TOBIAS/BINDdetect/_RUNX1_HUMAN.H11MO.0.A/beds/nofeeder_RUNX1only_unbound.bed"
nofeeder_PRDM1RUNX1_bound="/scratch/jibzhang/project/PRDM1_ATAC/TOBIAS/BINDdetect/_RUNX1_HUMAN.H11MO.0.A/beds/nofeeder_PRDM1RUNX1_bound.bed"
nofeeder_PRDM1RUNX1_unbound="/scratch/jibzhang/project/PRDM1_ATAC/TOBIAS/BINDdetect/_RUNX1_HUMAN.H11MO.0.A/beds/nofeeder_PRDM1RUNX1_unbound.bed"

feeder_ChIP_PRDM1="/scratch/jibzhang/project/PRDM1_NK/consensus_peak/feeder_d13/Homer/motif/feeder_d13_PRDM1uniq_peaks.bed"
feeder_ChIP_RUNX1="/scratch/jibzhang/project/PRDM1_NK/consensus_peak/feeder_d13/Homer/motif/feeder_d13_RUNX1uniq_peaks.bed"
feeder_ChIP_PRDM1RUNX1="/scratch/jibzhang/project/PRDM1_NK/consensus_peak/feeder_d13/Homer/motif/feeder_d13_PRDM1RUNX1_peaks.bed"
feeder_PRDM1_bound="/scratch/jibzhang/project/PRDM1_ATAC/TOBIAS/BINDdetect/_RUNX1_HUMAN.H11MO.0.A/beds/feeder_PRDM1only_bound.bed"
feeder_PRDM1_unbound="/scratch/jibzhang/project/PRDM1_ATAC/TOBIAS/BINDdetect/_RUNX1_HUMAN.H11MO.0.A/beds/feeder_PRDM1only_unbound.bed"
feeder_RUNX1_bound="/scratch/jibzhang/project/PRDM1_ATAC/TOBIAS/BINDdetect/_RUNX1_HUMAN.H11MO.0.A/beds/feeder_RUNX1only_bound.bed"
feeder_RUNX1_unbound="/scratch/jibzhang/project/PRDM1_ATAC/TOBIAS/BINDdetect/_RUNX1_HUMAN.H11MO.0.A/beds/feeder_RUNX1only_unbound.bed"
feeder_PRDM1RUNX1_bound="/scratch/jibzhang/project/PRDM1_ATAC/TOBIAS/BINDdetect/_RUNX1_HUMAN.H11MO.0.A/beds/feeder_PRDM1RUNX1_bound.bed"
feeder_PRDM1RUNX1_unbound="/scratch/jibzhang/project/PRDM1_ATAC/TOBIAS/BINDdetect/_RUNX1_HUMAN.H11MO.0.A/beds/feeder_PRDM1RUNX1_unbound.bed"

bedtools intersect -wa -a $ATAC_PRDM1 -b $nofeeder_ChIP_PRDM1 |uniq > $nofeeder_PRDM1_bound
bedtools intersect -v -a $ATAC_PRDM1 -b $nofeeder_ChIP_RUNX1 |uniq > $nofeeder_PRDM1_unbound

bedtools intersect -wa -a $ATAC_PRDM1 -b $feeder_ChIP_PRDM1 |uniq > $feeder_PRDM1_bound
bedtools intersect -v -a $ATAC_PRDM1 -b $feeder_ChIP_RUNX1 |uniq > $feeder_PRDM1_unbound

bedtools intersect -wa -a $ATAC_PRDM1 -b $nofeeder_ChIP_RUNX1 |uniq > $nofeeder_RUNX1_bound
bedtools intersect -v -a $ATAC_PRDM1 -b $nofeeder_ChIP_RUNX1 |uniq > $nofeeder_RUNX1_unbound

bedtools intersect -wa -a $ATAC_PRDM1 -b $feeder_ChIP_RUNX1 |uniq > $feeder_RUNX1_bound
bedtools intersect -v -a $ATAC_PRDM1 -b $feeder_ChIP_RUNX1 |uniq > $feeder_RUNX1_unbound

bedtools intersect -wa -a $ATAC_PRDM1 -b $nofeeder_ChIP_PRDM1RUNX1 |uniq > $nofeeder_PRDM1RUNX1_bound
bedtools intersect -v -a $ATAC_PRDM1 -b $nofeeder_ChIP_PRDM1RUNX1 |uniq > $nofeeder_PRDM1RUNX1_unbound

bedtools intersect -wa -a $ATAC_PRDM1 -b $feeder_ChIP_PRDM1RUNX1 |uniq > $feeder_PRDM1RUNX1_bound
bedtools intersect -v -a $ATAC_PRDM1 -b $feeder_ChIP_PRDM1RUNX1 |uniq > $feeder_PRDM1RUNX1_unbound

TOBIAS PlotHeatmap --TFBS $feeder_PRDM1_bound $feeder_PRDM1_unbound --TFBS-labels "bound" "unbound" \
--TFBS $nofeeder_PRDM1_bound $nofeeder_PRDM1_unbound --TFBS-labels "bound" "unbound" \
--signals $bw_feeder $bw_nofeeder --output /scratch/jibzhang/project/PRDM1_ATAC/TOBIAS/TFBS_feedervsnofeeder_PRDM1only_plots.heatmap.pdf \
--signal_labels feeder nofeeder --share_colorbar --sort_by -1 --title "PRDM1only_feedervsnofeeder_heatmap"

TOBIAS PlotHeatmap --TFBS $feeder_RUNX1_bound $feeder_RUNX1_unbound --TFBS-labels "bound" "unbound" \
--TFBS $nofeeder_RUNX1_bound $nofeeder_RUNX1_unbound --TFBS-labels "bound" "unbound" \
--signals $bw_feeder $bw_nofeeder --output /scratch/jibzhang/project/PRDM1_ATAC/TOBIAS/TFBS_feedervsnofeeder_RUNX1only_plots.heatmap.pdf \
--signal_labels feeder nofeeder --share_colorbar --sort_by -1 --title "RUNX1only_feedervsnofeeder_heatmap"

TOBIAS PlotHeatmap --TFBS $feeder_PRDM1RUNX1_bound $feeder_PRDM1RUNX1_unbound --TFBS-labels "bound" "unbound" \
--TFBS $nofeeder_PRDM1RUNX1_bound $nofeeder_PRDM1RUNX1_unbound --TFBS-labels "bound" "unbound" \
--signals $bw_feeder $bw_nofeeder --output /scratch/jibzhang/project/PRDM1_ATAC/TOBIAS/TFBS_feedervsnofeeder_PRDM1RUNX1_plots.heatmap.pdf \
--signal_labels feeder nofeeder --share_colorbar --sort_by -1 --title "PRDM1RUNX1_feedervsnofeeder_heatmap"