rule cnv_facets:
    input:
        bam = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam"
    output:
        vcf = dir_out + "/{sample}/" + squence_type + "/cnv_facets/{sample}.cnv_facets.csv.gz",
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.cnv_facets.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.cnv_facets.bmk"
    params:
        ref=config["ref_bwa"],
        out_prefix = dir_out + "/{sample}/" + squence_type + "/cnv_facets/{sample}.cnv_facets",
        dbsnp_common=config["dbsnp_cnvfacets"],
        facetscnv_gbuild=config["facetscnv_gbuild"]
    threads:
        config["threads_facets"]
    shell:
        '''
        cnv_facets.R -N {threads} -t {input.bam} -vcf {params.dbsnp_common} -o {params.out_prefix} \
        --gbuild {params.facetscnv_gbuild} --rnd-seed 10 &> {log}
        '''

#export PATH="/home/zgu_labs/anaconda3/envs/copywriteR/bin:$PATH"
rule copywriteR:
    input:
        bam=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam"
    output:
        tsv=dir_out + "/{sample}/" + squence_type + "/copywriteR/results_output.tsv"
    log:        dir_out + "/{sample}/" + squence_type + "/log/{sample}.copywriteR.log"
    benchmark:  dir_out + "/{sample}/" + squence_type + "/log/{sample}.copywriteR.bmk"
    params:
        helper=config["copywriteR_helper"],
        dir_out=dir_out + "/{sample}/" + squence_type + "/copywriteR"
    threads: 4
    shell:
        '''
        Rscript_copywriteR /home/jibzhang/bin/R_script/run_copywriteR_nopair.R \
        {threads} {params.helper} {input.bam} {params.dir_out}
        '''

rule cnmops_wes_germline:
    input:
	    input_bams=expand(dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam",sample=samplelist)
    output:
	    cnv_beds=expand(dir_out + "/cnmops_wgs_germline/{sample}.bed",sample=samplelist)
    params:
	    result_dir=dir_out + "/cnmops_wes_germline/"
    log:
        dir_out + "/cnmops_wes_germline/cnmops_wes_germline.log"
    benchmark:
        dir_out + "/cnmops_wes_germline/cnmops_wes_germline.bmk"
    threads: 4
    shell:
        '''
        Rscript /coh_labs/jochan/pipeline/WES/scripts/cnmops_wes.R {params.result_dir} 500 {threads} {input.input_bams} 2> {log}
        '''

rule ControlFREEC_estimate_gender:
    input:
        bam = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam",
        bai = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam.bai"
    output:
        chrX   = "temp/ControlFREEC/{sample}/mosdepth/chrX.mosdepth.global.dist.txt",
        chrY   = "temp/ControlFREEC/{sample}/mosdepth/chrY.mosdepth.global.dist.txt",
        gender = "temp/ControlFREEC/{sample}_gender.txt",
    params:
        prefixX = "temp/ControlFREEC/{sample}/mosdepth/chrX",
        prefixY = "temp/ControlFREEC/{sample}/mosdepth/chrY",
    threads: 1
    log:    "logs/ControlFREEC_estimate_gender/{sample}.log"
    benchmark: "benchmarks/ControlFREEC/{sample}.txt"
    shell:
        '''
        >&2 echo 'Estimate chrX coverage';
        mosdepth -n --fast-mode --chrom X {params.prefixX} {input.bam};
        >&2 echo 'Estimate chrY coverage';
        mosdepth -n --fast-mode --chrom Y {params.prefixY} {input.bam};
        XCOV=$(grep 'total\t10\t' {output.chrX} | cut -f 3);
        YCOV=$(grep 'total\t10\t' {output.chrY} | cut -f 3);
        if (( $(echo \$XCOV > 0.1 && ! $YCOV > 0.1\ | bc) )); then gender=XX; else gender=XY; fi; 
        echo $gender > {output.gender} 
        '''

rule ControlFREEC_calling:
    input:
        bam = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam",
        bai = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam.bai",
        gender = rules.ControlFREEC_estimate_gender.output.gender,
        chrdir = config["chromosome_sequences_dir"],
        chrlen = config["chromosome_lengths"],
    output:
        cnv   = temp("temp/ControlFREEC/{sample}/{sample}.bam_CNVs"),
        ratio = temp("temp/ControlFREEC/{sample}/{sample}.bam_ratio.txt"),
    params:
        outdir = "temp/ControlFREEC/{sample}",
    threads: 28
    log:    "logs/ControlFREEC_calling/{sample}.log"
    shell:
        "(module load tools ngs gcc/8.2.0 anaconda3/4.4.0 "
        "intel/perflibs/64/2019_update3 intel/redist/2019_update2 R/3.5.0 "
        "samtools/1.9 bedtools/2.27.1 sambamba/0.6.7 control-freec/11.5; "
        "tools/ControlFREEC/wrapper.py "
        "--input {input.bam} "
        "--sex $(cat {input.gender}) "
        "--outputDir {params.outdir} "
        "--maxThreads {threads} "
        "--sambamba $(which sambamba) "
        "--SambambaThreads 2 "
        "--coefficientOfVariation 0.05 "
        "--minCNAlength 1 "
        "--forceGCcontentNormalization 0 "
        "--ploidy 2 "
        "--chrFiles {input.chrdir} "
        "--chrLenFile {input.chrlen} "
        ">&2 "
        ") 2> {log}"

rule ControlFREEC_significance:
    input:
        cnv   = "temp/ControlFREEC/{sample}/{sample}.bam_CNVs",
        ratio = "temp/ControlFREEC/{sample}/{sample}.bam_ratio.txt",
    output:
        txt   = "temp/ControlFREEC/{sample}/{sample}.bam_CNVs.p.value.txt"
    threads: 1
    log: "logs/ControlFREEC_significance/{sample}.log"
    shell:
        "(module load tools ngs gcc/8.2.0 intel/perflibs/64/2019_update3 intel/redist/2019_update2 R/3.5.0; "
        # "module load module load rdxplorer/3.2; "
        "Rscript tools/ControlFREEC/assess_significance.R {input.cnv} {input.ratio} "
        ") 2> {log}"

rule ControlFREEC_filter:
    input:
        "temp/ControlFREEC/{sample}/{sample}.bam_CNVs.p.value.txt",
    output:
        p1 = "temp/ControlFREEC/{sample}/{sample}_1.txt",
        p5 = "temp/ControlFREEC/{sample}/{sample}_5.txt",
    threads: 1
    log: "logs/ControlFREEC_filter/{sample}.log"
    shell:
        "perl -ane 'print if $F[5] <= 0.05' < {input} > {output.p5}; "
        "perl -ane 'print if $F[5] <= 0.01' < {input} > {output.p1}; "

rule ControlFREEC_postprocess:
    input:
        txt = "temp/ControlFREEC/{sample}/{sample}.bam_CNVs.p.value.txt"
    output:
        bed = "results/ControlFREEC/{sample}.bed"
    log: "logs/ControlFREEC_postprocess/{sample}.log"
    shell:
        "("
        "tail -n +2 {input.txt} | cut -f 1-4,6 > {output.bed} "
        ") 2> {log}"

rule cnv_WES_nompair:
    input:
        copyWriteR = rules.copywriteR.output.tsv
    output:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.cnv_WES_nopair.done"
    shell:
        '''
        touch {output}
        '''