#cd /ref_genomes/genome_anno/human/gtf/V102
#python \
#/home/zgu_labs/anaconda3/pkgs/bioconductor-dexseq-1.36.0-r40_0/lib/R/library/DEXSeq/python_scripts/dexseq_prepare_annotation.py -r no \
#    Homo_sapiens.GRCh38.V102.withChr.gtf \
#    Homo_sapiens.GRCh38.V102.withChr.HTseqExon.gff

rule Transform:
    input:
        bmk=dir_out+"/temp/"+squence_type+"/log/star/{sample}_star.bmk",
        junction=dir_out+"/temp/"+squence_type+"/{sample}/star/SJ.out.tab"
    output:
        dir_out+"/temp/"+squence_type+"/star/{sample}_junctions.bed"
    log:        dir_out + "/{sample}/" + squence_type + "/log/{sample}.junctions.log"
    benchmark:  dir_out + "/{sample}/" + squence_type + "/log/{sample}.junctions.bmk"
    shell:
        '''
        awk 'BEGIN {{OFS="\t"}}{{print $1, $2-20-1, $3+20, "JUNCBJ"NR, $7, ($4 == 1)? "+":"-",$2-20-1, $3+20,"255,0,0", 2,"20,20","0,300"}}' {input.junction} > {output}
        '''

rule PSI:
    input:
        bed=rules.Transform.output,
        bam=dir_out+"/temp/"+squence_type+"/star/{sample}.bam"
    output:
        dir_out+"/temp/"+squence_type+"/star/{sample}_exonic_parts.psi"
    log:        dir_out + "/{sample}/" + squence_type + "/log/{sample}.PSI.log"
    benchmark:  dir_out + "/{sample}/" + squence_type + "/log/{sample}.PSI.bmk"
    params:
        gtf="/home/jibzhang/reference/Mouse_ensembl_107/Mus_musculus.GRCm39.107_reduced.gtf",
        prefix="{sample}",
        psi="{sample}_exonic_parts.psi"
    shell:
        '''
        /home/jibzhang/bin/PSI.sh StartPSI {params.gtf} 151 {input.bam} {input.bed} {params.prefix} 2> {log}
        mv {params.psi} {output}
        '''

rule HTSeq_count:
    input:
        bam=    dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam",
        bmk=    dir_out + "/{sample}/" + squence_type + "/log/{sample}_markDup.bmk",
    output:
        HTSeq_count=    dir_out + "/{sample}/" + squence_type + "/HTSeq/{sample}.HTSeq",
    log:        dir_out + "/{sample}/" + squence_type + "/log/{sample}.HTSeq_count.log"
    benchmark:  dir_out + "/{sample}/" + squence_type + "/log/{sample}.HTSeq_count.bmk"
    params:
        htseq_count=config["htseq-count"],
        python=config["python"],
        gtf=config["rat_gtf"]
    shell:
        '''
        {params.htseq_count} -m intersection-strict -f bam -r pos -s no -a 10 {input.bam} {params.gtf} 1> {output.HTSeq_count} 2> {log}
        '''

rule HTSeqDUX4:
    input:
        bed="/ref_genomes/genome_anno/human/gtf/V102/Homo_sapiens.GRCh38.V102.withChr.DUX4.bed",
        bam=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam",
        htseq=rules.HTSeq_count.output.HTSeq_count,
        bmk=dir_out + "/{sample}/" + squence_type + "/log/{sample}.HTSeq_count.bmk"
    output:
        HTSeq=dir_out + "/{sample}/" + squence_type + "/HTSeq_DUX4/{sample}.DUX4patched.HTSeq",
    log:        dir_out + "/{sample}/" + squence_type + "/log/{sample}.HTSeqDUX4.log"
    benchmark:  dir_out + "/{sample}/" + squence_type + "/log/{sample}.HTSeqDUX4.bmk"
    params:
        pl=""
    shell:
        '''
        HTSeqDUX4patch.pl -f {input.bed} -b {input.bam} -h {input.htseq} -o {output.HTSeq}
        '''
rule DEXSeq:
    input:
        bam=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam",
        bmk=dir_out + "/{sample}/" + squence_type + "/log/{sample}_markDup.bmk",
    output:
        DEXSeq=dir_out + "/{sample}/" + squence_type + "/DEXSeq/{sample}.DEXSeq",
    log:        dir_out + "/{sample}/" + squence_type + "/log/{sample}.DEXSeq.log"
    benchmark:  dir_out + "/{sample}/" + squence_type + "/log/{sample}.DEXSeq.bmk"
    params:
        dexseq="/home/zgu_labs/bin/dexseq_count.py",
        gff=config["exon_gff"],
    shell:
        '''
        python {params.dexseq} -p yes -r pos -s no -a 10 -f bam {params.gff} {input.bam} {output.DEXSeq} 2> {log}
        '''

rule DEXSeqDUX4:
    input:
        bed="/ref_genomes/genome_anno/human/gtf/V102/Homo_sapiens.GRCh38.V102.withChr.Exon.DUX4.bed",
        bam=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam",
        dexseq=rules.DEXSeq.output.DEXSeq,
        bmk=dir_out + "/{sample}/" + squence_type + "/log/{sample}.DEXSeq.bmk"
    output:
        DEXSeq= dir_out + "/{sample}/" + squence_type + "/DEXSeq_DUX4/{sample}.DUX4patched.DEXSeq",
    log:        dir_out + "/{sample}/" + squence_type + "/log/{sample}.DEXSeqDUX4.log"
    benchmark:  dir_out + "/{sample}/" + squence_type + "/log/{sample}.DEXSeqDUX4.bmk"
    params:
        pl=""
    shell:
        '''
        HTSeqDUX4patch.pl -f {input.bed} -b {input.bam} -h {input.dexseq} -o {output.DEXSeq}
        '''

rule DEXSeq_geneAggr:
    input:
        bam=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam",
        bmk=dir_out + "/{sample}/" + squence_type + "/log/{sample}_markDup.bmk",
    output:
        DEXSeq=dir_out + "/{sample}/" + squence_type + "/DEXSeq/{sample}.DEXSeq_geneAggr",
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.DEXSeq_geneAggr.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.DEXSeq_geneAggr.bmk"
    params:
        dexseq="/home/zgu_labs/bin/dexseq_count.py",
        gff=config["dexseq_gff"],
    shell:
        '''
        python {params.dexseq} -p yes -r pos -s no -a 10 -f bam {params.gff} {input.bam} {output.DEXSeq} 2> {log}
        '''

rule featureCounts:
    input:
        gtf=config["ref_gtf"],
        ref=config["ref_gatk"],
        bam=rules.samtools_sort.output.bam,
        bmk=rules.samtools_sort.benchmark
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
        featureCounts -p -F GTF -g gene_id -t exon -T {threads} -G {input.ref} -Q 10 --tmpDir {params.dir} \
        -a {input.gtf} -o {output.tab} {input.bam} &> {log}
        '''
