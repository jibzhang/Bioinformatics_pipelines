rule buildkmer:
    input:
        ref=config["ref_bwa"]
    output: "/coh_labs/jochan/reference_genome/hg38/hg38_pathseqkmer.hss.bfi"
    log:        dir_out + "/log/buildkmer.log"
    benchmark:  dir_out + "/log/buildkmer.bmk"
    shell:
        '''
        gatk.422 --java-options "-Xmx60G" PathSeqBuildKmers --reference {input.ref} --output {output} --bloom-false-positive-probability 0.001 --kmer-mask 16 --kmer-size 31 &> {log}
        '''

rule creatimage:
    input:
        human_ref = config["ref_bwa"],
        ebv_ref = config["EBV_genome"]
    output:
        human_ima = "/coh_labs/jochan/reference_genome/hg38/hg38.fasta.img",
        ebv_ima = "/coh_labs/jochan/reference_genome/EBV/ebv.fasta.img"
    log: dir_out + "/log/creatimage.log"
    benchmark:  dir_out + "/log/creatimage.bmk"
    shell:
        '''
        gatk.422 BwaMemIndexImageCreator -I {input.human_ref} -O {output.human_ima} &> {log}
        gatk.422 BwaMemIndexImageCreator -I {input.ebv_ref} -O {output.ebv_ima}
        '''

rule Taxonomy:
    input:
        ebv_ref = config["EBV_genome"],
        ebv_catalog = config["EBV_catalog"],
        ebv_taxdump = config["EBV_taxdump"]
    output:
        "/coh_labs/jochan/reference_genome/EBV/ebv.taxonomy.db"
    log: dir_out + "/log/Taxonomy.log"
    benchmark:  dir_out + "/log/Taxonomy.bmk"
    shell:
        '''
        gatk.422 PathSeqBuildReferenceTaxonomy --reference {input.ebv_ref} --output {output} --refseq-catalog {input.ebv_catalog} --tax-dump {input.ebv_taxdump} --min-non-virus-contig-length 2000 &> {log}
        '''

rule pathseq:
    input:
        bam =  dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam",
        human_image = rules.creatimage.output.human_ima,
        ebv_image = rules.creatimage.output.ebv_ima,
        ebv_dict = config["EBV_dict"],
        kmer = rules.buildkmer.output,
        ebv_ref= config["EBV_genome"],
        taxonomy = rules.Taxonomy.output
    output:
        bam = dir_out + "/{sample}/" + squence_type + "/pathoseq/{sample}.bam",
        score = dir_out + "/{sample}/" + squence_type + "/pathoseq/{sample}.score.txt",
        filter_metrics = dir_out + "/{sample}/" + squence_type + "/pathoseq/{sample}.filter_metrics.txt",
        score_metrics = dir_out + "/{sample}/" + squence_type + "/pathoseq/{sample}.score_metrics.txt"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}_pathoseq.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}_pathoseq.bmk"
    shell:
        '''
        gatk PathSeqPipelineSpark --input {input.bam} --kmer-file {input.kmer} --filter-bwa-image {input.human_image} --microbe-bwa-image {input.ebv_image} \
        --microbe-dict {input.ebv_dict} --taxonomy-file {input.taxonomy} --is-host-aligned true --min-clipped-read-length 50 --min-score-identity 0.90 --identity-margin 0.02 \
        --scores-output {output.score} --output {output.bam} --filter-metrics {output.filter_metrics} --score-metrics {output.score_metrics} &> {log}
        '''
