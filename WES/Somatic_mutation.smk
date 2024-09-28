rule Mutect2_tumor_only_mode:
    input:
        ref_fa=config["ref_gatk"],
        bam=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam",
        bai=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam.bai"
    output:
        vcf=dir_out + "/{sample}/" + squence_type + "/mutect2/{sample}.vcf.gz",
        f1r2=dir_out + "/{sample}/" + squence_type + "/mutect2/{sample}_f1r2.tar.gz",
        mutect2_stats=dir_out + "/{sample}/" + squence_type + "/mutect2/{sample}.vcf.gz.stats"
    log:
        dir_out+"/{sample}/"+squence_type+"/log/{sample}_mutect2.log"
    benchmark:
        dir_out+"/{sample}/"+squence_type+"/log/{sample}_mutect2.bmk"
    params:
        gatk=config["gatk"],
        germline_resource=config["mutect2_germline_resource"],
        pon=config["mutect2_pon"]
    shell:
        '''
        {params.gatk} Mutect2  \
        -R {input.ref_fa}  \
        -I {input.bam} \
        --germline-resource {params.germline_resource} \
        --panel-of-normals {params.pon}  \
        --max-mnp-distance 0 \
        --f1r2-tar-gz {output.f1r2} \
        -O {output.vcf} &> {log}
        '''

rule CalculateContamination:
    input:
        pileup_table=dir_out + "/{sample}/" + squence_type + "/getpileup/{sample}_pileups.table"
    output:
        seg_table=dir_out + "/{sample}/" + squence_type + "/variant_filter/{sample}_segments.table",
        contam_table=dir_out + "/{sample}/" + squence_type + "/variant_filter/{sample}_contamination.table"
    log:
        dir_out+"/{sample}/"+squence_type+"/log/{sample}_contamination.log"
    benchmark:
        dir_out+"/{sample}/"+squence_type+"/log/{sample}_contamination.bmk"
    params:
        gatk=config["gatk"]
    shell:
        '''
        {params.gatk} CalculateContamination \
        -I {input.pileup_table} \
        -tumor-segmentation {output.seg_table} \
        -O {output.contam_table} &> {log}
        '''

rule LearnReadOrientation:
    input:
        f1r2=rules.Mutect2_tumor_only_mode.output.f1r2
    output:
        readorientation=dir_out + "/{sample}/" + squence_type + "/variant_filter/{sample}_readorientation.tar.gz"
    log:
        dir_out+"/{sample}/"+squence_type+"/log/{sample}_readorientation.log"
    benchmark:
        dir_out+"/{sample}/"+squence_type+"/log/{sample}_learnreadorient.bmk"
    params:
        gatk=config["gatk"]
    shell:
        '''
        {params.gatk} LearnReadOrientationModel -I {input.f1r2} -O {output.readorientation} &> {log}
        '''

rule Filtermutectcall:
    input:
        vcf=rules.Mutect2_tumor_only_mode.output.vcf,
        ref_fa=config["ref_gatk"],
        seg_table=rules.CalculateContamination.output.seg_table,
        contam_table=rules.CalculateContamination.output.contam_table,
        readorientation=rules.LearnReadOrientation.output.readorientation,
        mutect2_stats=rules.Mutect2_tumor_only_mode.output.mutect2_stats
    output:
        vcf=dir_out + "/{sample}/" + squence_type + "/variant_filter/{sample}_filtered.vcf"
    log:
        dir_out+"/{sample}/"+squence_type+"/log/{sample}_filtermutect.log"
    benchmark:
        dir_out+"/{sample}/"+squence_type+"/log/{sample}_filtermutect.bmk"
    params:
        gatk=config["gatk"]
    shell:
        '''
        {params.gatk} FilterMutectCalls \
        -V {input.vcf} \
        -R {input.ref_fa} \
        --contamination-table {input.contam_table} \
        --tumor-segmentation {input.seg_table} \
        --stats {input.mutect2_stats} \
        --ob-priors {input.readorientation} \
        -O {output.vcf} &> {log}
        '''

rule filtersuperdup:
    input:
        vcf=rules.Filtermutectcall.output.vcf
    output:
        vcf = dir_out + "/{sample}/" + squence_type + "/variant_filter/{sample}_superdup.vcf"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}_superdup.bmk"
    params:
        superdup=config["superdup"]
    shell:
        '''
        bedtools intersect -header -v -a {input.vcf} -b {params.superdup} > {output.vcf}
        '''

rule strand_filter:
    input:
        vcf=rules.filtersuperdup.output.vcf
    output:
        vcf=dir_out + "/{sample}/" + squence_type + "/variant_filter/{sample}_strand_filter.vcf"
    log:
        dir_out+"/{sample}/"+squence_type+"/log/{sample}_filtermutation.log"
    benchmark:
        dir_out+"/{sample}/"+squence_type+"/log/{sample}_filtermutation.bmk"
    params:
        fltmut=config["filtmutation"]
    shell:
        '''
        {params.fltmut} -i {input.vcf} --min-coverage 10 --min-reads2 4 --strand-filter 1 --min-var-freq 0.03 --output-file {output.vcf} &> {log}
        '''

rule filterpass:
    input:
        vcf=rules.strand_filter.output.vcf,
        ref_fa=config["ref_gatk"]
    output:
        vcf = dir_out + "/{sample}/" + squence_type + "/variant_filter/{sample}_pass.vcf"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}_pass.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}_pass.bmk"
    params:
        gatk=config["gatk"]
    shell:
        '''
        {params.gatk} SelectVariants --exclude-filtered -R {input.ref_fa} -V {input.vcf} -O {output.vcf} &> {log}
        '''

rule Mutect2_split:
    input:
        vcf=rules.filterpass.output.vcf
    output:
        snv=    dir_out + "/{sample}/" + squence_type + "/mutect2_split/{sample}.snv.vcf",
        indel=  dir_out + "/{sample}/" + squence_type + "/mutect2_split/{sample}.indel.vcf",
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.Mutect2_split.bmk"
    shell:
        '''
        scripts/vcfSplit.sh -v {input.vcf} -s {output.snv} -i {output.indel}
        '''

rule get_bed:
    input:
        snv_in=    rules.Mutect2_split.output.snv,
        indel_in=  rules.Mutect2_split.output.indel
    output:
        snv_bed=dir_out + "/{sample}/" + squence_type + "/mutect2_split/{sample}.snv.bed",
        indel_bed=dir_out + "/{sample}/" + squence_type + "/mutect2_split/{sample}.indel.bed"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.get_bed.bmk"
    shell:
        '''
        bcftools query -f '%CHROM\\t%POS0\\t%END\\n' {input.snv_in} > {output.snv_bed}
        bcftools query -f '%CHROM\\t%POS0\\t%END\\n' {input.indel_in} > {output.indel_bed}
        '''

# rule get_PON:
#     input:
#         files_snv=expand(dir_out + "/{sample}/" + squence_type + "/mutect2_split/{sample}.snv.bed",sample=samplelist_PON),
#         files_indel=expand(dir_out + "/{sample}/" + squence_type + "/mutect2_split/{sample}.indel.bed",sample=samplelist_PON)
#     output:
#         bed_snv_f=dir_out + "/stats/" + squence_type + "/get_PON/PON_snv.forFiltering.bed",
#         bed_indel_f=dir_out + "/stats/" + squence_type + "/get_PON/PON_indel.forFiltering.bed"
#     benchmark:
#         dir_out + "/stats/" + squence_type + "/log/get_PON.bmk"
#     params:
#         bed_snv=dir_out + "/stats/" + squence_type + "/get_PON/PON_snv.bed",
#         bed_indel=dir_out + "/stats/" + squence_type + "/get_PON/PON_indel.bed"
#     shell:
#         '''
#         cat {input.files_snv}| sort -k1,1 -k2,2n|bedtools merge -i - -c 1 -o count >  {params.bed_snv}
#         cat {params.bed_snv} |awk '$4 >=2' > {output.bed_snv_f}
#
#         cat {input.files_indel}| sort -k1,1 -k2,2n|bedtools merge -i - -c 1 -o count >  {params.bed_indel}
#         cat {params.bed_indel} |awk '$4 >=2' > {output.bed_indel_f}
#         '''

rule filter_PON:
    input:
        snv_in = rules.Mutect2_split.output.snv,
        indel_in = rules.Mutect2_split.output.indel,
        PON_snv = dir_out + "/stats/EXOME/get_PON/PON_snv.forFiltering.bed",
        PON_indel = dir_out + "/stats/EXOME/get_PON/PON_indel.forFiltering.bed"
    output:
        snv_out = dir_out + "/{sample}/" + squence_type + "/filter_PON/{sample}.snv.vcf",
        indel_out = dir_out + "/{sample}/" + squence_type + "/filter_PON/{sample}.indel.vcf"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.filter_PON.bmk"
    shell:
        '''
        bedtools intersect -v -a {input.snv_in} -b {input.PON_snv} > {output.snv_out}
        bedtools intersect -v -a {input.indel_in} -b {input.PON_indel} > {output.indel_out}
        '''

rule vcf_QC_Mutect2:
    input:
        snv_in=    rules.filter_PON.output.snv_out,
        indel_in=   rules.filter_PON.output.indel_out
    output:
        snv_out = dir_out + "/{sample}/" + squence_type + "/vcf_QC_Mutect2/{sample}.snv.vcf",
        indel_out = dir_out + "/{sample}/" + squence_type + "/vcf_QC_Mutect2/{sample}.indel.vcf",
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.vcf_QC_Mutect2.bmk"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.vcf_QC_Mutect2.log"
    shell:
        '''
        scripts/1.vcf_QC_Mutect2.pl -i {input.snv_in} -n 10 -D 4 -h 0.05 -o {output.snv_out} &> {log}
        scripts/1.vcf_QC_Mutect2.pl -i {input.indel_in} -n 10 -D 4 -h 0.05 -o {output.indel_out} &>> {log}
        '''

rule DvsG_filter_noMatch:
    input:
        snv_in = rules.vcf_QC_Mutect2.output.snv_out,
        indel_in = rules.vcf_QC_Mutect2.output.indel_out,
        snv_control=dir_out + "/stats/EXOME/get_PON/PON_control.snv.vcf",
        indel_control=dir_out + "/stats/EXOME/get_PON/PON_control.indel.vcf"
    output:
        snv_out=dir_out + "/{sample}/" + squence_type + "/DvsG_filter_noMatch/{sample}.snv.vcf",
        indel_out=dir_out + "/{sample}/" + squence_type + "/DvsG_filter_noMatch/{sample}.indel.vcf"
    benchmark:  dir_out + "/{sample}/" + squence_type + "/log/{sample}.DvsG_filter_noMatch.bmk"
    log:        dir_out + "/{sample}/" + squence_type + "/log/{sample}.DvsG_filter_noMatch.log"
    shell:
        '''
        2.1.DvsG_filter.pl -a {input.snv_in} -c {input.snv_control} -o {output.snv_out} &> {log}
        2.1.DvsG_filter.pl -a {input.indel_in} -c {input.indel_control} -o {output.indel_out} &>> {log}
        '''

rule DvsG_cns_noMatch:
    input:
        snv_in=rules.DvsG_filter_noMatch.output.snv_out,
        indel_in=rules.DvsG_filter_noMatch.output.indel_out
    output:
        snv_out = dir_out + "/{sample}/" + squence_type + "/DvsG_cns_noMatch/{sample}.snv.vcf",
        indel_out = dir_out + "/{sample}/" + squence_type + "/DvsG_cns_noMatch/{sample}.indel.vcf",
    benchmark:  dir_out + "/{sample}/" + squence_type + "/log/{sample}.DvsG_cns_noMatch.bmk"
    log:        dir_out + "/{sample}/" + squence_type + "/log/{sample}.DvsG_cns_noMatch.log"
    params:
        fa=config["ref_gatk"],
        bam_control=dir_out + "/DLBCL_Song_N12/EXOME/bam/DLBCL_Song_N12.bam"
    shell:
        '''
        2.2.snv_DvsG_cns.pl    -f {params.fa} -a {input.snv_in}   -c {params.bam_control} -o {output.snv_out}   &> {log}
        2.2.indel_DvsG_cns.pl  -f {params.fa} -a {input.indel_in} -c {params.bam_control} -o {output.indel_out} &>> {log}
        '''

rule filter_public:
    input:
        snv=rules.DvsG_cns_noMatch.output.snv_out,
        indel=rules.DvsG_cns_noMatch.output.indel_out
    output:
        snv_out=dir_out + "/{sample}/" + squence_type + "/filter_public/{sample}.snv.vcf",
        indel_out=dir_out + "/{sample}/" + squence_type + "/filter_public/{sample}.indel.vcf",
    benchmark:  dir_out + "/{sample}/" + squence_type + "/log/{sample}.filter_public.bmk"
    log:        dir_out + "/{sample}/" + squence_type + "/log/{sample}.filter_public.log"
    params:
        ref_snv_HP="/home/zgu_labs/refData/NormalRNAseq/Normal.HC.GRCh38.snv",
        ref_indel_HP="/home/zgu_labs/refData/NormalRNAseq/Normal.HC.GRCh38.indel",
        N_snv_HP=3,N_indel_HP=3
    shell:
        '''
        rmgerm.pl -a {input.snv}     -n {params.N_snv_HP}    -c {params.ref_snv_HP}   -o {output.snv_out}  &> {log}      
        rmgerm.pl -a {input.indel}   -n {params.N_indel_HP}  -c {params.ref_indel_HP} -o {output.indel_out}  &>> {log}
        '''

rule filter_dbSNP_2:
    input:
        snv_in = rules.vcf_QC_Mutect2.output.snv_out,
        indel_in = rules.vcf_QC_Mutect2.output.indel_out
    output:
        snv_out = dir_out + "/{sample}/" + squence_type + "/filter_dbSNP/{sample}.snv.tsv",
        indel_out = dir_out + "/{sample}/" + squence_type + "/filter_dbSNP/{sample}.indel.tsv",
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.filter_dbSNP2.bmk"
    shell:
        '''
        scripts/3.snv_byDbsnpBB.pl -b /ref_genomes/dbSNP/human/153_hg38/dbSnp153.bb -i {input.snv_in} -o {output.snv_out} 
        scripts/3.indel_byDbsnpBB.pl -b /ref_genomes/dbSNP/human/153_hg38/dbSnp153.bb -i {input.indel_in} -o {output.indel_out} 
        '''

rule Annotation:
    input:
        snv=rules.filter_dbSNP_2.output.snv_out,
        indel=rules.filter_dbSNP_2.output.indel_out
    output:
        snv_out = dir_out + "/{sample}/" + squence_type + "/Annotation/{sample}.snv.tsv",
        indel_out = dir_out + "/{sample}/" + squence_type + "/Annotation/{sample}.indel.tsv",
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.Annotation.bmk"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.Annotation.log"
    params:
        snv=dir_out + "/{sample}/" + squence_type + "/Annotation/{sample}.snv.2.tsv",
        indel=dir_out + "/{sample}/" + squence_type + "/Annotation/{sample}.indel.2.tsv",
        N_snv=1,N_indel=1
    shell:
        '''
        scripts/AnnoSNV.pl -i {input.snv} -o {params.snv} &> {log}
        scripts/runVep.sh {params.snv} {output.snv_out} {params.N_snv} &>> {log}       
        scripts/AnnoINDEL.pl -i {input.indel} -o {params.indel}   &>> {log}
        scripts/runVep.sh {params.indel} {output.indel_out} {params.N_indel}  &>> {log}
        rm {params.snv} {params.indel}
        '''

rule keep_nonSlient_mutation:
    input:
        snv_in= rules.Annotation.output.snv_out,
        indel_in= rules.Annotation.output.indel_out
    output:
        snv_out=dir_out + "/{sample}/" + squence_type + "/keep_nonSlient_mutation/{sample}.snv.tsv",
        indel_out=dir_out + "/{sample}/" + squence_type + "/keep_nonSlient_mutation/{sample}.indel.tsv",
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.keep_nonSlient_mutation.bmk"
    shell:
        '''
        cat {input.snv_in} |awk 'NR==1 || $23~/exon|frame_shift_del|frame_shift_ins|in-frame_del|in-frame_ins|missense|nonsense|splice|intron-NS/' > {output.snv_out}
        cat {input.indel_in} |awk 'NR==1 || $23~/exon|frame_shift_del|frame_shift_ins|in-frame_del|in-frame_ins|missense|nonsense|splice|intron-NS/' > {output.indel_out}
        '''

rule filter_gnomAD:
    input:
        snv_in=rules.keep_nonSlient_mutation.output.snv_out,
        indel_in=rules.keep_nonSlient_mutation.output.indel_out
    output:
        snv_out = dir_out + "/{sample}/" + squence_type + "/filter_gnomAD/{sample}.snv.mini.tsv",
        indel_out = dir_out + "/{sample}/" + squence_type + "/filter_gnomAD/{sample}.indel.mini.tsv"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.filter_gnomAD.bmk"
    params:
        sample = "{sample}",
        snv=dir_out + "/{sample}/" + squence_type + "/filter_gnomAD/{sample}.snv.tsv",
        indel=dir_out + "/{sample}/" + squence_type + "/filter_gnomAD/{sample}.indel.tsv"
    shell:
        '''
        cat {input.snv_in} | awk 'NR==1 || $98 < 0.001' > {params.snv}
        cat {params.snv} | cut -f 1-16,21-23,25,26,29-34,44,98 |perl -lane '$,="\\t";$sample="{params.sample}";print @F,$sample' > {output.snv_out}
        cat {input.indel_in} | awk 'NR==1 || $98 < 0.001' > {params.indel}
        cat {params.indel} | cut -f 1-16,21-23,25,26,29-34,44,98 |perl -lane '$,="\\t";$sample="{params.sample}";print @F,$sample' > {output.indel_out}
        '''

rule Funcotator:
    input:
        vcf=rules.filterpass.output.vcf,
        ref_fa=config["ref_gatk"],
        data_source=config["funcotator_data_sources"]
    output:
        vcf=dir_out + "/{sample}/" + squence_type + "/funcotator/{sample}_funcotator.vcf"
    log:
        dir_out+"/{sample}/"+squence_type+"/log/{sample}_funcotator.log"
    benchmark:
        dir_out+"/{sample}/"+squence_type+"/log/{sample}_funcotator.bmk"
    params:
        gatk=config["gatk"]
    shell:
        '''
        {params.gatk} Funcotator \
        -R {input.ref_fa} \
        --variant {input.vcf} \
        --output-file-format VCF \
        --data-sources-path {input.data_source} \
        --ref-version hg38 \
        -O {output.vcf} &> {log}
        '''

rule vep:
    input:
        vcf=rules.Funcotator.output.vcf
    output:
        vcf=dir_out + "/{sample}/" + squence_type + "/funcotator/{sample}_vep.vcf",
        filter=dir_out + "/{sample}/" + squence_type + "/funcotator/{sample}_vep_filtered.vcf"
    log:
        dir_out+"/{sample}/"+squence_type+"/log/{sample}_vep.log"
    benchmark:
        dir_out+"/{sample}/"+squence_type+"/log/{sample}_vep.bmk"
    params:
        cache=config["vep_cache"],
        stats=dir_out + "/{sample}/" + squence_type + "/funcotator/{sample}_vep.html"
    shell:
        '''
        PATH="/home/zgu_labs/anaconda3/envs/vep/bin:$PATH"
        vep -cache --species homo_sapiens -a GRCh38 --force --pick -e --vcf --filter_common \
        --fork 1 --dir_cache {params.cache} -i {input.vcf} \
        --stats_file {params.stats} -o {output.vcf} &> {log}
        PATH="/home/zgu_labs/anaconda3/envs/vep/bin:$PATH" filter_vep -i {output.vcf} --force_overwrite --format vcf -filter "(MAX_AF < 0.01 or not MAX_AF) and (IMPACT != MODIFIER)" -o {output.filter}
        '''

rule variantotable:
    input:
        vcf=rules.vep.output.filter
    output:
        table=dir_out + "/{sample}/" + squence_type + "/funcotator/{sample}_SNV.table"
    log:
        dir_out+"/{sample}/"+squence_type+"/log/{sample}_fltmutotable.log"
    benchmark:
        dir_out+"/{sample}/"+squence_type+"/log/{sample}_fltmutotable.bmk"
    params:
        gatk=config["gatk"]
    shell:
        '''
        {params.gatk} VariantsToTable \
        -V {input.vcf} \
        -F CHROM -F POS -F REF -F ALT -F TYPE -F FUNCOTATION -F CSQ -F POPAF \
        -GF GT -GF DP -GF AD -GF AF -GF F1R2 -GF F2R1 -GF SB \
        -O {output.table} &> {log}
        '''

rule SOBDetector:
    input:
        tsv=rules.variantotable.output.table,
        bam=dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam"
    output:
        tsv=dir_out + "/{sample}/" + squence_type + "/funcotator/{sample}_SOBD.tsv"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}_SOBD.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}_SOBD.bmk"
    shell:
        '''
        java -jar scripts/SOBDetector_v1.0.4.jar --input-type Table --input-variants {input.tsv} --input-bam {input.bam} --output-variants {output.tsv} &> {log}
        '''

rule artifact_filter:
    input:
        tsv=rules.SOBDetector.output.tsv
    output:
        tsv = dir_out + "/{sample}/" + squence_type + "/funcotator/{sample}_artifact_filter.tsv"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}_artifact_filter.bmk"
    shell:
        '''
        cat {input.tsv}|grep -v "FUNCOTATION"|perl -lane '
        @funcotator=split("[=[|]",$F[5]);
        @funcotator=map {{ $_ eq "" ? "NA" : $_ }} @funcotator;
        @vep=split("[|]",$F[6]);
        @vep=map {{ $_ eq "" ? "NA" : $_ }} @vep;
        @nonsilent_mutation =("MISSENSE","NONSENSE","NONSTOP","SPLICE_SITE","IN_FRAME_DEL","IN_FRAME_INS","FRAME_SHIFT_INS",
        "FRAME_SHIFT_DEL","START_CODON_SNP","START_CODON_INS","START_CODON_DEL","DE_NOVO_START_IN_FRAME",
        "DE_NOVO_START_OUT_FRAME");
        $condition=grep(/$funcotator[6]/i, @nonsilent_mutation);
        $nonsilent_var = $condition>0 ? "Nonsilent" : "Silent";
        print join "\\t", $F[0],$F[1],$F[2],$F[3],$F[4],$funcotator[1],$funcotator[6],$funcotator[14],$funcotator[19],$F[9],$F[10],$F[11],$vep[2],$vep[21],$vep[36],$vep[37]
        if $nonsilent_var eq "Nonsilent" & $F[21] < 0.7' > {output.tsv}
        '''

rule perl_filter_Mutect2noNormal:
    input:
        tsv=rules.artifact_filter.output.tsv
    output:
        tsv=dir_out + "/{sample}/" + squence_type + "/perl_filter_Mutect2noNormal/{sample}.perl_filter_Mutect2noNormal.tsv"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.perl_filter_Mutect2noNormal.log"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.perl_filter_Mutect2noNormal.bmk"
    params:
        perl_params='''-F"\\t" -lane '$,="\\t"; $chr=$F[0]; $start=$F[1]; $refBase=$F[2]; $varBase=$F[3]; $depth=$F[9]; $MAF=$F[11]; $mutType=$F[4]; 
        $gene=$F[5]; $mutClass=$F[6]; $strand=$F[7]; $aac=$F[8]; $Impact=$F[12]; $SIFT=$F[14]; $Polyphen=$F[15];
        @ref=split("",$refBase); $ref_len=$#ref; $end=$F[1]+$ref_len; @AD_D=split(",",$F[10]); $varDepth=@AD_D[1];
        print $chr,$start,$refBase,$varBase,$end,$gene,$strand,$depth,$varDepth,$MAF,$mutType,$mutClass,$aac,$Impact,$SIFT,$Polyphen' '''
    shell:
        '''
        cat <(echo 'chr,start,refBase,varBase,end,Gene,strand,depth,varDepth,MAF,mutType,mutClass,AAC,Impact,SIFT,Polyphen'|sed 's/,/\\t/g') \
        <(cat {input.tsv} |perl {params.perl_params}|awk '$6>=10') > {output.tsv}
        '''

rule filter_dbSNP:
    input:
        tsv=rules.perl_filter_Mutect2noNormal.output.tsv
    output:
        snv=dir_out + "/{sample}/" + squence_type + "/perl_filter_Mutect2noNormal/{sample}.dbsnp_filter.snv.tsv",
        indel=dir_out + "/{sample}/" + squence_type + "/perl_filter_Mutect2noNormal/{sample}.dbsnp_filter.indel.tsv"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.filter_dbSNP.bmk"
    params:
        snv_tmp=dir_out + "/{sample}/" + squence_type + "/perl_filter_Mutect2noNormal/{sample}.snvtmp.tsv",
        indel_tmp=dir_out + "/{sample}/" + squence_type + "/perl_filter_Mutect2noNormal/{sample}.indeltmp.tsv"
    shell:
        '''
        cat {input.tsv} |awk -v OFS="\\t" '$11=="SNP" || NR==1' > {params.snv_tmp}
        cat {input.tsv} |awk -v OFS="\\t" '$11=="INDEL" || NR==1'  > {params.indel_tmp}

        scripts/3.snv_byDbsnpBB.pl -b /ref_genomes/dbSNP/human/153_hg38/dbSnp153.bb -i {params.snv_tmp} -o {output.snv} 
        scripts/3.indel_byDbsnpBB.pl -b /ref_genomes/dbSNP/human/153_hg38/dbSnp153.bb -i {params.indel_tmp} -o {output.indel} 
        
        rm {params.snv_tmp} {params.indel_tmp}
        '''

rule search_cosmic:
    input:
        snv=rules.filter_dbSNP.output.snv,
        indel=rules.filter_dbSNP.output.indel,
        strelka_snv= dir_out + "/{sample}/" + squence_type + "/Strelka/results/variants/filtered.somatic.snvs.txt",
        strelka_indel = dir_out + "/{sample}/" + squence_type + "/Strelka/results/variants/filtered.somatic.indel.txt",
        filter_snv = rules.filter_public.output.snv_out,
        filter_indel = rules.filter_public.output.indel_out
    output:
        all = dir_out + "/{sample}/" + squence_type + "/perl_filter_Mutect2noNormal/{sample}.allmut.tsv"
    benchmark:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.search_cosmic.bmk"
    log:
        dir_out + "/{sample}/" + squence_type + "/log/{sample}.search_cosmic.log"
    params:
        sample = "{sample}"
    shell:
        '''
        Rscript scripts/search_cosmic.R {input.snv} {input.indel} {params.sample} {input.strelka_snv} {input.strelka_indel} {input.filter_snv} {input.filter_indel} &> {log}
        '''

rule done_mutation:
    input:
        dbsnp_mutect=rules.search_cosmic.benchmark,
        #dbsnp_gnomad=rules.filter_gnomAD.benchmark
    output:
        done = dir_out + "/{sample}/" + squence_type + "/log/doneMutation.txt"
    shell:
        '''
        touch {output.done}
        '''

# rulemicrosec:
# input:
#     maf = rules.search_cosmic.output.all, \
#     bam = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam", \
#     bai = dir_out + "/{sample}/" + squence_type + "/bam/{sample}.bam.bai"
# output:
#     convert=dir_out + "/{sample}/" + squence_type + "/perl_filter_Mutect2noNormal/{sample}_Mutations_converted.xlsx"
# params:
#     dir=dir_out + "/FFPEfilter/{sample}",
#     maf=dir_out + "/FFPEfilter/{sample}/{sample}_Mutations_converted.xlsx",
#     bam=dir_out + "/FFPEfilter/{sample}/{sample}.T.bam",
#     bai=dir_out + "/FFPEfilter/{sample}/{sample}.T.bam.bai",
#     temp=dir_out + "/{sample}/" + squence_type + "/perl_filter_Mutect2noNormal/{sample}_Mutations.xlsx"
# benchmark:
#     dir_out + "/{sample}/" + squence_type + "/log/{sample}.microsec.bmk"
# log:
#     dir_out + "/{sample}/" + squence_type + "/log/{sample}.microsec.log"
# shell:
#     '''
#     if [ -d "{params.dir}" ]; then rm -rf {params.dir};fi
#     mkdir {params.dir}
#     ln -s {input.bam} {params.bam}
#     ln -s {input.bai} {params.bai}
#     Rscript ffpefiltering/convert_docker.R {input.maf} {params.temp} {output.convert}
#     ln -s {output.convert} {params.maf}
#     rm {params.temp}
#     '''

# rule ffpefilter:
#     input:
#         maf = expand(dir_out + "/{sample}/" + squence_type + "/log/{sample}.microsec.bmk",sample=samplelist)
#     benchmark:
#         dir_out + "/FFPEfilter/ffpefilter.bmk"
#     log:
#         dir_out + "/FFPEfilter/ffpefilter.log"
#     params:
#         dir=dir_out + "/FFPEfilter",
#         samples = dir_out + "/FFPEfilter/sample_info.tsv"
#     shell:
#         '''
#         Rscript scripts/MicroSEC.R {params.dir} {params.samples} Y &> {log}
#         '''

# rule search_dbsnp_gnomad:
#     input:
#         file_out_variants_sample=rules.artifact_filter.output.tsv,
#         dbsnp_common="/ref_genomes/dbSNP/human/153_hg38/dbSnp153Common.bb",
#         dbsnp_all="/ref_genomes/dbSNP/human/153_hg38/dbSnp153.bb",
#         gnomad_genome="/ref_genomes/gnomad/gnomad.genomes.r3.0.sites.vcf.bgz",
#         gnomad_exome= "/ref_genomes/gnomad/GRCh38/Exomes/r2.1.1/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz"
#     output:
#         tsv=dir_out+"/{sample}/" + squence_type + "/funcotator/{sample}.dbsnp_gnomad.tsv"
#     log:
#         dir_out+"/{sample}/"+squence_type+"/log/{sample}.dbsnp_gnomad.log"
#     benchmark:
#         dir_out+"/{sample}/"+squence_type+"/log/{sample}.dbsnp_gnomad.bmk"
#     params:
#         bigBedToBed="/home/zgu_labs/bin/software/bigbed/bigBedToBed",
#         bcftools="bcftools"
#     shell:
#         '''
#         Rscript scripts/variants_search_dbsnp_gnomad.R {input.file_out_variants_sample} {params.bigBedToBed} {params.bcftools} {input.dbsnp_common} {input.dbsnp_all} {input.gnomad_genome} {input.gnomad_exome} {output.tsv} 2>{log}
#         for i in {output.tsv}; do awk '{{var=FILENAME; split(var,a,/[\/.]/); print a[8]"\t"$0}}' $i > $i.bk; mv $i.bk $i; done
#         '''
#
# rule Mutect2_annotation_snv:
#     input:
#         snv= rules.Mutect2_split.output.snv
#     output:
#         snv=   dir_out + "/{sample}/" + squence_type + "/variant_filter/{sample}.Mutect2.snv.tsv",
#     benchmark:
#         dir_out + "/{sample}/" + squence_type + "/log/{sample}.Mutect2_annotation_snv.bmk"
#     params:
#         filtermutect = config["filtermutect"],
#         germSNP_snv=    "/home/zgu_labs/refData/NormalRNAseq/Normal.HC.GRCh38.snv",
#         snv1 =          dir_out + "/{sample}/" + squence_type + "/variant_filter/{sample}.1.snv",
#         snv2 =          dir_out + "/{sample}/" + squence_type + "/variant_filter/{sample}.2.snv",
#         snv4 =          dir_out + "/{sample}/" + squence_type + "/variant_filter/{sample}.4.snv",
#         N_qc=4,H_qc=0.01,N_depth=10,N_filter_Germ=3,N_vep=1
#     shell:
#         '''
#         PATH="/home/zgu_labs/anaconda3/envs/vep/bin:$PATH"
#         {params.filtermutect} -i {input.snv} -n {params.N_qc} -D {params.N_depth} -h {params.H_qc} -o {params.snv1}
#         rmgerm.pl -a {params.snv1} -n {params.N_filter_Germ} -c {params.germSNP_snv} -o {params.snv2}
#         AnnoSNV.pl -i {params.snv2} -o {params.snv4}
#         runVep.sh {params.snv4} {output.snv} {params.N_vep}
#         rm {params.snv1} {params.snv4}
#         '''
#
# rule Mutect2_annotation_indel:
#     input:
#         indel= rules.Mutect2_split.output.indel
#     output:
#         indel= dir_out + "/{sample}/" + squence_type + "/variant_filter/{sample}.Mutect2.indel.tsv",
#     benchmark:
#         dir_out + "/{sample}/" + squence_type + "/log/{sample}.Mutect2_annotation_indel.bmk"
#     params:
#         filtermutect = config["filtermutect"],
#         germSNP_indel=  "/home/zgu_labs/refData/NormalRNAseq/Normal.HC.GRCh38.indel",
#         indel1=dir_out + "/{sample}/" + squence_type + "/variant_filter/{sample}.1.indel",
#         indel2=dir_out + "/{sample}/" + squence_type + "/variant_filter/{sample}.2.indel",
#         indel4=dir_out + "/{sample}/" + squence_type + "/variant_filter/{sample}.4.indel",
#         N_qc=4,H_qc=0.01,N_depth=10,N_filter_Germ=3,N_vep=1
#     shell:
#         '''
#         {params.filtermutect} -i {input.indel} -n {params.N_qc} -D {params.N_depth} -h {params.H_qc} -o {params.indel1}
#         rmgerm.pl -a {params.indel1} -n {params.N_filter_Germ} -c {params.germSNP_indel} -o {params.indel2}
#         AnnoINDEL.pl -i {params.indel2} -o {params.indel4}
#         runVep.sh {params.indel4} {output.indel} {params.N_vep}
#         rm {params.indel1} {params.indel4}
#         '''
#
# rule done_mutect:
#     input:
#         Mutect2_snv = rules.Mutect2_annotation_snv.output.snv,
#         Mutect2_indel = rules.Mutect2_annotation_indel.output.indel
#     benchmark:
#         dir_out + "/{sample}/" + squence_type + "/log/{sample}.annotation_done.bmk"
#     params:
#         tsv = dir_out + "/{sample}/" + squence_type + "/Mutation/{sample}.mutect2_anno.tsv"
#     output:
#         tsv = dir_out + "/{sample}/" + squence_type + "/Mutation/{sample}.mutect2_mutect2_filter.tsv"
#     shell:
#         '''
#         cat {input.Mutect2_snv} <(cat {input.Mutect2_indel}|awk -v OFS="\\t" 'NR>1') > {params.tsv}
#         awk -F '\t' 'NR == 1 {{print;next}} NR>2 && $9>=4 && $11>=0.01 && ($16=="missense"||$16=="nonsense"||$16=="frame_shift_del"||$16=="frame_shift_ins"||$16=="in-frame_del"||$16=="in-frame_ins"||$16=="splice"||$16=="startloss"||$16=="stoploss")' {params.tsv} > {output.tsv}
#         '''
#
# rule perl_combine_Mutect2noNormal:
#     input:
#         funcotator_mut=expand(dir_out+"/{sample}/"+squence_type+"/log/{sample}.dbsnp_gnomad.bmk",sample=samplelist[0]),
#         perlfilter_mut=expand(dir_out + "/{sample}/" + squence_type + "/log/{sample}.search_dbsnp_human_Mutect2noNormal.bmk",sample=samplelist[0])
#     output:
#         tsv=dir_out + "/id/All_mutation_" + squence_type + ".tsv",
#         perlanno_tsv=dir_out + "/id/Perlanno_Allmutation_" + squence_type + ".tsv"
#     benchmark:
#         dir_out + "/id/All_mutation_" + squence_type + ".bmk"
#     params:
#         file_prefix=dir_out+"/*/" + squence_type + "/funcotator/*.dbsnp_gnomad.tsv",
#         file_header=dir_out+"/Jejoye/EXOME/funcotator/Jejoye.dbsnp_gnomad.tsv",
#         perlanno_file_prefix=dir_out + "/*/" + squence_type + "/search_dbsnp_human_Mutect2noNormal/*.search_dbsnp.tsv",
#         perlanno_file_header=dir_out+"/Jejoye/EXOME/search_dbsnp_human_Mutect2noNormal/Jejoye.search_dbsnp.tsv"
#     shell:
#         '''
#         cat <(head -n1 {params.file_header}|perl -lane  ' $,="\t";print "Sample",@F[1..$#F]') \
#         <(ls {params.file_prefix}|xargs -i echo "cat {{}}|grep -v varDepth"|bash) > {output.tsv}
#
#         cat <(head -n1 {params.perlanno_file_header}|perl -lane  ' $,="\t";print "Sample",@F[1..$#F]') \
#         <(ls {params.perlanno_file_prefix}|xargs -i echo "cat {{}}|grep -v varDepth"|bash) > {output.perlanno_tsv}
#         '''
#
# rule search_gnomAD_Mutect2noNormal:
#     input:
#         tsv_snv = rules.filter_dbSNP.output.snv,
#         tsv_indel = rules.filter_dbSNP.output.indel
#     output:
#         tsv_snv = dir_out + "/{sample}/" + squence_type + "/search_gnomAD_Mutect2noNormal/{sample}.search_gnomAD.snv.tsv",
#         tsv_indel = dir_out + "/{sample}/" + squence_type + "/search_gnomAD_Mutect2noNormal/{sample}.search_gnomAD.indel.tsv"
#     log:
#         dir_out + "/{sample}/" + squence_type + "/log/{sample}.search_gnomAD_Mutect2noNormal.log"
#     benchmark:
#         dir_out + "/{sample}/" + squence_type + "/log/{sample}.search_gnomAD_Mutect2noNormal.bmk"
#     params:
#         script="/coh_labs/jochan/pipeline/WES/scripts/search_gnomad.R",
#         vcf_dir_genome="/ref_genomes/gnomad/GRCh38/Genomes/v3.1/",
#         vcf_gnomADExome="/ref_genomes/gnomad/GRCh38/Exomes/r2.1.1/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz",
#         bcftools="bcftools"
#     shell:
#         '''
#         Rscript {params.script} {input.tsv_snv} {input.tsv_indel} {params.vcf_dir_genome} {params.vcf_gnomADExome} {params.bcftools} {output.tsv_snv} {output.tsv_indel} &> {log}
#         '''
#
# rule search_dbsnp_human_Mutect2noNormal:
#     input:
#         tsv_snv = rules.search_gnomAD_Mutect2noNormal.output.tsv_snv,
#         tsv_indel = rules.search_gnomAD_Mutect2noNormal.output.tsv_indel
#     output:
#         tsv_snv = dir_out + "/{sample}/" + squence_type + "/search_dbsnp_human_Mutect2noNormal/{sample}.search_dbsnp.snv.tsv",
#         tsv_indel = dir_out + "/{sample}/" + squence_type + "/search_dbsnp_human_Mutect2noNormal/{sample}.search_dbsnp.indel.tsv",
#         tsv = dir_out + "/{sample}/" + squence_type + "/search_dbsnp_human_Mutect2noNormal/{sample}.search_dbsnp.tsv"
#     log:
#         dir_out + "/{sample}/" + squence_type + "/log/{sample}.search_dbsnp_human_Mutect2noNormal.log"
#     benchmark:
#         dir_out + "/{sample}/" + squence_type + "/log/{sample}.search_dbsnp_human_Mutect2noNormal.bmk"
#     params:
#         script="/coh_labs/jochan/pipeline/WES/scripts/search_dbsnp_human.R",
#         bigBedToBed="/home/zgu_labs/anaconda3/bin/bigBedToBed",
#         df_dbsnp="/ref_genomes/dbSNP/human/153_hg38/dbSnp153.bb"
#     shell:
#         '''
#         Rscript {params.script} {input.tsv_snv} {input.tsv_indel} {params.df_dbsnp} {params.bigBedToBed} {output.tsv_snv} {output.tsv_indel} &> {log}
#
#         cat {output.tsv_snv} <(cat {output.tsv_indel}|awk -v OFS="\\t" 'NR>1') > {output.tsv}
#
#         for i in {output.tsv}; do awk '{{var=FILENAME; split(var,a,/[\/.]/); print a[9]"\t"$0}}' $i > $i.bk; mv $i.bk $i; done
#         '''