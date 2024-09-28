include: "common_rule.smk"

sample_1=['COH000337_M40D40WT','COH000339_M40D90WT','COH000344_F25D40WT',     'COH000345_F25D90WT',     'COH000337_M40D40WT','COH000344_F25D40WT','COH000338_M40D40KO','COH000346_F25D40KO_eRNP','COH000338_M40D40KO', 'COH000340_M40D90KO']
sample_2=['COH000338_M40D40KO','COH000340_M40D90KO','COH000346_F25D40KO_eRNP','COH000347_F25D90KO_eRNP','COH000339_M40D90WT','COH000345_F25D90WT','COH000340_M40D90KO','COH000347_F25D90KO_eRNP','COH000341_M40D230KO','COH000341_M40D230KO']
name_1=  ['M40D40WT',          'M40D90WT',          'F25D40WT',               'F25D90WT',               'M40D40WT',          'F25D40WT',          'M40D40KO',           'F25D40KO',              'M40D40KO',            'M40D90KO']
name_2=  ['M40D40KO',          'M40D90KO',          'F25D40KO',               'F25D90KO',               'M40D90WT',          'F25D90WT',          'M40D90KO',           'F25D90KO',              'M40D230KO',           'M40D230KO']
contrast=['M40D40_KO_vs_WT',   'M40D90_KO_vs_WT',   'F25eRNPD40_KO_vs_WT',    'F25eRNPD90_KO_vs_WT',    'M40WT_D90_vs_D40',  'F25WT_D90_vs_D40',  'M40KO_D90_vs_D40',   'F25eRNP_D90_vs_D40',    'M40KO_D230_vs_D40',   'M40KO_D230_vs_D90']
contra_dic = dict((z[0],list(z[1:])) for z in zip (contrast,sample_1, sample_2, name_1, name_2))

rule all:
    input:
        expand(dir_out+"/pairwise_contr/log/{contrast}.homer.bmk", contrast=contra_dic.keys())

rule contrast:
    input:
        rrbs_1= lambda wildcards: dir_out + "/methylkit/{s1}_pe.bismark.cov".format(s1=contra_dic[wildcards.contrast][0]),
        rrbs_2=lambda wildcards: dir_out + "/methylkit/{s2}_pe.bismark.cov".format(s2=contra_dic[wildcards.contrast][1])
    output:
        dir_out + "/pairwise_contr/{contrast}/methylkit/{contrast}_annotatr_annotation.gene.txt"
    log:
        dir_out+"/pairwise_contr/log/{contrast}.methylkit.log"
    benchmark:
        dir_out+"/pairwise_contr/log/{contrast}.methylkit.bmk"
    params:
        anno=dir_out + "/pairwise_contr/{contrast}/methylkit/{contrast}_annotatr_annotation.txt",
        name1=lambda wildcards: contra_dic[wildcards.contrast][2],
        name2=lambda wildcards: contra_dic[wildcards.contrast][3],
        out=dir_out + "/pairwise_contr/{contrast}/methylkit/",
        prefix="{contrast}",
        hmedip=dir_hmedip + "/macs2/bdgdiff_broad/{contrast}_c3.0_cond2.bed"
    shell:
        '''
        Rscript scripts/methylkit_pairwise_analysis.R {input.rrbs_1} {input.rrbs_2} {params.name1} {params.name2} {params.out} {params.prefix} {params.hmedip} &> {log}
        scripts/SEandGenes.pl {params.anno} > {output}
        '''

rule motif_tile:
    input:
        dir_out + "/pairwise_contr/{contrast}/methylkit/{contrast}_annotatr_annotation.txt"
    output:
        dir_out + "/pairwise_contr/{contrast}/homer/homerResults.html"
    log:
        dir_out + "/pairwise_contr/log/{contrast}.homer.log"
    benchmark:
        dir_out + "/pairwise_contr/log/{contrast}.homer.bmk"
    params:
        motifdir=dir_out + "/pairwise_contr/{contrast}/homer",
        bed=dir_out + "/pairwise_contr/{contrast}/methylkit/{contrast}_hypermethyl_promoter.bed"
    shell:
        '''
        awk '{{if (($10 == "hyper") && ($20 == "hg38_genes_promoters")) {{print $2,$3,$4}} }}' OFS="\t" {input} | uniq -u > {params.bed}
        findMotifsGenome.pl {params.bed} hg38 {params.motifdir} -size given -mask &> {log}
        '''

rule motif_single:
    input:
        dir_out + "/methylkit/{contrast}/{contrast}_Diffmethyl_hyper.bed"
    output:
        dir_out + "/methylkit/{contrast}/homer/homerResults.html"
    log:
        dir_out + "/methylkit/log/{contrast}.homer.log"
    benchmark:
        dir_out + "/methylkit/log/{contrast}.homer.bmk"
    params:
        motifdir=dir_out + "/methylkit/{contrast}/homer",
        anno=dir_out + "/methylkit/{contrast}/{contrast}_hypermethyl_annotation.txt",
        bed=dir_out + "/methylkit/{contrast}/{contrast}_hypermethyl_promoter.bed"
    shell:
        '''
        annotatePeaks.pl {input} hg38 > {params.anno}
        tail -n +2 {params.anno} | awk -v OFS="\t" '{{if ($8=="promoter-TSS") {{print $1,$2,$3,$4,$5}} }}' | sort -k 1,1 -k 2,2n > {params.bed}
        findMotifsGenome.pl {params.bed} hg38 {params.motifdir} -size 200 -mask &> {log}
        '''