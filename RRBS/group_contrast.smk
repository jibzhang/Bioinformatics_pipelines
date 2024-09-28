include: "common_rule.smk"

group1=['s4_TET2KO','s4_DKO','s4_DKO','s8_TET2KO','s8_DKO','s8_DKO']
group2=['s4_WT','s4_WT','s4_TET2KO','s8_WT','s8_WT','s8_TET2KO']
contrast = [f"{a}vs{b}" for a, b in zip(group1, group2)]
contra_dic = dict((z[0],list(z[1:])) for z in zip (contrast,group1,group2))

rule all:
    input:
        expand(dir_out+"/group_contr/log/{contrast}.Homer_motif.bmk", contrast=contra_dic.keys())

rule contrast:
    input:
        sample1 = lambda wildcards: expand(dir_out + '/{s1}/RRBS/bismark_methylextr/{s1}_pe.bismark.cov', s1=groups[contra_dic[wildcards.contrast][0]]),
        sample2 = lambda wildcards: expand(dir_out + '/{s2}/RRBS/bismark_methylextr/{s2}_pe.bismark.cov', s2=groups[contra_dic[wildcards.contrast][1]])
    output:
        hyper = dir_out + "/methylkit/{contrast}/{contrast}_Diffmethyl_25p_hyper.txt",
        hypo=dir_out+ "/methylkit/{contrast}/{contrast}_Diffmethyl_25p_hypo.txt"
    log:
        dir_out+"/group_contr/log/{contrast}.methylkit.log"
    benchmark:
        dir_out+"/group_contr/log/{contrast}.methylkit.bmk"
    params:
        name1=lambda wildcards: contra_dic[wildcards.contrast][0],
        name2=lambda wildcards: contra_dic[wildcards.contrast][1],
    shell:
        '''
        Rscript scripts/methylkit_group_analysis.R {params.name2} {params.name1} &> {log}
        '''

rule Homer:
    input:
        hyper = rules.contrast.output.hyper,
        hypo = rules.contrast.output.hypo
    output:
        hypo_pe = dir_out + "/methylkit/{contrast}/{contrast}_Diffmethyl_25p_hypo_promenhan.bed",
        hyper_pe = dir_out + "/methylkit/{contrast}/{contrast}_Diffmethyl_25p_hyper_promenhan.bed"
    log:
        dir_out + "/group_contr/log/{contrast}.Homer_motif.log"
    benchmark:
        dir_out + "/group_contr/log/{contrast}.Homer_motif.bmk"
    params:
        hyper_bed = dir_out + "/methylkit/{contrast}/{contrast}_Diffmethyl_25p_hyper.bed",
        hypo_bed = dir_out + "/methylkit/{contrast}/{contrast}_Diffmethyl_25p_hypo.bed",
        hyper_dir = dir_out + "/methylkit/{contrast}/Homer/Hyper",
        hypo_dir = dir_out+ "/methylkit/{contrast}/Homer/Hypo",
        promoter_enhancer = "/coh_labs/jochan/CD4_TET2KO_H3K27ac/consensus_peak/merged_enhancers_for_all_plus_promoters.bed"
    shell:
        '''
        cut -f2- {input.hyper} | tail -n +2 > {params.hyper_bed}
        bedtools intersect -a {params.hyper_bed} -b {params.promoter_enhancer} -wa | uniq > {output.hyper_pe}
        cut -f2- {input.hypo} | tail -n +2 > {params.hypo_bed}
        bedtools intersect -a {params.hypo_bed} -b {params.promoter_enhancer} -wa | uniq > {output.hypo_pe}
        if [ $(wc -l < {output.hyper_pe}) -gt 50 ]; then findMotifsGenome.pl {output.hyper_pe} hg38 {params.hyper_dir} -size given &> {log}; fi
        if [ $(wc -l < {output.hypo_pe}) -gt 50 ]; then findMotifsGenome.pl {output.hypo_pe} hg38 {params.hypo_dir} -size given &> {log}; fi
        '''
