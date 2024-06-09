rule mixcr_align:
    input:
        fq1 = dir_out + "/temp/fq_trimed/{sample}_1.fq.gz",
        fq2 = dir_out + "/temp/fq_trimed/{sample}_2.fq.gz"
    output:
        mixcr_align=    dir_out + "/{sample}/" + squence_type + "/mixcr/{sample}.mixcr_align.vdjca",
    log:        dir_out + "/{sample}/" + squence_type + "/log/mixcr/{sample}.mixcr_align.log"
    benchmark:  dir_out + "/{sample}/" + squence_type + "/log/mixcr/{sample}.mixcr_align.bmk"
    threads:
            4
    shell:
        '''
        mixcr align -t {threads} -f -p rna-seq -s hs -OallowPartialAlignments=true {input.fq1} {input.fq2} {output.mixcr_align} > {log}
        '''

rule mixcr_assemblePartial:
    input:
        mixcr_align =  rules.mixcr_align.output.mixcr_align
    output:
        mixcr_align_rescue= dir_out + "/{sample}/" + squence_type + "/mixcr/{sample}.mixcr_assemblePartial.vdjca",
    log:                    dir_out + "/{sample}/" + squence_type + "/log/mixcr/{sample}.mixcr_assemblePartial.log"
    benchmark:              dir_out + "/{sample}/" + squence_type + "/log/mixcr/{sample}.mixcr_assemblePartial.bmk"
    params:
        rescue1=dir_out + "/{sample}/" + squence_type + "/mixcr/{sample}.alignRescue1.vdjca",
        rescue2=dir_out + "/{sample}/" + squence_type + "/mixcr/{sample}.alignRescue2.vdjca",
    threads:
            1
    shell:
        '''
        mixcr assemblePartial -f {input.mixcr_align} {params.rescue1} > {log}
        mixcr assemblePartial -f {params.rescue1} {params.rescue2} >> {log}
        mixcr assemblePartial -f {params.rescue2} {output.mixcr_align_rescue} >> {log}
        '''

rule mixcr_extend:
    input:
        mixcr_align_rescue = rules.mixcr_assemblePartial.output.mixcr_align_rescue
    output:
        mixcr_extend=       dir_out + "/{sample}/" + squence_type + "/mixcr/{sample}.mixcr_extend.vdjca",
    log:                    dir_out + "/{sample}/" + squence_type + "/log/mixcr/{sample}.mixcr_extend.log"
    benchmark:              dir_out + "/{sample}/" + squence_type + "/log/mixcr/{sample}.mixcr_extend.bmk"
    threads:
            4
    shell:
        '''
        mixcr extend -f -t {threads} {input.mixcr_align_rescue} {output.mixcr_extend} > {log}
        '''

rule mixcr_assemble:
    input:
        mixcr_extend =  rules.mixcr_extend.output.mixcr_extend
    output:
        assemble=       dir_out + "/{sample}/" + squence_type + "/mixcr/{sample}.mixcr_assemble.clones.clns",
    log:                dir_out + "/{sample}/" + squence_type + "/log/mixcr/{sample}.mixcr_assemble.log"
    benchmark:          dir_out + "/{sample}/" + squence_type + "/log/mixcr/{sample}.mixcr_assemble.bmk"
    threads:
            4
    shell:
        '''
        mixcr assemble -t {threads} {input.mixcr_extend} {output.assemble} > {log}
        '''

rule mixcr_exportClones:
    input:
        assemble =  rules.mixcr_assemble.output.assemble
    output:
        exportTRAClones =  dir_out + "/{sample}/" + squence_type + "/mixcr/{sample}.mixcr_TRA_exportClones.tsv",
        exportTRBClones =  dir_out + "/{sample}/" + squence_type + "/mixcr/{sample}.mixcr_TRB_exportClones.tsv"
    log:                dir_out + "/{sample}/" + squence_type + "/log/mixcr/{sample}.mixcr_exportClones.log"
    benchmark:          dir_out + "/{sample}/" + squence_type + "/log/mixcr/{sample}.mixcr_exportClones.bmk"
    threads:
            1
    shell:
        '''
        mixcr exportClones -c TRA -f {input.assemble} {output.exportTRAClones} > {log}
        mixcr exportClones -c TRB -f {input.assemble} {output.exportTRBClones} 2> {log}
        '''


