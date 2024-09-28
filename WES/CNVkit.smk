rule CNVkit_coverage:
	input:
		bam=dir_out + "/{sample}/" + squence_type + "/temp/samtools_sort/{sample}.bam",
		bai=dir_out + "/{sample}/" + squence_type + "/temp/samtools_sort/{sample}.bai"
	output:
		target_coverage=dir_out + "/{sample}/" + squence_type + "/CNVkit/{sample}.targetcoverage.cnn",
		antitarget_coverage=dir_out + "/{sample}/" + squence_type + "/CNVkit/{sample}.antitargetcoverage.cnn"
	benchmark:
		dir_out+"/{sample}/"+squence_type+"/log/{sample}_CNVkit.bmk"
	params:
		target= config["cnvkit_target"],
		antitarget= config["cnvkit_antitarget"]
	threads: 24
	shell:
		'''
		cnvkit.py coverage {input.bam} {params.target} -p {threads} -o {output.target_coverage}
		cnvkit.py coverage {input.bam} {params.antitarget} -p {threads} -o {output.antitarget_coverage}
		'''

rule CNVkit_fix:
	input:
		target = rules.CNVkit_coverage.output.target_coverage,
		antitarget = rules.CNVkit_coverage.output.antitarget_coverage
	output:
		dir_out + "/{sample}/" + squence_type + "/CNVkit/{sample}.cnr"
	benchmark:
		dir_out+"/{sample}/"+squence_type+"/log/{sample}_cnvkitfix.bmk"
	params:
		cnvkit_ref=config["cnvkit_ref"]
	shell:
		'''
		cnvkit.py fix {input.target} {input.antitarget} {params.cnvkit_ref} -o {output}
		'''

rule CNVkit_fix_noPON:
	input:
		target = rules.CNVkit_coverage.output.target_coverage,
		antitarget = rules.CNVkit_coverage.output.antitarget_coverage
	output:
		dir_out + "/{sample}/" + squence_type + "/CNVkit_noPON/{sample}.cnr"
	benchmark:
		dir_out+"/{sample}/"+squence_type+"/log/{sample}_cnvkitfix_noPON.bmk"
	params:
		cnvkit_ref=config["cnvkit_ref_flat"]
	shell:
		'''
		cnvkit.py fix {input.target} {input.antitarget} {params.cnvkit_ref} -o {output}
		'''

rule CNVkit_get_cns:
	input:
		fixed = rules.CNVkit_fix.output,
		vcf = dir_out + "/{sample}/" + squence_type + "/variant_filter/{sample}_filtered.vcf"
	output:
		segment = dir_out + "/{sample}/" + squence_type + "/CNVkit/{sample}.cns",
		Int_CNV = dir_out + "/{sample}/" + squence_type + "/CNVkit/{sample}.call.cns"
	shell:
		'''
 		cnvkit.py segment {input.fixed} --drop-low-coverage --vcf {input.vcf} --min-variant-depth 10 -o {output.segment} 
 		cnvkit.py call {output.segment} -y -v {input.vcf} -m threshold -t=-1.1,-0.4,0.3,0.7 -o {output.Int_CNV} 
 		'''

rule CNVkit_get_cns_noPON:
	input:
		fixed = rules.CNVkit_fix_noPON.output,
		vcf = dir_out + "/{sample}/" + squence_type + "/variant_filter/{sample}_filtered.vcf"
	output:
		segment = dir_out + "/{sample}/" + squence_type + "/CNVkit_noPON/{sample}.cns"
	shell:
		'''
 		cnvkit.py segment {input.fixed} --drop-low-coverage --vcf {input.vcf} --min-variant-depth 10 -o {output.segment} 
 		'''

rule CNVkit_get_genemetrics:
	input:
		cnr = rules.CNVkit_fix.output,
		cns = rules.CNVkit_get_cns.output.segment
	output: dir_out + "/{sample}/" + squence_type + "/CNVkit/{sample}.genemetrics.txt"
	params: gender = lambda wildcards: sex[wildcards.sample]
	shell:
		'''
		cnvkit.py genemetrics {input.cnr} -s {input.cns} -m 20 -y -x {params.gender} -o {output}
		'''

rule CNVkit_noPON_get_genemetrics:
	input:
		cnr = rules.CNVkit_fix_noPON.output,
		cns = rules.CNVkit_get_cns_noPON.output.segment
	output: dir_out + "/{sample}/" + squence_type + "/CNVkit_noPON/{sample}.genemetrics.txt"
	params: gender = lambda wildcards: sex[wildcards.sample]
	shell:
		'''
		cnvkit.py genemetrics {input.cnr} -s {input.cns} -m 3 -x {params.gender} -o {output}
		'''

rule CNVkit_scatter:
	input:
		cnr=rules.CNVkit_fix.output,
		cns=rules.CNVkit_get_cns.output.segment
	output: dir_out + "/{sample}/" + squence_type + "/CNVkit/{sample}.scatter.pdf"
 	params:
 		title = "CNVkit_scatter_plot_for_{sample}"
	shell:
 		'''
 		cnvkit.py scatter {input.cnr} -s {input.cns} -o {output} --title {params.title}
		'''

rule CNVkit_scatter_noPON:
	input:
		cnr=rules.CNVkit_fix_noPON.output,
		cns=rules.CNVkit_get_cns_noPON.output.segment
	output: dir_out + "/{sample}/" + squence_type + "/CNVkit_noPON/{sample}.scatter.pdf"
 	params:
 		title = "CNVkit_scatter_plot_for_{sample}"
	shell:
 		'''
 		cnvkit.py scatter {input.cnr} -s {input.cns} -o {output} --title {params.title}
		'''

rule CNVkit_diagram:
 	input:
 		cnr=rules.CNVkit_fix.output,
 		cns=rules.CNVkit_get_cns.output.segment
 	output: dir_out + "/{sample}/" + squence_type + "/CNVkit/{sample}.diagram.pdf"
	params:
 		title = "CNVkit_diagram_plot_for_{sample}"
	shell:
		'''
		cnvkit.py diagram  --segment {input.cns} {input.cnr} -o {output} --title {params.title}
		'''

rule CNVkit_diagram_noPON:
 	input:
 		cnr=rules.CNVkit_fix_noPON.output,
 		cns=rules.CNVkit_get_cns_noPON.output.segment
 	output: dir_out + "/{sample}/" + squence_type + "/CNVkit_noPON/{sample}.diagram.pdf"
	params:
 		title = "CNVkit_diagram_plot_for_{sample}"
	shell:
		'''
		cnvkit.py diagram  --segment {input.cns} {input.cnr} -o {output} --title {params.title}
		'''

rule CNVkit_heatmap:
	input:
		lambda wildcards: [dir_out + "/"+ ind + "/" + squence_type + "/CNVkit/" + ind + ".call.cns" for ind in EBV_dict[wildcards.EBV]]
	output: dir_out + "/CNVkit/EBV_{EBV}/CNVkit_Heatmap_{EBV}.pdf"
 	params: title = "CNVkit_Heatmap_plot_{EBV}"
	shell:
		'''
		cnvkit.py heatmap {input} -d -v -o {output} --title {params.title}
		'''

rule CNVkit_metrics:
	input:
		cnr = lambda wildcards: [dir_out + "/" + ind + "/" + squence_type + "/CNVkit/" + ind + ".cnr" for ind in EBV_dict[wildcards.EBV]],
		cns = lambda wildcards: [dir_out + "/" + ind + "/" + squence_type + "/CNVkit/" + ind + ".cns" for ind in EBV_dict[wildcards.EBV]],
		genemetrics = lambda wildcards: [dir_out + "/" + ind + "/" + squence_type + "/CNVkit/" + ind + ".genemetrics.txt" for ind in EBV_dict[wildcards.EBV]]
	output: dir_out + "/CNVkit/EBV_{EBV}/CNVkit_metrics_{EBV}.txt"
	benchmark:
		dir_out + "/log/CNVkit_metrics_{EBV}.bmk"
	params:
		outdir = dir_out + "/CNVkit/EBV_{EBV}"
	shell:
		'''
		cnvkit.py metrics {input.cnr} -s {input.cns} --drop-low-coverage -o {output}
		ln -s {input.genemetrics} {params.outdir}
		'''

rule CNVkit_metrics_noPON:
	input:
		cnr = lambda wildcards: [dir_out + "/" + ind + "/" + squence_type + "/CNVkit_noPON/" + ind + ".cnr" for ind in EBV_dict[wildcards.EBV]],
		cns = lambda wildcards: [dir_out + "/" + ind + "/" + squence_type + "/CNVkit_noPON/" + ind + ".cns" for ind in EBV_dict[wildcards.EBV]],
		genemetrics = lambda wildcards: [dir_out + "/" + ind + "/" + squence_type + "/CNVkit_noPON/" + ind + ".genemetrics.txt" for ind in EBV_dict[wildcards.EBV]]
	output: dir_out + "/CNVkit_noPON/EBV_{EBV}/CNVkit_metrics_{EBV}.txt"
	benchmark:
		dir_out + "/log/CNVkit_noPON_metrics_{EBV}.bmk"
	params:
		outdir = dir_out + "/CNVkit_noPON/EBV_{EBV}"
	shell:
		'''
		cnvkit.py metrics {input.cnr} -s {input.cns} --drop-low-coverage -o {output}
		ln -s {input.genemetrics} {params.outdir}
		'''

rule CNVkit_gistic:
	input:
		cnr=lambda wildcards: [dir_out + "/" + ind + "/" + squence_type + "/CNVkit/" + ind + ".cnr" for ind in EBV_dict[wildcards.EBV]],
		cns=lambda wildcards: [dir_out + "/" + ind + "/" + squence_type + "/CNVkit/" + ind + ".cns" for ind in EBV_dict[wildcards.EBV]]
	benchmark:
		dir_out + "/log/CNVkit_gistic_{EBV}.bmk"
	output:
		seg = dir_out + "/CNVkit/EBV_{EBV}/CNVkit_gistic_{EBV}.seg",
		marker = dir_out + "/CNVkit/EBV_{EBV}/CNVkit_gistic_marker_{EBV}.txt"
	shell:
		'''
        cnvkit.py export seg {input.cns} -o {output.seg}
        cnvkit.py export gistic {input.cnr} -o {output.marker}
        '''

rule CNVkit_gistic_noPON:
	input:
		cnr=lambda wildcards: [dir_out + "/" + ind + "/" + squence_type + "/CNVkit_noPON/" + ind + ".cnr" for ind in EBV_dict[wildcards.EBV]],
		cns=lambda wildcards: [dir_out + "/" + ind + "/" + squence_type + "/CNVkit_noPON/" + ind + ".cns" for ind in EBV_dict[wildcards.EBV]]
	benchmark:
		dir_out + "/log/CNVkit_noPON_gistic_{EBV}.bmk"
	output:
		seg = dir_out + "/CNVkit_noPON/EBV_{EBV}/CNVkit_gistic_{EBV}.seg",
		marker = dir_out + "/CNVkit_noPON/EBV_{EBV}/CNVkit_gistic_marker_{EBV}.txt"
	shell:
		'''
        cnvkit.py export seg {input.cns} -o {output.seg}
        cnvkit.py export gistic {input.cnr} -o {output.marker}
        '''


rule done_CNVkit:
    input:
        diagram = rules.CNVkit_diagram.output,
        scatter = rules.CNVkit_scatter.output,
	genes = rules.CNVkit_get_genemetrics.output
    output:
        done = dir_out + "/{sample}/" + squence_type + "/log/doneCNVkit.txt"
    shell:
        '''
        touch {output.done}
        '''