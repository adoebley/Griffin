#griffin_GC_correction.snakefile
#Anna-Lisa Doebley
#Template made 2021-12-13
#Ha Lab
#Fred Hutchinson Cancer Research Center

"""
#before running snakemake, do in tmux terminal:
ml snakemake/5.19.2-foss-2019b-Python-3.7.4
ml Python/3.7.4-foss-2019b-fh1

#command to run snakemake (remove -np at end when done validating):
snakemake -s griffin_GC_and_mappability_correction.fred_hutch.snakefile --latency-wait 60 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} -J {cluster.JobName}" -j 40 -np

#command to run on restart
#need to add '-q restart-new' according to scicomp. This should be fixed at some point.
snakemake -s griffin_GC_and_mappability_correction.fred_hutch.snakefile --latency-wait 60 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p restart-new -q restart-new --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} -J {cluster.JobName} --requeue" -j 40 -np
"""

configfile: "config/samples.yaml"
configfile: "config/config.fred_hutch.yaml"
configfile: "config/cluster_slurm.yaml"

rule all:
	input: 
		expand("{out_dir}/mappability_bias/{samples}.mappability_bias.txt", samples=config["samples"], out_dir=config['out_dir']),
		expand("{out_dir}/mappability_plots/{samples}.mappability_bias.pdf", samples=config["samples"], out_dir=config['out_dir']),
		expand("{out_dir}/mappability_plots/{samples}.mappability_bias.read_coverage_distribution.pdf", samples=config["samples"], out_dir=config['out_dir']),

		expand("{out_dir}/GC_counts/{samples}.GC_counts.txt", samples=config["samples"], out_dir=config['out_dir']), #GC read counts
		#
		expand("{out_dir}/GC_bias/{samples}.GC_bias.txt", samples=config["samples"], out_dir=config['out_dir']),
		expand("{out_dir}/GC_plots/{samples}.GC_bias.summary.pdf", samples=config["samples"], out_dir=config['out_dir']),
		#
		expand("{out_dir}/samples.GC.yaml", out_dir=config['out_dir'])

rule mappability_bias:
	input:
		bam_file = lambda wildcards: config["samples"][wildcards.samples]
	output:
		mappability_bias = "{out_dir}/mappability_bias/{samples}.mappability_bias.txt",
		mappability_plot = "{out_dir}/mappability_plots/{samples}.mappability_bias.pdf",
		mappability_plot2 = "{out_dir}/mappability_plots/{samples}.mappability_bias.read_coverage_distribution.pdf",
		tmp_dir = temp(directory("{out_dir}/tmp_{samples}/"))
	params:
		sample_name = "{samples}",
		mappability_bw = config['mappability_bw'],
		encode_exclude = config['encode_exclude'],
		centromeres = config['centromeres'],
		gaps = config['gaps'],
		patches = config['patches'],
		alternative_haplotypes = config['alternative_haplotypes'],
		chrom_sizes = config['chrom_sizes'],
		map_q=config['map_quality'],
		CPU = config['mappability_bias']['ncpus'],

		griffin_mappability_script = config['griffin_scripts_dir']+'/griffin_mappability_correction.py'

	shell:
		"time {params.griffin_mappability_script} \
		--bam_file {input.bam_file} \
		--bam_file_name {params.sample_name} \
		--output {output.mappability_bias} \
		--output_plot {output.mappability_plot} \
		--mappability {params.mappability_bw} \
		--exclude_paths {params.encode_exclude} {params.centromeres} {params.gaps} {params.patches} {params.alternative_haplotypes} \
		--chrom_sizes {params.chrom_sizes} \
		--map_quality {params.map_q} \
		--CPU {params.CPU} \
		--tmp_dir {output.tmp_dir}"

rule GC_counts:
	input:
		bam_file = lambda wildcards: config["samples"][wildcards.samples]
	output:
		GC_counts_file = protected("{out_dir}/GC_counts/{samples}.GC_counts.txt"),
	params:
		sample_name = "{samples}",
		mappable_regions_path = config['mappable_regions'],
		ref_seq = config['reference_genome'],
		chrom_sizes = config['chrom_sizes'],
		map_q=config['map_quality'],
		size_range=config['GC_bias_size_range'],
		CPU = config['GC_counts']['ncpus'],

		griffin_GC_counts_script = config['griffin_scripts_dir']+'/griffin_GC_counts.py'

	shell:
		"time {params.griffin_GC_counts_script} \
		--bam_file {input.bam_file} \
		--bam_file_name {params.sample_name} \
		--mappable_regions_path {params.mappable_regions_path} \
		--ref_seq {params.ref_seq} \
		--chrom_sizes {params.chrom_sizes} \
		--out_dir {wildcards.out_dir} \
		--map_q {params.map_q} \
		--size_range {params.size_range} \
		--CPU {params.CPU}"

rule GC_bias:
	input:
		GC_counts_file = "{out_dir}/GC_counts/{samples}.GC_counts.txt"
	output:
		GC_bias_file = "{out_dir}/GC_bias/{samples}.GC_bias.txt",
		GC_plots_file = "{out_dir}/GC_plots/{samples}.GC_bias.summary.pdf"
	params:
		sample_name = "{samples}",
		mappable_name =  config['mappable_regions'].rsplit('/',1)[1].rsplit('.',1)[0],
		genome_GC_frequency = config['genome_GC_frequency'],
		size_range=config['GC_bias_size_range'],
		griffin_GC_bias_script = config['griffin_scripts_dir']+'/griffin_GC_bias.py'

	shell:
		"time {params.griffin_GC_bias_script} \
		--bam_file_name {params.sample_name} \
		--mappable_name {params.mappable_name} \
		--genome_GC_frequency {params.genome_GC_frequency} \
		--out_dir {wildcards.out_dir} \
		--size_range {params.size_range}"

rule make_samples_yaml:
	input:
		mappability_bias = expand("{out_dir}/mappability_bias/{samples}.mappability_bias.txt", samples=config["samples"], out_dir=config['out_dir']),
		GC_bias_file = expand("{out_dir}/GC_bias/{samples}.GC_bias.txt", samples=config["samples"], out_dir=config['out_dir'])#needs to exist for the abspath to work
	output:
		samples_yaml = "{out_dir}/samples.GC.yaml"
	params:
		sample_dict = config['samples'].keys(),
		mappable_name = config['mappable_regions'].rsplit('/',1)[1].rsplit('.',1)[0]
	run:
		import os
		with open(output['samples_yaml'], 'w') as f:
			f.write('samples:\n')
			for key in params['sample_dict']:
				f.write('  '+key+':\n')
				f.write('    bam: '+config['samples'][key]+'\n')
				GC_bias_path = config['out_dir']+'/GC_bias/'+key+'.GC_bias.txt'
				GC_bias_path = os.path.abspath(GC_bias_path)
				f.write('    GC_bias: '+GC_bias_path+'\n')
				mappability_bias_path = config['out_dir']+'/mappability_bias/'+key+'.mappability_bias.txt'
				mappability_bias_path = os.path.abspath(mappability_bias_path)
				f.write('    mappability_bias: '+mappability_bias_path+'\n')

