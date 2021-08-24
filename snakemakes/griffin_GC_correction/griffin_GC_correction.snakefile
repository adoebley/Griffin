#griffin_GC_correction.snakefile
#Anna-Lisa Doebley
#Template made 2021-04-06
#Ha Lab
#Fred Hutchinson Cancer Research Center

"""
#before running snakemake, do in tmux terminal:
ml snakemake/5.19.2-foss-2019b-Python-3.7.4
ml Python/3.7.4-foss-2019b-fh1 #this loads deeptools
PATH="$PATH:/fh/fast/ha_g/user/adoebley/projects/griffin_paper/Griffin/scripts/"

#command to run snakemake (remove -np at end when done validating):
snakemake -s griffin_GC_correction.snakefile --latency-wait 60 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} -J {cluster.JobName}" -j 40 -np

#command to run on restart
#need to add '-q restart-new' according to scicomp. This should be fixed at some point.
snakemake -s griffin_GC_correction.snakefile --latency-wait 60 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p restart-new -q restart-new --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} -J {cluster.JobName} --requeue" -j 40 -np
"""

configfile: "config/samples.yaml"
configfile: "config/config.yaml"
configfile: "config/cluster_slurm.yaml"

rule all:
	input: #commented out files are produced at the same time as other files
		expand("{out_dir}/{mapable_name}/GC_counts/{samples}.GC_counts.txt", samples=config["samples"], out_dir=config['out_dir'], mapable_name = config['mapable_regions'].rsplit('/',1)[1].rsplit('.',1)[0]), #GC read counts
		#
		expand("{out_dir}/{mapable_name}/GC_bias/{samples}.GC_bias.txt", samples=config["samples"], out_dir=config['out_dir'], mapable_name = config['mapable_regions'].rsplit('/',1)[1].rsplit('.',1)[0]),
		expand("{out_dir}/{mapable_name}/GC_plots/{samples}.GC_bias.summary.pdf", samples=config["samples"], out_dir=config['out_dir'], mapable_name = config['mapable_regions'].rsplit('/',1)[1].rsplit('.',1)[0]),
		#
		expand("{out_dir}/samples.GC.yaml", out_dir=config['out_dir'])


rule GC_counts:
	input:
		bam_file = lambda wildcards: config["samples"][wildcards.samples]
	output:
		GC_counts_file = protected("{out_dir}/{mapable_name}/GC_counts/{samples}.GC_counts.txt"),
	params:
		sample_name = "{samples}",
		mapable_regions = config['mapable_regions'],
		ref_seq = config['reference_genome'],
		chrom_sizes = config['chrom_sizes'],
		map_q=config['map_quality'],
		size_range=config['GC_bias_size_range'],
		CPU = config['GC_counts']['ncpus']
	shell:
		"time griffin_GC_counts.py --bam_file {input.bam_file} --bam_file_name {params.sample_name} \
		--mapable_regions {params.mapable_regions} --ref_seq {params.ref_seq} --chrom_sizes {params.chrom_sizes} \
		--out_dir {wildcards.out_dir} --map_q {params.map_q} --size_range {params.size_range} --CPU {params.CPU}"

rule GC_bias:
	input:
		GC_counts_file = "{out_dir}/{mapable_name}/GC_counts/{samples}.GC_counts.txt"
	output:
		GC_bias_file = "{out_dir}/{mapable_name}/GC_bias/{samples}.GC_bias.txt",
		GC_plots_file = "{out_dir}/{mapable_name}/GC_plots/{samples}.GC_bias.summary.pdf"
	params:
		sample_name = "{samples}",
		mapable_name =  config['mapable_regions'].rsplit('/',1)[1].rsplit('.',1)[0],
		genome_GC_frequency = config['genome_GC_frequency'],
		size_range=config['GC_bias_size_range']
	shell:
		"time griffin_GC_bias.py --bam_file_name {params.sample_name} --mapable_name {params.mapable_name} \
		--genome_GC_frequency {params.genome_GC_frequency} --out_dir {wildcards.out_dir} --size_range {params.size_range}"

rule make_samples_yaml:
	input:
		GC_bias_file = expand("{out_dir}/{mapable_name}/GC_bias/{samples}.GC_bias.txt", samples=config["samples"], out_dir=config['out_dir'], mapable_name = config['mapable_regions'].rsplit('/',1)[1].rsplit('.',1)[0])#needs to exist for the abspath to work
	output:
		samples_yaml = "{out_dir}/samples.GC.yaml"
	params:
		sample_dict = config['samples'].keys(),
		mapable_name = config['mapable_regions'].rsplit('/',1)[1].rsplit('.',1)[0]
	run:
		import os
		with open(output['samples_yaml'], 'w') as f:
			f.write('samples:\n')
			for key in params['sample_dict']:
				f.write('  '+key+':\n')
				f.write('    bam: '+config['samples'][key]+'\n')
				GC_bias_path = config['out_dir']+'/'+params['mapable_name']+'/GC_bias/'+key+'.GC_bias.txt'
				GC_bias_path = os.path.abspath(GC_bias_path)
				f.write('    GC_bias: '+GC_bias_path+'\n')
