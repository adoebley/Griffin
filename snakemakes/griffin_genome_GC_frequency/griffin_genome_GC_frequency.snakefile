#griffin_genome_GC_frequency.snakefile
#Anna-Lisa Doebley
#Template made 2021-12-13
#Ha Lab
#Fred Hutchinson Cancer Research Center

"""

#command to run snakemake (remove -np at end when done validating):
snakemake -s griffin_genome_GC_frequency.snakefile -np

#command to run snakemake on a slurm cluster (remove -np at end when done validating):
snakemake -s griffin_genome_GC_frequency.snakefile --latency-wait 60 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} -J {cluster.JobName}" -j 30 -np
"""


configfile: "config/config.yaml"

rule all:
	input:
		expand('{out_dir}/{mappable_name}.{length}bp.GC_frequency.txt', out_dir = config['out_dir'], mappable_name = config['mappable_regions'].rsplit('/',1)[1].rsplit('.',1)[0], length = [m for m in range(config['size_range'][0],config['size_range'][1]+1)])


rule calc_GC_frequency:
	input:
		mappable_path=config['mappable_regions']
	output:
		out_file = '{out_dir}/{mappable_name}.{length}bp.GC_frequency.txt'
	params:
		griffin_GC_frequency_script = config['griffin_scripts_dir']+'/griffin_calc_GC_frequency.py',
		mappable_regions = config['mappable_regions'],
		ref_seq = config['reference_genome'],
		chrom_sizes = config['chrom_sizes'],
		read_length = config['read_length']

	shell:
		'time {params.griffin_GC_frequency_script} \
		--mappable_regions_path {params.mappable_regions} \
		--ref_seq {params.ref_seq} \
		--chrom_sizes {params.chrom_sizes} \
		--out_dir {wildcards.out_dir} \
		--fragment_length {wildcards.length} \
		--read_length {params.read_length}'