#griffin_genome_GC_frequency.snakefile
#Anna-Lisa Doebley
#Template made 2021-04-06
#Ha Lab
#Fred Hutchinson Cancer Research Center

"""
#before running snakemake, do in tmux terminal:
ml snakemake/5.19.2-foss-2019b-Python-3.7.4
ml Python/3.7.4-foss-2019b-fh1
PATH="$PATH:/fh/fast/ha_g/user/adoebley/projects/griffin_paper/scripts/"


#command to run snakemake (remove -np at end when done validating):
snakemake -s griffin_genome_GC_frequency.snakefile --latency-wait 60 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} -J {cluster.JobName}" -j 30 -np
"""


configfile: "config/config.yaml"

rule all:
	input:
		expand('{out_dir}/{mapable_name}.{length}bp.GC_frequency.txt', out_dir = config['out_dir'], mapable_name = config['mapable_regions'].rsplit('/',1)[1].rsplit('.',1)[0], length = [m for m in range(config['size_range'][0],config['size_range'][1]+1)])


rule calc_GC_frequency:
	input:
		mapable_path=config['mapable_regions']
	output:
		out_file = '{out_dir}/{mapable_name}.{length}bp.GC_frequency.txt'
	params:
		mapable_file = config['mapable_regions'],
		ref_seq = config['reference_genome'],
		chrom_sizes = config['chrom_sizes']

	shell:
		'griffin_calc_GC_frequency.py --mapable_regions {params.mapable_file} --ref_seq {params.ref_seq} --chrom_sizes {params.chrom_sizes} \
		--out_dir {wildcards.out_dir} --fragment_length {wildcards.length}'