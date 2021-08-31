#griffin_filter_sites.snakefile
#Anna-Lisa Doebley
#Template made 2021-05-16
#Ha Lab
#Fred Hutchinson Cancer Research Center

"""
#before running snakemake, do in terminal:
PATH="$PATH:/path/to/scripts/"


#command to run snakemake (remove -np at end when done validating):
snakemake -s griffin_filter_sites.snakefile -np

#or

#command for running on a slurm cluster
snakemake -s griffin_filter_sites.snakefile --latency-wait 60 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} -J {cluster.JobName}" -j 40 -np
"""

configfile: "config/sites.yaml"
configfile: "config/config.yaml"

rule all:
	input:
		expand("{map_dir}/{site_files}.high_mapability.txt", site_files=config["site_files"], map_dir=config['split_by_mapability_dir']),
		expand("{map_dir}/{site_files}.low_mapability.txt", site_files=config["site_files"], map_dir=config['split_by_mapability_dir'])


rule split_by_mapability:
	input:
		site_files = lambda wildcards: config["site_files"][wildcards.site_files]
	output:
		high_mapability = "{map_dir}/{site_files}.high_mapability.txt",
		low_mapability = "{map_dir}/{site_files}.low_mapability.txt"
	params:
		file_name = '{site_files}',
		map_dir=config['split_by_mapability_dir'],
		mapability=config['mapability_track'],
		chrom_col=config['chrom_column'],
		start_col=config['start_column'],
		end_col=config['end_column'],
		strand_column=config['strand_column'],
		chroms = config['chroms'],
		norm_window_values=config['norm_window_values'],
        targeted_panel=config['targeted_panel'],
        targeted_window_columns=config['targeted_window_columns'],
		threshold=config['threshold']
	shell:
		"griffin_filter_sites.py -i {input.site_files} --name {params.file_name} -o {params.map_dir} -m {params.mapability} \
		-c {params.chrom_col} -s {params.start_col} -e {params.end_col} --strand_column {params.strand_column} --chroms {params.chroms} \
		--window_values {params.norm_window_values} --targeted_panel {params.targeted_panel} \
		--targeted_window_columns {params.targeted_window_columns} \
        --threshold {params.threshold}"

