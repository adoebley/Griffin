Demo for running the Griffin snakemakes on a small bam file. 


This demo uses a conda environment. You will need to install conda:
	conda 4.10.3 (installed with miniconda https://docs.conda.io/en/latest/miniconda.html)
	(other versions may work, this is just the one used for testing)


For all pipelines start by initializing the conda environment and installing the necessary packages and software within the environment:
	conda create --name griffin_demo python=3.7.4
	conda activate griffin_demo
	pip install snakemake==5.5.4
	pip install pandas==1.2.4
	pip install scipy==1.7.1
	pip install pysam==0.15.4
	pip install pyBigWig==0.3.17
	pip install matplotlib==3.4.1
	conda install samtools=1.13


#####################
griffin_filter_sites (total time to run this demo step 10-15 minutes)

1. If you haven't already activated the conda environment created above, activate it:
	conda activate griffin_demo

2. Copy snakemakes/griffin_filter_sites/ to a location where you would like to do the analysis (ex. a directory called run_demo)
	mkdir run_demo
	cp -r snakemakes/griffin_filter_sites run_demo

3. Copy the sites.yaml from demo/griffin_filter_sites_demo_files/config into the snakemake config folder
	cp demo/griffin_filter_sites_demo_files/config/sites.yaml run_demo/griffin_filter_sites/config/

4. Download the mapability track (1.2gb) from the link below (if you don't have wget, you can open the link in a browser) and put it in a location of your choice (ex. Ref)
	wget https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k50.Umap.MultiTrackMappability.bw
	mv k50.Umap.MultiTrackMappability.bw Ref/

5. Navigate to the scripts folder and add it to your path:
	cd scripts/
	export PATH=$PATH:$PWD

6. Navigate to the folder with the snakefile
	cd ../run_demo/griffin_filter_sites/

7. Open the config.yaml file (run_demo/griffin_filter_sites/config/config.yaml) and update the paths to the reference files to point to the location of these files on your computer (either relative or absolute path):
	mapability_track: ../../Ref/k50.Umap.MultiTrackMappability.bw

8. Open the sites.yaml (run_demo/griffin_filter_sites/config/sites.yaml) and update the path to the demo sites file:
	CTCF_demo: ../../demo/griffin_filter_sites_demo_files/input/CTCF.hg38.2000demo.txt

9. Run the snakemake (expected runtime 11 seconds):
	snakemake -s griffin_filter_sites.snakefile -np #dry run to print a list of jobs
	snakemake -s griffin_filter_sites.snakefile

10. The expected outputs are 1000 high mapability sites and 1000 low mapability sites.


#####################
griffin_GC_correction (total time to run the griffin_GC_correction demo ~45 minutes)

1. If you haven't already activated the conda environment created above, activate it:
	conda activate griffin_demo

2. Copy snakemakes/griffin_GC_correction/ to a location where you would like to do the analysis (ex. a directory called run_demo)
	mkdir run_demo
	cp -r snakemakes/griffin_GC_correction run_demo

3. Copy the samples.yaml from demo/griffin_GC_correction_demo_files/config into the snakemake config folder
	cp demo/griffin_GC_correction_demo_files/config/samples.yaml run_demo/griffin_GC_correction/config/

4. Download the reference genome from the link below (if you don't have wget, you can open the link in a browser), unzip it, and put it in a location of your choice (ex. Ref)
	wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
	gunzip hg38.fa.gz
	mv hg38.fa Ref/

5. Convert the demo cram file to bam and create an index (takes ~1 minute):
	samtools view -b -T Ref/hg38.fa -o demo/bam/Healthy_GSM1833219_downsampled.sorted.mini.bam demo/bam/Healthy_GSM1833219_downsampled.sorted.mini.cram
	samtools index demo/bam/Healthy_GSM1833219_downsampled.sorted.mini.bam

6. Navigate to the scripts folder and add this folder to your path:
	cd scripts/
	export PATH=$PATH:$PWD

7. Navigate to the folder with the snakefile
	cd ../run_demo/griffin_GC_correction/

8. Open the config.yaml file (run_demo/griffin_GC_correction/config/config.yaml) and update the paths to the reference files to point to the location of these files on your computer (either relative or absolute path):
	reference_genome: ../../Ref/hg38.fa
	chrom_sizes: ../../Ref/hg38.standard.chrom.sizes
	mapable_regions: ../../Ref/repeat_masker.mapable.k50.Umap.hg38.bedGraph
	genome_GC_frequency: ../../Ref/genome_GC_frequency

9. Open the samples.yaml (run_demo/griffin_GC_correction/config/samples.yaml) and update the path to the demo bam file:
	Healthy_demo: ../../demo/bam/Healthy_GSM1833219_downsampled.sorted.mini.bam

10. If you do not want the snakemake to use 8 cores, open the cluster_slurm.yaml (run_demo/griffin_GC_correction/config/cluster_slurm.yaml) and edit the ncpus parameter (line 18) for the GC_counts step (other parameters in this file are not used unless launching to a slurm cluster). Increasing the CPUs will make the snakemake run faster, as long as your computer has the CPUs available. 

11. Run the snakemake (expected runtime 30 minutes - 8 cores):
	snakemake -s griffin_GC_correction.snakefile -np #dry run to print a list of jobs
	snakemake -s griffin_filter_sites.snakefile

12. The outputs should be very similar to the expected outputs (plots should appear identical) but there might be small differences due to the way the GC correction handles unknown bases 

#####################
griffin_nucleosome_profiling (total time to run the griffin_nucleosome_profiling demo ~15 minutes)

1. If you haven't already activated the conda environment created above, activate it:
	conda activate griffin_demo

2. Copy snakemakes/griffin_nucleosome_profiling/ to a location where you would like to do the analysis (ex. a directory called run_demo)
	mkdir run_demo
	cp -r snakemakes/griffin_nucleosome_profiling run_demo

3. Copy the samples.GC.yaml and sites.yaml from demo/griffin_nucleosome_profiling_demo_files/config into the snakemake config folder
	cp demo/griffin_nucleosome_profiling_demo_files/config/* run_demo/griffin_nucleosome_profiling/config/

4. If you haven't already downloaded the reference genome, download it from the link below (if you don't have wget, you can open the link in a browser), unzip it, and put it in a location of your choice (ex. Ref)
	wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
	gunzip hg38.fa.gz
	mv hg38.fa Ref/

5. Convert the demo cram file to bam and create an index (takes ~1 minute):
	samtools view -b -T Ref/hg38.fa -o demo/bam/Healthy_GSM1833219_downsampled.sorted.mini.bam demo/bam/Healthy_GSM1833219_downsampled.sorted.mini.cram
	samtools index demo/bam/Healthy_GSM1833219_downsampled.sorted.mini.bam

6. Navigate to the scripts folder and add this folder to your path:
	cd scripts/
	export PATH=$PATH:$PWD

7. Navigate to the folder with the snakefile
	cd ../run_demo/griffin_nucleosome_profiling/

8. Open the config.yaml file (run_demo/griffin_nucleosome_profiling/config/config.yaml) and update the path to the reference genome to point to the location of this file on your computer (either relative or absolute path):
	reference_genome: ../../Ref/hg38.fa

9. Open the sites.yaml (run_demo/griffin_nucleosome_profiling/config/sites.yaml) and update the path to the demo sites file (if you have run the filter sites demo, you can use the path to your results instead):
  CTCF_demo: ../../demo/griffin_filter_sites_demo_files/expected_results/CTCF_demo.high_mapability.txt

10. Open the samples.GC.yaml (run_demo/griffin_nucleosome_profiling/config/samples.yaml) and update the path to the demo bam file and GC correction file (if you have run the GC correction demo, you can use the path to your results instead):
	samples:
	  Healthy_demo:
	    bam: ../../demo/bam/Healthy_GSM1833219_downsampled.sorted.mini.bam
	    GC_bias: ../../demo/griffin_GC_correction_demo_files/expected_results/Healthy_demo.GC_bias.txt

11. Run the snakemake (runtime ~1minute):
	snakemake -s griffin_nucleosome_profiling.snakefile -np #dry run to print a list of jobs
	snakemake -s griffin_nucleosome_profiling.snakefile

12. The results should be identical to the expected results in demo/griffin_nucleosome_profiling_demo_files/expected_results (if you used your own GC_correction results, they may be slightly different)



#####################
About the demo bam: 
The demo bam is a downsampled version of a healthy donor bam file from GEO. This sample is comprised of a mixture of cfDNA from multiple healthy donors. 

It was originally published in:
Snyder MW, Kircher M, Hill AJ, Daza RM, Shendure J. Cell-free DNA comprises an in vivo nucleosome footprint that informs its tissues-of-origin. Cell. 2016 Jan 14;164(1-2):57-68.

The file was initially downloaded from GEO where it is freely available:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1833219

The file was then converted to bam, realigned to hg38 (http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz), and downsampled to create a small file for this demo using picard downsample with a probability of P=0.004 for a coverage of 0.26x, to further shrink the file, it was converted to cram with cramtools 


