Demo for running the Griffin snakemakes on a small bam file. 

About the demo bam: 
The demo bam is a downsampled version of a healthy donor bam file from GEO. This sample is comprised of a mixture of cfDNA from multiple healthy donors. 

It was originally published in:
Snyder MW, Kircher M, Hill AJ, Daza RM, Shendure J. Cell-free DNA comprises an in vivo nucleosome footprint that informs its tissues-of-origin. Cell. 2016 Jan 14;164(1-2):57-68.

The file was initially downloaded from GEO:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1833219

The file was then converted to bam, realigned to hg38 (http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz), and downsampled to create a small file for this demo using picard downsample with a probability of P=0.004 for a coverage of 0.26x

Prior to running these demos you will need to install:

conda 4.10.3 (installed with miniconda https://docs.conda.io/en/latest/miniconda.html)
(other versions may work, this is just the one used for testing)

For all snakemakes you can start by initializing the conda environment and installing the necessary packages within the environment:
	conda create --name griffin_demo python=3.7.4
	conda activate griffin_demo
	pip install snakemake==5.5.4
	pip install pandas==1.2.4
	pip install scipy==1.7.1
	pip install pysam==0.15.4
	pip install pyBigWig==0.3.17


griffin_filter_sites

1. If you haven't already activated the conda environment created above, activate it:
	conda activate griffin_demo

1. Copy snakemakes/griffin_filter_sites/ to a location where you would like to do the analysis (ex. a directory called run_demo)
	mkdir run_demo
	cp -r snakemakes/griffin_filter_sites run_demo

2. Copy the config files from demo/griffin_filter_sites_demo_files/config into the snakemake config folder
	cp demo/griffin_filter_sites_demo_files/config/* run_demo/griffin_filter_sites/config/

3. Download the mapability track (1.2gb) from the link below and put it in a location of your choice (ex Ref)
	wget https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k50.Umap.MultiTrackMappability.bw
	mv k50.Umap.MultiTrackMappability.bw Ref/
	(you can also open the link in a browser and download the file that way)

4. Navigate to the scripts folder and add it to your path:
	cd scripts/
	export PATH=$PATH:$PWD

4. Navigate to the folder with the snakefile
	cd ../run_demo/griffin_filter_sites/

5. Open the config.yaml file (run_demo/griffin_filter_sites/config/config.yaml) and update the path to the mapability track:
	mapability_track: ../../Ref/k50.Umap.MultiTrackMappability.bw

6. Open the sites.yaml (run_demo/griffin_filter_sites/config/sites.yaml) and update the path to the demo sites file:
	CTCF_demo: ../../demo/griffin_filter_sites_demo_files/input/CTCF.hg38.2000demo.txt

7. Run the snakemake:
	snakemake -s griffin_filter_sites.snakefile -np #dry run to print a list of jobs
	snakemake -s griffin_filter_sites.snakefile

8. Check that the outputs (in the 'sites' directory) match the expected outputs (1000 high mapability sites and 1000 low mapability sites)







