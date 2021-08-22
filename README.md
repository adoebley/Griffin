## Griffin
A flexible framework for nucleosome profiling of cell-free DNA


## Description
To run Griffin, use the snakemakes in the the 'snakemakes' directory

1. griffin_genome_GC_frequncy
    Calculate the frequency of fragments with each GC content across the mappable regions of the reference genome
    For hg38, this step is already complete and results are in Ref/genome_GC_frequency
    Griffin has not been tested on genome builds other than hg38, but this snakemake is provided in case you would like to try a different genome build or different filter for mappable regions
    
2. griffin_GC_correction
    Calculate the GC bias for a given set of bam files
    To run this step:
      create a samples.yaml with your list of bam files and place it in config (see config/example_samples.yaml for format)
      edit config.yaml to provide the path to the reference genome (hg38)
      follow the directions at the top of griffin_GC_correction.snakemake to run the snakemake
      
    Outputs:
      repeat_masker.mapable.k50.Umap.hg38/GC_bias/<sample_name>.GC_bias.txt
        The GC bias of fragments with each length and GC content, this is used for GC correction
      repeat_masker.mapable.k50.Umap.hg38/GC_counts/<sample_name>.GC_counts.txt
        Intermediate file with the number of fragments with each length and GC content
      repeat_masker.mapable.k50.Umap.hg38/GC_plots/
        Assorted plots of the GC bias for each sample
      samples.GC.yaml
        A config file for use in the nucleosome profiling step
      
3. griffin_filter_sites
    If using a new set of sites (not previously filtered) you will need to filter them to remove low mappability sites.
    If you have your own strategy for removing low mappability sites, you can skip this step but will need to add a column with the header 'position' to your sites file for subsequent steps.
    To run this step:
      Create a sites.yaml with paths to your lists of sites and place it in config (see config/example_sites.yaml for format)
      Site lists must be tab separated with a header at the top. At a minimum they must contain columns with the chromosome and position
      Edit config.yaml to specify the location of your mappability track (k50.Umap.MultiTrackMappability.hg38.bw can be downloaded from: https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k50.Umap.MultiTrackMappability.bw)
      Edit config.yaml to specify the name of the column with the chromosome and position or beginning and end of an interval containing the site. 
      If 'position' is a column in your input:
        chrom_column: Chrom
        start_column: position
        end_column: position
      If you only have an interval start and end:
        chrom_column: Chrom
        start_column: Start
        end_column: End   
      Follow the directions at the top of griffin_filter_sites.snakefile to run the snakemake
      
    Outputs:
      sites/<site_list_name>.counts.txt
        Summary of the number of low and high mappability sites
      sites/<site_list_name>.high_mapability.txt
        high mappability sites to be used in subsequent steps
      sites/<site_list_name>.low_mapability.txt
        low mappability sites
      
 4. griffin_nucleosome_profiling
      Run nucleosome profiling for a given set of site lists and a given set of bam files
      To run this step:
        Copy the samples.GC.yaml from the griffin_GC_correction step into the config directory
        Make a sites.yaml containing paths to the high mappability output files from griffin_filter_sites (see config/example_sites.yaml for format)
        Edit config.yaml to provide the path to the reference genome (hg38)
        Edit other config settings as needed
        Follow the directions at the top of griffin_nucleosome_profiling.snakefile to run the snakemake
        
      Outputs:
        results/coverage/all_site/<sample_name>.all_sites.coverage.txt
          nucleosome profiles and metadata for each site list. 
          Both GC corrected and non-GC corrected profiles are in this file and must be separated for downstream analysis (GC_correction column). Coverage profile data is labeled with the start coordinate of the bin. For instance, the column labeled -15 contains the coverage information for -15bp to 0bp relative to the site location.
          results/coverage/<site_name>/<sample_name>.<site_name>.coverage.txt
            These folders contain intermediate files with the coverage profiles for individual site lists. These have been concatenated into results/coverage/all_site/<sample_name>.all_sites.coverage.txt
                  
