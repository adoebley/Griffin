# Download Reference genome
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
# Unzip reference genome
gunzip hg38.fa.gz
# Move the reference genome to Ref directory
mv hg38.fa Ref/

# Download the mappability track
wget https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k100.Umap.MultiTrackMappability.bw
# Move the mappability track to Ref directory
mv k100.Umap.MultiTrackMappability.bw Ref/

# Convert the demo cram file to bam 
samtools view -b -T Ref/hg38.fa -o demo/bam/Healthy_GSM1833219_downsampled.sorted.mini.bam demo/bam/Healthy_GSM1833219_downsampled.sorted.mini.cram
# Index the bam file
samtools index demo/bam/Healthy_GSM1833219_downsampled.sorted.mini.bam

# Navigate to the folder with the griffin_GC_and_mappability_correction snakefile
cd run_demo/griffin_GC_and_mappability_correction/

# Dry run GC and mappability correction
snakemake -s griffin_GC_and_mappability_correction.snakefile --cores 8 -np 

# Navigate to the folder with the griffin_nucleosome_profiling snakefile
cd ../griffin_nucleosome_profiling/
# Dry run Nucleosome profiling
snakemake -s griffin_nucleosome_profiling.snakefile --cores 8 -np