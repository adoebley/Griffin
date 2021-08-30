Demo for running the Griffin snakemakes on a small bam file. 

About the demo bam: 
The demo bam is a downsampled version of a healthy donor bam file from GEO. This sample is comprised of a mixture of cfDNA from multiple healthy donors. 

It was originally published in:
Snyder MW, Kircher M, Hill AJ, Daza RM, Shendure J. Cell-free DNA comprises an in vivo nucleosome footprint that informs its tissues-of-origin. Cell. 2016 Jan 14;164(1-2):57-68.

The file was initially downloaded from GEO:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1833219

The file was then converted to bam, realigned to hg38 (http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz), and downsampled to create a small file for this demo using picard downsample with a probability of P=0.004 for a coverage of 0.26x

#downsample took 116 minutes and resulted in a 1.21gb file
time /app/software/Java/1.8.0_181/bin/java -jar $EBROOTPICARD/picard.jar DownsampleSam I=/fh/scratch/delete90/ha_g/realigned_bams/cfDNA_cancer_Snyder_hg38/results/Healthy_GSM1833219/Healthy_GSM1833219_recalibrated.bam O=downsample/Healthy_GSM1833219_downsampled.bam P=0.004
