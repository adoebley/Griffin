{\rtf1\ansi\ansicpg1252\cocoartf2639
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fmodern\fcharset0 Courier;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\deftab720
\pard\pardeftab720\partightenfactor0

\f0\fs26 \cf0 \expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 `encode_unified_GRCh38_exclusion_list.bed`  \
This is the encode unified GRCh38 exclusion list (https://www.encodeproject.org/files/ENCFF356LFX/)\
\
`hg38.standard.chrom.sizes`  \
Chromosome sizes for hg38\
\
`hg38_alternative_haplotypes.bed`  \
`hg38_centromeres.bed`  \
`hg38_fix_patches.bed`  \
`hg38_gaps.bed` \
Downloaded from UCSC table browser (https://genome.ucsc.edu/cgi-bin/hgTables)\
\
\
`k100_minus_exclusion_lists.mappable_regions.hg38.bed`  \
A bed file of all regions to be used for GC correction.  \
Contains all intervals with perfect mappability (mappability=1) for 100bp reads obtained from the UCSC genome browser (https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k100.Umap.MultiTrackMappability.bw)\
Additionally, all regions in `encode_unified_GRCh38_exclusion_list.bed` `hg38_alternative_haplotypes.bed` `hg38_centromeres.bed` `hg38_fix_patches.bed` and `hg38_gaps.bed` have been removed\
}