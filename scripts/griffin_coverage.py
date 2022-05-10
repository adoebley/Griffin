#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import sys
import argparse
import pandas as pd
import pysam
import pybedtools
import pyBigWig
import numpy as np
import time
import yaml 
from multiprocessing import Pool


# In[ ]:


# from matplotlib import pyplot as plt
# %matplotlib inline

# %load_ext autoreload
# %autoreload 2

# ##sample specific params for testing
# sample_name = 'MBC_1041_1_ULP'
# bam_path = '../../../../griffin_revisions_1/MBC_copy_bams/bam_file_copies/MBC_1041_1_ULP_recalibrated.bam'
# GC_bias_path = '../../../../griffin_revisions_1/GC_correction/MBC_ULP_GC_and_mappability_correction/results/GC_bias/MBC_1041_1_ULP.GC_bias.txt'

# mappability_bias_path = '../../../../griffin_revisions_1/GC_correction/MBC_ULP_GC_and_mappability_correction/results/mappability_bias/MBC_1041_1_ULP.mappability_bias.txt'
# mappability_correction = 'True'

# # mappability_bias_path = 'none'
# # mappability_correction = 'False'

# tmp_dir = 'tmp'

# ref_seq_path = '/fh/fast/ha_g/grp/reference/GRCh38/GRCh38.fa'
# mappability_bw='../../../../griffin_revisions_1/genome/k100.Umap.MultiTrackMappability.hg38.bw'
# chrom_sizes_path = '/fh/fast/ha_g/grp/reference/GRCh38/hg38.standard.chrom.sizes'

# # #additional params for testing
# sites_yaml = '/fh/fast/ha_g/user/adoebley/projects/griffin_revisions_1/MBC/CNA_correction_100kb_ATAC_np/config/sites.yaml'
# griffin_scripts_dir = '../'

# chrom_column = 'Chrom'
# position_column = 'position'
# strand_column = 'Strand'
# chroms = ['chr'+str(m) for m in np.arange(1,23)]

# norm_window = [-5000, 5000] 
# # norm_window = [-5000, 5000] 
# sz_range = [100, 200]
# map_q = 20

# number_of_sites = 500
# sort_by = 'Chrom'
# #sort_by = 'peak.count'
# ascending = 'False'

# # number_of_sites = 'none'
# # sort_by = 'none'
# # ascending = 'none'

# CPU = 6


# In[ ]:


parser = argparse.ArgumentParser()

parser.add_argument('--sample_name', help='name of sample', required=True)
parser.add_argument('--bam', help='bam file', required=True)
parser.add_argument('--GC_bias', help='GC bias info from griffin_GC_bias', required=True)
parser.add_argument('--mappability_bias', help='mappability bias info from griffin_mappability_bias', required=True)
parser.add_argument('--mappability_correction', help='True/False; whether to perform mappability correction', required=True)

parser.add_argument('--tmp_dir', help = 'directory for temporary outputs (may be large)', required=True)

parser.add_argument('--reference_genome',help = 'path to the reference genome',required=True)
parser.add_argument('--mappability_bw',help = 'bigWig file of genome wide mappability scores',required=True)
parser.add_argument('--chrom_sizes_path', help='path to chrom sizes file', required=True)

parser.add_argument('--sites_yaml', help='.bed file of sites', required=True)
parser.add_argument('--griffin_scripts_dir', help='path/to/scripts/', required=True)

parser.add_argument('--chrom_column',help='name of column containing chromosome number', default='Chrom')
parser.add_argument('--position_column',help='name of column containing chromosome position', default='Chrom')
parser.add_argument('--strand_column',help='name of column containing the strand (+ or -)', default='Strand')
parser.add_argument('--chroms', help='chroms to include when selecting sites', nargs='*', default=['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22'])

parser.add_argument('--norm_window',help='start and end of the window to be used for normalization',nargs=2, type=int, default=(-5000,5000))
parser.add_argument('--size_range',help='acceptable size range for fragments (to filter out genomic contamination)',nargs=2, type=int, default=(0,500))
parser.add_argument('--map_quality',help='minimum mapping quality', type=int, default=60)

parser.add_argument('--number_of_sites',help='number of sites to analyze', default='NA')
parser.add_argument('--sort_by',help='how to select the sites to analyze', default='none')
parser.add_argument('--ascending',help='whether to sort in ascending or descending order when selecting sites', default='NA')

parser.add_argument('--CPU',help='cpu available for parallelizing', type = int, required = True)


args = parser.parse_args()


sample_name = args.sample_name
bam_path = args.bam
GC_bias_path = args.GC_bias
mappability_bias_path = args.mappability_bias
mappability_correction = args.mappability_correction

tmp_dir = args.tmp_dir 

ref_seq_path = args.reference_genome
mappability_bw = args.mappability_bw
chrom_sizes_path = args.chrom_sizes_path

sites_yaml=args.sites_yaml
griffin_scripts_dir = args.griffin_scripts_dir

chrom_column = args.chrom_column
position_column=args.position_column
strand_column=args.strand_column
chroms = args.chroms

norm_window =args.norm_window
sz_range=args.size_range
map_q=args.map_quality

number_of_sites=args.number_of_sites
sort_by=args.sort_by
ascending=args.ascending

CPU = args.CPU


# In[ ]:


#print parameters for easy troubleshooting
if ascending.lower()=='false':
    ascending=False
elif ascending.lower()=='true':
    ascending=True
else:
    ascending='none'
    
print('\nparameters:')

print('\tsample_name = "'+sample_name+'"')
print('\tbam_path = "'+bam_path+'"')
print('\tGC_bias_path = "'+GC_bias_path+'"')
print('\tmappability_bias_path = "'+mappability_bias_path+'"')

print('\ttmp_dir = "'+tmp_dir+'"')

print('\tref_seq_path = "'+ref_seq_path+'"')
print('\tmappability_bw = "'+mappability_bw+'"')
print('\tchrom_sizes_path = "'+chrom_sizes_path+'"')

print('\tsites_yaml = "'+os.path.abspath(sites_yaml)+'"')
print('\tgriffin_scripts_dir = "'+griffin_scripts_dir+'"')

print('\tchrom_column = "'+chrom_column+'"')
print('\tposition_column = "'+position_column+'"')
print('\tstrand_column = "'+strand_column+'"')
print('\tchroms = ',chroms)

print('\tnorm_window = ', norm_window)
print('\tsz_range =',sz_range)
print('\tmap_q =',map_q)

print('\tnumber_of_sites = "'+str(number_of_sites)+'"')
print('\tsort_by = "'+sort_by+'"')
print('\tascending = "'+str(ascending)+'"')

print('\tCPU =',CPU)
print('\n')
sys.stdout.flush()


# In[ ]:


overall_start_time = time.time()


# In[ ]:


if not mappability_correction.lower() == 'true':
    del(mappability_bw,mappability_bias_path)
    print('Skipping mappability correction')
else:
    print('Correcting mappability')


# In[ ]:


#define global parameters and open global files
########################################
#GET GC BIAS
########################################
#open the GC_bias file 
GC_bias = pd.read_csv(GC_bias_path, sep='\t')

#get rid of extremely low GC bias values
#these fragments will now be excluded 
#these fragments are extremely rare so it is difficult to get a good estimate of GC bias
GC_bias['smoothed_GC_bias'] = np.where(GC_bias['smoothed_GC_bias']<0.05,np.nan,GC_bias['smoothed_GC_bias'])

GC_bias = GC_bias[['length','num_GC','smoothed_GC_bias']]
GC_bias = GC_bias.set_index(['num_GC','length']).unstack()

#convert to a dictionary
GC_bias = GC_bias.to_dict()

#get rid of values where the num_GC is greater than the length (included due to the way I made the dict)
GC_bias2 = {}
for key in GC_bias.keys():
    length = key[1]
    GC_bias2[length] = {}
    for num_GC in range(0,length+1):
        bias = GC_bias[key][num_GC]
        GC_bias2[length][num_GC]=bias
GC_bias = GC_bias2 
del(GC_bias2)

if mappability_correction.lower() == 'true':
    ########################################
    #GET mappability bias
    ########################################
    mappability_bias = pd.read_csv(mappability_bias_path, sep='\t')
    mappability_bias['mappable_percent'] = np.round(mappability_bias['mappable_percent']).astype(int) #convert indexes to integers
    mappability_bias = mappability_bias[['mappable_percent','smoothed_map_bias']].set_index('mappable_percent').to_dict()['smoothed_map_bias']
    print('using smoothed map bias')
    def closest_key(dictionary, i):
        sorted_keys=np.array(sorted(dictionary.keys()))
        idx = (np.abs(sorted_keys - i)).argmin() 
        closest_key=sorted_keys[idx]
        return(closest_key)

    #depending on the read length, some values might be missing, add them
    for i in range(0,101):
        if i in mappability_bias:
            pass
        else:
            mappability_bias[i] = mappability_bias[closest_key(mappability_bias,i)]
            print('adding mappabilty',i,'to dict')


# In[ ]:


#snakemake should create these folders, but if not using the snakemake, this is needed
tmp_sample_dir = tmp_dir+'/'+sample_name
if not os.path.exists(tmp_sample_dir): 
    os.mkdir(tmp_sample_dir)

tmp_pybedtools = tmp_sample_dir+'/tmp_pybedtools'
if not os.path.exists(tmp_pybedtools): 
    os.mkdir(tmp_pybedtools)
pybedtools.set_tempdir(tmp_pybedtools)

tmp_bigWig = tmp_sample_dir+'/tmp_bigWig'
if not os.path.exists(tmp_bigWig): 
    os.mkdir(tmp_bigWig)


# In[ ]:


#import the griffin scripts
sys.path.insert(0, griffin_scripts_dir)
import griffin_functions


# In[ ]:


#import the site_lists
with open(sites_yaml,'r') as f:
    sites = yaml.safe_load(f)
sites = sites['site_lists']

all_sites = pd.DataFrame()
for site_name in sites.keys():
    site_file = sites[site_name]
    current_sites = griffin_functions.import_and_filter_sites(site_name,site_file,strand_column,chrom_column,position_column,chroms,ascending,sort_by,number_of_sites)
    all_sites = all_sites.append(current_sites, ignore_index=True).copy()
    sys.stdout.flush()


# In[ ]:


#number of bp to fetch upstream and downstream of the site
upstream_bp = norm_window[0]-sz_range[0] #this should be negative
downstream_bp = norm_window[1]+sz_range[0] #this should be positive
all_sites = griffin_functions.define_fetch_interval('Total sites',all_sites,chrom_column,position_column,
                                                    chroms,chrom_sizes_path,upstream_bp,downstream_bp)


# In[ ]:


#convert to pybedtools and merge overlapping segments
start_time = time.time()
all_sites_bed = pybedtools.BedTool.from_dataframe(all_sites[[chrom_column,'fetch_start','fetch_end']])
all_sites_bed = all_sites_bed.sort()
all_sites_bed = all_sites_bed.merge()
print('Intervals to fetch:\t'+str(len(all_sites_bed)))
print('Total bp to fetch:\t'+str(all_sites_bed.total_coverage()))
sys.stdout.flush()

#split the long intervals 
to_fetch = all_sites_bed.to_dataframe()
to_fetch['length'] = to_fetch['end']-to_fetch['start']

print('Max fetch length: '+str(to_fetch['length'].max())+' bp')

to_fetch = to_fetch[['chrom','start','end']]
to_fetch = to_fetch.sort_values(by=['chrom','start']).reset_index(drop=True)
to_fetch = to_fetch.reset_index() #add an index column
split_len = pybedtools.BedTool.from_dataframe(to_fetch[['chrom','start','end']]).total_coverage()

# print('Intervals to fetch (after splitting long intervals):\t'+str(len(to_fetch)))
# print('Total bp to fetch (after splitting long intervals):\t'+str(split_len))
sys.stdout.flush()
pybedtools.cleanup(verbose=False, remove_all=True)


# In[ ]:


def collect_fragments(input_list):
    i,chrom,start,end = input_list
    #open the bam file for each pool worker (otherwise individual pool workers can close it)
    bam_file = pysam.AlignmentFile(bam_path)
    
    #open the ref seq
    ref_seq=pysam.FastaFile(ref_seq_path)

    #make dicts to hold the fetched positions
    columns = np.arange(start,end,1)
    cov_dict={m:0 for m in columns} 
    GC_cov_dict={m:0 for m in columns} 
    
    if mappability_correction.lower() == 'true':
        mappability = pyBigWig.open(mappability_bw)    
        GC_map_cov_dict={m:0 for m in columns}   
        
    #fetch reads
    fetched=bam_file.fetch(contig=chrom, start=start, stop=end) #fetch reads that map to the region of interest
                
    ########################
    #count coverage
    ########################
    for read in fetched:
        #filter out reads
        if abs(read.template_length)>=sz_range[0] and abs(read.template_length)<=sz_range[1]         and read.is_paired==True and read.mapping_quality>=map_q and read.is_duplicate==False and read.is_qcfail==False:
            #only use fw reads with positive fragment lengths (negative indicates an abnormal pair)
            #all paired end reads have a fw and rv read so we don't need the rv read to find the midpoint.
            if read.is_reverse==False and read.template_length>0:
                fragment_start = read.reference_start #for fw read, read start is fragment start
                fragment_end = read.reference_start+read.template_length
                midpoint = int(np.floor((fragment_start+fragment_end)/2))
                                
                #count the GC content
                fragment_seq = ref_seq.fetch(read.reference_name,fragment_start,fragment_end)
                fragment_seq = np.array(list(fragment_seq.upper()))
                fragment_seq[np.isin(fragment_seq, ['A','T','W'])] = 0
                fragment_seq[np.isin(fragment_seq, ['C','G','S'])] = 1
                rng = np.random.default_rng(fragment_start)
                fragment_seq[np.isin(fragment_seq, ['N','R','Y','K','M','B','D','H','V'])] = rng.integers(2, size=len(fragment_seq[np.isin(fragment_seq, ['N','R','Y','K','M','B','D','H','V'])])) #random integer in range(2) (i.e. 0 or 1)
                fragment_seq = fragment_seq.astype(float)
    
                if mappability_correction.lower() == 'true':
                    #find the two read locations for mappability correction
                    fw_read_map = mappability.values(chrom,read.reference_start,read.reference_start+read.reference_length)
                    fw_read_map = np.mean(np.nan_to_num(fw_read_map)) #replace any nan with zero and take the mean

                    rv_read_map = mappability.values(chrom,read.reference_start+read.template_length-read.reference_length,read.reference_start+read.template_length)
                    rv_read_map = np.mean(np.nan_to_num(rv_read_map)) #replace any nan with zero and take the mean
                                
                #check that the site is in the window          
                if midpoint>=start and midpoint<end:
                    #count the fragment
                    cov_dict[midpoint]+=1

                    ##get the GC bias
                    read_GC_content = sum(fragment_seq)
                    read_GC_bias = GC_bias[abs(read.template_length)][read_GC_content]

                    #count the fragment weighted by GC bias
                    if not np.isnan(read_GC_bias):
                        GC_cov_dict[midpoint]+=(1/read_GC_bias)
                        
                    if mappability_correction.lower() == 'true':
                        #get the mappability bias
                        read_map = np.int32(np.round(100*(fw_read_map+rv_read_map)/2))
                        read_map_bias = mappability_bias[read_map]
                        GC_map_cov_dict[midpoint]+=(1/read_GC_bias)*(1/read_map_bias)
                        
                    #print(read_GC_bias,read_map,read_map_bias)
                    
                else: #if fragment doesn't fully overlap
                    continue
                    
                del(read,midpoint,fragment_seq)
                
            else:
                #print('reverse',read.is_reverse)
                continue
                
    output = pd.DataFrame(pd.Series(cov_dict, name = 'uncorrected'))
    output['GC_corrected'] = pd.Series(GC_cov_dict)
    if mappability_correction.lower() == 'true':
        output['GC_map_corrected'] = pd.Series(GC_map_cov_dict)
    output = output[output['uncorrected']>0] #don't waste memory on positions with no coverage
    output['chrom'] = chrom
    output = output.reset_index().rename(columns = {'index':'position'})
    output['uncorrected'] = output['uncorrected'].astype(float)
    output['GC_corrected'] = np.round(output['GC_corrected'],5)
    if mappability_correction.lower() == 'true':
        output['GC_map_corrected'] = np.round(output['GC_map_corrected'],5)
        mappability.close()

    bam_file.close()
    ref_seq.close()
    
    
    if (i+1)%1000==0:
        printout = griffin_functions.progress_report([chrom,start,end],'intervals',start_time,time.time(),i,len(to_fetch))
        print(printout)
        sys.stdout.flush()
    
    return(output)


# In[ ]:


#run the analysis 
print('Starting fetch')
sys.stdout.flush()
start_time = time.time()

p = Pool(processes=CPU) #use the specified number of processes
results = p.map(collect_fragments, to_fetch.values, 1) #Send only one interval to each processor at a time.

elapsed_time = time.time()-overall_start_time
print('Done with fetch '+str(int(np.floor(elapsed_time/60)))+' min '+str(int(np.round(elapsed_time%60)))+' sec')
del(elapsed_time)
sys.stdout.flush()


# In[ ]:


print('Starting export')
sys.stdout.flush()
start_time = time.time()

chrom_sizes = pd.read_csv(chrom_sizes_path, sep='\t', header=None)
chrom_sizes = chrom_sizes[chrom_sizes[0].isin(chroms)]

uncorrected_bw = pyBigWig.open(tmp_bigWig+"/"+sample_name+".uncorrected.bw", "w")
GC_bw = pyBigWig.open(tmp_bigWig+"/"+sample_name+".GC_corrected.bw", "w")

uncorrected_bw.addHeader([(a,b) for a,b in chrom_sizes.values])
GC_bw.addHeader([(a,b) for a,b in chrom_sizes.values])

if mappability_correction.lower() == 'true':
    GC_map_bw = pyBigWig.open(tmp_bigWig+"/"+sample_name+".GC_map_corrected.bw", "w")
    GC_map_bw.addHeader([(a,b) for a,b in chrom_sizes.values])

for i in range(len(results)):
    current = results[i]
    if len(current)>0:
        if np.nansum(current['uncorrected'])>0:
            uncorrected_bw.addEntries(list(current['chrom']), list(current['position']), ends = list(current['position']+1), values = list(current['uncorrected']))  
        else:
            print('no uncorrected reads:')
            print(current)
            
        if np.nansum(current['GC_corrected'])>0:
            GC_bw.addEntries(list(current['chrom']), list(current['position']), ends = list(current['position']+1), values = list(current['GC_corrected'])) 
        else:
            print('no GC corrected reads:')
            print(current)
        
        if mappability_correction.lower() == 'true':
            if np.nansum(current['GC_map_corrected'])>0:   
                GC_map_bw.addEntries(list(current['chrom']), list(current['position']), ends = list(current['position']+1), values = list(current['GC_map_corrected']))  
            else:
                print('no GC map corrected reads:')
                print(current)
            
    if (i+1)%10000==0 and len(current)>0:
        printout = griffin_functions.progress_report(list(current.iloc[0][['chrom','position']]),'intervals',start_time,time.time(),i,len(results))
        print(printout)
        sys.stdout.flush()

if mappability_correction.lower() == 'true':
    GC_map_bw.close()     
    
uncorrected_bw.close()
GC_bw.close()


elapsed_time = time.time()-overall_start_time
print('Done with export '+str(int(np.floor(elapsed_time/60)))+' min '+str(int(np.round(elapsed_time%60)))+' sec')
sys.stdout.flush()


# In[ ]:




