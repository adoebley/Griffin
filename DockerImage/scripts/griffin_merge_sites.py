#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import sys
import argparse
import pandas as pd
import pysam
import pyBigWig
import numpy as np
import time
from scipy.signal import savgol_filter
from scipy.stats import zscore
import yaml 
from multiprocessing import Pool
import pybedtools


# In[ ]:


# from matplotlib import pyplot as plt
# %matplotlib inline

# #sample specific params for testing
# sample_name = 'MBC_1041_1_ULP'
# uncorrected_bw_path = 'tmp/MBC_1041_1_ULP/tmp_bigWig/MBC_1041_1_ULP.uncorrected.bw'
# GC_corrected_bw_path = 'tmp/MBC_1041_1_ULP/tmp_bigWig/MBC_1041_1_ULP.GC_corrected.bw'
# GC_map_corrected_bw_path = 'tmp/MBC_1041_1_ULP/tmp_bigWig/MBC_1041_1_ULP.GC_map_corrected.bw'
# mappability_correction = 'True'

# # GC_map_corrected_bw_path = 'none'
# # mappability_correction = 'False'

# # sample_name = 'HD45.ctDNA.WGS.FC19269447'
# # uncorrected_bw_path = 'tmp/HD45.ctDNA.WGS.FC19269447/tmp_bigWig/HD45.ctDNA.WGS.FC19269447.uncorrected.bw'
# # GC_corrected_bw_path = 'tmp/HD45.ctDNA.WGS.FC19269447/tmp_bigWig/HD45.ctDNA.WGS.FC19269447.GC_corrected.bw'
# # GC_map_corrected_bw_path = 'tmp/HD45.ctDNA.WGS.FC19269447/tmp_bigWig/HD45.ctDNA.WGS.FC19269447.GC_map_corrected.bw'

# tmp_dir = 'tmp'
# results_dir = 'results'

# mappability_bw='../../../../griffin_revisions_1/genome/k100.Umap.MultiTrackMappability.hg38.bw'
# chrom_sizes_path = '/fh/fast/ha_g/grp/reference/GRCh38/hg38.standard.chrom.sizes'

# # #additional params for testing
# sites_yaml = '/fh/fast/ha_g/user/adoebley/projects/griffin_revisions_1/MBC/CNA_correction_100kb_ATAC_np/config/sites.yaml'
# # sites_yaml = '/fh/fast/ha_g/user/adoebley/projects/nucleosome_profiling_griffin/add_mappability_3/sites/test_sites.yaml'

# griffin_scripts_dir = '../'

# chrom_column = 'Chrom'
# position_column = 'position'
# strand_column = 'Strand'
# chroms = ['chr'+str(m) for m in np.arange(1,23)]

# # norm_window = [-250000, 250000]
# norm_window = [-5000, 5000] #for testing

# save_window = [-1000, 1000]#for testing
# center_window = [-30,30] #define the center of the interval for feature calculation
# fft_window = [-960,960]
# fft_index = 10
# smoothing_length = 167 #fragment_length

# encode_exclude = '../../../../griffin_revisions_1/genome/encode_unified_GRCh38_exclusion_list.bed'
# centromere_path = '../../../../griffin_revisions_1/genome/hg38_centromeres.bed'
# gap_path = '../../../../griffin_revisions_1/genome/hg38_gaps.bed'
# patch_path = '../../../../griffin_revisions_1/genome/hg38_fix_patches.bed'
# alternative_haplotype_path = '../../../../griffin_revisions_1/genome/hg38_alternative_haplotypes.bed'

# exclude_paths = [encode_exclude,centromere_path,gap_path,patch_path,alternative_haplotype_path]
# del(encode_exclude,centromere_path,gap_path,patch_path)

# # exclude_paths = ['none']

# step = 15

# # CNA_normalization = 'False'
# CNA_normalization = 'False'

# individual = 'False'
# smoothing = 'True'

# exclude_outliers_parameter = 'True'
# exclude_zero_mappability_parameter = 'True'

# # number_of_sites = 'none'
# # sort_by = 'none'
# # ascending = 'none'

# number_of_sites = 500
# sort_by = 'Chrom'
# #sort_by = 'peak.count'
# ascending = 'False'

# CPU = 1


# In[ ]:


parser = argparse.ArgumentParser()

parser.add_argument('--sample_name', help='name of sample', required=True)
parser.add_argument('--uncorrected_bw_path', help='uncorrected bigWig from griffin_coverage', required=True)
parser.add_argument('--GC_corrected_bw_path', help='GC_corrected bigWig from griffin_coverage', required=True)
parser.add_argument('--GC_map_corrected_bw_path', help='GC_corrected bigWig from griffin_coverage', required=True)
parser.add_argument('--mappability_correction', help='True/False; whether to perform mappability correction', required=True)

parser.add_argument('--tmp_dir', help = 'directory for temporary outputs (may be large)', required=True)
parser.add_argument('--results_dir', help = 'directory for results', required=True)

parser.add_argument('--mappability_bw',help = 'bigWig file of genome wide mappability scores',required=True)
parser.add_argument('--chrom_sizes_path', help='path to chrom sizes file', required=True)

parser.add_argument('--sites_yaml', help='.bed file of sites', required=True)
parser.add_argument('--griffin_scripts_dir', help='path/to/scripts/', required=True)

parser.add_argument('--chrom_column',help='name of column containing chromosome number', default='Chrom')
parser.add_argument('--position_column',help='name of column containing chromosome position', default='Chrom')
parser.add_argument('--strand_column',help='name of column containing the strand (+ or -)', default='Strand')
parser.add_argument('--chroms', help='chroms to include when selecting sites', nargs='*', default=['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22'])

parser.add_argument('--norm_window',help='start and end of the window to be used for normalization',nargs=2, type=int, default=(-5000,5000))
parser.add_argument('--save_window',help='start and end of the window to be saved in the outputs',nargs=2, type=int, default=(-1000,1000))
parser.add_argument('--center_window',help='start and end of the window to be used for calculating the central coverage feature',nargs=2, type=int, default=(-1000,1000))
parser.add_argument('--fft_window',help='start and end of the window to be used for calculating the amplitude feature (default is for bin size 15)',nargs=2, type=int, default=(-960,960))
parser.add_argument('--fft_index',help='index of the fft component to be saved as the amplitude feature (default is for 15bp bins and -960 to 960 fft_window)', type = int, default=10)
parser.add_argument('--smoothing_length',help='window length for the Savitzky-Golay smoothing function (should be approximately the mean fragment length)', type = int, default=165)

parser.add_argument('--exclude_paths', help='path to bed files of regions to filter out (excluded regions, centromeres, gaps, patches, alternative haplotypes), or "none" to not exclude any regions', required=True, nargs = '*')

parser.add_argument('--step',help='step size when calculating coverage', type=int, default=5)

parser.add_argument('--CNA_normalization',help='whether to normalize each site individually to the copy number within the normalization window',default='False', required = True)
parser.add_argument('--individual',help='save individual site coverage. TRUE WILL RESULT IN HUGE OUTPUT FILES. (True/False)',default='False', required = True)
parser.add_argument('--smoothing',help='whether to use a savgol filter to smooth sites (True/False)', default = 'True', required = True)

parser.add_argument('--exclude_outliers',help='whether to exclude bins with extreme outlier coverage >10SD above the mean (True/False)', default='True', required = True)
parser.add_argument('--exclude_zero_mappability',help='whether to exclude bins with zero mappability (True/False)', default='True', required = True)

parser.add_argument('--number_of_sites',help='number of sites to analyze', default='none')
parser.add_argument('--sort_by',help='how to select the sites to analyze', default='none')
parser.add_argument('--ascending',help='whether to sort in ascending or descending order when selecting sites', default='none')

parser.add_argument('--CPU',help='CPU available for parallelizing', type = int, required = True)

args = parser.parse_args()


sample_name = args.sample_name
uncorrected_bw_path = args.uncorrected_bw_path
GC_corrected_bw_path = args.GC_corrected_bw_path
GC_map_corrected_bw_path = args.GC_map_corrected_bw_path
mappability_correction = args.mappability_correction

tmp_dir = args.tmp_dir
results_dir = args.results_dir

mappability_bw = args.mappability_bw
chrom_sizes_path = args.chrom_sizes_path

sites_yaml = args.sites_yaml
griffin_scripts_dir = args.griffin_scripts_dir
chrom_column = args.chrom_column
position_column = args.position_column
strand_column = args.strand_column
chroms = args.chroms

norm_window = args.norm_window
save_window = args.save_window
center_window = args.center_window
fft_window = args.fft_window
fft_index = args.fft_index
smoothing_length = args.smoothing_length

exclude_paths = args.exclude_paths

step = args.step

CNA_normalization = args.CNA_normalization
individual = args.individual
smoothing = args.smoothing

exclude_outliers_parameter = args.exclude_outliers
exclude_zero_mappability_parameter = args.exclude_zero_mappability

number_of_sites = args.number_of_sites
sort_by = args.sort_by
ascending = args.ascending

CPU = args.CPU


# In[ ]:


overall_start_time = time.time()


# In[ ]:


if not mappability_correction.lower() == 'true':
    del(GC_map_corrected_bw_path)
    print('Skipping mappability correction')
    
    results_dict_template = {'uncorrected': {'input_path':uncorrected_bw_path},
                'GC_corrected': {'input_path':GC_corrected_bw_path}}
else:
    print('Correcting mappability')
    
    results_dict_template = {'uncorrected': {'input_path':uncorrected_bw_path},
                    'GC_corrected': {'input_path':GC_corrected_bw_path},
                    'GC_map_corrected': {'input_path':GC_map_corrected_bw_path}}


# In[ ]:


norm_window=[int(np.ceil(norm_window[0]/step)*step),int(np.floor(norm_window[1]/step)*step)] #round to the nearest step inside the window
save_window=[int(np.ceil(save_window[0]/step)*step),int(np.floor(save_window[1]/step)*step)] #round to the nearest step inside the window
center_window=[int(np.ceil(center_window[0]/step)*step),int(np.floor(center_window[1]/step)*step)] #round to the nearest step inside the window
fft_window=[int(np.ceil(fft_window[0]/step)*step),int(np.floor(fft_window[1]/step)*step)] #round to the nearest step inside the window

all_positions = np.arange(norm_window[0],norm_window[1])
norm_columns = np.arange(norm_window[0],norm_window[1],step)
save_columns = np.arange(save_window[0],save_window[1],step)
center_columns = np.arange(center_window[0],center_window[1],step)
fft_columns = np.arange(fft_window[0],fft_window[1],step)

smoothing_length=int(np.round(smoothing_length/step)*step) #round fragment length to the nearest step

if ascending.lower()=='false':
    ascending=False
elif ascending.lower()=='true':
    ascending=True
else:
    ascending='none'
    
print('Normalization window (rounded down to step):', norm_window)
print('Saving window (rounded down to step):', save_window)
print('Center window (rounded down to step):', center_window)
print('FFT window (rounded down to step):', fft_window)
print('Savgol filter smoothing window:', smoothing_length)
print('Ascending is:',ascending)
sys.stdout.flush()


# In[ ]:


print('Excluding regions:',exclude_paths)
print('Excluding bins with coverage outliers:',exclude_outliers_parameter)
print('Excluding bins with zero mappability:',exclude_zero_mappability_parameter)


# In[ ]:


#snakemake should create these folders, but if not using the snakemake, this is needed
tmp_sample_dir = tmp_dir+'/'+sample_name
if not os.path.exists(tmp_sample_dir): 
    os.mkdir(tmp_sample_dir)

tmp_pybedtools = tmp_sample_dir+'/tmp_pybedtools2'
if not os.path.exists(tmp_pybedtools): 
    os.mkdir(tmp_pybedtools)
pybedtools.set_tempdir(tmp_pybedtools)

tmp_bigWig = tmp_sample_dir+'/tmp_bigWig'
if not os.path.exists(tmp_bigWig): 
    os.mkdir(tmp_bigWig)

#make results dir
results_sample_dir = results_dir+'/'+sample_name
if not os.path.exists(results_sample_dir): 
    os.mkdir(results_sample_dir)


# In[ ]:


#import the griffin scripts
sys.path.insert(0, griffin_scripts_dir)
import griffin_functions


# In[ ]:


if exclude_paths==['none']:
    print('No excluded regions.')

else: #if there are regions to be excluded
    #get the excluded regions
    merged_exclude_regions = pybedtools.BedTool('\n', from_string=True)

    #create an empty bed file
    excluded_regions_bw = pyBigWig.open(tmp_bigWig+"/excluded_regions.bw", "w")
    chrom_sizes = pd.read_csv(chrom_sizes_path, sep='\t', header=None)
    chrom_sizes = chrom_sizes[chrom_sizes[0].isin(chroms)]
    excluded_regions_bw.addHeader([(a,b) for a,b in chrom_sizes.values])

    for path in exclude_paths:
        print('excluding:',path)
        current_regions = pybedtools.BedTool(path)
        merged_exclude_regions = merged_exclude_regions.cat(current_regions)    
        del(current_regions)
    merged_exclude_regions = merged_exclude_regions.to_dataframe()
    merged_exclude_regions = merged_exclude_regions[merged_exclude_regions['chrom'].isin(chroms)]
    pybedtools.cleanup()
    excluded_regions_bw.addEntries(list(merged_exclude_regions['chrom']), list(merged_exclude_regions['start']), ends = list(merged_exclude_regions['end']), values = [1.0 for m in range(len(merged_exclude_regions))])  

    excluded_regions_bw.close()


# In[ ]:


#import the site_lists
with open(sites_yaml,'r') as f:
    sites = yaml.safe_load(f)
sites = sites['site_lists']
print('Analyzing '+str(len(sites))+' site lists')


# In[ ]:


def fetch_bw_values(bw_path,current_sites,site_name,name):
    elapsed_time = time.time()-overall_start_time
    print(site_name+' '+name+' starting fetch '+str(int(np.floor(elapsed_time/60)))+' min '+str(int(np.round(elapsed_time%60)))+' sec')
    del(elapsed_time)
    sys.stdout.flush()
    
    bw = pyBigWig.open(bw_path)
    
    fw_markers = ['+',1,'1']
    rv_markers = ['-',-1,'-1']

    results = pd.DataFrame(np.zeros([len(current_sites),norm_window[1]-norm_window[0]]))
    start_time = time.time()
        
    for i in range(len(current_sites)):
        chrom,start,end,strand = current_sites.iloc[i][[chrom_column,'fetch_start','fetch_end',strand_column]]

        values = bw.values(chrom, start, end, numpy=True)
        values = np.nan_to_num(values) #turn nan into zero because bw doesn't store zero
        
        #if the window extends beyond the end of the chromosome add np.nan to fill the gap
        if len(values)<(norm_window[1]-norm_window[0]):
            ###################
            position = current_sites.iloc[i][position_column]
            array_start = position + norm_window[0]
            array_end = position + norm_window[1]

            temp_series = pd.Series(np.full(norm_window[1]-norm_window[0], np.nan), index = np.arange(array_start,array_end))
            temp_series[np.arange(start,end)] = values
            
            #print('too_short',i,len(values),norm_window[1]-norm_window[0])
            #print(chrom,start,end,strand)
            #print(values)
            
            values = temp_series.values
            del(temp_series, position, array_start,array_end)
            ###################
        if strand in rv_markers:
            values = values[::-1]
        results.iloc[i] =  pd.Series(values)
        
        if (i+1)%10000==0:
            printout = griffin_functions.progress_report([site_name,name,chrom,start,end],'intervals',start_time,time.time(),i,len(current_sites))
            print(printout,', size',np.round(sys.getsizeof(results)/(1024**3),2),'GB')
            sys.stdout.flush()
    results.columns = np.arange(norm_window[0],norm_window[1])
    
    printout = griffin_functions.progress_report([site_name,name,'fetch_complete'],'intervals',start_time,time.time(),i,len(current_sites))
    print(printout,', size',np.round(sys.getsizeof(results)/(1024**3),2),'GB')
    sys.stdout.flush()
    return(results)


# In[ ]:


def sum_bins(results,name):
    #If a bin has an np.nan value, the whole bin will become np.nan
    summed = np.sum(results.values.reshape(len(results),int(len(results.columns)/step), step), axis = 2)    
    summed = pd.DataFrame(summed)
    summed.columns = norm_columns
     
    return(summed)


# In[ ]:


def exclude_regions(results,excluded_regions,name):
    results = np.where(excluded_regions>0,np.nan,results)#if any bp in the bin were excluded (1) exclude the bin
    results = pd.DataFrame(results)
    results.columns = norm_columns
    return(results)


# In[ ]:


def exclude_zero_mappability(results,mappability_values,name):
    results = np.where(mappability_values>0,results,np.nan)#only retain positions where the mappability is >0
    results = pd.DataFrame(results)
    results.columns = norm_columns
    return(results)


# In[ ]:


def make_outlier_mask(results,site_name):
    max_value = results.max().max()
    min_cutoff = 2 #minimum coverage that must be retained even if it is an outlier
    print(site_name,'max_bin_coverage is',max_value, 'midpoints')
    
    scores = pd.DataFrame(zscore(results.values, axis = None, nan_policy='omit'))
    outlier_mask = pd.DataFrame(np.where(scores<10,1,np.nan))
    outlier_mask.columns = norm_columns
    
    
    if (results*outlier_mask).max().max()<min_cutoff:
        print(site_name, 'low coverage, resetting the outlier cutoff to '+str(min_cutoff))
        outlier_mask = pd.DataFrame(np.where(results<=min_cutoff,1,np.nan))
        outlier_mask.columns = norm_columns
        
    outlier_cutoff = (results*outlier_mask).max().max()
    print(site_name,'masking sites with  >', outlier_cutoff, 'midpoints')

    return(outlier_mask,outlier_cutoff)


# In[ ]:


def normalize_and_smooth(results,site_name,name):
    #get the mean midpoints per valid position in each site
    mean_reads_per_bp_in_normalization_window = np.nanmean(results[norm_columns],axis = 1)/step
    mean_reads_per_bp_in_saved_window = np.nanmean(results[save_columns],axis = 1)/step
        
    #normalize individual sites to 1 to remove CNA
    if CNA_normalization.lower() == 'true':
        print(site_name,name,'normalizing CNAs')
        mean_data = np.nanmean(results.values,axis = 1, keepdims=True)
        #replace zero with nan so there aren't any infinities in the output
        mean_data = np.where(mean_data==0,np.nan,mean_data)
        results[norm_columns] = results[norm_columns]/mean_data  
        
    #take the mean of all sites
    if not individual.lower()=='true':
        print(site_name,name,'averaging sites')
        results = pd.DataFrame(pd.Series(np.nanmean(results[norm_columns], axis = 0), index=norm_columns)).T
        results.columns = norm_columns
        mean_reads_per_bp_in_normalization_window = np.nanmean(mean_reads_per_bp_in_normalization_window)
        mean_reads_per_bp_in_saved_window = np.nanmean(mean_reads_per_bp_in_saved_window)
        
    #smooth the sites
    if smoothing.lower()=='true':
        print(site_name,name,'smoothing')
        #savgol window should be approx one fragment length but it must be odd
        savgol_window=np.floor(smoothing_length/step)
        if savgol_window%2==0:
            savgol_window=savgol_window+1
        savgol_window=int(savgol_window)

        results[norm_columns] = savgol_filter(results[norm_columns], savgol_window, 3)
    
    #normalize the average site to 1
    print(site_name,name,'correcting for read depth')
    mean_value = np.nanmean(results[norm_columns])
    results[norm_columns] = results[norm_columns]/mean_value
    
    #save only plot columns 
    results = results[save_columns].copy()
    
    results['mean_reads_per_bp_in_normalization_window'] = mean_reads_per_bp_in_normalization_window
    results['mean_reads_per_bp_in_saved_window'] = mean_reads_per_bp_in_saved_window
    
    return(results)


# In[ ]:


def calculate_features(results):
    results['mean_coverage'] = results[save_columns].mean(axis = 1, skipna=False) #pandas skips na by default
    results['central_coverage'] = results[center_columns].mean(axis = 1, skipna=False)
    fft_res = np.fft.fft(results[fft_columns])
    results['amplitude'] = np.abs(fft_res[:,fft_index])
    
    return(results)


# In[ ]:


def merge_sites(input_list):
    site_name,site_file = input_list
    #get the site lists and define the fetch interval
    current_sites = griffin_functions.import_and_filter_sites(site_name,site_file,strand_column,chrom_column,position_column,chroms,ascending,sort_by,number_of_sites)   
    current_sites = griffin_functions.define_fetch_interval(site_name,current_sites,chrom_column,position_column,chroms,chrom_sizes_path,norm_window[0],norm_window[1])

    sys.stdout.flush()

    #dict to hold results 
    results_dict = results_dict_template.copy()

    #fetch coverage and sum into bins of length step
    for key in results_dict.keys():
        results_dict[key]['coverage'] = fetch_bw_values(results_dict[key]['input_path'],current_sites,site_name,key)
        results_dict[key]['coverage'] = sum_bins(results_dict[key]['coverage'],key)
    
    #exclude specified regions
    if not exclude_paths == ['none']:
        print(site_name+' - excluding specified regions.')
        regions_to_exclude = fetch_bw_values(tmp_bigWig+"/excluded_regions.bw",current_sites,site_name,'to_exclude')
        regions_to_exclude = sum_bins(regions_to_exclude,'to_exclude')
        for key in results_dict.keys():
            results_dict[key]['coverage'] = exclude_regions(results_dict[key]['coverage'],regions_to_exclude,key)
        del(regions_to_exclude)

    #exclude zero mappability
    if exclude_zero_mappability_parameter.lower()=='true':
        print(site_name+' - excluding zero mappability.')
        #fetch excluded_regions
        mappability_values = fetch_bw_values(mappability_bw,current_sites,site_name,'mappabilty')  
        #replace zero with np.nan for mappability
        #when summing bins, any bin with one or more zeros will now be np.nan
        mappability_values[all_positions] = np.where(mappability_values[all_positions]==0,np.nan,mappability_values[all_positions])
        mappability_values = sum_bins(mappability_values,'mappability')
        
        for key in results_dict.keys():
            results_dict[key]['coverage'] = exclude_zero_mappability(results_dict[key]['coverage'],mappability_values,key)
        del(mappability_values)

    #mask out bins with coverage >10 SD above the mean
    if exclude_outliers_parameter.lower()=='true':
        print(site_name+' - excluding outliers.')
        outlier_mask,outlier_cutoff = make_outlier_mask(results_dict['uncorrected']['coverage'],site_name)
        sys.stdout.flush()

        for key in results_dict.keys():
            results_dict[key]['coverage'] = results_dict[key]['coverage']*outlier_mask
    else:
        outlier_cutoff='NA'

    #normalize to a mean of 1 and smooth
    for key in results_dict.keys():
        results_dict[key]['coverage'] = normalize_and_smooth(results_dict[key]['coverage'],site_name,key)
    sys.stdout.flush()

    #get features
    for key in results_dict.keys():
        results_dict[key]['coverage'] = calculate_features(results_dict[key]['coverage'])

    #get metadata
    for key in results_dict.keys():
        results_dict[key]['coverage']['outlier_cutoff']=outlier_cutoff
        results_dict[key]['coverage']['exclude_zero_mappability']=exclude_zero_mappability_parameter 
        results_dict[key]['coverage']['correction'] = key
        results_dict[key]['coverage']['number_of_sites']=len(current_sites)
        results_dict[key]['coverage']['site_name']=site_name
        results_dict[key]['coverage']['smoothing']=smoothing
        results_dict[key]['coverage']['CNA_normalization']=CNA_normalization
        results_dict[key]['coverage']['sample']=sample_name    
        results_dict[key]['coverage'] = results_dict[key]['coverage'].copy()
        
    #if saving individual sites, keep the locations
    if individual.lower()=='true':
        current_sites = current_sites.drop(columns = ['site_name'])
        for key in results_dict.keys():
            results_dict[key]['coverage'] = results_dict[key]['coverage'].merge(current_sites, left_index=True, right_index=True, validate = 'one_to_one')
        
    elapsed_time = time.time()-overall_start_time
    print(site_name+' merge complete '+str(int(np.floor(elapsed_time/60)))+' min '+str(int(np.round(elapsed_time%60)))+' sec')
    del(elapsed_time)
    sys.stdout.flush()

    return(results_dict)


# In[ ]:


# # for testing
# # for site_name in sites.keys()[1]:
# for site_name in [list(sites.keys())[1]]:
#     site_name = site_name
#     site_file = sites[site_name]
#     results = merge_sites([site_name,site_file])
#     break


# In[ ]:


#run the analysis 
to_do_list = [[key,sites[key]] for key in sites.keys()]

p = Pool(processes=CPU) #use the specified number of processes
results = p.map(merge_sites, to_do_list, 1) #Send only one interval to each processor at a time.

elapsed_time = time.time()-overall_start_time
print('Done_calculating profiles '+str(int(np.floor(elapsed_time/60)))+' min '+str(int(np.round(elapsed_time%60)))+' sec')
del(elapsed_time)
sys.stdout.flush()


# In[ ]:


for key in results_dict_template.keys():
    current_results = pd.DataFrame()
    for i in range(len(results)):
        current_results = current_results.append(results[i][key]['coverage'])
    current_out_path = results_sample_dir+'/'+sample_name+'.'+key+'.coverage.tsv'
    current_results.to_csv(current_out_path,sep='\t', index = False, float_format='%.5f')


# In[ ]:


# # plt.plot(save_columns, results[0]['uncorrected']['coverage'][save_columns].mean())
# # plt.plot(save_columns, results[0]['GC_corrected']['coverage'][save_columns].mean())
# # plt.plot(save_columns, results[0]['GC_map_corrected']['coverage'][save_columns].mean())

# for i in range(len(results)):
#     df = results[i]['GC_map_corrected']['coverage']
#     site_name = df['site_name'].values[0]
#     plt.plot(save_columns, df[save_columns].mean(), label= site_name)
# plt.legend()
# plt.title(sample_name);


# In[ ]:


if not exclude_paths==['none']:
    os.remove(tmp_bigWig+"/excluded_regions.bw")


# In[ ]:





# In[ ]:





# In[ ]:




