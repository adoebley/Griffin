#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#import stuff 
import pandas as pd
import numpy as np
import math
import time
import pyBigWig
import sys
import argparse


# In[ ]:


# #TSS panel test
# # #define paths and parameters
# in_file='/fh/fast/ha_g/user/adoebley/data/SCLC_targeted_panel_sites/all_sites/TSS_targets.bed'
# in_file_name = 'TSS_targets'
# out_dir='./'
# mapability_file='../../downloads/genome/k50.Umap.MultiTrackMappability.hg38.bw'

# # #define the columns in the TFBS files
# chrom_col='Chrom'
# start_col='TSS'
# end_col='TSS'
# strand_col='Strand'
# chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22','chrX','chrY']

# window_values=(-100,255) #norm window
# targeted_window_columns=('window_start','window_end')
# targeted_panel = 'True'

# threshold = 0.95


# In[ ]:





# In[ ]:


parser = argparse.ArgumentParser()

parser.add_argument('-i','--in_file', help='file of TFBS to annotate', required=True)
parser.add_argument('--name', help='in file name', required=True)
parser.add_argument('-o','--out_dir', help='directory for output files', required=True)
parser.add_argument('-m','--mapability_file',help='.bw file with mapability values e.g. UCSC hg38 k50.Umap.MultiTrackMappability.bw', required=True)

parser.add_argument('-c','--chrom_column',help='name of column containing chromosome number', default='Chrom')
parser.add_argument('-s','--start_column',help='name of column containing the start of the TFBS (this will be averaged with the end to identify the TFBS center)', default='Start')
parser.add_argument('-e','--end_column',help='name of column containing the end of the TFBS (this will be averaged with the start to identify the TFBS center)', default='End')
parser.add_argument('--strand_column',help='name of column containing the strand (+ or -)', default='Strand')
parser.add_argument('--chroms', help='chromosomes to include when selecting sites', nargs='*')
parser.add_argument('--window_values',help='start and end of window to be analyzed around each TFBS',nargs=2, type=int, required=True)
parser.add_argument('--targeted_panel',help="whether the sites are from a targeted panel",default='False')
parser.add_argument('--targeted_window_columns',help='column names that specify the start and end of the window for a targeted region',nargs=2,default=('NA','NA'))

parser.add_argument('--threshold',help='define cutoff for high mapability', default=0.95, type=float)

args = parser.parse_args()

in_file=args.in_file
in_file_name = args.name
out_dir=args.out_dir.rstrip('/')+'/'
mapability_file=args.mapability_file

chrom_col=args.chrom_column
start_col=args.start_column
end_col=args.end_column
strand_col = args.strand_column
chroms = args.chroms

window_values=args.window_values
targeted_panel=args.targeted_panel
targeted_window_columns=args.targeted_window_columns

threshold=args.threshold


# In[ ]:


#set up parameters
print('\narguments provided:')
print('\tin_file = "'+in_file+'"')
print('\tin_file_name = "'+in_file_name+'"')
print('\tout_dir = "'+out_dir+'"')
print('\tmapability_file = "'+mapability_file+'"')

print('\tchrom_col = "'+chrom_col+'"')
print('\tstart_col = "'+start_col+'"')
print('\tend_col = "'+end_col+'"')
print('\tstrand_col = "'+strand_col+'"')
print('\tchroms = ',chroms)

print('\twindow_values=',window_values)
print('\ttargeted_panel = "'+targeted_panel+'"')

window_start_column=targeted_window_columns[0]
window_end_column=targeted_window_columns[1]
print('\ttargeted_window_columns = ',targeted_window_columns)

print('\tthreshold = '+str(threshold))

print('\n')


# In[ ]:


#import info about the sites
sites=pd.read_csv(in_file, sep='\t')
print('number_of_sites:',len(sites)) 
sys.stdout.flush()

#identify the TF position
sites['position'] = (sites[start_col]+sites[end_col])/2
sites['position'] = np.floor(sites['position'])
sites['position'] = sites['position'].astype(int)

#identify the window around each site to be used for mapability
sites['norm_window_start']=sites['position']+window_values[0]
sites['norm_window_end']=sites['position']+window_values[1]

#identify the window around reverse sites if direction is specified
if strand_col in sites.columns:
    print('flipping_reverse_sites')
    rv_sites = sites[sites[strand_col]=='-'].copy()
    rv_sites['norm_window_start']=rv_sites['position']-window_values[1]
    rv_sites['norm_window_end']=rv_sites['position']-window_values[0]
    sites[sites[strand_col]=='-'] = rv_sites.copy() #replace the rv sites with the flipped sites
    del(rv_sites)
    
     
#drop any sites that don't span the full window
if targeted_panel.lower()=='true': 
    
    sites['relative_window_start']=sites[window_start_column]-sites['position']
    sites['relative_window_end']=sites[window_end_column]-sites['position']
    
    if strand_col in sites.columns: #flip the orientation of the window for reverse sites
        print('flipping_reverse_sites')
        rv_sites = sites[sites[strand_col]=='-'].copy()
        rv_sites[['relative_window_start','relative_window_end']]=rv_sites[['relative_window_end','relative_window_start']]*-1
        sites[sites[strand_col]=='-'] = rv_sites #replace the rv sites with the flipped sites
        
    sites=sites[(sites['relative_window_start']<=window_values[0]) & (sites['relative_window_end']>=window_values[1])]
    #drop the new columns as they are no longer needed
    sites = sites.drop(columns=['relative_window_start','relative_window_end'])    

    print('sites that span the window:',len(sites))
    


# In[ ]:


if strand_col in sites.columns: 
    print('fw_sites:',len(sites[sites[strand_col]=='+']))
    print('rv_sites:',len(sites[sites[strand_col]=='-']))


# In[ ]:


sites = sites[sites[chrom_col].isin(chroms)]
print('sites_after_removing_non_specified_chroms:',len(sites))


# In[ ]:


mapability = pyBigWig.open(mapability_file)
chroms =  sites[chrom_col].unique()
for chrom in chroms:
    try:
        mapability.values(chrom,100000,100010)
    except:
        print('###\n###\nChromosome names dont match chromosomes in mapability file, check chromosome formatting\n###\n###\n')
        sys.exit(1)
        


# In[ ]:


#run analysis
start_time=time.time()

print(in_file_name)
sys.stdout.flush()

#import the mapability data
mapability = pyBigWig.open(mapability_file)

#make list to hold mean values
mean_values=[]

#get 1% intervals for tracking progress
one_percent=int(len(sites)/100)
if one_percent==0: #if there are less than 100 sites
    one_percent=1

window_len = window_values[1]-window_values[0]

for i in range(len(sites)):
    if i%one_percent==0:
        print(i,time.time()-start_time)
        sys.stdout.flush()
    #get the location of the site (chromosome and center of site)
    chrom = sites.iloc[i][chrom_col]
    
    position=sites.iloc[i]['position']
    start=sites.iloc[i]['norm_window_start']
    end=sites.iloc[i]['norm_window_end']
    
#     print(sites.iloc[i][strand_col],start-position)
    #fetch the values for that site
    try:
        fetched=mapability.values(chrom, start, end)
        #convert nan values in the data to zeros and convert fetched to np array
        fetched=np.nan_to_num(fetched) 
        
    except: #if the site can't be fetched
        fetched=[0 for m in range(window_len)]
        #convert the  data to an np array
        fetched=np.array(fetched)

    #reshape the data for adding to the array for all sites
    try:
        fetched=fetched.reshape(1,window_len)

    except: #if the full site wasn't fetched
        fetched=np.array([0 for m in range(window_len)])
        fetched=fetched.reshape(1,window_len)

    #add the fetched data to the array of data for all sites
    mean_values.append(fetched.mean())

    del (fetched,position,chrom)

#drop the start and end columns because these will be recalculated based on step in future analyses
sites = sites.drop(columns=['norm_window_start','norm_window_end'])    
print(time.time()-start_time)
sys.stdout.flush()


# In[ ]:


#calculate the mean value per site
sites['mean_mapability']=mean_values

overall_TF_mapability={'total_sites':len(sites)}

#split_list_of_sites for export and count list lengths

high_sites=sites[(sites['mean_mapability']>=threshold)]
low_sites=sites[(sites['mean_mapability']<threshold)]
high_sites.to_csv(out_dir+in_file_name+'.high_mapability'+'.txt', sep='\t', index=False)
low_sites.to_csv(out_dir+in_file_name+'.low_mapability'+'.txt', sep='\t', index=False)

overall_TF_mapability['high']=len(high_sites)
overall_TF_mapability['low']=len(low_sites)

overall_TF_mapability['mean_mapability']=sites['mean_mapability'].mean()

overall_TF_mapability=pd.DataFrame(overall_TF_mapability, index=[0])
print(overall_TF_mapability)


overall_TF_mapability.to_csv(out_dir+in_file_name+'.counts.txt', sep='\t', index=False)


# In[ ]:





# In[ ]:




