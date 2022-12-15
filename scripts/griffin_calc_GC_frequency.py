#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pysam
import os
import pandas as pd
import numpy as np
import time
import argparse
import sys
import pybedtools


# In[ ]:


parser = argparse.ArgumentParser()

parser.add_argument('--mappable_regions_path', help='highly mappable regions to be used in GC correction, bed or bedGraph format', required=True)
parser.add_argument('--ref_seq',help='reference sequence (fasta format)',required=True)
parser.add_argument('--chrom_sizes',help='path to chromosome sizes for the reference seq',required=True)
parser.add_argument('--out_dir',help='folder for results',required=True)
parser.add_argument('--fragment_length',help='length of fragment (in bp) for which GC will be calculated',type=int, required=True)
parser.add_argument('--read_length',help='length of read (in bp)',type=int, required=True)

args = parser.parse_args()

mappable_regions_path=args.mappable_regions_path
ref_seq_path = args.ref_seq
chrom_sizes_path = args.chrom_sizes
out_dir = args.out_dir
fragment_length = args.fragment_length
read_length = args.read_length


# In[ ]:


mappable_name = mappable_regions_path.rsplit('/',1)[1].rsplit('.',1)[0]
out_file = out_dir+'/'+mappable_name+'.'+str(fragment_length)+'bp.GC_frequency.txt'


# In[ ]:


#keep autosomes only
chroms = ['chr'+str(m) for m in range(1,23)]


# In[ ]:


if read_length>fragment_length:
    read_length = fragment_length 


# In[ ]:


print('mappable_regions_path:',mappable_regions_path)
print('out_file:',out_file)
print('fragment_length:',fragment_length)
print('read_length:',read_length)


# In[ ]:


#import filter
mappable_intervals = pd.read_csv(mappable_regions_path, sep='\t', header=None)

mappable_intervals = mappable_intervals[mappable_intervals[0].isin(chroms)]

print('chroms:', chroms)
print('number_of_intervals:',len(mappable_intervals))
sys.stdout.flush()


# In[ ]:


#get chrom sizes info
chrom_sizes = pd.read_csv(chrom_sizes_path, sep='\t', header=None)

#also keep as a dict
chrom_size_dict = chrom_sizes.set_index(0).to_dict()[1]


# In[ ]:


#import the ref_seq
ref_seq=pysam.FastaFile(ref_seq_path)


# In[ ]:


#count the GC content of all fragments where the forward read overlaps the specified regions
start_time = time.time()
print('Calculating forward read frequency')

#create the GC frequencies dict
fw_GC_dict = {}
for num_GC in range(0,fragment_length+1):
    fw_GC_dict[num_GC]=0
    
for i in range(len(mappable_intervals)):
    chrom = mappable_intervals.iloc[i][0]
    start = mappable_intervals.iloc[i][1]+1
    end = mappable_intervals.iloc[i][2]-1
    if i%5000==0:
        print('interval',i,':',chrom,start,end,'seconds:',np.round(time.time()-start_time))
        sys.stdout.flush()
    
    #count up all possible fw reads that overlap the interval
    #adjust the start and end so it includes all possible fragment that overlap the interval 
    adjusted_start = start-read_length
    adjusted_end = end+fragment_length
    
    if adjusted_start<0:
        adjusted_start = 0
    if adjusted_end>chrom_size_dict[chrom]:
        adjusted_end = chrom_size_dict[chrom]
        print(chrom,chrom_size_dict[chrom],'modifying_end_to_end_of_chromosome')
    
    #count the GC content
    fetched = ref_seq.fetch(chrom,adjusted_start,adjusted_end)
    fetched = np.array(list(fetched.upper()))
    fetched[np.isin(fetched, ['A','T','W'])] = 0
    fetched[np.isin(fetched, ['C','G','S'])] = 1
    rng = np.random.default_rng(start)
    fetched[np.isin(fetched, ['N','R','Y','K','M','B','D','H','V'])] = rng.integers(2, size=len(fetched[np.isin(fetched, ['N','R','Y','K','M','B','D','H','V'])])) #random integer in range(2) (i.e. 0 or 1)
    fetched = fetched.astype(float)
    
    n=len(fetched)
    
    window_sum = int(sum(fetched[:fragment_length]))
    #print(k,len(fetched[:k]),window_sum)

    fw_GC_dict[window_sum]+=1
    for i in range(n-fragment_length):
        window_sum = int(window_sum - fetched[i] + fetched[i+fragment_length])
        fw_GC_dict[window_sum]+=1


# In[ ]:


#count the GC content of all fragments where the reverse read overlaps the specified regions
print('Calculating reverse read frequency')

#create the GC frequencies dict
rv_GC_dict = {}
for num_GC in range(0,fragment_length+1):
    rv_GC_dict[num_GC]=0
    
for i in range(len(mappable_intervals)):
    chrom = mappable_intervals.iloc[i][0]
    start = mappable_intervals.iloc[i][1]+1 #skip the first and last positions because these reads aren't fetched by pysam
    end = mappable_intervals.iloc[i][2]-1
    if i%5000==0:
        print('interval',i,':',chrom,start,end,'seconds:',np.round(time.time()-start_time))
        sys.stdout.flush()
    
    #count up all possible rv reads that overlap the interval
    #adjust the start and end so it includes all possible fragment that overlap the interval 
    adjusted_start = start-fragment_length
    adjusted_end = end+read_length
    
    if adjusted_start<0:
        adjusted_start = 0
    if adjusted_end>chrom_size_dict[chrom]:
        adjusted_end = chrom_size_dict[chrom]
        print(chrom,chrom_size_dict[chrom],'modifying_end_to_end_of_chromosome')

    #count the GC content
    fetched = ref_seq.fetch(chrom,adjusted_start,adjusted_end)
    fetched = np.array(list(fetched.upper()))
    fetched[np.isin(fetched, ['A','T','W'])] = 0
    fetched[np.isin(fetched, ['C','G','S'])] = 1
    rng = np.random.default_rng(start)
    fetched[np.isin(fetched, ['N','R','Y','K','M','B','D','H','V'])] = rng.integers(2, size=len(fetched[np.isin(fetched, ['N','R','Y','K','M','B','D','H','V'])])) #random integer in range(2) (i.e. 0 or 1)
    fetched = fetched.astype(float)
    
    n=len(fetched)

    window_sum = int(sum(fetched[:fragment_length]))
    #print(k,len(fetched[:k]),window_sum)

    rv_GC_dict[window_sum]+=1
    for i in range(n-fragment_length):
        window_sum = int(window_sum - fetched[i] + fetched[i+fragment_length])
        rv_GC_dict[window_sum]+=1


# In[ ]:


#convert to df and export
GC_df = pd.DataFrame()
#save GC dict
current = (pd.Series(rv_GC_dict)+pd.Series(fw_GC_dict)).reset_index()
current = current.rename(columns={'index':'num_GC',0:'number_of_fragments'})
current['length']=fragment_length
current = current[['length','num_GC','number_of_fragments']]
GC_df = GC_df.append(current, ignore_index=True)
GC_df.to_csv(out_file,sep='\t',index=False)


# In[ ]:


print('done')


