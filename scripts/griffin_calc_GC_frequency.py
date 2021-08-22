#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pysam
import os
import pandas as pd
import numpy as np
import time
import argparse
import sys


# In[2]:


#This script calculates the frequency of each GC content for fragments that overlap the non-blacklisted areas
#This is performed for each fragment size in the range specified
#this only needs to be performed once for each filter


# In[3]:


# #arguments for testing 
# mapable_path = '/fh/fast/ha_g/user/adoebley/projects/griffin_paper/downloads/genome/repeat_masker.mapable.k50.Umap.hg38.bedGraph'

# ref_seq_path = '/fh/fast/ha_g/grp/reference/GRCh38/GRCh38.fa'
# chrom_sizes_path = '/fh/fast/ha_g/grp/reference/GRCh38/hg38.standard.chrom.sizes'
# out_dir = './tmp'

# # step = 1
# length = 50 #fragment length


# In[4]:


parser = argparse.ArgumentParser()

parser.add_argument('--mapable_regions', help='highly mapable regions to be used in GC correction, bed or bedGraph format', required=True)
parser.add_argument('--ref_seq',help='reference sequence (fasta format)',required=True)
parser.add_argument('--chrom_sizes',help='path to chromosome sizes for the reference seq',required=True)
parser.add_argument('--out_dir',help='folder for results',required=True)
parser.add_argument('--fragment_length',help='length of fragment (in bp) for which GC will be calculated',type=int, required=True)

args = parser.parse_args()

mapable_path=args.mapable_regions
ref_seq_path = args.ref_seq
chrom_sizes_path = args.chrom_sizes
out_dir = args.out_dir
length = args.fragment_length


# In[5]:


print('arguments provided:')

print('\tmapable_path = "'+mapable_path+'"')
print('\tref_seq_path = "'+ref_seq_path+'"')
print('\tchrom_sizes_path = "'+chrom_sizes_path+'"')
print('\tout_dir = "'+out_dir+'"')
print('\tlength = '+str(length))


# In[6]:


mapable_name = mapable_path.rsplit('/',1)[1].rsplit('.',1)[0]
out_file = out_dir+'/'+mapable_name+'.'+str(length)+'bp.GC_frequency.txt'
print('output path:',out_file)

if not os.path.exists(out_dir):
    os.mkdir(out_dir)


# In[7]:


sys.stdout.flush()


# In[8]:


#import filter
mapable_intervals = pd.read_csv(mapable_path, sep='\t', header=None)

#keep autosomes only
chroms = ['chr'+str(m) for m in range(1,23)]
mapable_intervals = mapable_intervals[mapable_intervals[0].isin(chroms)]

print('chroms:', chroms)
print('number_of_intervals:',len(mapable_intervals))
sys.stdout.flush()


# In[9]:


#get chrom sizes info
chrom_sizes = pd.read_csv(chrom_sizes_path, sep='\t', header=None)

#also keep as a dict
chrom_size_dict = chrom_sizes.set_index(0).to_dict()[1]


# In[10]:


#import the ref_seq
ref_seq=pysam.FastaFile(ref_seq_path)


# In[11]:


#create the GC frequencies dict
GC_dict = {}

GC_dict={}
for num_GC in range(0,length+1):
    GC_dict[num_GC]=0


# In[12]:


start_time = time.time()

k = length #just keeping this compatable with the previous version

for i in range(len(mapable_intervals)):
    chrom = mapable_intervals.iloc[i][0]
    start = mapable_intervals.iloc[i][1]
    end = mapable_intervals.iloc[i][2]
    if i%5000==0:
        print('interval',i,':',chrom,start,end,'seconds:',np.round(time.time()-start_time))
        sys.stdout.flush()
    #adjust the start and end so it includes all fragments that overlap the interval 
    adjusted_start = start-k
    adjusted_end = end+k
    
    if adjusted_start<0:
        adjusted_start = 0
    if adjusted_end>chrom_size_dict[chrom]:
        adjusted_end = chrom_sizes_dict[chrom]
        print(chrom,chrom_sizes_dict[chrom],'adjusting_end')

    fetched = ref_seq.fetch(chrom,adjusted_start,adjusted_end)
    fetched = fetched.replace('g','G').replace('c','C').replace('a','A').replace('t','T').replace('n','N')
    fetched = np.array(list(fetched.replace('G','1').replace('C','1').replace('A','0').replace('T','0').replace('N','2')),dtype=float)

    #swap the 2 for a random 1 or 0 #there has to be a better way to do this but I can't figure it out
    #the 0 or 1 is required because the sliding window sum algorithm only does integers
    #unknown nucleotides should be quite rare if the filter is done correctly
    fetched[fetched==2]=np.random.randint(2) #random integer in range(2) (i.e. 0 or 1)

    n=len(fetched)

    window_sum = int(sum(fetched[:k]))
    #print(k,len(fetched[:k]),window_sum)

    GC_dict[window_sum]+=1
    for i in range(n-k):
        window_sum = int(window_sum - fetched[i] + fetched[i+k])
        #print(k,window_sum)
        GC_dict[window_sum]+=1


# In[13]:


#convert to df and export
GC_df = pd.DataFrame()
#save GC dict
current = pd.Series(GC_dict).reset_index()
current = current.rename(columns={'index':'num_GC',0:'number_of_fragments'})
current['length']=length
current = current[['length','num_GC','number_of_fragments']]
GC_df = GC_df.append(current, ignore_index=True)
GC_df.to_csv(out_file,sep='\t',index=False)


# In[14]:


print('done')


# In[ ]:




