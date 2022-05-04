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

from multiprocessing import Pool



# In[ ]:


parser = argparse.ArgumentParser()

parser.add_argument('--bam_file', help='sample_bam_file', required=True)
parser.add_argument('--bam_file_name', help='sample name (does not need to match actual file name)', required=True)
parser.add_argument('--mappable_regions_path', help='highly mappable regions to be used in GC correction, bedGraph or bed foramt', required=True)

parser.add_argument('--ref_seq',help='reference sequence (fasta format)',required=True)
parser.add_argument('--chrom_sizes',help='path to chromosome sizes for the reference seq',required=True)

parser.add_argument('--out_dir',help='folder for GC bias results',required=True)

parser.add_argument('--map_q',help='minimum mapping quality for reads to be considered',type=int,required=True)
parser.add_argument('--size_range',help='range of read sizes to be included',nargs=2, type=int, required=True)

parser.add_argument('--CPU',help='number of CPU for parallelizing', type=int, required=True)

args = parser.parse_args()

bam_file_path = args.bam_file
bam_file_name = args.bam_file_name
mappable_regions_path=args.mappable_regions_path

ref_seq_path = args.ref_seq
chrom_sizes_path = args.chrom_sizes
out_dir = args.out_dir

map_q = args.map_q
size_range = args.size_range
CPU = args.CPU


# In[ ]:


print('arguments provided:')

print('\tbam_file_path = "'+bam_file_path+'"')
print('\tbam_file_name = "'+bam_file_name+'"')
print('\tmappable_regions_path = "'+mappable_regions_path+'"')

print('\tref_seq_path = "'+ref_seq_path+'"')
print('\tchrom_sizes_path = "'+chrom_sizes_path+'"')
print('\tout_dir = "'+out_dir+'"')

print('\tmap_q = '+str(map_q))
print('\tsize_range = '+str(size_range))
print('\tCPU = '+str(CPU))


# In[ ]:


out_file = out_dir +'/GC_counts/'+ bam_file_name+'.GC_counts.txt'

print('out_file',out_file)

#create a directory for the GC data
if not os.path.exists(out_dir +'/GC_counts/'):
    os.mkdir(out_dir +'/GC_counts/')


# In[ ]:


#import filter
mappable_intervals = pd.read_csv(mappable_regions_path, sep='\t', header=None)

#remove non standard chromosomes and X and Y
chroms = ['chr'+str(m) for m in range(1,23)]
mappable_intervals = mappable_intervals[mappable_intervals[0].isin(chroms)]

print('chroms:', chroms)
print('number_of_intervals:',len(mappable_intervals))

sys.stdout.flush()


# In[ ]:


def collect_reads(sublist):
    #create a dict for holding the frequency of each read length and GC content
    GC_dict = {}
    for length in range(size_range[0],size_range[1]+1):
        GC_dict[length]={}
        for num_GC in range(0,length+1):
            GC_dict[length][num_GC]=0
        
    #import the bam file
    #this needs to be done within the loop otherwise it gives a truncated file warning
    bam_file = pysam.AlignmentFile(bam_file_path, "rb")
    print('sublist intervals:',len(sublist))
    
    #this might also need to be in the loop
    #import the ref_seq
    ref_seq=pysam.FastaFile(ref_seq_path)
    
    for i in range(len(sublist)):
        chrom = sublist.iloc[i][0]
        start = sublist.iloc[i][1]
        end = sublist.iloc[i][2]
        if i%5000==0:
            print('interval',i,':',chrom,start,end,'seconds:',np.round(time.time()-start_time))
            sys.stdout.flush()
        #fetch any read that overlaps the inteterval 
        fetched = bam_file.fetch(chrom,start,end)
        for read in fetched:
            #use both fw (positive template length) and rv (negative template length) reads
            if (read.is_reverse==False and read.template_length>=size_range[0] and read.template_length<=size_range[1]) or             (read.is_reverse==True and -read.template_length>=size_range[0] and -read.template_length<=size_range[1]):
                #qc filters, some longer fragments are considered 'improper pairs' but I would like to keep these
                if read.is_paired==True and read.mapping_quality>=map_q and read.is_duplicate==False and read.is_qcfail==False:
                    if read.is_reverse==False:
                        fragment_start = read.reference_start
                        fragment_end = read.reference_start+read.template_length
                    elif read.is_reverse==True:
                        fragment_end = read.reference_start + read.reference_length
                        fragment_start = fragment_end + read.template_length
    
                    #count the GC content
                    fragment_seq = ref_seq.fetch(read.reference_name,fragment_start,fragment_end)
                    fragment_seq = np.array(list(fragment_seq.upper()))
                    fragment_seq[np.isin(fragment_seq, ['A','T','W'])] = 0
                    fragment_seq[np.isin(fragment_seq, ['C','G','S'])] = 1
                    rng = np.random.default_rng(fragment_start)
                    fragment_seq[np.isin(fragment_seq, ['N','R','Y','K','M','B','D','H','V'])] = rng.integers(2, size=len(fragment_seq[np.isin(fragment_seq, ['N','R','Y','K','M','B','D','H','V'])])) #random integer in range(2) (i.e. 0 or 1)
                    fragment_seq = fragment_seq.astype(float)
                    
                    
                    num_GC = int(fragment_seq.sum())
                    GC_dict[abs(read.template_length)][num_GC]+=1

    print('done')
    return(GC_dict)


# In[ ]:


start_time = time.time()
p = Pool(processes=CPU) #use the available CPU
sublists = np.array_split(mappable_intervals,CPU) #split the list into sublists, one per CPU

GC_dict_list = p.map(collect_reads, sublists, 1)


# In[ ]:


all_GC_df = pd.DataFrame()
for i,GC_dict in enumerate(GC_dict_list):
    GC_df = pd.DataFrame()
    for length in GC_dict.keys():
        current = pd.Series(GC_dict[length]).reset_index()
        current = current.rename(columns={'index':'num_GC',0:'number_of_fragments'})
        current['length']=length
        current = current[['length','num_GC','number_of_fragments']]
        GC_df = GC_df.append(current, ignore_index=True)
    GC_df = GC_df.set_index(['length','num_GC'])
    all_GC_df[i] = GC_df['number_of_fragments']
    del(GC_df,GC_dict)
    
all_GC_df = all_GC_df.sum(axis=1)
all_GC_df = pd.DataFrame(all_GC_df).rename(columns = {0:'number_of_fragments'})
all_GC_df = all_GC_df.reset_index()
all_GC_df.to_csv(out_file,sep='\t',index=False)


# In[ ]:


print('done')


# In[ ]:





# In[ ]:




