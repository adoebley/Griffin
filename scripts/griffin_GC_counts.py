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


# ##arguments for testing 

# bam_file_path = '/fh/scratch/delete90/ha_g/realigned_bams/cfDNA_MBC_ULP_hg38/realign_bam_paired_snakemake-master/results/MBC_1041_1_ULP/MBC_1041_1_ULP_recalibrated.bam'
# bam_file_name = 'MBC_1041_1_ULP'
# mapable_path = '../../downloads/genome/repeat_masker.mapable.k50.Umap.hg38.bedGraph'

# ref_seq_path = '/fh/fast/ha_g/grp/reference/GRCh38/GRCh38.fa'
# chrom_sizes_path = '/fh/fast/ha_g/grp/reference/GRCh38/hg38.standard.chrom.sizes'

# out_dir = './tmp/'

# map_q = 20
# size_range = [15,500]

# CPU = 4


# In[ ]:


parser = argparse.ArgumentParser()

parser.add_argument('--bam_file', help='sample_bam_file', required=True)
parser.add_argument('--bam_file_name', help='sample name (does not need to match actual file name)', required=True)
parser.add_argument('--mapable_regions', help='highly mapable regions to be used in GC correction, bedGraph or bed foramt', required=True)

parser.add_argument('--ref_seq',help='reference sequence (fasta format)',required=True)
parser.add_argument('--chrom_sizes',help='path to chromosome sizes for the reference seq',required=True)

parser.add_argument('--out_dir',help='folder for GC bias results',required=True)

parser.add_argument('--map_q',help='minimum mapping quality for reads to be considered',type=int,required=True)
parser.add_argument('--size_range',help='range of read sizes to be included',nargs=2, type=int, required=True)

parser.add_argument('--CPU',help='number of CPU for parallelizing', type=int, required=True)

args = parser.parse_args()

bam_file_path = args.bam_file
bam_file_name = args.bam_file_name
mapable_path=args.mapable_regions

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
print('\tmapable_regions = "'+mapable_path+'"')

print('\tref_seq_path = "'+ref_seq_path+'"')
print('\tchrom_sizes_path = "'+chrom_sizes_path+'"')
print('\tout_dir = "'+out_dir+'"')

print('\tmap_q = '+str(map_q))
print('\tsize_range = '+str(size_range))
print('\tCPU = '+str(CPU))


# In[ ]:


mapable_name = mapable_path.rsplit('/',1)[1].rsplit('.',1)[0]
out_file = out_dir +'/'+mapable_name+'/GC_counts/'+ bam_file_name+'.GC_counts.txt'

print('out_file',out_file)


# In[ ]:


#create a directory for the GC data
if not os.path.exists(out_dir +'/'+mapable_name):
    os.mkdir(out_dir +'/'+mapable_name)
if not os.path.exists(out_dir +'/'+mapable_name+'/GC_counts/'):
    os.mkdir(out_dir +'/'+mapable_name+'/GC_counts/')


# In[ ]:


#import filter
mapable_intervals = pd.read_csv(mapable_path, sep='\t', header=None)

#remove non standard chromosomes and X and Y
chroms = ['chr'+str(m) for m in range(1,23)]
mapable_intervals = mapable_intervals[mapable_intervals[0].isin(chroms)]

print('chroms:', chroms)
print('number_of_intervals:',len(mapable_intervals))

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
        #fetch any read that overlaps the inteterval (don't need to extend the interval because the fetch function does this automatically)
        fetched = bam_file.fetch(chrom,start,end)
        for read in fetched:
            #use both fw (positive template length) and rv (negative template length) reads
            if (read.is_reverse==False and read.template_length>=size_range[0] and read.template_length<=size_range[1]) or             (read.is_reverse==True and -read.template_length>=size_range[0] and -read.template_length<=size_range[1]):
                #qc filters, some longer fragments are considered 'improper pairs' but I would like to keep these
                if read.is_paired==True and read.mapping_quality>=map_q and read.is_duplicate==False and read.is_qcfail==False:
                    if read.is_reverse==False:
                        read_start = read.reference_start
                        read_end = read.reference_start+read.template_length
                    elif read.is_reverse==True:
                        read_end = read.reference_start + read.reference_length
                        read_start = read_end + read.template_length

                    fragment_seq = ref_seq.fetch(read.reference_name,read_start,read_end)
                    #tally up the GC content
                    fragment_seq=fragment_seq.replace('g','G').replace('c','C').replace('a','A').replace('t','T').replace('n','N')

    #                 #################
    #                 ##logic check####
    #                 #################
    #                 if read.is_reverse==False:
    #                     if fragment_seq[0:read.reference_length]==read.query_sequence and len(fragment_seq)==read.template_length:
    #                         print('fw match',read.reference_length)
    #                     else:
    #                         print(fragment_seq[0:read.reference_length],read.reference_length,'fw')
    #                         print(read.query_sequence,len(read.query_sequence),'fw')
    #                         print(len(fragment_seq),read.template_length)
    #                         print('\n')
    #                 elif read.is_reverse==True:
    #                     if fragment_seq[-read.reference_length:]==read.query_sequence and len(fragment_seq)==-read.template_length:
    #                         print('rv match',read.reference_length)
    #                     else:
    #                         print(fragment_seq[-read.reference_length:],read.reference_length,'rv')
    #                         print(read.query_sequence,len(read.query_sequence),'rv')
    #                         print(len(fragment_seq),read.template_length)
    #                         print('\n')                        
    #                 #################

                    #split and convert to numpy array
                    fragment_seq = np.array(list(fragment_seq))
                    #replace with values
                    fragment_seq[(fragment_seq=='G') | (fragment_seq=='C')]=1
                    fragment_seq[(fragment_seq=='A') | (fragment_seq=='T')]=0
                    fragment_seq[(fragment_seq=='N')]=np.random.randint(2) #choose a random 0 or 1 for N (so that you always get an integer) #should be very rare if the filter is done right
                    fragment_seq = fragment_seq.astype(int)

                    num_GC = int(fragment_seq.sum())
                    GC_dict[abs(read.template_length)][num_GC]+=1

    print('done')
    return(GC_dict)


# In[ ]:


start_time = time.time()
p = Pool(processes=CPU) #use the available CPU
sublists = np.array_split(mapable_intervals,CPU) #split the list into sublists, one per CPU

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





# In[ ]:




