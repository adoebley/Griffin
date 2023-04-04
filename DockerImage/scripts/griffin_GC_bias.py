#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pysam
import os
#import pybedtools #not used
import pandas as pd
import numpy as np
import time
import argparse
import sys
from matplotlib import pyplot as plt


# In[ ]:


parser = argparse.ArgumentParser()

parser.add_argument('--bam_file_name', help='sample name (does not need to match actual file name)', required=True)

parser.add_argument('--mappable_name', help='name of mappable regions file (with .bed removed)', required=True)

parser.add_argument('--genome_GC_frequency',help='folder containing GC counts in the reference sequence (made by generate_reference_files.snakemake)',required=True)

parser.add_argument('--out_dir',help='folder for GC bias results',required=True)

parser.add_argument('--size_range',help='range of read sizes to be included',nargs=2, type=int, required=True)

args = parser.parse_args()

bam_file_name = args.bam_file_name
mappable_name=args.mappable_name
genome_GC_frequency = args.genome_GC_frequency
out_dir = args.out_dir
size_range = args.size_range


# In[3]:


print('arguments provided:')

print('\tbam_file_name = "'+bam_file_name+'"')
print('\tmappable_name = "'+mappable_name+'"')

print('\tgenome_GC_frequency = "'+genome_GC_frequency+'"')
out_dir = out_dir.rstrip('/')
print('\tout_dir = "'+out_dir+'"')

print('\tsize_range = '+str(size_range))


# In[4]:


#For now I'm going to keep the smoothing bin size as a set variable
GC_smoothing_step = 20


# In[5]:


#input is the out_file from the previous step
in_file = out_dir +'/GC_counts/'+ bam_file_name+'.GC_counts.txt'
print('in_file:',in_file)

#output is smoothed version
smoothed_out_file = out_dir +'/GC_bias/'+ bam_file_name+'.GC_bias.txt'

#plot files
plot_file1 = out_dir +'/GC_plots/'+ bam_file_name+'.GC_bias.pdf'
plot_file2 = out_dir +'/GC_plots/'+ bam_file_name+'.GC_bias.summary.pdf'
plot_file3 = out_dir +'/GC_plots/'+ bam_file_name+'.GC_bias.key_lengths.pdf'

print('out_file:',smoothed_out_file)
sys.stdout.flush()


# In[6]:


#create output folders if needed
if not os.path.exists(out_dir +'/GC_plots/'):
    os.mkdir(out_dir +'/GC_plots/')
if not os.path.exists(out_dir +'/GC_bias/'):
    os.mkdir(out_dir +'/GC_bias/')


# In[7]:


#import the GC info from the genome
frequency_prefix = genome_GC_frequency+'/'+mappable_name+'.'

GC_freq = pd.DataFrame()
for i in range(size_range[0],size_range[1]+1):
    current_path = frequency_prefix+str(i)+'bp.GC_frequency.txt'
    current_data = pd.read_csv(current_path,sep='\t')
    GC_freq = GC_freq.append(current_data, ignore_index=True)
    
GC_freq['GC_content']=GC_freq['num_GC']/GC_freq['length']
GC_freq = GC_freq.sort_values(by=['GC_content','length']).reset_index(drop=True)


# In[8]:


#import GC counts from the sample
GC_df = pd.read_csv(in_file, sep='\t')

GC_df['GC_content']=GC_df['num_GC']/GC_df['length']
GC_df = GC_df.sort_values(by=['GC_content','length']).reset_index(drop=True)


# In[9]:


#calculate the GC_bias
new_df = pd.DataFrame()
for length in range(size_range[0],size_range[1]+1):
    current = GC_df[GC_df['length']==length].copy().reset_index(drop=True)
    current_freq = GC_freq[GC_freq['length']==length].copy().reset_index(drop=True)
    
    #save the frequency of each GC content in the genome
    current['number_of_positions']=current_freq['number_of_fragments']
    
    #calculate the GC bias
    current_bias = current['number_of_fragments']/current['number_of_positions']    
    current['GC_bias'] = current_bias

    #normalize to a mean of 1 for each fragment length(compute GC bias does this same thing)
    current['GC_bias'] = current['GC_bias']/np.nanmean(current['GC_bias'])
    new_df = new_df.append(current, ignore_index=True)
    
    #print(length,len(current['GC_bias']),np.nanmean(current['GC_bias']))
    
new_df = new_df.sort_values(by=['GC_content','length']).reset_index(drop=True)


# In[10]:


def median_smoothing(current,fraction):
    bin_size=int(len(current)*fraction)
    if bin_size<50:
        bin_size=50
    medians = []

    for i in range(len(current)):
        start = int(i-bin_size/2)
        end = int(i+bin_size/2)
        #if the bin starts before the beginning, just take the first bin
        if start<0:
            start=0
            end=bin_size
        #if the bin extends beyond the end, take the last bin
        if end>=len(current):
            start=len(current)-bin_size
            end=len(current)
        current_median = np.nanmedian(current['GC_bias'].iloc[start:end])
        medians.append(current_median)
    return(medians)


# In[ ]:





# In[11]:


#smooth GC bias by size bin

start_time = time.time()

new_df2 = pd.DataFrame()
for length in new_df['length'].unique():
    if length%20==0:
        print(length, time.time()-start_time)
        sys.stdout.flush()
        
    #get a bin of similar sized fragments
    min_len = int(length - (GC_smoothing_step/2))
    max_len = int(length + (GC_smoothing_step/2))
    
    current = new_df[(new_df['length']>=min_len) & (new_df['length']<=max_len)].copy()

    #perform smoothing
    fit = median_smoothing(current,.05)  
    current['smoothed_GC_bias']=fit
    
    #only keep smoothed values for the selected length
    current = current[current['length']==length]
    
    #get rid of values for GC contents that are never observed
    current['smoothed_GC_bias'] = np.where(current['number_of_positions']==0,np.nan,current['smoothed_GC_bias'])
    
    #normalize to a mean of 1
    current['smoothed_GC_bias'] = current['smoothed_GC_bias']/np.nanmean(current['smoothed_GC_bias'])
    
    new_df2 = new_df2.append(current,ignore_index=True)
    
    #print(length,len(current),np.nanmean(current['smoothed_GC_bias']))
    
new_df = new_df2


# In[12]:


#export results
new_df2.to_csv(smoothed_out_file,sep='\t',index=False)


# In[13]:


#generate one plot per size bin

#set up a figure for plotting
plot_indexes = np.arange(size_range[0]+GC_smoothing_step,size_range[1]+GC_smoothing_step,GC_smoothing_step)
lengths_to_plot = plot_indexes
x_dim = 6
y_dim = int(np.ceil(len(plot_indexes)/6))
empty_plots = int(x_dim*y_dim - len(plot_indexes))
plot_indexes = np.append(plot_indexes,[np.nan for m in range(empty_plots)])
plot_indexes = np.reshape(plot_indexes,(y_dim,x_dim))
fig, axes = plt.subplots(y_dim,x_dim, figsize = (5*x_dim,3.5*y_dim), sharex = True, sharey = True)
axes = axes.reshape(y_dim,x_dim) #make sure the axes array is two dimension (just in case it has less than 7 value)

#do the plotting
min_len = 0 
for max_len in lengths_to_plot:
    if max_len%20==0:
        print(max_len)
        
    #pick the axis
    current_index = np.where(plot_indexes==max_len)
    current_index = (current_index[0][0],current_index[1][0])
    current_ax = axes[current_index]

    #pick the data
    current1 = new_df2[(new_df2['length']>min_len) & (new_df2['length']<=max_len)].copy()
    
    #plot the smoothed data over top
    for length2 in current1['length'].unique():
        current2 = current1[current1['length']==length2]       
        current_ax.plot(current2['GC_content'],current2['smoothed_GC_bias'], label=str(length2)+'bp')
    
    current_ax.set_title(str(min_len) + 'bp to '+str(max_len)+'bp')
    current_ax.legend(ncol = 2)
    
    min_len = max_len
    
for i in range(x_dim):
    axes[y_dim-1,i].set_xlabel('GC content')
    
for i in range(y_dim):
    axes[i,0].set_ylabel('coverage bias')

ylim = axes[0,0].get_ylim()

old_title = axes[0,0].get_title()
axes[0,0].set_title(bam_file_name+'\n'+mappable_name + '\n' + old_title)

fig.tight_layout()

plt.savefig(plot_file1)

plt.close('all')


# In[14]:


#key lengths
selected_lengths = np.arange(100,201,GC_smoothing_step)

fig,ax = plt.subplots(1)

# for_color = len(selected_lengths)-1
# color = (1-(i/for_color),.5*(1-(i/for_color)), i/for_color)

for i,length in enumerate(selected_lengths):
    current = new_df2[new_df2['length']==length]
    ax.plot(current['GC_content'],current['smoothed_GC_bias'], label = str(length)+'bp')
    
ax.legend(ncol = 2, bbox_to_anchor = [1,1], loc = 'upper left')

ax.set_xlabel('GC content')
ax.set_ylabel('coverage bias')
ax.set_title(bam_file_name+'\n'+mappable_name)

fig.tight_layout()
fig.savefig(plot_file3)
plt.close('all')


# In[15]:


#summary figure
selected_lengths = np.arange(size_range[0],size_range[1],GC_smoothing_step)

fig,ax = plt.subplots(1)

for length in selected_lengths:
    current = new_df2[new_df2['length']==length]
    ax.plot(current['GC_content'],current['smoothed_GC_bias'], label = str(length)+'bp')
ax.legend(ncol = 2, bbox_to_anchor = [1,1], loc = 'upper left')

ax.set_xlabel('GC content')
ax.set_ylabel('coverage bias')
ax.set_title(bam_file_name+'\n'+mappable_name)

fig.tight_layout()
fig.savefig(plot_file2)
plt.close('all')


# In[ ]:





# In[ ]:


# plot_file4 = out_dir +'/'+mappable_name+'/GC_plots/'+ bam_file_name+'.GC_bias.test.pdf'

# selected_lengths = np.arange(size_range[0],size_range[1],GC_smoothing_step)

# fig,ax = plt.subplots(1)

# for length in selected_lengths:
#     current = new_df2[new_df2['length']==length]
    
#     ax.plot(current['GC_content'],current['GC_bias'],alpha=.2,marker='.')
        
# #reset the color cycle
# # for Matplotlib version >= 1.5
# plt.gca().set_prop_cycle(None)
    
    
# for length in selected_lengths:
#     current = new_df2[new_df2['length']==length]
    
#     ax.plot(current['GC_content'],current['smoothed_GC_bias'], label = length)
    

# ax.legend(ncol = 2, bbox_to_anchor = [1,1])

# ax.set_xlabel('GC content')
# ax.set_ylabel('coverage bias')
# ax.set_title(bam_file_name+'\n'+mappable_name)
# ax.set_ylim(-.1,new_df2['smoothed_GC_bias'].max()+.1)

# fig.tight_layout()
# fig.savefig(plot_file4)


# In[ ]:


# #generate one plot per size bin
# #raw_data
# plot_file4 = out_dir +'/'+mappable_name+'/GC_plots/'+ bam_file_name+'.GC_bias.test.pdf'

# #set up a figure for plotting
# plot_indexes = np.arange(size_range[0]+GC_smoothing_step,size_range[1]+GC_smoothing_step,GC_smoothing_step)
# lengths_to_plot = plot_indexes
# x_dim = 6
# y_dim = int(np.ceil(len(plot_indexes)/6))
# empty_plots = int(x_dim*y_dim - len(plot_indexes))
# plot_indexes = np.append(plot_indexes,[np.nan for m in range(empty_plots)])
# plot_indexes = np.reshape(plot_indexes,(y_dim,x_dim))
# fig, axes = plt.subplots(y_dim,x_dim, figsize = (5*x_dim,3.5*y_dim), sharex = True, sharey = True)
# axes = axes.reshape(y_dim,x_dim) #make sure the axes array is two dimension (just in case it has less than 7 value)

# #do the plotting
# min_len = 0 
# for max_len in lengths_to_plot:
#     if max_len%20==0:
#         print(max_len)
        
#     #pick the axis
#     current_index = np.where(plot_indexes==max_len)
#     current_index = (current_index[0][0],current_index[1][0])
#     current_ax = axes[current_index]

#     #pick the data
#     current1 = new_df2[(new_df2['length']>min_len) & (new_df2['length']<=max_len)].copy()
    
#     #plot the raw data
#     for length2 in current1['length'].unique():
#         current2 = current1[current1['length']==length2]
#         current_ax.plot(current2['GC_content'],current2['GC_bias'],alpha=.2,marker='.')
        
#     #reset the color cycle
#     # for Matplotlib version >= 1.5
#     plt.gca().set_prop_cycle(None)
    
#     #plot the smoothed data over top
#     for length2 in current1['length'].unique():
#         current2 = current1[current1['length']==length2]       
#         current_ax.plot(current2['GC_content'],current2['smoothed_GC_bias'], label=length2)
    
#     current_ax.set_title(str(min_len) + 'bp to '+str(max_len)+'bp')
#     current_ax.legend(ncol = 2)
    
#     min_len = max_len
    
# for i in range(x_dim):
#     axes[y_dim-1,i].set_xlabel('GC content')
    
# for i in range(y_dim):
#     axes[i,0].set_ylabel('coverage bias')

# axes[0,0].set_ylim(ylim)

# old_title = axes[0,0].get_title()
# axes[0,0].set_title(bam_file_name+'\n'+mappable_name + '\n' + old_title)

# fig.tight_layout()

# plt.savefig(plot_file4)
# plt.close('all')


# In[ ]:




