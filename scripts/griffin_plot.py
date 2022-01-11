#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import time
import sys
import argparse


# In[ ]:


# #for ipynb
# %matplotlib inline

# in_dir = '../../../demo/griffin_nucleosome_profiling/results/coverage/all_sites/'

# in_files = []
# for file in os.listdir(in_dir):
#     in_files.append(in_dir+'/'+file)
    
# del(in_dir)

# #actual arguments
# in_files = in_files

# plot_window = [-500,500]
# step = 15

# individual = 'True'
# out_dir = './tmp'


# In[ ]:


parser = argparse.ArgumentParser()

parser.add_argument('--in_files', help='coverge files for all samples', nargs = '*', required=True)

parser.add_argument('--plot_window',help='start and end of window to be plotted',nargs=2, type=int, default=(-1000,1000))
parser.add_argument('--step',help='step size when calculating coverage', type=int, default=5)

parser.add_argument('--individual',help='if individual sites were saved in previous steps. (True/False)',default='False')
parser.add_argument('--out_dir',help='folder for results (new plots folder will be generated inside)',required=True)

args = parser.parse_args()

in_files = args.in_files

plot_window=args.plot_window
step = args.step

individual = args.individual
out_dir = args.out_dir


# In[ ]:


print('\nArguments provided:')
print('\tin_files = ',in_files)

print('\tplot_window = '+str(plot_window))
plot_window=[int(np.ceil(plot_window[0]/step)*step),int(np.floor(plot_window[1]/step)*step)] #round to the nearest step inside the window
print('\t#plot_window rounded to step:',plot_window)

print('\tstep =',step)
print('\tindividual = "'+individual+'"')
print('\tout_dir = "'+out_dir+'"')


# In[ ]:


#set up global variables
plot_columns = np.arange(plot_window[0],plot_window[1],step)
str_plot_columns = [str(m) for m in plot_columns]

print(plot_columns)


# In[ ]:


start_time = time.time()
data = pd.DataFrame()
for file in in_files:
    new_file = pd.read_csv(file,sep='\t')
    if individual.lower()=='true':
        print('taking means',file)
        for site_name in new_file['site_name'].unique():
            print(site_name, time.time()-start_time)
            sys.stdout.flush()
            for normalization in ['none','GC_corrected']:
                current = new_file[(new_file['site_name']==site_name) & (new_file['GC_correction']==normalization)]
                current = current[str_plot_columns].mean()
                current['site_name']=site_name
                current['GC_correction']=normalization
                current['sample']=new_file['sample'].iloc[0]
                data = data.append(current, ignore_index=True)
    else:
        current=new_file
        data = data.append(current, ignore_index=True)


# In[ ]:


#generate plots
for site_name in data['site_name'].unique():
    fig,axes = plt.subplots(1,2,figsize=(10,3.5), sharey = 'row')
    for i,normalization in enumerate(['none','GC_corrected']):
        ax = axes[i]
        for sample in data['sample'].unique():
            current = data[(data['sample']==sample) & (data['site_name']==site_name) & (data['GC_correction']==normalization)]
            ax.plot(plot_columns, current[str_plot_columns].T, label=sample)
            ax.tick_params(labelleft=True)
        ax.set_title(site_name+' '+normalization)

    axes[0].set_ylabel('normalized coverage')
    axes[0].set_xlabel('distance from site')
    axes[1].set_xlabel('distance from site')

    if len(data['sample'].unique())<15:
        axes[1].legend(bbox_to_anchor=[1,1],loc = 'upper left')
    else:
        axes[1].legend(bbox_to_anchor=[1,1],loc = 'upper left',ncol=2)

    fig.tight_layout()
    plt.savefig(out_dir+'/plots/'+site_name+'.pdf')
    plt.close('all')


# In[ ]:


#if info about individual sites was kept, the averaging process can take quite a while. Save for later use. 
if individual.lower()=='true':
    data.to_csv(out_dir+'/plots/'+'mean_data_for_all_sites.txt', sep='\t', index=False)
    


# In[ ]:





# In[ ]:




