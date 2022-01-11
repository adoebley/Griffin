#!/usr/bin/env python
# coding: utf-8

# In[ ]:


def import_and_filter_sites(site_name,site_file,strand_column,chrom_column,position_column,chroms,ascending,sort_by,number_of_sites):
    import pandas as pd
    current_sites = pd.read_csv(site_file,sep='\t')
    if strand_column not in current_sites.columns:
        current_sites[strand_column]=0 
            
    #throw out sites that aren't on the selected chroms
    current_sites = current_sites[current_sites[chrom_column].isin(chroms)]
    
    #select the sites to use if specified
    if sort_by.lower()=='none': #if using all sites 
        print(site_name,'processing all '+str(len(current_sites))+' sites')
    
    else: #othewise sort by the specified column
        current_sites=current_sites.sort_values(by=sort_by,ascending=ascending).reset_index(drop=True)#sort and reset index
        current_sites=current_sites.iloc[0:int(number_of_sites)]
        print(site_name+'\tprocessing',len(current_sites),'sites\t('+
              str(sort_by),'range after sorting: ',min(current_sites[sort_by]),'to',
              str(max(current_sites[sort_by]))+')')

    current_sites = current_sites[[chrom_column,position_column,strand_column]]
    current_sites['site_name']=site_name
    return(current_sites)


# In[2]:


def define_fetch_interval(name_to_print,sites,chrom_column,position_column,chroms,chrom_sizes_path,upstream_bp,downstream_bp):
    import pandas as pd
    import numpy as np
    #separate fw and reverse sites
    fw_markers = ['+',1,'1']
    rv_markers = ['-',-1,'-1']
    fw_sites = sites[sites['Strand'].isin(fw_markers)].copy()
    rv_sites = sites[sites['Strand'].isin(rv_markers)].copy()

    undirected_sites = sites[~(sites['Strand'].isin(fw_markers+rv_markers))].copy()

    if len(rv_sites)+len(fw_sites)+len(undirected_sites)==len(sites):
        print(name_to_print+' (fw/rv/undirected/total): '+
              str(len(fw_sites))+'/'+
              str(len(rv_sites))+'/'+
              str(len(undirected_sites))+'/'+
              str(len(sites)))
    else: #I don't think this should ever happen...
        print('total fw sites:\t\t'+str(len(fw_sites)))
        print('total rv sites:\t\t'+str(len(rv_sites)))
        print('total undirected sites:'+'\t'+str(len(undirected_sites)))
        print('total sites:\t\t'+str(len(sites)))
        sys.exit('Problem with strand column')

    #set up to fetch a window extending across the desired window
    fw_sites['fetch_start'] = fw_sites[position_column]+upstream_bp
    fw_sites['fetch_end'] = fw_sites[position_column]+downstream_bp

    undirected_sites['fetch_start'] = undirected_sites[position_column]+upstream_bp
    undirected_sites['fetch_end'] = undirected_sites[position_column]+downstream_bp

    #for reverse sites, flip the window
    rv_sites['fetch_start'] = rv_sites[position_column]-downstream_bp
    rv_sites['fetch_end'] = rv_sites[position_column]-upstream_bp
    
    #merge fw and reverse back together and sort them back into the original order
    sites = fw_sites.append(rv_sites).append(undirected_sites).sort_index()
    sites = sites.sort_values(by = [chrom_column,position_column]).reset_index(drop=True)

    chrom_sizes = pd.read_csv(chrom_sizes_path, sep='\t', header=None)
    chrom_sizes = chrom_sizes[chrom_sizes[0].isin(chroms)]
    chrom_sizes = chrom_sizes.set_index(0)

    adjusted_ends_df = pd.DataFrame()

    for chrom in chroms:
        length = chrom_sizes.loc[chrom][1]
        current = sites[sites[chrom_column]==chrom].copy()
        current['fetch_start'] = np.where(current['fetch_start']<0,0,current['fetch_start'])
        current['fetch_end'] = np.where(current['fetch_end']>length,length,current['fetch_end'])    
        adjusted_ends_df = adjusted_ends_df.append(current)
    adjusted_ends_df = adjusted_ends_df.sort_values(by = [chrom_column,position_column]).reset_index(drop=True)
    adjusted_ends_df = adjusted_ends_df.copy()

    return(adjusted_ends_df)


# In[3]:


def progress_report(name,unit,start_time,current_time,item_index,total_items):
    import numpy as np
    time_elapsed = current_time-start_time
    elapsed_min = int(np.floor(time_elapsed/60))
    elapsed_sec = int(np.floor(time_elapsed%60))
    
    expected_time = (time_elapsed/(item_index+1))*total_items
    
    remaining_time = expected_time-time_elapsed
    remaining_min = int(np.floor(remaining_time/60))
    remaining_sec = int(np.floor(remaining_time%60))
    
    name = [str(m) for m in name]
    name_str = '_'.join(name)
    printout = name_str+': '+ str(item_index+1) +' of '+ str(total_items)+' '+str(unit)+' done in '+           str(elapsed_min)+' min '+str(elapsed_sec)+' sec, '+str(remaining_min)+' min '+str(remaining_sec)+' sec remaining'
    return(printout)


# In[ ]:




