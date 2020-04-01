# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 28 21:15:00 2018

@author: Jhnnyboy
"""

#------------------------------------------------------------------------------
#                             load/ convert data
#------------------------------------------------------------------------------
import os
import scipy.io as sio
session = 'X:/03. Lab Procedures and Protocols/Python/John/data'
os.chdir(session) # cd equivalent

# load pre-made dataset of subregions
sample_data = sio.loadmat('sample_data')

# erase irrelevant keys in the dictionary - could modify for loop later 
# on instead
del sample_data['__globals__']
del sample_data['__header__']
del sample_data['__version__']

# replace nans
import numpy as np
for key, value in sample_data.items():
    sample_data[key][np.isnan(sample_data[key])]=0
    sample_data[key] = np.ndarray.transpose(sample_data[key])
    
# calculate mean firing rate and create dicts for each subregion
acc = {"cho_tj": sample_data['ACC_choiceT'].mean(axis=1),
       "sam_tj": sample_data['ACC_sampleT'].mean(axis=1),
       "delay" : sample_data['ACC_delay'].mean(axis=1),
       "iti"   : sample_data['ACC_iti'].mean(axis=1)
       }

prl = {"cho_tj": sample_data['PrL_choiceT'].mean(axis=1),
       "sam_tj": sample_data['PrL_sampleT'].mean(axis=1),
       "delay" : sample_data['PrL_delay'].mean(axis=1),
       "iti"   : sample_data['PrL_iti'].mean(axis=1)
       }

vo = {"cho_tj": sample_data['VO_choiceT'].mean(axis=1),
      "sam_tj": sample_data['VO_sampleT'].mean(axis=1),
      "delay" : sample_data['VO_delay'].mean(axis=1),
      "iti"   : sample_data['VO_iti'].mean(axis=1)
      }

#------------------------------------------------------------------------------
#                       format data and dictionary
#------------------------------------------------------------------------------
   
import pandas as pd
data = np.array(pd.concat(
        [pd.DataFrame(
        (np.concatenate((prl['cho_tj'],acc['cho_tj'],vo['cho_tj']),axis=0))),
        pd.DataFrame(
        (np.concatenate((prl['sam_tj'],acc['sam_tj'],vo['sam_tj']),axis=0))),
        pd.DataFrame(
        (np.concatenate((prl['delay'],acc['delay'],vo['delay']),axis=0))), 
        pd.DataFrame(
        (np.concatenate((prl['iti'],acc['iti'],vo['iti']),axis=0)))],axis=1))

# create labels
twos = (np.full((np.size(vo['cho_tj']), 1), 2))
twos.resize((np.size(vo['cho_tj']),))
twos = twos.astype(float)

labels = np.concatenate((np.zeros(np.size(prl['cho_tj'])),
                   np.ones(np.size(acc['cho_tj'])),
                   twos))
labels.resize((np.size(labels),1))
labels = labels.astype(int)

# create feature names for easy reference
feature_names = ['choice t-junct', 
                 'sample t-junct',
                 'delay period', 
                 'intertrial-interval']

target_names  = np.asarray(['PrL','ACC', 'VO'])
 
subregion_dataset = {"data"         : data,
                     "feature_names": feature_names,
                     "target"       : labels,
                     "target_names" : target_names}

print("Target names: {}".format(subregion_dataset['target_names']))
print("Feature names: {}".format(subregion_dataset['feature_names']))

#------------------------------------------------------------------------------
#                                 Plot data
#------------------------------------------------------------------------------

# plot all data
sample_dataframe = pd.DataFrame(subregion_dataset['data'],
                                columns=subregion_dataset['feature_names'])
# convert to int32
labels.resize((np.size(labels),))

# plot data
colors=['green','blue','red']
axarr = pd.plotting.scatter_matrix(sample_dataframe, c=colors,
                           marker='o')

pd.plotting.scatter_matrix()

