# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 13:31:58 2019

@author: jstout
"""

#------------------------------------------------------------------------------
#                      convert and load variables
#------------------------------------------------------------------------------
import scipy.io as sio
import os
sessions_folder = 'X:/01.Experiments/John n Andrew/Dual optogenetics w mPFC recordings/All Subjects - DNMP/Good performance/Medial Prefrontal Cortex'
session = '/Baby Groot 9-16-18'
os.chdir(sessions_folder+session) # cd equivalent

# load Int and TimeVariables dictionaries
IntVars       = sio.loadmat('Int_file.mat')
TimeVariables = sio.loadmat('VT1.mat')

# If you want to isolate an element from the dictionary
Int        = IntVars["Int"]
ExtractedY = TimeVariables["ExtractedY"]
TimeStamps = TimeVariables["TimeStamps_VT"]

# index sample and choice
sam = Int[::2]
cho = Int[1::2]

# Load .txt files into one variable like in matlab
import glob
path     = sessions_folder+session+'/TT*.txt'
clusters = glob.glob(path)

#------------------------------------------------------------------------------
#                          calculate firing rate
#------------------------------------------------------------------------------
import numpy as np
import pandas as pd

# designate locations of maze
mazei = [0, 5]
    
# initialize some variables
mean_sample  = np.array([])
mean_choice  = np.array([])
mean_delay   = np.array([])
sample_array = pd.DataFrame(np.array([]))
choice_array = pd.DataFrame(np.array([]))
delay_array  = pd.DataFrame(np.array([]))

# loop across clusters
for i in range(0,np.size(clusters)):
    spkTimes = ()
    spkTimes = np.loadtxt(clusters[i])
    
    # loop over trials
    fr_sample = np.array([])
    fr_choice = np.array([])
    fr_delay  = np.array([])
    
    for ii in range(0,np.size(sam,0)):
        
        # firing rate during stem of sample trials
        fr_s = (np.size(np.where(np.logical_and(spkTimes>=sam[ii,mazei[0]], 
             spkTimes<=sam[ii,mazei[1]]))))/((sam[ii,mazei[1]]-
             sam[ii,mazei[0]])/1e6)  
 
        # firing rate during stem of choice trials
        fr_c = (np.size(np.where(np.logical_and(spkTimes>=cho[ii,mazei[0]], 
             spkTimes<=cho[ii,mazei[1]]))))/((cho[ii,mazei[1]]-
             cho[ii,mazei[0]])/1e6)
 
        # firing rate during startbox of delay
        fr_d = (np.size(np.where(np.logical_and(spkTimes>=sam[ii,7], 
             spkTimes<=cho[ii,0]))))/((cho[ii,0]-sam[ii,7])/1e6)     

        # store trial data
        fr_sample = np.append(fr_sample, fr_s)
        fr_choice = np.append(fr_choice, fr_c)
        fr_delay  = np.append(fr_delay, fr_d)
        
        # calculate mean
        mean_temp1 = np.mean(fr_sample)
        mean_temp2 = np.mean(fr_choice)
        mean_temp3 = np.mean(fr_delay)
    
    # store all cells mean firing rate data
    mean_sample = np.append(mean_sample,mean_temp1)
    mean_choice = np.append(mean_choice,mean_temp2)
    mean_delay  = np.append(mean_delay,mean_temp3)
    
    # store trial firing rate
    fr_sample    = pd.DataFrame(fr_sample)
    fr_choice    = pd.DataFrame(fr_choice)
    fr_delay     = pd.DataFrame(fr_delay)
    sample_array = pd.concat([sample_array, fr_sample], axis=1)
    choice_array = pd.concat([choice_array, fr_choice], axis=1)
    delay_array  = pd.concat([delay_array, fr_delay], axis=1)
    
# change to numpy
sample_array = np.array(sample_array)
choice_array = np.array(choice_array)
delay_array  = np.array(delay_array)

# remove nans
sample_array = np.delete(sample_array, 0, 1)
choice_array = np.delete(choice_array, 0, 1)
delay_array  = np.delete(delay_array, 0, 1)

# concatenate data
fr_data = {"mean_sample" : mean_sample,
           "mean_choice" : mean_choice,
           "mean_delay"  : mean_delay,
           "array_sample": sample_array,
           "array_choice": choice_array,
           "array_delay" : delay_array}
    
