# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 09:57:52 2020

@author: uggriffin
"""

# libraries for opening matlab files
import scipy.io as sio
import os

# define a variable representing the session directory
sessions_folder = 'C:/Users/uggriffin/Documents/GitHub/CodingClub/'
session         = '/Data'
os.chdir(sessions_folder+session) # cd equivalent

# load Int and TimeVariables dictionaries
TimeVariables = sio.loadmat('VT1.mat')
IntVars       = sio.loadmat('Int_file.mat')

# If you want to isolate an element from the dictionary
ExtractedX = TimeVariables["ExtractedX"]
ExtractedY = TimeVariables["ExtractedY"]
TimeStamps = TimeVariables["TimeStamps"]
Int        = IntVars["Int"]

# convert data
import numpy as np                             # import numPy
numSamples = np.size(ExtractedX)               # how many samples total?
X_data     = np.resize(ExtractedX,(numSamples,))   # reformat
Y_data     = np.resize(ExtractedY,(numSamples,))
TimeStamps = np.resize(TimeStamps,(numSamples,))

"""
Remove instances where data is missing
"""
# find instances where x is zero
idxZero = []
idxZero = np.where(X_data == 0)

# remove the data
X_data = np.delete(X_data,idxZero)
Y_data = np.delete(Y_data,idxZero)
TimeStamps = np.delete(TimeStamps,idxZero)

# find instances where x is zero
idxZero = []
idxZero = np.where(Y_data == 0)

# remove the data
X_data = np.delete(X_data,idxZero)
Y_data = np.delete(Y_data,idxZero)
TimeStamps = np.delete(TimeStamps,idxZero)

"""
Get data around T-junction
"""
# plot data on a trial-by-trial basis
numTrials = np.size(Int,0)

# remember that python uses 0 indexing
ts_pos = [] # initialize variable
x_pos  = []
y_pos  = []

# for loop - get x, y, and times for every trial
for i in range(0,numTrials): # notice that syntax may change, but concept doesn't
    idx = np.where(np.logical_and(TimeStamps>Int[i,4], TimeStamps<Int[i,1]))
    ts_pos.append(TimeStamps[idx])
    x_pos.append(X_data[idx])
    y_pos.append(Y_data[idx])

# load in matplotlib library
import matplotlib.pyplot as plt # a plotting library

# sample trials
sample_trials = array=np.arange(0,numTrials,2)

# choice trials
choice_trials = array=np.arange(1,numTrials,2)

# plot
i = choice_trials[17]
plt.plot(x_pos[i],y_pos[i])

import scipy.stats as stat

# calc_IdPhi(x_pos[1], y_pos[1])


""" 
calculate IdPhi - 
the absolute, integrated, change in angular motion
"""

diff_x = np.diff(x_pos[1])
diff_y = np.diff(y_pos[1])

plt.plot(diff_x,diff_y)

turns = 0
# angular velocity of motion (Phi)
o_motion = np.arctan2(diff_y, diff_x)
# change in orientation of motion (dPhi)
dphi = np.diff(o_motion)
# integrated absolute dPhi
abs_dphi = np.absolute(dphi)
abs_int_dphi = np.sum(abs_dphi)
for i in abs_dphi:
  if i > 0.5:
    turns += 1 
return abs_int_dphi, stat.circstd(dphi), turns