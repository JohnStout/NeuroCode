# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 09:20:06 2020

@author: Jesse from Mizumori lab
"""

def calc_IdPhi(x_pos, y_pos):
    ''' 
    calculate IdPhi - 
    the absolute, integrated, change in angular motion
    '''
    diff_x = np.diff(x_pos)
    diff_y = np.diff(y_pos)
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