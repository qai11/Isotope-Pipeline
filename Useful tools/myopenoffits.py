#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 12:53:45 2023

@author: qai11
"""

from astropy.io import fits
import os
#%%
star_name = 'J0354019.fit'
path = '/home/users/qai11/data/MATLAB/Version_1_5_Karen_Linux/RAW_DATA/20240213/'+star_name

hdul = fits.open(path)
list1 = hdul[0].header
list1.keys
#%%
'''edits the header and saves it'''

list1['HERCEXPT']=' Stellar'
hdul.writeto(star_name,overwrite=True)
hdul.close()

    

#%%
'''Reads certain parts of the header'''
path2 = '/media/astro/MJArchive/octans/20240224'



fitsfiles=[f for f in os.listdir(path2) if f.endswith('.fit')]
for i in fitsfiles:
    fitspath=os.path.join(path2, i)
    header = fits.getheader(fitspath, ignore_missing_end=True)
    print(header['Filter'])
    print(header['EXPTIME'])
keys = hdul[0].header.keys()

#%%

hdul = fits.open(path)
list1 = hdul[0].header
print(list1)
print(list1['OBJECT'])

#%%
'''script for fixing megara fits files'''

star='hd_'+'100407'

path3 = '/home/users/qai11/Documents/megara_pipeline/Reduced_Data/'+star+'/'


fitsfiles=[f for f in os.listdir(path3) if f.endswith('reduced.fits')]
for i in fitsfiles:
    fitspath=os.path.join(path3, i)
    hdul=fits.open(fitspath)
    header = fits.getheader(fitspath, ignore_missing_end=True)
    start_time = header['EXP_TIME']#saves the data that needs to be swapped
    exp_time = header['START']
    header['EXP_TIME']= exp_time #changes the data that needs to be swapped
    header['START'] = start_time
    hdul.writeto(i,overwrite=True)#writes all of it to file
    hdul.close()

