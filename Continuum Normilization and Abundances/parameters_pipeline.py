"""
Title: parameters_pipeline.py
Author: Quin Aicken Davies
Date: 06/10/24

Description: This runs the code that does rv correction, merging, adjusts the continuum 
placement,generates the parameters then generates the metallities of all listed stars. 
This will allow direct input to interpolate atmospheres for each star and allow for 
isotopes to be estimated without having to manually input the parameters. 

Appropriate iSpec edits will be listed on the Isotopes wiki page in the near future.
"""
#%%
import os
import sys
import numpy as np
import logging
import multiprocessing
from multiprocessing import Pool
from matplotlib import pyplot as plt
import pandas as pd
from scipy.stats import chisquare
import scipy.optimize as opt
from copy import deepcopy
import time
from functools import partial

#%%
#This must be run in the correct order or it wont do what is needed
'''Note: The linelist must contain isotopes and the relevant molecules
for the abundances and isotopes to be found and be consistent.
Most things should be on parallel pool so should work quickly if all things go well.'''
#First Deianira must be run to completion first
#Then rv_combine.py
#Then continuum_adjust.py
#Then find_params.py
#Then find_abund.py
#Then interp_atmos.py
#Then the data can be run through to find the isotopes

'''At the moment I have the stars all hard coded in,
in order to change this we will have to remove that and add in a loop
once a star list is loaded in, an easy change which I can do if requested but
not needed at this moment.'''

'''Files are saved as they go so don't need to  be output as named variables'''

#%%
# Import the scripts to be run
from rv_combine import rv_combine
from continuum_adjust import continuum_adjust
from find_params import find_params
# from find_abund import find_abund
from line_by_line_abunds import find_abundance
from interp_atmos import interp_atmos

#--- iSpec directory -------------------------------------------------------------
if os.path.exists('/home/users/qai11/iSpec_v20201001'):
    "Location of the files on Uni computer"
    ispec_dir = '/home/users/qai11/iSpec_v20201001'
else:
    "location of data on Mac"
    ispec_dir = '/Users/quin/Desktop/2024_Data/iSpec_v20230804'
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec


# star = ['hd_128620','hd_157244','hd_102870']
# star = ['hd_18884']
# star = ['hd_102870']
# star = ['hd_45588']
# star = ['hd_128620']
# star = ['hd_45588','hd_100407','hd_102870','hd_128620','hd_128621','hd_146233','hd_157244','hd_160691','moon']
# star = ['hd_2151','hd_11695','hd_18907','hd_10700','hd_23249','hd_22049','hd_18884','hd_165499','hd_156098']
# star =['hd_45588','hd_100407','hd_102870','hd_128620','hd_128621','hd_146233','hd_157244','hd_160691','moon',
#        'hd_2151','hd_11695','hd_18907','hd_10700','hd_23249','hd_22049','hd_18884','hd_165499','hd_156098']

# star =['hd_100407','hd_102870','hd_128620','hd_128621','hd_146233','hd_157244','hd_160691','moon',
    #    'hd_2151','hd_18907','hd_10700','hd_23249','hd_22049','hd_18884']
  
# star = ['hd_156098']
  
# # star = ['hd_2151']
# star = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
#     'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
#     'hd_102870','hd_45588','hd_156098']
# element = ["Mg", "Si", "Ca", "Ti", "Sc","V","Cr","Mn","Co", "Ni", "Y", "Ba", "La", "Nd", "Eu", "Sr", "Zr"]#Rb isnt included
# element = ["Mg","Eu","Ba"]
star = ['hd_10700']
element = ["Mg"]
#%%
# start = time.time()
'''finish setting up parallel pools for each thing'''
#Run the rv correction and Merging 
# rv_combine(star)
#Run the continuum adjustment
# try:
#     pool = Pool(os.cpu_count()-1)
#     pool.map(continuum_adjust, star)
# finally:
#     pool.close()
#     pool.join()
#Run the parameter finding
# try:
#     pool = Pool(os.cpu_count()-1)
#     pool.map(find_params, star)
# finally:
#     pool.close()
#     pool.join()
# # Run the abundance finding
try:
    pool = Pool(os.cpu_count()-1)
    pool.map(partial(find_abundance, elements=element), star)
finally:
    pool.close()
    pool.join()
#Run the interpolation of the atmospheres
# interp_atmos(star)


#%%
# star = ['hd_157244']
# Load and save the stars with a cut wavelength range for MgH fitting
# for star_name in star:
#     try:
#         #Uni PC
#         #  spectrum = ispec.read_spectrum(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/{star_name}_adjusted.fits')
#         spectrum = ispec.read_spectrum(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/rv_corrected/median_spectrum_{star_name}.txt')
#     except:
#         #MAC
#         spectrum = ispec.read_spectrum(f'/Users/quin/Desktop/2024_Data/Fixed_fits_files/{star_name}/{star_name}_adjusted.fits')
        
#     spectrum['waveobs'] = spectrum['waveobs']*10
#     #Cut the spectrum to be between 5100 and 5200
#     wfilter = ispec.create_wavelength_filter(spectrum, wave_base=5100, wave_top=5200)
#     cutted_star_spectrum = spectrum[wfilter]
#     #Save the cut spectrum
#     try:
#         #Uni PC
#         ispec.write_spectrum(cutted_star_spectrum,f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/{star_name}_5100-5200.txt')
#     except:
#         #MAC
#         ispec.write_spectrum(cutted_star_spectrum,f'/Users/quin/Desktop/2024_Data/Fixed_fits_files/{star_name}/{star_name}_5100-5200.txt')
#     print(f"{star_name} saved successfully.")
    
# end = time.time()

# print(f'Time taken: {(end - start)/3600} Hrs')
# print('Pipeline Complete done!')


# %%
