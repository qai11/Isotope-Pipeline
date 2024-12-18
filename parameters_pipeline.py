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
# star = ['hd_157244']
# star = ['hd_102870']
# star = ['hd_45588']
# star = ['hd_128620']
# star = ['hd_45588','hd_100407','hd_102870','hd_128620','hd_128621','hd_146233','hd_157244','hd_160691','moon']
# star = ['hd_2151','hd_11695','hd_18907','hd_10700','hd_23249','hd_22049','hd_18884','hd_165499','hd_156098']
star =['hd_45588','hd_100407','hd_102870','hd_128620','hd_128621','hd_146233','hd_157244','hd_160691','moon','hd_2151','hd_11695','hd_18907','hd_10700','hd_23249','hd_22049','hd_18884','hd_165499','hd_156098']
# star = ['hd_2151']
element = ["Mg", "Si", "Ca", "Ti", "Sc","V","Cr","Mn","Co", "Ni", "Y", "Ba", "La", "Nd", "Eu", "Sr", "Zr"]#Rb isnt included
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
# try:
#     pool = Pool(os.cpu_count()-1)
#     pool.map(partial(find_abundance, elements=element), star)
# finally:
#     pool.close()
#     pool.join()
#Run the interpolation of the atmospheres
# interp_atmos(star)

#%%
# '''Caculate the final abundances from the line by line abundances'''
# import pandas as pd

# # Define the list of stars and elements
# # stars = ['hd_45588']

# # Create an empty dataframe with the desired columns
# new_abund_df = pd.DataFrame(columns=['element', 'code', 'Abund', 'A(X)', '[X/H]', '[X/Fe]', 'eAbund', 'eA(X)', 'e[X/H]', 'e[X/Fe]'])

# # Iterate over each star and element
# for star_name in star:
#     for el in element:
#         # Read the data for the current star and element
#         line_abund = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/{star_name}/line_by_line_{el}.txt', delimiter=' ')
        
#         # If there are at least 3 data points, use median and mad, otherwise use mean and mad
#         if len(line_abund) >= 3:
#             median_XH = line_abund['[X/H]'].median()
#             mad_XH = line_abund['[X/H]'].mad()
#             median_XFe = line_abund['[X/Fe]'].median()
#             mad_XFe = line_abund['[X/Fe]'].mad()
#             median_Abund = line_abund['Abund'].median()
#             mad_Abund = line_abund['Abund'].mad()
#             median_A = line_abund['A(X)'].median()
#             mad_A = line_abund['A(X)'].mad()
#         else:
#             median_XH = line_abund['[X/H]'].mean()
#             mad_XH = line_abund['[X/H]'].mad()
#             median_XFe = line_abund['[X/Fe]'].mean()
#             mad_XFe = line_abund['[X/Fe]'].mad()
#             median_Abund = line_abund['Abund'].mean()
#             mad_Abund = line_abund['Abund'].mad()
#             median_A = line_abund['A(X)'].mean()
#             mad_A = line_abund['A(X)'].mad()
        
#         # Extract the code for the current element
#         code = line_abund['code'].iloc[0]
        
#         # Append the calculated values to the dataframe
#         new_abund_df = new_abund_df.append({
#             'element': el, 
#             'code': code, 
#             'Abund': median_Abund, 
#             'A(X)': median_A, 
#             '[X/H]': median_XH, 
#             '[X/Fe]': median_XFe, 
#             'eAbund': mad_Abund, 
#             'eA(X)': mad_A, 
#             'e[X/H]': mad_XH, 
#             'e[X/Fe]': mad_XFe
#         }, ignore_index=True)

#     # Save the dataframe to a txt file for the current star
#     df.to_csv(f'/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/{star_name}/summary_abundances_{star_name}.txt', sep=' ', index=False)

#     print(f"Data for {star_name} saved successfully.")

# star = ['hd_157244']
#Load and save the stars with a cut wavelength range for MgH fitting
for star_name in star:
    try:
        #Uni PC
         spectrum = ispec.read_spectrum(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/{star_name}_adjusted.fits')
    except:
        #MAC
        spectrum = ispec.read_spectrum(f'/Users/quin/Desktop/2024_Data/Fixed_fits_files/{star_name}/{star_name}_adjusted.fits')
        
    spectrum['waveobs'] = spectrum['waveobs']
    #Cut the spectrum to be between 5100 and 5200
    wfilter = ispec.create_wavelength_filter(spectrum, wave_base=510, wave_top=520)
    cutted_star_spectrum = spectrum[wfilter]
    #Save the cut spectrum
    try:
        #Uni PC
        ispec.write_spectrum(cutted_star_spectrum,f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/{star_name}_510-520.txt')
    except:
        #MAC
        ispec.write_spectrum(cutted_star_spectrum,f'/Users/quin/Desktop/2024_Data/Fixed_fits_files/{star_name}/{star_name}_510-520.txt')
    print(f"{star_name} saved successfully.")
    
# end = time.time()

# print(f'Time taken: {(end - start)/3600} Hrs')
# print('Pipeline Complete done!')


# %%
