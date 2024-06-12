"""
Title: spectrum_plot.py
Author: Quin Aicken Davies
Date: 13/05/24

Description: This model inputs the .txt files from ispec and 
plots the spectra for the different stellar parameters.
"""
#%%
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
import numpy as np
import logging
import multiprocessing
from multiprocessing import Pool


#%%
#Define the path to the data and star information
# data_information = pd.read_csv(f'/home/users/qai11/Documents/quin-masters-code/Cool_stars.csv',sep=' ')
# data_information = pd.read_csv(f'/home/users/qai11/Documents/quin-masters-code/Hot_stars.csv',sep=' ')
#data_information = pd.read_csv(f'/home/users/qai11/Documents/quin-masters-code/Masters_stars.csv',sep=',')
data_information = pd.read_csv(f'/home/users/qai11/Documents/quin-masters-code/Sun.csv',sep=' ')

star_name = []
for i,name in enumerate(data_information['ID2']):
    name = name.lower()
    name = name[:2] + '_' + name[2:]
    star_name = np.append(star_name,name)

# star_name = ['Sun']
# star_name = ['Sun','hd_100407','hd_128621']
# star_name = ['Sun','hd_45588','hd_102870','hd_128620']
#%%
#prints every star with selected stellar parameters
indexes = []
plt.figure(figsize=(13,7))
for i,name in enumerate(star_name):
    try:
        path = f'/home/users/qai11/Documents/Fixed_fits_files/{name}/{name}.txt'  
        # Define the data
        data = pd.read_csv(path, sep='	')
        #Defining the vectors
        wavelength_vector = data['waveobs']
        intensity_vector = data['flux']
        plt.plot(wavelength_vector, intensity_vector,label =f"{data_information['ID2'][i]}, TEFF:{data_information['TEFF'][i]}, LOGG:{data_information['LOGG'][i]}, Fe/H:{data_information['FEH'][i]}, Spectral Type:{data_information['SPT'][i]}")
        indexes = np.append(indexes,i)
    except:
        continue

#plots the entire spectra
plt.legend(loc='lower left')
plt.xlabel('Wavelength (nm)',size=12)
plt.ylabel('Intensity',size=12)
plt.title('Spectra of stars with selected stellar parameters',size=15)   
#plt.xlim(513.4,514.1)
#plt.ylim(0.78,1.01)
#plots the vertical lines for the Mg isotope lines
# line_names = pd.read_csv('/home/users/qai11/iSpec_v20201001/input/regions/Mg_isotope_Quin.txt',sep='	') 
# for j,name in enumerate(line_names['wave_peak']):
#     plt.axvline(x=name,ymin=0,ymax=1, color='gold', linestyle='-')
#     plt.axvspan(xmin=line_names['wave_base'][j],xmax=line_names['wave_top'][j], color='gold', alpha=0.5)
plt.savefig('/home/users/qai11/Documents/RASNZ_Talks/Sun_spectra_full.png',dpi=300)
plt.show()



# %%
