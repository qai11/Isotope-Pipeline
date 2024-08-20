"""
Title: rv_combine.py
Author: Ethan Bull and Quin Aicken Davies
Date: 30/07/2024

Description: Uses the iSpec to rv correct spectrum and combine them
"""
#%%
import numpy as np
import glob
import os
from astropy.io import fits
import pandas as pd
import time
import scipy as sp
import sys

#--- iSpec directory -------------------------------------------------------------
#ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/"
ispec_dir = '/home/users/qai11/iSpec_v20201001/'
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec

start = time.time()

# spectrum = ispec.read_spectrum(f'{data_path}/{star_name}/median_spectrum_{star_name}.txt')
# star = ['hd_45588','hd_100407','hd_102870','hd_128620','hd_128621','hd_11695','hd_146233','hd_156098','hd_157244','hd_160691','moon']
star = ['hd_11695','hd_146233','hd_156098','hd_157244','hd_160691','moon']

for star_name in star:
    # Open spectra one at a time in a folder
    folder_path = f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}'  # Replace with the actual folder path
    files = glob.glob(folder_path + '/*.fits')

    start = time.time()
    rad_vel = []
    all_spectra = []
    for file in files:
        # Load the data
        spectrum = ispec.read_spectrum(file)
        # Radial velocity shift the spectrum for further analysis
        template = ispec.read_spectrum(ispec_dir + '/input/spectra/templates/NARVAL.Sun.370_1048nm/template.txt.gz')
        #Model the spectrum ccf
        model, ccf = ispec.cross_correlate_with_template(spectrum, template, lower_velocity_limit=-200, upper_velocity_limit=200, velocity_step=1.0, fourier=False)
        #Use the shift for the ccf to calculate the shift for the spectrum
        rv = np.round(model[0].mu(),2)
        #Add to a list for later
        rad_vel = np.append(rad_vel, rv)
        #Print the rv for the spectrum
        print(f'Radial Velocity correction: {rv} km/s')
        #Save the error
        rv_err = np.round(model[0].emu(),2)
        #Shift the Spectrum
        rv_shift_median_spectrum = ispec.correct_velocity(spectrum, rv)
        #Save the rv corrected spectrum
        #Make a new folder for the rv corrected spectra
        if not os.path.exists(folder_path + '/rv_corrected'):
            os.makedirs(folder_path + '/rv_corrected')
        #Make a dataframe for the rv corrected spectrum
        rv_shifted_data = pd.DataFrame({'waveobs': rv_shift_median_spectrum['waveobs'], 'flux': spectrum['flux'], 'err': spectrum['err']})
        #Save the shifted spectrum
        # rv_shifted_data.to_csv(folder_path + '/rv_corrected/corrected_' + file[-13:-5] + '.txt', sep=' ',index=False)
        
        '''MERGE THE DATA'''
        # Load the data into data and wave for median addition
        data = fits.getdata(file)
        all_spectra.append(spectrum['flux'])
    
    #median the fluxes
    median_flux = np.median(all_spectra, axis=0)
    # rebuild data for ispec
    median_data = pd.DataFrame({'waveobs': rv_shift_median_spectrum['waveobs'], 'flux': median_flux, 'err': np.zeros_like(median_flux)})
    # Save the median spectrum
    median_data.to_csv(folder_path + '/rv_corrected/median_spectrum_' + star_name + '.txt', sep=' ',index=False)
    
end = time.time()

print(f'Time taken: {end - start}')




# %%
