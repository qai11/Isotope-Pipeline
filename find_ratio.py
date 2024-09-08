"""
Title: find_ratio.py
Author: Quin Aicken Davies
Date: 22/08/24

Description: This is the main script for the find_ratio project. 
This uses parameters found from find_params.py and find_abund.py
to calculate the ratio of the isotopes of interest. It requires 
the input of the fits files from the previous scripts.
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

#--- iSpec directory -------------------------------------------------------------
ispec_dir = '/home/users/qai11/iSpec_v20201001'
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec

#%%

def make_temp_file(filename):
    '''This function will create a temporary file for use in MOOG outputs'''
    f  = open(filename, "a+") 
    f.write('')
    f.close() 

def create_moog_input(spectrum, in_file, out_file, wave_region):
    '''This function will create the input file for MOOG'''
    standard_output = 'out1'
    summary_output = 'out2'
    
    make_temp_file(standard_output)
    make_temp_file(summary_output)
    make_temp_file(out_file)

    par_string = "synth\n" +\
    "standard_out   '" + standard_output +"'\n"                    + \
    "summary_out    '" + summary_output +"'\n"                     + \
    "smoothed_out   '" + out_file +"'\n"                    + \
    "model_in       'qaimodel.moog'\n"                          + \
    "lines_in       'quinlist.MgH'\n"                           + \
    "observed_in    '" + spectrum +"'\n"               + \
    "atmosphere    1\n"                                         + \
    "molecules     2\n"                                         + \
    "lines         2\n"                                         + \
    "flux/int      0\n"                                         + \
    "plotpars      1\n"                                         + \
    wave_region + " 0.15 1.05\n"                          + \
    str(par['rv']) + "      0.000   0.000    1.00\n"                   + \
    "d          0.047 0.0 0.0 "+ str(par['s']) +" 0.0\n"        + \
    "abundances   4    1\n"                                     + \
    "6            0.0500000\n"                                  + \
    "12           " + str(par['mg']) + "\n"                     + \
    "24           0.50000\n"                                    + \
    "26           0.20000\n"                                    + \
    "isotopes      4    1\n"                                    + \
    "606.01212     5.0\n"                                       + \
    "112.00124     "+ str(par['i_24']) +"\n"                    + \
    "112.00125     "+ str(par['i_25']) +"\n"                    + \
    "112.00126     "+ str(par['i_26']) +"\n"                    + \
    "obspectrum 5\n"                                            + \
    "synlimits\n"                                               + \
    wave_region + " 0.01 5.0\n"                           + \
    "plot 2\n"                                                  + \
    "damping 2\n"


    # writing that string to a file 
    par_file  = open(in_file, "w+") 
    par_file.write(par_string)
    par_file.close() 
    return in_file, out_file



start = time.time()

star = ['hd_45588','hd_100407','hd_102870','hd_128620','hd_128621','hd_11695','hd_146233','hd_156098','hd_157244','hd_160691','moon']

#--- Read in the fits files -----------------------------------------------------
# Define the path to the fits files
fits_path = '/home/users/qai11/Documents/Fixed_fits_files/hd_160691/'

# Read in the fits files
fits_files = os.listdir(fits_path)
fits_files = [fits_path + file for file in fits_files]

# Read in the fits files
data = ispec.read_fits_files(fits_files)

#--- Read in the data -----------------------------------------------------------
# Define the path to the data
data_path = '/home/users/qai11/Documents/Fixed_fits_files/hd_160691/'

# Read in the data
data = ispec.read_data(data_path)

#--- Define the isotopes of interest -------------------------------------------
# Define the isotopes of interest
isotopes = ['Mg24', 'Mg25', 'Mg26']

#--- Define the lines of interest ---------------------------------------------
# Define the lines of interest
lines = [4703.0, 5528.4]

#--- Define the parameters ----------------------------------------------------
# Define the parameters

for star_name in star:
    star_spectrum = ispec.read_spectrum(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/{star_name}_adjusted.fits')
    parameters = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/iteration_parameters/{star_name}_iter_param.txt')
    teff = parameters['teff'][len(parameters)-1]
    logg = parameters['logg'][len(parameters)-1]
    MH = parameters['MH'][len(parameters)-1]
    Vmic = parameters['Vmic'][len(parameters)-1]
    Vmac = parameters['Vmac'][len(parameters)-1]
    alpha = parameters['alpha'][len(parameters)-1]
    vsini = parameters['Vsini'][len(parameters)-1]
element_name = 'Mg'
wave_base = 480
wave_top = 680
resolution = 82000
code = "moog"

#--- Calculate the ratio -------------------------------------------------------
# Calculate the ratio
# ratio = ispec.determine_abundances_using_synth_spectra(data, teff, logg, MH, Vmic, Vmac, alpha, vsini, max_iterations, element_name, wave_base, wave_top, resolution, code)

#--- Plot the ratio ------------------------------------------------------------
# Plot the ratio
plt.plot(ratio)
plt.xlabel('Wavelength')
plt.ylabel('Ratio')
plt.title('Ratio of Isotopes') 
plt.show()

#--- Save the ratio ------------------------------------------------------------
# Save the ratio
# ratio_path = '/home/users/qai11/Documents/Fixed_fits_files/hd_160691/ratio.csv'
ratio.to_csv(ratio_path)

#--- End of script ------------------------------------------------------------
end = time.time()

print(f'Time taken in sec: {end - start}')
#%%