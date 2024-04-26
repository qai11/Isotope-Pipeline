#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 11:35:17 2024

@author: Quin Aicken Davies
"""
#%%
from astropy.io import fits
import numpy as np
from scipy.io import loadmat
import matplotlib.pyplot as plt
import glob
import os
#%%
'''For opening the data and showing the arrangement of the headers'''

# current_path = os.path.dirname(os.path.abspath(__file__))

star_name = 'hd_45588'
if star_name == None:
    star_name = input('Enter the name of the star: ')

'''For opening the data and showing the arrangement of the headers'''

path = '/home/users/qai11/Documents/megara_pipeline/REDUCED_DATA/' + star_name +'/'


# Find all files in the path that match the pattern *_reduced.fits
file_pattern = path + 'J*.mat'
print(file_pattern)
reduced_files = glob.glob(file_pattern)

# Print the list of matching files         

print(reduced_files)

#%%
for file in reduced_files:
    '''Creates a new fits file for adding in the header information'''
    new_hdul = fits.HDUList()
    new_hdul.append(fits.ImageHDU())

    hdu1=fits.PrimaryHDU()
    new_hudl=fits.HDUList([hdu1])

    '''Load directly from the unmerged orders in .mat'''

    # Load the .mat file
    mat_data = loadmat(file)

    # Extract the data from the 'J J0353028' key
    j_data = mat_data[file[-12:-4]]

    # Determine the maximum length along axis 1
    max_length = max([j_data[key].shape[1] for key in j_data.dtype.names])

    # Create a list to store the padded arrays
    padded_arrays = []

    # Pad each array to match the maximum length along axis 1
    for key in j_data.dtype.names:
        array = j_data[key]
        pad_length = max_length - array.shape[1]
        padded_array = np.pad(array, ((0, 0), (0, pad_length)), mode='constant')
        padded_arrays.append(padded_array)

    # Concatenate the padded arrays along axis 1
    concatenated_data = np.concatenate(padded_arrays[:-9], axis=1)

    '''Puts the padded arrays (matlab data in a cube form) into two columns per cell'''
    column_list = []
    for i in range(len(padded_arrays[:-9])):#this should be 94 when running fully
        array_data = concatenated_data[0, i]
        if array_data.shape[1] >= 2:  # Check if array has at least 2 columns
            # Extract the first two columns
            first_column = array_data[:, 0]
            second_column = array_data[:, 1]
            # Combine the first and second columns into a single array
            combined_columns = np.column_stack((first_column, second_column)) 
            column_list.append(combined_columns)
            
    #Takes the column list, find where the spectra overlap and cut them both in between.
    for i, order in enumerate(column_list):
        try:
            next_order = column_list[i+1]#Defining the second order
            wave_max = order[-1,0]
            wave_min = next_order[0,0]
            wave_split = (wave_max - wave_min)/2 #Finds the middle of the overlap
            column_list[i] = order[order[:,0] < wave_min + wave_split]#Cuts the end 
            column_list[i+1] = next_order[next_order[:,0] > wave_max - wave_split]#Cuts the beginning 
        except:
            None

    # Concatenate all arrays into one long array
    long_array = np.concatenate(column_list)
    #Find and add in the barycentric correction into the dataset
    c_light = 299792  # Speed of light in km/s
    barycorr = float(padded_arrays[-1])  # Barycentric correction in km/s
    d_lambda = barycorr * (long_array[:, 0] / c_light)
    corrected_long_array = long_array[:, 0] + d_lambda

    # Convert the list of corrected wavelengths to a NumPy array
    corrected_long_array = corrected_long_array.T
    
    # Create a new table HDU with columns 'col1' and 'col2' to hold the array data
    cols = []
    cols.append(fits.Column(name='wave', format='D', array=corrected_long_array))
    cols.append(fits.Column(name='flux', format='D', array=long_array[:, 1]))
    coldefs = fits.ColDefs(cols)
    new_table_hdu = fits.BinTableHDU.from_columns(coldefs)

    # Add the new table HDU to the HDUList
    new_hdul.append(new_table_hdu)
    try:
        new_folder = os.mkdir('/home/users/qai11/Documents/Fixed_fits_files/' + star_name)
    except:
        None
    '''Location to save the file'''
    split_file_name = '/home/users/qai11/Documents/Fixed_fits_files/' + star_name + '/' + file[-12:-4] + '.fits'

    '''Write to file'''
    new_hdul.writeto(split_file_name,overwrite=True)

    new_hdul.close()
    # %%
