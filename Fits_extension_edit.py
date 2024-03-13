#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 11:35:17 2024

@author: qai11
"""

from astropy.io import fits
import numpy as np
from scipy.io import loadmat
#%%
'''For opening the data and showing the arrangement of the headers'''

path = '/home/users/qai11/Documents/megara_pipeline/Reduced_Data/hd_100407/'

file_name = 'J0353028_reduced.fits'

hdul = fits.open(path + file_name)
list1 = hdul[0].header#shows the information in the header

hdul.info()
#%%

data = hdul[0].data

'''Creates a new fits file for adding in the header information'''
new_hdul = fits.HDUList()
new_hdul.append(fits.ImageHDU())

hdu1=fits.PrimaryHDU()
new_hudl=fits.HDUList([hdu1])

#%%
'''do it directly from the unmerged orders'''

# Load the .mat file
mat_data = loadmat(path +'J0353028.mat')

# Extract the data from the 'J J0353028' key
j_data = mat_data['J0353028']

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

column_list = []

for i in range(94):
    array_data = concatenated_data[0, i]
    if array_data.shape[1] >= 2:  # Check if array has at least 2 columns
        # Extract the first two columns
        first_column = array_data[:, 0]
        second_column = array_data[:, 1]
        # Combine the first and second columns into a single array
        combined_columns = np.column_stack((first_column, second_column))
        column_list.append(combined_columns)

# Concatenate all arrays into one long array
long_array = np.concatenate(column_list)


#%%

# Create a new table HDU with columns 'col1' and 'col2' to hold the array data
cols = []
cols.append(fits.Column(name='wave', format='D', array=long_array[:, 0]))
cols.append(fits.Column(name='flux', format='D', array=long_array[:, 1]))
coldefs = fits.ColDefs(cols)
new_table_hdu = fits.BinTableHDU.from_columns(coldefs)

# Add the new table HDU to the HDUList
new_hdul.append(new_table_hdu)

#%%
'''Location to save the file'''
split_file_name = '/home/users/qai11/Documents/Fixed_fits_files/' + file_name

'''Write to file'''
new_hdul.writeto(split_file_name,overwrite=True)


#%%
hdul.close()