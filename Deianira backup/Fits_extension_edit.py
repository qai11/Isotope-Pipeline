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
import scipy
import matplotlib.pyplot as plt
import glob
import os
from Barycentric_correction_iSpec import calculate_barycentric_velocity_correction
import matplotlib.pyplot as plt
import pandas as pd
#%%

# current_path = os.path.dirname(os.path.abspath(__file__))

'''For opening the data and showing the arrangement of the headers'''

#Checks if the star name is already defined
star_name = 'hd_100407'
if star_name == None:
    star_name = input('Enter the name of the star: ')
    
#Checks if the reduced data path is already defined    
reduced_data = '/home/users/qai11/Documents/megara_pipeline/REDUCED_DATA/'
if reduced_data == None:
    reduced_data = input('Enter the path to the reduced data: ')

#Checks if the post reduced data path is already defined
post_reduced_data = '/home/users/qai11/Documents/Fixed_fits_files/'
if post_reduced_data == None:
    post_reduced_data = input('Enter the path to the reduced data: ')
path = reduced_data + star_name +'/'


# Find all files in the path that match the pattern *_reduced.fits
file_pattern = path + 'J*.mat'
print(file_pattern)
reduced_files = glob.glob(file_pattern)

# Print the list of matching files         

print(reduced_files)

#Defining the Sigmoid function for weighted average
def sigmoid(x): #x is an array
    delta_x = max(x)-min(x)
    new_x = (x-min(x))-delta_x/2
    return (1 / (1 + np.exp(-new_x)))

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

    # Create a list to store the matlab input for accessing later
    padded_arrays = []

    # Pad each array to match the maximum length along axis 1
    for key in j_data.dtype.names:
        array = j_data[key]
        pad_length = max_length - array.shape[1]
        padded_array = np.pad(array, ((0, 0), (0, pad_length)), mode='constant')
        padded_arrays.append(padded_array)

        # Extract the names from the .mats for headers later
        names = [key for key in j_data.dtype.names]

        #Create an empty dictionary to store the arrays along with their names
        arrays_with_names = {}

        # Iterate over the padded arrays and their corresponding names
        for array, name in zip(padded_arrays[94:], names[94:]):
            # Add the array as a value to the dictionary with the name as the key
            arrays_with_names[name] = array
        
    # Concatenate the padded arrays along axis 1
    concatenated_data = np.concatenate(padded_arrays[0:94], axis=1)

    '''Puts the padded arrays (matlab data in a cube form) into two columns per cell'''
    column_list = []
    for i in range(10,15):
    # for i in range(len(padded_arrays[0:94])):#this should be 94 when running fully
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
            order_1_overlap = order[order[:,0] >= wave_min] #values from Order 1 that overlap order 2
            order_2_overlap = next_order[next_order[:,0] <= wave_max] #values from Order 2 that overlap order 1
            #Currently this isnt working so comment out for now
            # # Calculate the signal-to-noise ratio (S/N) in wavelength bins through the overlap region
            # bin_size = 0.1  # Define the size of each wavelength bin
            # overlap_bins = np.arange(wave_min, wave_max, bin_size)  # Create an array of bin edges
            # snr_bins = []  # List to store the S/N values for each bin

            # # Iterate over the wavelength bins
            # for bin_start in overlap_bins:
            #     bin_end = bin_start + bin_size
            #     bin_data = order_1_overlap[(order_1_overlap[:, 0] >= bin_start) & (order_1_overlap[:, 0] < bin_end)]
            #     if len(bin_data) > 0:
            #         snr = np.mean(bin_data[:, 1]) / np.std(order_2_overlap[:, 1])
            #         snr_bins.append(snr)

            # # Print the S/N values for each bin
            # print("Signal-to-Noise Ratio (S/N) in wavelength bins:")
            # print(snr_bins)
            
            # # Calculate the weights using the formula
            # weights = 1 / (np.square(order_1_overlap[:, 1]) * np.square(sigmoid(order_1_overlap[:, 0])))

            # # Apply the weights to the intensity values
            # weighted_intensity = np.multiply(weights, order_1_overlap[:, 1])

            # # Calculate the weighted average
            # weighted_average = np.sum(weighted_intensity) / np.sum(weights)

            # # Print the weighted average
            # print("Weighted Average:", weighted_average)
            
        
            # Calculate the weighted average using sigmoid function of our y values
            weighted_average = np.multiply(sigmoid(order_1_overlap[:, 0]), order_1_overlap[:, 1])
            column_list[i+1] = next_order[next_order[:,0] > wave_max]
            smooth_order = order[order[:,0] < wave_min]
            #separates the wavelength and intensity columns for adding onto
            smooth_wavelength = smooth_order[:,0]
            smooth_intensity = smooth_order[:,1]
            #adds on the weighted average and the overlapping region of order 1
            smooth_intensity = np.append(smooth_intensity, weighted_average)
            smooth_wavelength = np.append(smooth_wavelength, order_1_overlap[0])
            #combines them all back into the same file
            smooth_order = np.array([smooth_wavelength, smooth_intensity])
            #replaces the column list with the new spectrum
            column_list[i] = smooth_order
        except:
            None      
   

    # Concatenate all arrays into one long array
    long_array = np.concatenate(column_list)
    
    if star_name == "Cal_Dark":
        """If dark file, no barycentric correction is needed."""
        barycorr = 0
        c_light = 299792
        d_lambda = barycorr * (long_array[:, 0] / c_light)
        corrected_long_array = long_array[:, 0] + d_lambda

        # Convert the list of corrected wavelengths to a NumPy array
        corrected_long_array = corrected_long_array.T  
    elif not arrays_with_names['bcorr'].ravel()[0].item():
        'Calls Barycentric correction calculation if not in the .mat file'             
        datetime = [arrays_with_names['expdate'].ravel()[0].item()]
        ra = arrays_with_names['RA_j2000'].ravel()[0].item()
        dec = arrays_with_names['DEC_j2000'].ravel()[0].item()
        coordinates = [j for k in [[float(x) for x in ra.split()], [float(x) for x in dec.split()]] for j in k] #this prob wont work check later lol
        barycorr = calculate_barycentric_velocity_correction(datetime, coordinates, deq=0)
        print(barycorr)
        c_light = 299792  # Speed of light in km/s
        d_lambda = barycorr * (long_array[:, 0] / c_light)
        corrected_long_array = long_array[:, 0] + d_lambda

        # Convert the list of corrected wavelengths to a NumPy array
        corrected_long_array = corrected_long_array.T
    else:
        """If the barycentric correction is in the .mat file, use that."""
        #Find and add in the barycentric correction into the dataset
        c_light = 299792  # Speed of light in km/s
        barycorr = arrays_with_names['bcorr'].ravel()[0].item()  # Barycentric correction in km/s
        d_lambda = barycorr * (long_array[:, 0] / c_light)
        corrected_long_array = long_array[:, 0] + d_lambda

        # Convert the list of corrected wavelengths to a NumPy array
        corrected_long_array = corrected_long_array.T
            
    '''create the fits files'''
    
    '''Create the primary extention header'''
    header_keywords = {
        'JD': ('jd', 'Julian Date'),
        'SNR': ('signal_to_noise', 'Signal-to-Noise Ratio'),
        'OBJECT': ('starname', 'Star Name'),
        'FWMT': ('FWMT', 'Flux Weighted Mean Time'),
        'HERCFTC': ('counts', 'Flux meter total counts'),
        'FIBRE': ('fibre', 'HERCULES fibre position'),
        'EXPTIME': ('exptime', 'Exposure Time in seconds'),
        'HERCEXPT': ('exptype', 'Exposure Type'),
        'DATE-OBS': ('expdate', 'Exposure Date yyyy-mm-dd'),
        'TIME-OBS': ('rec_start', 'Recorded Start Time in UTC hh:mm:ss.s'), 
        'MIDPOINT': ('midtime_utc','Midtime of exposure in UTC hh:mm:ss.sss'),
        'MIDTT': ('midtime_tt','Midtime of exposure in Terrestrial Time hh:mm:ss.sss'),
        'FILENAME': ('sname','File name'),
        'OBJCTRA': ('RA_j2000', 'Right Ascension (J2000) hh:mm:ss.s'),
        'OBJCTDEC': ('DEC_j2000', 'Declination (J2000) dd:mm:ss'),
        'BLUECHOP': ('blue_data_chop', 'Blue Data Chop'),
        'BCORR': ('bcorr', 'Barycentric Correction in km/s')
    }
    
    for keyword, name in header_keywords.items():
        try:
            value = arrays_with_names[name].ravel()[0].item()
            if keyword == 'expdate' or keyword == 'utc_mid' or keyword == 'tt_mid':
                value = value[1:-1]
            new_hdul[0].header.set(keyword, value)
        except:
            None

    '''Creating the secondary header datatable'''
    # Create a new table HDU with columns 'col1' and 'col2' to hold the array data
    cols = []
    cols.append(fits.Column(name='wave', format='D', array=corrected_long_array/10))
    cols.append(fits.Column(name='flux', format='D', array=long_array[:, 1]))
    coldefs = fits.ColDefs(cols)
    new_table_hdu = fits.BinTableHDU.from_columns(coldefs)

    # Add the new table HDU to the HDUList (This is for adding the table into he secondary header)
    new_hdul.append(new_table_hdu)
    '''Try and make a directory for them to be saved in'''
    save_folder = post_reduced_data + star_name + '/'
    try:
        new_folder = os.mkdir(save_folder)
    except:
        None
    '''Location to save the file'''
    split_file_name = save_folder + '/' + file[-12:-4] + '.fits'

    '''Write to file'''
    new_hdul.writeto(split_file_name,overwrite=True)

    new_hdul.close()
    # %%
