"""
Title: create_cube.py
Author: Quin Aicken Davies, Heather Sinclair-Wentworth, Ethan Bull
Date: 11/03/24

Description: This module reshapes the data from the '.mat' files into a cube form for use in the
Deianira pipeline.
"""
from astropy.io import fits
import numpy as np
import scipy
import matplotlib.pyplot as plt
import glob
import os
import warnings
import matplotlib.pyplot as plt
import pandas as pd
from scipy.io import loadmat

def create_cube(file, j_data, max_length):
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
            
    return column_list, arrays_with_names