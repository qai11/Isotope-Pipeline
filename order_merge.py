"""
Title: order_merge.py
Author: Quin Aicken Davies, Heather Sinclair-Wentworth, Ethan Bull
Date: 11/03/24

Description: This module takes column_list, a list of numpy arrays, each containing 
a wavelength and intensity column, and merges the overlapping regions of the spectra. 
The overlapping regions are weighted by the sigmoid function of the wavelength values, 
and the weighted average is calculated. The weighted average is then added to the smoothed 
spectrum. The module returns a single numpy array containing the merged spectra.
"""
from astropy.io import fits
import numpy as np
from scipy.io import loadmat
import scipy
import matplotlib.pyplot as plt
import glob
import os
import warnings
import matplotlib.pyplot as plt
import pandas as pd

#Defining the Sigmoid function for weighted average
def sigmoid(x): #x is an array
    delta_x = max(x)-min(x)
    new_x = (x-min(x))-delta_x/2
    return (1 / (1 + np.exp(-new_x)))

def order_merge(column_list):
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
    
    return long_array