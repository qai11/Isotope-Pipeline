"""
Title: bcorr.py
Author: Quin Aicken Davies, Heather Sinclair-Wentworth, Ethan Bull
Date: 11/03/24

Description: Takes the '.mat' files from Megara and loads them into the pipeline. 
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
from Barycentric_correction_iSpec import calculate_barycentric_velocity_correction

def bcorr(long_array, arrays_with_names): 
    if not arrays_with_names['bcorr'].ravel()[0].item():
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
        
    return corrected_long_array