"""
Title: open_mats.py
Author: Quin Aicken Davies, Heather Sinclair-Wentworth, Ethan Bull
Date: 11/03/24

Description: Takes the '.mat' files from Megara and loads them into the pipeline. 
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

'''For opening the data and showing the arrangement of the headers'''
def open_mats(file):
    '''Load directly from the unmerged orders in .mat'''

    # Load the .mat file
    mat_data = loadmat(file)

    # Extract the data from the 'J J0353028' key
    j_data = mat_data[file[-12:-4]]

    # Determine the maximum length along axis 1
    max_length = max([j_data[key].shape[1] for key in j_data.dtype.names])
    
    return j_data, max_length