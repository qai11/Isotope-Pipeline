#%%
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
from create_cube import create_cube
from bcorr import bcorr
from create_fits import create_fits
from order_merge import order_merge
from open_mats import open_mats


warnings.filterwarnings("ignore")
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
reduced_files = glob.glob(file_pattern)

for file in reduced_files:
    '''Load directly from the unmerged orders in .mat'''
    j_data, max_length = open_mats(file)
    
    '''Create the cube'''
    column_list, arrays_with_names = create_cube(file, j_data, max_length)
    
    '''Merge the orders'''
    long_array = order_merge(column_list)
    
    '''Barycentric correction'''
    corrected_long_array = bcorr(long_array, arrays_with_names)
    
    '''Create the fits files'''
    create_fits(arrays_with_names,corrected_long_array, long_array, post_reduced_data, star_name, file)
    




# %%
