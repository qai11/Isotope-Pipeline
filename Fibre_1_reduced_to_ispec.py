#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 11:35:17 2024

@author: qai11
NOTE: THIS IS FROM THE hello COMIT ON THE GITLAB
"""
#%%
from astropy.io import fits
import numpy as np
from scipy.io import loadmat
#%%
'''For opening the data and showing the arrangement of the headers'''

path = '/home/users/qai11/Documents/megara_pipeline/REDUCED_DATA/hd_88158/'

file_name = 'J0006029_reduced.fits'

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

# Create a new table HDU with columns 'col1' and 'col2' to hold the array data
cols = []
cols.append(fits.Column(name='wave', format='D', array=data[:, 0]))
cols.append(fits.Column(name='flux', format='D', array=data[:, 1]))
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
