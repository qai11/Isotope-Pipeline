"""
Title: create_fits.py
Author: Quin Aicken Davies, Heather Sinclair-Wentworth, Ethan Bull
Date: 11/03/24

Description: This module is used to create the fits files from the data that has been reduced
and reconfigured. It takes the header information from the '.mat' files and adds them to
the primary header. The data is then added to the secondary header and saved as a '.fits' file.
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

def create_fits(arrays_with_names,corrected_long_array, long_array, post_reduced_data, star_name, file):
    '''create the fits files''' 
    '''Creates a new fits file for adding in the header information'''
    new_hdul = fits.HDUList()
    new_hdul.append(fits.ImageHDU())

    hdu1=fits.PrimaryHDU()
    new_hudl=fits.HDUList([hdu1])
    
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