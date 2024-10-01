"""add_error_spectra.py"""

#%%
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy.signal import find_peaks
import os
import glob
import pandas as pd
from astropy.io import fits

# raw_files_loc = '/home/users/ebu36/data/ASTR690/MEGARA/RAW_DATA/' # GBS Stars
raw_files_loc = '/home/users/qai11/data/MATLAB/Version_1_5_Karen_Linux/RAW_DATA/' # SAUCE Stars
reduced_files_loc = '/home/users/qai11/Documents/Fixed_fits_files/' #Change this and lines above

# make census of all raw files to refer to later
raw_files = glob.glob(raw_files_loc + '**/' + '*.fit*')
print(f'Total raw files: {len(raw_files)}')

# set up dataframe to track star name and location of each file
census = pd.DataFrame(columns=['star_name', 'file_name', 'file_loc'])
for file in raw_files:
   with fits.open(file) as hdul:
       star_name = hdul[0].header['OBJECT']
   census = pd.concat([census,pd.DataFrame({'star_name': star_name, 'file_name': file[file.rindex('/')+1:], 'file_loc': file}, index = [0])], ignore_index=True)


# get reduced files for a specified star; this should be case insensitive

def get_reduced_files(star):
   reduced_files = glob.glob(reduced_files_loc + star + '/*.fits')
   return reduced_files

def add_err_to_fits(star):
   star = star.lower()
   reduced_files = get_reduced_files(star)
   if len(reduced_files) == 0:
       reduced_files = get_reduced_files(star.upper())

   reduced_files = [x[x.rindex('/')+1:-1] for x in reduced_files]
   print(f'Total reduced files: {len(reduced_files)}')

   # open a raw file for each reduced file to measure the peaks of each order
   for file in reduced_files:
       raw_file = census[census['file_name'] == file]
       try:
            raw_file = raw_file['file_loc'].values[0]
       except:
           continue
       with fits.open(raw_file) as hdul:
           data = hdul[0].data
       # get the peaks of each order
       peaks_c, heights_c = find_peaks(data[2000, 1200:3900], height=data[0,0]*3, distance = 10, prominence=100) 
       peaks_c += 1200
       dark_peaks_c = []
       for i, peak in enumerate(peaks_c):
           try:
               dark_peaks_c.append(int((peak + peaks_c[i+1])/2))
           except:
               dark_peaks_c.append(int(peak+10))
       peaks_l, heights_l = find_peaks(data[1600, 1200:3900], height=data[0,0]*3, distance = 10, prominence=100) 
       peaks_l += 1200
       dark_peaks_l = []
       for i, peak in enumerate(peaks_l):
           try:
               dark_peaks_l.append(int((peak + peaks_l[i+1])/2))
           except:
               dark_peaks_l.append(int(peak+10))
       peaks_r, heights_r = find_peaks(data[2400, 1200:3900], height=data[0,0]*3, distance = 10, prominence=100) 
       peaks_r += 1200
       dark_peaks_r = []
       for i, peak in enumerate(peaks_r):
           try:
               dark_peaks_r.append(int((peak + peaks_r[i+1])/2))
           except:
               dark_peaks_r.append(int(peak+10))
           
       if len(peaks_c) != len(peaks_l) or len(peaks_c) != len(peaks_r) or len(peaks_l) != len(peaks_r):
           # find shortest and make all the same length
           peaks = [peaks_c, peaks_l, peaks_r]
           min_peaks = min(peaks, key=len)
           peaks_l = peaks_l[-len(min_peaks):]
           peaks_r = peaks_r[-len(min_peaks):]
           peaks_c = peaks_c[-len(min_peaks):]
           dark_peaks_l = dark_peaks_l[-len(min_peaks):]
           dark_peaks_r = dark_peaks_r[-len(min_peaks):]
           dark_peaks_c = dark_peaks_c[-len(min_peaks):]
           
       # figure to test
       # plt.figure()
       # plt.imshow(data)
       # plt.scatter(peaks_c, [2000]*len(peaks_c), c='r', s=3)
       # plt.scatter(dark_peaks_c, [2000]*len(dark_peaks_c), c='b', s=3)
       # plt.scatter(peaks_l, [1600]*len(peaks_l), c='r', s=3)
       # plt.scatter(dark_peaks_l, [1600]*len(dark_peaks_l), c='b', s=3)
       # plt.scatter(peaks_r, [2400]*len(peaks_r), c='r', s=3)
       # plt.scatter(dark_peaks_r, [2400]*len(dark_peaks_r), c='b', s=3)
       
       # calculate errors
       def error_func(peaks, dark_peaks):
           # print(peaks, dark_peaks)
           bias = 21
           error = np.sqrt(np.sum([[x - dark_peaks[i] + bias for i, x in enumerate(peaks)], np.square([x + bias - x for x in dark_peaks])], axis=0))
           error = np.ones(len(error))/error
           # print(error)
           return error
       
       errors_c = error_func(data[2000, peaks_c], data[2000,dark_peaks_c])
       errors_l = error_func(data[1600, peaks_l], data[1600,dark_peaks_l])
       errors_r = error_func(data[2400, peaks_r], data[2400,dark_peaks_r])
       errors = np.mean([errors_c, errors_l, errors_r], axis=0)
       error = np.mean(errors)
       
       # multiply flux values by error to get error column to save in fits
       # open reduced fits file
       try:
           with fits.open(reduced_files_loc + star + '/' + file + 's') as hdul:
               data = hdul[1].data
       except:
           with fits.open(reduced_files_loc + star.upper() + '/' + file + 's') as hdul:
               data = hdul[1].data
       wave = data['wave']
       flux = data['flux']
       err = flux*error
       
       cols = []
       cols.append(fits.Column(name='wave', format='D', array=wave))
       cols.append(fits.Column(name='flux', format='D', array=flux))
       cols.append(fits.Column(name='err', format='D', array=err))
       coldefs = fits.ColDefs(cols)
       new_table_hdu = fits.BinTableHDU.from_columns(coldefs)
       
       try:
           with fits.open(reduced_files_loc + star + '/' + file + 's') as hdul:
               hdul[1].data = new_table_hdu.data
               hdul.writeto(reduced_files_loc + star + '/' + file + 's', overwrite=True)
       except:
           with fits.open(reduced_files_loc + star.upper() + '/' + file + 's') as hdul:
               hdul[1].data = new_table_hdu.data
               hdul.writeto(reduced_files_loc + star.upper() + '/' + file + 's', overwrite=True)

       # plt.figure()
       # plt.plot(errors)
       # plt.plot(errors_c)
       # plt.plot(errors_l)
       # plt.plot(errors_r)
       # plt.axhline(error, c='r')
       # plt.show()

include = ['HD', 'moon', 'hd']
# include = ['moon']
   
reduced_stars_folders = os.walk(reduced_files_loc)
reduced_stars_folders = [x[0] for x in reduced_stars_folders]
star_names = [x for x in reduced_stars_folders if any([y in x for y in include])]
star_names = [x[x.rindex('/')+1:] for x in star_names]
star_names = list(set(star_names))

for star in star_names:
   add_err_to_fits(star)
   print(f'Added error to {star}')

# %%
