"""
Title: print_params.py
Author: Quin Aicken Davies
Date: 08/08/24

Description: Prints the parameters and errors for each star in the final_params.txt and final_errors.txt files.
"""
#%%
import numpy as np
import pandas as pd

# %%
star = ['hd_45588','hd_100407','hd_102870','hd_128620','hd_128621']
iteration = 0
params = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/parameters/final_params.txt')
errors = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/parameters/final_errors.txt')
for i in star:
    print('\n')
    print(star[iteration])
    print(f'teff = {params["teff"][iteration]:.2f} pm {errors["teff"][iteration]:.2f}')
    print(f'logg = {params["logg"][iteration]:.2f} pm {errors["logg"][iteration]:.2f}')
    print(f'MH = {params["MH"][iteration]:.2f} pm {errors["MH"][iteration]:.2f}')
    print(f'vsini = {params["vsini"][iteration]:.2f} pm {errors["vsini"][iteration]:.2f}')
    iteration +=1
# %%
