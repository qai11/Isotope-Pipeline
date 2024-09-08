"""
Title: interp_atmos.py
Author: Quin Aicken Davies
Date: 03/09/24

Description: This script interpolates the atmosphere grids 
in iSpec to the parameters of a given spectra. This will allow moog
to run the correct grid and abundances for finding isotopes.
"""
#%%
import os
import sys
import numpy as np
import logging
import multiprocessing
from multiprocessing import Pool
from matplotlib import pyplot as plt
import pandas as pd
from scipy.stats import chisquare
import scipy.optimize as opt
from copy import deepcopy
import time

#--- iSpec directory -------------------------------------------------------------
if os.path.exists('/home/users/qai11/iSpec_v20201001'):
    "Location of the files on Uni computer"
    ispec_dir = '/home/users/qai11/iSpec_v20201001'
else:
    "location of data on Mac"
    ispec_dir = '/Users/quin/Desktop/2024_Data/iSpec_v20230804'
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec

#%%
"""TEST: using hd_45588 as a test star, Manual input of 102870"""
star_name = 'hd_102870'

# params = pd.read_csv(f'/Users/quin/Desktop/2024_Data/parameters/hd_45588_final_params.txt')
save_path = '/Users/quin/Desktop/2024_Data/'
def interpolate(model, params, code='moog'):
    '''This function will interpolate the atmosphere layers for a given set of parameters'''
    # teff = params['teff'][0]
    teff = 6083.0
    # logg = params['logg'][0]
    logg = 4.1
    # MH = params['MH'][0]
    MH = 0.24
    # alpha = params['alpha'][0]
    alpha = 0.0
    # Load model atmospheres
    modeled_layers_pack = ispec.load_modeled_layers_pack(model)
    atmosphere_layers = ispec.interpolate_atmosphere_layers(modeled_layers_pack, {'teff':teff, 'logg':logg, 'MH':MH, 'alpha':alpha}, code=code)
    atmosphere_layers_file = save_path + star_name + "_atmosphere.moog"
    atmosphere_layers_file = ispec.write_atmosphere(atmosphere_layers, teff, logg, MH, atmosphere_filename=atmosphere_layers_file, code=code)
    

model = ispec_dir + "/input/atmospheres/MARCS.GES/"
interpolate(model, params, code='moog')

# %%
