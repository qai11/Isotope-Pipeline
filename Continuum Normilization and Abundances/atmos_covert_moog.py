"""
Title: atmos_convert_moog.py
Author: Quin Aicken Davies
Date: 03/07/2024

Description: Converts MARCS model atmospheres to MOOG format.
"""
#%%
import os
import subprocess
import numpy as np
import re
import pandas as pd
import pymoog
from scipy.spatial import Delaunay
import pickle
from matplotlib import pyplot as plt
from pymoog.model import read_marcs_model, save_marcs_model, marcs2moog
# from pymoog import linelist
import sys

#--- iSpec directory -------------------------------------------------------------
#ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/"
ispec_dir = '/home/users/qai11/iSpec_v20201001/'
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec


#%%
"""Interpolate the atmosphere using iSpec, this will allow me to use the atmosphere for moog and adjust for 
the correct parameters"""
def interpolate_atmosphere(code="spectrum"):
    #--- Synthesizing spectrum -----------------------------------------------------
    # Parameters
    teff = 4777.0
    logg = 4.44
    MH = 0.00
    alpha = 0.00

    # Selected model amtosphere, linelist and solar abundances
    model = ispec_dir + "/input/atmospheres/MARCS.GES/"

    # Load model atmospheres
    modeled_layers_pack = ispec.load_modeled_layers_pack(model)

    # Validate parameters
    if not ispec.valid_atmosphere_target(modeled_layers_pack, {'teff':teff, 'logg':logg, 'MH':MH, 'alpha':alpha}):
        msg = "The specified effective temperature, gravity (log g) and metallicity [M/H] \
                fall out of theatmospheric models."
        print(msg)

    # Prepare atmosphere model
    atmosphere_layers = ispec.interpolate_atmosphere_layers(modeled_layers_pack, {'teff':teff, 'logg':logg, 'MH':MH, 'alpha':alpha}, code=code)
    atmosphere_layers_file = "example_atmosphere.txt"
    atmosphere_layers_file = ispec.write_atmosphere(atmosphere_layers, teff, logg, MH, atmosphere_filename=atmosphere_layers_file, code=code)

#%%

linelist_file = '/home/users/qai11/quin-masters-code/quinlist.MgH'

ll_data = pymoog.line_data.read_linelist(linelist_file)
#%%
# s = pymoog.synth.synth(6083, 4.1,    0.24,       5100,     5200,          82000,line_list=linelist_file)
# #                      Teff, logg, [Fe/H], wav_start(A), wav_end(A), resolution 
# s.prepare_file()
# s.run_moog()
# s.read_spectra()

#%%
# Plot the synthesized spectra
# plt.figure()
# plt.plot(s.wav, s.flux)
# plt.show()
# %%
'''Testing coversion of MARCS model to MOOG format'''

original_file = 'quin_model.mod'
path = '/home/users/qai11/Documents/Fixed_fits_files/hd_102870/'
save_file = 'quin_model.moog'

molecules = ['606.0', '106.0', '607.0', '608.0', '107.0', '108.0', '112.0', '707.0', '708.0', '808.0', '12.1', '60808.0', '10108.0', '101.0', '6.1', '7.1', '8.1', '822.0', '22.1']
# abun_change = {0: -0.27}

model_dict = read_marcs_model(path + original_file)
#this changes the overall metalicity abundance of the model
model_dict['[M/H]']=0.5 #cant do this you need to interpolate probably using ispec
moog_model_dict = marcs2moog(model_dict, path + save_file, abun_change=None, molecules_include=molecules)

# %%
