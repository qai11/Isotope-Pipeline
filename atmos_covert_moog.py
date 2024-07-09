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
#%%

linelist_file = '/home/users/qai11/quin-masters-code/quinlist.MgH'

ll_data = pymoog.line_data.read_linelist(linelist_file)
#%%
s = pymoog.synth.synth(6083, 4.1,    0.24,       5100,     5200,          82000,line_list=linelist_file)
#                      Teff, logg, [Fe/H], wav_start(A), wav_end(A), resolution 
s.prepare_file()
s.run_moog()
s.read_spectra()

#%%
# Plot the synthesized spectra
plt.figure()
plt.plot(s.wav, s.flux)
plt.show()
# %%
'''Testing coversion of MARCS model to MOOG format'''

original_file = 'quin_model.mod'
path = '/home/users/qai11/Documents/Fixed_fits_files/hd_102870/'
save_file = 'quin_model.moog'

molecules = ['101.0', '106.0', '107.0', '108.0', '606.0', '607.0', '608.0', '707.0', '708.0', '808.0', '10108.0', '60808.0']
# abun_change = {0: -0.27}

model_dict = read_marcs_model(path + original_file)
#this changes the overall metalicity abundance of the model
model_dict['[M/H]']=-0.27
moog_model_dict = marcs2moog(model_dict, path + save_file, abun_change=None, molecules_include=molecules)

# %%
