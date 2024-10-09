"""
Title: parameters_pipeline.py
Author: Quin Aicken Davies
Date: 06/10/24

Description: This runs the code that does rv correction, merging, adjusts the continuum 
placement,generates the parameters then generates the metallities of all listed stars. 
This will allow direct input to interpolate atmospheres for each star and allow for 
isotopes to be estimated without having to manually input the parameters. 

Appropriate iSpec edits will be listed on the Isotopes wiki page in the near future.
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

#%%
#This must be run in the correct order or it wont do what is needed
'''Note: The linelist must contain isotopes and the relevant molecules
for the abundances and isotopes to be found and be consistent.
Most things should be on parallel pool so should work quickly if all things go well.'''
#First Deianira must be run to completion first
#Then rv_combine.py
#Then continuum_adjust.py
#Then find_params.py
#Then find_abund.py
#Then interp_atmos.py
#Then the data can be run through to find the isotopes

'''At the moment I have the stars all hard coded in,
in order to change this we will have to remove that and add in a loop
once a star list is loaded in, an easy change which I can do if requested but
not needed at this moment.'''

'''Files are saved as they go so don't need to  be output as named variables'''

#%%
# Import the scripts to be run
from rv_combine import rv_combine
from continuum_adjust import continuum_adjust
from find_params import find_params
# from find_abund import find_abund
from line_by_line_abunds import run_abunds
from interp_atmos import interp_atmos


star = ['hd_45588','hd_100407','hd_102870','hd_128620','hd_128621','hd_11695','hd_146233','hd_156098','hd_157244','hd_160691','moon']
element = ["Mg", "Si", "Ca", "Ti", "Sc","V","Cr","Mn","Co", "Ni", "Y", "Ba", "La", "Nd", "Eu", "Sr", "Zr","Rb"]
#%%
start = time.time()
'''finish setting up parallel pools for each thing'''
#Run the rv correction and Merging
rv_combine(star)
#Run the continuum adjustment
try:
    pool = Pool(os.cpu_count()-1)
    pool.map(continuum_adjust, star)
finally:
    pool.close()
    pool.join()
#Run the parameter finding
try:
    pool = Pool(os.cpu_count()-1)
    pool.map(find_params, star)
finally:
    pool.close()
    pool.join()
#Run the abundance finding
try:
    pool = Pool(os.cpu_count()-1)
    pool.map(run_abunds, star, element)
finally:
    pool.close()
    pool.join()
#Run the interpolation of the atmospheres
interp_atmos(star)

end = time.time()

print(f'Time taken: {(end - start)/3600} Hrs')
print('Pipeline Complete done!')


# %%
