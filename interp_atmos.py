"""
Title: interp_atmos.py
Author: Quin Aicken Davies
Date: 03/08/24

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
# from pymoog.model import read_marcs_model, save_marcs_model, marcs2moog, read_Kurucz_model, kurucz2moog

#--- iSpec directory -------------------------------------------------------------
if os.path.exists('/home/users/qai11/iSpec_v20201001'):
    "Location of the files on Uni computer"
    ispec_dir = '/home/users/qai11/iSpec_v20201001'
else:
    "location of data on Mac"
    ispec_dir = '/Users/quin/Desktop/2024_Data/iSpec_v20230804'
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec

def interp_atmos(star):
# star = ['hd_45588','hd_100407','hd_102870','hd_128620','hd_128621','hd_11695','hd_146233','hd_156098','hd_157244','hd_160691','moon']

    '''Function for running in parameters_pipeline file'''

    """TEST: using hd_45588 as a test star, Manual input of 102870"""
    # star_name = 'hd_102870'
    # star = ['hd_45588','hd_100407','hd_102870','hd_128620','hd_128621','hd_11695','hd_146233','hd_156098','hd_157244','hd_160691','moon']
    # vmicro = ['1.33','1.40','1.33','1.20','0.99','1.75','1.07','1.00','1.88','1.17','1.06']
    star_num = 0
    # def interpolate(model, star_name, save_path, params=None, code='moog'):
    #         '''This function will interpolate the atmosphere layers for a given set of parameters'''
    #         teff = params['teff'][0]
    #         print(teff)
    #         # teff = 6083.0
    #         logg = params['logg'][0]
    #         # logg = 4.1
    #         MH = params['MH'][0]
    #         # MH = 0.24
    #         alpha = params['alpha'][0]
    #         # alpha = 0.0
    #         # Load model atmospheres
    #         modeled_layers_pack = ispec.load_modeled_layers_pack(model)
    #         atmosphere_layers = ispec.interpolate_atmosphere_layers(modeled_layers_pack, {'teff':teff, 'logg':logg, 'MH':MH, 'alpha':alpha}, code=code)
    #         atmosphere_layers_file = save_path + star_name + "_atmosphere.moog"
    #         atmosphere_layers_file = ispec.write_atmosphere(atmosphere_layers, teff, logg, MH, atmosphere_filename=atmosphere_layers_file, code=code)
    #         print("Atmosphere interpolated successfully!")
            # ,teff,logg,MH,alpha,Vmic,Vmac,Vsini,limb_darkening_coeff,R,CHI-SQUARE,DOF
    for star_name in star:
        try:
            #UNI
            params = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/iteration_parameters/{star_name}_iter_param.txt',header=None)
            save_path = '/home/users/qai11/Documents/Fixed_fits_files/' + star_name + '/'
        except:
            #MAC
            # params = pd.read_csv(f'/Users/quin/quin-masters-code/Masters_stars.csv',header=None, delimiter=',')
            params = pd.read_csv(f'/Users/quin/Desktop/2024_Data/Fixed_fits_files/{star_name}/{star_name}_params.txt',header=None)
            save_path = '/Users/quin/Desktop/2024_Data/Fixed_fits_files/' + star_name + '/'
        
        code='moog'
        
        # model = ispec_dir + "/input/atmospheres/MARCS.GES/"
        model = ispec_dir + "/input/atmospheres/ATLAS9.Castelli/"
        # interpolate(model, star_name, save_path, params=None, code='moog')
        '''This function will interpolate the atmosphere layers for a given set of parameters'''
        teff = float(params[1][-1:].values[0])
        # print(teff)
        # teff = 6093
        logg = float(params[2][-1:].values[0])
        # print(logg)
        # logg = 	4.08
        Mh = float(params[3][-1:].values[0])
        # print(MH)
        # MH = 0.13
        #its alpha enhancement only so check and if less than zero set to 0.0
        alpha = float(params[4][-1:].values[0])
        if alpha < 0:
            alpha = 0.0
        # print(alpha)
        vmic = 1.5
        # Load model atmospheres
        modeled_layers_pack = ispec.load_modeled_layers_pack(model)
        atmosphere_layers = ispec.interpolate_atmosphere_layers(modeled_layers_pack, {'teff':teff, 'logg':logg, 'MH':Mh, 'alpha':alpha}, code=code)
        atmosphere_layers_file = save_path + star_name + "_atmosphere.moog"
        atmosphere_layers_file = ispec.write_atmosphere(atmosphere_layers, teff, logg, Mh, atmosphere_filename=atmosphere_layers_file, code=code)
        print(f"{star_name}Atmosphere interpolated successfully!")

        '''Add the molecules to the atmosphere file by reopeneing the file and adding the molecules'''
        """FROM moog.py lines 134-158 in iSpec"""
        # Append microturbulence, solar abundances and metallicity
        moog_atmosphere = open(save_path + star_name + "_atmosphere.moog", "a")
        # atom_abundances = abundances[np.logical_and(abundances['code'] > 1, abundances['code'] <= 92)] # Don't update hydrogen or helium abundances
        # moog_atmosphere.write("  %.2f\n" % (float(params[5][-1:].values[0])))
        moog_atmosphere.write("  %.2f\n" % (vmic))
        moog_atmosphere.write("NATOMS=   %i %.2f\n" % (0, float(params[3][-1:].values[0])))

        # Molecule list as used by Jorge Melendez (private communication) from iSpec
        moog_atmosphere.write("NMOL      28\n")
        moog_atmosphere.write("  101.0   106.0   107.0   108.0   112.0  126.0\n")
        moog_atmosphere.write("  606.0   607.0   608.0\n")
        moog_atmosphere.write("  707.0   708.0\n")
        moog_atmosphere.write("  808.0   812.0   822.0   823.0   840.0\n")
        moog_atmosphere.write("  10108.0 10820.0 60808.0\n")
        moog_atmosphere.write("  6.1     7.1     8.1   12.1  20.1  22.1  23.1  26.1  40.1\n")
        moog_atmosphere.close()
        
        print(f"{star_name}Text appended successfully!")

# %%
