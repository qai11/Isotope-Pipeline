"""
Title: find_abund.py
Author: Quin Aicken Davies
Date: 19/08/24

Description: This will take the values from find_params.py and set the values
for finding abundance values for elements in a list.
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
#ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/"
ispec_dir = '/home/users/qai11/iSpec_v20201001'
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec

def find_abund(star, element):
    '''Function for running in parameters_pipeline file'''
    #%%
    def determine_abundances_using_synth_spectra(star_spectrum, teff, logg, MH, Vmic, Vmac, alpha, vsini, max_iterations, element_name, wave_base=480, wave_top=680, resolution=82000, code="moog"):
        # star_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
        star_spectrum = star_spectrum
        wave_base=wave_base
        wave_top=wave_top
        
    #--- Continuum fit -------------------------------------------------------------
        model = "Fixed value" # 
        degree = 1
        nknots = None # Automatic: 1 spline every 5 nm
        from_resolution = resolution

        # Strategy: Filter first median values and secondly MAXIMUMs in order to find the continuum
        order='median+max'
        median_wave_range=0.05
        max_wave_range=1.0

        star_continuum_model = ispec.fit_continuum(star_spectrum, from_resolution=from_resolution, \
                                    nknots=nknots, degree=degree, \
                                    median_wave_range=median_wave_range, \
                                    max_wave_range=max_wave_range, \
                                    model=model, order=order, \
                                    automatic_strong_line_detection=True, \
                                    strong_line_probability=0.5, \
                                    use_errors_for_fitting=True,fixed_value=1.0)
        #--- Normalize -------------------------------------------------------------
        normalized_star_spectrum = ispec.normalize_spectrum(star_spectrum, star_continuum_model, consider_continuum_errors=False)
        # Use a fixed value because the spectrum is already normalized
        star_continuum_model = ispec.fit_continuum(star_spectrum, fixed_value=1.0, model="Fixed value")
        #--- Model spectra ----------------------------------------------------------
        # Parameters
        initial_teff = teff
        initial_logg = logg
        initial_MH = MH
        initial_alpha = alpha #ispec.determine_abundance_enchancements(initial_MH)
        initial_vmic = Vmic #ispec.estimate_vmic(initial_teff, initial_logg, initial_MH)
        initial_vmac = Vmac #ispec.estimate_vmac(initial_teff, initial_logg, initial_MH)
        initial_vsini = vsini
        initial_limb_darkening_coeff = 0.6
        initial_R = resolution
        initial_vrad = 0
        max_iterations = max_iterations

        # Selected model amtosphere, linelist and solar abundances
        model = ispec_dir + "/input/atmospheres/MARCS.GES/"
        atomic_linelist_file = ispec_dir + "/input/linelists/transitions/GBS_LIST.480_800nm/atomic_lines.tsv"
        
        if "ATLAS" in model:
            solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.1998/stdatom.dat"
        else:
            # MARCS
            solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.2007/stdatom.dat"

        isotope_file = ispec_dir + "/input/isotopes/SPECTRUM.lst"

        # Load chemical information and linelist
        atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, wave_base=np.min(star_spectrum['waveobs']), wave_top=np.max(star_spectrum['waveobs']))
        atomic_linelist = atomic_linelist[atomic_linelist['theoretical_depth'] >= 0.01] # Select lines that have some minimal contribution in the sun

        isotopes = ispec.read_isotope_data(isotope_file)

        # Load model atmospheres
        modeled_layers_pack = ispec.load_modeled_layers_pack(model)

        # Load SPECTRUM abundances
        solar_abundances = ispec.read_solar_abundances(solar_abundances_file)

        # Free parameters
        #free_params = ["teff", "logg", "MH", "vmic", "vmac", "vsini", "R", "vrad", "limb_darkening_coeff"]
        # free_params = ["vrad"]
        free_params = []

        # Free individual element abundance (WARNING: it should be coherent with the selected line regions!)
        chemical_elements_file = ispec_dir + "/input/abundances/chemical_elements_symbols.dat"
        chemical_elements = ispec.read_chemical_elements(chemical_elements_file)

        element_name = element_name
        #, "Si", "Ca", "Ti", "Sc","V","Cr","Mn","Co", "Ni", "Y", "Ba", "La", "Nd", "Eu"
        free_abundances = ispec.create_free_abundances_structure([element_name], chemical_elements, solar_abundances)
        free_abundances['Abund'] += initial_MH # Scale to metallicity

        linelist_free_loggf = None

        # Line regions
        line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_GES/{}_synth_good_for_params_all.txt".format(code))
        #line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_GES/{}_synth_good_for_params_all_extended.txt".format(code))
        #line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_VALD/{}_synth_good_for_params_all.txt".format(code))
        #line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_VALD/{}_synth_good_for_params_all_extended.txt".format(code))
        # Select only the lines to get abundances from
        # line_regions = line_regions[np.logical_or(line_regions['note'] == element_name+' 1', line_regions['note'] == element_name+' 2')]
        # line_regions = ispec.adjust_linemasks(normalized_star_spectrum, line_regions, max_margin=0.5)
        
        # line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/fe_lines.txt".format(code))
        line_regions = None
        
        # Read segments if we have them or...
        #segments = ispec.read_segment_regions(ispec_dir + "/input/regions/fe_lines_segments.txt")
        # ... or we can create the segments on the fly:
        # segments = ispec.create_segments_around_lines(line_regions, margin=0.25)
        segments = ispec.read_segment_regions(ispec_dir + "/input/regions/Quin_segments.txt")


        obs_spec, modeled_synth_spectrum, params, errors, abundances_found, loggf_found, status, stats_linemasks = \
                ispec.model_spectrum(normalized_star_spectrum, star_continuum_model, \
                modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, free_abundances, linelist_free_loggf, initial_teff, \
                initial_logg, initial_MH, initial_alpha, initial_vmic, initial_vmac, initial_vsini, \
                initial_limb_darkening_coeff, initial_R, initial_vrad, free_params, segments=segments, \
                linemasks=line_regions, \
                enhance_abundances=True, \
                use_errors = True, \
                vmic_from_empirical_relation = False, \
                vmac_from_empirical_relation = False, \
                max_iterations=max_iterations, \
                tmp_dir = None, \
                code=code)

        # Selected model amtosphere, linelist and solar abundances
        model = ispec_dir + "/input/atmospheres/MARCS.GES/"
        
        return abundances_found

    #%%

    start = time.time()
    # star = ['hd_45588','hd_100407','hd_102870','hd_128620','hd_128621','hd_11695','hd_146233','hd_156098','hd_157244','hd_160691','moon']
    # star = ['hd_102870']
    # element = ["Mg", "Si", "Ca", "Ti", "Sc","V","Cr","Mn","Co", "Ni", "Y", "Ba", "La", "Nd", "Eu", "Sr", "Zr","Rb"]
    # element = ["Mg", "Si"]
    star_number = 0
    # errors_df = pd.DataFrame()
    # params_df = pd.DataFrame()
    for star_name in star:
        abund_params_df=pd.DataFrame()
        element_number = 0
        for element_name in element:
            print(f"element number  = {element_number}")
            star_spectrum = ispec.read_spectrum(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/{star_name}_adjusted.fits')
            parameters = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/iteration_parameters/{star_name}_iter_param.txt')
            teff = parameters['teff'][len(parameters)-1]
            logg = parameters['logg'][len(parameters)-1]
            MH = parameters['MH'][len(parameters)-1]
            Vmic = parameters['Vmic'][len(parameters)-1]
            Vmac = parameters['Vmac'][len(parameters)-1]
            alpha = parameters['alpha'][len(parameters)-1]
            vsini = parameters['Vsini'][len(parameters)-1]
            max_iterations = 20
            abundances_found = determine_abundances_using_synth_spectra(star_spectrum, teff, logg, MH, Vmic, Vmac, alpha, vsini, max_iterations, element_name=element_name, wave_base=480, wave_top=680, resolution=82000, code="moog")
            #Save the parameters in a dataframe
            #Add the abundances to a pandas dataframe
            abund_params_df = pd.concat([abund_params_df,pd.DataFrame(abundances_found, index=[f'{element_number}'])])
            # os.rename('/home/users/qai11/Documents/Fixed_fits_files/iteration_parameters/iter_param.txt',f'/home/users/qai11/Documents/Fixed_fits_files/iteration_parameters/{star_name}_iter_param.txt')
            element_number+=1
        star_number+=1

            
        #Save all the results to one file
        abund_params_df.to_csv(f'/home/users/qai11/Documents/Fixed_fits_files/parameters/{star_name}_final_indiv_abund.txt')
        print('Files saved')
        
    end = time.time()

    print(f'Time taken in sec: {end - start}')
    # %%
