"""
Title: line_by_line_abunds.py
Author: Heather Sinclair-Wentworth, edits by Quin Aicken Davies, Some extra code for saving from Ethan Bull
Date: 06/10/24

Description: Find the abundances of elements by line-by-line analysis using iSpec. Saves all the synth regions 
for each line for plotting later and checking condition of fit.
"""
#%%
import os
import sys
import numpy as np
import logging
import multiprocessing
from multiprocessing import Pool
import numpy as np
from astropy.io import fits
import pandas as pd
import matplotlib.pyplot as plt
import time
import ast

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
start=time.time()


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import time
import os, sys
import ispec
import glob
from multiprocessing import Pool
from functools import partial
#--- iSpec directory -------------------------------------------------------------
if os.path.exists('/home/users/qai11/iSpec_v20201001'):
    "Location of the files on Uni computer"
    ispec_dir = '/home/users/qai11/iSpec_v20201001'
else:
    "location of data on Mac"
    ispec_dir = '/Users/quin/Desktop/2024_Data/iSpec_v20230804'
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec

time_start = time.time()

def find_abundance(star, elements=None, overwrite_line_regions=False):
    spectrum = ispec.read_spectrum(f'/home/users/qai11/Documents/Fixed_fits_files/{star}/{star}_adjusted.fits')
    parameters = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/parameters/{star}_final_params.txt', sep=',',index_col=0)
    if elements is None:
        raise ValueError('Element must be specified')

    star_continuum_model = ispec.fit_continuum(spectrum, fixed_value=1.0, model="Fixed value")
    code = 'MOOG'
    free_params = []
    
    # Set up the parameters
    teff = parameters['teff'].values[0]
    logg = parameters['logg'].values[0]
    m_h = parameters['MH'].values[0]
    alpha = parameters['alpha'].values[0]
    vmic = parameters['vmic'].values[0]
    vmac = parameters['vmac'].values[0]
    vsini = parameters['vsini'].values[0]
    limb_darkening_coeff = parameters['limb_darkening_coeff'].values[0]
    resolving_power = parameters['R'].values[0]
    vrad = 0
    max_iterations = 20
   
#    wave_base = wave_min
#    wave_top = wave_max 
    for element in elements:    
        # atomic_linelist_file = ispec_dir + "/input/linelists/transitions/Quin_GES_LIST.420_920nm/atomic_lines.tsv"
        atomic_linelist_file = ispec_dir + "/input/linelists/transitions/VALD.300_1100nm/atomic_lines.tsv"
        atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, wave_base=np.min(spectrum['waveobs']), wave_top=np.max(spectrum['waveobs']))
        atomic_linelist = atomic_linelist[atomic_linelist['theoretical_depth'] >= 0.01]  # Select lines that have some minimal contribution in the sun
        
        model = ispec_dir + '/input/atmospheres/MARCS.GES/'
        model_atmosphere = ispec.load_modeled_layers_pack(model)
        
        solar_abundance = ispec_dir + '/input/abundances/Grevesse.2007/stdatom.dat'
        solar_abundances = ispec.read_solar_abundances(solar_abundance)
        
        isotope_file = ispec_dir + "/input/isotopes/SPECTRUM.lst"
        isotopes = ispec.read_isotope_data(isotope_file)
        
        chemical_elements_file = ispec_dir + "/input/abundances/chemical_elements_symbols.dat"
        chemical_elements = ispec.read_chemical_elements(chemical_elements_file)
        
        #    #Build line regions file for specified element from atomic linelist, if it doesn't exist
        line_region_loc = '/home/users/qai11/Documents/quin-masters-code/Linelists/'
        line_regions = line_region_loc + f'{element}_lines.csv'
        
        #Saves the line regions 
        if not os.path.exists(line_regions) or overwrite_line_regions:
            full_line_regions = pd.read_csv(ispec_dir + '/input/regions/47000_GES/limited_but_with_missing_elements_moog_synth_good_for_abundances_all_extended.txt', sep='\t')
            element_line_regions = full_line_regions[full_line_regions['element'].str.split().str[0] == element]
            # print(full_line_regions['element'].str.split().str[0])
            element_line_regions = element_line_regions[['wave_peak', 'wave_base', 'wave_top', 'note']]
            element_line_regions.to_csv(line_regions, sep='\t', index=False)
            
        line_regions = ispec.read_line_regions(line_regions)
        
        output_dir = f'/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/{star}/'
        ispec.mkdir_p(output_dir)
        
        i = 0
        for a, line in enumerate(line_regions):
            element_name = '_'.join(line['note'].split())
            filename = f'{output_dir}'+ 'line_by_line' + f'_{element}'
            
            try:
                free_abundances = ispec.create_free_abundances_structure([element], chemical_elements, solar_abundances)
            except:
                raise ValueError(f'Element {element} not found in chemical_elements')
            
            free_abundances['Abund'] += m_h #Scale to metallicity
            
            linelist_free_loggf = None
            
            individual_line_regions = line_regions[a:a+1] # Keep recarray structure

            # Segment
            segments = ispec.create_segments_around_lines(individual_line_regions, margin=0.25)
            wfilter = ispec.create_wavelength_filter(spectrum, regions=segments) # Only use the segment

            if len(spectrum[wfilter]) == 0 or np.any(spectrum['flux'][wfilter] <= 0):
                continue
            
            try:
                obs_spec, modeled_synth_spectrum, params, errors, abundances_found, loggf_found, status, stats_linemasks = \
                        ispec.model_spectrum(spectrum[wfilter], star_continuum_model, \
                        model_atmosphere, atomic_linelist, isotopes, solar_abundances, free_abundances, linelist_free_loggf, teff, \
                        logg, m_h, alpha, vmic, vmac, vsini, \
                        limb_darkening_coeff, resolving_power, vrad, free_params, segments=segments, \
                        linemasks=individual_line_regions, \
                        enhance_abundances=True, \
                        use_errors = True, \
                        vmic_from_empirical_relation = False, \
                        vmac_from_empirical_relation = False, \
                        max_iterations=max_iterations, \
                        tmp_dir = None, \
                        code=code)
            except:
                continue
            synth = pd.DataFrame(modeled_synth_spectrum[['waveobs', 'flux']], columns=['waveobs', 'flux'])
            
            abundances_found = pd.DataFrame(abundances_found, index=[i])
            abundances_found['element'] = element_name
            abundances_found = abundances_found[['element', 'code', 'Abund', 'A(X)', '[X/H]', '[X/Fe]', 'eAbund', 'eA(X)', 'e[X/H]', 'e[X/Fe]']]
            abundances_found = pd.concat([abundances_found, pd.DataFrame({'wave_peak': line['wave_peak']}, index=[i])], axis=1)
            # print('concatenated')
            
            #Define folder for saving. If it doesn't exist, create it
            synth_folder = output_dir + 'synth_regions_' + element + '/'
            if not os.path.exists(synth_folder):
                os.mkdir(synth_folder)
            filename_synth = f'{synth_folder}'+ f'synth_{element}' + '_' + str(line['wave_peak']) + '.txt'
            synth.to_csv(filename_synth, sep=' ', index=False, mode='w', header=True)
            
            
            if i == 0:
                abundances_found.to_csv(filename + '.txt', sep=' ', index=True, mode='w', header=True)
            else:
                abundances_found.to_csv(filename + '.txt', sep=' ', index=True, mode='a', header=False)
            i += 1


# def mad(data):
#    return np.median(np.abs(data - np.median(data)))

# new_abunds = True
# elements = ['Ca','Fe','Mg','Si','Ti','C','O']
# if type(elements) == list:
#    for element in elements:
#        print(f'Finding abundances for {element}')
#        if new_abunds:
#            try:
#                pool = Pool(os.cpu_count() - 1)
#                print(f'Opened pools on {os.cpu_count()-1} threads.')
#                if len(star_names) > 0:
#                    fit = pool.map(partial(find_abundance, element = element), star_names)
#            finally:
#                pool.close()
#                pool.join()
#        print(f'Finished finding abundances for {element}')
#        for star in star_names:
#            abundances = pd.read_csv(f'/home/users/ebu36/data/ASTR690/SAUCE/echelle_fits/SAUCE_stars/abunds/{star}/line_by_line_{element}.txt', sep=' ')
#            # Errors from MAD
#            median_abund_values = pd.DataFrame({'Abund': np.median(abundances['Abund'].values), 'A(X)': np.median(abundances['A(X)'].values), '[X/H]': np.median(abundances['[X/H]'].values), '[X/Fe]': np.median(abundances['[X/Fe]'].values),\
#                                                'eAbund': mad(abundances['eAbund']), 'eA(X)': mad(abundances['eA(X)']), 'e[X/H]': mad(abundances['[X/H]']), 'e[X/Fe]': mad(abundances['[X/Fe]'])}, index=[0])
#            median_abund_values.to_csv(f'/home/users/ebu36/data/ASTR690/SAUCE/echelle_fits/SAUCE_stars/abunds/{star}/median_abund_{element}.txt', sep=' ', index=False)
# else:
#    raise ValueError('Elements must be a list of strings')




end = time.time()
print(f'It took {end-start} seconds to run.')


# %%
