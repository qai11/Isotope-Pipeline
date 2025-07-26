"""
Title: 1_iter_continuum_adjust.py
Author: Quin Aicken Davies
Date: 10/07/2024

Description: Uses the iSpec and my own code to adjust the continuum level for a stellar spectrum.
Only does it once with hope of making the isotope fit better.
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

def synthesize_spectrum(teff,logg,MH,vsini,wave_base, wave_top, code="moog",wave_step=0.001):
    #--- Synthesizing spectrum -----------------------------------------------------
    # Parameters
    # teff = 6083
    # logg = 4.1
    # MH = 0.27
    alpha = ispec.determine_abundance_enchancements(MH)
    microturbulence_vel = ispec.estimate_vmic(teff, logg, MH) # 1.07
    macroturbulence = ispec.estimate_vmac(teff, logg, MH) # 4.21
    # vsini = 2.0
    limb_darkening_coeff = 0.6
    resolution = 82000

    # Wavelengths to synthesis
    #regions = ispec.read_segment_regions(ispec_dir + "/input/regions/fe_lines_segments.txt")
    regions = None
    # wave_base = wave_min
    # wave_top = wave_max


    # Selected model amtosphere, linelist and solar abundances
    model = ispec_dir + "/input/atmospheres/MARCS.GES/"
    
    atomic_linelist_file = ispec_dir + "/input/linelists/transitions/Quin_GES_LIST.420_920nm/atomic_lines.tsv"

    isotope_file = ispec_dir + "/input/isotopes/SPECTRUM.lst"

    # Load chemical information and linelist
    atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, wave_base=wave_base, wave_top=wave_top)
    atomic_linelist = atomic_linelist[atomic_linelist['theoretical_depth'] >= 0.01] # Select lines that have some minimal contribution in the sun

    isotopes = ispec.read_isotope_data(isotope_file)

    if "ATLAS" in model:
        solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.1998/stdatom.dat"
    else:
        # MARCS
        solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.2007/stdatom.dat"

    # Load model atmospheres
    modeled_layers_pack = ispec.load_modeled_layers_pack(model)
    # Load SPECTRUM abundances
    solar_abundances = ispec.read_solar_abundances(solar_abundances_file)

    ## Custom fixed abundances
    #fixed_abundances = ispec.create_free_abundances_structure(["C", "N", "O"], chemical_elements, solar_abundances)
    #fixed_abundances['Abund'] = [-3.49, -3.71, -3.54] # Abundances in SPECTRUM scale (i.e., x - 12.0 - 0.036) and in the same order ["C", "N", "O"]
    ## No fixed abundances
    fixed_abundances = None

    # Validate parameters
    if not ispec.valid_atmosphere_target(modeled_layers_pack, {'teff':teff, 'logg':logg, 'MH':MH, 'alpha':alpha}):
        msg = "The specified effective temperature, gravity (log g) and metallicity [M/H] \
                fall out of theatmospheric models."
        print(msg)

    # Prepare atmosphere model
    atmosphere_layers = ispec.interpolate_atmosphere_layers(modeled_layers_pack, {'teff':teff, 'logg':logg, 'MH':MH, 'alpha':alpha}, code=code)

    # Synthesis
    synth_spectrum = ispec.create_spectrum_structure(np.arange(wave_base, wave_top, wave_step))
    synth_spectrum['flux'] = ispec.generate_spectrum(synth_spectrum['waveobs'], \
            atmosphere_layers, teff, logg, MH, alpha, atomic_linelist, isotopes, solar_abundances, \
            fixed_abundances, microturbulence_vel = microturbulence_vel, \
            macroturbulence=macroturbulence, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, \
            R=resolution, regions=regions, verbose=1,
            code=code)
    # ##--- Save spectrum ------------------------------------------------------------
    logging.info("Saving spectrum...")
    synth_filename = "hd_102870_synth_%s.fits" % (code)
    ispec.write_spectrum(synth_spectrum, synth_filename)
    return synth_spectrum['waveobs'], synth_spectrum['flux'], synth_spectrum

def normalize_whole_spectrum_with_template(star_spectrum, synthetic_spectrum):
    """
    Use a template to normalize the whole spectrum
    """
    # star_spectrum = ispec.read_spectrum("/home/users/qai11/Documents/Fixed_fits_files/hd_102870/J0354027.txt")
    # synth_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/templates/Synth.Sun.300_1100nm/template.txt.gz")

    #--- Continuum fit -------------------------------------------------------------
    model = "Template"
    nknots = None # Automatic: 1 spline every 5 nm (in this case, used to apply a gaussian filter)
    from_resolution = 82000
    median_wave_range=1

    #strong_lines = ispec.read_line_regions(ispec_dir + "/input/regions/strong_lines/absorption_lines.txt")
    # strong_lines = ispec.read_line_regions(ispec_dir + "/input/regions/relevant/relevant_line_masks.txt")
    strong_lines = None
    star_continuum_model = ispec.fit_continuum(star_spectrum, from_resolution=from_resolution, \
                                ignore=strong_lines, \
                                nknots=nknots, \
                                median_wave_range=median_wave_range, \
                                model=model, \
                                template=synthetic_spectrum)

    #--- Continuum normalization ---------------------------------------------------
    logging.info("Continuum normalization...")
    normalized_star_spectrum = ispec.normalize_spectrum(star_spectrum, star_continuum_model, consider_continuum_errors=False)
    # Use a fixed value because the spectrum is already normalized
    # star_continuum_model = ispec.fit_continuum(star_spectrum, fixed_value=1.0, model="Fixed value")
    return normalized_star_spectrum, star_continuum_model

#%%

"""Read in the star spectrum and the synthetic spectrum for fitting later"""
# star = ['hd_45588','hd_100407','hd_102870','hd_128620','hd_128621']
# star = ['hd_45588','hd_100407','hd_102870','hd_128620','hd_128621','hd_11695','hd_146233','hd_156098','hd_157244','hd_160691','moon']
# star = ['hd_11695','hd_146233','hd_156098','hd_157244','hd_160691','moon']

star = ['hd_128620']
# star = ['hd_157244']
# star = ['hd_102870']

# def continuum_adjust(star_name):
for star_name in star:
    try:
        #Uni computer
        star_spectrum = ispec.read_spectrum(f"/home/users/qai11/Documents/Fixed_fits_files/{star_name}/rv_corrected/median_spectrum_{star_name}.txt")
    except:
        #Mac
        star_spectrum = ispec.read_spectrum(f"/Users/quin/Desktop/2024_Data/Fixed_fits_files/{star_name}/rv_corrected/median_spectrum_{star_name}.txt")

    #open up the file information and copy the spectrum for correction
    star_wave = star_spectrum['waveobs']
    star_flux = star_spectrum['flux']
    star_flux_err = star_spectrum['err']
    wave_step = 0.00017#np.round(star_wave[1]-star_wave[0],13)
    spectrum_copy = deepcopy(star_spectrum)
    
    #set the parameters for the synthetic spectrum
    if star_name == 'hd_128620':
        modeled_synth_wave, modeled_synth_flux, modeled_synth_spectrum = synthesize_spectrum(5792,4.31,0.26,1.9, wave_base=480, wave_top=680,wave_step=wave_step)
    elif star_name == 'hd_157244':
        modeled_synth_wave, modeled_synth_flux, modeled_synth_spectrum = synthesize_spectrum(4197,1.05,-0.05,5.4, wave_base=480, wave_top=680,wave_step=wave_step)
    elif star_name == 'hd_102870':
        modeled_synth_wave, modeled_synth_flux, modeled_synth_spectrum = synthesize_spectrum(6083,4.1,0.27,2.0, wave_base=480, wave_top=680,wave_step=wave_step)
    
    normalized_spectrum, star_continuum_model = normalize_whole_spectrum_with_template(spectrum_copy, modeled_synth_spectrum)
    
    #Load file path and save
    try:
        #Uni computer
        star_filename = f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/'+ f"{star_name}_adjusted.fits" 
        ispec.write_spectrum(normalized_spectrum, star_filename)
    except:
        #Mac
        star_filename = f'/Users/quin/Desktop/2024_Data/Fixed_fits_files/{star_name}/'+ f"{star_name}_adjusted.fits"
        ispec.write_spectrum(normalized_spectrum, star_filename)
    print(f'{star_name} file saved')

#plot for visual inspection
plt.figure()
plt.plot(star_wave,star_flux)
plt.plot(modeled_synth_wave,modeled_synth_flux)
plt.plot(normalized_spectrum['waveobs'],normalized_spectrum['flux'])
plt.xlim(513.38, 513.55)
plt.ylim(0.8,1.05)
plt.legend(['Observed','Synthetic','Normalized'])
plt.show()
# %%

#Resave the original spectrum as a fits file
star = ['hd_128620','hd_157244','hd_102870']
for star_name in star:
    try:
        #Uni computer
        star_spectrum = ispec.read_spectrum(f"/home/users/qai11/Documents/Fixed_fits_files/{star_name}/rv_corrected/median_spectrum_{star_name}.txt")
    except:
        #Mac
        star_spectrum = ispec.read_spectrum(f"/Users/quin/Desktop/2024_Data/Fixed_fits_files/{star_name}/rv_corrected/median_spectrum_{star_name}.txt")

    try:
        #Uni computer
        star_filename = f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/'+ f"{star_name}_adjusted.fits" 
        ispec.write_spectrum(star_spectrum, star_filename)
        print(f'{star_name} file saved')
    except:
        #Mac
        star_filename = f'/Users/quin/Desktop/2024_Data/Fixed_fits_files/{star_name}/'+ f"{star_name}_adjusted.fits"
        ispec.write_spectrum(star_spectrum, star_filename)
        print(f'{star_name} file saved')
# %%
