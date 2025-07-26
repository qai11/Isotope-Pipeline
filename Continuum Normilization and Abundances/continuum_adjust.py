"""
Title: continuum_adjust.py
Author: Quin Aicken Davies
Date: 09/07/2024

Description: Uses the iSpec and my own code to adjust the continuum level for a stellar spectrum.
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


def continuum_adjust(star_name):
    '''Function for running in parameters_pipeline file'''
    start = time.time()
    #%%
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
    #%%
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
    def determine_astrophysical_parameters_using_synth_spectra(star_spectrum, teff, logg,MH, vsini, max_iterations,loop_iteration,wave_base=480, wave_top=680, resolution=82000, code="moog",wave_step=0.001):
        # star_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
        star_spectrum = star_spectrum
        wave_base=wave_base
        wave_top=wave_top
        # #--- Radial Velocity determination with template -------------------------------
        # logging.info("Radial velocity determination with template...")
        # # - Read synthetic template
        # #template = ispec.read_spectrum(ispec_dir + "/input/spectra/templates/Atlas.Arcturus.372_926nm/template.txt.gz")
        # #template = ispec.read_spectrum(ispec_dir + "/input/spectra/templates/Atlas.Sun.372_926nm/template.txt.gz")
        # template = ispec.read_spectrum(ispec_dir + "/input/spectra/templates/NARVAL.Sun.370_1048nm/template.txt.gz")
        # #template = ispec.read_spectrum(ispec_dir + "/input/spectra/templates/Synth.Sun.300_1100nm/template.txt.gz")

        # models, ccf = ispec.cross_correlate_with_template(star_spectrum, template, \
        #                         lower_velocity_limit=-200, upper_velocity_limit=200, \
        #                         velocity_step=1.0, fourier=False)

        # # Number of models represent the number of components
        # components = len(models)
        # # First component:
        # rv = np.round(models[0].mu(), 2) # km/s
        # rv_err = np.round(models[0].emu(), 2) # km/s
        #--- Radial Velocity correction ------------------------------------------------
        # logging.info("Radial velocity correction... %.2f +/- %.2f" % (rv, rv_err))
        # star_spectrum = ispec.correct_velocity(star_spectrum, rv)
        #--- Resolution degradation ----------------------------------------------------
        # NOTE: The line selection was built based on a solar spectrum with R ~ 47,000 and GES/VALD atomic linelist.
        # from_resolution = 80000
        # to_resolution = 47000
        # star_spectrum = ispec.convolve_spectrum(star_spectrum, to_resolution, from_resolution)
        #--- Continuum fit -------------------------------------------------------------
        model = "Fixed value" # "Polynomy"
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
        initial_alpha = ispec.determine_abundance_enchancements(initial_MH)
        initial_vmic = ispec.estimate_vmic(initial_teff, initial_logg, initial_MH)
        initial_vmac = ispec.estimate_vmac(initial_teff, initial_logg, initial_MH)
        initial_vsini = vsini
        initial_limb_darkening_coeff = 0.6
        initial_R = resolution
        initial_vrad = 0
        max_iterations = max_iterations

        # Selected model amtosphere, linelist and solar abundances
        model = ispec_dir + "/input/atmospheres/MARCS.GES/"

        atomic_linelist_file = ispec_dir + "/input/linelists/transitions/VALD.300_1100nm/atomic_lines.tsv"
        #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/VALD.1100_2400nm/atomic_lines.tsv"
        # atomic_linelist_file = ispec_dir + "/input/linelists/transitions/GESv6_atom_hfs_iso.420_920nm/atomic_lines.tsv"
        #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/GESv6_atom_nohfs_noiso.420_920nm/atomic_lines.tsv"
        # atomic_linelist_file = ispec_dir + "/input/linelists/transitions/Quin_GES_LIST.420_920nm/atomic_lines.tsv"
        

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
        free_params = ["teff", "logg", "MH", "vmic", "vmac", "vsini","alpha"]

        # Free individual element abundance
        free_abundances = None
        linelist_free_loggf = None

        # Line regions
        # line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_GES/{}_synth_good_for_params_all.txt".format(code))
        #line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_GES/{}_synth_good_for_params_all_extended.txt".format(code))
        #line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_VALD/{}_synth_good_for_params_all.txt".format(code))
        #line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_VALD/{}_synth_good_for_params_all_extended.txt".format(code))
        ## Select only some lines to speed up the execution (in a real analysis it is better not to do this)
        # line_regions = line_regions[np.logical_or(line_regions['note'] == 'Ti 1', line_regions['note'] == 'Ti 2')]
        # line_regions = ispec.adjust_linemasks(normalized_star_spectrum, line_regions, max_margin=0.5)
        # Read segments if we have them or...
        #segments = ispec.read_segment_regions(ispec_dir + "/input/regions/fe_lines_segments.txt")
        # ... or we can create the segments on the fly:
        line_regions = None
        # segments = ispec.create_segments_around_lines(line_regions, margin=0.25)
        # segments = None
        # line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/GBS_optical_list.tsv")
        # segments = ispec.create_segments_around_lines(line_regions, margin=0.25)
        segments = ispec.read_segment_regions(ispec_dir + "/input/regions/Quin_segments.txt")
        
        ### Add also regions from the wings of strong lines:
        ## H beta
        #hbeta_lines = ispec.read_line_regions(ispec_dir + "input/regions/wings_Hbeta.txt")
        #hbeta_segments = ispec.read_segment_regions(ispec_dir + "input/regions/wings_Hbeta_segments.txt")
        #line_regions = np.hstack((line_regions, hbeta_lines))
        #segments = np.hstack((segments, hbeta_segments))
        ## H alpha
        #halpha_lines = ispec.read_line_regions(ispec_dir + "input/regions/wings_Halpha.txt")
        #halpha_segments = ispec.read_segment_regions(ispec_dir + "input/regions/wings_Halpha_segments.txt")
        #line_regions = np.hstack((line_regions, halpha_lines))
        #segments = np.hstack((segments, halpha_segments))
        ## Magnesium triplet
        #mgtriplet_lines = ispec.read_line_regions(ispec_dir + "input/regions/wings_MgTriplet.txt")
        #mgtriplet_segments = ispec.read_segment_regions(ispec_dir + "input/regions/wings_MgTriplet_segments.txt")
        #line_regions = np.hstack((line_regions, mgtriplet_lines))
        #segments = np.hstack((segments, mgtriplet_segments))

        obs_spec, modeled_synth_spectrum, params, errors, abundances_found, loggf_found, status, stats_linemasks = \
                ispec.model_spectrum(normalized_star_spectrum, star_continuum_model, \
                modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, free_abundances, linelist_free_loggf, initial_teff, \
                initial_logg, initial_MH, initial_alpha, initial_vmic, initial_vmac, initial_vsini, \
                initial_limb_darkening_coeff, initial_R, initial_vrad, free_params, segments=segments, \
                linemasks=line_regions, \
                enhance_abundances=True, \
                use_errors = True, \
                vmic_from_empirical_relation = False, \
                vmac_from_empirical_relation = True, \
                max_iterations=max_iterations, \
                tmp_dir = None, \
                code=code)
                
        ##--- Save results -------------------------------------------------------------
        # logging.info("Saving results...")
        # dump_file = "example_results_synth_%s.dump" % (code)
        # logging.info("Saving results...")
        # ispec.save_results(dump_file, (params, errors, abundances_found, loggf_found, status, stats_linemasks))
        # # If we need to restore the results from another script:
        # params, errors, abundances_found, loggf_found, status, stats_linemasks = ispec.restore_results(dump_file)

        # logging.info("Saving synthetic spectrum...")
        # synth_filename = "example_modeled_synth_%s.fits" % (code)
        # ispec.write_spectrum(modeled_synth_spectrum, synth_filename)

        return obs_spec, modeled_synth_spectrum, params, errors, abundances_found, loggf_found 

    #%%
    """Read in the star spectrum and the synthetic spectrum for fitting later"""
    # star = ['hd_45588','hd_100407','hd_102870','hd_128620','hd_128621']
    # star = ['hd_45588','hd_100407','hd_102870','hd_128620','hd_128621','hd_11695','hd_146233','hd_156098','hd_157244','hd_160691','moon']
    # star = ['hd_11695','hd_146233','hd_156098','hd_157244','hd_160691','moon']

    # def continuum_adjust(star_name):
    try:
        #Uni computer
        star_spectrum = ispec.read_spectrum(f"/home/users/qai11/Documents/Fixed_fits_files/{star_name}/rv_corrected/median_spectrum_{star_name}.txt")
    except:
        #Mac
        star_spectrum = ispec.read_spectrum(f"/Users/quin/Desktop/2024_Data/Fixed_fits_files/{star_name}/rv_corrected/median_spectrum_{star_name}.txt")
    star_wave = star_spectrum['waveobs']
    star_flux = star_spectrum['flux']
    star_flux_err = star_spectrum['err']
    wave_step = 0.00017#np.round(star_wave[1]-star_wave[0],13)
   
    
    """New code for normalizing the spectrum"""
    iteration_number = 0
    chi2_df = pd.DataFrame(columns=['iteration_number','chi2'])
    loop_spectrum = deepcopy(star_spectrum)
    errors_df = pd.DataFrame()
    params_df = pd.DataFrame()
    first_loop=True
    while True:
        print(star_name[iteration_number])
        # Perform the normalization iteratively until the chi-squared value is minimized
        loop_spectrum, modeled_synth_spectrum, params, errors, abundances_found, loggf_found = determine_astrophysical_parameters_using_synth_spectra(loop_spectrum, teff=5792,logg=4.31,MH=0.26,vsini=2.0, max_iterations=5, loop_iteration=iteration_number,wave_base=480, wave_top=680, resolution=82000, code="moog",wave_step=0.001)
        # loop_spectrum, modeled_synth_spectrum, params, errors, abundances_found, loggf_found = determine_astrophysical_parameters_using_synth_spectra(loop_spectrum, teff=6080,logg=4.1,MH=0.24,vsini=2.0, max_iterations=1, loop_iteration=iteration_number,wave_base=480, wave_top=680, resolution=82000, code="moog",wave_step=0.001)
        if first_loop:
            errors_df = pd.DataFrame(errors, index=np.arange(0,1,1))
            params_df = pd.DataFrame(params, index=np.arange(0,1,1))
        else:
            #Add the errors to a pandas dataframe
            errors_df = pd.concat([errors_df,pd.DataFrame(errors, index=list(f'{iteration_number}'))])
            #Add the parameters to a pandas dataframe
            params_df = pd.concat([params_df,pd.DataFrame(params, index=list(f'{iteration_number}'))])
        
        #Normalise the star spectrum using the template from the synthetic spectrum
        loop_spectrum, star_continuum_model = normalize_whole_spectrum_with_template(loop_spectrum, modeled_synth_spectrum)
        iteration_number += 1
        
        #Calculate the chi squared value and add it to a list
        chi2 = np.sum(((loop_spectrum['flux'][1:-1] - modeled_synth_spectrum['flux'][1:-1])**2)/ modeled_synth_spectrum['flux'][1:-1])
        chi2_df.loc[len(chi2_df)] = [iteration_number, chi2]
        
        if first_loop:
            first_loop=False
            continue
            
        # Check if the new chi-squared value is smaller than the initial value, and if the iteration number is less than 15
        if ((chi2_df.loc[len(chi2_df)-1,'chi2'] - chi2_df.loc[len(chi2_df)-2,'chi2']) < 5) and (iteration_number < 2):
            #prints the chi2 value for the iteration
            print(chi2_df.loc[len(chi2_df)-1,'chi2'])
            continue
        else:
            print(chi2_df)
            print(errors_df)
            print(params_df)
            break
    
    # chi2_df.to_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/{star_name}_chi2.txt')
    # errors_df.to_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/{star_name}_errors.txt')
    # params_df.to_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/{star_name}_params.txt')
    try:
        #Uni computer
        star_filename = f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/'+ f"{star_name}_adjusted.fits" 
    except:
        #Mac
        star_filename = f'/Users/quin/Desktop/2024_Data/Fixed_fits_files/{star_name}/'+ f"{star_name}_adjusted.fits"
    ispec.write_spectrum(loop_spectrum, star_filename)  
    print('Files saved')

    # try:
    #     pool=Pool(os.cpu_count()-1)
    #     pool.map(continuum_adjust, star)
    # finally:
    #     pool.close()
    #     pool.join()
    
    end = time.time()

    print(f'Time taken: {end - start}')
#%%
