"""
Title: lbl_check.py
Author: Ethan Bull, Quin Aicken Davies
Date: 06/10/24

Description: This code plots all of the lines in the line by line files for a given star.
This is for visual inspection of the lines to check if the lines are being fit correctly.
"""
#%%
import numpy as np
from astropy.io import fits
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt
# from closest_factor import calc_closest_factors as factor
# from gen_even_dist import dist_stars
import os
import sys
import glob
#--- iSpec directory -------------------------------------------------------------
if os.path.exists('/home/users/qai11/iSpec_v20201001'):
    "Location of the files on Uni computer"
    ispec_dir = '/home/users/qai11/iSpec_v20201001'
else:
    "location of data on Mac"
    ispec_dir = '/Users/quin/Desktop/2024_Data/iSpec_v20230804'
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec

abunds_loc = '/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/'
obs_spectra_loc = '/home/users/qai11/Documents/Fixed_fits_files/' #final fitting spectrum
synth_spectra_loc = abunds_loc # for now
line_regions_loc = '/home/users/qai11/Documents/quin-masters-code/Linelists/'
save_loc = '/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/LBL_CHECK_NEW'

# First, open abundance lbl file to check wavelengths
# element = 'Ti' # put element name here
element = ["Eu","Ba","Mg"]#Rb isnt included
#"Nd", "Eu", "Sr", "Zr"
for element in element:
    # if not os.path.exists(save_loc + element):
    #     os.mkdir(save_loc + element)
    #LIST OF STAR NAMES
    stars =['hd_45588','hd_100407','hd_102870','hd_128620','hd_128621','hd_146233','hd_157244'
            ,'hd_160691','moon','hd_2151','hd_11695','hd_18907','hd_10700','hd_23249','hd_22049'
            ,'hd_18884','hd_165499','hd_156098']
    for star_name in stars:
        # star_name = 'hd_134606'
        abunds = pd.read_csv(abunds_loc + star_name + '/good_lbl/line_by_line_' + element + '.txt', sep=' ')
        #read in the wavelengths of the lines
        abunds_wave = pd.read_csv(abunds_loc + star_name + '/good_lbl/line_by_line_' + element + '.txt', sep=' ',dtype={'wave_peak':str})
        # also open the linelist to get more information
        line_regions = pd.read_csv(line_regions_loc + element + '_lines.csv', sep='\t')
        # combine abunds and line_regions
        abunds = pd.merge(abunds, line_regions, on='wave_peak')
            
        # open the observed_spectrum and synthetic spectrum
        obs_spec = ispec.read_spectrum(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/{star_name}_adjusted.fits')
        
        y_sub = 4
        x_sub = np.ceil(len(abunds) / y_sub).astype(int)
        #generate plot with subplots for each line
        fig, axs = plt.subplots(x_sub, 4, figsize=(4*y_sub, 3*x_sub))
        axs=axs.flatten()
        for i, line in abunds.iterrows():
            # define wavelength to plot
            wave_base = line['wave_base']
            wave_top = line['wave_top']
            #Read in the wavepeak from the abunds_wave file
            wave_peak = abunds_wave['wave_peak'][i]
            #YOU NEED TO ROUND THE STRINGS TO 5 DECIMAL PLACES as it wont match otherwise. Do this later
            synth_spec = glob.glob(synth_spectra_loc + star_name + '/synth_regions_' + element 
                                    + '/' + 'synth_' + element + '_' + wave_peak + '*')[0]
            synth_spec = pd.read_csv(synth_spec, sep=' ')
            
            extra_tolerance = 0.1 # extra tolerance to plot around the line
            obs_spec_region = obs_spec[(obs_spec['waveobs'] > wave_base - extra_tolerance) & (obs_spec['waveobs'] < wave_top + extra_tolerance)]
            synth_spec_region = synth_spec[(synth_spec['waveobs'] > wave_base - extra_tolerance) & (synth_spec['waveobs'] < wave_top + extra_tolerance)]
            # plot the observed spectrum
            axs[i].scatter(obs_spec_region['waveobs'], obs_spec_region['flux'], label='Observed Spectrum', c='k', s=1)
            # plot the synthetic spectrum
            axs[i].plot(synth_spec_region['waveobs'], synth_spec_region['flux'], label='Synthetic Spectrum', c='r')
            # plot the line region
            axs[i].axvline(x=line['wave_peak'], c='b', ls='--')
            axs[i].axvspan(line['wave_base'], line['wave_top'], color='y', alpha=0.4)
            axs[i].set_xlim(wave_base - extra_tolerance, wave_top + extra_tolerance)
            axs[i].set_title(f'{element} {line["wave_peak"]:.2f} nm, [X/H] = {line["[X/H]"]:.3f}')
            axs[i].set_xlabel('Wavelength (nm)')
            axs[i].set_ylabel('Flux')
        axs[0].legend()
        fig.tight_layout()
        diff_axes = len(axs) - len(abunds)
        for i in range(diff_axes):
            axs[-int(i+1)].set_axis_off()
        if not os.path.exists(f'{save_loc}/{element}'):
            os.mkdir(f'{save_loc}/{element}/')
        plt.savefig(f'{save_loc}/{element}/{star_name}.png', dpi = 300)
    plt.close()
    
    #%%
"""EDITS WITH THE GOOD LBL LISTS FOR TAKING THE AVERAGE FOR THE FINAL ABUNDANCE"""
import numpy as np
from astropy.io import fits
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt
# from closest_factor import calc_closest_factors as factor
# from gen_even_dist import dist_stars
import os
import sys
import glob
#--- iSpec directory -------------------------------------------------------------
if os.path.exists('/home/users/qai11/iSpec_v20201001'):
    "Location of the files on Uni computer"
    ispec_dir = '/home/users/qai11/iSpec_v20201001'
else:
    "location of data on Mac"
    ispec_dir = '/Users/quin/Desktop/2024_Data/iSpec_v20230804'
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec

abunds_loc = '/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/'
obs_spectra_loc = '/home/users/qai11/Documents/Fixed_fits_files/' #final fitting spectrum
synth_spectra_loc = abunds_loc # for now
line_regions_loc = '/home/users/qai11/Documents/quin-masters-code/Linelists/'
save_loc = '/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/LBL_CHECK'

# First, open abundance lbl file to check wavelengths
# element = 'Ti' # put element name here
element = ["Eu", "Sr", "Zr"]#Rb isnt included
#"Nd", "Eu", "Sr", "Zr"
for element in element:
    # if not os.path.exists(save_loc + element):
    #     os.mkdir(save_loc + element)
    #LIST OF STAR NAMES
    stars =['hd_45588','hd_100407','hd_102870','hd_128620','hd_128621','hd_146233','hd_157244'
            ,'hd_160691','moon','hd_2151','hd_11695','hd_18907','hd_10700','hd_23249','hd_22049'
            ,'hd_18884','hd_165499','hd_156098']
    for star_name in stars:
        # star_name = 'hd_134606'
        abunds = pd.read_csv(abunds_loc + star_name + '/good_lbl/line_by_line_' + element + '.txt', sep=' ')
        #read in the wavelengths of the lines
        abunds_wave = pd.read_csv(abunds_loc + star_name + '/good_lbl/line_by_line_' + element + '.txt', sep=' ',dtype={'wave_peak':str})
        # also open the linelist to get more information
        line_regions = pd.read_csv(line_regions_loc + element + '_lines.csv', sep='\t')
        # combine abunds and line_regions
        abunds = pd.merge(abunds, line_regions, on='wave_peak')
            
        # open the observed_spectrum and synthetic spectrum
        obs_spec = ispec.read_spectrum(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/{star_name}_adjusted.fits')
        
        y_sub = 4
        x_sub = np.ceil(len(abunds) / y_sub).astype(int)
        #generate plot with subplots for each line
        fig, axs = plt.subplots(x_sub, 4, figsize=(4*y_sub, 3*x_sub))
        axs=axs.flatten()
        for i, line in abunds.iterrows():
            # define wavelength to plot
            wave_base = line['wave_base']
            wave_top = line['wave_top']
            #Read in the wavepeak from the abunds_wave file
            wave_peak = abunds_wave['wave_peak'][i]
            #YOU NEED TO ROUND THE STRINGS TO 5 DECIMAL PLACES as it wont match otherwise. Do this later
            synth_spec = glob.glob(synth_spectra_loc + star_name + '/synth_regions_' + element 
                                    + '/' + 'synth_' + element + '_' + wave_peak + '*')[0]
            synth_spec = pd.read_csv(synth_spec, sep=' ')
            
            extra_tolerance = 0.1 # extra tolerance to plot around the line
            obs_spec_region = obs_spec[(obs_spec['waveobs'] > wave_base - extra_tolerance) & (obs_spec['waveobs'] < wave_top + extra_tolerance)]
            synth_spec_region = synth_spec[(synth_spec['waveobs'] > wave_base - extra_tolerance) & (synth_spec['waveobs'] < wave_top + extra_tolerance)]
            # plot the observed spectrum
            axs[i].scatter(obs_spec_region['waveobs'], obs_spec_region['flux'], label='Observed Spectrum', c='k', s=1)
            # plot the synthetic spectrum
            axs[i].plot(synth_spec_region['waveobs'], synth_spec_region['flux'], label='Synthetic Spectrum', c='r')
            # plot the line region
            axs[i].axvline(x=line['wave_peak'], c='b', ls='--')
            axs[i].axvspan(line['wave_base'], line['wave_top'], color='y', alpha=0.4)
            axs[i].set_xlim(wave_base - extra_tolerance, wave_top + extra_tolerance)
            axs[i].set_title(f'{element} {line["wave_peak"]:.2f} nm, [X/H] = {line["[X/H]"]:.3f}')
            axs[i].set_xlabel('Wavelength (nm)')
            axs[i].set_ylabel('Flux')
        axs[0].legend()
        fig.tight_layout()
        diff_axes = len(axs) - len(abunds)
        for i in range(diff_axes):
            axs[-int(i+1)].set_axis_off()
        if not os.path.exists(f'{save_loc}/{element}'):
            os.mkdir(f'{save_loc}/{element}/')
        plt.savefig(f'{save_loc}/{element}/{star_name}.png', dpi = 300)
    plt.close()