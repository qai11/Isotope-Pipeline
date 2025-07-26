"""
Title: moog_abundance.py
Author: Quin Aicken Davies
Date: 25/06/2024

Description: Uses the iSpec code to generate a spectrum using moog.
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

#--- iSpec directory -------------------------------------------------------------
#ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/"
ispec_dir = '/home/users/qai11/iSpec_v20201001/'
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec

#%%

'''determine_abundances needs; atmosphere_layers, teff, logg, MH, alpha, 
linemasks, abundances, microturbulence_vel = 2.0, ignore=None, verbose=0, 
gui_queue=None, timeout=1800, isotopes=None, code="spectrum", tmp_dir=None'''

#if we set isotopes to not none it might work

#%%

determine_abundances = ispec.determine_abundances(atmosphere_layers=None, teff=6083, logg=4.1, MH=0.24, alpha=0.0, linemasks=None, abundances=None, microturbulence_vel=2.0, ignore=None, verbose=0, gui_queue=None, timeout=1800, isotopes=None, code="moog", tmp_dir=None)


#%%
'''Code from Heather for generating synthetic spectra and comparing to observed spectra'''

def compare_to_synth_spec(target_name, params):
   
   star_spectrum = ispec.read_spectrum("/home/users/qai11/Documents/Fixed_fits_files/" + target_name + '/hd_102870_25-06_avg.txt')
   #star_spectrum = rv_shift(target_name)
   wfilter = ispec.create_wavelength_filter(star_spectrum, wave_base=500.0, wave_top=675.0)
   star_spectrum = star_spectrum[wfilter]
   #star_spectrum = star_spectrum.reset_index(drop = True)
   #--- Synthesizing spectrum -----------------------------------------------------
   # Parameters
   teff = params[0] #5771.0
   logg = params[1] #4.44
   MH = params[2] #0.00
   alpha = params[3] #ispec.determine_abundance_enchancements(MH)
   microturbulence_vel = params[4] #ispec.estimate_vmic(teff, logg, MH) # 1.07
   macroturbulence = params[5]#ispec.estimate_vmac(teff, logg, MH) # 4.21
   vsini = params[6] #1.60 # Sun
   limb_darkening_coeff = 0.6
   resolution = 82000
   wave_step = (675 - 500) / len(star_spectrum)
   regions = None
   code = "MOOG"

   # Wavelengths to synthesis
   segments = np.recarray((1,),  dtype=[('wave_base', float), ('wave_top', float)])
   segments['wave_base'][0] = 500
   segments['wave_top'][0] = 675


   # Selected model amtosphere, linelist and solar abundances
   model = ispec_dir + "/input/atmospheres/MARCS.GES/"

   atomic_linelist_file = ispec_dir + "/input/linelists/transitions/GBS_LIST.480_800nm/atomic_lines.tsv"

   isotope_file = ispec_dir + "/input/isotopes/quinlist.MgH"

   # Load chemical information and linelist
   atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, wave_base=500, wave_top=675)
   atomic_linelist = atomic_linelist[atomic_linelist['theoretical_depth'] >= 0.01] # Select lines that have some minimal contribution in the sun

   isotopes = ispec.read_isotope_data(isotope_file)

   if "ATLAS" in model:
       solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.1998/stdatom.dat"
   else:
       # MARCS
       #solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.2007/stdatom.dat"
       solar_abundances_file = ispec_dir + "/input/abundances/Kurucz/stdatom.dat"

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
   synth_spectrum = ispec.create_spectrum_structure(np.arange(500, 675, wave_step))
   synth_spectrum['flux'] = ispec.generate_spectrum(synth_spectrum['waveobs'], \
           atmosphere_layers, teff, logg, MH, alpha, atomic_linelist, isotopes, solar_abundances, \
           fixed_abundances, microturbulence_vel = microturbulence_vel, \
           macroturbulence=macroturbulence, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, \
           R=resolution, regions=regions, verbose=1,
           code=code)
   
   plt.figure(1)
   plt.plot(star_spectrum["waveobs"], star_spectrum["flux"], label = "Star Spectrum")
   plt.plot(synth_spectrum["waveobs"], synth_spectrum["flux"], label = "Synthetic Spectrum")
   plt.legend()
   
   plt.figure(2)
   plt.plot(star_spectrum["waveobs"][1:-1], star_spectrum["flux"][1:-1] - synth_spectrum["flux"][1:-1], label = "Difference")
   mean = np.median(star_spectrum["flux"][1:-1] - synth_spectrum["flux"][1:-1])
   standard_deviation = np.std(star_spectrum["flux"][1:-1] - synth_spectrum["flux"][1:-1])
   plt.axhline(y = 0.1, c = "r", label = "Noise")
   plt.axhspan(mean-standard_deviation/2, mean+standard_deviation/2,  color = "g", label = "Zero", alpha = 0.5)
   plt.legend()
   plt.show()
   
   diff_spec = pd.DataFrame()
   diff_spec["waveobs"] = synth_spectrum["waveobs"][0:-1]
   diff_spec["flux"] = star_spectrum["flux"][0:-1] - synth_spectrum["flux"][0:-1]
   print(star_spectrum["flux"])
   # synth_spectrum = pd.DataFrame(synth_spectrum)
   print(synth_spectrum["flux"])
   diff_spec["err"] = synth_spectrum["err"][0:-1]
   print(diff_spec["flux"])
   diff_spec.to_csv("/home/users/qai11/Documents/quin-masters-code/Test_Fits/" + target_name + "_diff.txt", index = False, sep = "\t")
   
#%%
   
'Parameters order: teff, logg, MH, alpha, microturbulence_vel, macroturbulence, vsini'
output = compare_to_synth_spec('hd_102870', [6083, 4.1, 0.24, 0.0, 1.07, 4.21, 2])
# %%
