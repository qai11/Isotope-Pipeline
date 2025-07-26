"""
Title: moog_ispec.py
Author: Quin Aicken Davies
Date: 03/07/2024

Description: converts moog line list for ispec
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
# Define the file paths
moog_linelist_file = '/home/users/qai11/Documents/quin-masters-code/quinlist.txt'
ispec_linelist_file = '/home/users/qai11/Documents/quin-masters-code/ispec_linelist.tsv'
#%%
# Read the MOOG linelist
moog_linelist = pd.read_csv(moog_linelist_file, delim_whitespace=True, header=None)
moog_linelist.columns = ['wave_A', 'species', 'excitation_potential', 'loggf', 'damping_potential','reference_code']

# Convert element to species notation required by iSpec
# For neutral elements, iSpec uses element + ".0", for singly ionized elements, element + ".1"
moog_linelist['element'] = moog_linelist['species'].apply(lambda x: f"{int(x)}.0" if x.is_integer() else f"{int(x)}.1")
moog_linelist['wave_nm'] = moog_linelist['wave_A'] / 10  # Convert wavelength from Å to nm

# Reorder and rename columns to match iSpec format
ispec_linelist = moog_linelist[['element','wave_A','wave_nm', 'loggf', 'excitation_potential', 'species']]
ispec_linelist.columns = ['element','wave_A','wave_nm', 'loggf', 'EP', 'species']

# Add damping parameters if needed
ispec_linelist['gamma'] = 0.0  # Van der Waals broadening (optional)
ispec_linelist['gamma6'] = 0.0  # Stark broadening (optional)
ispec_linelist['vdW'] = 0.0  # Natural broadening (optional)

# Save the converted linelist to CSV
ispec_linelist.to_csv(ispec_linelist_file, index=False, sep='\t')
print(f"Converted linelist saved to {ispec_linelist_file}")

# %%


def synthesize_spectrum(code="spectrum"):
    #--- Synthesizing spectrum -----------------------------------------------------
    # Parameters
    teff = 6083
    logg = 4.1
    MH = 0.24
    alpha = ispec.determine_abundance_enchancements(MH)
    microturbulence_vel = ispec.estimate_vmic(teff, logg, MH) # 1.07
    macroturbulence = ispec.estimate_vmac(teff, logg, MH) # 4.21
    vsini = 2 
    limb_darkening_coeff = 0.6
    resolution = 82000
    wave_step = 0.001

    # Wavelengths to synthesis
    #regions = ispec.read_segment_regions(ispec_dir + "/input/regions/fe_lines_segments.txt")
    regions = None
    wave_base = 515.0 # Magnesium triplet region
    wave_top = 525.0


    # Selected model amtosphere, linelist and solar abundances
    #model = ispec_dir + "/input/atmospheres/MARCS/"
    model = ispec_dir + "/input/atmospheres/MARCS.GES/"
    #model = ispec_dir + "/input/atmospheres/MARCS.APOGEE/"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.APOGEE/"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Castelli/"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kurucz/"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kirby/"

    #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/VALD.300_1100nm/atomic_lines.tsv"
    #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/VALD.1100_2400nm/atomic_lines.tsv"
    # atomic_linelist_file = ispec_dir + "/input/linelists/transitions/GESv6_atom_hfs_iso.420_920nm/atomic_lines.tsv"
    #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/GESv6_atom_nohfs_noiso.420_920nm/atomic_lines.tsv"
    atomic_linelist_file = ispec_linelist_file
    
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
    #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2005/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2009/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Anders.1989/stdatom.dat"

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
    ##--- Save spectrum ------------------------------------------------------------
    logging.info("Saving spectrum...")
    synth_filename = "example_synth_%s.fits" % (code)
    ispec.write_spectrum(synth_spectrum, synth_filename)

# # Plot the synthetic spectrum
# import matplotlib.pyplot as plt
# plt.plot(spectrum['wave'], spectrum['flux'])
# plt.xlabel('Wavelength (Å)')
# plt.ylabel('Flux')
# plt.show()

synthesize_spectrum(code="moog")

# %%
