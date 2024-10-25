"""
Title: write_ispec_to_moog.py
Author: Quin Aicken Davies
Date: 01/10/2024

Description: Uses the iSpec to convert the atomic linelist to a MOOG readable format.
"""

#%%
import numpy as np
import glob
import os
from astropy.io import fits
import pandas as pd
import time
import scipy as sp
import sys

#--- iSpec directory -------------------------------------------------------------
#ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/"
ispec_dir = '/home/users/qai11/iSpec_v20201001/'
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec

from ispec.lines import __moog_write_atomic_linelist
from ispec.lines import write_atomic_linelist
from ispec.lines import read_atomic_linelist

#%%

#--- Read atomic linelist -------------------------------------------------------

linelist_save_loc = '/home/users/qai11/Documents/Fixed_fits_files/quinlinelist.in'
linelist_filename = '/home/users/qai11/iSpec_v20201001/input/linelists/transitions/Quin_GES_LIST.420_920nm/atomic_lines.tsv'
linelist = read_atomic_linelist(linelist_filename, wave_base=None, wave_top=None)

line_list = write_atomic_linelist(linelist, linelist_filename=linelist_save_loc, code='moog', tmp_dir=None)

# __moog_write_atomic_linelist(linelist, linelist_filename=None, tmp_dir=None)


# %%

def synthesize_spectrum(atomic_linelist_file=None, code="spectrum"):
    #--- Synthesizing spectrum -----------------------------------------------------
    # Parameters
    teff = 5771.0
    logg = 4.44
    MH = 0.0
    alpha = ispec.determine_abundance_enchancements(MH)
    microturbulence_vel = ispec.estimate_vmic(teff, logg, MH) # 1.07
    macroturbulence = ispec.estimate_vmac(teff, logg, MH) # 4.21
    vsini = 1.6 
    limb_darkening_coeff = 0.6
    resolution = 82000
    wave_step = 0.001

    # Wavelengths to synthesis
    #regions = ispec.read_segment_regions(ispec_dir + "/input/regions/fe_lines_segments.txt")
    regions = None
    wave_base = 510.0 # Magnesium triplet region
    wave_top = 540.0


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
    # atomic_linelist_file = ispec_linelist_file
    atomic_linelist_file = ispec_dir + "/input/linelists/transitions" + atomic_linelist_file
    
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
    # logging.info("Saving spectrum...")
    # synth_filename = "example_synth_%s.fits" % (code)
    # ispec.write_spectrum(synth_spectrum, synth_filename)
    
    return synth_spectrum

# synth_spectrum1 = synthesize_spectrum(atomic_linelist_file = '/GESv6_atom_nohfs_noiso.420_920nm/atomic_lines.tsv',code="moog")

synth_spectrum2 = synthesize_spectrum(atomic_linelist_file = '/Quin_GES_LIST.420_920nm/atomic_lines.tsv', code="moog")
#%%

# Plot the synthetic spectrum
import matplotlib.pyplot as plt
plt.plot(synth_spectrum1['waveobs'], synth_spectrum1['flux'])
plt.plot(synth_spectrum2['waveobs'], synth_spectrum2['flux'])
plt.xlim(510, 540.35)
plt.ylim(0.8,1)
plt.legend(['GES', 'Quin'])
plt.xlabel('Wavelength (Å)')
plt.ylabel('Flux')
# plt.xlim(510.5, 510.6)
# plt.savefig('/home/users/qai11/Documents/Masters_Figures/Method/GESV6_VS_GES_MOL.png')
#THERE Appears the be a small difference between the Cu1 lines in the two linelists
#one had molecules in this region and one does not at approximately 510.519-510.521
plt.show()
# %%
