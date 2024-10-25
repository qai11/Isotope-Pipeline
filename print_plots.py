"""
Title: print_plots.py
Author: Quin Aicken Davies
Date: 14/10/2024

Description: Plots any extra plots I need for descriptions
"""
#%%
from scipy.stats import chisquare
from scipy import interpolate
from os import listdir
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from astropy.io import fits
#%%
def velocity_correction(wavelength, rv):
    ''' 
    Transforming the wavelengths in velocity space

    Lets hope my algebra is right. Should take the velocity and scale it depending on the rv_correction array that should be defined at the top of the notebook.

    Parameters
    ----------
    wavelength : array
        The array you want to scale
    rv : int
        The radial velocity correction

    Returns
    -------
    arr
        The wavelength in the rest frame
    '''

    return wavelength * (1+(-rv/300000))

def call_pymoogi(filename):
    #https://github.com/madamow/pymoogi/blob/master/README.md
    os.system('echo q | pymoogi ' + filename)
    
def get_region(r):
    if r == 0:
        lw = 5134.42
        uw = 5140.46
    elif r == 1:
        lw = 5134.42
        uw = 5134.85
    elif r == 2:
        lw = 5138.59
        uw = 5138.95
    elif r == 3:
        lw = 5140.04
        uw = 5140.46
    else:
        print('wavelength region error')
        lw = 0
        uw = 1
    return lw, uw

def get_lines(r):
    if r == 0:
        wl  = [5134.569, 5134.653, 5134.734,5138.710, 5138.768, 5138.785, 
        5138.823, 5138.860, 5140.202, 5140.251, 5140.286, 5140.302, 5140.358]
        iso = [24, 25, 26, 24, 25, 25, 26, 26, 24, 25, 25, 26, 26]
    elif r == 1:
        wl  = [5134.569, 5134.653, 5134.734]
        iso = [24, 25, 26]
    elif r == 2:
        wl  = [5138.710, 5138.768, 5138.785, 5138.823, 5138.860]
        iso = [24, 25, 25, 26, 26]
    elif r == 3:
        wl  = [5140.202, 5140.251, 5140.286, 5140.302, 5140.358]
        iso = [24, 25, 25, 26, 26]
    else:
        print('line region error')
        wl = [0]
        iso = [0]
    return wl, iso

def interp_smooth(raw, smooth):

    # perform a cubic spline on the data to make the wavelengths line up with eachother
    tck = interpolate.splrep(smooth.wavelength, smooth.flux, s=0)

    # evaluate the value of the spline at the wavelength points of the original spectra
    new_flux = interpolate.splev(raw.wavelength, tck, der=0)
    
    # add the model flux to the dataframe
    raw['model_flux'] = pd.Series(new_flux, index=raw.index)
    
    # return the dataframe with the new interpolated column in it 
    return raw
    
def make_model_plots(raw, smooth, output_filename, region, rv):

    # fig = plt.figure(constrained_layout=True, figsize = (8, 3))
    # gs = fig.add_gridspec(2, 1, height_ratios = [1, 0.3])
    # ax1 = fig.add_subplot(gs[0])
    # ax2 = fig.add_subplot(gs[1], sharex=ax1)

    # shift things backwards so it looks correct whe you plot things
    wavelength_shifted = velocity_correction(raw.wavelength, -rv)

    # Getting the plot bounds
    lw, uw = get_region(region)
    cropped_flux = raw[(raw.wavelength > lw) & (raw.wavelength < uw)].flux
    min_flux = cropped_flux.min()
    max_flux = cropped_flux.max()

    #getting the positions to plot the lines
    wl, iso = get_lines(region)
    text_24 = False
    text_25 = False
    text_26 = False

    # square surrounding the fitting region
    # ax1.fill_between([lw, uw], min_flux - 0.01, 1, facecolor = '#CCDBFD', alpha = 0.3)

    col24 = '#73454E'
    col25 = '#995C68'
    col26 = '#BA8C95'
    
    plt.figure(figsize=(10, 5))
    
    # for plotting the vertical lines of the isotopes
    for i in range(len(wl)):
        # if iso[i] == 24:
        #     ax1.plot([wl[i], wl[i]], [min_flux - 0.06, max_flux + 0.05], color = col24, linestyle = '--', dashes=(5, 1))
        #     if not text_24:
        #         ax1.text(wl[i] - 0.0, min_flux - 0.09, r'$^{24}\rm{MgH}$', color = col24, fontsize=10)
        #         text_24 = True
        # if iso[i] == 25:
        #     ax1.plot([wl[i], wl[i]], [min_flux - 0.04, max_flux + 0.05], color = col25,linestyle = '--', dashes=(3, 1))
        #     if not text_25:
        #         ax1.text(wl[i] - 0.0, min_flux - 0.07, r'$^{25}\rm{MgH}$', color = col25, fontsize=10)
        #         text_25 = True
        if iso[i] == 26:
            plt.plot([wl[i], wl[i]], [min_flux - 0.02, max_flux + 0.05], color = col26,linestyle = '--', dashes=(1, 1))
            if not text_26:
                plt.text(wl[i] - 0.0, min_flux - 0.05, r'$^{26}\rm{MgH}$', color = col26, fontsize=10)
                text_26 = True

    plt.plot(wavelength_shifted, raw.model_flux, color ='#E1A4A7', label = 'model')
    # plt.plot(smooth.wavelength, smooth.flux, color ='#eca1a6', label = 'model')
    plt.plot(wavelength_shifted, raw.flux, '.', color ='#2a9d8f', label = 'spectra')
    #ax1.plot(smooth.wavelength, smooth.flux, '.', color ='#d6cbd3', label = 'model')

    # ax2.plot(wavelength_shifted, raw.flux - raw.model_flux, color ='#ada397')

    plt.legend(frameon=False)
    plt.set_xlim(lw - 0.2, uw + 0.2)
    plt.set_xlabel(r'Wavelength ($\AA$)')

    plt.set_ylim(min_flux - 0.1, max_flux + 0.03)
    plt.set_ylabel('Norm. flux')

    # Limits for the residual plot
    # ax2.set_ylim(-0.024,0.024)
    

    # ax1.ticklabel_format(useOffset=False)
    # plt.setp(ax1.get_xticklabels(), visible=False)
    # ax2.ticklabel_format(useOffset=False)

    # ax1.tick_params(direction='in', axis='both', which='both', bottom=True,top=True, left=True, right=True)
    # ax1.xaxis.set_minor_locator(AutoMinorLocator())
    # ax1.yaxis.set_minor_locator(AutoMinorLocator())
    # ax2.tick_params(direction='in', axis='both', which='both', bottom=True,top=True, left=True, right=True)
    # ax2.xaxis.set_minor_locator(AutoMinorLocator())
    # ax2.yaxis.set_minor_locator(AutoMinorLocator())
    # gs.update(wspace=0, hspace=0, left=0, right=1, bottom=0, top=1)
    

    # plt.savefig('./plots/'+ output_filename, facecolor='white', bbox_inches='tight', dpi = 300)
    # plt.close()
    
def make_temp_file(filename):
    # will need these files in your directory - wont make them apparently...
    f = open(filename, "a+") 
    f.write('')
    f.close() 

def generate_parameter_string(raw_spec_filename, in_filename, out_filename, wavelength_region, par,star_name):
    # I doubt I'll ever want to change these so initialise them here
    standard_out = 'out1'
    summary_out  = 'out2'

    # will need these files in your directory - wont make them apparently...
    make_temp_file(standard_out)
    make_temp_file(summary_out)
    make_temp_file(out_filename)
    #par inputs the percentage ratios of the isotopes
    # print(str(par['i_24']))
    # print(str(par['i_25']))
    # print(str(par['i_26']))
    # print(str(par['mg']))
    # print(str(par['rv']))
    # print(wavelength_region)
    #t4070g040m18.newmod
    par_string = "synth\n" +\
    "standard_out   '" + standard_out +"'\n"                    + \
    "summary_out    '" + summary_out +"'\n"                     + \
    "smoothed_out   '" + out_filename +"'\n"                    + \
    f"model_in       '{star_name}_atmosphere.moog'\n"                          + \
    "lines_in       'quinlinelist.in'\n"                           + \
    "observed_in    '" + raw_spec_filename +"'\n"               + \
    "atmosphere    1\n"                                         + \
    "molecules     2\n"                                         + \
    "lines         2\n"                                         + \
    "flux/int      0\n"                                         + \
    "plotpars      1\n"                                         + \
    wavelength_region + " 0.15 1.05\n"                          + \
    str(par['rv']) + "      0.000   0.000    1.00\n"                   + \
    "d          0.047 0.0 0.0 "+ str(par['s']) +" 0.0\n"        + \
    "abundances   4    1\n"                                     + \
    "6            0.0500000\n"                                  + \
    "12           " + str(par['mg']) + "\n"                     + \
    "24           0.50000\n"                                    + \
    "26           0.20000\n"                                    + \
    "isotopes      4    1\n"                                    + \
    "606.01212     5.0\n"                                       + \
    "112.00124     "+ str(par['i_24']) +"\n"                    + \
    "112.00125     "+ str(par['i_25']) +"\n"                    + \
    "112.00126     "+ str(par['i_26']) +"\n"                    + \
    "obspectrum 5\n"                                            + \
    "synlimits\n"                                               + \
    wavelength_region + " 0.01 5.0\n"                           + \
    "plot 2\n"                                                  + \
    "damping 2\n"


    # writing that string to a file 
    par_file  = open(in_filename, "w+") 
    par_file.write(par_string)
    par_file.close() 
    return in_filename, out_filename

def read_raw_spectra(filename):
    return pd.read_table(filename, sep="\s+", usecols=[0,1], 
                         header=0, names = ['wavelength', 'flux'])

def read_smoothed_spectra(filename, rv):
    # different to reading raw spectra because we have to skip some headder rows
    smooth = pd.read_table(filename, sep="\s+", header=None, skiprows = [0,1],
                         names = ['wavelength', 'flux'])
    # run_interpolation for the values of the raw spectra wavelength
    smooth.wavelength = velocity_correction(smooth.wavelength, rv)
    return smooth

def make_filenames(par, prefix):
    str_s = str(round(par['s'],   2)).replace('.', '')
    str_mg = str(round(par['mg'],   3)).replace('.', '')
    str_24 = str(round(par['i_24'], 3)).replace('.', '')
    str_25 = str(round(par['i_25'], 3)).replace('.', '')
    str_26 = str(round(par['i_26'], 3)).replace('.', '')
    str_rv = str(round(par['rv'],   2)).replace('.', '')

    return prefix + '_s'+ str_s +'_mg'+ str_mg + '_i' \
     + str_24 + '_' + str_25  + '_' + str_26 + '_rv' + str_rv
     
def optimise_model_fit(raw_spec_filename, raw_spectra, region, wavelength_region, guess,star_name):

    # creating the in and out filenames based on the guess parameters
    in_filename  = make_filenames(guess, 'in')
    out_filename = make_filenames(guess, 'out')


    # creates a parameter string in the directory that moog can read
    generate_parameter_string(raw_spec_filename, in_filename, out_filename, wavelength_region, guess,star_name)

    # create the smoothed spectra by calling pymoogi
    smoothed_spectrum = call_pymoogi(in_filename)
    print(smoothed_spectrum)
    call_pymoogi(in_filename)

    # read in the smoothed model spectra and calculate the chi squared value
    # cs = get_chi_squared(raw_spectra, out_filename, region, guess, make_plot = True)
    cs = None
    
    # return a dataframe with a single row (to be added to a larger df later)
    return pd.DataFrame({'filename'   : out_filename, 
                         'chi_squared': cs, 
                         's'          : guess['s'],
                         'mg'         : guess['mg'],
                         'i_24'       : guess['i_24'],
                         'i_25'       : guess['i_25'],
                         'i_26'       : guess['i_26'],
                         'rv'         : guess['rv'],
                        #  'ratio'      : calc_ratio(guess['i_24'], guess['i_25'], guess['i_26'])
                         }, index=[1])
    
def initial_guess():

    s = 8.41 #has to be here for some reason or it breaks, seems to break move below 1.4
    mg = 0.055
    i_24 = 1
    i_25 = 6
    i_26 = 9
    rv = 0

    # return the guess as a dictionary
    return {'s'    : s,
            'mg'   : mg, 
            'i_24' : i_24, 
            'i_25' : i_25, 
            'i_26' : i_26, 
            'rv'   : rv}
#0.9 for plots
mg = initial_guess()

def get_wavelength_region(raw_wavelength):
    lower_wavelength = raw_wavelength[0]
    upper_wavelength = raw_wavelength[len(raw_wavelength)-1] # -1 isnt working for some reason
    return str(np.round(lower_wavelength, 2)) + ' ' + str(np.round(upper_wavelength, 2)) 
    
def model_finder(star_name):
    
    data_path = f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests/'
    region = 1
    os.chdir(data_path)
    os.system('mkdir plots')
    # initial guesses as a dictionary
    guess = initial_guess()

    raw_spec_filename = f'{star_name}_5100-5200.txt'
    raw_spectra       = read_raw_spectra(raw_spec_filename)
    wavelength_region = get_wavelength_region(raw_spectra.wavelength)

    # add the first chi_squyared value to the dataframe
    chi_df = optimise_model_fit(raw_spec_filename, raw_spectra, 
                                region, wavelength_region, guess,star_name)

star_name = 'hd_157244'
model_finder(star_name)
# %%
'''idk what happened above but it made the output file for some reason so ill print that'''
mg24 = str(mg['i_24']).replace('.', '')
mg25 = str(mg['i_25']).replace('.', '')
mg26 = str(mg['i_26']).replace('.', '')
mg_all= str(mg['mg']).replace('.', '')
s_all = str(mg['s']).replace('.', '')

#HD 102870
# smoothed = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/hd_102870/out_s841_mg{mg_all}_i{mg24}_{mg25}_{mg26}_rv0', sep="     ", header=None, skiprows = [0,1])
# raw = pd.read_csv('/home/users/qai11/Documents/Fixed_fits_files/hd_102870/hd_102870_5100-5200.txt', sep="	", header=None)
#HD128620
smoothed = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests/out_s{s_all}_mg{mg_all}_i{mg24}_{mg25}_{mg26}_rv0', sep="     ", header=None, skiprows = [0,1])
raw = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests/{star_name}_5100-5200.txt', sep="	", header=None)

plt.figure(figsize=(8, 4))
plt.plot(smoothed[0], smoothed[1])
plt.plot(raw[0], raw[1])
plt.legend(['Model', 'Observed'])
# plt.xlim(5133.8, 5135.5) #Region 1
# plt.xlim(5133, 5135)
plt.xlim(5130,5140)
# plt.xlim(5138.59, 5138.95) #Region 2
# plt.xlim(5140.04, 5140.46) #Region 3
# plt.xlim(5138.59,5138.95)
#HD102870
# plt.ylim(0.86, 1) #Region 1
# plt.ylim(0.92, 1) #Region 2
# plt.ylim(0.94, 1) #Region 3
#HD128620
# plt.ylim(0.84, 1.05) #Region 1
# plt.ylim(0.8, 1) #Region 2
# plt.ylim(0.85, 1) #Region 3

plt.xlabel('Wavelength ($\AA$)',fontsize=14)
plt.ylabel('Norm. Flux',fontsize=14)
# if mg['i_26'] == 1:
#     print("Mg26")
#     plt.text(5134.376, 0.87, 'Mg26 \n5134.376', fontsize=12, color='black',horizontalalignment='center')
#     plt.text(5134.734, 0.87, 'Mg26 \n5134.734', fontsize=12, color='black',horizontalalignment='center')
#     plt.annotate('', xy=(5134.376, 0.97), xytext=(5134.376, 0.885), arrowprops=dict(arrowstyle='->', color='black'))
#     plt.annotate('', xy=(5134.734, 0.96), xytext=(5134.734, 0.885), arrowprops=dict(arrowstyle='->', color='black'))
#     plt.savefig('/home/users/qai11/Documents/Masters_Figures/Method/Model_Mg26.png', dpi=300, bbox_inches='tight')
# elif mg['i_25'] == 1:
#     print("Mg25")
#     plt.text(5134.295, 0.87, 'Mg25 \n5134.295', fontsize=12, color='black',horizontalalignment='center')
#     plt.text(5134.656, 0.87, 'Mg25 \n5134.656', fontsize=12, color='black',horizontalalignment='center')
#     plt.text(5135.160, 0.87, 'Mg25 \n5135.160', fontsize=12, color='black',horizontalalignment='center')
#     plt.annotate('', xy=(5134.295, 0.969), xytext=(5134.295, 0.885), arrowprops=dict(arrowstyle='->', color='black'))
#     plt.annotate('', xy=(5134.656, 0.958), xytext=(5134.656, 0.885), arrowprops=dict(arrowstyle='->', color='black'))
#     plt.annotate('', xy=(5135.160, 0.973), xytext=(5135.160, 0.885), arrowprops=dict(arrowstyle='->', color='black'))
#     plt.savefig('/home/users/qai11/Documents/Masters_Figures/Method/Model_Mg25.png', dpi=300, bbox_inches='tight')
# elif mg['i_24'] == 1:
#     print("Mg24")
#     plt.text(5134.208, 0.87, 'Mg24 \n5134.208', fontsize=12, color='black',horizontalalignment='center')
#     plt.text(5134.570, 0.87, 'Mg24 \n5134.570', fontsize=12, color='black',horizontalalignment='center')
#     plt.text(5135.111, 0.87, 'Mg24 \n5135.111', fontsize=12, color='black',horizontalalignment='center')
#     plt.annotate('', xy=(5134.208, 0.963), xytext=(5134.208, 0.885), arrowprops=dict(arrowstyle='->', color='black'))
#     plt.annotate('', xy=(5134.570, 0.96), xytext=(5134.570, 0.885), arrowprops=dict(arrowstyle='->', color='black'))
#     plt.annotate('', xy=(5135.111, 0.973), xytext=(5135.111, 0.885), arrowprops=dict(arrowstyle='->', color='black'))
#     plt.savefig('/home/users/qai11/Documents/Masters_Figures/Method/Model_Mg24.png', dpi=300, bbox_inches='tight')
# # else:
# print("some cobination of Mg24, Mg25 and Mg26")
# plt.savefig('/home/users/qai11/Documents/Masters_Figures/Method/Model_Mg24_Mg25_Mg26.png', dpi=300, bbox_inches='tight')
plt.show()
#%%
def calc_ratio(i_24, i_25, i_26):
    i24_percentage=1/(0.01*i_24)
    i25_percentage=1/(0.01*i_25)
    i26_percentage=1/(0.01*i_26)

    isotope_sum = i24_percentage + i25_percentage + i26_percentage

    i24_ratio = (i24_percentage/isotope_sum) * 100
    i25_ratio = (i25_percentage/isotope_sum) * 100
    i26_ratio = (i26_percentage/isotope_sum) * 100
    
    return str(round(i24_ratio,2)) + '_' + str(round(i25_ratio,2)) + '_' + str(round(i26_ratio,2))

print(calc_ratio(mg['i_24'], mg['i_25'], mg['i_26']))

# %%

import astropy.io.fits as fits

def make_hr_diagram(gbs_filename, obs_filename):
    '''Makes an HR digram with inverted axes, Teff on x and Logg on y. Coloured by MH'''
    #Open the files
    gbs = fits.open(gbs_filename)
    obs = pd.read_csv(obs_filename, delimiter=',')
    #Print the files
    plt.figure(figsize=(8, 4))

    # Convert star names to lists
    obs_star_names = obs['ID2'].tolist()
    for i in range(len(obs_star_names)):
        obs_star_names[i] = obs_star_names[i].replace('hd_', 'HD')
    gbs_star_names = gbs[1].data['HD'].tolist()

    # Find common star names
    common_star_names = list(set(obs_star_names) & set(gbs_star_names))
    # Find stars that are in the obs list but not in the gbs list
    missing_star_names = list(set(obs_star_names) - set(gbs_star_names))
    # Find stars that are in the gbs list but not in the obs list
    gbs_only_star_names = list(set(gbs_star_names) - set(obs_star_names))
    
    
    # Plot stars that are in the gbs list but not in the obs list (gbs-only stars in blue)
    gbs_only_mask = [star in gbs_only_star_names for star in gbs_star_names]
    plt.scatter(gbs[1].data['TEFF'][gbs_only_mask], gbs[1].data['LOGG'][gbs_only_mask], c='#377eb8', s=10, lw=0,label='GBS-only Stars',)
    
    # Plot stars that are in the gbs list (common stars in red)
    common_mask = [star in common_star_names for star in gbs_star_names]
    plt.scatter(gbs[1].data['TEFF'][common_mask], gbs[1].data['LOGG'][common_mask], c='#ff7f00', s=20, label='Common Stars')

    # Plot stars that are in the obs list but not in the gbs list (missing stars in yellow)
    obs_missing_mask = [star in missing_star_names for star in obs_star_names]
    plt.scatter(obs['TEFF'][obs_missing_mask], obs['LOGG'][obs_missing_mask], c='#ff7f00', s=20, label='Missing Stars')
    
    # plt.scatter()
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    plt.xlabel('Teff (K)',fontsize=14)
    plt.ylabel('Logg (cm/s$^{2}$)',fontsize=14)
    plt.legend(['GBS-only Stars', 'Observed Stars'])
    plt.savefig('/home/users/qai11/Documents/Masters_Figures/Method/HR_diagram.png', dpi=300,bbox_inches='tight')
    # Show the plot
    plt.show()
    gbs.close()
    
make_hr_diagram('/home/users/qai11/Documents/quin-masters-code/gbs_stars.fit', '/home/users/qai11/Documents/quin-masters-code/Masters_stars.csv')

# %%
import astropy.io.fits as fits
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



"""Plot the synthetic spectrum the observed spectrum and adjusted spectrum"""
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
        # regions = None
        regions = ispec.read_segment_regions(ispec_dir + "/input/regions/Quin_segments.txt")
        # wave_base = wave_min
        # wave_top = wave_max


        # Selected model amtosphere, linelist and solar abundances
        model = ispec_dir + "/input/atmospheres/MARCS.GES/"
        
        atomic_linelist_file = ispec_dir + "/input/linelists/transitions/Quin_GES_LIST.420_920nm/atomic_lines.tsv"

        isotope_file = ispec_dir + "/input/isotopes/SPECTRUM.lst"

        # Load chemical information and linelist
        atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, wave_base=wave_base, wave_top=wave_top)
        # print(atomic_linelist)
        # atomic_linelist = atomic_linelist[atomic_linelist['theoretical_depth'] >= 0.01] # Select lines that have some minimal contribution in the sun
        # print(atomic_linelist['theoretical_depth'])
        print(len(atomic_linelist))
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
        # logging.info("Saving spectrum...")
        # synth_filename = "hd_102870_synth_%s.fits" % (code)
        # ispec.write_spectrum(synth_spectrum, synth_filename)
        return synth_spectrum

synth = synthesize_spectrum(teff=6080,logg=4.1,MH=0.24,vsini=2.0,wave_base=509.9, wave_top=520.1, code="moog",wave_step=0.001)
#%%
adjusted = ispec.read_spectrum('/home/users/qai11/Documents/Fixed_fits_files/hd_102870/hd_102870_adjusted.fits')
observed = ispec.read_spectrum('/home/users/qai11/Documents/Fixed_fits_files/hd_102870/J0354027.txt')
#%%
'''continuum adjust figure'''
# plt.figure(figsize=(10, 5))
# plt.plot(synth['waveobs'], synth['flux'])
# plt.plot(adjusted['waveobs'], adjusted['flux'])
# plt.plot(observed['waveobs'], observed['flux'])
# plt.legend(['Synthetic Spectrum', 'Adjusted Spectrum', 'Observed Spectrum'])
# plt.xlim(516.8, 517.5)
# plt.ylim(0, 1.1)
# plt.xlabel('Wavelength (nm)',fontsize=14)
# plt.ylabel('Norm. Flux',fontsize=14)
# plt.savefig('/home/users/qai11/Documents/Masters_Figures/Method/Continuum_adjust.png', dpi=150,bbox_inches='tight')
# plt.show()


# %%
'''Plot the parameters of the stars'''
# Import necessary libraries if not already imported
import pandas as pd
import matplotlib.pyplot as plt

# All stars including special stars
star_list = ['hd_100407', 'hd_102870', 'hd_128620', 'hd_128621', 'hd_11695', 
             'hd_146233', 'hd_157244', 'hd_160691', 'moon', 'hd_45588', 'hd_156098']

# Special stars you want to highlight
star_list_special = ['hd_45588', 'hd_156098']  # These stars will be highlighted in the plot

# Create empty DataFrames to store parameters and errors separately
all_params = pd.DataFrame()
all_errors = pd.DataFrame()

# Loop through each star and load parameters and errors
for star_name in star_list:
    # Load parameters and errors for each star
    params = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/parameters/{star_name}_final_params.txt', sep=',', index_col=None)
    errors = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/parameters/{star_name}_final_errors.txt', sep=',', index_col=None)
    
    # Add a column for the star name
    params['star'] = star_name
    errors['star'] = star_name

    # Append the data for this star to all_params and all_errors DataFrames
    all_params = pd.concat([all_params, params], ignore_index=True)
    all_errors = pd.concat([all_errors, errors], ignore_index=True)

# Open the benchmark values from the literature
benchmark_params = pd.read_csv('/home/users/qai11/Documents/quin-masters-code/Masters_stars.csv')

# Filter the data to only include stars that are in both DataFrames
merged_data = pd.merge(all_params, benchmark_params, left_on='star', right_on='ID2')

# Define the parameters to plot
parameters = ['teff', 'logg', 'MH']  # Parameters from all_params
benchmark_params_columns = ['TEFF', 'LOGG', 'FEH']  # Corresponding columns in benchmark_params

# Create a figure with three stacked plots
fig, axes = plt.subplots(3, 1, figsize=(6, 9))  # 3 rows, 1 column

# Loop through parameters and plot them
for i, (param, benchmark_param) in enumerate(zip(parameters, benchmark_params_columns)):
    # Calculate the difference between parameters and benchmark values
    y_values = merged_data[param] - merged_data[benchmark_param]
    x_values = merged_data['LOGG']

    # Extract corresponding uncertainties for the merged stars
    star_errors = merged_data['star'].map(all_errors.set_index('star')[param])

    # Plot regular stars with error bars
    regular_stars_mask = ~merged_data['star'].isin(star_list_special)
    axes[i].errorbar(x_values[regular_stars_mask], y_values[regular_stars_mask],
                     yerr=star_errors[regular_stars_mask], fmt='o', color='#377eb8', label='Benchmark Stars', markersize=2)

    # Plot special stars with error bars
    special_stars_mask = merged_data['star'].isin(star_list_special)
    axes[i].errorbar(x_values[special_stars_mask], y_values[special_stars_mask],
                     yerr=star_errors[special_stars_mask], fmt='o', color='#ff7f00', label='Unknown Stars', markersize=2)

    axes[i].set_xlabel('log($\it{g}$) (cm/s$^{2}$)', fontsize=12)

    # Manually modify the y-axis label
    if param == 'teff':
        axes[i].set_ylabel('T$_{\it{eff}}$ Difference (K)', fontsize=12)
    elif param == 'logg':
        axes[i].set_ylabel('log($\it{g}$) Difference (cm/s$^{2}$)', fontsize=12)
    elif param == 'MH':
        axes[i].set_ylabel('[M/H] Difference LTE (dex)', fontsize=12)

    # Add legend
    axes[i].legend(loc='best', fontsize=10)

# Adjust layout and show plot
plt.tight_layout()
# plt.savefig('/home/users/qai11/Documents/Masters_Figures/Method/Parameter_comparison.png', dpi=300, bbox_inches='tight')
plt.show()

# %%
'''Plot the MgH feature for all stars at 5133.8 to 5135.5'''
import pandas as pd
import matplotlib.pyplot as plt

# Open the benchmark values from the literature
benchmark_params = pd.read_csv('/home/users/qai11/Documents/quin-masters-code/Masters_stars.csv')

# Filter the data to only include stars that are in both DataFrames
merged_data = pd.merge(all_params, benchmark_params, left_on='star', right_on='ID2')

for star_name in merged_data['star']:
    # Read the spectrum file for the current star
    spectrum = ispec.read_spectrum(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/{star_name}_adjusted.fits')
    
    # Filter the benchmark_params to get the corresponding SPT, Teff, Logg, and Mh for the current star
    star_params = benchmark_params[benchmark_params['ID2'] == star_name]
    
    # Extract the values (if present)
    star_spt = star_params['SPT'].values
    star_teff = star_params['TEFF'].values
    star_logg = star_params['LOGG'].values
    star_mh = star_params['FEH'].values
    
    # Plot the spectrum
    plt.figure(figsize=(8, 4))
    plt.plot(spectrum['waveobs'] * 10, spectrum['flux'])
    plt.xlim(5133.8, 5135.5)
    plt.ylim(0.6, 1.1)
    
    # Create the legend with SPT, Teff, Logg, and Mh
    if len(star_spt) > 0 and len(star_teff) > 0 and len(star_logg) > 0 and len(star_mh) > 0:
        plt.legend([f"Star: {star_name}, SPT: {star_spt[0]}, Teff: {star_teff[0]}, Logg: {star_logg[0]}, [M/H]: {star_mh[0]}"])
    else:
        plt.legend(["SPT: N/A, Teff: N/A, Logg: N/A, [M/H]: N/A"])
    
    # plt.show()




# %%
