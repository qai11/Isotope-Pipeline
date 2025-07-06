"""
Title: Isotope_pipeline_V2.py
Author: Quin Aicken Davies
Date: 17/06/2025

Description: With Failure of Isotope_pipeline_V1.py, this script is a rework of the original pipeline.
This pipeline reworks the fitting routine to use a more robust fitting method, 
and includes additional features such as: Giant or Dwarf star Selection.
"""

#%%
from scipy.stats import chisquare
from scipy import interpolate
from scipy import optimize
from scipy.optimize import least_squares
from os import listdir
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import sys
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

from ispec.spectrum import __convolve_spectrum

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
        lw = 5138.55
        uw = 5138.95
    elif r == 3:
        lw = 5140.04
        uw = 5140.46
    elif r == 4:
        lw = 5134.0
        uw = 5134.4
    elif r == 5:
        lw = 5134.9
        uw = 5135.3
    elif r == 6:
        lw = 5135.9
        uw = 5136.3
    elif r == 7:
        lw = 5136.2
        uw = 5136.6
    elif r == 8:
        lw = 5138.2
        uw = 5138.6
    elif r == 9:
        lw = 5141.0
        uw = 5141.45
    elif r == 10:
        lw = 5133.0
        uw = 5133.4
    else:
        print('wavelength region error')
        lw = 0
        uw = 1
    return lw, uw


def interp_smooth(obs_flux, smooth):
    # perform a cubic spline on the data to make the wavelengths line up with eachother
    tck = interpolate.splrep(smooth.wavelength, smooth.flux, s=0)

    # evaluate the value of the spline at the wavelength points of the original spectra
    new_flux = interpolate.splev(raw.wavelength, tck, der=0)
    
    # add the model flux to the dataframe
    raw['model_flux'] = pd.Series(new_flux, index=raw.index)
    
    # return the dataframe with the new interpolated column in it 
    return obs_flux


def make_temp_file(filename):
    # will need these files in your directory - wont make them apparently...
    f  = open(filename, "a+") 
    f.write('')
    f.close() 

def generate_parameter_string(raw_spec_filename, in_filename, out_filename, wavelength_region, par,star_name,linelist,vsini,Fe,CN,CC):
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
    # par_string = "synth\n" +\
    # "standard_out   '" + standard_out +"'\n"                    + \
    # "summary_out    '" + summary_out +"'\n"                     + \
    # "smoothed_out   '" + out_filename +"'\n"                    + \
    # f"model_in       '{star_name}_atmosphere.moog'\n"                          + \
    # f"lines_in       '{linelist}'\n"                           + \
    # "observed_in    '" + raw_spec_filename +"'\n"               + \
    # "atmosphere    1\n"                                         + \
    # "molecules     2\n"                                         + \
    # "lines         2\n"                                         + \
    # "flux/int      0\n"                                         + \
    # "plotpars      1\n"                                         + \
    # wavelength_region + " 0.15 1.05\n"                          + \
    # str(par['rv']) + "      0.000   0.000    1.00\n"                   + \
    # "d          0.047 0.0 0.0 "+ str(par['s']) +" 0.0\n"        + \
    # "abundances   3    1\n"                                     + \
    # "6            0.000001\n"                                  + \
    # "12           " + str(par['mg']) + "\n"                     + \
    # "22           0.20000\n"                                    + \
    # "isotopes      5    1\n"                                    + \
    # "607.01214     1.0\n"                                       + \
    # "606.01212     3.0\n"                                       + \
    # "112.00124     "+ str(par['i_24']) +"\n"                    + \
    # "112.00125     "+ str(par['i_25']) +"\n"                    + \
    # "112.00126     "+ str(par['i_26']) +"\n"                    + \
    # "obspectrum 5\n"                                            + \
    # "synlimits\n"                                               + \
    # wavelength_region + " 0.01 5.0\n"                           + \
    # "plot 2\n"                                                  + \
    # "damping 2\n"
    par_string = "synth\n" +\
    "standard_out   '" + standard_out +"'\n"                    + \
    "summary_out    '" + summary_out +"'\n"                     + \
    "smoothed_out   '" + out_filename +"'\n"                    + \
    f"model_in       '{star_name}_atmosphere.moog'\n"                          + \
    f"lines_in       '{linelist}'\n"                           + \
    "observed_in    '" + raw_spec_filename +"'\n"               + \
    "atmosphere    1\n"                                         + \
    "molecules     2\n"                                         + \
    "lines         2\n"                                         + \
    "flux/int      0\n"                                         + \
    "plotpars      1\n"                                         + \
    wavelength_region + " 0.15 1.05\n"                          + \
    str(par['rv']) + "      0.000   0.000    1.00\n"                   + \
    "d          0.06 "+str(vsini)+" 0.6 "+ str(par['s']) +" 0.0\n"        + \
    "abundances   5    1\n"                                     + \
    "6            0.200000\n"                                  + \
    "12           " + str(par['mg']) + "\n"                     + \
    "22           0.20000\n"                                    + \
    "24           0.10000\n"                                    + \
    "26           " + str(Fe) + "\n"                            + \
    "isotopes      5    1\n"                                    + \
    "607.01214     " + str(CN) + "\n"                            + \
    "606.01212     " + str(CC) + "\n"                            + \
    "112.00124     "+ str(par['i_24']) +"\n"                    + \
    "112.00125     "+ str(par['i_25']) +"\n"                    + \
    "112.00126     "+ str(par['i_26']) +"\n"                    + \
    "obspectrum 5\n"                                            + \
    "synlimits\n"                                               + \
    wavelength_region + " 0.01 5.0\n"                           + \
    "plot 1\n"                                                  + \
    "damping 0\n"

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


def get_chi_squared(obs_flux, out_filename, region, guess, vsini, make_plot=False):
    
    # read in the smoothed data
    model_flux = read_smoothed_spectra(out_filename, guess['rv'])
    obs_flux = interp_smooth(obs_flux, smooth)
    residual = obs_flux[region[0]:region[1]] - model_flux[region[0]:region[1]]
    return np.sum(residual**2)

def make_filenames(par, prefix):
    str_s = str(round(par['s'],   2)).replace('.', '')
    str_mg = str(round(par['mg'],   3)).replace('.', '')
    str_24 = str(round(par['i_24'], 3)).replace('.', '')
    str_25 = str(round(par['i_25'], 3)).replace('.', '')
    str_26 = str(round(par['i_26'], 3)).replace('.', '')
    str_rv = str(round(par['rv'],   2)).replace('.', '')

    return prefix + '_s'+ str_s +'_mg'+ str_mg + '_i' \
     + str_24 + '_' + str_25  + '_' + str_26 + '_rv' + str_rv

def get_wavelength_region(region):
    '''Try cutting out the range'''
    # lower_wavelength = raw_wavelength[0]
    # upper_wavelength = raw_wavelength[len(raw_wavelength)-1] # -1 isnt working for some reason
    if region == 1:
        '''region 1'''
        lower_wavelength = 5131
        upper_wavelength = 5138
    if region == 2:
        '''region 2'''
        lower_wavelength =  5135
        upper_wavelength = 5142
    if region == 3:
        '''region 3'''
        lower_wavelength = 5136
        upper_wavelength = 5143
    if region == 4:
        '''region 1'''
        lower_wavelength = 5131
        upper_wavelength = 5138
    if region == 5:
        '''region 2'''
        lower_wavelength =  5131
        upper_wavelength = 5138
    if region == 6:
        '''region 3'''
        lower_wavelength = 5133
        upper_wavelength = 5139
    if region == 7:
        '''region 1'''
        lower_wavelength = 5133
        upper_wavelength = 5139
    if region == 8:
        '''region 2'''
        lower_wavelength =  5135
        upper_wavelength = 5142
    if region == 9:
        '''region 3'''
        lower_wavelength = 5136
        upper_wavelength = 5143
    if region == 10:
        '''region 3'''
        lower_wavelength = 5131
        upper_wavelength = 5138 
    
    # print(str(np.round(lower_wavelength, 2)) + ' ' + str(np.round(upper_wavelength, 2)) )
    return str(np.round(lower_wavelength, 2)) + ' ' + str(np.round(upper_wavelength, 2)) 


def optimise_model_fit(raw_spec_filename, raw_spectra, region, wavelength_region, guess, star_name, linelist, vsini, Fe, CN, CC):

    # Initial guess: s, i_24, i_25, i_26
    x0 = [guess['s'], guess['i_24'], guess['i_25'], guess['i_26']]
    print(x0)

    # Optional: set bounds if needed
    bounds = ([4, 0.1, 1, 1], [9, 8, 8, 8])  # bounds

    result = least_squares(
        residuals_to_minimize,
        x0=x0,
        bounds=bounds,
        args=(raw_spec_filename, raw_spectra, region, wavelength_region, guess, star_name, linelist, vsini, Fe, CN, CC),
        loss='soft_l1'  # Robust to outliers; or 'linear' if clean
    )

    # Use final parameters to recompute everything
    guess['s'], guess['i_24'], guess['i_25'], guess['i_26'] = result.x
    in_filename  = make_filenames(guess, 'in')
    out_filename = make_filenames(guess, 'out')
    generate_parameter_string(raw_spec_filename, in_filename, out_filename, wavelength_region, guess, star_name, linelist, vsini, Fe, CN, CC)
    #Call pymoogi with the generated input file
    call_pymoogi(in_filename)

    # Get final chi-squared
    cs = get_chi_squared(raw_spectra, out_filename, region, guess, vsini, make_plot=False)

    return pd.DataFrame({
        'filename'   : out_filename, 
        'chi_squared': cs, 
        's'          : guess['s'],
        'mg'         : guess['mg'],
        'i_24'       : guess['i_24'],
        'i_25'       : guess['i_25'],
        'i_26'       : guess['i_26'],
        'rv'         : guess['rv'],
        'ratio'      : calc_ratio(guess['i_24'], guess['i_25'], guess['i_26'])
    }, index=[1])

def residuals_to_minimize(x, raw_spec_filename, raw_spectra, region, wavelength_region, guess, star_name, linelist, vsini, Fe, CN, CC):
    # Update guess dictionary with current parameters
    guess = guess.copy()
    guess['s']     = x[0]
    guess['i_24']  = x[1]
    guess['i_25']  = x[2]
    guess['i_26']  = x[3]

    # Generate filenames
    in_filename  = make_filenames(guess, 'in')
    out_filename = make_filenames(guess, 'out')

    # Generate MOOG input + run pymoogi
    generate_parameter_string(raw_spec_filename, in_filename, out_filename, wavelength_region, guess, star_name, linelist, vsini, Fe, CN, CC)
    call_pymoogi(in_filename)

    # Get the smoothed spectrum and residuals
    residuals = get_residuals(raw_spectra, out_filename, region)  # <-- You will need to write or already have this
    return residuals

def get_residuals(obs_flux, syn_filename, region):
    # Load model spectrum
    model_flux = read_spectrum(syn_filename)  # or however you access it
    # Mask/select the region
    observed = obs_flux[region[0]:region[1]]
    model = model_flux[region[0]:region[1]]
    return observed - model

def model_finder(star_name,linelist,region,vsini,MgH,Fe,CN,CC):
    data_path = f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests_paper/'
    'change wavelength range'
    region = region
    os.chdir(data_path)
    os.system('mkdir plots')
    # initial guesses as a dictionary
    guess = initial_guess(MgH)
    #Open the raw spectra file
    raw_spec_filename = data_path + f'{star_name}_5100-5200.txt'
    raw_spectra       = read_raw_spectra(raw_spec_filename)
    #Open the wavelength region based on what region is selected
    wavelength_region = get_wavelength_region(region)

    # add the first chi_squyared value to the dataframe
    chi_df = optimise_model_fit(raw_spec_filename, raw_spectra, 
                                region, wavelength_region, guess,star_name,linelist,vsini,Fe,CN,CC)



    #THE STUFF BELOW MIGHT NEEED TO BE REMOVED
    best_guess_chi = chi_df.chi_squared.iloc[0] # should only be 1 thing in the df atm

    # currently this is the best guess we have
    best_guess = guess
    while len(chi_df) < 300:
        # add the neighbours to the dataframe
        chi_df = find_minimum_neighbour(raw_spec_filename, raw_spectra, 
                        wavelength_region, region, best_guess, chi_df,star_name,linelist,vsini,Fe,CN,CC)
        
        # get the best chi-squared fit
        chi_df = chi_df.sort_values(by = ['chi_squared'])
        minimum_model = chi_df.iloc[0]
        print(chi_df)
        print('minimum model identified: ', minimum_model)

        new_best_guess_chi = minimum_model.chi_squared
        new_best_guess = reconstruct_min_chi(minimum_model)
        print('current best guess: ', best_guess_chi)
        print('new best guess: ', new_best_guess_chi)
        
        if best_guess_chi > new_best_guess_chi:
            best_guess = new_best_guess
            best_guess_chi = new_best_guess_chi
            print('now using ', best_guess)
        elif best_guess_chi == new_best_guess_chi:
            print('could not find a better option - exiting')
            break

    return chi_df
        
def calc_ratio(i_24, i_25, i_26):
    i24_percentage=1/(0.01*i_24)
    i25_percentage=1/(0.01*i_25)
    i26_percentage=1/(0.01*i_26)

    isotope_sum = i24_percentage + i25_percentage + i26_percentage

    i24_ratio = (i24_percentage/isotope_sum) * 100
    i25_ratio = (i25_percentage/isotope_sum) * 100
    i26_ratio = (i26_percentage/isotope_sum) * 100
    
    return str(round(i24_ratio,2)) + '_' + str(round(i25_ratio,2)) + '_' + str(round(i26_ratio,2))


def calc_moog(r_24, r_25, r_26):
    i24=1/(0.01*r_24)
    i25=1/(0.01*r_24)
    i26=1/(0.01*r_24)
    return [i24, i25, i26]

def calc_moog_string(r_24, r_25, r_26):
    i24=1/(0.01*r_24)
    i25=1/(0.01*r_24)
    i26=1/(0.01*r_24)
    return str(round(i24,2)) + '_' + str(round(i25,2)) + '_' + str(round(i26,2))

def initial_guess(MgH):
    # initial guess for the parameters
    # s = 7.5
    # mg = MgH
    # i_24 = 2
    # i_25 = 15
    # i_26 = 13
    # rv = 0
    # New guess after first pass.
    s = 8.9
    mg = MgH
    i_24 = 0.8
    i_25 = 15
    i_26 = 3.6
    rv = 0
    #region best fit first pass
    # s = 8.9
    # mg = 0.03
    # i_24 = 3.7
    # i_25 = 5.5
    # i_26 = 15
    # rv = 0
    # return the guess as a dictionary
    return {'s'    : s, 
            'mg'   : mg, 
            'i_24' : i_24, 
            'i_25' : i_25, 
            'i_26' : i_26, 
            'rv'   : rv}

#%%
# star_list = ['hd_11695','hd_18884']
# star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
#     'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
#     'hd_102870','hd_45588','hd_156098']
# star_list = ['moon','hd_18907']
star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
    'hd_10700','hd_100407']
vpass = 5
linelist = 'quinlinelist.in'
for star_name in star_list:
    #open masters stars csv
    star_info = pd.read_csv(f'/home/users/qai11/Documents/quin-masters-code/Masters_stars.csv', sep=',')
    #get the star regions
    regions = star_info[star_info['ID2'] == star_name]['regions'].apply(ast.literal_eval).values[0]
    #extract the vsini
    vsini = star_info[star_info['ID2'] == star_name]['VSINI'].values[0]
    Fe = star_info[star_info['ID2'] == star_name]['Fe'].values[0]
    CN = star_info[star_info['ID2'] == star_name]['CN'].values[0]
    CC = star_info[star_info['ID2'] == star_name]['CC'].values[0]
    
    # if star_info['LOGG'].values[0] >= 3.5: Implement
    #     type_of_star = 'dwarf'
    # else:
    #     type_of_star = 'giant'
        
    #Open summary abundances file for Mg abundance
    summary_abundances = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/{star_name}/good_lbl/summary_abundances_{star_name}.txt', sep='\s+', engine='python')
    #Extract the Mg [X/H] and error
    MgH = summary_abundances.loc[summary_abundances['element']=='Mg',['[X/H]','e[X/H]']]
    MgH = MgH['[X/H]'].values[0]
    for region in regions:
        csv_out = model_finder(star_name,linelist,region,vsini,MgH,Fe,CN,CC)
        csv_out.to_csv(f'all_fits_region_{region}_pass_{vpass}.csv')


# %%
