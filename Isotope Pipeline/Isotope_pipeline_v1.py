#################
# Rapid         #
# AuTomatic     #
# Isotope       #
# Optimisation  #
#################
# From https://github.com/madeleine-mckenzie/RAtIO/blob/main/RAtIO.py
"""
Title: Isotope_pipeline_V1.py
Author: Quin Aicken Davies
Date: May 2024

Description: With Failure of Isotope_pipeline_V1.py, this script is a rework of the original pipeline.
Initially based on RATIO.py by Madeline Mckenzie, This work improves on the original with addition of automation
and extention to more isotope regions."""
#%%
#from my_imports import *
from scipy.stats import chisquare
from scipy import interpolate
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

def calc_chi(raw, r):
    """Calculate the chi squared value for the raw spectrum against the model flux."""
    # hard coded wavelength bounds
    lw, uw = get_region(r)
    # get the wavelength bounds for the region
    raw_flux = raw[(raw.wavelength > lw) & (raw.wavelength < uw)]
    # Calculate the chi squared value manually, could be used to compare to the scipy version
    # In my use case this is the only thing that worked.
    chisquare = np.sum(((raw_flux.flux - raw_flux.model_flux)**2)/ raw_flux.model_flux)
    
    return chisquare 

def interp_smooth(raw, smooth):
    """Perform an interpolation to align the smoothed model flux with the raw spectrum."""
    # perform a cubic spline on the data to make the wavelengths line up with eachother
    tck = interpolate.splrep(smooth.wavelength, smooth.flux, s=0)

    # evaluate the value of the spline at the wavelength points of the original spectra
    new_flux = interpolate.splev(raw.wavelength, tck, der=0)

    # add the model flux to the dataframe
    raw['model_flux'] = pd.Series(new_flux, index=raw.index)
    
    # return the dataframe with the new interpolated column in it 
    return raw

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
    "6            0.100000\n"                                  + \
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

def change_s(d, increase=True):
    change_by = 0.5
    ll = 4 # lower limit 
    ul = 9 # upper limit

    # if changing the values is within the limits for that parameter
    if increase and d['s'] + change_by <= ul:
        d['s'] += change_by
        return d
    elif not increase and d['s'] - change_by >= ll:
        d['s'] -= change_by
        return d
    else:
        return None

#Disabled for lbl Mg use instead
# def change_mg(d, increase=True):
#     change_by = 0.02
#     ll = -1 # lower limit 
#     ul = 1.5 # upper limit

#     # if changing the values is within the limits for that parameter
#     if increase and d['mg'] + change_by <= ul:
#         d['mg'] += change_by
#         d['mg'] = round(d['mg'], 2)
#         return d
#     elif not increase and d['mg'] - change_by >= ll:
#         d['mg'] -= change_by
#         d['mg'] = round(d['mg'], 2)
#         return d
#     else:
#         return None

def change_24(d, Mg24_step, increase=True):
    change_by = Mg24_step
    ll = 0.1 # lower limit 
    ul = 8 # upper limit

    # if changing the values is within the limits for that parameter
    if increase and d['i_24'] + change_by <= ul:
        d['i_24'] += change_by
        return d
    elif not increase and d['i_24'] - change_by >= ll:
        d['i_24'] -= change_by
        return d
    else:
        return None

def change_25(d, increase=True):
    change_by = 0.2
    ll = 1 # lower limit 
    ul = 8 # upper limit

    # if changing the values is within the limits for that parameter
    if increase and d['i_25'] + change_by <= ul:
        d['i_25'] += change_by
        d['i_25'] = round(d['i_25'],3)
        return d
    elif not increase and d['i_25'] - change_by >= ll:
        d['i_25'] -= change_by
        d['i_25'] = round(d['i_25'],3)
        return d
    else:
        return None

def change_26(d, increase=True):
    change_by = 0.4
    ll = 1 # lower limit 
    ul = 8 # upper limit

    # if changing the values is within the limits for that parameter
    if increase and d['i_26'] + change_by <= ul:
        d['i_26'] += change_by
        d['i_26'] = round(d['i_26'],3)
        return d
    elif not increase and d['i_26'] - change_by >= ll:
        d['i_26'] -= change_by
        d['i_26'] = round(d['i_26'],3)
        return d
    else:
        return None
    
#Disabled as corrected in parameters pipeline.
# def change_rv(d, increase=True):
#     change_by = 1
#     ll =-5 # lower limit 
#     ul = 5 # upper limit

#     # if changing the values is within the limits for that parameter
#     if increase and d['rv'] + change_by <= ul:
#         d['rv'] += change_by
#         d['rv'] = round(d['rv'], 2)
#         return d
#     elif not increase and d['rv'] - change_by >= ll:
#         d['rv'] -= change_by
#         d['rv'] = round(d['rv'], 2)
#         return d
#     else:
#         return None

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
    
def get_chi_squared(raw, out_filename, region, guess,vsini, make_plot = True):
    """Grab the Chi squared to compare to previous models."""
    # read in the smoothed data
    smooth = read_smoothed_spectra(out_filename, guess['rv'])
    
    raw = interp_smooth(raw, smooth)
    
    # # make a plot of the model
    # if make_plot:
    #     make_model_plots(raw, smooth, out_filename, region, guess['rv'])

    # return the chi quared value over the line
    return calc_chi(raw, region)

def make_filenames(par, prefix):
    """Creation of the filenames based on the parameters. From MM"""
    str_s = str(round(par['s'],   2)).replace('.', '')
    str_mg = str(round(par['mg'],   3)).replace('.', '')
    str_24 = str(round(par['i_24'], 3)).replace('.', '')
    str_25 = str(round(par['i_25'], 3)).replace('.', '')
    str_26 = str(round(par['i_26'], 3)).replace('.', '')
    str_rv = str(round(par['rv'],   2)).replace('.', '')

    return prefix + '_s'+ str_s +'_mg'+ str_mg + '_i' \
     + str_24 + '_' + str_25  + '_' + str_26 + '_rv' + str_rv

def get_wavelength_region(raw_wavelength,region):
    '''Try cutting out the range'''
    # lower_wavelength = raw_wavelength[0]
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
    return str(np.round(lower_wavelength, 2)) + ' ' + str(np.round(upper_wavelength, 2)) 

def optimise_model_fit(raw_spec_filename, raw_spectra, region, wavelength_region, guess,star_name,linelist,vsini,Fe,CN,CC):
    """Code that collates and runs everything then outputs the dataframe with the fit"""
    # creating the in and out filenames based on the guess parameters
    in_filename  = make_filenames(guess, 'in')
    out_filename = make_filenames(guess, 'out')

    # creates a parameter string in the directory that moog can read
    generate_parameter_string(raw_spec_filename, in_filename, out_filename, wavelength_region, guess,star_name,linelist,vsini,Fe,CN,CC)

    # create the smoothed spectra by calling pymoogi
    call_pymoogi(in_filename)

    # read in the smoothed model spectra and calculate the chi squared value
    cs = get_chi_squared(raw_spectra, out_filename, region, guess,vsini, make_plot = False)
    
    # return a dataframe with a single row (to be added to a larger df later)
    return pd.DataFrame({'filename'   : out_filename, 
                         'chi_squared': cs, 
                         's'          : guess['s'],
                         'mg'         : guess['mg'],
                         'i_24'       : guess['i_24'],
                         'i_25'       : guess['i_25'],
                         'i_26'       : guess['i_26'],
                         'rv'         : guess['rv'],
                         'ratio'      : calc_ratio(guess['i_24'], guess['i_25'], guess['i_26'])
                         }, index=[1])

def generate_neighbours(guess, region,Mg24_step):
    """Generates the neighbours of the guess parameters either side of the given parameters."""

    # a list of dictionaries
    new_guesses = []
    # Note: must pass a copy of the dictionary!!!
    new_guesses.append(change_s(guess.copy(), True))  # increase 
    new_guesses.append(change_s(guess.copy(), False)) # decrease 

    #replaced by lbl Mg abundance use. Can reenable: WARNING: doesn't always work.
    # new_guesses.append(change_mg(guess.copy(), True))  # increase 
    # new_guesses.append(change_mg(guess.copy(), False)) # decrease 

    # only optimise for these for the individual regions, not the whole thing
    if not region == -1:
        new_guesses.append(change_24(guess.copy(), Mg24_step, True))  # increase 
        new_guesses.append(change_24(guess.copy(), Mg24_step, False)) # decrease 

        new_guesses.append(change_25(guess.copy(), True))  # increase 
        new_guesses.append(change_25(guess.copy(), False)) # decrease 

        new_guesses.append(change_26(guess.copy(), True))  # increase 
        new_guesses.append(change_26(guess.copy(), False)) # decrease 

    #new_guesses.append(change_rv(guess.copy(), True))  # increase 
    #new_guesses.append(change_rv(guess.copy(), False)) # decrease 

    # Strip the none values
    new_guesses = [i for i in new_guesses if i != None]
    
    return new_guesses

def filter_guesses(guess_arr, chi_df):
    """Check to make sure you arent running a model you have already done"""
    # compare the strings of making the file name to what is currently
    no_duplicates_arr = []
    
    for dict in guess_arr:
        # make the string
        dict_name = make_filenames(dict, 'out')
        have_found = False # we havent found it yet

        for file in chi_df.filename:
            if dict_name == file:
                have_found = True
        
        # if you get throught the loop and you havent found it:
        if not have_found:
            no_duplicates_arr.append(dict)

    return no_duplicates_arr


def reconstruct_min_chi(min):
    """Reconstruct the minimum chi squared model from the dataframe row. Extra check could be removed and done elsewhere"""
    return      {'s'    : min.s, 
                 'mg'   : round(min.mg, 2), 
                 'i_24' : min.i_24, 
                 'i_25' : min.i_25, 
                 'i_26' : min.i_26, 
                 'rv'   : min.rv}

def find_minimum_neighbour(raw_spec_filename, raw_spectra, wavelength_region, region, guess, chi_df,star_name,linelist,vsini,Fe,CN,CC,Mg24_step):
    """Runs the neightbours of the guess parameters to find the minimum chi squared value. 
    Using the generate neighbours function and the filter function for checking."""
    # generate neighbours close to the guess (that havent already been run)
    guess_arr = generate_neighbours(guess, region,Mg24_step)
    print('The length of the guess array before filtering is: ', len(guess_arr))
    guess_arr = filter_guesses(guess_arr, chi_df)
    print('The length of the guess array after filtering is: ', len(guess_arr))
    
    # if there are no neighbours that we can run, this must be our minimum
    if len(guess_arr) == 0:
        return chi_df

    # run optimise_model_fit on the neighbours
    for par in guess_arr:
        # add the new chi squared values to the df
        chi_of_model = optimise_model_fit(raw_spec_filename, raw_spectra, region, wavelength_region, par,star_name,linelist,vsini,Fe,CN,CC)
        chi_df = pd.concat([chi_df,chi_of_model])
    
    # return chi_df with the results of the new models
    return chi_df

def model_finder(star_name,linelist,region,vsini,MgH,Fe,CN,CC,Mg24_step):
    """Find the model that best fits the data for a given star and region."""
    data_path = f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests_paper/'
    # change wavelength range
    region = region
    os.chdir(data_path)
    # os.system('mkdir plots')
    # initial guesses as a dictionary
    guess = initial_guess(MgH,Mg24_step)
    
    # read in the raw spectra
    raw_spec_filename = data_path + f'{star_name}_5100-5200.txt'
    raw_spectra       = read_raw_spectra(raw_spec_filename)
    wavelength_region = get_wavelength_region(raw_spectra.wavelength,region)

    # add the first chi_squyared value to the dataframe
    chi_df = optimise_model_fit(raw_spec_filename, raw_spectra, 
                                region, wavelength_region, guess,star_name,linelist,vsini,Fe,CN,CC)

    best_guess_chi = chi_df.chi_squared.iloc[0] # should only be 1 thing in the df atm

    # currently this is the best guess we have
    best_guess = guess
    while len(chi_df) < 300:
        # add the neighbours to the dataframe
        chi_df = find_minimum_neighbour(raw_spec_filename, raw_spectra, 
                        wavelength_region, region, best_guess, chi_df,star_name,linelist,vsini,Fe,CN,CC,Mg24_step)
        
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

def initial_guess(MgH,Mg24_step):
    s = 0
    mg = 0
    i_24 = 0
    i_25 = 0
    i_26 = 0
    rv = 0
    
    if Mg24_step == 0.5:
        # initial guess for the parameters
        s = 7.5
        mg = MgH
        i_24 = 2
        i_25 = 15
        i_26 = 13
        rv = 0
        print('Using the coarse pass as the initial guess: ', 
              f's = {s}, mg = {mg}, i_24 = {i_24}, i_25 = {i_25}, i_26 = {i_26}, rv = {rv}')
    else:
        fit_pass = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests_paper/all_fits_region_{region}_pass_{vpass}_coarse.csv', sep=',')
        # except:
        #     fit_pass = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests_paper/all_fits_region_{region}.csv', sep=',')
        
        #Create a dataframe with the name of the best fit file
        best_fit = fit_pass.loc[fit_pass['chi_squared'].idxmin()]['filename']
        
        # Extract The macroturbulence value from the filename
        s_block = best_fit.split('_')[1]  # 's41'
        s_value = float(s_block.lstrip('s')[:-1] + '.' + s_block[-1])
        
        #Extract the mg isotope values from the filename
        i_block = best_fit.split('_')[3:6]  # ['i20', '64', '130']

        # Remove 'i' from the first element and apply decimal for input back into the IS
        # All numbers get a decimal 1 digit from the right
        def convert(s):
            s = s.lstrip('i')  # remove 'i' if present
            return float(s[:-1] + '.' + s[-1])  # insert decimal one digit from end

        # Apply to all three parts
        converted = [convert(val) for val in i_block]
        
        # initial guess for the parameters
        s = s_value
        mg = MgH
        i_24 = converted[0]
        i_25 = converted[1]
        i_26 = converted[2]
        rv = 0
        print('Using the first pass best fit as the initial guess for fine: ', 
              f's = {s}, mg = {mg}, i_24 = {i_24}, i_25 = {i_25}, i_26 = {i_26}, rv = {rv}')
    # New guess after first pass.
    # s = 8.9
    # mg = MgH
    # i_24 = 0.8
    # i_25 = 15
    # i_26 = 3.6
    # rv = 0
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

#OLD INPUTS
# vpass = 1
# star_name = 'hd_157244'
# linelist = 'quinlinelist.in'
# # region = 4
# regions = [1,3,4,5,10]
# vsini = 5.4
# for region in regions:
#     csv_out = model_finder(star_name,linelist,region,vsini)
#     print(csv_out)

#     csv_out.to_csv(f'all_fits_region_{region}_pass_{vpass}.csv')

#%%
# star_list = ['hd_11695','hd_18884']
# star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
#     'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
#     'hd_102870','hd_45588','hd_156098']
# star_list = ['moon','hd_18907']

'''all stars below 5300K'''
# star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
#     'hd_10700','hd_100407'] 
'''giants which play up'''
star_list = ['hd_18884','hd_157244','hd_23249'] 
'''Test star'''
# star_list = ['hd_10700']
#Pass for testing purposes
vpass = 20
#name of the linelist it should look for
linelist = 'quinlinelist.in'
for star_name in star_list:
    #open masters stars csv which is a list of stars with regions, abudnaces and vsini
    star_info = pd.read_csv(f'/home/users/qai11/Documents/Isotope-Pipeline/Masters_stars.csv', sep=',')
    #get the star regions
    regions = star_info[star_info['ID2'] == star_name]['regions'].apply(ast.literal_eval).values[0]
    #extract the vsini
    vsini = star_info[star_info['ID2'] == star_name]['VSINI'].values[0]
    Fe = star_info[star_info['ID2'] == star_name]['Fe'].values[0]
    CN = star_info[star_info['ID2'] == star_name]['CN'].values[0]
    CC = star_info[star_info['ID2'] == star_name]['CC'].values[0]
     
    #Open summary abundances file for Mg abundance(This is a line by line magnesium abundance)
    summary_abundances = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/{star_name}/good_lbl/summary_abundances_{star_name}.txt', sep='\s+', engine='python')
    #Extract the Mg [X/H] and error
    MgH = summary_abundances.loc[summary_abundances['element']=='Mg',['[X/H]','e[X/H]']]
    MgH = MgH['[X/H]'].values[0]
    # for region in regions:
    #     csv_out = model_finder(star_name,linelist,region,vsini,MgH,Fe,CN,CC)
    #     csv_out.to_csv(f'all_fits_region_{region}_pass_{vpass}.csv')
    #Run a coarse and then fine search for the best fit
    # regions = [1]
    for region in regions:
        """Add a coarse search to find the best fit for the region. Then run the fine search after wards using the best fit coarse fit.
        This will allow for a more accurate fit to the data but will require a new variable for stepsize for Mg24."""
        # coarse search
        csv_out = model_finder(star_name,linelist,region,vsini,MgH,Fe,CN,CC,Mg24_step=0.5)
        csv_out.to_csv(f'all_fits_region_{region}_pass_{vpass}_coarse.csv')
        print('Finished coarse search for region: ', region)
    for region in regions:
        # fine search
        csv_out = model_finder(star_name,linelist,region,vsini,MgH,Fe,CN,CC,Mg24_step=0.1)
        csv_out.to_csv(f'all_fits_region_{region}_pass_{vpass}_fine.csv')     
        print('Finished fine search for region: ', region) 
# %%
#%% --------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------



