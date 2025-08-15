"""
Title: Weighted_avg_models.py
Author: Quin Aicken Davies
Date: 13/08/25

Description: Generates just the weighted average models for plotting using the isotope pipeline
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
import sys 
import os

#--- iSpec directory -------------------------------------------------------------
if os.path.exists('/home/users/qai11/iSpec_v20201001'):
    "Location of the files on Uni computer"
    ispec_dir = '/home/users/qai11/iSpec_v20201001'
else:
    "location of data on Mac"
    ispec_dir = '/Users/quin/Desktop/2024_Data/iSpec_v20230804'
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec

# from ispec import apply_post_fundamental_effects

def call_pymoogi(filename):
    #https://github.com/madamow/pymoogi/blob/master/README.md
    os.system('echo q | pymoogi ' + filename)

def interp_smooth(raw, smooth):

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
    f = open(filename, "a+") 
    f.write('')
    f.close() 

def generate_parameter_string(raw_spec_filename, in_filename, out_filename, wavelength_region, par,star_name,linelist,vsini,Fe,CN,CC,stronglines):
    # I doubt I'll ever want to change these so initialise them here
    standard_out = 'out1'
    summary_out  = 'out2'

    # will need these files in your directory - wont make them apparently...
    make_temp_file(standard_out)
    make_temp_file(summary_out)
    make_temp_file(out_filename)
    #par inputs the percentage ratios of the isotopes
    if stronglines != None:
        # print(stronglines)
        par_string = "synth\n" +\
        "standard_out   '" + standard_out +"'\n"                    + \
        "summary_out    '" + summary_out +"'\n"                     + \
        "smoothed_out   '" + out_filename +"'\n"                    + \
        f"model_in       '{star_name}_atmosphere.moog'\n"                          + \
        f"lines_in       '{linelist}'\n"                           + \
        f"stronglines_in '{stronglines}'\n"                           + \
        "observed_in    '" + raw_spec_filename +"'\n"               + \
        "atmosphere    1\n"                                         + \
        "molecules     2\n"                                         + \
        "lines         2\n"                                         + \
        "strong        1\n"                                         + \
        "flux/int      0\n"                                         + \
        "plotpars      1\n"                                         + \
        wavelength_region + " 0.15 1.05\n"                          + \
        str(par['rv']) + "      0.000   0.000    1.00\n"                   + \
        "d          0.06 "+str(vsini)+" 0.6 "+ str(par['s']) +" 0.0\n"        + \
        "abundances   4    1\n"                                     + \
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
        
        #  smooth-type  FWHM-Gauss  vsini     LimbDarkeningCoeff    FWHM-Macro     FWHM-Loren
        
    elif stronglines == None:
        # print(stronglines)
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
        # "abundances   4    1\n"                                     + \
        # "6            0.000001\n"                                  + \
        # "12           " + str(par['mg']) + "\n"                     + \
        # "22           0.20000\n"                                    + \
        # "24           0.10000\n"                                    + \
        # "isotopes      5    1\n"                                    + \
        # "607.01214     8.0\n"                                       + \
        # "606.01212     2.0\n"                                       + \
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
        "isotopes      4    1\n"                                    + \
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

def make_filenames(par, prefix):
    str_s = str(round(par['s'],   2)).replace('.', '')
    str_mg = str(round(par['mg'],   3)).replace('.', '')
    str_24 = str(round(par['i_24'], 3)).replace('.', '')
    str_25 = str(round(par['i_25'], 3)).replace('.', '')
    str_26 = str(round(par['i_26'], 3)).replace('.', '')
    str_rv = str(round(par['rv'],   2)).replace('.', '')

    return prefix + '_s'+ str_s +'_mg'+ str_mg + '_i' \
     + str_24 + '_' + str_25  + '_' + str_26 + '_rv' + str_rv
     
def optimise_model_fit(raw_spec_filename, raw_spectra, region, wavelength_region, guess,star_name,linelist,vsini,Fe,CN,CC,stronglines):

    # creating the in and out filenames based on the guess parameters
    in_filename  = make_filenames(guess, 'in')
    out_filename = make_filenames(guess, 'out')
    
    # creates a parameter string in the directory that moog can read
    # generate_parameter_string(raw_spec_filename, in_filename, out_filename, wavelength_region, guess,star_name,linelist,stronglines,vsini)
    generate_parameter_string(raw_spec_filename, in_filename, out_filename, wavelength_region, guess,star_name,linelist,vsini,Fe,CN,CC,stronglines)
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
    
def initial_guess(vpass,star_name,MgH):
    # #Current star
    # s = 9
    # mg = -0.09
    # i_24 = 1.16
    # i_25 = 14.99
    # i_26 = 12.96
    # rv = 0
    # import the values from the weighted average file
    # Open the weighted average file
    weighted_avg_file = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/weighted_avg_iso_abund_paper_vpass_{vpass}.csv', sep=',')
    
    # initial guess for the parameters by matching the star and values
    s = weighted_avg_file[weighted_avg_file['Unnamed: 0'] == star_name]['s'].values[0]
    mg = MgH
    i_24 = weighted_avg_file[weighted_avg_file['Unnamed: 0'] == star_name]['i_24'].values[0]
    i_25 = weighted_avg_file[weighted_avg_file['Unnamed: 0'] == star_name]['i_25'].values[0]
    i_26 = weighted_avg_file[weighted_avg_file['Unnamed: 0'] == star_name]['i_26'].values[0]
    rv = 0
    # return the guess as a dictionary
    return {'s'    : s,
            'mg'   : mg, 
            'i_24' : i_24, 
            'i_25' : i_25, 
            'i_26' : i_26, 
            'rv'   : rv}
    
#0.9 for plots
# mg = initial_guess()

def get_wavelength_region(r, as_string=False, synth_width=8.0):
    """
    Returns a tuple of (synth_lower, synth_upper), or string if as_string=True.
    Synth range is widened around the central region to ensure MOOG can synthesize properly.
    """
    # Original analysis regions (used for residuals)
    regions = {
        0: (5128.0, 5146.45),
        1: (5134.42, 5134.85),
        2: (5138.55, 5138.95),
        3: (5140.04, 5140.46),
        4: (5134.0,  5134.4),
        5: (5134.9,  5135.3),
        6: (5135.9,  5136.3),
        7: (5136.2,  5136.6),
        8: (5138.2,  5138.6),
        9: (5141.0,  5141.45),
        10: (5133.0, 5133.4),
    }

    if r not in regions:
        print('wavelength region error')
        lw, uw = 0.0, 1.0
    else:
        lw, uw = regions[r]

    # Calculate center of original region
    center = (lw + uw) / 2.0
    half_synth = synth_width / 2.0

    # Define synthesis region
    synth_lw = center - half_synth
    synth_uw = center + half_synth

    if as_string:
        return f"{synth_lw:.2f} {synth_uw:.2f}"
    else:
        return synth_lw, synth_uw

def model_finder(star_name,linelist,region, stronglines,vsini, MgH,Fe,CN,CC,vpass):
    try:
        #Uni computer
        # data_path = f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests/'
        data_path = f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests_paper/'
        os.chdir(data_path)
    except:
        #MAC
        if os.path.exists(f'/Users/quin/Desktop/2024_Data/Fixed_fits_files/{star_name}/moog_tests_paper/'):
            data_path = f'/Users/quin/Desktop/2024_Data/Fixed_fits_files/{star_name}/moog_tests_paper/'
            os.chdir(data_path)
        else:
            os.mkdir(f'/Users/quin/Desktop/2024_Data/Fixed_fits_files/{star_name}/moog_tests_paper/')
            data_path = f'/Users/quin/Desktop/2024_Data/Fixed_fits_files/{star_name}/moog_tests_paper/'
            os.chdir(data_path)
    region = region
    
    # os.system('mkdir plots')
    # initial guesses as a dictionary
    guess = initial_guess(vpass,star_name,MgH)

    raw_spec_filename = data_path + f'{star_name}_5100-5200.txt'
    raw_spectra       = pd.read_table(raw_spec_filename, sep="\s+", usecols=[0,1], 
                         header=0, names = ['wavelength', 'flux'])
    wavelength_region = get_wavelength_region(region, as_string=True)
    

    # add the first chi_squyared value to the dataframe
    chi_df = optimise_model_fit(raw_spec_filename, raw_spectra, region, wavelength_region, guess,star_name,linelist,vsini,Fe,CN,CC,stronglines)
    
    return chi_df['filename'].values[0]
    # make_model_plots(raw, smooth, out_filename, region, guess['rv'])
#%%
# stronglines= None
# region = 1
# vsini = 0.4
# # linelist = 'quinlist.MgH'
# model_finder(star_name,linelist,region, stronglines=None,vsini)

#%%
import ast
stronglines= None
'''all stars below 5300K'''
star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
    'hd_10700','hd_100407'] 
# star_list = ['hd_10700']
vpass = 24
linelist = 'quinlinelist.in'
#define the weighted average df for plotting
w_avg_files = pd.DataFrame(columns=['star_name','filename'])
#Open the weighted average file to inset the ratio
for star_name in star_list:
    #open masters stars csv which is a list of stars with regions, abudnaces and vsini
    star_info = pd.read_csv(f'/home/users/qai11/Documents/Isotope-Pipeline/Masters_stars.csv', sep=',')
    #get the star regions 
    # regions = star_info[star_info['ID2'] == star_name]['regions'].apply(ast.literal_eval).values[0]
    region = 0
    #extract the vsini
    vsini = star_info[star_info['ID2'] == star_name]['VSINI'].values[0]
    Fe = star_info[star_info['ID2'] == star_name]['Fe'].values[0]
    CN = star_info[star_info['ID2'] == star_name]['CN'].values[0]
    CC = star_info[star_info['ID2'] == star_name]['CC'].values[0]

    #Open summary abundances file for Mg abundance(This is a line by line magnesium abundance)
    summary_abundances = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/{star_name}/good_lbl/summary_abundances_{star_name}.txt', sep='\s+', engine='python')
    #Extract the Mg [X/H] and error
    MgH = summary_abundances.loc[summary_abundances['element']=='Mg',['[X/H]','e[X/H]']]
    #The solar stuff: https://www.aanda.org/articles/aa/pdf/2021/09/aa40445-21.pdf#page=21.70
    MgH = MgH['[X/H]'].values[0] 

    #open the 
    
    """"region is a list of regions to run the model finder on"""
    file_out = model_finder(star_name,linelist,region, stronglines,vsini, MgH,Fe,CN,CC,vpass)
    
    #appen each loop to the w_avg_files dataframe
    w_avg_files = w_avg_files.append({'star_name':star_name,
                                    'filename': file_out,
                                    }, ignore_index=True)
    
    print(f"Star {star_name} region {region} done")
    #Create a df which adds the outfile name to the csv_out dataframe
    
#save the df 
w_avg_files.to_csv(f'/home/users/qai11/Documents/Fixed_fits_files/w_avg_models_vpass_{vpass}.csv', index=False)

        
print(f'avg files {w_avg_files}')
        
# %%
