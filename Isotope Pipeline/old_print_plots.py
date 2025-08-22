"""Old printing functions for compatibility."""

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
    
def make_temp_file(filename):
    # will need these files in your directory - wont make them apparently...
    f = open(filename, "a+") 
    f.write('')
    f.close() 

def generate_parameter_string(raw_spec_filename, in_filename, out_filename, wavelength_region, par,star_name,linelist,stronglines,vsini):
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
    if stronglines != None:
        print(stronglines)
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
        "6            0.200000\n"                                  + \
        "12           " + str(par['mg']) + "\n"                     + \
        "22           0.20000\n"                                    + \
        "24           0.10000\n"                                    + \
        "isotopes      5    1\n"                                    + \
        "607.01214     0.2\n"                                       + \
        "606.01212     2.0\n"                                       + \
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
        print(stronglines)
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
        "6            0.200000\n"                                  + \
        "12           " + str(par['mg']) + "\n"                     + \
        "22           0.20000\n"                                    + \
        "24           0.10000\n"                                    + \
        "26           0.18000\n"                                    + \
        "isotopes      4    1\n"                                    + \
        "606.01212     5.0\n"                                       + \
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

def make_filenames(par, prefix):
    str_s = str(round(par['s'],   2)).replace('.', '')
    str_mg = str(round(par['mg'],   3)).replace('.', '')
    str_24 = str(round(par['i_24'], 3)).replace('.', '')
    str_25 = str(round(par['i_25'], 3)).replace('.', '')
    str_26 = str(round(par['i_26'], 3)).replace('.', '')
    str_rv = str(round(par['rv'],   2)).replace('.', '')

    return prefix + '_s'+ str_s +'_mg'+ str_mg + '_i' \
     + str_24 + '_' + str_25  + '_' + str_26 + '_rv' + str_rv
     
def optimise_model_fit(raw_spec_filename, raw_spectra, region, wavelength_region, guess,star_name,linelist,stronglines,vsini):

    # creating the in and out filenames based on the guess parameters
    in_filename  = make_filenames(guess, 'in')
    out_filename = make_filenames(guess, 'out')


    # creates a parameter string in the directory that moog can read
    generate_parameter_string(raw_spec_filename, in_filename, out_filename, wavelength_region, guess,star_name,linelist,stronglines,vsini)

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
    #SUN
    # s = 8.9 #has to be here for some reason or it breaks, seems to break move below 1.4
    # mg = 0.5
    # i_24 = 3.1
    # i_25 = 15
    # i_26 = 16.5
    # rv = 0
    #Current star
    s = 8.9
    mg = 0.16
    i_24 = 8.5
    i_25 = 35
    i_26 = 16
    rv = 0

    # s = 0.0 #has to be here for some reason or it breaks, seems to break move below 1.4
    # mg = 0.0
    # i_24 = 3.1
    # i_25 = 15
    # i_26 = 18.5
    # rv = 0

    # return the guess as a dictionary
    return {'s'    : s,
            'mg'   : mg, 
            'i_24' : i_24, 
            'i_25' : i_25, 
            'i_26' : i_26, 
            'rv'   : rv}
#0.9 for plots
mg = initial_guess()

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
    # upper_wavelength = raw_wavelength[len(raw_wavelength)-1] # -1 isnt working for some reason
    # print(str(np.round(lower_wavelength, 2)) + ' ' + str(np.round(upper_wavelength, 2)) )
    return str(np.round(lower_wavelength, 2)) + ' ' + str(np.round(upper_wavelength, 2)) 

def model_finder(star_name,linelist,region,stronglines,vsini):
    try:
        #Uni computer
        # data_path = f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests/'
        data_path = f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests/'
        os.chdir(data_path)
    except:
        #MAC
        if os.path.exists(f'/Users/quin/Desktop/2024_Data/Fixed_fits_files/{star_name}/moog_tests/'):
            data_path = f'/Users/quin/Desktop/2024_Data/Fixed_fits_files/{star_name}/moog_tests/'
            os.chdir(data_path)
        else:
            os.mkdir(f'/Users/quin/Desktop/2024_Data/Fixed_fits_files/{star_name}/moog_tests/')
            data_path = f'/Users/quin/Desktop/2024_Data/Fixed_fits_files/{star_name}/moog_tests/'
            os.chdir(data_path)
    region = region
    
    # os.system('mkdir plots')
    # initial guesses as a dictionary
    guess = initial_guess()

    raw_spec_filename = data_path + f'{star_name}_5100-5200.txt'
    raw_spectra       = read_raw_spectra(raw_spec_filename)
    wavelength_region = get_wavelength_region(raw_spectra.wavelength,region)

    # add the first chi_squyared value to the dataframe
    chi_df = optimise_model_fit(raw_spec_filename, raw_spectra, 
                                region, wavelength_region, guess,star_name,linelist,stronglines,vsini)
    
    # make_model_plots(raw, smooth, out_filename, region, guess['rv'])

star_name = 'hd_157244'
linelist = 'quinlinelist.in'
# stronglines = 'quinstronglines.in'
# stronglines = 'quinbarklem.in'
stronglines= None
region = 3
vsini = 5.4
# linelist = 'quinlist.MgH'
model_finder(star_name,linelist,region, stronglines,vsini)
# %%
'''idk what happened above but it made the output file for some reason so ill print that'''
mg24 = str(mg['i_24']).replace('.', '')
mg25 = str(mg['i_25']).replace('.', '')
mg26 = str(mg['i_26']).replace('.', '')
mg_all= str(mg['mg']).replace('.', '')
s_all = str(mg['s']).replace('.', '')
if -0.001 < mg['mg'] < 0.001 :
    mg_all = '00'
#HD 102870
# smoothed = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/hd_102870/out_s841_mg{mg_all}_i{mg24}_{mg25}_{mg26}_rv0', sep="     ", header=None, skiprows = [0,1])
# raw = pd.read_csv('/home/users/qai11/Documents/Fixed_fits_files/hd_102870/hd_102870_5100-5200.txt', sep="	", header=None)
#HD128620
try:
    #Uni Computer paths
    smoothed = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests/out_s{s_all}_mg{mg_all}_i{mg24}_{mg25}_{mg26}_rv0', sep="     ", header=None, skiprows = [0,1])
    raw = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests/{star_name}_5100-5200.txt', sep="	", header=None)
except:
    #Mac computer paths
    smoothed = pd.read_csv(f'/Users/quin/Desktop/2024_Data/Fixed_fits_files/{star_name}/moog_tests/out_s{s_all}_mg{mg_all}_i{mg24}_{mg25}_{mg26}_rv0', sep="     ", header=None, skiprows = [0,1])
    raw = pd.read_csv(f'/Users/quin/Desktop/2024_Data/Fixed_fits_files/{star_name}/moog_tests/{star_name}_5100-5200.txt', sep="	", header=None)

print(f'out_s{s_all}_mg{mg_all}_i{mg24}_{mg25}_{mg26}_rv0')
#%%
'''Interpolating the model to make it fit the observed data better'''
# apply_post_fundamental_effects(waveobs, fluxes, segments, macroturbulence = 3.0, vsini = 2.0, limb_darkening_coeff = 0.60, R=500000, vrad=(0,), verbose=0)

# segments = ispec_dir + '/Quin_segments.txt'
# fluxes = ispec.apply_post_fundamental_effects(smoothed[0], smoothed[1], segments, macroturbulence = 4.21, vsini = 1.6, limb_darkening_coeff = 0.60, R=82000, vrad=(0,), verbose=0)

'''Plotting the model and observed spectra'''
plt.figure(figsize=(8, 4))


plt.xlabel('Wavelength ($\AA$)',fontsize=14)
plt.ylabel('Norm. Flux',fontsize=14)

save = True
# save = False
# region = 10
'''Region 1,9,10'''
if region == 1:     
    plt.plot(smoothed[0], smoothed[1]-0.02)
    plt.plot(raw[0], raw[1])
    plt.legend(['Model', 'Observed'])
    # plt.xlim(5134, 5135.3)
    # plt.xlim(5130,5140)
    # plt.ylim(0.05, 1.01)
    # plt.ylim(0.75,1.01)
    #Find min wavelength
    lw, uw = get_region(region)
    cropped_flux = raw[(raw[0] > lw) & (raw[0] < uw)][1]
    min_flux = cropped_flux.min()
    max_flux = cropped_flux.max()
    plt.ylim(min_flux-0.05,1.01)
    plt.xlim(lw - 0.4, uw + 0.5)
    #Plot the box where the fitting region is
    plt.fill_between([lw, uw], min_flux - 0.01, 1, facecolor = '#CCDBFD', alpha = 0.3)
    #plot mg24 lines
    plt.axvline(x=5134.208, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
    plt.axvline(x=5134.570, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
    plt.axvline(x=5135.111, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
    #plot the mg25 lines
    plt.axvline(x=5134.295, ymin=0, linestyle='--', color='black',lw=1,alpha=0.5)
    plt.axvline(x=5134.656, ymin=0, linestyle='--', color='black',lw=1,alpha=0.5)
    plt.axvline(x=5135.160, ymin=0, linestyle='--', color='black',lw=1,alpha=0.5)
    #plot the mg26 lines
    plt.axvline(x=5134.734, ymin=0, color='black',lw=1,alpha=0.5)
    plt.axvline(x=5134.376, ymin=0, color='black',lw=1,alpha=0.5)
    plt.axvline(x=5135.24, ymin=0, color='black',lw=1,alpha=0.5)
    
'''Region 2'''
if region == 2:
    plt.plot(smoothed[0], smoothed[1])
    plt.plot(raw[0], raw[1])
    plt.legend(['Model', 'Observed'])    
    # plt.xlim(5138.59, 5138.95)
    # plt.xlim(5130,5140)
    # plt.xlim(5138,5140)
    # plt.ylim(0.89, 0.98)
    
    #Find min wavelength
    lw, uw = get_region(region)
    cropped_flux = raw[(raw[0] > lw) & (raw[0] < uw)][1]
    min_flux = cropped_flux.min()
    max_flux = cropped_flux.max()
    plt.ylim(min_flux-0.05,1.01)
    plt.xlim(lw - 0.4, uw + 0.4)
    #Plot the box where the fitting region is
    plt.fill_between([lw, uw], min_flux - 0.01, 1, facecolor = '#CCDBFD', alpha = 0.3)
    
    #mg24
    plt.axvline(x=5138.710, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
    #mg25
    plt.axvline(x=5138.768, ymin=0, linestyle='--', color='black',lw=1,alpha=0.5)
    plt.axvline(x=5138.785, ymin=0, linestyle='--', color='black',lw=1,alpha=0.5)
    #mg26
    plt.axvline(x=5138.826, ymin=0, color='black',lw=1,alpha=0.5)
    plt.axvline(x=5138.862, ymin=0, color='black',lw=1,alpha=0.5)

'''Region 3'''
if region == 3:
    plt.plot(smoothed[0], smoothed[1])
    plt.plot(raw[0], raw[1])
    plt.legend(['Model', 'Observed'])
    # plt.xlim(5140.04, 5140.46)
    # plt.xlim(5139.04, 5141.46)
    # plt.ylim(0.90, 1.01)
    
    #Find min wavelength
    lw, uw = get_region(region)
    cropped_flux = raw[(raw[0] > lw) & (raw[0] < uw)][1]
    min_flux = cropped_flux.min()
    max_flux = cropped_flux.max()
    plt.ylim(min_flux-0.05,1.01)
    plt.xlim(lw - 0.4, uw + 0.4)
    #Plot the box where the fitting region is
    plt.fill_between([lw, uw], min_flux - 0.01, 1, facecolor = '#CCDBFD', alpha = 0.3)
    
    #mg24
    # plt.axvline(x=5140.206, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
    plt.axvline(x=5140.229, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
    #mg25
    # plt.axvline(x=5140.253, ymin=0, linestyle='--', color='black',lw=1,alpha=0.5)
    plt.axvline(x=5140.286, ymin=0, linestyle='--', color='black',lw=1,alpha=0.5)
    #mg26
    # plt.axvline(x=5140.302, ymin=0, color='black',lw=1,alpha=0.5)
    plt.axvline(x=5140.359, ymin=0, color='black',lw=1,alpha=0.5)

'''Region 4'''
if region == 4:
    plt.plot(smoothed[0], smoothed[1])
    plt.plot(raw[0], raw[1])
    plt.legend(['Model', 'Observed'])
    # plt.xlim(5134.04, 5134.4)
    # plt.ylim(0.90, 1.01)
    # plt.ylim(0.4,1.01)
    
    #Find min wavelength
    lw, uw = get_region(region)
    cropped_flux = raw[(raw[0] > lw) & (raw[0] < uw)][1]
    min_flux = cropped_flux.min()
    max_flux = cropped_flux.max()
    plt.ylim(min_flux-0.05,1.01)
    plt.xlim(lw - 0.4, uw + 0.4)
    #Plot the box where the fitting region is
    plt.fill_between([lw, uw], min_flux - 0.01, 1, facecolor = '#CCDBFD', alpha = 0.3)
    
    #mg24
    plt.axvline(x=5134.208, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
    #mg25
    plt.axvline(x=5134.295, ymin=0, linestyle='--', color='black',lw=1,alpha=0.5)
    #mg26
    plt.axvline(x=5134.376, ymin=0, color='black',lw=1,alpha=0.5)
    
'''Region 5'''
if region == 5:
    plt.plot(smoothed[0], smoothed[1])
    plt.plot(raw[0], raw[1])
    plt.legend(['Model', 'Observed'])
    # plt.xlim(5134.95,5135.32)
    # plt.ylim(0.93, 1.01)
    # plt.ylim(0.4,1.01)
    
    #Find min wavelength
    lw, uw = get_region(region)
    cropped_flux = raw[(raw[0] > lw) & (raw[0] < uw)][1]
    min_flux = cropped_flux.min()
    max_flux = cropped_flux.max()
    plt.ylim(min_flux-0.05,1.01)
    plt.xlim(lw - 0.4, uw + 0.4)
    #Plot the box where the fitting region is
    plt.fill_between([lw, uw], min_flux - 0.01, 1, facecolor = '#CCDBFD', alpha = 0.3)
    
    #mg24
    # plt.axvline(x=5135.072, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
    plt.axvline(x=5135.111, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
    # plt.axvline(x=5135.178, ymin=0, linestyle='--', color='black',lw=1,alpha=0.5)
    #mg25
    plt.axvline(x=5135.160, ymin=0, linestyle='--', color='black',lw=1,alpha=0.5)
    #mg26
    plt.axvline(x=5135.240, ymin=0, color='black',lw=1,alpha=0.5)

'''Region 6'''
if region == 6:
    plt.plot(smoothed[0], smoothed[1])
    plt.plot(raw[0], raw[1])
    plt.legend(['Model', 'Observed'])
    # plt.xlim(5135.97, 5136.2)
    # plt.ylim(0.8, 1.01)
    # plt.ylim(0.4,1.01)
    
    #Find min wavelength
    lw, uw = get_region(region)
    cropped_flux = raw[(raw[0] > lw) & (raw[0] < uw)][1]
    min_flux = cropped_flux.min()
    max_flux = cropped_flux.max()
    plt.ylim(min_flux-0.05,1.01)
    plt.xlim(lw - 0.4, uw + 0.4)
    #Plot the box where the fitting region is
    plt.fill_between([lw, uw], min_flux - 0.01, 1, facecolor = '#CCDBFD', alpha = 0.3)
    
    #mg24
    # plt.axvline(x=5136.025, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
    # plt.axvline(x=5136.080, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
    plt.axvline(x=5136.123, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
    # plt.axvline(x=5136.229, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
    # plt.axvline(x=5136.258, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
    # plt.axvline(x=5136.294, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
    #mg25
    plt.axvline(x=5136.087, ymin=0, linestyle='--', color='black',lw=1,alpha=0.5)
    # plt.axvline(x=5136.318, ymin=0, linestyle='--', color='black',lw=1,alpha=0.5)
    #mg26
    plt.axvline(x=5136.144 , ymin=0, color='black',lw=1,alpha=0.5)

'''Region 7'''
if region == 7:
    plt.plot(smoothed[0], smoothed[1])
    plt.plot(raw[0], raw[1])
    plt.legend(['Model', 'Observed'])
    # plt.xlim(5136.33,5136.55)
    # plt.ylim(0.90, 1.01)
    # plt.ylim(0.4,1.01)
    
    #Find min wavelength
    lw, uw = get_region(region)
    cropped_flux = raw[(raw[0] > lw) & (raw[0] < uw)][1]
    min_flux = cropped_flux.min()
    max_flux = cropped_flux.max()
    plt.ylim(min_flux-0.05,1.01)
    plt.xlim(lw - 0.4, uw + 0.4)
    #Plot the box where the fitting region is
    plt.fill_between([lw, uw], min_flux - 0.01, 1, facecolor = '#CCDBFD', alpha = 0.3)
    43
    #mg24
    # plt.axvline(x=5136.229, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
    # plt.axvline(x=5136.258, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
    # plt.axvline(x=5136.294, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
    plt.axvline(x=5136.439, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
    #mg25
    plt.axvline(x=5136.502, ymin=0, linestyle='--', color='black',lw=1,alpha=0.5)
    #mg26
    plt.axvline(x=5136.560, ymin=0, color='black',lw=1,alpha=0.5)

'''Region 8'''
if region == 8:
    plt.plot(smoothed[0], smoothed[1])
    plt.plot(raw[0], raw[1])
    plt.legend(['Model', 'Observed'])
    # plt.xlim(5138.2,5138.62)
    # plt.ylim(0.90, 1.01)
    # plt.ylim(0.4,1.01)
    
    #Find min wavelength
    lw, uw = get_region(region)
    cropped_flux = raw[(raw[0] > lw) & (raw[0] < uw)][1]
    min_flux = cropped_flux.min()
    max_flux = cropped_flux.max()
    plt.ylim(min_flux-0.05,1.01)
    plt.xlim(lw - 0.4, uw + 0.4)
    #Plot the box where the fitting region is
    plt.fill_between([lw, uw], min_flux - 0.01, 1, facecolor = '#CCDBFD', alpha = 0.3)
    
    #mg24
    plt.axvline(x=5138.486, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
    #mg25
    plt.axvline(x=5138.427, ymin=0, linestyle='--', color='black',lw=1,alpha=0.5)
    #mg26
    plt.axvline(x=5138.501, ymin=0, color='black',lw=1,alpha=0.5)


'''Region 9'''
if region == 9:
    plt.plot(smoothed[0], smoothed[1])
    plt.plot(raw[0], raw[1])
    plt.legend(['Model', 'Observed'])
    # plt.xlim(5140.8,5141.2)
    # plt.ylim(0.80, 1.01)
    # plt.ylim(0.4,1.01)
    
    #Find min wavelength
    lw, uw = get_region(region)
    cropped_flux = raw[(raw[0] > lw) & (raw[0] < uw)][1]
    min_flux = cropped_flux.min()
    max_flux = cropped_flux.max()
    plt.ylim(min_flux-0.05,1.01)
    plt.xlim(lw - 0.4, uw + 0.4)
    #Plot the box where the fitting region is
    plt.fill_between([lw, uw], min_flux - 0.01, 1, facecolor = '#CCDBFD', alpha = 0.3)
    
    #mg24
    plt.axvline(x=5141.234, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
    #mg25
    plt.axvline(x=5141.288, ymin=0, linestyle='--', color='black',lw=1,alpha=0.5)
    #mg26
    plt.axvline(x=5141.338, ymin=0, color='black',lw=1,alpha=0.5)

'''Region 10'''
if region == 10:
    plt.plot(smoothed[0], smoothed[1])
    plt.plot(raw[0], raw[1])
    plt.legend(['Model', 'Observed'])
    # plt.xlim(5133.0,5133.4)
    # plt.ylim(0.80, 1.01)
    # plt.ylim(0.4,1.01)
    
    #Find min wavelength
    lw, uw = get_region(region)
    cropped_flux = raw[(raw[0] > lw) & (raw[0] < uw)][1]
    min_flux = cropped_flux.min()
    max_flux = cropped_flux.max()
    plt.ylim(min_flux-0.05,1.01)
    plt.xlim(lw - 0.4, uw + 0.4)
    #Plot the box where the fitting region is
    plt.fill_between([lw, uw], min_flux - 0.01, 1, facecolor = '#CCDBFD', alpha = 0.3)
    
    #mg24
    plt.axvline(x=5133.174, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
    #mg25
    plt.axvline(x=5133.231, ymin=0, linestyle='--', color='black',lw=1,alpha=0.5)
    #mg26
    plt.axvline(x=5133.292, ymin=0, color='black',lw=1,alpha=0.5)

if save == True:
    def calc_ratio(i_24, i_25, i_26):
        i24_percentage=1/(0.01*i_24)
        i25_percentage=1/(0.01*i_25)
        i26_percentage=1/(0.01*i_26)

        isotope_sum = i24_percentage + i25_percentage + i26_percentage
        print(f"sum {isotope_sum}")

        i24_ratio = (i24_percentage/isotope_sum) * 100
        i25_ratio = (i25_percentage/isotope_sum) * 100
        i26_ratio = (i26_percentage/isotope_sum) * 100

        return str(round(i24_ratio,2)) + '_' + str(round(i25_ratio,2)) + '_' + str(round(i26_ratio,2))

    mg_ratio = calc_ratio(mg['i_24'], mg['i_25'], mg['i_26'])
    print(mg_ratio)

# %%
