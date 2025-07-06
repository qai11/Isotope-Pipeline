"""
Title: Ratio_uncertainties.py
Author: nQuin Aicken Davies
Date: 05/02/25

Description: This code takes some of the ratio calulations and calculated the uncertainties.

"""
#%%
#Regenerate the spectra to get the chi squared values for the hessian

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

def generate_parameter_string(raw_spec_filename, in_filename, out_filename, wavelength_region, par,star_name,linelist,vsini,Fe,CN,CC):
    # I doubt I'll ever want to change these so initialise them here
    standard_out = 'out1'
    summary_out  = 'out2'

    # will need these files in your directory - wont make them apparently...
    make_temp_file(standard_out)
    make_temp_file(summary_out)
    make_temp_file(out_filename)
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

def get_chi_squared(raw, out_filename, region, guess,vsini, make_plot = True):
    """Calculate chi square for uncertainties"""
    # read in the smoothed data
    smooth = read_smoothed_spectra(out_filename, guess['rv'])

    raw = interp_smooth(raw, smooth)

    # return the chi squared value over the line
    return calc_chi(raw, region)

def calc_chi(raw, r):
    """Equation for calculating chi squared"""
    # hard coded wavelength bounds
    lw, uw = get_region(r)

    raw_flux = raw[(raw.wavelength > lw) & (raw.wavelength < uw)]
    
    chisquare = np.sum(((raw_flux.flux - raw_flux.model_flux)**2)/ raw_flux.model_flux)
    
    return chisquare 

def make_filenames(par, prefix):
    str_s = str(round(par['s'],   2)).replace('.', '')
    str_mg = str(round(par['mg'],   3)).replace('.', '')
    str_24 = str(round(par['i_24'], 3)).replace('.', '')
    str_25 = str(round(par['i_25'], 3)).replace('.', '')
    str_26 = str(round(par['i_26'], 3)).replace('.', '')
    str_rv = str(round(par['rv'],   2)).replace('.', '')

    return prefix + '_s'+ str_s +'_mg'+ str_mg + '_i' \
     + str_24 + '_' + str_25  + '_' + str_26 + '_rv' + str_rv
     
def optimise_model_fit(raw_spec_filename, raw_spectra, region, wavelength_region, guess,star_name,linelist,vsini,Fe,CN,CC):

    # creating the in and out filenames based on the guess parameters
    in_filename  = make_filenames(guess, 'in')
    out_filename = make_filenames(guess, 'out')


    # creates a parameter string in the directory that moog can read
    generate_parameter_string(raw_spec_filename, in_filename, out_filename, wavelength_region, guess,star_name,linelist,vsini,Fe,CN,CC)

    # create the smoothed spectra by calling pymoogi
    smoothed_spectrum = call_pymoogi(in_filename)
    print(smoothed_spectrum)
    call_pymoogi(in_filename)

    # read in the smoothed model spectra and calculate the chi squared value
    cs = get_chi_squared(raw_spectra, out_filename, region, guess,vsini, make_plot = True)
    # cs = None
    
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
    
def initial_guess(guess_params,MgH):

    #Current star
    s = guess_params['s']
    mg = MgH
    i_24 = guess_params['i_24']
    i_25 = guess_params['i_25']
    i_26 = guess_params['i_26']
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

def model_finder(star_name,linelist,region,vsini,guess_params,Fe,CN,CC,MgH):
    try:
        #Uni computer
        # data_path = f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests_paper/'
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
    guess = initial_guess(guess_params,MgH)

    raw_spec_filename = data_path + f'{star_name}_5100-5200.txt'
    raw_spectra       = read_raw_spectra(raw_spec_filename)
    wavelength_region = get_wavelength_region(raw_spectra.wavelength,region)

    # add the first chi_squyared value to the dataframe
    chi_df = optimise_model_fit(raw_spec_filename, raw_spectra, 
                                region, wavelength_region, guess,star_name,linelist,vsini,Fe,CN,CC)
    
    # make_model_plots(raw, smooth, out_filename, region, guess['rv'])
    best_guess_chi = chi_df.chi_squared.iloc[0] # should only be 1 thing in the df atm
    
    df_save = pd.DataFrame(columns=['filename', 'chi_squared', 's', 'mg', 'i_24', 'i_25', 'i_26', 'rv'])
    df_save = df_save.append(chi_df)
    print(df_save)
    return df_save

# star_name = 'hd_156098'
# linelist = 'quinlinelist.in'
# # stronglines = 'quinstronglines.in'
# # stronglines = 'quinbarklem.in'
# stronglines= None
# region = 6
# vsini = 9.2
# linelist = 'quinlist.MgH'

def compute_hessian_from_saved(params, saved_models, step_sizes, chi_sq,Fe,CN,CC,MgH, max_step_size=0.01):
    """
    Compute the Hessian matrix using pre-saved chi-square values for different parameter sets.
    
    Parameters:
    - params: List or array of 5 current parameter values [s, mg, i24, i25, i26]
    - saved_models: Dictionary mapping parameter tuples to chi-square values
    - step_sizes: List (or dictionary) specifying step sizes for each parameter
    - max_step_size: Default step size if not provided for a parameter
    
    Returns:
    - Hessian matrix: A 5x5 NumPy array
    """
    if len(params) != 4:
        raise ValueError("Expected 4 parameters (s, i24, i25, i26), but received {}".format(len(params)))

    num_params = 4
    hessian = np.zeros((num_params, num_params))

    # Ensure step_sizes is a list of length 5
    if isinstance(step_sizes, dict):
        step_sizes = [step_sizes.get(i, max_step_size) for i in range(num_params)]

    # Compute second derivatives numerically
    for i in range(num_params):
        for j in range(i, num_params):  # Use symmetry (Hessian is symmetric)
            step_i = step_sizes[i]
            step_j = step_sizes[j]

            # Base chi-square value
            chi2_base = saved_models.get(tuple(params), None)
            if chi2_base is None:
                continue  # Skip if base chi-square is not available

            # Increment all parameters
            params_up = params.copy()
            params_up[i] += step_i
            params_up[j] += step_j
            #generate the model for the chi squrared value
            model_up = model_finder(star_name,linelist,region,vsini, params_up,Fe,CN,CC,MgH)
            #Add the chi squared value to the model
            chi2_up = model_up['chi_squared'].iloc[0]

            # Decrement all parameters
            params_down = params.copy()
            params_down[i] -= step_i
            params_down[j] -= step_j
            #generate the model for the chi squrared value
            model_down = model_finder(star_name,linelist,region,vsini, params_down,Fe,CN,CC,MgH)
            #Add the chi squared value to the model
            chi2_down = model_down['chi_squared'].iloc[0]
            # chi2_down = saved_models.get(tuple(params_down), None)

            # Increment only the first parameter
            params_i_up = params.copy()
            params_i_up[i] += step_i
            #generate the model for the chi squrared value
            model_i_up = model_finder(star_name,linelist,region,vsini, params_i_up,Fe,CN,CC,MgH)
            #Add the chi squared value to the model
            chi2_i_up = model_i_up['chi_squared'].iloc[0]
            # chi2_i_up = saved_models.get(tuple(params_i_up), None)


            #Increment only the second parameter
            params_j_up = params.copy()
            params_j_up[j] += step_j
            chi2_j_up = saved_models.get(tuple(params_j_up), None)
            #Generate the model for the chi squrared value
            model_j_up = model_finder(star_name,linelist,region,vsini, params_j_up,Fe,CN,CC,MgH)
            #Add the chi squared value to the model
            chi2_j_up = model_j_up['chi_squared'].iloc[0]

            
            # Compute second derivative if all required chi-square values exist
            if None not in (chi2_up, chi2_down, chi2_i_up, chi2_j_up):
                hessian[i, j] = (chi2_up - chi2_i_up - chi2_j_up + chi2_down) / (step_i * step_j)
                # hessian[j, i] = chi_sq*hessian[i, j]  # Symmetric matrix with chisq
                hessian[j, i] = hessian[i, j]  # Symmetric matrix
    
    return hessian

def compute_covariance_matrix(hessian):
    # Invert the Hessian matrix to get the covariance matrix
    covariance_matrix = np.linalg.inv(hessian)
    #if the matrix has negative entries, convert them to positive
    covariance_matrix = np.abs(covariance_matrix)
    parameter_uncertainties = np.sqrt(np.diag(covariance_matrix))  # Extract 1σ uncertainties
    
    return parameter_uncertainties


def calc_ratio(i_24, i_25, i_26, sigma_i24, sigma_i25, sigma_i26):
    # Calculate isotope percentages
    i24_percentage = 1 / (0.01 * i_24)
    i25_percentage = 1 / (0.01 * i_25)
    i26_percentage = 1 / (0.01 * i_26)

    # Calculate total isotope sum
    isotope_sum = i24_percentage + i25_percentage + i26_percentage

    # Calculate uncertainty in each percentage (simple approximation based on the inverse relationship)
    sigma_i24_percentage = (1 / (0.01 * i_24))**2 * sigma_i24
    sigma_i25_percentage = (1 / (0.01 * i_25))**2 * sigma_i25
    sigma_i26_percentage = (1 / (0.01 * i_26))**2 * sigma_i26

    # Calculate uncertainty in total sum (sum of uncertainties)
    sigma_sum = np.sqrt(sigma_i24_percentage**2 + sigma_i25_percentage**2 + sigma_i26_percentage**2)

    # Calculate each ratio
    i24_ratio = (i24_percentage / isotope_sum) * 100
    i25_ratio = (i25_percentage / isotope_sum) * 100
    i26_ratio = (i26_percentage / isotope_sum) * 100

    # Calculate uncertainty in each ratio
    sigma_i24_ratio = (i24_percentage / isotope_sum) * np.sqrt(
        (sigma_i24_percentage / isotope_sum)**2 + (sigma_sum / isotope_sum**2)**2)
    sigma_i25_ratio = (i25_percentage / isotope_sum) * np.sqrt(
        (sigma_i25_percentage / isotope_sum)**2 + (sigma_sum / isotope_sum**2)**2)
    sigma_i26_ratio = (i26_percentage / isotope_sum) * np.sqrt(
        (sigma_i26_percentage / isotope_sum)**2 + (sigma_sum / isotope_sum**2)**2)

    # Format results as a string
    result = (f"Mg-24 Ratio: {i24_ratio:.2f}% ± {sigma_i24_ratio:.2f}%\n"
              f"Mg-25 Ratio: {i25_ratio:.2f}% ± {sigma_i25_ratio:.2f}%\n"
              f"Mg-26 Ratio: {i26_ratio:.2f}% ± {sigma_i26_ratio:.2f}%")

    #Put the ratios into a list with the uncertainties
    ratios = [sigma_i24_ratio, sigma_i25_ratio, sigma_i26_ratio]
    
    return ratios


# step_sizes = [0.0001, 0.01, 0.1, 0.001, 0.01]  # Default step sizes for each parameter (s, mg, i24, i25, i26)
# regions = [9]
# vsini = 1.6
# star_name = 'moon'
# vpass = 1

# linelist = 'quinlinelist.in'

stronglines= None

step_sizes = [0.0001, 0.05, 0.1, 0.1]  # step sizes for each parameter (mg, i24, i25, i26)
# star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
#     'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
#     'hd_102870','hd_45588','hd_156098']
# star_list = ['hd_2151','hd_102870','hd_45588','hd_156098']
# star_list = ['hd_102870','hd_45588','hd_156098']


# star_list = ['moon','hd_18907']
star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
    'hd_10700','hd_100407']
vpass = 6


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
    
    #Open summary abundances file for Mg abundance
    summary_abundances = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/{star_name}/good_lbl/summary_abundances_{star_name}.txt', sep='\s+', engine='python')
    #Extract the Mg [X/H] and error
    MgH = summary_abundances.loc[summary_abundances['element']=='Mg',['[X/H]','e[X/H]']]
    MgH = MgH['[X/H]'].values[0]
        
        
# for star in star_name:
#Setup df to save the best fit uncertainties
    par_unc_df = pd.DataFrame(columns=['s','mg','d_s', 'i_24','d_i_24','i_25', 'd_i_25', 'i_26','d_i_26','d_R_24', 'd_R_25', 'd_R_26'])
    for region in regions:
        #Load saved chi-square values
        #Load the best fit values for the region
        fit_pass = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests_paper/all_fits_region_{region}_pass_{vpass}.csv', sep=',')
        # fit_pass = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests_paper/all_fits_region_{region}.csv', sep=',')
        #Create a dictionary mapping parameter tuples to chi-square values
        saved_models = dict(zip(fit_pass[['s', 'i_24', 'i_25', 'i_26']].itertuples(index=False), fit_pass['chi_squared']))
        #Create a dataframe with the best fit values in it
        best_params = fit_pass.loc[fit_pass['chi_squared'].idxmin()][['s', 'i_24', 'i_25', 'i_26']]
        best_ratio = fit_pass.loc[fit_pass['chi_squared'].idxmin()][['ratio']]
        best_mg = fit_pass.loc[fit_pass['chi_squared'].idxmin()]['mg']
        #split the ratio into the individual isotopes
        best_ratio = best_ratio['ratio'].split('_')
        # print(best_params['mg'])
        
        #Compute the Hessian matrix

        hessian = compute_hessian_from_saved(best_params, saved_models, step_sizes, fit_pass['chi_squared'].min(),Fe,CN,CC,MgH)

        par_unc = compute_covariance_matrix(hessian)
        
        #Calculate the ratios and uncertainties
        mg_ratio_unc = calc_ratio(best_params['i_24'], best_params['i_25'], best_params['i_26'], par_unc[1], par_unc[2], par_unc[3])
        
        #Add the best fit parameters from S, Mg, i24, i25, i26 and the mg ratio to the dataframe
        par_unc_df = par_unc_df.append({
            'region': region,
            'pass': vpass,
            's': best_params['s'],
            'd_s': par_unc[0],
            'mg': best_mg,
            # 'd_mg': par_unc[1],
            'i_24': best_params['i_24'],
            'd_i_24': par_unc[1],
            'i_25': best_params['i_25'],
            'd_i_25': par_unc[2],
            'i_26': best_params['i_26'],
            'd_i_26': par_unc[3],
            'R_24': best_ratio[0],
            'd_R_24': mg_ratio_unc[0],
            'R_25': best_ratio[1],
            'd_R_25': mg_ratio_unc[1],
            'R_26': best_ratio[2],
            'd_R_26': mg_ratio_unc[2]
        }, ignore_index=True)
    #Save the df to a csv file
    par_unc_df.to_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests_paper/par_unc_{star_name}_paper_vpass_{vpass}_2.csv')


# mg_ratio = calc_ratio(par_unc[2], par_unc[3], par_unc[4])
# print(mg_ratio)


# %%

# #Calculation for testing
"""This one doesnt work how i want it to"""
# # Raw abundances (e.g., i_24, i_25, i_26) and their uncertainties (e.g., sigma_i24, sigma_i25, sigma_i26)
# i_24, i_25, i_26 = best_params['i_24'], best_params['i_25'],best_params['i_26']  # example raw abundances
# sigma_i24, sigma_i25, sigma_i26 = par_unc[2], par_unc[3], par_unc[4]  # example uncertainties in raw abundances

# # Calculate isotope percentages
# i24_percentage = 1 / (0.01 * i_24)
# i25_percentage = 1 / (0.01 * i_25)
# i26_percentage = 1 / (0.01 * i_26)

# # Calculate total isotope sum
# isotope_sum = i24_percentage + i25_percentage + i26_percentage

# # Calculate uncertainty in each percentage (simple approximation based on the inverse relationship)
# sigma_i24_percentage = (1 / (0.01 * i_24))**2 * sigma_i24
# sigma_i25_percentage = (1 / (0.01 * i_25))**2 * sigma_i25
# sigma_i26_percentage = (1 / (0.01 * i_26))**2 * sigma_i26

# # Calculate uncertainty in total sum (sum of uncertainties)
# sigma_sum = np.sqrt(sigma_i24_percentage**2 + sigma_i25_percentage**2 + sigma_i26_percentage**2)

# # Calculate each ratio and its uncertainty
# i24_ratio = (i24_percentage / isotope_sum) * 100
# i25_ratio = (i25_percentage / isotope_sum) * 100
# i26_ratio = (i26_percentage / isotope_sum) * 100

# # Calculate uncertainty in each ratio
# sigma_i24_ratio = (i24_percentage / isotope_sum) * np.sqrt(
#     (sigma_i24_percentage / isotope_sum)**2 + (sigma_sum / isotope_sum**2)**2)
# sigma_i25_ratio = (i25_percentage / isotope_sum) * np.sqrt(
#     (sigma_i25_percentage / isotope_sum)**2 + (sigma_sum / isotope_sum**2)**2)
# sigma_i26_ratio = (i26_percentage / isotope_sum) * np.sqrt(
#     (sigma_i26_percentage / isotope_sum)**2 + (sigma_sum / isotope_sum**2)**2)

# # Print results
# print(f"Mg-24 Ratio: {i24_ratio:.2f}% ± {sigma_i24_ratio:.2f}%")
# print(f"Mg-25 Ratio: {i25_ratio:.2f}% ± {sigma_i25_ratio:.2f}%")
# print(f"Mg-26 Ratio: {i26_ratio:.2f}% ± {sigma_i26_ratio:.2f}%")

#%%
#further testing
def convert_uncertainties_to_percentage(sigma, i_24, i_25, i_26, i24_ratio, i25_ratio, i26_ratio):
    # Extract uncertainties from covariance matrix (square root of diagonal)
    sigma_values = np.sqrt(np.diag(sigma))

    # Convert uncertainties to percentage form
    sigma_i24_percentage = (sigma_values[0] / i_24) * i24_ratio
    sigma_i25_percentage = (sigma_values[1] / i_25) * i25_ratio
    sigma_i26_percentage = (sigma_values[2] / i_26) * i26_ratio

    return sigma_i24_percentage, sigma_i25_percentage, sigma_i26_percentage

# Example input (you should replace these with your actual values)
sigma = np.linalg.inv(hessian) # Example covariance matrix
sigma = np.abs(sigma)
sigma = sigma[-3:, -3:]
i_24, i_25, i_26 = best_params['i_24'], best_params['i_25'],best_params['i_26']  # Raw values
i24_ratio, i25_ratio, i26_ratio = 85.23, 3.41, 11.36  # Already computed ratios in percentage

# Convert uncertainties
sigma_i24_perc, sigma_i25_perc, sigma_i26_perc = convert_uncertainties_to_percentage(
    sigma, i_24, i_25, i_26, i24_ratio, i25_ratio, i26_ratio)

print(f"Mg-24 Ratio Uncertainty: {sigma_i24_perc:.2f}%")
print(f"Mg-25 Ratio Uncertainty: {sigma_i25_perc:.2f}%")
print(f"Mg-26 Ratio Uncertainty: {sigma_i26_perc:.2f}%")
# %%
