"""
Title: Voigt_abundance.py
Author: Quin Aicken Davies
Date: 10/06/2024

Description: Tests for calculation of abundance using ispec extensions using a voight model.
"""

#%%
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
import numpy as np
import logging
import multiprocessing
from multiprocessing import Pool
from scipy.optimize import curve_fit
import lmfit
from lmfit.models import LorentzianModel
from specutils import Spectrum1D
import astropy.units as u
from lmfit.models import VoigtModel

#%%
# #Define the path to the data and star information
# # data_information = pd.read_csv(f'/home/users/qai11/Documents/quin-masters-code/Cool_stars.csv',sep=' ')
# # data_information = pd.read_csv(f'/home/users/qai11/Documents/quin-masters-code/Hot_stars.csv',sep=' ')
# #data_information = pd.read_csv(f'/home/users/qai11/Documents/quin-masters-code/Masters_stars.csv',sep=',')
# data_information = pd.read_csv(f'/home/users/qai11/Documents/quin-masters-code/Sun.csv',sep=' ')

# #Takes all stellar parameters from the csv file into a data frame
# # star_name = []
# # for i,name in enumerate(data_information['ID2']):
# #     name = name.lower()
# #     name = name[:2] + '_' + name[2:]
# #     star_name = np.append(star_name,name)


#%%
name = 'hd_102870'

#import the data from the .txt files
path = f'/home/users/qai11/Documents/Fixed_fits_files/{name}/J0354027.txt'  
# Define the data
data = pd.read_csv(path, sep='	')

#Defining the vectors and putting into a spectrum array
wavelength_vector = data['waveobs'].to_numpy()
intensity_vector = data['flux'].to_numpy()
#Define the spectrum with units nano metres and Jansky (standard flux unit)
spec = Spectrum1D(spectral_axis=wavelength_vector * u.nm, flux=intensity_vector * u.Jy)

#%%
def voigt_model(norm_spec, fit_region, sigma, y_offset=1):
    '''Defines the voigt model for the fit region using limfit. Then fits to the 
    given spectrum
    
    Parameters
    ----------
    norm_spec : Spectrum1D
    The normalized spectrum to be fitted
    
    fit_region : array
    The region of the spectrum to be fitted
    
    sigma : float
    The width of the line to be fitted
    
    y_offset : float
    The offset of the line to be fitted
    
    Returns
    -------
    fit_result : lmfit.model.ModelResult
    The result of the fit to the spectrum using the lorentzian model
    '''
    #Define the normalized spectrum
    spec_line = Spectrum1D(spectral_axis=norm_spec.spectral_axis[fit_region], flux=norm_spec.flux[fit_region])
    
    #Define the lorentzian model
    voigt_mod = VoigtModel()
    
    #Define inital guess from the observed specturm
    pars = voigt_mod.guess(y_offset - spec_line.flux.value, x=spec_line.spectral_axis.value)
    pars['sigma'].set(value=sigma, vary=True, expr='')
    
    #Fit the model to the spectrum
    fit_result = voigt_mod.fit(y_offset - spec_line.flux.value, pars, x=spec_line.spectral_axis.value)

    return fit_result

def voigt_fit(spec_norm, fit_region, y_fit_offset=1, plot=False):
    '''Uses the voigt model to fit a given spectral line using varying sigma values using lmfit
    
    PARAMETERS
    ----------
    sigma_vals : array
    sigma is the width of the line, and is varied to find the best fit
    
    chisqr_vals : array
    chisqr values for each sigma value Used to find the best fit for sigma
    
    Returns
    -------
    fit_final : lmfit.model.ModelResult
    The result of the fit to the spectrum using the lorentzian model
    '''
    sigma_vals = np.arange(start=0, stop=1.1, step=0.01)
    chisqr_vals = np.zeros(shape=sigma_vals.shape)
    for j, sigma in enumerate(sigma_vals):
        fit_init = voigt_model(norm_spec=spec_norm, fit_region=fit_region, y_offset=y_fit_offset,
                               sigma=sigma)

        # SAVING CHISQR RESULTS TO EMPTY ARRAY
        chisqr_vals[j] = fit_init.chisqr

    # TAKING BEST CHI-SQUARED VALUE CORRESPONDING TO THE sigma VALUE
    best_chisqr = min(chisqr_vals)
    best_sigma = sigma_vals[np.where(chisqr_vals == best_chisqr)[0][0]]

    # REDOING MODEL-FIT, BUT WITH THE BEST sigma VALUE
    fit_final = voigt_model(norm_spec=spec_norm, fit_region=fit_region, y_offset=y_fit_offset,
                            sigma=best_sigma)
    
    print(fit_final.fit_report())
    
    if plot == True:
        plt.figure()
        plt.plot(spec_norm.spectral_axis[fit_region], spec_norm.flux[fit_region], label="Observed")
        plt.plot(spec_norm.spectral_axis[fit_region], 1-fit_final.best_fit, label="Lorentz Fit")
        plt.xlabel("Wavelength (Å)", fontsize=11, fontweight="bold")
        plt.ylabel("Normalized Flux", fontsize=11, fontweight="bold")
        plt.minorticks_on()
        plt.legend()
        plt.show()
        
    return fit_final


#TEST ON MG LINE at 518.3nm    
wave_518nm = 518.3604
xlim_518nm = 518.2
xmax_518nm = 518.5

region_518nm = np.where(
        (spec.spectral_axis.value >= xlim_518nm) & (spec.spectral_axis.value <= xmax_518nm))
fit_518nm = voigt_fit(spec_norm=spec, fit_region=region_518nm, plot=True)

fwhm_518nm = fit_518nm.params['fwhm'].value
 

# %%
# ---COMBINING DATA INTO ARRAYS--- #
# width_L = np.array([fwhm_4303AA,fwhm_4520AA]) # All fwhm for all lines with large z 
# wave_L = np.array([wave_4303AA, wave_4520AA]) # All rest wavelengths for all lines with large z
# z_L = np.array([z_4303AA,z_4520AA]) # All large z values

# width_S = np.array([fwhm_4491AA])
# wave_S = np.array([wave_4491AA])
# z_S = np.array([z_4491AA])