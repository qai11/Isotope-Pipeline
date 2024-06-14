"""
Title: Lorentz_abundance.py
Author: Quin Aicken Davies
Date: 10/06/2024

Description: Tests for calculation of abundance using ispec extensions and a lorentz model.
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
def lorentz_model(norm_spec, fit_region, sigma, y_offset=1):
    '''Defines the lorentzian model for the fit region using limfit. Then fits to the 
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
    lorentz_mod = LorentzianModel()
    
    #Define inital guess from the observed specturm
    pars = lorentz_mod.guess(y_offset - spec_line.flux.value, x=spec_line.spectral_axis.value)
    pars['sigma'].set(value=sigma, vary=True, expr='')
    
    #Fit the model to the spectrum
    fit_result = lorentz_mod.fit(y_offset - spec_line.flux.value, pars, x=spec_line.spectral_axis.value)

    return fit_result

def lorentz_fit(spec_norm, fit_region, y_fit_offset=1, plot=False):
    '''Uses the lorentz model to fit a given spectral line using varying sigma values using lmfit
    
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
        fit_init = lorentz_model(norm_spec=spec_norm, fit_region=fit_region, y_offset=y_fit_offset,
                               sigma=sigma)

        # SAVING CHISQR RESULTS TO EMPTY ARRAY
        chisqr_vals[j] = fit_init.chisqr

    # TAKING BEST CHI-SQUARED VALUE CORRESPONDING TO THE sigma VALUE
    best_chisqr = min(chisqr_vals)
    best_sigma = sigma_vals[np.where(chisqr_vals == best_chisqr)[0][0]]

    # REDOING MODEL-FIT, BUT WITH THE BEST sigma VALUE
    fit_final = lorentz_model(norm_spec=spec_norm, fit_region=fit_region, y_offset=y_fit_offset,
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

#%%
#TEST ON MG LINE at 516.73nm    
wave_5167AA = 516.73
xmin_5167AA = 516.657
xmax_5167AA = 516.8

region_5167AA = np.where(
        (spec.spectral_axis.value >= xmin_5167AA) & (spec.spectral_axis.value <= xmax_5167AA))
fit_5167AA = lorentz_fit(spec_norm=spec, fit_region=region_5167AA, plot=True)

fwhm_5167AA = fit_5167AA.params['fwhm'].value

#%%
#TEST ON MG LINE at 516.9nm    
wave_5169AA = 516.90
xmin_5169AA = 516.84
xmax_5169AA = 516.96

region_5169AA = np.where(
        (spec.spectral_axis.value >= xmin_5169AA) & (spec.spectral_axis.value <= xmax_5169AA))
fit_5169AA = lorentz_fit(spec_norm=spec, fit_region=region_5169AA, plot=True)

fwhm_5169AA = fit_5169AA.params['fwhm'].value

#%%
#TEST ON MG LINE at 517.269nm    
wave_5172AA = 517.269
xmin_5172AA = 517.18
xmax_5172AA = 517.35

region_5172AA = np.where(
        (spec.spectral_axis.value >= xmin_5172AA) & (spec.spectral_axis.value <= xmax_5172AA))
fit_5172AA = lorentz_fit(spec_norm=spec, fit_region=region_5172AA, plot=True)

fwhm_5172AA = fit_5172AA.params['fwhm'].value


#%%
#TEST ON MG LINE at 518.3nm    
wave_5183AA = 518.3604
xmin_5183AA = 518.2
xmax_5183AA = 518.5

region_5183AA = np.where(
        (spec.spectral_axis.value >= xmin_5183AA) & (spec.spectral_axis.value <= xmax_5183AA))
fit_5183AA = lorentz_fit(spec_norm=spec, fit_region=region_5183AA, plot=True)

fwhm_5183AA = fit_5183AA.params['fwhm'].value

#%%
#TEST ON MG LINE at 513.4nm    
wave_5134AA = 513.458
xmin_5134AA = 513.43
xmax_5134AA = 513.47

region_5134AA = np.where(
        (spec.spectral_axis.value >= xmin_5134AA) & (spec.spectral_axis.value <= xmax_5134AA))
fit_5134AA = lorentz_fit(spec_norm=spec, fit_region=region_5134AA, plot=True)

fwhm_5134AA = fit_5134AA.params['fwhm'].value
 

# %%
# ---COMBINING DATA INTO ARRAYS--- #
# width_L = np.array([fwhm_4303AA,fwhm_4520AA]) # All fwhm for all lines with large z 
# wave_L = np.array([wave_4303AA, wave_4520AA]) # All rest wavelengths for all lines with large z
# z_L = np.array([z_4303AA,z_4520AA]) # All large z values

# width_S = np.array([fwhm_4491AA])
# wave_S = np.array([wave_4491AA])
# z_S = np.array([z_4491AA])