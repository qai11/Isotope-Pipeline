"""
Title: abundance.py
Author: Quin Aicken Davies
Date: 10/06/2024

Description: Tests for calculation of abundance using ispec extensions.
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
#Define the path to the data and star information
# data_information = pd.read_csv(f'/home/users/qai11/Documents/quin-masters-code/Cool_stars.csv',sep=' ')
# data_information = pd.read_csv(f'/home/users/qai11/Documents/quin-masters-code/Hot_stars.csv',sep=' ')
#data_information = pd.read_csv(f'/home/users/qai11/Documents/quin-masters-code/Masters_stars.csv',sep=',')
data_information = pd.read_csv(f'/home/users/qai11/Documents/quin-masters-code/Sun.csv',sep=' ')

#Takes all stellar parameters from the csv file into a data frame
# star_name = []
# for i,name in enumerate(data_information['ID2']):
#     name = name.lower()
#     name = name[:2] + '_' + name[2:]
#     star_name = np.append(star_name,name)

name = 'hd_102870'
#%%
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
#Fit a lorentzian functions to the data?
# def lorentzian(x, x0, gamma, A):
#     '''Gamma is the width of the peak, x0 is the center of the peak, A is the amplitude of the peak.'''
#     return (A/np.pi) * ((gamma/2)**2 / ((x - x0)**2 + (gamma/2)**2))

# # Example usage:
# x_data = np.linspace(-10, 10, 100)
# y_data = lorentz(x_data, 0, 1, 1) + np.random.normal(0, 0.1, 100)

# indexes = np.where((wavelength_vector > 518.18) & (wavelength_vector <= 518.58))[0]
# x_data = wavelength_vector[indexes].to_numpy()
# y_data = intensity_vector[indexes].to_numpy()



#%%
# sigma = 0.04
# amplitude = 0.099999

# fit_params = lmfit.Parameters()
# fit_params.add('center', value=518.362)
# fit_params.add('sigma', value=sigma)
# fit_params.add('amplitude', value=amplitude)
# fit_params.add('fwhm', value=2.0000000*sigma)
# fit_params.add('height', value=0.3183099*amplitude/max(1e-15, sigma))

# def residual(params, x, data):
#     x0 = params['x0'].value
#     gamma = params['gamma'].value   
#     A = params['A'].value
#     model = lorentz(x, x0, gamma, A)
#     return data - model


# lorentz_mod = LorentzianModel()

# y_offset = 1

# pars = lorentz_mod.guess(y_offset - y_data, x=x_data)
# # pars['amplitude'].set(value=0.1)
# out = lorentz_mod.fit(y_data, pars, x=x_data)

# out

def lorentz_model(norm_spec, fit_region, gamma, y_offset=1):
    '''Defines the lorentzian model for the fit region using limfit. Then fits to the 
    given spectrum'''
    #Define the normalized spectrum
    spec_line = Spectrum1D(spectral_axis=norm_spec.spectral_axis[fit_region], flux=norm_spec.flux[fit_region])
    
    #Define the lorentzian model
    lorentz_mod = LorentzianModel()
    
    #Define inital guess from the observed specturm
    pars = lorentz_mod.guess(y_offset - spec_line.flux.value, x=spec_line.spectral_axis.value)
    #might need to define gamma just making a note here
    # pars['amplitude'].set(value=0.1)
    
    #Fit the model to the spectrum
    fit_result = lorentz_mod.fit(y_offset - spec_line.flux.value, pars, x=spec_line.spectral_axis.value)

    return fit_result

def lorentz_fit(spec_norm, fit_region, y_fit_offset=1, plot=False):
    '''Uses the lorentz model to fit a given spectral line using varying gamma values using lmfit'''
    gamma_vals = np.arange(start=0, stop=1.1, step=0.01)
    chisqr_vals = np.zeros(shape=gamma_vals.shape)
    for j, gamma in enumerate(gamma_vals):
        fit_init = lorentz_model(norm_spec=spec_norm, fit_region=fit_region, y_offset=y_fit_offset,
                               gamma=gamma)

        # SAVING CHISQR RESULTS TO EMPTY ARRAY
        chisqr_vals[j] = fit_init.chisqr

    # TAKING BEST CHI-SQUARED VALUE CORRESPONDING TO THE GAMMA VALUE
    best_chisqr = min(chisqr_vals)
    best_gamma = gamma_vals[np.where(chisqr_vals == best_chisqr)[0][0]]

    # REDOING MODEL-FIT, BUT WITH THE BEST GAMMA VALUE
    fit_final = lorentz_model(norm_spec=spec_norm, fit_region=fit_region, y_offset=y_fit_offset,
                            gamma=best_gamma)
    
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

    
wave_518nm = 518.362
xlim_518nm = 518.1
xmax_518nm = 518.6

region_518nm = np.where(
        (spec.spectral_axis.value >= xlim_518nm) & (spec.spectral_axis.value <= xmax_518nm))
fit_518nm = lorentz_fit(spec_norm=spec, fit_region=region_518nm, plot=True)

fwhm_518nm = fit_518nm.params['fwhm'].value
 
    
    
    
    
    
#%%
# result = lmfit.minimize(residual, fit_params, args=(x_data, y_data))
#%%
# x_fit = np.linspace(518.0, 518.75, len(x_data))
# y_fit = lorentz(x_data, result.params['x0'].value, result.params['gamma'].value, result.params['A'].value)
#%%
# plt.plot(x_data, y_data, label='Data')
# plt.plot(x_data, out.init_fit, 'r-', label='Fit')
# plt.plot(x_data, out.best_fit, 'g-', label='Best Fit')
# plt.plot(x_data,lorentzian(x_data, 518.362, 0.06, -2.8)+1,label='Lorentzian')
# # plt.xlim(518.48, 518.49)
# plt.legend()
# plt.show()
# plt.plot
# #Takes the stars and performs the abundance calculation
# # %%
