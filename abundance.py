"""
Title: abundance.py
Author: Quin Aicken Davies
Date: 10/06/2024

Description: Tests for calculation of abundance using ispec extensions.
"""
import numpy as np

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
#Defining the vectors
wavelength_vector = data['waveobs']
intensity_vector = data['flux']
#%%
#Fit a lorentzian functions to the data?
# def lorentz(x, x0, gamma, A):
#     '''Gamma is the width of the peak, x0 is the center of the peak, A is the amplitude of the peak.'''
#     return (A/np.pi) * ((gamma/2)**2 / ((x - x0)**2 + (gamma/2)**2))

# # Example usage:
# x_data = np.linspace(-10, 10, 100)
# y_data = lorentz(x_data, 0, 1, 1) + np.random.normal(0, 0.1, 100)
indexes = np.where((wavelength_vector > 518.0) & (wavelength_vector <= 518.75))[0]
x_data = wavelength_vector[indexes].to_numpy()
y_data = intensity_vector[indexes].to_numpy()

#%%

fit_params = lmfit.Parameters()
fit_params.add('x0', value=518.35)
fit_params.add('gamma', value=0.1)
fit_params.add('A', value=0.8)

# def residual(params, x, data):
#     x0 = params['x0'].value
#     gamma = params['gamma'].value   
#     A = params['A'].value
#     model = lorentz(x, x0, gamma, A)
#     return data - model

lorentz = LorentzianModel(x_data, 518.35, 0.1, 0.8)

#%%
result = lmfit.minimize(lorentz, fit_params, args=(x_data, y_data))

result

x_fit = np.linspace(518.0, 518.75, len(x_data))
y_fit = lorentz(x_data, result.params['x0'].value, result.params['gamma'].value, result.params['A'].value)
#%%
plt.plot(x_data, y_data, label='Data')
plt.plot(x_fit, y_fit, 'r-', label='Fit')
# plt.xlim(518.0, 518.75)
plt.legend()
plt.show()
#Takes the stars and performs the abundance calculation
# %%
