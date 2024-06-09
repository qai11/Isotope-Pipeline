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

#%%
#Define the path to the data and star information
# data_information = pd.read_csv(f'/home/users/qai11/Documents/quin-masters-code/Cool_stars.csv',sep=' ')
# data_information = pd.read_csv(f'/home/users/qai11/Documents/quin-masters-code/Hot_stars.csv',sep=' ')
#data_information = pd.read_csv(f'/home/users/qai11/Documents/quin-masters-code/Masters_stars.csv',sep=',')
data_information = pd.read_csv(f'/home/users/qai11/Documents/quin-masters-code/Sun.csv',sep=' ')

#Takes all stellar parameters from the csv file into a data frame
star_name = []
for i,name in enumerate(data_information['ID2']):
    name = name.lower()
    name = name[:2] + '_' + name[2:]
    star_name = np.append(star_name,name)

#%%
#Fit a lorentzian functions to the data?

#Takes the stars and performs the abundance calculation