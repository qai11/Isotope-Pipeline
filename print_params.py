"""
Title: print_params.py
Author: Quin Aicken Davies
Date: 08/08/24

Description: Prints the parameters and errors for each star in the final_params.txt and final_errors.txt files.
"""
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# star = ['hd_45588','hd_100407','hd_102870','hd_128620','hd_128621']
star = ['hd_45588','hd_100407','hd_102870','hd_128620','hd_128621','hd_11695','hd_146233','hd_156098','hd_157244','hd_160691','moon']
# star = ['hd_11695','hd_146233','hd_156098','hd_157244','hd_160691','moon']
params = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/parameters/final_params.txt')
errors = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/parameters/final_errors.txt')

# %%
iteration = 0
for i in star:
    print('\n')
    print(star[iteration])
    print(f'teff = {params["teff"][iteration]:.0f} pm {errors["teff"][iteration]:.0f}')
    print(f'logg = {params["logg"][iteration]:.2f} pm {errors["logg"][iteration]:.2f}')
    print(f'MH = {params["MH"][iteration]:.2f} pm {errors["MH"][iteration]:.2f}')
    print(f'vmac = {params["vmac"][iteration]:.2f} pm {errors["vmac"][iteration]:.2f}')
    print(f'vsini = {params["vsini"][iteration]:.2f} pm {errors["vsini"][iteration]:.2f}')
    iteration +=1
    
# %%
#plot the convergence of the parameters to the final values for each star
for i,star_name in enumerate(star):
    plot_params = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/iteration_parameters/{star_name}_iter_param.txt')
    # plot_params = plot_params.append(params.loc[i],ignore_index=True)
    plt.figure()
    plt.plot(plot_params['logg'],label='logg')
    plt.xlabel('Iteration')
    plt.ylabel('logg')
    plt.legend()
    plt.figure()
    plt.plot(plot_params['MH'],label='MH')
    plt.xlabel('Iteration')
    plt.ylabel('[M/H]')
    plt.legend()
    plt.figure()
    plt.plot(plot_params['teff'],label='teff')
    plt.xlabel('Iteration')
    plt.ylabel('Teff (K)')
    plt.legend()
    plt.figure()
    plt.plot(plot_params['CHI-SQUARE'],label='$\chi^2$')
    plt.xlabel('Iteration')
    plt.ylabel('$\chi^2$')
    plt.legend()

#%%
#plot the convergence of the parameters to the final values for each star onto sublplots
star = ['hd_45588']
for i,star_name in enumerate(star):
    plot_params = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/iteration_parameters/{star_name}_iter_param.txt')
    fig, axs = plt.subplots(2, 2, figsize=(8, 8))
    fig.tight_layout(pad=3)  # Add padding between subplots
    axs[0, 0].plot(plot_params['teff'],label='teff')
    axs[0, 0].set_ylabel('Teff (K)')
    
    axs[0, 1].plot(plot_params['logg'],label='logg')
    axs[0, 1].set_ylabel('logg')
    
    axs[1, 0].plot(plot_params['MH'],label='MH')
    axs[1, 0].set_xlabel('Iteration')
    axs[1, 0].set_ylabel('[M/H]')
    
    axs[1, 1].plot(plot_params['CHI-SQUARE'],label='$\chi^2$')
    axs[1, 1].set_xlabel('Iteration')
    axs[1, 1].set_ylabel('$\chi^2$')  
    
plt.savefig(f'/home/users/qai11/Documents/Masters_Figures/Method/{star_name}_convergence.png',dpi=150)
    
    
    
#%%
'''Plots all the convergence values onto one plot'''
# plt.figure()
# for i,star_name in enumerate(star):
#     plot_params = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/iteration_parameters/{star_name}_iter_param.txt')
#     plt.plot(plot_params['logg'],label='logg')
#     plt.xlabel('Iteration')
#     plt.ylabel('logg')
#     plt.legend(star)    
 
# plt.figure()
# for i,star_name in enumerate(star):
#     plot_params = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/iteration_parameters/{star_name}_iter_param.txt')
#     plt.plot(plot_params['MH'],label='MH')
#     plt.xlabel('Iteration')
#     plt.ylabel('[M/H]')
#     plt.legend(star)  
    
# plt.figure()
# for i,star_name in enumerate(star):
#     plot_params = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/iteration_parameters/{star_name}_iter_param.txt')
#     plt.plot(plot_params['teff'],label='teff')
#     plt.xlabel('Iteration')
#     plt.ylabel('Teff (K)')
#     plt.legend(star)  
    
# plt.figure()
# for i,star_name in enumerate(star):
#     plot_params = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/iteration_parameters/{star_name}_iter_param.txt')
#     plt.plot(plot_params['CHI-SQUARE'],label='$\chi^2$')
#     plt.xlabel('Iteration')
#     plt.ylabel('$\chi^2$')
#     plt.legend(star)     
        

# %%
